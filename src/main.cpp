// The simulator: sightline scan -> star draws -> light curves -> detection -> Fisher forecast -> output.
#include "common.h"
#include "types.h"
#include "run/config.h"
#include "run/inputs.h"
#include "run/sightlines.h"
#include "run/outputs.h"
#include "run/summary.h"
#include "surveys/footprints.h"
#include "surveys/schedule.h"
#include "sim/state.h"
#include "sim/sightline.h"
#include "sim/draw.h"
#include "sim/observe.h"
#include "sim/characterize.h"
#include "sim/record.h"

time_t _timeNow;
unsigned int _randVal;
unsigned int _dummyVal;
FILE * _randStream;

///==============================================================//
///                                                              //                                                    /
///                  Main program                                //
///                                                              //
///==============================================================//

int main(int argc, char** argv) {
    // NOTE: srand(time(0)) used to be called here. Nothing in this project ever
    // calls rand() -- the RNG is the seeded mt19937_64 in Bulge.h -- so it did
    // nothing except make the run look clock-seeded, which is the opposite of
    // the reproducibility the provenance block below is for.

    RunConfig cfg;
    int exitCode = 0;

    // ---- Stage 1: command line ----
    if (!parseCommandLine(argc, argv, cfg, exitCode)) return exitCode;

    // ---- Stage 2: the GBTDS detector layout, and the grid steps it bounds ----
    const GbtdsLayout gl = loadGbtdsLayout();
    GridSteps steps;
    if (int rc = resolveGridSteps(cfg, gl, steps)) return rc;

    // --------------------- Allocate objects ------------------------
    auto s  = std::make_unique<source>();
    auto l  = std::make_unique<lens>();
    auto as = std::make_unique<astromet>();
    auto cm = std::make_unique<CMD>();
//    auto ga = std::make_unique<galactic>();
    auto ex = std::make_unique<extin>();
    auto ls = std::make_unique<lsst>();
    auto ro = std::make_unique<roman>();
    auto co = std::make_unique<covarian>();
    // Step H3's second forecast: the same event with the satellite offset zeroed. Allocated
    // once beside `co` rather than per event -- covarian owns several vectors, and building
    // one per detection would cost more than the Fisher call it serves.
    auto coNS = std::make_unique<covarian>();

    // Step H1: Roman's observer position. satScale multiplies L2_OFFSET_AU inside
    // lightcurve(), so 0 puts Roman back at the centre of the Earth -- the pre-H1 behaviour,
    // and the "off" run of Step H3's experiment.
    as->satScale = cfg.noSatPar ? 0.0 : 1.0;

    // ---- Stage 3: read the input files ----
    RomanSchedule sched;
    if (int rc = loadInputs(*ls, *ro, *ex, *cm, sched)) return rc;

    // ---- Stage 4: open the output files and write their headers ----
    RunOutputs outs;
    if (int rc = openOutputs(cfg, outs)) return rc;

    // ---- Stage 5: the sightline list ----
    SightlineGrid grid;
    if (int rc = buildSightlines(cfg, gl, steps, *ls, *ro, grid)) return rc;

    // ---- Stage 6: run provenance, then (for --dry-run) stop before drawing a star ----
    if (int rc = writeRunProvenance(cfg, gl, steps, grid, sched, *ex)) return rc;

    if (cfg.dryRun) {
        printDryRunStrata(grid, steps);
        return 0;
    }

    // ---- Stage 7: the Monte Carlo over sightlines ----
    // One sightline at a time: set it up, draw events until its budget is met, aggregate it.
    // One event at a time: draw -> pre-select -> light curve -> detect + characterise -> record.
    // The stages live in src/sim/; the state they share is in include/sim/state.h.
    RunTotals totals;
    totals.NDetClassTE.resize(GG + 1);
    SimContext ctx{cfg, gl, sched, *s, *l, *as, *cm, *ex, *ls, *ro, *co, *coNS, outs, totals};
    SightlineState st;
    st.records.reserve(1000);   // rough upper bound on icon per field
    EventState ev;

    for (const auto& sightline : grid.scan) {
        const SightlineStart go = setupSightline(ctx, st, sightline);
        if (go == SightlineStart::Skip) continue;
        if (go == SightlineStart::Stop)  break;

        do {
            drawEvent(ctx, st, ev);
            if (preselectEvent(ctx, st)) simulateLightCurve(ctx, st, ev);
            characterizeEvent(ctx, st, ev);
            recordEvent(ctx, st, ev);
        } while ((st.icon < cfg.iconTarget or st.nlens < cfg.nlensTarget or st.nerr < cfg.nerrTarget)
                 and st.nsim < cfg.maxDraws);

        finishSightline(ctx, st);
    }

    // ---- Stage 8: end-of-run report ----
    writeRunSummary(cfg, totals, *l);

    return(0);
}
