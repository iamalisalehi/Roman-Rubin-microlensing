// The simulator: sightline scan -> star draws -> light curves -> detection -> Fisher forecast -> output.
#include "common.h"
#include "types.h"
#include "util/random.h"
#include "run/config.h"
#include "run/inputs.h"
#include "run/sightlines.h"
#include "run/outputs.h"
#include "run/summary.h"
#include "run/sample_dump.h"
#include "run/histograms.h"
#include "galaxy/density.h"
#include "galaxy/kinematics.h"
#include "events/source.h"
#include "events/lens.h"
#include "events/lightcurve.h"
#include "surveys/visits.h"
#include "surveys/footprints.h"
#include "surveys/schedule.h"
#include "surveys/noise.h"
#include "fisher/fisher.h"
#include "fisher/linalg.h"
#include <utility>   // std::pair, for the Roman field list

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
    
    std::vector<EventRecord> records;
    records.reserve(1000); // rough upper bound on icon per field

    // ---- Stage 3: read the input files ----
    RomanSchedule sched;
    if (int rc = loadInputs(*ls, *ro, *ex, *cm, sched)) return rc;

    // ---- The sightline loop's working state: per-event and per-field variables ----
    int    dclsEvent = DET_NONE;
    // Bin indices for the seven detection-efficiency axes. gg (tE) was the only one ever
    // computed; the other six were commented out here and at their call site, which is why
    // every efficiency column but tE has been a column of zeros (Deviation 46).
    int    save = 0, flagL, gg = -1, ss = 0, qq = 0, ww = 0, vv = 0, zz = 0, pp = 0;
    int    nri = -1, nde = -1, icon;
    int    nlens;// hh; // nde1, nri1,
    std::array<int, NDETCLASS> nDetClass{}; // per-field detection-taxonomy counts (DetClass, Bulge.h)
    // Run-wide totals. The joint-only class -- events neither telescope can find alone but the
    // combination can -- is the strongest evidence for the joint fit and is expected to be rare,
    // so it needs statistics pooled over every field, not per-field counts that are individually
    // too small to quote. Broken down by tE bin as well, since the whole science case is that
    // the gain is tE-dependent.
    std::array<long, NDETCLASS> NDetClassTot{};
    long nDchiMismatch = 0;   // Deviation 79: events where dchiL != dchiL_L + dchiL_R (should stay 0)
    std::vector<std::array<long, NDETCLASS>> NDetClassTE(GG + 1);
    long nSimTot = 0;
    // ndw MUST start at 0: it doubles as the bound of the per-event buffer clear
    // below, which runs before ndw is reset and therefore reads the PREVIOUS
    // event's value. On the very first event there is no previous value.
    int    gi,       ndw = 0, sq, ndd;
    int    giR,      sqR, nddR; // Roman-side cursor/count, parallel to gi/sq/ndd
    int    flag_det; // nml = 0;
    int    ndw_L, ndw_R;           // per-instrument epoch counts (ndw stays the joint/shared total)
    // Step R1. Per-event, per-survey tallies of epochs at which the two images were both
    // detectable and far enough apart. Reset with ndw_L/ndw_R below -- a counter that leaks
    // across events is the same bug the `ndw` note above guards against.
    long   nres5_L, nres20_L, nresPSF_L, nres5_R, nres20_R, nresPSF_R;
    double dsepMax_L, dsepMax_R;   // largest separation while both were detectable [mas]
    int    flag_det_L, flag_det_R; // per-instrument run-test result (flag_det stays the joint one)
    int    detL, detR, detJ;       // per-instrument / joint detection booleans; FFG[0] = detL or detR or detJ
    int    flagf,  fi; // datf1, datf2;
    double errs,   errg, minc, cade; // fel, , mind
    double fdetRubin, testL, testR; // Step B2: per-survey pre-selection (fdet retired)
    bool   rubinDetectable, romanDetectable, acceptRubin, acceptRoman;
    double mincR, cadeR, errgR; // Roman-side cadence/error tracking, parallel to minc/cade/errg
    double errsR, magnioR; // Roman's per-exposure astrometric error (errRomanA) / noisy magnitude
    double magnio, test, deltaA; // dist,
    double Astar0, As1,  As0;
    double initial;
    double trajm, trajp; //  ddf;
    double chi1,  chi2,   chi3, chi1a, chi2a, chi3a, sil, sil2;
    // Per-instrument duplicates. chi1/chi2/chi3/chi1a/chi2a/chi3a above are the JOINT
    // accumulators (fed by both branches); _L/_R are Rubin-only/Roman-only respectively.
    double chi1_L, chi2_L, chi3_L, chi1a_L, chi2a_L, chi3a_L;
    double chi1_R, chi2_R, chi3_R, chi1a_R, chi2a_R, chi3a_R, silR, sil2R;
    double dchiL, dchiP,  dchiA;
    double dchiL_L, dchiP_L, dchiA_L, dchiL_R, dchiP_R, dchiA_R;
    double flag0, flag1,  flag2;
    double flag0_L, flag1_L, flag2_L, flag0_R, flag1_R, flag2_R;
    double vs1,   vs2,    def1p,  def2p, vsave, dt,    Mpeak;
    double ErtE,  ErpiE,  ErtetE, Erml,  Erdl,  Ermul, Ermus, Eru0, Erfb, nsim;
    double nErAvg; // events the precision means are actually averaged over
    double mbase, fblend, Gamma,  Neven, EFF,   EffiD, EffiL, nerr=0.0;
//    double shib,  Efi;
   
    std::array<int, 3> FFG;
    std::array<double, M> magni, magni0;
    std::array<double, 2> tE, RE, piE, tetE, Vt, u0, Ml, opd, Dl, Ds, vl, vs;
    std::array<double, 2> numd, Struc, murel, vsn, DelT, mbs, fb, fwhm, Map, nbl, Ext;
   
    //fil=fopen("./files/MONTLMC/files/EfLMC2B.dat","r");
    //for(int i=0; i<=GG; ++i){
    //fscanf(fil,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",//18
    //&l.NdtE[i],&l.NstE[i], &l.NdMl[i],&l.NsMl[i], &l.Ndpi[i],&l.Nspi[i],   &l.Ndu0[i], &l.Nsu0[i], &l.Ndmb[i],&l.Nsmb[i], 
    //&l.Ndfb[i],&l.Nsfb[i], &l.Ndmu[i],&l.Nsmu[i], &l.Nhalo[1],&l.Nhalo[0], &l.Nself[1],&l.Nself[0]);}
    //fclose(fil); 
    //cout<<"**** File extinctionf.txt was read ****"<<endl;    

    // ---- Stage 4: open the output files and write their headers ----
    RunOutputs outs;
    if (int rc = openOutputs(cfg, outs)) return rc;

    save = 0;

    // ---- Stage 5: the sightline list ----
    SightlineGrid grid;
    if (int rc = buildSightlines(cfg, gl, steps, *ls, *ro, grid)) return rc;

    // ---- Stage 6: run provenance, then (for --dry-run) stop before drawing a star ----
    if (int rc = writeRunProvenance(cfg, gl, steps, grid, sched, *ex)) return rc;

    if (cfg.dryRun) {
        printDryRunStrata(grid, steps);
        return 0;
    }

    // Names the sightline loop below uses for what the stages above produced.
    const std::vector<Sightline>& scan = grid.scan;
    std::ofstream&     fil2  = outs.fil2;
    std::ofstream&     fil2b = outs.fil2b;
    std::ofstream&     fil3  = outs.fil3;
    std::ofstream&     filg_in = outs.filg_in;   // opened in append mode per event
    const std::string& fnLDt = outs.fnLDt;
    const std::string& testf = outs.testf;
    const std::string& fnPair = outs.fnPair;
    SampleSpec&             dumpSpec = outs.dumpSpec;
    std::vector<DumpEpoch>& dumpBuf  = outs.dumpBuf;
    long&                   dumpSeq  = outs.dumpSeq;

    // ---- Stage 7: the Monte Carlo over sightlines ----
///HHHHHHHHHHHHHHHHHHHHH Monte Carlo Simulation HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    // Index-driven rather than accumulate-and-compare. The old loop did
    // `for (lon = 0.5; lon <= 0.6; lon += dd)` and lost its last row AND last column
    // to floating-point drift -- 0.02 is not representable in binary, so after five
    // additions the accumulator sits a few ulp above the bound and the comparison
    // fails. That scan advertised 36 sightlines and ran 25. Computing each position
    // as lonMin + i*gridStep keeps the count exact and makes it predictable ahead of
    // the run, which the provenance block and the coverage guard both rely on.

    // Sightline bookkeeping. The scan region reaches 2*FoV + the field reach past the Roman
    // field centres (Deviation 69), so it deliberately extends to sky only a Roman-overlapping
    // Rubin pointing's far edge can catch; part of the grid may have no coverage at all.
    // These counters are what make that visible: a density in deg^-2 computed downstream
    // must know how much of the scanned area yielded nothing, and why.
    int nSkipNoCoverage = 0; //no Rubin AND no Roman epochs -- never entered the star loop
    int nSkipBarren     = 0; //had epochs, but produced no characterised event to aggregate
    int nCapped         = 0; //stopped by --maxdraws with its budget unmet
    int nAggregated     = 0; //reached the per-sightline aggregation block
    // Areas rather than counts, because under stratification the sightlines no longer
    // stand for equal pieces of sky and a count is not an area (Step E1).
    double areaAggregated = 0.0, areaNoCoverage = 0.0, areaBarren = 0.0;
    int lastCol = -1;
    long iScan = -1;
    for (const auto& sightline : scan) {
        iScan += 1;
        s->lon = sightline.lon;
        s->lat = sightline.lat;
        // The sky area THIS sightline stands for. Written into every event row it produces;
        // any statistic pooled across sightlines has to weight by it (see the list build).
        const double wArea = sightline.area;
        // nri/nde keep their old meaning -- longitude column index, and index within that
        // column -- which is what the map file and the event table record. Assigned rather
        // than incremented: the old loop bumped nri once per iLon whether or not any
        // sightline in that column survived the corner cut, so an incrementing counter here
        // would silently renumber the columns the moment a fully-cut column existed. Under
        // stratification the column index is a fine-grid one, so it steps by kSub between
        // coarse columns.
        nri = sightline.col;
        if (sightline.col != lastCol) { nde = -1; lastCol = sightline.col; }
        nde += 1;

        // --start-index. Resuming an interrupted run. The nri/nde bookkeeping above runs for
        // skipped sightlines too, deliberately: nde counts position within a longitude column
        // and is written into the map file, so skipping it would renumber the first partial
        // column and a resumed run would not concatenate onto the interrupted one. Only the
        // simulation is skipped, not the numbering.
        //
        // The index to resume AT is the number of sightlines entered, which is the number of
        // "NEW STEP" lines in the interrupted run's log -- this print sits above the
        // no-coverage and barren skips, so it counts every sightline the scan reached. The map
        // file counts only sightlines that produced something to aggregate and is therefore a
        // smaller number; resuming at it would redo the difference and duplicate rows.
        //
        // This is safe precisely because `scan` is built deterministically from --stride,
        // --stride-roman and --stub and does not depend on the RNG. Note that the random
        // sequence IS advanced by the sightlines a full run would have simulated, so a
        // resumed run is not bit-identical to an uninterrupted one -- it is a valid
        // continuation with a different draw sequence, not a replay.
        if (iScan < cfg.startIndex) continue;
        if (cfg.endIndex >= 0 and iScan >= cfg.endIndex) break;   // Deviation 77
        // Deviation 77: this sightline's own random stream (see sightlineSeed in Bulge.h). The note
        // above, that a resumed run is "not a replay", no longer applies: it IS one.
        rng.seed(sightlineSeed(cfg.seedBase, iScan));
        {
            cout << ">>>>>>>>>>> NEW STEP " << nde << " <<<<<<<<\t nri:  " << nri << endl;
            cout << "longtitude: " << s->lon << "\t latitude: " << s->lat << endl;

            cade = 0.0;

            // Deviation 80: on LSSTCam's active silicon for this visit's pointing and rotation
            // (was: within a 1.75-deg circle). A cheap (l, b) distance cut first.
            double slRA, slDec;
            galToIcrs(s->lon, s->lat, slRA, slDec);
            auto rubinCovers = [&](int i) {
                const double dl = s->lon - ls->l[i], db = s->lat - ls->b[i];
                if (dl * dl + db * db > (RUBIN_MAX_RADIUS + 0.1) * (RUBIN_MAX_RADIUS + 0.1)) return false;
                return onLsstCam(slRA, slDec, ls->RA[i], ls->DEC[i], ls->rot[i]);
            };
            auto romanCovers = [&](int i) {
                return inDetector(gl, ro->layout[i], s->lon - ro->l[i], s->lat - ro->b[i]);
            };
            ndd  = matchVisibleEpochs("LSST",  rubinCovers, ls->tim, Nl,      ls->ct, minc);
            nddR = matchVisibleEpochs("Roman", romanCovers, ro->tim, NlRoman, ro->ct, mincR);

            // Rubin's depth per band at this sightline (Deviation 73): the median of its matched
            // visits' own 5-sigma depths, for the pre-selection below. A band with no visit here
            // gets -inf, so it cannot count toward "detectable in >= 2 bands".
            std::array<double, 6> rubinDepthMed;
            {
                std::array<std::vector<double>, 6> d5;
                for (int k = 0; k < ndd; ++k) {
                    const int v = int(ls->ct[k]);
                    d5[ls->filter[v]].push_back(ls->sig5[v]);
                }
                for (int b = 0; b < 6; ++b) {
                    if (d5[b].empty()) { rubinDepthMed[b] = -std::numeric_limits<double>::infinity(); continue; }
                    std::nth_element(d5[b].begin(), d5[b].begin() + d5[b].size() / 2, d5[b].end());
                    rubinDepthMed[b] = d5[b][d5[b].size() / 2];
                }
            }
            cout << "ndd (LSST): "  << ndd  << "\t minc (LSST): "  << minc  << endl;
            cout << "ndd (Roman): " << nddR << "\t minc (Roman): " << mincR << endl;

            // Empty sky. Neither survey visits this sightline, so no light curve can ever
            // have a datum, no event can be detected, and nerr can never advance -- the
            // do/while below would spin forever. Skipping is not an approximation: a
            // sightline with no epochs contributes exactly zero events to the yield.
            //
            // It has to be a `continue` rather than a run that finds nothing, because the
            // aggregation block at the end of this loop asserts CHECK(numd[1] != 0.0) and
            // CHECK(nerr != 0.0) -- reaching it with an empty field aborts the whole run.
            //
            // The skipped area is NOT written to the map file, so anything converting the
            // per-sightline Neven density into a total count must use the aggregated
            // sightline count reported at the end of the run, not the grid size.
            if (ndd == 0 and nddR == 0) {
                nSkipNoCoverage += 1;
                areaNoCoverage  += wArea;
                cout << "  SKIP: no Rubin and no Roman epochs at this sightline" << endl;
                continue;
            }
 
            icon  = 0;
            nlens = 0;
            nDetClass.fill(0);
            nsim  = 0.0;
            nerr  = 0.0;
            for (int i = 0; i < Num; ++i) { s->nssim[i] = 0.0;  s->nsdet[i] = 0.0; }
            for (int i = 0; i <= GG; ++i) { l->nstE[i]  = 0.0;  l->ndtE[i]  = 0.0; }

            s->TET = (360.0 - s->lon) / RAa;///radian s.lon/RA;//
            s->FI  = s->lat / RAa;
            
            Disk_model(*s, 1);
            int sightlineIdx = nearestSightline(*ex, s->lon, s->lat);

            records.clear();
            do { //Start of visible star
                nsim += 1.0;
                func_source(*s, *cm, *ex, sightlineIdx);
                func_lens(*l, *s, *ex, sightlineIdx);
//                std::cerr << "nsim=" << nsim << "  Ds=" << s->Ds << "  mass=" << s->mass
//                          << "  nums=" << s->nums << "  Ml=" << l->Ml << "  u0=" << l->u0 << "\n";
                optical_depth(*s);

                // tE-histogram bin for THIS event, computed for every draw rather than only
                // for detected ones: nstE is the denominator of the detection efficiency, so
                // it has to count everything simulated. Previously gg was only evaluated
                // inside the detection branch and neither counter was ever incremented, so
                // ndtE stayed identically zero and EFF, Gamma and Neven with it (the run then
                // aborted on CHECK(EFF > 0.0) as soon as a field managed to complete).
                gg = FunctE(*l);
                l->nstE[gg] += 1.0;
                nSimTot += 1;

                // The other six efficiency axes, on the same principle as tE: the DENOMINATOR
                // has to count everything simulated, so the bin is computed here, before the
                // detection test, and the numerator is incremented in the detection branch
                // using these same indices. Every input is already set for this draw --
                // func_source filled s->Map and s->blend, func_lens filled u0, pirel and
                // murel -- so this is six array lookups and no new physics.
                //
                // Unlike the tE pair there is no per-sightline lowercase pair for these; the
                // N* arrays accumulate over the whole run, which is what the EfLMC writer
                // reports and why nothing resets them per sightline.
                ss = FuncMl(*l);                    // lens mass
                qq = FuncPi(*l);                    // log10 relative parallax
                ww = Funcu0(*l);                    // impact parameter
                vv = FuncMu(*l);                    // relative proper motion
                zz = FuncMb(*l, s->Map[2]);         // source baseline magnitude, r band
                pp = FuncFb(*l, s->blend[2]);       // blend fraction, r band
                l->NsMl[ss] += 1.0;
                l->Nspi[qq] += 1.0;
                l->Nsu0[ww] += 1.0;
                l->Nsmu[vv] += 1.0;
                l->Nsmb[zz] += 1.0;
                l->Nsfb[pp] += 1.0;

                s->nssim[s->nums] += 1.0;
                flagf   = 0;
                dumpBuf.clear(); //Step S1: this draw's epoch buffer. See the note on `ndw`.
                dclsEvent = DET_NONE; //DetClass for this draw; stays NONE if no light curve
                initial = 0.0;
                // (Step B2: the old single `test = RandR(0.0,1.0)` draw consumed here by
                // `test <= s->blend[2]` is gone — testL/testR are now drawn fresh right
                // before the per-survey pre-selection check, below.)

                // Clear only the prefix the PREVIOUS event dirtied -- `ndw` is not reset
                // until a few lines below, so it still holds that count here. Clearing all
                // `coun` slots (Nl + NlRoman = 306,092, times seven arrays = 2.1M writes)
                // to reset the ~2,000 an event actually uses was ~150x of pure waste per
                // draw, and got 16x more expensive when NlRoman became the real visit count.
                //
                // Safe because every slot in [0, ndw) is fully written before it is read --
                // both fill branches write all seven arrays at index ndw before incrementing
                // it -- and FisherM reads only [0, ndw). Slots past the previous ndw are
                // therefore untouched since construction, i.e. already zero.
                //
                // Kept as a prefix clear rather than deleted outright so the clean-slate
                // invariant survives: if a future edit ever advances ndw without filling
                // every array, that shows up as a zero instead of as the previous event's
                // photometry silently entering this event's Fisher matrix.
                //
                // NOTE: the untouched tail of tele[] is 0 (its constructed value), not the
                // -1 the old full clear wrote. Unobservable today -- tele is only read at
                // [0, ndw) -- but relevant if anything ever scans the whole array.
                for (int i = 0; i < ndw; ++i) {
                    l->timn[i] = 0.0;  l->magn[i] = 0.0; l->soux[i] = 0.0;  l->souy[i] = 0.0;
                    l->errm[i] = 0.0;  l->erra[i] = 0.0; l->tele[i] = -1;
                }

                ndw     = 0;   flag_det = 0;
                ndw_L   = 0;   ndw_R    = 0;
                nres5_L = 0; nres20_L = 0; nresPSF_L = 0; dsepMax_L = -1.0;
                nres5_R = 0; nres20_R = 0; nresPSF_R = 0; dsepMax_R = -1.0;
                flag_det_L = 0; flag_det_R = 0;
                flag0   = 0.0; flag1    = 0.0; flag2 = 0.0;
                flag0_L = 0.0; flag1_L  = 0.0; flag2_L = 0.0;
                flag0_R = 0.0; flag1_R  = 0.0; flag2_R = 0.0;
                chi1    = 0.0; chi2     = 0.0; chi3  = 0.0;
                chi1_L  = 0.0; chi2_L   = 0.0; chi3_L = 0.0;
                chi1_R  = 0.0; chi2_R   = 0.0; chi3_R = 0.0;
                chi1a   = 0.0; chi2a    = 0.0; chi3a = 0.0;
                chi1a_L = 0.0; chi2a_L  = 0.0; chi3a_L = 0.0;
                chi1a_R = 0.0; chi2a_R  = 0.0; chi3a_R = 0.0;
                def1p   = 0.0; s->def1c = 0.0; vsave = 0.0;
                def2p   = 0.0; s->def2c = 0.0;
                s->errM = 0.0; s->errA  = 0.0;
                dchiL   = 0.0; dchiP    = 0.0; dchiA = 0.0;
                dchiL_L = 0.0; dchiP_L  = 0.0; dchiA_L = 0.0;
                dchiL_R = 0.0; dchiP_R  = 0.0; dchiA_R = 0.0;

                dt   = 60.0;///days
                fdetRubin = 0.0;
                romanDetectable = false;

                for (int i = 0; i < M; ++i) {
                    Mpeak = s->magb[i] - 2.5 * std::log10(l->A0 * s->blend[i] + 1.0 - s->blend[i]);
//                        cout << "i=" << i << "  Mab=" << s->Mab[i] << "  Map=" << s->Map[i]
//                             << "  blend=" << s->blend[i] << "  Mpeak=" << Mpeak << endl;
                    if (i < 6) { // LSST ugrizy
                        if (Mpeak <= rubinDepthMed[i] and s->magb[i] > rubinDepthMed[i] - RUBIN_SATU_BELOW_M5)
                            fdetRubin += 1.0;
                    } else {     // i == 6, Roman F146 — single band, no ">=2 filters" bar applies
                        if (Mpeak <= thre[i] and s->magb[i] > satu[i])    romanDetectable = true;
                    }
                }
                rubinDetectable = (fdetRubin > 1.0); // at least 2 of the 6 LSST bands

                testL = RandR(0.0, 1.0);
                testR = RandR(0.0, 1.0);
                // Independent draws per survey — reusing one draw for both would correlate
                // the Rubin-accept and Roman-accept decisions for no physical reason. Each
                // draw is weighted by that survey's OWN blend fraction (Step B2), replacing
                // the old single test <= s->blend[2] (LSST r-band only, for both surveys).
                acceptRubin = rubinDetectable and (testL <= s->blend[2]);
                acceptRoman = romanDetectable and (testR <= s->blend[6]);

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
                if (acceptRubin or acceptRoman) { //detectable by Rubin (>=2 LSST bands) or by Roman (F146)
                    cout << "************** DETECTABLE!!!!!! ********" << endl;
                    s->nsdet[s->nums] += 1.0;
                    flagf = 1;
                    // (The legacy magC/datC/BHLSSTMONTS demo dump that was gated here could never
                    // fire -- save < 0 with save = 0 -- and was removed in Deviation 78, with the
                    // random draw that fed it.)
    
                    gi = 0;
                    gi = 0; giR = 0;
                    for (double tim = float(0.0 * year - 100.0 - initial);  tim < float(10.0 * year + 100.0 + initial); tim = tim + dt) {
                        // Rubin's geometry (observer on Earth). The quantities derived below and
                        // shared across both branches -- Astar0/As0/As1, vs1/vs2, def1p/def2p,
                        // trajm/trajp, magni[] -- are all in this frame; the Roman branch
                        // recomputes the ones it needs in its own frame (Step H1).
                        lightcurve(*s, *l, *as, tim, 0);
                        Astar0   = double(s->ut0 * s->ut0 + 2.0) / std::sqrt(s->ut0 * s->ut0 * (s->ut0 * s->ut0 + 4.0)); //MAgnification equation
                        s->Astar = double(s->ut  * s->ut  + 2.0) / std::sqrt(s->ut  * s->ut  * (s->ut  * s->ut  + 4.0)); //MAgnification equation
                        As0      = double(Astar0   * s->blend[2] + 1.0 - s->blend[2]);
                        As1      = double(s->Astar * s->blend[2] + 1.0 - s->blend[2]); //LSST r-band
                        vs1      = double(s->mus1 + (s->def1c - def1p) / dt); //[mas/days]
                        vs2      = double(s->mus2 + (s->def2c - def2p) / dt); //[mas/days]
                        def1p    = s->def1c; //pervious
                        def2p    = s->def2c;

                        // Computed once per timestep (not just inside the LSST branch) since both
                        // instruments' astrometric chi-square terms need it, and it only depends on
                        // s->pos1b/pos2b/pos1c/pos2c, which lightcurve() already refreshed above.
                        trajm = std::sqrt(s->pos1b * s->pos1b + s->pos2b * s->pos2b); //stright + parallax
                        trajp = std::sqrt(s->pos1c * s->pos1c + s->pos2c * s->pos2c); //stright + parallax+lensing

    
                        for (int i = 0; i < M; ++i) {
                            magni0[i] = s->magb[i] - 2.5 * std::log10(Astar0   * s->blend[i] + 1.0 - s->blend[i]);
                            magni[i]  = s->magb[i] - 2.5 * std::log10(s->Astar * s->blend[i] + 1.0 - s->blend[i]);
                        }
                        // Rubin's representative-band model magnitude (RUBIN_REF_BANDS, config/parameters.h) --
                        // replaces the old hardcoded magni[2] (r-band) at the two use sites below.
                        // Reduces to exactly magni[2] when RUBIN_REF_BANDS = {2} (the default), since
                        // s->mbs[0]/s->fb[0] were built from the same combination in func_source.
                        double magniRubinRef = s->mbs[0] - 2.5 * std::log10(s->Astar * s->fb[0] + 1.0 - s->fb[0]);
                        sq = int(ls->ct[gi]);
                        sq  = int(ls->ct[gi]);
                        sqR = (nddR > 0 and giR < nddR) ? int(ro->ct[giR]) : -1;

                        // ---------------- LSST (ugrizy) ----------------
    
                        if (tim >= 0.0 and tim <= Tobs and tim >= ls->tim[int(ls->ct[0])] and tim <= ls->tim[int(ls->ct[ndd - 1])] and
                            gi < ndd and sq >= 0 and sq <= static_cast<int>(Nl) and tim >= ls->tim[sq]) {
    
                            fi = int(ls->filter[sq]);
    
                            // Deviation 73: this visit's own depth and saturation, not the SRD
                            // minimum -- the same depth that sets errg below.
                            const double m5v   = double(ls->sig5[sq]);
                            const double satuv = m5v - RUBIN_SATU_BELOW_M5;
                            if (magni[fi] >= satuv and magni[fi] <= m5v) {
                                errg = errlsstM(magni[fi], int(fi), m5v); //[mag]
                                errs = errlsstA(*ls, magniRubinRef); ///[mas]

                                // Step R1. Could Rubin have told the two images apart at THIS
                                // epoch? Inside the magnitude gate on purpose: the paper's
                                // criterion counts recorded data points, and an epoch the
                                // survey never recorded is not one.
                                {
                                    const ImagePair ip = imagePair(s->ut, l->tetE, s->magb[fi],
                                                                   s->blend[fi], m5v, satuv);
                                    if (ip.bothDetectable) {
                                        if (ip.sep >= RESOLVE_D_FAINT  * errs)          nres5_L   += 1;
                                        if (ip.sep >= RESOLVE_D_BRIGHT * errs)          nres20_L  += 1;
                                        if (ip.sep >= FWHM[fi] * ARCSEC_TO_MAS)         nresPSF_L += 1;
                                        if (ip.sep >  dsepMax_L)                        dsepMax_L  = ip.sep;
                                    }
                                }

                                deltaA = std::fabs(std::pow(10.0, -0.4 * errg) - 1.0) * (s->blend[fi] * s->Astar + 1.0 - s->blend[fi]);
                                magnio = magni[fi] + RandN(errg, NOISE_TRUNC_NSIGMA);
    
                                chi1 += std::fabs((magnio -   magni[fi]) * (magnio -   magni[fi]) / (errg * errg)); //real
                                chi2 += std::fabs((magnio -  magni0[fi]) * (magnio -  magni0[fi]) / (errg * errg)); //real no parallax
                                chi3 += std::fabs((magnio - s->magb[fi]) * (magnio - s->magb[fi]) / (errg * errg)); //baseline
                                chi1_L += std::fabs((magnio -   magni[fi]) * (magnio -   magni[fi]) / (errg * errg));
                                chi2_L += std::fabs((magnio -  magni0[fi]) * (magnio -  magni0[fi]) / (errg * errg));
                                chi3_L += std::fabs((magnio - s->magb[fi]) * (magnio - s->magb[fi]) / (errg * errg));

                                // erra is the per-coordinate sigma (Deviation 71; was drawn and
                                // divided as if each coordinate had variance 2 errs^2).
                                sil  = RandN(errs, NOISE_TRUNC_NSIGMA);
                                sil2 = RandN(errs, NOISE_TRUNC_NSIGMA);

                                chi1a += std::fabs((trajp + sil - trajp) * (trajp + sil - trajp) / (errs * errs)); //real trajectory
                                chi2a += std::fabs((trajp + sil - trajm) * (trajp + sil - trajm) / (errs * errs)); //real without lensing
                                chi3a += std::fabs( sil2  * sil2 / (errs * errs));
                                chi1a_L += std::fabs((trajp + sil - trajp) * (trajp + sil - trajp) / (errs * errs));
                                chi2a_L += std::fabs((trajp + sil - trajm) * (trajp + sil - trajm) / (errs * errs));
                                chi3a_L += std::fabs( sil2  * sil2 / (errs * errs));

                                CHECK(ndw < coun); // array-bounds guard — see note on `coun` sizing
                                l->timn[ndw] = tim;
                                l->magn[ndw] = magniRubinRef; // Rubin representative-band model value (RUBIN_REF_BANDS)
                                l->errm[ndw] = errg;
                                l->soux[ndw] = s->pos1c;
                                l->souy[ndw] = s->pos2c;
                                l->erra[ndw] = errs;
                                l->tele[ndw] = 0; // 0 = LSST
                                l->rseas[ndw] = -1; l->rroll[ndw] = -1;

                                // Step S1. Everything here was computed above for the
                                // detection test; nothing new is drawn. magnio in
                                // particular is the noisy datum chi1/chi2/chi3 just used.
                                if (dumpSpec.on)
                                    dumpBuf.push_back(DumpEpoch{
                                        tim, 0, int(fi),
                                        magnio, magni[fi], magni0[fi], errg,
                                        s->ut, s->ut0, s->Astar, Astar0,
                                        s->def1c, s->def2c, s->def1a, s->def2a,
                                        s->pos1b, s->pos2b, s->pos1c, s->pos2c,
                                        l->pos1,  l->pos2,  errs});
    
                                flag2 = 0.0;
                                if (std::fabs(magnio - s->magb[fi]) > std::fabs(OUTLIER_FLAG_NSIGMA * errg))    flag2 = 1.0;
                                if (ndw_L > 2 and float(flag0 + flag1 + flag2) > OUTLIER_RUN_THRESHOLD)   flag_det = 1;
                                flag0 = flag1;
                                flag1 = flag2;
                                flag2_L = flag2; // identical test; kept as an explicit Rubin-labeled copy
                                if (ndw_L > 2 and float(flag0_L + flag1_L + flag2_L) > OUTLIER_RUN_THRESHOLD)   flag_det_L = 1;
                                flag0_L = flag1_L;
                                flag1_L = flag2_L;
    
    
                                CHECK(sq >= 0);
                                CHECK(sq <= int(Nl - 1));
                                CHECK(fi >= 0);
                                CHECK(fi <= 5);
                                CHECK(errs >= 0.0);
                                CHECK(errg >= 0.0);
                                CHECK(deltaA >= 0.0);
                                CHECK(gi <= ndd);
    
                                s->errM += deltaA;
                                s->errA += errs;
                                vsave += std::sqrt(vs1 * vs1 + vs2 * vs2);

                                ndw += 1;
                                ndw_L += 1;
                            }//magnitude limit
                            gi += 1;
                        }
                        // ---------------- Roman (F146) — new ----------------
                        // Roman now feeds its own chi1_R/chi2_R/chi3_R (+ astrometric _R,
                        // placeholder error — see errsR below) and the joint chi1/chi2/chi3,
                        // alongside the Rubin-only _L versions built in the LSST branch above.
                        // See Step B1 of JOINT_FIT_REFACTOR_PLAN.md. s->errM/s->errA remain
                        // Rubin-only by design (ORIENTATION.md: they're specifically the
                        // LSST-only running-error accumulators, unrelated to detection).
                        if (nddR > 0 and tim >= 0.0 and tim <= Tobs and
                            tim >= ro->tim[int(ro->ct[0])] and tim <= ro->tim[int(ro->ct[nddR - 1])] and
                            giR < nddR and sqR >= 0 and sqR <= static_cast<int>(NlRoman) and tim >= ro->tim[sqR]) {

                            constexpr int fiR = 6; // F146 — the only Roman band modeled so far

                            // ---- Step H1: switch to Roman's observer position ----
                            // Everything above was computed with the observer on Earth, which is
                            // right for Rubin and wrong for Roman by the L2 offset. Roman sees a
                            // slightly different trajectory, so the magnification, the F146
                            // magnitude and the astrometric positions all have to be rebuilt here
                            // before any of them is recorded. `magni`/`magni0`/`trajm`/`trajp` are
                            // per-timestep scratch that the Rubin branch above has already
                            // consumed and that the next timestep overwrites, so they are
                            // rewritten in place rather than shadowed.
                            //
                            // No RNG is consumed by this call, so the random stream is untouched
                            // and L2_OFFSET_AU = 0 reproduces the pre-H1 run exactly.
                            lightcurve(*s, *l, *as, tim, 1);
                            {
                                const double Astar0R = double(s->ut0 * s->ut0 + 2.0)
                                                     / std::sqrt(s->ut0 * s->ut0 * (s->ut0 * s->ut0 + 4.0));
                                s->Astar = double(s->ut * s->ut + 2.0)
                                         / std::sqrt(s->ut * s->ut * (s->ut * s->ut + 4.0));
                                magni0[fiR] = s->magb[fiR] - 2.5 * std::log10(Astar0R  * s->blend[fiR] + 1.0 - s->blend[fiR]);
                                magni[fiR]  = s->magb[fiR] - 2.5 * std::log10(s->Astar * s->blend[fiR] + 1.0 - s->blend[fiR]);
                                trajm = std::sqrt(s->pos1b * s->pos1b + s->pos2b * s->pos2b);
                                trajp = std::sqrt(s->pos1c * s->pos1c + s->pos2c * s->pos2c);
                            }

                            if (magni[fiR] >= satu[fiR] and magni[fiR] <= thre[fiR]) {
                                errgR = errRomanM(*ro, magni[fiR]); //[mag] (Deviation 72)

                                magnioR = magni[fiR] + RandN(errgR, NOISE_TRUNC_NSIGMA);
                                chi1 += std::fabs((magnioR -   magni[fiR]) * (magnioR -   magni[fiR]) / (errgR * errgR));
                                chi2 += std::fabs((magnioR -  magni0[fiR]) * (magnioR -  magni0[fiR]) / (errgR * errgR));
                                chi3 += std::fabs((magnioR - s->magb[fiR]) * (magnioR - s->magb[fiR]) / (errgR * errgR));
                                chi1_R += std::fabs((magnioR -   magni[fiR]) * (magnioR -   magni[fiR]) / (errgR * errgR));
                                chi2_R += std::fabs((magnioR -  magni0[fiR]) * (magnioR -  magni0[fiR]) / (errgR * errgR));
                                chi3_R += std::fabs((magnioR - s->magb[fiR]) * (magnioR - s->magb[fiR]) / (errgR * errgR));

                                // Step H4: Roman's own per-exposure astrometric error, replacing the
                                // errlsstA() placeholder (Rubin's curve at Roman's magnitude, which had
                                // no reason to be right). Constants and sources in config/parameters.h; the model
                                // is per EXPOSURE, which is what one row of RomanBaseline.dat is.
                                errsR = errRomanA(magni[fiR]); //[mas]

                                // Step R1, Roman side. s->ut is Roman's OWN impact parameter
                                // here: lightcurve(..., 1) rebuilt the trajectory in the L2
                                // frame above, so this is not the Rubin value reused.
                                {
                                    const ImagePair ip = imagePair(s->ut, l->tetE, s->magb[fiR],
                                                                   s->blend[fiR], thre[fiR], satu[fiR]);
                                    if (ip.bothDetectable) {
                                        if (ip.sep >= RESOLVE_D_FAINT  * errsR)         nres5_R   += 1;
                                        if (ip.sep >= RESOLVE_D_BRIGHT * errsR)         nres20_R  += 1;
                                        if (ip.sep >= FWHM[fiR] * ARCSEC_TO_MAS)        nresPSF_R += 1;
                                        if (ip.sep >  dsepMax_R)                        dsepMax_R  = ip.sep;
                                    }
                                }

                                silR  = RandN(errsR, NOISE_TRUNC_NSIGMA);
                                sil2R = RandN(errsR, NOISE_TRUNC_NSIGMA);
                                chi1a += std::fabs((trajp + silR - trajp) * (trajp + silR - trajp) / (errsR * errsR));
                                chi2a += std::fabs((trajp + silR - trajm) * (trajp + silR - trajm) / (errsR * errsR));
                                chi3a += std::fabs( sil2R * sil2R / (errsR * errsR));
                                chi1a_R += std::fabs((trajp + silR - trajp) * (trajp + silR - trajp) / (errsR * errsR));
                                chi2a_R += std::fabs((trajp + silR - trajm) * (trajp + silR - trajm) / (errsR * errsR));
                                chi3a_R += std::fabs( sil2R * sil2R / (errsR * errsR));

                                flag2_R = 0.0;
                                if (std::fabs(magnioR - s->magb[fiR]) > std::fabs(OUTLIER_FLAG_NSIGMA * errgR))    flag2_R = 1.0;
                                if (ndw_R > 2 and float(flag0_R + flag1_R + flag2_R) > OUTLIER_RUN_THRESHOLD)   flag_det_R = 1;
                                flag0_R = flag1_R;
                                flag1_R = flag2_R;

                                CHECK(ndw < coun); // array-bounds guard — see note on `coun` sizing
                                l->timn[ndw] = tim;
                                l->magn[ndw] = magni[fiR]; // F146 magnitude — NOT magni[2]; FisherM's
                                                            // tt==1 branch compares against s.mbs[1]/s.fb[1],
                                                            // which are F146-based, so the reference point
                                                            // recorded here must be F146 too.
                                l->errm[ndw] = errgR;
                                l->soux[ndw] = s->pos1c;
                                l->souy[ndw] = s->pos2c;
                                // Step H4. This used to store `errs` -- the RUBIN astrometric error,
                                // left over from whichever Rubin epoch last set it, and in general from
                                // a different timestep entirely. So the astrometric Fisher matrix was
                                // being weighted by a stale value from the other telescope, not even by
                                // the errlsstA(magni[fiR]) the comment above it described: errsR was
                                // computed for the chi-squared terms and then thrown away. Now Roman's
                                // own per-exposure error is both used and stored.
                                l->erra[ndw] = errsR;
                                l->tele[ndw] = 1; // 1 = Roman/F146
                                // Season and roll of this exposure, for the astrometric noise
                                // variants' day blocks and frame groups (Deviation 71).
                                l->rseas[ndw] = sched.seasonOf(ro->tim[sqR]);
                                l->rroll[ndw] = ro->layout[sqR];
                                CHECK(l->rseas[ndw] >= 0);

                                // Step S1, Roman side. s->ut, s->def* and s->pos* are the
                                // L2-frame values lightcurve(..., 1) rebuilt above, not the
                                // Rubin ones from earlier in this timestep. Astar0 is
                                // recomputed from s->ut0 because the Roman branch's own
                                // Astar0R went out of scope before the magnitude gate.
                                if (dumpSpec.on)
                                    dumpBuf.push_back(DumpEpoch{
                                        tim, 1, fiR,
                                        magnioR, magni[fiR], magni0[fiR], errgR,
                                        s->ut, s->ut0, s->Astar, magnifOf(s->ut0),
                                        s->def1c, s->def2c, s->def1a, s->def2a,
                                        s->pos1b, s->pos2b, s->pos1c, s->pos2c,
                                        l->pos1,  l->pos2,  errsR});

                                CHECK(sqR >= 0);
                                CHECK(sqR <= int(NlRoman - 1));
                                CHECK(errgR >= 0.0);
                                CHECK(giR <= nddR);

                                ndw += 1;
                                ndw_R += 1;
                            }//magnitude limit
                            giR += 1;
                        }

                        // ---------------- Adaptive step size ----------------
                        // dt must respect whichever instrument's *next* epoch is sooner.
                        // Sizing dt to LSST's cadence alone (the original behavior) would
                        // silently step right over Roman's much denser epochs.
                        if ( tim < -50.0 or tim > (10.0 * year + 50.0)) {
                            dt = 60.0; //days
                        }
    
                        else {
                            bool lsstInWindow  = (tim >= ls->tim[int(ls->ct[0])] and tim <= ls->tim[int(ls->ct[ndd - 1])]);
                            bool romanInWindow = (nddR > 0 and tim >= ro->tim[int(ro->ct[0])] and tim <= ro->tim[int(ro->ct[nddR - 1])]);

                            if (lsstInWindow) {
                                if (gi > 0 and gi < ndd) cade = float(ls->tim[int(ls->ct[gi])] - ls->tim[int(ls->ct[gi - 1])]); //days
                                else cade = minc;
                            } else {
                                cade = TIME_STEP_OUTSIDE_LSST_DAYS; //days — matches the original "outside LSST window" fallback
                            }

                            if (romanInWindow) {
                                if (giR > 0 and giR < nddR) cadeR = float(ro->tim[int(ro->ct[giR])] - ro->tim[int(ro->ct[giR - 1])]); //days
                                else cadeR = mincR;
                            } else {
                                cadeR = cade; // Roman not active right now — don't let it constrain dt
                            }

                            dt = double(std::min(cade, cadeR));
                        }
    
                    }//end of loop time
                }// end of visible star

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH 
                for (int i = 0; i < nq; ++i) co->resu[i] = -1.0;
                // Same sentinel discipline for the per-survey results: an event that never
                // reaches FisherM must not inherit the previous event's sigmas.
                for (int q = 0; q < NSURV; ++q) {
                    co->okA[q] = 0; co->okB[q] = 0; co->nepochA[q] = 0;
                    co->flagi = 0;   // Deviation 78: was stale on uncharacterised events
                    co->condA[q] = -1.0; co->condB[q] = -1.0;
                    co->relMl[q] = -1.0;
                    for (int k = 0; k < Nx; ++k) co->Era[q][k] = -1.0;
                    for (int k = 0; k < Ny; ++k) co->Erb[q][k] = -1.0;
                    for (int v = 0; v < NAVAR; ++v) {
                        co->okBV[v][q] = 0; co->condBV[v][q] = -1.0; co->relMlV[v][q] = -1.0;
                        for (int k = 0; k < Ny; ++k) co->ErbV[v][q][k] = -1.0;
                    }
                }
    
                FFG[0] = 0;
                FFG[1] = 0;
                FFG[2] = 0;
                detL = 0; detR = 0; detJ = 0;

                //cout << "flagf: " << flagf << "ndw: " << ndw << endl;

                if (flagf == 0 or ndw <= 2) {
                    errg    = errlsstM(s->magb[2], 2, double(RUBIN_R_DEPTH5_FALLBACK)); //r-band
                    s->errA = errlsstA(*ls, s->magb[2]); //r-band
                    s->errM = std::fabs(std::pow(10.0, - 0.4 * errg) - 1.0); //r-band
                    dchiL = 0.0;
                    dchiP = 0.0;
                    dchiA = 0.0;
                    dchiL_L = 0.0; dchiP_L = 0.0; dchiA_L = 0.0;
                    dchiL_R = 0.0; dchiP_R = 0.0; dchiA_R = 0.0;
                    vsave = s->mus;
                }

                // Step H3's second forecast for this event. Declared one scope out from the
                // detection block below, because the row is written further down, where the
                // satellite observable and the coverage counts have been computed.
                // Step H3's no-satellite forecast for this same event. Declared here, before
                // the detection block, so the row write further down can see them whether or
                // not the Fisher step ran. -1.0 is the not-measured sentinel used everywhere
                // else in this code: a sigma that never inverted is not a large sigma.
                double sigtE_ns   = -1.0, sigpiE_ns  = -1.0, sigpiER_ns = -1.0;
                double sigtetE_ns = -1.0, sigpiEb_ns = -1.0, relMl_ns   = -1.0;
                double condA_ns   = -1.0, condB_ns   = -1.0;
                int    okNS       = 0,    okNSb      = 0;

                if (flagf > 0 and ndw > 2) { //if star is visible
                    cout << "************** DETECTABLE!!!!!! ********" << endl;
                    icon +=1;
                    vsave   = double(vsave   / (ndw + 0.000065645));
                    s->errM = double(s->errM / (ndw + 0.000065645));
                    s->errA = double(s->errA / (ndw + 0.000065645));
                    // Step H7. The lensing statistic is SIGNED: chi3 is the flat-baseline
                    // residual and chi1 the lensing-model residual, so chi3 - chi1 > 0 means
                    // the lensing model fits better and is the only direction that can count
                    // as evidence of lensing. The previous fabs() meant a lensing model that
                    // fitted WORSE than the baseline by more than the bar would have been
                    // reported as a detection. That never bit in practice, because the data
                    // are lensed and chi3 > chi1 for any event with real signal, but it also
                    // breaks the monotonicity H7 exists to establish: with fabs, an instrument
                    // whose epochs all fall outside the event contributes a small NEGATIVE
                    // difference from noise, and |dchi_L + dchi_R| can then fall below
                    // max(|dchi_L|, |dchi_R|). Signed, the sum is exact and monotone.
                    //
                    // dchiP and dchiA are not detection tests (nothing thresholds them); since
                    // Deviation 78 they are signed like dchiL, so a negative value (the simpler
                    // model fitting better) is visible rather than folded into a "size".
                    dchiL   = chi3  - chi1;             //lensing_effect (signed)
                    // Deviation 78: signed like dchiL (positive = the full model fits better).
                    dchiP   = chi2  - chi1;             //parallax_effect (signed)
                    dchiA   = chi2a - chi1a;            //deflection_effect (signed)
                    dchiL_L = chi3_L - chi1_L;
                    dchiP_L = (chi2_L  - chi1_L);   // signed (Deviation 78)
                    dchiA_L = (chi2a_L - chi1a_L);   // signed (Deviation 78)
                    dchiL_R = chi3_R - chi1_R;
                    dchiP_R = (chi2_R  - chi1_R);   // signed (Deviation 78)
                    dchiA_R = (chi2a_R - chi1a_R);   // signed (Deviation 78)
                    // flag_det is the JOINT run-test flag the table reports (Deviation 78): it was set in
                    // the Rubin branch only, so a Roman-only persistent signal left it at 0.
                    flag_det = (flag_det_L > 0 or flag_det_R > 0) ? 1 : 0;

                    // Three independent detection tests. detL/detR use each instrument's own
                    // epoch count (ndw_L/ndw_R) and run-test result — mixing in the joint ndw
                    // here would let Roman's dense epochs silently raise the bar for a purely
                    // Rubin-driven signal (and vice versa). detJ requires a persistent run in
                    // EITHER instrument's own cadence (flag_det_L or flag_det_R) rather than an
                    // interleaved run across two different cadences, which wouldn't mean anything
                    // (see Step B1 teaching brief). FFG[0] — the boolean that gates FisherM below,
                    // preserving the existing detection-then-Fisher ordering (plan §0.4) — is the
                    // union of all three, so an event Rubin alone clearly detects is never dropped
                    // from characterization just because the joint-scaled threshold happens to miss.
                    if (s->FWHM < Tobs and dchiL_L > cfg.dchiDet and flag_det_L > 0 and ndw_L > 10) detL = 1;
                    if (s->FWHM < Tobs and dchiL_R > cfg.dchiDet and flag_det_R > 0 and ndw_R > 10) detR = 1;
                    if (s->FWHM < Tobs and dchiL   > cfg.dchiDet and (flag_det_L > 0 or flag_det_R > 0) and ndw > 10) detJ = 1;

                    // Monotonicity. Adding data to an analysis cannot destroy information, so
                    // if either survey alone clears its detection bar, the combined stream --
                    // which contains that survey's data in full, plus more -- must clear its
                    // own.
                    //
                    // Since Step H7 this holds BY CONSTRUCTION rather than by patch. chi1 and
                    // chi3 are accumulated over both instruments' epochs in the same loops that
                    // fill chi1_L/chi1_R and chi3_L/chi3_R, so
                    //
                    //     chi1 = chi1_L + chi1_R   and   chi3 = chi3_L + chi3_R
                    //
                    // exactly, hence dchiL = dchiL_L + dchiL_R with the signed statistic above.
                    // All three tests now use the SAME fixed bar, so if either survey's own
                    // dchi clears it, the sum -- which is that dchi plus a quantity that is
                    // positive whenever the other survey sees any signal at all -- clears it
                    // too. The auxiliary gates agree: flag_det_L > 0 implies (flag_det_L or
                    // flag_det_R), and ndw_L > 10 implies ndw = ndw_L + ndw_R > 10.
                    //
                    // The old bar was 2*ndw, a threshold on the MEAN per-epoch chi-squared
                    // improvement rather than on total significance, which made pooling a
                    // survey with many low-signal epochs RAISE the joint bar without adding
                    // signal. In Roman's footprint that lifted the joint bar ~4,728 above
                    // Roman's own and vetoed 21.3% of all detections there.
                    //
                    // detJ_raw and the monotone patch are retained deliberately (H7 acceptance
                    // criterion 4) so the anomaly rate stays MEASURABLE and can be shown to
                    // have gone to zero rather than been hidden. If DET_ANOMALY is ever
                    // non-zero again, the construction above has been broken and the counter
                    // is how that gets noticed.
                    const int detJ_raw = detJ;
                    if (detL or detR) detJ = 1;

                    if (detL or detR or detJ) { //lensing — detected by Rubin, Roman, or the joint test
                        FFG[0] = 1; //Lensing
                        nlens += FFG[0];
                        FisherM(*s, *l, *as, *co, ndw);

                        if (co->flagi > 0) {
                            // flagi is always +1: FisherM's F*F^-1 checks are commented out
                            // (Deviation 40). The conditioning test that works is okA.
                            nerr += co->okA[SJOINT] ? 1.0 : 0.0;
                            ErrorCal(*co, *l, *s);

                            std::ofstream fil0_append(fnLDt, std::ios::app);
                            fil0_append << std::fixed << std::setprecision(5)
                                        << l->Ml       << " " << l->Dl       << " " << l->mul       << " " << s->fb[0]     << " " << s->mbs[0]   << " "
          <<std::setprecision(7)           << l->tE       << " " << l->murel    << " " << l->u0        << " " << s->lon       << " " << s->lat      << " "
                                        << l->piE      << " " << l->tetE     << " " << co->resu[1]  << " " << co->resu[2]  << " " << co->resu[3] << " "
                                        << co->resu[5] << " " << co->resu[9] << " " << co->resu[10] << " " << co->resu[13] << "\n";
                            fil0_append.close();
                        }

                        // -----------------------------------------------------------------
                        // Step H3. The satellite-parallax experiment as a CONTROLLED
                        // comparison on one event, instead of a difference between two runs.
                        //
                        // Why it must be done here. The obvious design -- two runs differing
                        // only in L2_OFFSET_AU, matched by row -- cannot work. The RNG is one
                        // stream and the per-event path draws conditionally on `acceptRubin
                        // or acceptRoman`, so moving the observer changes what is detectable,
                        // the streams fork at the first footprint sightline, and the two
                        // runs' DETECTED populations then differ. Measured on the v3 pair,
                        // that selection difference shifts the median sigma_tE -- which the
                        // satellite cannot physically touch -- by 15%, several times the
                        // effect being looked for. DEVIATIONS.md 35.
                        //
                        // Here the event is fixed: same draw, same epochs, same photometry,
                        // only the observer moves. satScale is the single knob carrying
                        // L2_OFFSET_AU into lightcurve(), and FisherM rebuilds the light
                        // curve from it, so flipping it and recomputing gives the forecast
                        // THIS event would have had with Roman at Earth.
                        //
                        // Re-evaluating is legitimate for an event whose data were generated
                        // with the offset on: the Fisher matrix is built from model
                        // derivatives at the true parameters, not from realised noise. The
                        // question is what each observing geometry can constrain.
                        if (cfg.pairSat) {
                            const double keepScale = as->satScale;
                            as->satScale = 0.0;
                            FisherM(*s, *l, *as, *coNS, ndw);
                            if (coNS->flagi > 0) ErrorCal(*coNS, *l, *s);
                            as->satScale = keepScale;

                            // Both halves of the forecast are recorded, because moving the
                            // observer acts on both and they carry different physics.
                            //
                            //   PHOTOMETRIC (okA, Era): the magnification depends on |u|, and
                            //   the two observers see different |u| at the same instant. That
                            //   difference is the simultaneous parallax baseline, and it lands
                            //   on piE.
                            //
                            //   ASTROMETRIC (okB, Erb): the centroid deflection is a vector,
                            //   theta_E * u/(u^2+2), so moving the observer changes its
                            //   direction as well as its size. Erb[0] is sigma(theta_E) and
                            //   Erb[3] the astrometric route to piE, which is independent of
                            //   the photometric one.
                            //
                            // relMl combines them: Ml = theta_E/(kappa piE) is the quantity a
                            // microlensing survey actually wants, and it needs one observable
                            // from each matrix, so it is the only place where a change in
                            // either shows up as a change in the science.
                            //
                            // SROMAN as well as SJOINT: the L2 offset is Roman's geometry, so
                            // Roman's own matrix is where any effect must appear first and
                            // undiluted by Rubin's epochs.
                            //
                            // Condition numbers for both matrices and both observers, so the
                            // conditioning question can be answered from the data instead of
                            // hypothesised. This is the diagnostic OPEN_ITEMS asked for.
                            okNS  = coNS->okA[SJOINT];
                            okNSb = coNS->okB[SJOINT];
                            sigtE_ns   = okNS  ? coNS->Era[SJOINT][1] : -1.0;
                            sigpiE_ns  = okNS  ? coNS->Era[SJOINT][3] : -1.0;
                            sigpiER_ns = coNS->okA[SROMAN] ? coNS->Era[SROMAN][3] : -1.0;
                            sigtetE_ns = okNSb ? coNS->Erb[SJOINT][0] : -1.0;
                            sigpiEb_ns = okNSb ? coNS->Erb[SJOINT][3] : -1.0;
                            relMl_ns   = coNS->relMl[SJOINT];
                            condA_ns   = coNS->condA[SJOINT];
                            condB_ns   = coNS->condB[SJOINT];
                        }
                    }

                    // Detection taxonomy (DetClass, Bulge.h). The joint fit is what makes a
                    // detection meaningful -- it sees strictly more data than either survey
                    // alone -- so the classes are distinguished by which telescopes ALSO
                    // detect the event unaided. DET_ANOMALY catches the case that should be
                    // impossible, a single-telescope detection the joint test misses.
                    dclsEvent = detClass(detL, detR, detJ);
                    const int dcls = dclsEvent;
                    nDetClass[dcls]    += 1;
                    NDetClassTot[dcls] += 1;
                    NDetClassTE[gg][dcls] += 1;
                    // Diagnostic only: how often the raw joint test contradicted a
                    // single-survey detection, before monotonicity was imposed.
                    // Both counters, and the run total is the one that was missing: it is
                    // read at the end of the run to decide whether to explain DET_ANOMALY,
                    // but was never incremented, so every run has reported exactly zero
                    // anomalies since the counter was introduced -- including the 2026-08-30
                    // production run, whose "0 ANOMALY" is this bug and not a measurement.
                    // The per-sightline counter was always right, which is why the anomaly
                    // was visible in the log all along and invisible in every summary.
                    // Deviation 79: the joint delta-chi2 is accumulated separately from the two
                    // per-survey ones; it must equal their sum on every event.
                    if (std::fabs((dchiL_L + dchiL_R) - dchiL) > 1e-6 * std::max(1.0, std::fabs(dchiL))) {
                        if (nDchiMismatch < 10)
                            std::cerr << "DCHI_MISMATCH sightline " << iScan << " joint " << dchiL
                                      << " L+R " << dchiL_L + dchiL_R << "\n";
                        nDchiMismatch += 1;
                    }
                    if (!detJ_raw and (detL or detR)) {
                        nDetClass[DET_ANOMALY]    += 1;
                        NDetClassTot[DET_ANOMALY] += 1;
                        // Deviation 79: log everything needed to explain it. The joint delta-chi2
                        // is the SUM of the two surveys' (each accumulated over its own epochs), and
                        // either term can be negative through noise -- a survey whose data happen to
                        // fit a flat baseline slightly better than the true model -- so a single-
                        // survey detection just over the bar can leave the sum just under it. That
                        // is physics, not a sign bug; "sum - joint" must be ~0 if the bookkeeping is
                        // right, and a non-zero value there WOULD be a bug.
                        std::cerr << std::setprecision(10)
                                  << "DET_ANOMALY_DETAIL sightline " << iScan << " lon " << s->lon
                                  << " lat " << s->lat << " | dchiL_L " << dchiL_L << " dchiL_R "
                                  << dchiL_R << " sum " << dchiL_L + dchiL_R << " joint " << dchiL
                                  << " (sum - joint " << (dchiL_L + dchiL_R) - dchiL << ") bar "
                                  << cfg.dchiDet << " | ndw_L " << ndw_L << " ndw_R " << ndw_R
                                  << " ndw " << ndw << " | flag_det_L " << flag_det_L
                                  << " flag_det_R " << flag_det_R << " | detL " << detL << " detR "
                                  << detR << " | tE " << l->tE << " u0 " << l->u0 << " t0 " << l->t0
                                  << "\n";
                    }

                    // Detected-event count for this tE bin. The joint detection is the one
                    // that defines "detected" here, per the taxonomy above; the per-class
                    // breakdown lives in nDetClass. gg was computed for every simulated event
                    // before the detection test, so nstE (the denominator) counts all draws
                    // and this counts the numerator -- which is what makes EFF an efficiency.
                    if (detJ) l->ndtE[gg] += 1.0;

                    // The numerators for the other six axes, using the bin indices computed
                    // for this draw before the detection test. Same definition of "detected"
                    // as the tE curve -- the joint test -- so all seven efficiencies describe
                    // the same thing and can be read side by side.
                    if (detJ) {
                        l->NdMl[ss] += 1.0;
                        l->Ndpi[qq] += 1.0;
                        l->Ndu0[ww] += 1.0;
                        l->Ndmu[vv] += 1.0;
                        l->Ndmb[zz] += 1.0;
                        l->Ndfb[pp] += 1.0;

                        // Nhalo and Nself were declared, written out, and never counted, so
                        // the EfLMC writer divided 0 by 0 and every row of every EfLMC file
                        // ever produced ends in two "-nan" columns. Counted here, with the
                        // definitions stated rather than guessed:
                        //   Nhalo[1]/Nhalo[0]  fraction of detections whose LENS is a halo
                        //                      star -- the population a MACHO search cares
                        //                      about, and negligible toward the bulge.
                        //   Nself[1]/Nself[0]  fraction that are bulge self-lensing: bulge
                        //                      lens AND bulge source, the dominant channel
                        //                      here and the one whose kinematics set tE.
                        l->Nhalo[0] += 1.0;
                        l->Nself[0] += 1.0;
                        if (l->struc == GalacticComponent::HALO) l->Nhalo[1] += 1.0;
                        if (l->struc == GalacticComponent::BULGE and
                            s->struc == GalacticComponent::BULGE) l->Nself[1] += 1.0;
                    }
                }
//            CHECK(gg >= 0);

            // Deviation 76: the observed peak, and the gap geometry measured from it.
            const auto [t0obs, uminObs] = observedPeak(*s, *l, *as);

            records.push_back(EventRecord{
                icon, static_cast<int>(FFG[0]),
                l->tE, l->RE/AU, l->piE, l->tetE, l->Vt, l->u0, l->Ml,
                s->opt*1.0e6, l->Dl, s->Ds, l->vl, s->vs,
                s->mbs[0], s->fb[0],
                gg, static_cast<int>(l->struc),
                s->FWHM/year, vsave/s->mus, l->DeltaT/s->errA, l->murel*year,
                co->resu[0], co->resu[1], co->resu[2], co->resu[3], co->resu[5],
                co->resu[9], co->resu[10], co->resu[13], co->resu[14],
                s->Map[2], s->nsbl[2],
                co->flagi, s->Ai[2],
                ndw_L, ndw_R,
                detL, detR, detJ,
                co->okA[SJOINT], co->okA[SRUBIN], co->okA[SROMAN],
                co->Era[SJOINT][1], co->Era[SRUBIN][1], co->Era[SROMAN][1],   //sigma(tE)
                co->Era[SJOINT][3], co->Era[SRUBIN][3], co->Era[SROMAN][3],   //sigma(piE)
                co->Erb[SJOINT][0], co->Erb[SRUBIN][0], co->Erb[SROMAN][0],   //sigma(tetE)
                dclsEvent,
                synergyClass(*co),
                co->condA[SJOINT], co->condA[SRUBIN], co->condA[SROMAN],
                // ---- the rest of the row (Step D1) ----
                l->t0, s->xi, s->lon, s->lat,
                s->mbs[1], s->fb[1],
                {s->magb[0], s->magb[1], s->magb[2], s->magb[3], s->magb[4], s->magb[5], s->magb[6]},
                {s->blend[0], s->blend[1], s->blend[2], s->blend[3], s->blend[4], s->blend[5], s->blend[6]},
                co->relMl[SJOINT], co->relMl[SRUBIN], co->relMl[SROMAN],
                co->okB[SJOINT], co->okB[SRUBIN], co->okB[SROMAN],
                co->condB[SJOINT], co->condB[SRUBIN], co->condB[SROMAN],
                sched.dtToSeasonEdge(t0obs), sched.zone(t0obs),
                nres5_L, nres20_L, nresPSF_L, nres5_R, nres20_R, nresPSF_R,
                dsepMax_L, dsepMax_R
            });
   
            // ------------------------------------------------------------------------------
            // Step H2. Two quantities that make the satellite-parallax effect visible in the
            // table instead of only implicit in the Fisher matrix.
            //
            // du_sat  the observer separation in Einstein radii at t0, piE * D_perp / AU.
            //         This is the amplitude of the satellite effect for this event and the
            //         natural x-axis of every Step H3 figure. It is NOT L2_OFFSET_AU * piE:
            //         D_perp is the separation projected perpendicular to the line of sight
            //         and runs ~0.87-0.99 of the full L2 offset around the year, so the
            //         simple product is a ceiling (DEVIATIONS.md 28.2). Computed by asking
            //         lightcurve() for both observers rather than re-deriving the projection
            //         here, so the two can never drift apart.
            //
            // nepL_pk, nepR_pk  epochs from each survey within +-2 tE of t0, i.e. while the
            //         event is actually magnified. Satellite parallax needs CONTEMPORANEOUS
            //         coverage: an event Roman saw in season 3 and Rubin saw in season 7 has
            //         none of it, however large ndw_L and ndw_R are. The plan asked for a
            //         single `nep_both` flag; two counts are the same cost and strictly more
            //         informative, and the flag is just (nepL_pk > 0 and nepR_pk > 0).
            //
            // Safe to call lightcurve() here: FisherM has already run, and the only state it
            // touches (s->ux/uy, s->def*, s->pos*, as->ue_n*) is not read by the row written
            // below and is recomputed from scratch by the next draw.
            // ------------------------------------------------------------------------------
            double duSat = 0.0;
            {
                lightcurve(*s, *l, *as, l->t0, 0);
                const double r1 = as->ue_n1, r2 = as->ue_n2;
                lightcurve(*s, *l, *as, l->t0, 1);
                const double dn1 = as->ue_n1 - r1, dn2 = as->ue_n2 - r2;
                duSat = l->piE * std::sqrt(dn1 * dn1 + dn2 * dn2);
            }
            int nepLpk = 0, nepRpk = 0;
            {
                const double win = 2.0 * l->tE;
                for (int i = 0; i < ndw; ++i) {
                    if (std::fabs(l->timn[i] - t0obs) > win) continue;   // Deviation 76
                    if (int(l->tele[i]) == 1) nepRpk += 1;
                    else                      nepLpk += 1;
                }
            }

            filg_in.open(testf, std::ios::app);
            filg_in << icon           << " " << FFG[0]         << " " << l->tE                      << " "
                    << l->RE / AU     << " " << l->piE         << " " << l->tetE                    << " "
                    << l->Vt          << " " << l->u0          << " " << l->Ml                      << " " 
                    << s->opt * 1.0e6 << " " << l->Dl          << " " << s->Ds                      << " "
                    << l->vl          << " " << s->vs          << " " << s->mbs[0]                  << " "
                    << s->fb[0]       << " " << gg             << " " << static_cast<int>(l->struc) << " "
                    << s->FWHM / year << " " << vsave / s->mus << " " << l->DeltaT / s->errA        << " " << l->murel * year << " "
                    << co->resu[0]    << " " << co->resu[1]    << " " << co->resu[2]                << " " 
                    << co->resu[3]    << " " << co->resu[5]    << " "
                    << co->resu[9]    << " " << co->resu[10]   << " " << co->resu[13]               << " " << co->resu[14]    << " "
                    << s->Map[2]      << " " << s->nsbl[2]     << " " << co->flagi                  << " " << s->Ai[2]        << " "
                    // per-survey bookkeeping (Step C5), appended so existing column indices hold
                    << ndw_L << " " << ndw_R << " "
                    << detL  << " " << detR  << " " << detJ << " "
                    << co->okA[SJOINT] << " " << co->okA[SRUBIN] << " " << co->okA[SROMAN] << " "
                    << co->Era[SJOINT][1] << " " << co->Era[SRUBIN][1] << " " << co->Era[SROMAN][1] << " "
                    << co->Era[SJOINT][3] << " " << co->Era[SRUBIN][3] << " " << co->Era[SROMAN][3] << " "
                    << co->Erb[SJOINT][0] << " " << co->Erb[SRUBIN][0] << " " << co->Erb[SROMAN][0] << " "
                    << dclsEvent << " " << synergyClass(*co) << " "
                    << co->condA[SJOINT] << " " << co->condA[SRUBIN] << " " << co->condA[SROMAN] << " "
                    // the rest of the row (Step D1) -- appended, so columns 1-57 keep their indices
                    << l->t0 << " " << s->xi << " " << s->lon << " " << s->lat << " "
                    << s->mbs[1] << " " << s->fb[1] << " ";
            for (int i = 0; i < M; ++i) filg_in << s->magb[i]  << " ";
            for (int i = 0; i < M; ++i) filg_in << s->blend[i] << " ";
            filg_in << co->relMl[SJOINT] << " " << co->relMl[SRUBIN] << " " << co->relMl[SROMAN] << " "
                    << co->okB[SJOINT]   << " " << co->okB[SRUBIN]   << " " << co->okB[SROMAN]   << " "
                    << co->condB[SJOINT] << " " << co->condB[SRUBIN] << " " << co->condB[SROMAN] << " "
                    // Gap geometry. dt_edge is NEGATIVE when t0 fell inside a Roman season;
                    // t0zone distinguishes a mid-mission gap (1) from before-launch/after-end
                    // (2), which must never be pooled -- only the former is gap-filling.
                    << sched.dtToSeasonEdge(t0obs) << " " << sched.zone(t0obs) << " "
                    // Sky area this event's sightline stands for, deg^2 (Step E1). Constant
                    // across an unstratified run; NOT constant once --stride-roman is used,
                    // and then any statistic pooled over sightlines must weight by it.
                    << wArea << " "
                    // Step H2: the satellite-parallax observable and contemporaneous coverage.
                    << duSat << " " << nepLpk << " " << nepRpk << " "
                    // Step R1: resolving the two images. Counts of qualifying epochs per
                    // survey, then the largest separation reached while both were detectable
                    // (-1 = never). The three bars differ only in what counts as "resolved".
                    << nres5_L << " " << nres20_L << " " << nresPSF_L << " " << dsepMax_L << " "
                    << nres5_R << " " << nres20_R << " " << nresPSF_R << " " << dsepMax_R << " "
                    // Deviation 71: astrometric noise variants N and P (joint, Roman).
                    << co->ErbV[AV_N][SJOINT][0] << " " << co->ErbV[AV_N][SROMAN][0] << " "
                    << co->ErbV[AV_P][SJOINT][0] << " " << co->ErbV[AV_P][SROMAN][0] << " "
                    << co->relMlV[AV_N][SJOINT] << " " << co->relMlV[AV_N][SROMAN] << " "
                    << co->relMlV[AV_P][SJOINT] << " " << co->relMlV[AV_P][SROMAN] << " "
                    << co->okBV[AV_N][SJOINT] << " " << co->okBV[AV_N][SROMAN] << " "
                    << co->okBV[AV_P][SJOINT] << " " << co->okBV[AV_P][SROMAN] << " "
                    << int(l->luminous) << " " << s->fLens[0] << " " << s->fLens[1] << " "
                    << t0obs << " " << uminObs << "\n";
            filg_in.close();

            // ------------------------------------------------------------------------------
            // Step S1. Commit this draw's buffered light curve if it fills a requested
            // sample class. Here, and not in the time loop, because detL/detR/detJ and the
            // Fisher sigmas -- which is what the classes are defined in terms of -- do not
            // exist until now.
            //
            // First match wins: an event is written once, under the first class in the spec
            // file whose quota is still open. Otherwise a `both` event would also land in
            // `any` and the same light curve would be drawn twice in one figure.
            // ------------------------------------------------------------------------------
            if (dumpSpec.on and not dumpBuf.empty()) {
                double maxShift = 0.0;
                for (const DumpEpoch& e : dumpBuf)
                    maxShift = std::max(maxShift,
                                        std::sqrt(e.def1c * e.def1c + e.def2c * e.def2c));

                const SampleFacts facts{detL, detR, detJ, l->tE, sched.zone(t0obs),
                                        nepLpk, nepRpk, ndw_L, ndw_R,
                                        co->okB[SROMAN], maxShift};

                for (SampleClass& c : dumpSpec.classes) {
                    if (c.kept >= c.quota)        continue;
                    if (!sampleMatch(c, facts))   continue;

                    std::ostringstream pr;
                    pr << std::setprecision(10)
                       << "# Sample event for class '" << c.name << "'.\n"
                       << "# One key per line. Angles in deg, times in days, masses in Msun,\n"
                       << "# distances in kpc, angular scales in mas. A sigma of -1 means the\n"
                       << "# parameter was NOT measured by that survey (inactive or no epochs),\n"
                       << "# never that it was measured to be -1.\n"
                       << "class "        << c.name        << "\n"
                       << "population "   << gPop->name    << "\n"
                       << "lon "          << s->lon        << "\n"
                       << "lat "          << s->lat        << "\n"
                       << "tE "           << l->tE         << "\n"
                       << "t0 "           << l->t0         << "\n"
                       << "u0 "           << l->u0         << "\n"
                       << "xi "           << s->xi         << "\n"
                       << "piE "          << l->piE        << "\n"
                       << "tetE "         << l->tetE       << "\n"
                       << "Ml "           << l->Ml         << "\n"
                       << "Dl "           << l->Dl         << "\n"
                       << "Ds "           << s->Ds         << "\n"
                       << "Vt "           << l->Vt         << "\n"
                       << "murel_yr "     << l->murel * year << "\n"
                       << "mus1 "         << s->mus1       << "\n"
                       << "mus2 "         << s->mus2       << "\n"
                       << "mul1 "         << l->mul1       << "\n"
                       << "mul2 "         << l->mul2       << "\n"
                       << "lens_struc "   << int(l->struc) << "\n"
                       << "mbs0 "         << s->mbs[0]     << "\n"
                       << "fb0 "          << s->fb[0]      << "\n"
                       << "mbs1 "         << s->mbs[1]     << "\n"
                       << "fb1 "          << s->fb[1]      << "\n";
                    pr << "magb";  for (int i = 0; i < M; ++i) pr << " " << s->magb[i];
                    pr << "\nblend"; for (int i = 0; i < M; ++i) pr << " " << s->blend[i];
                    pr << "\n"
                       << "ndw_L "        << ndw_L         << "\n"
                       << "ndw_R "        << ndw_R         << "\n"
                       << "nep_pk_L "     << nepLpk        << "\n"
                       << "nep_pk_R "     << nepRpk        << "\n"
                       << "detL "         << detL          << "\n"
                       << "detR "         << detR          << "\n"
                       << "detJ "         << detJ          << "\n"
                       << "dt_edge "      << sched.dtToSeasonEdge(t0obs) << "\n"
                       << "t0zone "       << sched.zone(t0obs)           << "\n"
                       << "du_sat "       << duSat         << "\n"
                       << "max_shift "    << maxShift      << "\n"
                       << "dsep_max_L "   << dsepMax_L     << "\n"
                       << "dsep_max_R "   << dsepMax_R     << "\n"
                       << "okA_J "  << co->okA[SJOINT] << " okA_L " << co->okA[SRUBIN]
                       << " okA_R " << co->okA[SROMAN] << "\n"
                       << "okB_J "  << co->okB[SJOINT] << " okB_L " << co->okB[SRUBIN]
                       << " okB_R " << co->okB[SROMAN] << "\n"
                       << "sigtE_J "   << co->Era[SJOINT][1] << " sigtE_L " << co->Era[SRUBIN][1]
                       << " sigtE_R "  << co->Era[SROMAN][1] << "\n"
                       << "sigpiE_J "  << co->Era[SJOINT][3] << " sigpiE_L " << co->Era[SRUBIN][3]
                       << " sigpiE_R " << co->Era[SROMAN][3] << "\n"
                       << "sigtetE_J "  << co->Erb[SJOINT][0] << " sigtetE_L " << co->Erb[SRUBIN][0]
                       << " sigtetE_R " << co->Erb[SROMAN][0] << "\n"
                       << "relMl_J "  << co->relMl[SJOINT] << " relMl_L " << co->relMl[SRUBIN]
                       << " relMl_R " << co->relMl[SROMAN] << "\n";

                    std::ostringstream idss;
                    idss << std::setw(3) << std::setfill('0') << (++dumpSeq);
                    writeSampleEvent(dumpSpec, c.name, idss.str(), dumpBuf, pr.str(),
                                     *s, *l, *as);
                    c.kept += 1;
                    std::cout << "  [sample] " << c.name << " " << idss.str()
                              << "  tE=" << l->tE << "d  Ml=" << l->Ml << "Msun"
                              << "  ndw_L=" << ndw_L << " ndw_R=" << ndw_R
                              << "  (" << c.kept << "/" << c.quota << ")" << std::endl;
                    break;
                }
            }

            // Step H3. One row per DETECTED event, carrying both forecasts for that same
            // event. Written here rather than beside the Fisher call because duSat and the
            // contemporaneous-coverage counts are computed above, and they are the axes
            // every H3 figure uses.
            if (cfg.pairSat and (detL or detR or detJ)) {
                std::ofstream fpair(fnPair, std::ios::app);
                const int okAs = co->okA[SJOINT], okBs = co->okB[SJOINT];
                fpair << std::setprecision(7)
                      << s->lon << " " << s->lat << " "
                      << l->tE  << " " << l->u0  << " " << l->piE << " " << l->tetE << " "
                      << duSat << " "
                      << okAs << " " << okNS << " " << okBs << " " << okNSb << " "
                      << (okAs ? co->Era[SJOINT][1] : -1.0) << " " << sigtE_ns   << " "
                      << (okAs ? co->Era[SJOINT][3] : -1.0) << " " << sigpiE_ns  << " "
                      << (co->okA[SROMAN] ? co->Era[SROMAN][3] : -1.0) << " " << sigpiER_ns << " "
                      << (okBs ? co->Erb[SJOINT][0] : -1.0) << " " << sigtetE_ns << " "
                      << (okBs ? co->Erb[SJOINT][3] : -1.0) << " " << sigpiEb_ns << " "
                      << co->relMl[SJOINT] << " " << relMl_ns << " "
                      << co->condA[SJOINT] << " " << condA_ns << " "
                      << co->condB[SJOINT] << " " << condB_ns << " "
                      << nepLpk << " " << nepRpk << " " << wArea << " "
                      << l->Ml << " " << l->Dl << " " << s->Ds << " " << l->Vt << "\n";
                fpair.close();
            }
//          


            //cout << "** End of saving in the file *********" << save << endl;
            //cout << "icon: " << icon << "\tnlens: " << nlens << "\tnerr: " << nerr << endl;
            // Per-sightline event budget -- see RunConfig. Was hardcoded to a stub
            // 20/5/1.0 with the production 850/150/2.0 commented out beside it.
            } while ((icon < cfg.iconTarget or nlens < cfg.nlensTarget or nerr < cfg.nerrTarget)
                     and nsim < cfg.maxDraws);

            // Did the sightline actually meet its budget, or did the cap stop it? The
            // distinction matters: a capped sightline's Poisson precision is whatever it
            // reached, not what was asked for, and averaging it in as though it were a
            // full sample would understate the error bars.
            const bool budgetMet = (icon  >= cfg.iconTarget and
                                    nlens >= cfg.nlensTarget and
                                    nerr  >= cfg.nerrTarget);
            if (!budgetMet) {
                nCapped += 1;
                cout << "  CAP: stopped at nsim = " << nsim << " with icon = " << icon
                     << "/" << cfg.iconTarget << ", nlens = " << nlens << "/" << cfg.nlensTarget
                     << ", nerr = " << nerr << "/" << cfg.nerrTarget << endl;
            }

            // Some epochs existed but nothing survived to be aggregated. Same reasoning as
            // the no-coverage skip above: the CHECK block below requires at least one
            // detected AND one characterised event, so this must not fall through.
            if (nlens < 1 or nerr <= 0.0) {
                nSkipBarren += 1;
                areaBarren  += wArea;
                cout << "  BARREN: " << nlens << " detected, " << nerr
                     << " characterised -- nothing to aggregate at this sightline" << endl;
                continue;
            }
            nAggregated    += 1;
            areaAggregated += wArea;
   
            for (int i = 0; i <= GG; ++i) {
                l->NstE[i] += l->nstE[i];
                l->NdtE[i] += l->ndtE[i];
                l->ndtE[i]  = double(l->ndtE[i] / (l->nstE[i] + eps)); // [0.0, 1.0]

                fil2 << std::fixed << std::setprecision(4)
                     << l->tEs[i]  << " "
                     << double(l->NdtE[i]  * 100.0 / (l->NstE[i] + eps)) << " "
                     << l->Mls[i]  << " "
                     << double(l->NdMl[i]  * 100.0 / (l->NsMl[i] + eps)) << " "
                     << l->pis[i]  << " "
                     << double(l->Ndpi[i]  * 100.0 / (l->Nspi[i] + eps)) << " "
                     << l->u0s[i]  << " "
                     << double(l->Ndu0[i]  * 100.0 / (l->Nsu0[i] + eps)) << " "
                     << l->mbs[i]  << " "
                     << double(l->Ndmb[i]  * 100.0 / (l->Nsmb[i] + eps)) << " "
                     << l->fbs[i]  << " "
                     << double(l->Ndfb[i]  * 100.0 / (l->Nsfb[i] + eps)) << " "
                     << l->mus[i]  << " "
                     << double(l->Ndmu[i]  * 100.0 / (l->Nsmu[i] + eps)) << " "
                     // The + eps every other column has, and these two lacked: before the
                     // counters above existed this was 0/0 and printed "-nan" on every row.
                     << double(l->Nhalo[1] * 100.0 / (l->Nhalo[0] + eps)) << " "
                     << double(l->Nself[1] * 100.0 / (l->Nself[0] + eps)) << "\n";

                fil2b << std::fixed  << std::setprecision(1)
                      << l->NdtE[i]  << " " << l->NstE[i]  << " "
                      << l->NdMl[i]  << " " << l->NsMl[i]  << " "
                      << l->Ndpi[i]  << " " << l->Nspi[i]  << " "
                      << l->Ndu0[i]  << " " << l->Nsu0[i]  << " "
                      << l->Ndmb[i]  << " " << l->Nsmb[i]  << " "
                      << l->Ndfb[i]  << " " << l->Nsfb[i]  << " "
                      << l->Ndmu[i]  << " " << l->Nsmu[i]  << " "
                      << l->Nhalo[1] << " " << l->Nhalo[0] << " "
                      << l->Nself[1] << " " << l->Nself[0] << "\n";
            }
   
            s->nstart = 0.0;
            for (int i = 0; i < Num; ++i) {
                test = double(s->nsdet[i]   / (s->nssim[i] + eps));
                s->nstart  += s->Rostari[i] * (s->Nstart   / s->Rostart) * test;
            } //Number/deg^{2}

            for (int i = 0; i < 2; ++i){
                tE[i]   = 0.0, RE[i]  = 0.0, piE[i]  = 0.0, tetE[i]  = 0.0; Vt[i]    = 0.0;
                u0[i]   = 0.0, Ml[i]  = 0.0, opd[i]  = 0.0; Dl[i]    = 0.0, Ds[i]    = 0.0;
                vl[i]   = 0.0, vs[i]  = 0.0, mbs[i]  = 0.0, fb[i]    = 0.0; numd[i]  = 0.0;
                fwhm[i] = 0.0; vsn[i] = 0.0; DelT[i] = 0.0; Struc[i] = 0.0; murel[i] = 0.0;
                Map[i]  = 0.0; nbl[i] = 0.0; Ext[i]  = 0.0;
            }

            EffiL = 0.0; EFF   = 0.0; Gamma  = 0.0; Neven = 0.0; EffiD = 0.0;
            ErtE  = 0.0; ErpiE = 0.0; ErtetE = 0.0; Erml  = 0.0; Erfb  = 0.0; nErAvg = 0.0;
            Erdl  = 0.0; Ermul = 0.0; Ermus  = 0.0; Eru0  = 0.0;
  
    int tempStruc;
    for (const auto& r : records) {
//        counter = r.counter;
        flagL = r.flagL;
        l->tE = r.tE;   l->RE = r.RE;   l->piE = r.piE;  l->tetE = r.tetE;
        l->Vt = r.Vt;   l->u0 = r.u0;   l->Ml  = r.Ml;
        s->opt = r.opt; l->Dl = r.Dl;   s->Ds  = r.Ds;   l->vl = r.vl; s->vs = r.vs;
        mbase = r.mbase; fblend = r.fblend; gg = r.gg; tempStruc = r.struc;
        s->FWHM = r.FWHM; vsave = r.vsave; l->DeltaT = r.DeltaT; l->murel = r.murel;
        co->resu[0]=r.resu0;  co->resu[1]=r.resu1;   co->resu[2]=r.resu2;
        co->resu[3]=r.resu3;  co->resu[5]=r.resu5;   co->resu[9]=r.resu9;
        co->resu[10]=r.resu10; co->resu[13]=r.resu13; co->resu[14]=r.resu14;
        s->Map[2]=r.Map2; s->nsbl[2]=r.nsbl2; co->flagi=r.flagi; s->Ai[2]=r.Ai2;
        // Must be replayed too: this loop runs after the whole field, so co->okA
        // otherwise holds whatever the LAST FisherM call left, not this event's.
        co->okA[SJOINT] = r.okJoint;

        l->struc = static_cast<GalacticComponent>(tempStruc);

        // ------------------ Update cumulative sums ------------------
        Struc[0] += 1.0;
        tE[0]    += l->tE;    RE[0]  += l->RE;     piE[0]  += l->piE;
        tetE[0]  += l->tetE;  Vt[0]  += l->Vt;     u0[0]   += l->u0;
        Ml[0]    += l->Ml;    opd[0] += s->opt;    Dl[0]   += l->Dl;
        Ds[0]    += s->Ds;    vl[0]  += l->vl;     vs[0]   += s->vs;
        mbs[0]   += mbase;    fb[0]  += fblend;    numd[0] += 1.0;
        fwhm[0]  += s->FWHM;  vsn[0] += vsave;     DelT[0] += l->DeltaT;
        murel[0] += l->murel; Map[0] += s->Map[2]; nbl[0]  += s->nsbl[2];
        Ext[0]   += s->Ai[2];

        // ------------------ Flagged case ---------------------------
        if (flagL > 0) {
            Struc[1] += 1.0;
            tE[1]    += l->tE;    RE[1]  += l->RE;     piE[1]  += l->piE;
            tetE[1]  += l->tetE;  Vt[1]  += l->Vt;     u0[1]   += l->u0;
            Ml[1]    += l->Ml;    opd[1] += s->opt;    Dl[1]   += l->Dl;
            Ds[1]    += s->Ds;    vl[1]  += l->vl;     vs[1]   += s->vs;
            mbs[1]   += mbase;    fb[1]  += fblend;    numd[1] += 1.0;
            fwhm[1]  += s->FWHM;  vsn[1] += vsave;     DelT[1] += l->DeltaT;
            murel[1] += l->murel; Map[1] += s->Map[2]; nbl[1]  += s->nsbl[2];
            Ext[1]   += s->Ai[2];

            EFF += static_cast<double>(l->ndtE[gg] / (l->tE / year)); // 1/years

            // Average only over events the joint fit could actually characterize. Step C4
            // reports sigma = -1 as an explicit "not characterizable" sentinel, and ErrorCal
            // divides it by the parameter value, so an unguarded sum pulls in large negative
            // numbers: one event in this field contributed resu[3] = -697, dragging the mean
            // fractional piE error negative and tripping CHECK(ErpiE > 0.0).
            //
            // Skipping those events is not the selection bias DEVIATIONS entry 8 warns about.
            // That rule is about never dropping an event from the joint-vs-single RATIO
            // statistics, where a missing single-survey sigma is itself the result. Here we are
            // forming a mean precision, and an event with no measurement has no precision to
            // average -- including it would be averaging a sentinel. The count of events the
            // mean is actually over is tracked separately so the denominator is honest.
            //
            // okA[SJOINT] is NOT sufficient on its own. It says the joint photometric
            // matrix inverted -- not that every parameter was in the fit. Since the joint
            // refactor gave each survey partition its own active parameter subset
            // (activePhotParams in Bulge.h), fb0 and mbs0 enter the joint fit only when the
            // event has Rubin epochs, and fb1/mbs1 only when it has Roman ones. An event
            // detected by Roman with no Rubin data therefore has a perfectly valid joint
            // fit in which Era[2] is still the -1.0 sentinel, and ErrorCal divides that by
            // fb0 regardless: one such event contributed resu[2] = -8499 at l=0.881,
            // b=-0.94 and dragged a 102-event mean to -83, tripping CHECK(Erfb > 0.0).
            //
            // So test the values themselves. An event is averaged only if all nine are real
            // measurements, which keeps every mean over the same event set and keeps
            // nErAvg meaningful as a single denominator. The cost is small and measured:
            // across 48,959 characterised events in the 2026-08-29 run exactly ONE was
            // excluded by this, and no column other than resu[2] was ever negative.
            const bool allMeasured = (co->flagi > 0 and co->okA[SJOINT]
                                      and co->resu[0]  >= 0.0 and co->resu[1]  >= 0.0
                                      and co->resu[2]  >= 0.0 and co->resu[3]  >= 0.0
                                      and co->resu[5]  >= 0.0 and co->resu[9]  >= 0.0
                                      and co->resu[10] >= 0.0 and co->resu[13] >= 0.0
                                      and co->resu[14] >= 0.0);
            if (allMeasured) {
                Eru0  += co->resu[0];  ErtE   += co->resu[1];  Erfb  += co->resu[2];
                ErpiE += co->resu[3];  ErtetE += co->resu[5];  Erml  += co->resu[9];
                Erdl  += co->resu[10]; Ermul  += co->resu[13]; Ermus += co->resu[14];
                nErAvg += 1.0;
            }
        }
    }

    for (int i = 0; i < 2; ++i) { // What is the purpose of this block?
        tE[i]    = double(tE[i]            / (numd[i] + eps)) / year; //[years]  
        RE[i]    = double(RE[i]            / (numd[i] + eps));//[AU]  
        piE[i]   = double(piE[i]           / (numd[i] + eps));//[]     
        tetE[i]  = double(tetE[i]          / (numd[i] + eps));//[mas]  
        Vt[i]    = double(Vt[i]            / (numd[i] + eps));//[km/s]  
        u0[i]    = double(u0[i]            / (numd[i] + eps));//[]  
        Ml[i]    = double(Ml[i]            / (numd[i] + eps));//[Msun]  
        opd[i]   = double(opd[i]           / (numd[i] + eps)) * u0m * u0m;//[] x10^{6}
        Dl[i]    = double(Dl[i]            / (numd[i] + eps));//[kpc]  
        Ds[i]    = double(Ds[i]            / (numd[i] + eps));//[kpc]  
        vl[i]    = double(vl[i]            / (numd[i] + eps));//[km/s]  
        vs[i]    = double(vs[i]            / (numd[i] + eps));//[km/s]  
        mbs[i]   = double(mbs[i]           / (numd[i] + eps));//[mag]  
        fb[i]    = double(fb[i]            / (numd[i] + eps));//[]
        fwhm[i]  = double(fwhm[i]          / (numd[i] + eps));//[years]
        vsn[i]   = double(vsn[i]           / (numd[i] + eps));//[]
        DelT[i]  = double(DelT[i]          / (numd[i] + eps));//[]
        murel[i] = double(murel[i]         / (numd[i] + eps));//[mas/year]
        Map[i]   = double(Map[i]           / (numd[i] + eps));//[mag]
        nbl[i]   = double(nbl[i]           / (numd[i] + eps));//[]
        Ext[i]   = double(Ext[i]           / (numd[i] + eps));//[mag]
        Struc[i] = double(Struc[i] * 100.0 / (numd[i] + eps));
    }

    EFF = double(EFF / (numd[1] + eps));//1/[years]
    Gamma = double(2.0 / M_PI * opd[0] * 1.0e-6 * EFF) / u0m; // 1/[star*year]  
    EffiL = double(numd[1] * 100.0 / (numd[0] + eps)); // probability of detecting lensing  
    EffiD = double(icon * 100.0 / (nsim    + eps)); // % of drawn stars that are visible (Deviation 78; was numd[0], i.e. 100% by construction)
    Neven = double(s->nstart * Gamma * 10.0);//deg^{-2}
   
    Eru0   = double(Eru0   / (nErAvg + eps));  
    Erfb   = double(Erfb   / (nErAvg + eps));  
    ErtE   = double(ErtE   / (nErAvg + eps));  
    ErpiE  = double(ErpiE  / (nErAvg + eps));  
    ErtetE = double(ErtetE / (nErAvg + eps));  
    Erml   = double(Erml   / (nErAvg + eps));  
    Erdl   = double(Erdl   / (nErAvg + eps));  
    Ermul  = double(Ermul  / (nErAvg + eps));  
    Ermus  = double(Ermus  / (nErAvg + eps));  
   
    //Maps
    fil3 << std::fixed << std::setprecision(6)
         << tE[0]    << " " << tE[1]    << " "
         << RE[0]    << " " << RE[1]    << " "
         << piE[0]   << " " << piE[1]   << " "
         << tetE[0]  << " " << tetE[1]  << " "
         << Vt[0]    << " " << Vt[1]    << " "
         << u0[0]    << " " << u0[1]    << " "
         << Ml[0]    << " " << Ml[1]    << " "
         << opd[0]   << " " << opd[1]   << " "
         << Dl[0]    << " " << Dl[1]    << " "
         << Ds[0]    << " " << Ds[1]    << " "
         << vl[0]    << " " << vl[1]    << " "
         << vs[0]    << " " << vs[1]    << " "
         << mbs[0]   << " " << mbs[1]   << " "
         << fb[0]    << " " << fb[1]    << " "
         << fwhm[0]  << " " << fwhm[1]  << " "
         << vsn[0]   << " " << vsn[1]   << " "
         << DelT[0]  << " " << DelT[1]  << " "
         << Struc[0] << " " << Struc[1] << " "
         << murel[0] << " " << murel[1] << " "
         << Map[0]   << " " << Map[1]   << " "
         << nbl[0]   << " " << nbl[1]   << " "
         << Ext[0]   << " " << Ext[1]   << " ";

    fil3 << std::setprecision(8)
         << EffiD      << " " << EffiL        << " "
         << std::log10(EFF) << " " << std::log10(Gamma) << " " << std::log10(Neven) << " "
         << Eru0       << " " << ErtE         << " " << Erfb         << " " << ErpiE << " " << ErtetE << " "
         << Erml       << " " << Erdl         << " " << Ermul        << " " << Ermus << " ";

    fil3 << std::setprecision(1)
         << nsim              << " " << numd[0]          << " " << numd[1] << " "
         << nerr              << " " << nri              << " " << nde     << " "
         << std::log10(s->Rostart) << " " << std::log10(s->Nstart) << " " << std::log10(s->nstart)
         // Step E1. Three columns appended, in this order:
         //   w_area  deg^2 of sky this sightline stands for -- no longer a run-wide constant
         //   lon,lat where it is. The map file had NO position column at all, so a row in it
         //           could not be tied to the events it produced, and the draw count `nsim`
         //           it records -- the denominator any pooled yield needs -- was unreachable
         //           from the event table. With these, an event joins its sightline on
         //           (lon, lat) and the correct pooled weight, w_area/nsim, is computable.
         << " " << std::setprecision(8) << wArea
         << " " << std::setprecision(6) << s->lon << " " << s->lat
         << "\n";

    // Flush the per-sightline outputs now rather than when the stream is destroyed. Every
    // production pause so far has been a kill, and a kill discards whatever is still
    // buffered: the 2026-09-06 chunk-1 stop lost the map rows of six completed sightlines
    // and left a half-written seventh, onto which the resuming run's first row was then
    // appended -- one 122-field line that made the whole file unreadable until the Python
    // reader learned to skip it (Deviation 41). The same stop cost EfLMC5/EfLMC5B the same
    // six blocks, which is why fil2/fil2b are flushed here too. A sightline costs minutes
    // of CPU, so three flushes per sightline are free, and what reaches disk is then what
    // the log says was finished.
    fil3.flush();
    fil2.flush();
    fil2b.flush();

    cout << "nsim:  "  << nsim    << "\t Ndetected:  " << icon    << "\t Nlensing:  " << nlens << "\t NError:  " << nerr << endl;
    cout << "Detection classes:";
    for (int c = 0; c < NDETCLASS; ++c)
        cout << "  " << detClassName(c) << ": " << nDetClass[c];
    cout << endl;
    if (nDetClass[DET_ANOMALY] > 0) {
        cout << "  WARNING: " << nDetClass[DET_ANOMALY] << " event(s) detected by one telescope "
             << "but NOT by the joint test. Adding data cannot destroy signal, so this is a "
             << "threshold inconsistency -- see DetClass in include/fisher/fisher.h." << endl;
    }
    cout << "numd0:  " << numd[0] << "\t numd1:  "     << numd[1] << endl;
    cout << "EFF:  "   << EFF     << "\t Gamma:  "     << Gamma   << "\t Neven:  "    << Neven << endl;
    cout << "l.tE:  "  << l->tE   << "\t l.Ml:  "      << l->Ml   << endl;
   
    CHECK(EFF > 0.0);
    CHECK(Gamma > 0.0);
    CHECK(EffiL > 0.0);
    CHECK(EffiD > 0.0);
    CHECK(Neven > 0.0);
    
    CHECK(DelT[0] > 0.0);
    CHECK(DelT[1] > 0.0);
    CHECK(vsn[0] > 0.0);
    CHECK(vsn[1] > 0.0);
    CHECK(fwhm[0] > 0.0);
    CHECK(fwhm[1] > 0.0);
    
    CHECK(fb[0] > 0.0);
    CHECK(fb[0] <= 1.0);
    CHECK(fb[1] > 0.0);
    CHECK(fb[1] <= 1.0);
    
    CHECK(mbs[0] > 0.0);
    CHECK(mbs[1] > 0.0);
    
    CHECK(vs[0] > 0.0);
    CHECK(vs[1] > 0.0);
    CHECK(vl[0] > 0.0);
    CHECK(vl[1] > 0.0);
    
    CHECK(Ds[0] > 0.0);
    CHECK(Ds[1] > 0.0);
    CHECK(Dl[0] > 0.0);
    CHECK(Dl[1] >= 0.0);
    
    CHECK(piE[1] >= 0.0);
    CHECK(tetE[1] > 0.0);
    
    CHECK(murel[0] > 0.0);
    CHECK(murel[1] > 0.0);
    
    if (nErAvg > 0.0) CHECK(ErtE > 0.0);
    if (nErAvg > 0.0) CHECK(ErpiE > 0.0);
    if (nErAvg > 0.0) CHECK(ErtetE > 0.0);
    if (nErAvg > 0.0) CHECK(Erml > 0.0);
    if (nErAvg > 0.0) CHECK(Erdl > 0.0);
    if (nErAvg > 0.0) CHECK(Ermul > 0.0);
    if (nErAvg > 0.0) CHECK(Ermus > 0.0);
    if (nErAvg > 0.0) CHECK(Eru0 > 0.0);
    if (nErAvg > 0.0) CHECK(Erfb > 0.0);
    
    CHECK(nerr != 0.0);
    CHECK(numd[0] != 0.0);
    CHECK(numd[1] != 0.0);
    CHECK(nsim != 0.0);
    
    // NOT an equality. `icon` counts stars that were OBSERVABLE (flagf > 0 and ndw > 2);
    // `numd[0]` counts every record pushed, and the push site sits OUTSIDE that gate, so a
    // star that was drawn but never observable still gets a record. The two are equal only
    // where every draw is observable, which is true in the dense stub patch that every run
    // before commit 81a6b04 used and false as soon as the scan reaches sparse sky -- at
    // l=-3.499, b=-1.98 it is 5 observable out of 15 drawn.
    //
    // Loosening this assertion changes no computed value. It does expose a real question
    // about what the per-sightline denominators mean -- EffiD = numd[0]/nsim is 100% by
    // construction, and the [0] means average over drawn rather than observed stars. That
    // is a science decision, recorded in OPEN_ITEMS.md, not something to change here.
    CHECK(icon <= numd[0]);
    CHECK(numd[1] == nlens);
     
    cout << "==============================================================" << endl;  
  
   }}//end of the sightline loop (was: right_accention and declinaton)

    // ---- Stage 8: end-of-run report ----
    RunTotals totals;
    totals.nDchiMismatch   = nDchiMismatch;
    totals.nAggregated     = nAggregated;
    totals.nSkipNoCoverage = nSkipNoCoverage;
    totals.nSkipBarren     = nSkipBarren;
    totals.nCapped         = nCapped;
    totals.areaAggregated  = areaAggregated;
    totals.areaNoCoverage  = areaNoCoverage;
    totals.areaBarren      = areaBarren;
    totals.nSimTot         = nSimTot;
    totals.NDetClassTot    = NDetClassTot;
    totals.NDetClassTE     = NDetClassTE;
    writeRunSummary(cfg, totals, *l);

    return(0);
}
