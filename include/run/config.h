// The run configuration (command-line flags) and the usage text.
#ifndef ROMAN_RUN_CONFIG_H
#define ROMAN_RUN_CONFIG_H

#include "common.h"

// Run configuration. Every field is a cost or coverage lever, so each is a runtime flag
// rather than a constant: a smoke test and a production run differ only in these numbers.
struct RunConfig {
    // Sightline grid. The scan steps by `stride * dd` degrees: stride=1 is the native
    // 0.02 deg grid (166,397 sightlines, not runnable), stride=10 is 0.20 deg (1,706).
    int    stride      = 10;

    // Delta-chi2 a lensing model must beat a flat baseline by for an event to count as
    // detected (detL, detR, detJ). A fixed bar, not scaled by the epoch count; see
    // DCHI_DET_DEFAULT in config/parameters.h. Recorded in run_provenance.txt.
    double dchiDet     = DCHI_DET_DEFAULT;

    // Day-shared astrometric error of the N and P noise variants [mas]; --ast-sigc N,P overrides the
    // AST_SIGC defaults, e.g. for a scan of theta_E precision against it. Recorded in run_provenance.txt.
    double astSigcN    = AST_SIGC_N;
    double astSigcP    = AST_SIGC_P;

    // Sightline grid inside Roman's footprint, in the same units. Roman covers ~2.6% of the
    // scan region, so a uniform grid spends most of its time where the joint fit is just
    // Rubin's. Sightlines within a GBTDS field's detector outline are visited on a
    // `strideRoman * dd` grid; everything else stays on the coarse `stride * dd` grid, and
    // every sightline carries the sky area it stands for so survey-wide totals are recoverable.
    // 0 means "same as --stride", which reproduces the unstratified scan exactly.
    int    strideRoman = 0;

    // Per-sightline event budget: the do/while over stars stops once ALL three are met.
    // icon counts detected events, nlens those also Fisher-characterised, nerr an
    // accumulated Fisher-error weight.
    int    iconTarget  = DEFAULT_EVENTS_TARGET;
    int    nlensTarget = DEFAULT_LENSES_TARGET;
    double nerrTarget  = DEFAULT_NERR_TARGET;

    // Hard cap on stars drawn at one sightline, regardless of the budgets above. The targets
    // are combined with AND and nerr only advances when FisherM succeeds, so a sightline
    // where no event is ever characterisable would never exit. Sightlines with no coverage
    // at all are skipped before the loop; the cap handles the partial case. The run reports
    // how many sightlines hit it. The default is ~59x the draws a well-covered bulge
    // sightline needs to meet the full budget.
    double maxDraws    = DEFAULT_MAXDRAWS;

    // Restrict the scan to a hardcoded 0.1x0.1 deg patch instead of the full region.
    bool   stubPatch   = false;
    // Index into the (deterministic) sightline scan vector at which to begin (resume).
    long   startIndex  = 0;
    // One past the last sightline to simulate (-1 = to the end), and the base seed from which
    // every sightline's own seed is derived. With --start-index they let a scan be cut into
    // independent chunks that reproduce the whole run.
    long   endIndex    = -1;
    unsigned long long seedBase = seed;

    // Put Roman at the centre of the Earth, removing the Earth-L2 spatial baseline while
    // leaving the timing untouched. With this flag the output is identical to a run without
    // the satellite-parallax term.
    bool   noSatPar    = false;
    // Characterise every detected event twice, once with Roman at L2 and once with the offset
    // zeroed, and write both forecasts to a side file. Two separate runs cannot isolate the
    // effect, because moving the observer changes which events are detected.
    bool   pairSat     = false;

    // Path to a sample-event dump spec, or empty for no dump. The dump consumes no RNG, so a
    // run with it is event-for-event the same as one without.
    std::string dumpSpec;

    // Build the sightline grid, write the provenance block, report the strata, and stop before
    // drawing a star. Gives the sightline counts before committing to a long run: refining the
    // footprint stratum by k multiplies its sightlines by k^2.
    bool   dryRun      = false;
};

void printUsage(const char* prog);

// Parse and validate argv into `cfg` (the population choice lands in gPop). Returns true when the
// run should go on; false means main() must return `exitCode` -- 0 after --help, the code of
// the error message already printed otherwise.
bool parseCommandLine(int argc, char** argv, RunConfig& cfg, int& exitCode);

// The sightline-grid steps, derived from the flags and the GBTDS detector size.
struct GridSteps {
    double gridStep;   // coarse grid step [deg] = stride * dd
    int    kSub;       // fine cells per coarse cell, per axis
    double fineStep;   // footprint grid step [deg] = strideRoman * dd
};

struct GbtdsLayout;
// Resolve --stride-roman against the detector size, fill `steps`, and
// refuse a grid that steps over whole detectors. May set cfg.strideRoman. Returns 0, or the exit
// code of the error printed.
int resolveGridSteps(RunConfig& cfg, const GbtdsLayout& gl, GridSteps& steps);

#endif // ROMAN_RUN_CONFIG_H
