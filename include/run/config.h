// The run configuration (command-line flags) and the usage text.
#ifndef ROMAN_RUN_CONFIG_H
#define ROMAN_RUN_CONFIG_H

#include "common.h"

// ---------------------------------------------------------------------------
// Run configuration. Every field here is a COST or COVERAGE lever, which is why
// they are runtime flags rather than constants: a smoke test and a production
// run differ only in these numbers, and having to edit and recompile to switch
// between them is how a run ends up with no record of what produced it.
// ---------------------------------------------------------------------------
struct RunConfig {
    // Sightline grid. The scan steps by `stride * dd` degrees, so stride=1 is the
    // native 0.02 deg grid (166,397 sightlines -- not runnable) and stride=10 is
    // 0.20 deg (1,706 sightlines, ~15 h at the production budget).
    int    stride      = 10;

    // Step H7. The delta-chi2 a lensing model must beat a flat baseline by before the event
    // counts as detected, for each of detL, detR and detJ. A FIXED bar, not one scaled by the
    // epoch count: see the derivation at DCHI_DET_DEFAULT in config/parameters.h. Exposed because every
    // yield this project reports is conditioned on it, so it belongs in run_provenance.txt
    // and has to be variable for a sensitivity test.
    double dchiDet     = DCHI_DET_DEFAULT;

    // Sightline grid INSIDE Roman's footprint, in the same units (Step E1). Roman covers
    // ~2.6% of the scan region, so a uniform grid spends 97% of its wall clock on sky where
    // the joint fit is Rubin's matrix and nothing can be learned about the combination. This
    // is the second stride: sightlines within a GBTDS field's detector outline are visited on
    // a `strideRoman * dd` grid, everything else stays on the coarse `stride * dd` grid, and
    // every sightline carries the sky area it stands for so survey-wide totals are recoverable.
    //
    // 0 means "same as --stride", which reproduces the unstratified scan exactly -- same
    // sightline positions, same order, same RNG stream, same areas. Nothing changes unless
    // asked for.
    int    strideRoman = 0;

    // Per-sightline event budget: the do/while over stars stops once ALL three are
    // met. icon counts detected events, nlens those also Fisher-characterised,
    // nerr an accumulated Fisher-error weight. These set the Poisson precision of
    // every per-sightline quantity.
    int    iconTarget  = DEFAULT_EVENTS_TARGET;
    int    nlensTarget = DEFAULT_LENSES_TARGET;
    double nerrTarget  = DEFAULT_NERR_TARGET;

    // Hard cap on stars drawn at ONE sightline, regardless of the budgets above.
    //
    // The three targets are combined with AND: the loop runs until icon, nlens AND nerr
    // are all met. nerr only advances when FisherM succeeds, so a sightline where no
    // event is ever characterisable cannot satisfy it and the loop never exits. That is
    // not hypothetical -- the 2026-08-29 production attempt drew 331,931 events at
    // sightline 0 with ndw_L = ndw_R = 0 on every one of them and had to be killed.
    //
    // Sightlines with NO coverage at all are skipped outright before the loop starts, so
    // this cap is for the partial case: a few epochs exist, events are occasionally
    // detected, but Fisher almost never succeeds. Such a sightline would still run far
    // past any useful precision. The cap bounds it and the run reports how many sightlines
    // hit it, so a cap set too low announces itself rather than silently truncating.
    //
    // Default sized from measurement, not taste: a well-covered bulge sightline (l=0.5,
    // b=-1.0, 2412 Rubin and 50401 Roman epochs) meets the full 850/150/2.0 budget in
    // exactly 850 draws -- every draw there is observable. 5e4 is ~59x that, so the cap
    // cannot bite a sightline that is merely unlucky. Draw cost scales with the epoch
    // count, so the sightlines that could approach the cap are the sparse ones, where a
    // draw is ~1 ms (measured: 331,931 draws in 5.5 min at a zero-epoch sightline) and
    // 5e4 draws costs under a minute.
    double maxDraws    = DEFAULT_MAXDRAWS;

    // Restrict the scan to the old hardcoded 0.1x0.1 deg patch instead of the full
    // region. Kept only so a run can be compared against the pre-Step-4 numbers.
    bool   stubPatch   = false;
    // Step: resume. Index into the (deterministic) sightline scan vector at which to begin.
    // 0 means "start from the beginning", which is what every non-resumed run wants.
    long   startIndex  = 0;
    // Deviation 77: one past the last sightline to simulate (-1 = to the end), and the base seed
    // from which every sightline's own seed is derived. Together with --start-index they let a
    // scan be cut into independent chunks that reproduce the whole run.
    long   endIndex    = -1;
    unsigned long long seedBase = seed;

    // Put Roman back at the centre of the Earth, killing the Earth-L2 spatial baseline while
    // leaving the timing untouched (Step H1). This is the "off" half of Step H3's
    // satellite-parallax experiment, and it is also how the H1 regression proves the new term
    // is a clean no-op when disabled: with this flag the run must reproduce the pre-H1 output
    // byte for byte.
    bool   noSatPar    = false;
    // Step H3. Characterise every detected event TWICE -- once with Roman at L2, once with
    // the offset zeroed -- and write both forecasts to a side file. Two separate runs cannot
    // answer this: moving the observer changes which events are detected, so the detected
    // populations differ by more than the effect (DEVIATIONS.md 35).
    bool   pairSat     = false;

    // Step S1. Path to a sample-event dump spec, or empty for "do not dump". The dump is
    // off by default and consumes no RNG when on, so a run with it is event-for-event the
    // same run as one without -- see the block above main().
    std::string dumpSpec;

    // Build the sightline grid, write the provenance block, report the strata, and stop
    // before drawing a single star. The point of stratifying (Step E1) is to decide how to
    // spend wall clock, and that decision needs the sightline counts BEFORE committing to a
    // multi-hour run: the footprint stratum is the expensive one (a GBTDS sightline carries
    // ~50,000 Roman epochs against ~2,400 Rubin ones), so refining it by k multiplies its
    // sightlines by k^2 and the run time by rather more than that.
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
// Step E1: resolve --stride-roman against the detector size (Deviation 69), fill `steps`, and
// refuse a grid that steps over whole detectors. May set cfg.strideRoman. Returns 0, or the exit
// code of the error printed.
int resolveGridSteps(RunConfig& cfg, const GbtdsLayout& gl, GridSteps& steps);

#endif // ROMAN_RUN_CONFIG_H
