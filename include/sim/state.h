// The state the per-event pipeline (src/sim/) shares: the long-lived objects, the per-sightline
// accumulators, and the small value structs that the stages of one event pass to each other.
#ifndef ROMAN_SIM_STATE_H
#define ROMAN_SIM_STATE_H

#include "common.h"
#include "types.h"
#include "run/config.h"
#include "run/outputs.h"
#include "run/summary.h"
#include "run/sightlines.h"
#include "galaxy/catalogue.h"
#include "galaxy/extinction.h"
#include "surveys/visits.h"
#include "surveys/footprints.h"
#include "surveys/schedule.h"
#include "fisher/fisher.h"

// Everything the pipeline reads or writes that outlives one event, by reference. Built once in
// main(); passed to every stage.
struct SimContext {
    const RunConfig&      cfg;
    const GbtdsLayout&    gl;
    const RomanSchedule&  sched;
    source&   s;                // the source star of the current draw
    lens&     l;                // the lens, and the detection-efficiency histograms
    astromet& as;               // astrometric scratch + Roman's observer offset (satScale)
    CMD&      cm;
    extin&    ex;
    lsst&     ls;               // Rubin visit list
    roman&    ro;               // Roman visit list
    covarian& co;               // Fisher results of the current event
    covarian& coNS;             // Step H3: the same event with the satellite offset zeroed
    RunOutputs& outs;           // output streams, file names, Step S1 dump state
    RunTotals&  run;            // run-wide counters, printed by writeRunSummary
};

// Per-sightline state: filled by setupSightline, grown by the draws, consumed by finishSightline.
struct SightlineState {
    // Position in the scan. iScan and lastCol persist across sightlines (resume bookkeeping).
    long   iScan   = -1;
    int    lastCol = -1;
    int    nri = -1, nde = -1;     // longitude column index, and index within that column
    double wArea = 0.0;            // sky area this sightline stands for [deg^2]
    int    sightlineIdx = 0;       // nearest extinction sightline
    // Epoch matching: how many visits of each survey cover this sightline, and the shortest gap.
    int    ndd = 0, nddR = 0;      // Rubin / Roman matched visits
    double minc = 0.0, mincR = 0.0;
    std::array<double, 6> rubinDepthMed{};   // median 5-sigma depth per LSST band (Deviation 73)
    // Event budget bookkeeping for the do/while over draws.
    int    icon = 0, nlens = 0;
    double nsim = 0.0, nerr = 0.0;
    std::array<int, NDETCLASS> nDetClass{};  // per-sightline detection-taxonomy counts
    std::vector<EventRecord> records;        // one row per draw, replayed by finishSightline
};

// The per-event values, one small struct per producing stage. Each is a plain value that the next stages
// take by const&, so a signature says what a stage reads. A default-constructed struct is the "nothing
// happened" value (no light curve, no detection, nothing measured).

// drawEvent: the bin indices for the seven detection-efficiency axes (gg = tE; Deviation 46).
struct EfficiencyBins { int gg = -1, ss = 0, qq = 0, ww = 0, vv = 0, zz = 0, pp = 0; };

// simulateLightCurve: what the time loop accumulated. All-zero (flagf = 0) if no light curve was generated.
struct LightCurveStats {
    int flagf = 0;                   // 1 if the light curve was generated at all
    // Epoch counts (ndw = joint total = ndw_L + ndw_R) and the run-test results.
    int ndw = 0, ndw_L = 0, ndw_R = 0;
    int flag_det_L = 0, flag_det_R = 0;
    // Step R1: epochs at which the two images were both detectable and far enough apart.
    long   nres5_L = 0, nres20_L = 0, nresPSF_L = 0, nres5_R = 0, nres20_R = 0, nresPSF_R = 0;
    double dsepMax_L = -1.0, dsepMax_R = -1.0;   // largest separation while both detectable [mas]
    // Chi-squared accumulators: chi1 = lensing model, chi2 = no-parallax model, chi3 = baseline;
    // the "a" versions are astrometric. Unsuffixed = joint, _L = Rubin only, _R = Roman only.
    double chi1 = 0, chi2 = 0, chi3 = 0, chi1a = 0, chi2a = 0, chi3a = 0;
    double chi1_L = 0, chi2_L = 0, chi3_L = 0, chi1a_L = 0, chi2a_L = 0, chi3a_L = 0;
    double chi1_R = 0, chi2_R = 0, chi3_R = 0, chi1a_R = 0, chi2a_R = 0, chi3a_R = 0;
    double vsave = 0.0;              // sum over the epochs of the source proper-motion speed
};

// detectEvent (verdicts, delta-chi-squared) and tallyDetection (dclsEvent).
struct Detection {
    double dchiL = 0, dchiP = 0, dchiA = 0;                 // signed statistics, joint
    double dchiL_L = 0, dchiP_L = 0, dchiA_L = 0;           // Rubin only
    double dchiL_R = 0, dchiP_R = 0, dchiA_R = 0;           // Roman only
    int detL = 0, detR = 0, detJ = 0, detJ_raw = 0;         // Rubin alone, Roman alone, joint (as patched), joint (raw)
    int dclsEvent = DET_NONE;        // DetClass of this draw; stays NONE if never tallied
};

// characterizeEvent, Step H3: the no-satellite forecast of the same event. -1 = not measured.
struct SatellitePair {
    double sigtE_ns = -1.0, sigpiE_ns = -1.0, sigpiER_ns = -1.0;
    double sigtetE_ns = -1.0, sigpiEb_ns = -1.0, relMl_ns = -1.0;
    double condA_ns = -1.0, condB_ns = -1.0;
    int    okNS = 0, okNSb = 0;
};

// characterizeEvent: the detection verdicts, whether FisherM ran, the Step H3 pair.
struct Characterization {
    Detection     det;
    SatellitePair pair;
    int    fisher = 0;               // 1 if the event was detected and FisherM ran
    double vMean = 0.0;              // mean source proper-motion speed (s.mus if the event is not visible)
};

// recordEvent (Step H2, Deviation 76): the observed peak, the satellite observable and the peak-window
// coverage; handed to commitSampleDump and writeSatellitePair.
struct PeakCoverage {
    double t0obs = 0.0, uminObs = 0.0;   // observed peak time and impact parameter
    double duSat = 0.0;                  // observer separation in Einstein radii at t0
    int    nepLpk = 0, nepRpk = 0;       // Rubin / Roman epochs within +-2 tE of the peak
};

#endif // ROMAN_SIM_STATE_H
