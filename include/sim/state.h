// The state the per-event pipeline (src/sim/) shares: the long-lived objects, the per-sightline
// accumulators and the per-draw values that main()'s sightline loop used to keep as ~100 locals.
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

// Per-draw state. One instance lives for the whole run: ndw in particular must survive from one
// draw to the next, because drawEvent clears only the prefix of the light-curve buffers that the
// PREVIOUS event dirtied (and from one sightline to the next, for the same reason).
struct EventState {
    // Bin indices for the seven detection-efficiency axes (gg = tE; Deviation 46).
    int gg = -1, ss = 0, qq = 0, ww = 0, vv = 0, zz = 0, pp = 0;
    int    flagf = 0;                // 1 if the light curve was generated at all
    int    dclsEvent = DET_NONE;     // DetClass of this draw
    double initial = 0.0;            // time-window padding of the light-curve loop [days]

    // Epoch counts (ndw = joint total = ndw_L + ndw_R) and the run-test results.
    int ndw = 0, ndw_L = 0, ndw_R = 0;
    int flag_det = 0, flag_det_L = 0, flag_det_R = 0;
    // Step R1: epochs at which the two images were both detectable and far enough apart.
    long   nres5_L = 0, nres20_L = 0, nresPSF_L = 0, nres5_R = 0, nres20_R = 0, nresPSF_R = 0;
    double dsepMax_L = -1.0, dsepMax_R = -1.0;   // largest separation while both detectable [mas]

    // Chi-squared accumulators: chi1 = lensing model, chi2 = no-parallax model, chi3 = baseline;
    // the "a" versions are astrometric. Unsuffixed = joint, _L = Rubin only, _R = Roman only.
    double chi1 = 0, chi2 = 0, chi3 = 0, chi1a = 0, chi2a = 0, chi3a = 0;
    double chi1_L = 0, chi2_L = 0, chi3_L = 0, chi1a_L = 0, chi2a_L = 0, chi3a_L = 0;
    double chi1_R = 0, chi2_R = 0, chi3_R = 0, chi1a_R = 0, chi2a_R = 0, chi3a_R = 0;
    double dchiL = 0, dchiP = 0, dchiA = 0;                 // signed statistics, joint
    double dchiL_L = 0, dchiP_L = 0, dchiA_L = 0;           // Rubin only
    double dchiL_R = 0, dchiP_R = 0, dchiA_R = 0;           // Roman only
    double vsave = 0.0;              // running proper-motion-speed sum, then its mean

    // Detection verdicts (Rubin alone, Roman alone, joint); FFG[0] gates the Fisher call.
    int detL = 0, detR = 0, detJ = 0, detJ_raw = 0;
    std::array<int, 3> FFG{};

    // Step H3: the no-satellite forecast of the same event. -1 = not measured.
    double sigtE_ns = -1.0, sigpiE_ns = -1.0, sigpiER_ns = -1.0;
    double sigtetE_ns = -1.0, sigpiEb_ns = -1.0, relMl_ns = -1.0;
    double condA_ns = -1.0, condB_ns = -1.0;
    int    okNS = 0, okNSb = 0;

    // Step H2 and Deviation 76: satellite observable, peak-window coverage, observed peak.
    double duSat = 0.0;
    int    nepLpk = 0, nepRpk = 0;
    double t0obs = 0.0, uminObs = 0.0;
};

#endif // ROMAN_SIM_STATE_H
