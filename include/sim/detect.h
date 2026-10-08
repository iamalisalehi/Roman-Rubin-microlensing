// Stage 3a of an event: the detection test and its bookkeeping.
#ifndef ROMAN_SIM_DETECT_H
#define ROMAN_SIM_DETECT_H

#include "sim/state.h"

// Signed delta-chi-squared per survey and jointly, the three detection verdicts
// (detL, detR, detJ) and the monotone patch; detJ_raw is the joint verdict before the patch.
Detection detectEvent(SimContext& ctx, const LightCurveStats& lc);

// After the Fisher call: the DetClass taxonomy, the delta-chi mismatch and anomaly diagnostics, and
// the numerators of the seven detection-efficiency axes. Sets det.dclsEvent.
void tallyDetection(SimContext& ctx, SightlineState& st, const EfficiencyBins& bins, const LightCurveStats& lc,
                    Detection& det);

#endif // ROMAN_SIM_DETECT_H
