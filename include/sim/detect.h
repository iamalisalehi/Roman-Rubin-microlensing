// Stage 3a of an event: the detection test and its bookkeeping.
#ifndef ROMAN_SIM_DETECT_H
#define ROMAN_SIM_DETECT_H

#include "sim/state.h"

// Step H7: signed delta-chi-squared per survey and jointly, the three detection verdicts
// (ev.detL, ev.detR, ev.detJ) and the monotone patch. Sets ev.detJ_raw.
void detectEvent(SimContext& ctx, EventState& ev);

// After the Fisher call: the DetClass taxonomy, the delta-chi mismatch and anomaly diagnostics, and
// the numerators of the seven detection-efficiency axes.
void tallyDetection(SimContext& ctx, SightlineState& st, EventState& ev);

#endif // ROMAN_SIM_DETECT_H
