// Stage 1 of an event: draw the source and lens, bin them for the efficiency tables, reset the
// per-event accumulators, and run the per-survey pre-selection.
#ifndef ROMAN_SIM_DRAW_H
#define ROMAN_SIM_DRAW_H

#include "sim/state.h"

// func_source, func_lens, optical_depth, the seven efficiency-axis bins and their denominators, then
// the per-event resets. Draws random numbers. `prevNdw` is the epoch count of the PREVIOUS draw's light
// curve (LightCurveStats::ndw): the number of light-curve slots in ctx.l that it dirtied and this call clears.
EfficiencyBins drawEvent(SimContext& ctx, SightlineState& st, int prevNdw);

// The Step B2 pre-selection: is the event bright enough to be seen by Rubin (>= 2 bands) or by
// Roman, then one accept draw per survey. Returns true if either accepts, i.e. a light curve is
// to be generated. Draws random numbers (testL, testR).
bool preselectEvent(SimContext& ctx, const SightlineState& st);

#endif // ROMAN_SIM_DRAW_H
