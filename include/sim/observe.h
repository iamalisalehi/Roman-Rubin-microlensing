// Stage 2 of an event: the light curve as each survey records it -- the time loop over epochs with its
// Rubin (ugrizy) and Roman (F146) blocks. This is the hot path.
#ifndef ROMAN_SIM_OBSERVE_H
#define ROMAN_SIM_OBSERVE_H

#include "sim/state.h"

// Steps the event through time, taking Rubin and Roman epochs as they fall due, drawing the noisy
// photometry and astrometry and accumulating the chi-squared sums, epoch counts, resolution tallies
// (Step R1) and the Step S1 buffer. Fills the light-curve arrays of `l` up to the returned ndw. Draws random numbers.
LightCurveStats simulateLightCurve(SimContext& ctx, const SightlineState& st);

#endif // ROMAN_SIM_OBSERVE_H
