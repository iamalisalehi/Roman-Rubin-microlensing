// Stage 2 of an event: the light curve as each survey records it -- the time loop over epochs with its
// Rubin (ugrizy) and Roman (F146) blocks. This is the hot path.
#ifndef ROMAN_SIM_OBSERVE_H
#define ROMAN_SIM_OBSERVE_H

#include "sim/state.h"

// Steps the event through time, taking Rubin and Roman epochs as they fall due, drawing the noisy
// photometry and astrometry and accumulating the chi-squared sums, epoch counts, resolution tallies
// and the sample-dump buffer. Fills the light-curve arrays of `l` up to the returned ndw. Draws random numbers.
LightCurveStats simulateLightCurve(SimContext& ctx, const SightlineState& st);

// The epochs of the time loop on this sightline, strictly increasing [simulation days]: the times simulateLightCurve
// lays out for every event (the grid depends only on what setupSightline matched for the sightline, not on the event),
// and the times the Fisher matrices then see (they are copies of them). Needs ls and ro matched to the sightline.
// setupSightline hands them to the light-curve model (vbmSetObserverGrid).
std::vector<double> sightlineEpochTimes(const SightlineState& st, const lsst& ls, const roman& ro);

#endif // ROMAN_SIM_OBSERVE_H
