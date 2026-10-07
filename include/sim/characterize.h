// Stage 3b of an event: Fisher-matrix characterisation of a detected event, and the Step H3
// no-satellite second forecast.
#ifndef ROMAN_SIM_CHARACTERIZE_H
#define ROMAN_SIM_CHARACTERIZE_H

#include "sim/state.h"

// Resets the Fisher results to the not-measured sentinel, then, for an observable event, runs the
// detection test (detectEvent), FisherM/ErrorCal, the LpLMC row, the Step H3 no-satellite forecast
// (if --pair-satellite) and tallyDetection. Returns the verdicts, the pair and the mean speed.
// Draws no random numbers.
Characterization characterizeEvent(SimContext& ctx, SightlineState& st, const EfficiencyBins& bins,
                                   const LightCurveStats& lc);

// Step H3: one row per detected event with both forecasts, to the side file. Needs the Step H2
// quantities (PeakCoverage) that recordEvent computes first.
void writeSatellitePair(SimContext& ctx, const SightlineState& st, const Characterization& ch,
                        const PeakCoverage& pk);

#endif // ROMAN_SIM_CHARACTERIZE_H
