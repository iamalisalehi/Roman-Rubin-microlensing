// Stage 4 of an event: the row it leaves behind -- the EventRecord, the per-event table line, and the
// sample light-curve dump.
#ifndef ROMAN_SIM_RECORD_H
#define ROMAN_SIM_RECORD_H

#include "sim/state.h"

// Pushes the EventRecord, computes the satellite observable and peak-window coverage, writes the
// per-event table row, then commits the sample dump and the satellite-pair row. Draws no random numbers.
void recordEvent(SimContext& ctx, SightlineState& st, const EfficiencyBins& bins, const LightCurveStats& lc,
                 const Characterization& ch);

// Write this draw's buffered light curve if it fills a requested sample class.
void commitSampleDump(SimContext& ctx, const LightCurveStats& lc, const Characterization& ch,
                      const PeakCoverage& pk);

#endif // ROMAN_SIM_RECORD_H
