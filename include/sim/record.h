// Stage 4 of an event: the row it leaves behind -- the EventRecord, the per-event table line, and the
// Step S1 sample light-curve dump.
#ifndef ROMAN_SIM_RECORD_H
#define ROMAN_SIM_RECORD_H

#include "sim/state.h"

// Pushes the EventRecord, computes the Step H2 satellite observable and peak-window coverage, writes the
// per-event table row, then commits the Step S1 dump and the Step H3 pair row. Draws no random numbers.
void recordEvent(SimContext& ctx, SightlineState& st, EventState& ev);

// Step S1: write this draw's buffered light curve if it fills a requested sample class.
void commitSampleDump(SimContext& ctx, const EventState& ev);

#endif // ROMAN_SIM_RECORD_H
