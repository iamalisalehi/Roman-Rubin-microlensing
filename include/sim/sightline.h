// Per-sightline stages of the Monte Carlo: set one sightline up (position, epoch matching, resets), and
// aggregate and write it out once its draws are done.
#ifndef ROMAN_SIM_SIGHTLINE_H
#define ROMAN_SIM_SIGHTLINE_H

#include "sim/state.h"

// What setupSightline decided: run the sightline, skip it (before --start-index, or no coverage), or
// stop the scan (--end-index reached).
enum class SightlineStart { Simulate, Skip, Stop };

// Position, the --start-index/--end-index bookkeeping, this sightline's RNG stream, epoch matching for
// Rubin and Roman, the per-band depths, and the resets the draws accumulate into.
SightlineStart setupSightline(SimContext& ctx, SightlineState& st, const Sightline& sightline);

// After the draws: budget check, the empty-sightline (barren) skip, the efficiency and map rows, and
// the per-sightline log. Returns early, writing nothing, for a barren sightline.
void finishSightline(SimContext& ctx, SightlineState& st);

#endif // ROMAN_SIM_SIGHTLINE_H
