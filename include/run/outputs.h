// The per-event output table: its column header.
#ifndef ROMAN_RUN_OUTPUTS_H
#define ROMAN_RUN_OUTPUTS_H

#include "common.h"

// Column names of the per-event table, in write order. Kept next to the writer's data so
// the two cannot drift: a header that disagrees with its columns is worse than none.
const char* eventTableHeader();

#endif // ROMAN_RUN_OUTPUTS_H
