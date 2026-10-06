// Drawing a source star from the CMD and extinction data.
#ifndef ROMAN_EVENTS_SOURCE_H
#define ROMAN_EVENTS_SOURCE_H

#include "common.h"
#include "types.h"
#include "galaxy/catalogue.h"
#include "galaxy/extinction.h"

void   func_source(source & s, CMD & cm, const extin& ex, int sightlineIdx);

#endif // ROMAN_EVENTS_SOURCE_H
