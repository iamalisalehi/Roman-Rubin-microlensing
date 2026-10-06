// Bin finders for the per-event histograms: which grid cell of lens::tEs/Mls/pis/... a drawn value falls in.
#ifndef ROMAN_RUN_HISTOGRAMS_H
#define ROMAN_RUN_HISTOGRAMS_H

#include "common.h"
#include "types.h"

int    Funcu0(lens & l);
int    FunctE(lens & l);
int    FuncMl(lens & l);
int    FuncPi(lens & l);
int    FuncMu(lens & l);
int    FuncMb(lens & l, double);
int    FuncFb(lens & l, double);

#endif // ROMAN_RUN_HISTOGRAMS_H
