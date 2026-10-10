// Per-exposure measurement-error models: LSST photometry/astrometry, Roman photometry/astrometry.
#ifndef ROMAN_SURVEYS_NOISE_H
#define ROMAN_SURVEYS_NOISE_H

#include "common.h"
#include "surveys/visits.h"

double errlsstM(double,int,double);
double errlsstA(double errPhotMag, double fwhmArcsec);   //Rubin per-visit astrometric error [mas]
double errRomanA(double errPhotMag);   //Roman WFI per-exposure astrometric error [mas]

// TODO: confirm what sigma_roman.txt represents (a fixed mag-vs-error lookup, or a per-visit-depth
// formula like errlsstM). The signature assumes the simpler case (no per-visit depth).
double errRomanM(const roman & ro, double mag);

#endif // ROMAN_SURVEYS_NOISE_H
