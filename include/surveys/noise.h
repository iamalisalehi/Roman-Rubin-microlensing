// Per-exposure measurement-error models: LSST photometry/astrometry, Roman photometry/astrometry.
#ifndef ROMAN_SURVEYS_NOISE_H
#define ROMAN_SURVEYS_NOISE_H

#include "common.h"
#include "surveys/visits.h"

double errlsstM(double,int,double);
double errlsstA(lsst & ls,  double);
double errRomanA(double magF146);   //Roman WFI per-exposure astrometric error [mas] (Step H4)

//double errELT(lsst & ls,double,int);
// TODO(Ali): wire this to whatever sigma_roman.txt actually represents (a fixed
// mag-vs-error lookup, or a per-visit-depth-dependent formula like errlsstM). Signature
// below assumes the simpler case (no per-visit depth); adjust if you need `sig5`.
double errRomanM(const roman & ro, double mag);

#endif // ROMAN_SURVEYS_NOISE_H
