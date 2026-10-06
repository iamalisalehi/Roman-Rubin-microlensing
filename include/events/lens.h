// Drawing a lens (position, mass, kinematics) and the lens mass functions.
#ifndef ROMAN_EVENTS_LENS_H
#define ROMAN_EVENTS_LENS_H

#include "common.h"
#include "types.h"
#include "galaxy/extinction.h"

void   func_lens( lens & l, source & s, const extin & ex, int sightlineIdx);

double drawKroupaInitialMass();
double remnantMass(double initialMass);
double drawLogUniformMass(double lo, double hi);
double drawNeutronStarMass();
double drawPowerLawMass(double lo, double hi, double alpha);
double drawLensMass(bool* luminous = nullptr);   // *luminous: a living star (Deviation 74)

#endif // ROMAN_EVENTS_LENS_H
