// Drawing a lens (position, mass, kinematics) and the lens mass functions.
#ifndef ROMAN_EVENTS_LENS_H
#define ROMAN_EVENTS_LENS_H

#include "common.h"
#include "types.h"
#include "galaxy/extinction.h"
#include "galaxy/catalogue.h"

void   func_lens( lens & l, source & s, const CMD & cm, const extin & ex, int sightlineIdx);

double drawKroupaInitialMass();
double remnantMass(double initialMass);
double drawLogUniformMass(double lo, double hi);
double drawNeutronStarMass();
double drawPowerLawMass(double lo, double hi, double alpha);
double drawCatalogueLens(const CMD& cm, GalacticComponent comp, bool* luminous, std::array<double, 7>& mab);
double drawLensMass(bool* luminous = nullptr);   // *luminous: a living star

#endif // ROMAN_EVENTS_LENS_H
