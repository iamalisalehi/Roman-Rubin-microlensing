// Lens-source relative velocity and proper motion.
#ifndef ROMAN_GALAXY_KINEMATICS_H
#define ROMAN_GALAXY_KINEMATICS_H

#include "common.h"
#include "types.h"

void   vrel(source & s, lens & l);
// Asymmetric drift [km/s] of thin-disc age bin / thick-disc episode idx at radius R [kpc]; 0 for other components.
double discAsymmetricDrift(GalacticComponent comp, int idx, double R);

#endif // ROMAN_GALAXY_KINEMATICS_H
