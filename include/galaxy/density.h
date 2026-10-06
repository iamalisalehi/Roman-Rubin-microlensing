// Galactic stellar density model (thin disk, bulge, thick disk, halo) and the optical depth it implies.
#ifndef ROMAN_GALAXY_DENSITY_H
#define ROMAN_GALAXY_DENSITY_H

#include "common.h"
#include "types.h"

void   optical_depth(source & s);
void   Disk_model(source & s, int);

#endif // ROMAN_GALAXY_DENSITY_H
