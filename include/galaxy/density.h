// Galactic stellar density model (thin disk, bulge, thick disk, halo) and the optical depth it implies.
#ifndef ROMAN_GALAXY_DENSITY_H
#define ROMAN_GALAXY_DENSITY_H

#include "common.h"
#include "types.h"

void   optical_depth(source & s);
void   Disk_model(source & s, int);
double thinBinDensity(int ii, double Rb, double zb);  // Msun/pc^3 of age bin ii at (R, z) [kpc]
double thickEpisodeDensity(int k, double Rb, double zb);  // Msun/pc^3 of thick-disc episode k (0 young, 1 old)
double thickDensity(double Rb, double zb);            // Msun/pc^3 of the thick disc (both episodes)
double haloDensity(double Rb, double zb);             // Msun/pc^3 of the stellar halo

#endif // ROMAN_GALAXY_DENSITY_H
