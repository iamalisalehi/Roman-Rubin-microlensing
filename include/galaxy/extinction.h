// 3D dust extinction: the A_V(distance) tables, sightline lookup, and the CCM89 reddening law.
#ifndef ROMAN_GALAXY_EXTINCTION_H
#define ROMAN_GALAXY_EXTINCTION_H

#include "common.h"

// The extinction tables (Deviation 70): files/ext/ext_tables.dat, built by maps.py from DECaPS and
// Marshall through dustmaps. One shared distance grid; per table position its (l, b) and A_V on that
// grid, row-major in `ext`. Sizes come from the file's header at run time -- no NFILES/NROWS to keep
// in step -- and A_V is stored as float (the maps' own errors are tenths of a magnitude), so the
// ~16,000 positions x 399 distances take ~25 MB, against 148 MB for the 2,518 per-pointing tables
// with their own distance columns before.
struct extin
{
    std::vector<double> dist;    // n_dist, kpc, strictly increasing
    std::vector<double> l, b;    // n_tables
    std::vector<float>  ext;     // n_tables * n_dist, A_V [mag]
    int nTables = 0, nDist = 0;
    std::string built, k;        // provenance, from the header ("# ext_tables ...", "# k ...")
};

int    nearestSightline(const extin& ex, double lon, double lat);
void   readExtinction(extin& ex, const std::string& path); // files/ext/ext_tables.dat; exits on bad input
double interpExtinctionAlongSightline(const extin& ex, int k, double dist);

double CCM89_a(double lambda_um);
double CCM89_b(double lambda_um);
double AlAv(double lambda_um, double Rv);

#endif // ROMAN_GALAXY_EXTINCTION_H
