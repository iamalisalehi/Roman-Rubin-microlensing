// Stellar catalogues: the Besancon/CMD component tables and the luminous-lens mass-luminosity table.
#ifndef ROMAN_GALAXY_CATALOGUE_H
#define ROMAN_GALAXY_CATALOGUE_H

#include "common.h"

struct CMD {
    // Thin disk
    std::vector<double> logT_thin;
    std::vector<double> mass_thin;
    std::vector<std::array<double, M>> Mab_thin; // M × N1
    std::vector<double> typ_thin;
    std::vector<double> cl_thin;
    std::vector<double> age_thin;

    // Bulge
    std::vector<double> logT_bulge;
    std::vector<double> mass_bulge;
    std::vector<std::array<double, M>> Mab_bulge; // M × N2
    std::vector<double> typ_bulge; 
    std::vector<double> cl_bulge;
    std::vector<double> age_bulge;

    // Thick disk
    std::vector<double> logT_thick;
    std::vector<double> mass_thick;
    std::vector<std::array<double, M>> Mab_thick; // M × N3
    std::vector<double> typ_thick;
    std::vector<double> cl_thick;
    std::vector<double> age_thick;

    // Halo
    std::vector<double> logT_halo;
    std::vector<double> mass_halo;
    std::vector<std::array<double, M>> Mab_halo; // M × N4
    std::vector<double> typ_halo;
    std::vector<double> cl_halo;
    std::vector<double> age_halo;

    // Constructor
    CMD()
        : logT_thin(N1),  mass_thin(N1),  Mab_thin(N1),  typ_thin(N1),  cl_thin(N1),  age_thin(N1),
          logT_bulge(N2), mass_bulge(N2), Mab_bulge(N2), typ_bulge(N2), cl_bulge(N2), age_bulge(N2),
          logT_thick(N3), mass_thick(N3), Mab_thick(N3), typ_thick(N3), cl_thick(N3), age_thick(N3),
          logT_halo(N4),  mass_halo(N4),  Mab_halo(N4),  typ_halo(N4),  cl_halo(N4),  age_halo(N4)
    {}
};

void   read_cmd(CMD & cm);

// Luminous lenses: main-sequence mass -> absolute magnitude (AB, ugrizy + F146) per
// Galactic component, from CMD/components/lens_ml.dat (CMD/lens_ml_table.py). Bins of 0.02 Msun
// over 0.08-1.00 Msun; lighter lenses are brown dwarfs and treated as dark.
struct LensML {
    std::array<std::vector<double>, 4> mmid;
    std::array<std::vector<std::array<double, 7>>, 4> mab;
};
inline LensML gLensML;
void   readLensML(const std::string& path);
bool   lensAbsMag(int comp, double mass, std::array<double, 7>& mab);   // false if dark

#endif // ROMAN_GALAXY_CATALOGUE_H
