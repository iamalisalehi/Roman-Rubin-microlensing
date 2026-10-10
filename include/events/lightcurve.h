// The microlensing light curve and astrometric track at one epoch and over many (evaluateModel), the observed peak,
// and the two-image resolution test.
#ifndef ROMAN_EVENTS_LIGHTCURVE_H
#define ROMAN_EVENTS_LIGHTCURVE_H

#include "common.h"
#include "types.h"

// The two images of one source at impact parameter u, with their own (unblended) magnitudes.
// `magBase`/`blendFrac` are the BLENDED baseline magnitude and the source's flux share in the
// filter being observed, i.e. s.magb[i] and s.blend[i] -- per FILTER, not per telescope.
struct ImagePair {
    double sep;            //angular separation of the two images [mas]; -1 if undefined
    double magPlus;        //apparent magnitude of the major image [mag]
    double magMinus;       //apparent magnitude of the minor image [mag]
    bool   bothDetectable; //both within [saturation, single-visit depth] in this filter
};

inline ImagePair imagePair(double u, double tetE, double magBase, double blendFrac,
                           double thrMag, double satMag)
{
    ImagePair ip{-1.0, 99.0, 99.0, false};
    if (!(u > 0.0) or !(tetE > 0.0) or !(blendFrac > 0.0)) return ip;

    const double root = std::sqrt(u * u + 4.0);
    ip.sep = tetE * root;                                   //[mas]

    const double Aplus  = (u * u + 2.0) / (2.0 * u * root) + 0.5;
    const double Aminus = (u * u + 2.0) / (2.0 * u * root) - 0.5;
    if (!(Aminus > 0.0) or !(Aplus > 0.0)) return ip;       //minor image formally extinguished

    // Source magnitude from the blended baseline: m_source = m_base - 2.5 log10(blendFrac).
    // Blend light is deliberately not added back, so a faint minor image is not made to look
    // detectable on light that is not its own.
    const double magSource = magBase - 2.5 * std::log10(blendFrac);
    ip.magPlus  = magSource - 2.5 * std::log10(Aplus);
    ip.magMinus = magSource - 2.5 * std::log10(Aminus);

    ip.bothDetectable = (ip.magPlus  <= thrMag and ip.magPlus  >= satMag and
                         ip.magMinus <= thrMag and ip.magMinus >= satMag);
    return ip;
}

// A(u) = (u^2 + 2) / (u * sqrt(u^2 + 4)), the point-source point-lens magnification.
// Written in the form used inside the time loop so dumped values match the simulation bit for bit.
inline double magnifOf(double u)
{
    return double(u * u + 2.0) / std::sqrt(u * u * (u * u + 4.0));
}

void   lightcurve(source & s, lens & l, astromet & as, double, int tele);

// What the Fisher matrices need from the model at one epoch.
struct ModelPoint {
    double A;             //magnification, magnifOf(u)
    double u;             //impact parameter, s.ut
    double pos1c, pos2c;  //measured light centroid [mas], s.pos1c / s.pos2c
};

// The current event (s, l, as as they stand) seen by observer `tele` at the n epochs t[0..n-1]:
// out[i] is the model at t[i]. The one place FisherM asks for the model. See lightcurve.cpp.
void   evaluateModel(source & s, lens & l, astromet & as, int tele, const double* t, int n,
                     ModelPoint* out);

// Finds the OBSERVED peak: the Earth-frame, parallax-bent closest approach.
// Returns {time of the peak, impact parameter there}. See the comment in lightcurve.cpp.
std::pair<double, double> observedPeak(source& s, lens& l, astromet& as);

#endif // ROMAN_EVENTS_LIGHTCURVE_H
