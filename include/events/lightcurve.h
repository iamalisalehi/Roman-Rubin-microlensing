// The microlensing light curve and astrometric track at one epoch, the observed peak, and the two-image resolution test.
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

    // The SOURCE's own magnitude, recovered from the blended baseline: blendFrac is the
    // source's share of the aperture flux, so m_source = m_base - 2.5 log10(blendFrac) and is
    // always the FAINTER of the two (blendFrac <= 1). Each image then carries its own
    // magnification. The blend light is deliberately NOT added back: an image that is resolved
    // from its twin is resolved from the neighbours too, and re-adding the full blend would
    // make a faint minor image look detectable on light that is not its own.
    const double magSource = magBase - 2.5 * std::log10(blendFrac);
    ip.magPlus  = magSource - 2.5 * std::log10(Aplus);
    ip.magMinus = magSource - 2.5 * std::log10(Aminus);

    ip.bothDetectable = (ip.magPlus  <= thrMag and ip.magPlus  >= satMag and
                         ip.magMinus <= thrMag and ip.magMinus >= satMag);
    return ip;
}

// A(u) = (u^2 + 2) / (u * sqrt(u^2 + 4)), the point-source point-lens magnification.
// Written in exactly the form the two call sites inside the time loop use, so the
// dumped magnification is what the simulation computed and not an algebraically
// equal rearrangement that could round differently.
inline double magnifOf(double u)
{
    return double(u * u + 2.0) / std::sqrt(u * u * (u * u + 4.0));
}

void   lightcurve(source & s, lens & l, astromet & as, double, int tele);

// Finds the OBSERVED peak (Deviation 76): the Earth-frame, parallax-bent closest approach.
// Returns {time of the peak, impact parameter there}. See the comment in lightcurve.cpp.
std::pair<double, double> observedPeak(source& s, lens& l, astromet& as);

#endif // ROMAN_EVENTS_LIGHTCURVE_H
