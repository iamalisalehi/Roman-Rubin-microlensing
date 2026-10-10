// The microlensing light curve and astrometric track over many epochs of one observer (evaluateModel,
// evaluateTrack), the observed peak, and the two-image resolution test.
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

// The single-epoch kernel behind the two ways of asking for the model below. Only lightcurve.cpp calls it.
void   lightcurve(source & s, lens & l, astromet & as, double, int tele);

// The code asks for the model in exactly two ways, and both take a whole light curve: one observer, many
// epochs, one call. A library that computes light curves (rather than single points) can replace the
// kernel behind them without any caller changing.
//
//   evaluateModel  the slim path for the Fisher matrices, which differentiate it hundreds of times per
//                  event and need only the four numbers of ModelPoint.
//   evaluateTrack  everything else (the time loop, the sample dump, the observed peak, the satellite
//                  offset), which need the whole geometry and are called once per light curve.

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

// What every other caller reads from the model at one epoch: the lens-source geometry with and without
// the microlensing parallax, the magnification, the centroid deflection, the source, centroid and lens
// positions, and the observer's offset.
struct TrackPoint {
    double u, uNoPlx;       //impact parameter WITH / WITHOUT the microlensing parallax [Einstein radii]
    double A, ANoPlx;       //magnifOf(u), magnifOf(uNoPlx)
    double def1c, def2c;    //centroid deflection WITH parallax [mas]
    double def1a, def2a;    //centroid deflection WITHOUT parallax [mas]
    double pos1b, pos2b;    //unlensed source position: proper motion + its own parallax [mas]
    double pos1c, pos2c;    //measured light centroid [mas]
    double lens1, lens2;    //lens position: its proper motion + its own parallax [mas]
    double ue1, ue2;        //the observer's sky-projected offset from Earth's position at t = 0 [AU] (the
                            //legacy as.ue_n1 / as.ue_n2); only the satellite separation du_sat reads it
};

// The current event seen by observer `tele` (0 = Rubin on Earth, 1 = Roman at L2) at the n epochs
// t[0..n-1]: out[i] is the model at t[i]. Uses no random numbers. See lightcurve.cpp.
void   evaluateTrack(source & s, lens & l, astromet & as, int tele, const double* t, int n,
                     TrackPoint* out);

// Finds the OBSERVED peak: the Earth-frame, parallax-bent closest approach.
// Returns {time of the peak, impact parameter there}. See the comment in lightcurve.cpp.
std::pair<double, double> observedPeak(source& s, lens& l, astromet& as);

#endif // ROMAN_EVENTS_LIGHTCURVE_H
