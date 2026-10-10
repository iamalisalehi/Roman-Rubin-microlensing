// Light curve and astrometric track of one event at one epoch and over many (evaluateModel), and
// the observed-peak finder.
#include "events/lightcurve.h"

// tele: which observatory is asking -- 0 = Rubin (on Earth), 1 = Roman (at Sun-Earth L2).
void lightcurve(source & s, lens & l, astromet & as, double timh, int tele)
{
    double dvex, dvey, Ve_x, tt;
    double abs1 = 0.0, abs2 = 0.0; //undifferenced projection at t = timh
    double pis = double(1.0/s.Ds); //[mas] source parallax
    double pil = double(1.0/l.Dl); //[mas] lens parallax
    double int1 = 0.0, int2 = 0.0;


    for (int ig = 0; ig < 2; ++ig) {
        if(ig == 0) tt = timh;
        if(ig == 1) tt = 0.0;
        
        dvex = +vearth * std::sin(omegae * tt + M_PI / 2.0) / omegae;//km/s
        dvey = -vearth * std::cos(omegae * tt + M_PI / 2.0) / omegae;//km/s
        as.Ve_n1 =  std::cos(tetp) * dvex * std::sin(l.deltao) - dvey * std::cos(l.deltao);
        Ve_x     = -std::cos(tetp) * dvex * std::cos(l.deltao) - dvey * std::sin(l.deltao);
        as.Ve_n2 = -std::sin(s.FI) * Ve_x + std::cos(s.FI) * std::sin(tetp) * dvex;

        if(ig == 0) { int1  = as.Ve_n1; int2  = as.Ve_n2;
                      abs1  = as.Ve_n1; abs2  = as.Ve_n2; } //kept UNdifferenced
        if(ig == 1) { int1 -= as.Ve_n1; int2 -= as.Ve_n2; }
    }

    // Satellite-parallax term. The loop above computes ue(t) = P(X_E(t)) - P(X_E(0)): the observer's
    // projected displacement measured from Earth's position at t = 0, which fixes the meaning of u0
    // and t0. Roman sits on the Sun-Earth line, X_R(t) = (1 + f) X_E(t) with f = L2_OFFSET_AU ~ 0.01,
    // and P is linear in (dvex, dvey), so with the same single origin for both observers
    //
    //     ue_Roman(t) = ue_Rubin(t) + f * P(X_E(t)),
    //
    // i.e. the differenced term plus f times the UNdifferenced projection (abs1/abs2).
    // Letting each observer subtract its own t = 0 position instead would only rescale Earth's
    // parallax ellipse by (1+f) and drop the constant inter-observer offset, which is the satellite
    // parallax. Check: at t = 0 the two observers must differ, |delta u| ~ f * piE ~ 1e-3.
    // Rubin has fSat = 0.
    const double fSat = (tele == 1) ? L2_OFFSET_AU * as.satScale : 0.0;

    as.ue_n1 = int1 + fSat * abs1; //[radian] without dimention
    as.ue_n2 = int2 + fSat * abs2; //[radian] without dimention

    s.ux=-l.u0*std::sin(s.xi) + (timh-l.t0)*std::cos(s.xi)/l.tE + l.piE*as.ue_n1;//[] +parallax
    s.uy= l.u0*std::cos(s.xi) + (timh-l.t0)*std::sin(s.xi)/l.tE + l.piE*as.ue_n2;//[] +parallax

    s.ut0=std::sqrt((s.ux-l.piE*as.ue_n1)*(s.ux-l.piE*as.ue_n1) + (s.uy-l.piE*as.ue_n2)*(s.uy-l.piE*as.ue_n2));
    s.ut= std::sqrt(s.ux*s.ux+ s.uy*s.uy);

    if(s.ut==0.0)   s.ut=1.0e-50;
    if(s.ut0==0.0)  s.ut0=1.0e-50;


    s.def1a = (s.ux - l.piE * as.ue_n1) * l.tetE / (s.ut0 * s.ut0 + 2.0); //x-deflection, without parallax[mas]
    s.def2a = (s.uy - l.piE * as.ue_n2) * l.tetE / (s.ut0 * s.ut0 + 2.0); //y-deflection, without parallax[mas]

    s.def1c = s.ux * l.tetE / (s.ut * s.ut + 2.0); //x-deflection[mas]
    s.def2c = s.uy * l.tetE / (s.ut * s.ut + 2.0); //y-deflection[mas]

    s.pos1b = -l.u0 * l.tetE * std::sin(s.xi) + s.mus1 * (timh - l.t0) - as.ue_n1 * pis; //x-source trajectory+parallax[mas]
    s.pos2b = +l.u0 * l.tetE * std::cos(s.xi) + s.mus2 * (timh - l.t0) - as.ue_n2 * pis; //y-source trajectory+parallax[mas]

    // The MEASURED centroid: the light-weighted position of everything in the PSF. The lensed source
    // (flux fb*A, at its unlensed position + the deflection), the lens's own light (fLens, at the
    // lens: -u*thetaE from the source) and the unresolved neighbours (the rest of the baseline flux,
    // at their light centroid blendOff, which moves with the source). As the source brightens its
    // share grows and the centroid slides from the neighbours toward it. Fractions are of this
    // telescope's baseline flux.
    {
        const int    tt  = (tele == 1) ? 1 : 0;
        const double fs  = s.fb[tt], fL = s.fLens[tt];
        const double fn  = std::max(0.0, 1.0 - fs - fL);
        const double u2  = s.ut * s.ut;
        const double A   = (u2 + 2.0) / std::sqrt(u2 * (u2 + 4.0));
        const double den = fs * A + 1.0 - fs;
        s.pos1c = s.pos1b + (fs * A * s.def1c - fL * l.tetE * s.ux + fn * s.blendOff[tt][0]) / den; //x-centroid[mas]
        s.pos2c = s.pos2b + (fs * A * s.def2c - fL * l.tetE * s.uy + fn * s.blendOff[tt][1]) / den; //y-centroid[mas]
    }

    l.pos1  = l.mul1 * (timh - l.t0) - as.ue_n1 * pil ;//x-lens trajectory && parallax[mas]
    l.pos2  = l.mul2 * (timh - l.t0) - as.ue_n2 * pil ;//y-lens trajectory && parallax[mas]
}

// The model for one parameter set at all of one telescope's epochs: the Fisher matrices ask for it
// through here and nowhere else, one call per parameter set and telescope. A whole light curve per
// call is the shape of a library that computes light curves (one observer, many epochs) rather than
// single points, which can replace the loop below without FisherM changing. For now it is
// lightcurve() epoch by epoch, and leaves s/as holding the state at t[n-1].
void evaluateModel(source & s, lens & l, astromet & as, int tele, const double* t, int n,
                   ModelPoint* out)
{
    for (int i = 0; i < n; ++i) {
        lightcurve(s, l, as, t[i], tele);
        out[i] = ModelPoint{magnifOf(s.ut), s.ut, s.pos1c, s.pos2c};
    }
}

// The OBSERVED peak. The table's t0 and u0 are the closest approach of the straight line in
// lightcurve()'s gauge (parallax referenced to Earth's position at t = 0); the event an observer
// sees is the parallax-bent trajectory, whose closest approach comes at another time (median offset
// ~0.12 tE on a sample of events). Found for the Earth observer on a grid of 600 points over
// t0 +- 3 tE, then golden-section refinement. Uses no random numbers; leaves s/as holding the state
// at the returned time (callers recompute what they need).
std::pair<double, double> observedPeak(source& s, lens& l, astromet& as)
{
    auto uAt = [&](double t) { lightcurve(s, l, as, t, 0); return s.ut; };
    const double lo = l.t0 - 3.0 * l.tE, hi = l.t0 + 3.0 * l.tE;
    const int    N  = 600;
    const double h  = (hi - lo) / N;
    int best = 0; double ubest = uAt(lo);
    for (int k = 1; k <= N; ++k) { const double u = uAt(lo + k * h); if (u < ubest) { ubest = u; best = k; } }
    double a = lo + std::max(best - 1, 0) * h, b = lo + std::min(best + 1, N) * h;
    const double g = 0.5 * (std::sqrt(5.0) - 1.0);
    double c = b - g * (b - a), d = a + g * (b - a), fc = uAt(c), fd = uAt(d);
    for (int it = 0; it < 60 and (b - a) > 1e-6 * l.tE; ++it) {
        if (fc < fd) { b = d; d = c; fd = fc; c = b - g * (b - a); fc = uAt(c); }
        else         { a = c; c = d; fc = fd; d = a + g * (b - a); fd = uAt(d); }
    }
    const double t = 0.5 * (a + b);
    return {t, uAt(t)};
}
