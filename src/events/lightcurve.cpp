// Light curve and astrometric track of one event over many epochs (evaluateModel, evaluateTrack), and
// the observed-peak finder. The model itself is the library's, behind events/vbm_model.h.
#include "events/lightcurve.h"
#include "events/vbm_model.h"

namespace {

// The event as the model takes it, from the simulator's structs: the source's size is s.rho, and its model
// (a point or a uniform disc) is the run's (vbmSourceModel).
ModelEvent eventOf(const source & s, const lens & l)
{
    return ModelEvent{l.u0, l.tE, l.t0, l.piE, s.xi, l.tetE, s.mus1, s.mus2, 1.0 / s.Ds, s.rho, vbmSourceModel()};
}

// Whether observer `tele` (0 = Rubin on Earth, 1 = Roman) sits at Roman's L2 orbit. as.satScale is a
// switch, not a scale: 1 = Roman at L2 (physical), 0 = Roman at the centre of the Earth, for the
// satellite-parallax comparison (--no-satellite-parallax, --pair-satellite).
bool atL2(const astromet & as, int tele)
{
    CHECK(tele == 0 or tele == 1);
    CHECK(as.satScale == 0.0 or as.satScale == 1.0);
    return tele == 1 and as.satScale == 1.0;
}

// The model's epochs, kept between calls so the hot path does not allocate: [0] the photometric call's,
// [1] the astrometric call's (evaluateTrack with positions).
std::vector<ModelEpoch>& epochScratch(int n, int which = 0)
{
    static thread_local std::vector<ModelEpoch> v[2];
    if (int(v[which].size()) < n) v[which].resize(n);
    return v[which];
}

// The MEASURED centroid: the light-weighted position of everything in the PSF. The lensed source
// (flux fb*A, at its unlensed position + the deflection), the lens's own light (fLens, at the
// lens: -u*thetaE from the source) and the unresolved neighbours (the rest of the baseline flux,
// at their light centroid blendOff, which moves with the source). As the source brightens its
// share grows and the centroid slides from the neighbours toward it. Fractions are of this
// telescope's baseline flux.
// The offsets from the unlensed source are what the light-weighting acts on: def = (image centroid) -
// (unlensed source), and -thetaE (ux, uy) = (lens) - (source).
void measuredCentroid(const source & s, const lens & l, int tele, const ModelEpoch & e,
                      double& def1c, double& def2c, double& pos1c, double& pos2c)
{
    const int    tt  = (tele == 1) ? 1 : 0;
    const double fs  = s.fb[tt], fL = s.fLens[tt];
    const double fn  = std::max(0.0, 1.0 - fs - fL);
    const double den = fs * e.A + 1.0 - fs;
    def1c = e.imgx - e.srcx;
    def2c = e.imgy - e.srcy;
    pos1c = e.srcx + (fs * e.A * def1c - fL * l.tetE * e.ux + fn * s.blendOff[tt][0]) / den; //x-centroid[mas]
    pos2c = e.srcy + (fs * e.A * def2c - fL * l.tetE * e.uy + fn * s.blendOff[tt][1]) / den; //y-centroid[mas]
}

}  // namespace

// The model for one parameter set at all of one telescope's epochs: the Fisher matrices ask for it
// through here and nowhere else, one call per parameter set and telescope -- one call to the library
// for the whole array. The photometric Fisher pass needs only the magnification, which the library
// computes without the lens's motion and the image centroids; the astrometric pass asks for them.
void evaluateModel(const source & s, const lens & l, const astromet & as, int tele, const double* t, int n,
                   bool centroids, ModelPoint* out)
{
    std::vector<ModelEpoch>& e = epochScratch(n);
    vbmEvaluate(eventOf(s, l), atL2(as, tele), centroids, t, n, e.data());
    for (int i = 0; i < n; ++i) {
        out[i] = ModelPoint{e[i].A, e[i].u, 0.0, 0.0};
        if (centroids) {
            double d1, d2;
            measuredCentroid(s, l, tele, e[i], d1, d2, out[i].pos1c, out[i].pos2c);
        }
    }
}

// The same light curve for everything that is not the Fisher matrix, with the whole geometry per
// epoch instead of four numbers: the time loop, the sample dump, the observed peak and the satellite
// offset ask for it through here and nowhere else. Same shape as evaluateModel above.
//
// u, A, ux, uy and everything computed from them are the photometric call's, whether or not positions are
// asked for: the library's astrometric call perturbs the parallax vector at every epoch of an array
// (vbm_model.cpp), which would make a run depend on whether the positions are written. The positions (the
// source, the lens and the centroid, whose difference is the deflection def1c, def2c) come from a second
// call over the same epochs, the astrometric one, and cost about as much as the first. Without them the
// positions are NaN (so that a position read by mistake is loud) and the deflection is the point-lens
// expression from the source-lens vector: the images' centroid lies thetaE u/(u^2 + 2) from the source,
// along (ux, uy), exact for a point lens. Under a finite source the library's centroid differs from it while
// the source is within a few radii of the lens (at the peak of an event that crosses the disc by half the
// shift and in the opposite direction).
//
// The "no parallax" values are the geocentric straight line the light curve is parametrised by,
//     u_np = u0 (-sin xi, cos xi) + tau (cos xi, sin xi),   tau = (t - t0) / tE,
// whose length is sqrt(u0^2 + tau^2). In the geocentric frame it is the tangent to the true track at t0
// (the Earth at rest there), the natural model without parallax; the observer's departure from it, the
// Earth's orbit curving away from its velocity at t0, is the annual parallax, and Roman's offset from the
// Earth is added on top. Its deflection is thetaE u_np / (u_np^2 + 2).
void evaluateTrack(const source & s, const lens & l, const astromet & as, int tele, const double* t, int n,
                   bool positions, TrackPoint* out)
{
    std::vector<ModelEpoch>& e = epochScratch(n);
    const ModelEvent ev = eventOf(s, l);
    vbmEvaluate(ev, atL2(as, tele), false, t, n, e.data());
    const ModelEpoch* ea = nullptr;     // the astrometric call's epochs
    if (positions) {
        vbmEvaluate(ev, atL2(as, tele), true, t, n, epochScratch(n, 1).data());
        ea = epochScratch(n, 1).data();
    }
    const double sx = std::sin(s.xi), cx = std::cos(s.xi);
    constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
    for (int i = 0; i < n; ++i) {
        const ModelEpoch & m = e[i];
        TrackPoint & p = out[i];
        p.u = m.u;
        p.A = m.A;
        p.ux = m.ux;
        p.uy = m.uy;

        const double tau = (t[i] - l.t0) / l.tE;
        const double vx = -l.u0 * sx + tau * cx, vy = l.u0 * cx + tau * sx;
        double uNp = std::sqrt(l.u0 * l.u0 + tau * tau);
        if (uNp == 0.0) uNp = 1.0e-50;
        p.uNoPlx = uNp;
        p.ANoPlx = vbmMagnification(ev, uNp);   // the same source, on the straight line
        p.def1a = vx * l.tetE / (uNp * uNp + 2.0);   //x-deflection, without parallax[mas]
        p.def2a = vy * l.tetE / (uNp * uNp + 2.0);   //y-deflection, without parallax[mas]

        if (positions) {
            // the centroid, in the astrometric call's own frame and with its own A, ux, uy
            measuredCentroid(s, l, tele, ea[i], p.def1c, p.def2c, p.pos1c, p.pos2c);
            p.pos1b = ea[i].srcx;   //x-source trajectory+parallax[mas]
            p.pos2b = ea[i].srcy;
            p.lens1 = ea[i].lensx;  //x-lens trajectory && parallax[mas]
            p.lens2 = ea[i].lensy;
        } else {
            const double k = l.tetE / (m.u * m.u + 2.0);
            p.def1c = k * m.ux;
            p.def2c = k * m.uy;
            p.pos1b = p.pos2b = p.pos1c = p.pos2c = p.lens1 = p.lens2 = NaN;
        }
    }
}

// The OBSERVED peak. The table's t0 and u0 are the closest approach of the geocentric straight line;
// the event an observer sees is the parallax-bent trajectory, whose closest approach comes at another
// time (the Earth's orbit curves away from the line it has at t0). Found for the Earth observer on a grid
// of 600 points over t0 +- 3 tE, then golden-section refinement. The window is cut to the days the Earth's
// ephemeris covers (the library would extrapolate it): it matters only for events many years long, whose
// peak is not defined by the survey anyway. Uses no random numbers; changes nothing.
std::pair<double, double> observedPeak(const source& s, const lens& l, const astromet& as)
{
    // u(t) at one time; the golden-section search below asks for one point at a time.
    auto uAt = [&](double t) { ModelPoint p; evaluateModel(s, l, as, 0, &t, 1, false, &p); return p.u; };
    double first, last;
    vbmEarthSpan(&first, &last);
    const double lo = std::max(l.t0 - 3.0 * l.tE, first), hi = std::min(l.t0 + 3.0 * l.tE, last);
    const int    N  = 600;
    const double h  = (hi - lo) / N;
    // The grid, as one light curve.
    std::vector<double>     tg(N + 1);
    std::vector<ModelPoint> pg(N + 1);
    tg[0] = lo;
    for (int k = 1; k <= N; ++k) tg[k] = lo + k * h;
    evaluateModel(s, l, as, 0, tg.data(), N + 1, false, pg.data());
    int best = 0; double ubest = pg[0].u;
    for (int k = 1; k <= N; ++k) { const double u = pg[k].u; if (u < ubest) { ubest = u; best = k; } }
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
