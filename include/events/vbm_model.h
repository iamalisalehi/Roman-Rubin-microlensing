// The light-curve model's interface to the vendored VBMicrolensing library (external/VBMicrolensing/):
// the point-lens magnification, of a point source or of a uniform disc, along the source-lens track with the
// parallax of the observer, and the astrometric centroid.
//
// This module is the only code that includes the library; its types stay out of this header. It owns one
// VBMicrolensing object, set up once (vbmInit), pointed at one sightline at a time (vbmSetSightline). The
// simulator reaches the model through evaluateModel and evaluateTrack (events/lightcurve.h), which call
// vbmEvaluate once per observer and light curve.
//
// What the library does with the sky, as the rest of the code must read it:
//   - Times are JD - 2450000. A simulation day t (days from the first Rubin bulge visit) is
//     t + SIM_DAY0_JD2450000 on that clock.
//   - Parallax light curves are computed in the GEOCENTRIC frame at t0_par = t0 (Gould 2004): the observer
//     moves on the Earth's orbit minus its straight-line motion at t0, so at t0 the Earth is at rest. The
//     parallax vector piE and the source-lens proper motion are then geocentric, and the heliocentric
//     proper motions the Galactic model draws differ from them by the Earth's projected velocity at t0
//     times the relative parallax (earthVelocityNE gives the velocity).
//   - The sky axes are celestial: piE = (piE_N, piE_E) on (North, East), and the library works with the
//     Earth's offset on (South, West), the opposite directions.
//   - Roman is "satellite 1" (the Earth is satellite 0): its geocentric position, a daily table, is added
//     to the observer's offset, which is how the library does space parallax.
//   - The source is a point (PSPL) or a uniform disc of radius rho Einstein radii (ESPL, finite source with a
//     point lens; no limb darkening). The disc's magnification and centroid come from the library's table of
//     the uniform disc near the lens, and from a quadrupole expansion in rho far from it (the library
//     switches where the correction falls below Tol/10, Tol = VBM_TOL). The table covers rho from 1e-4 to
//     100: below 1e-4 it reads the 1e-4 row and is as accurate (1.5e-5 at worst against the exact value,
//     for rho 1e-6 to 5e-5).
#ifndef ROMAN_EVENTS_VBM_MODEL_H
#define ROMAN_EVENTS_VBM_MODEL_H

#include "common.h"

// One sightline's sky frame.
struct SkyFrame {
    double lDeg, bDeg;       // Galactic longitude and latitude [deg]
    double raDeg, decDeg;    // the same direction in ICRS [deg]
    // The 2x2 orthogonal matrix taking a vector's components on (e_l, e_b), the unit vectors toward increasing
    // Galactic longitude and latitude at (l, b), to its components on (North, East) there:
    //     v_N = toNE[0][0] v_l + toNE[0][1] v_b,    v_E = toNE[1][0] v_l + toNE[1][1] v_b.
    // Orthogonal with determinant -1: (e_l, e_b) is east-then-north on the sky, like (E, N), so against
    // (N, E) it is a rotation by the angle between Galactic north and celestial north at the sightline
    // followed by the swap of the axes.
    double toNE[2][2];
};

// Loads the library's tables (finite-source magnification, the Sun's geocentric ephemeris) and Roman's
// ephemeris, and sets the accuracy and the time convention. Idempotent. Throws if a file is missing or does
// not cover the simulation, naming it; run from the repository root, like everything that reads data.
void vbmInit();

// Points the library at the sightline (l, b) [deg] and returns its sky frame. Calls vbmInit if needed.
// The frame stays available to the rest of the code through vbmSkyFrame until the next call.
SkyFrame vbmSetSightline(double lDeg, double bDeg);

// The current sightline's sky frame; throws if no sightline has been set.
const SkyFrame& vbmSkyFrame();

// The simulation days over which the library's ephemerides can be evaluated for the Earth: the span of the
// Sun table, less the rows the library reads beyond the requested time. A light curve asked for outside it is
// extrapolated silently by the library, so callers that search over a long window clamp to this span.
void vbmEarthSpan(double* firstDay, double* lastDay);

// The Earth's velocity projected on the sky of the current sightline, at simulation day tSimDay, on
// (North, East) [AU/day]: heliocentric, from the library's Sun table. The velocity v_E,perp in
// mu_geo = mu_hel - v_E,perp * pi_rel (pi_rel in rad, mu in rad/day).
void earthVelocityNE(double tSimDay, double v[2]);

// The parallax vector piE_vec on (North, East) of an event whose geocentric source-minus-lens relative motion
// points along xi on (e_l, e_b) [rad]. By the convention of the library (Gould 2004) piE_vec is parallel to
// the LENS-minus-source relative proper motion, the opposite of the direction the source travels, so it is
// piE * toNE * (-cos xi, -sin xi).
void parallaxVectorNE(double piE, double xi, double v[2]);

// An event's relative motion in the geocentric frame at t0, from the heliocentric proper motions the Galactic
// model draws. The observer's own velocity shifts every proper motion by v_E,perp * parallax, so the
// lens-source relative proper motion, seen from an Earth at rest in the frame of the Earth at t0, is
//     mu_geo = mu_hel - v_E,perp(t0) * pi_rel          (lens minus source; mas/day with pi_rel in mas)
// and the event lasts tE = thetaE / |mu_geo| instead of thetaE / |mu_hel|. xi is the direction of the
// geocentric SOURCE-minus-lens motion on (e_l, e_b), the straight-line direction of the light curve.
struct GeocentricMotion {
    double tE;     // Einstein crossing time [days]
    double xi;     // [rad], in [0, 2 pi)
    double muGeo;  // |mu_geo| [mas/day]
};
// t0 on the simulation clock [days]; tetE [mas]; pirel = 1/Dl - 1/Ds [mas]; mus*, mul* the heliocentric proper
// motions of the source and the lens on (e_l, e_b) [mas/day]. Needs a sightline set.
GeocentricMotion geocentricMotion(double t0, double tetE, double pirel, double mus1, double mus2,
                                  double mul1, double mul2);

// The lens's own motion on the sky, from the parameters a light curve and the astrometry measure. The light curve
// gives the geocentric lens-source relative proper motion: its size is thetaE / tE, its direction the parallax
// vector's, i.e. xi (parallaxVectorNE: piE_vec is along lens minus source, so the source moves along +xi).
// geocentricMotion's inverse, then, with the source's own heliocentric motion added:
//     mu_rel,geo (lens minus source, e_l e_b) = -(thetaE / tE) (cos xi, sin xi)
//     mu_rel,hel = mu_rel,geo + v_E,perp(t0) * pi_rel,   pi_rel = thetaE piE     (N, E; Gould 2004)
//     mu_L,hel   = mu_S + mu_rel,hel                                           (mu_S = (mus1, mus2) from astrometry)
//     D_L = 1 / (pi_rel + pi_S),   v_L = D_L mu_L,hel      (pi_S the source parallax, taken as known)
// At the true parameters it returns the draw's (mul1, mul2) and v_L = mul * Dl * AU / (1000 * 86400): the identity
// that pins every sign and frame, checked on every characterised event (ErrorCal's lmCheck).
struct LensMotion {
    double mu[2];  // mu_L,hel on (e_l, e_b) [mas/day]
    double v[2];   // v_L = D_L mu_L,hel on (e_l, e_b) [km/s], heliocentric
    double DL;     // [kpc]
};
// What the function takes from the event and the sightline beyond the fitted parameters: Earth's velocity on
// (N, E) at t0, the frame matrix, the source parallax.
struct MotionContext {
    double vE[2];       // Earth's velocity on (North, East) at t0 [AU/day]
    double toNE[2][2];  // the sightline's SkyFrame::toNE
    double piS;         // source parallax [mas]
};
// t0 on the simulation clock [days]; piS [mas]. Needs a sightline set.
MotionContext motionContext(double t0, double piS);
// tE [days] and xi [rad] geocentric, piE [], tetE [mas], mus1/mus2 heliocentric on (e_l, e_b) [mas/day].
LensMotion lensMotion(const MotionContext& c, double tE, double xi, double piE, double tetE, double mus1, double mus2);

// The source of the light-curve model: a point (PSPL) or a uniform disc (ESPL).
enum class SourceModel { PSPL, ESPL };

// The source model of the run, which evaluateModel and evaluateTrack read (events/lightcurve.cpp) when they
// build the ModelEvent of the simulator's source and lens. Set once, before the first light curve; the
// default is ESPL.
void        vbmSetSourceModel(SourceModel m);
SourceModel vbmSourceModel();

// The model of one event, in the parameters the simulator carries.
struct ModelEvent {
    double u0;           // impact parameter at t0, > 0 [Einstein radii]
    double tE;           // Einstein crossing time, geocentric [days]
    double t0;           // time of closest approach of the geocentric straight line [simulation days]
    double piE;          // microlensing parallax |piE_vec| = pi_rel / thetaE []
    double xi;           // direction of the geocentric source-minus-lens relative motion on (e_l, e_b) [rad]
    double tetE;         // angular Einstein radius [mas]
    double mus1, mus2;   // source proper motion, HELIOCENTRIC, on (e_l, e_b) [mas/day]
    double piS;          // source parallax [mas]
    double rho;          // source radius in Einstein radii, theta* / thetaE >= 0; read only by ESPL
    SourceModel source;  // ESPL with rho = 0 is the point source
};

// The model at one epoch, on (e_l, e_b).
struct ModelEpoch {
    double A;                // magnification of the model: magnifOf(u) for a point source (bit for bit)
    double u;                // |source - lens| [Einstein radii]
    double ux, uy;           // source - lens, WITH the parallax of the observer [Einstein radii]
    // Positions [mas] on the sky from an origin at the lens at t0; set only when astrometry is asked for.
    double imgx, imgy;       // centroid of the lensed images
    double lensx, lensy;     // the lens
    double srcx, srcy;       // the source as it would be unlensed: lens + thetaE (ux, uy)
};

// The event seen by one observer, the Earth (atL2 = false) or Roman at L2, at the n epochs t[0..n-1]
// [simulation days]: out[i] is the model at t[i]. One call of the library for the whole array. Without
// astrometry the library computes the photometric track only (the positions are not set). Throws if the
// library had to extrapolate an ephemeris, i.e. an epoch is outside vbmEarthSpan (or, for Roman, outside
// the span of its table).
//
// u0 > 0 and the line u0 (-sin xi, cos xi) + tau (cos xi, sin xi) are the legacy straight-line geometry
// of the source relative to the lens on (e_l, e_b); only the parallax is the library's.
//
// The astrometric call perturbs the parallax vector at every epoch (see vbm_model.cpp), which the
// photometric one does not: the two tracks differ by up to ~2e-10 x (epoch index) Einstein radii, 2e-7 at the
// end of a 1024-epoch array (and A by that over u).
void vbmEvaluate(const ModelEvent& ev, bool atL2, bool astrometry, const double* t, int n, ModelEpoch* out);

// The magnification of the event's source at the impact parameter u > 0 [Einstein radii], as vbmEvaluate
// computes it along the track, but for a given u and so without any parallax: the straight line the
// light curve is parametrised by. magnifOf(u) for a point source, the library's uniform-disc ESPLMag2 for a
// finite one.
double vbmMagnification(const ModelEvent& ev, double u);

// TEST HELPER, not used by the simulation: how far the track of the source relative to the lens, in
// Einstein radii, moves when the observer is Roman at L2 instead of the Earth, at simulation day tSimDay,
// for the parallax vector (piEN, piEE). dy = (y1 - y1_Earth, y2 - y2_Earth), y1 along the motion and y2
// across it, as the library defines them. Returns the library's extrapolation flag: 0 when both
// ephemerides cover tSimDay, 1 for the Sun table, 2 for Roman's table.
int vbmTestSatelliteShift(double tSimDay, double piEN, double piEE, double dy[2]);

#endif // ROMAN_EVENTS_VBM_MODEL_H
