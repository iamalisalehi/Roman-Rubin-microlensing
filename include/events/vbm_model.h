// The light-curve model's interface to the vendored VBMicrolensing library (external/VBMicrolensing/).
//
// This module is the only code that includes the library; its types stay out of this header. It owns one
// VBMicrolensing object, set up once (vbmInit), pointed at one sightline at a time (vbmSetSightline).
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
SkyFrame vbmSetSightline(double lDeg, double bDeg);

// The Earth's velocity projected on the sky of the current sightline, at simulation day tSimDay, on
// (North, East) [AU/day]: heliocentric, from the library's Sun table. The velocity v_E,perp in
// mu_geo = mu_hel - v_E,perp * pi_rel (pi_rel in rad, mu in rad/day).
void earthVelocityNE(double tSimDay, double v[2]);

// TEST HELPER, not used by the simulation: how far the track of the source relative to the lens, in
// Einstein radii, moves when the observer is Roman at L2 instead of the Earth, at simulation day tSimDay,
// for the parallax vector (piEN, piEE). dy = (y1 - y1_Earth, y2 - y2_Earth), y1 along the motion and y2
// across it, as the library defines them. Returns the library's extrapolation flag: 0 when both
// ephemerides cover tSimDay, 1 for the Sun table, 2 for Roman's table.
int vbmTestSatelliteShift(double tSimDay, double piEN, double piEE, double dy[2]);

#endif // ROMAN_EVENTS_VBM_MODEL_H
