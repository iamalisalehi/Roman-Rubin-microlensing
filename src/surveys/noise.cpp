// Photometric and astrometric error models of the two surveys (surveys/noise.h).
#include "surveys/noise.h"

///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Error LSST calculations                    //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
double errlsstM(double mag, int fi, double sig5){ //LSST Photometric Error 

    double x, Delta1 = 0.0;
    x = std::pow(10.0, 0.4 * (mag - sig5));
    Delta1 = std::sqrt(std::fabs((LSST_ERR_C04 - gama[fi]) * x + gama[fi] * x * x));
    
    if (Delta1 < 0.0001)   Delta1 = 0.0001;

    CHECK(sig5 >= 0.0);
    CHECK(sig5 <= 40.0);
    CHECK(Delta1 >= 0.00001);
    CHECK(x > 0.0);
   
    return std::sqrt(delta2 * delta2 + Delta1 * Delta1);
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double errlsstA(lsst & ls, double ghadr){ //LSST Astrometric Error  //Change it!!

    double error = -1.0, shib = 0.0;

    if (ghadr < ls.mag[0] or  ghadr == ls.mag[0]) error = double(ls.err[0]);

    else if (ghadr > ls.mag[Na-1] or ghadr == ls.mag[Na-1]) {
        shib = (ls.err[Na-1] - ls.err[Na-2]) / (ls.mag[Na-1] - ls.mag[Na-2]);
        error = double(ls.err[Na-1] + shib * (ghadr - ls.mag[Na-1]));
    }

    else {
        for (int i = 1; i < Na; ++i) {
            if (double((ghadr - ls.mag[i]) * (ghadr - ls.mag[i-1])) < 0.0 or ghadr == ls.mag[i-1]) {
                shib = (ls.err[i] - ls.err[i-1]) / (ls.mag[i] - ls.mag[i-1]);
                error = double(ls.err[i-1] + shib * (ghadr - ls.mag[i-1]));
                break;
            }
        }
    }

    CHECK(error > 0.0);
    CHECK(error >= ls.err[0]);
    CHECK(ghadr >= 0.0);

    // Renormalise the shipped mission-averaged curve to a PER-VISIT error, which is what
    // l.erra[] means and what FisherM assumes. See the LSST_AST_* block in config/parameters.h for the
    // two independent checks that fix the factor at 26.74. Applied here rather than by
    // editing files/sigmaA_LSST.txt so the input data stay as delivered and the correction
    // is visible in the code that depends on it.
    error *= LSST_AST_RENORM;

    CHECK(error >= LSST_AST_FLOOR * 0.999);
    CHECK(std::isfinite(error));

    return(error);
}


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
// Roman WFI astrometric error for one F146 exposure, in milliarcseconds (Step H4).
//
// Replaces the errlsstA() placeholder that stood in for Roman -- Rubin's astrometric error
// curve evaluated at Roman's magnitude, which had no reason to be right and was flagged in
// OPEN_ITEMS.md. Constants, their sources and the per-exposure caveat are in config/parameters.h.
//
// Three regimes:
//   m <= 20.62   1.1 mas       centroiding floor, 1% of the 110 mas pixel. A systematic,
//                              not photon noise, so it does NOT improve for brighter stars.
//   20.62 - 23.5 rises at 0.333/mag   interpolation between the two published anchors.
//   m >  23.5    rises at 0.4/mag     background dominated, SNR ~ counts.
//
// Unlike errlsstA this reads no data file and needs no instrument struct, so it takes the
// magnitude alone. Its photometric sibling errRomanM(), defined below, takes the roman struct,
// because it needs the lookup table.
double errRomanA(double magF146AB){

    // The anchors below are VEGA magnitudes; the simulator's are AB (Deviation 72).
    const double magF146 = magF146AB - F146_AB_MINUS_VEGA;
    double error = ROMAN_AST_FLOOR;

    if (magF146 > ROMAN_AST_MBKG) {
        error = ROMAN_AST_SBKG * std::pow(10.0, ROMAN_AST_SLOPE_BKG * (magF146 - ROMAN_AST_MBKG));
    }
    else if (magF146 > ROMAN_AST_MFLR) {
        error = ROMAN_AST_FLOOR * std::pow(10.0, ROMAN_AST_SLOPE_SRC * (magF146 - ROMAN_AST_MFLR));
    }

    // The floor is a floor: the source-dominated branch is continuous with it at
    // ROMAN_AST_MFLR by construction, but clamp anyway so no future edit to the constants can
    // silently return a precision better than Roman can centroid.
    if (error < ROMAN_AST_FLOOR) error = ROMAN_AST_FLOOR;

    CHECK(error >= ROMAN_AST_FLOOR);
    CHECK(std::isfinite(error));

    return(error);
}

// Roman's per-exposure F146 photometric error [mag] at AB magnitude `mag` (Deviation 72): the
// Penny et al. 2019 curve, anchored at load time to the 66-s 5-sigma depth (see config/parameters.h), with
// log(err) interpolated linearly in magnitude. Brighter than the table: its first value (the 1 mmag
// floor dominates there); fainter: extrapolated along the last segment (the recording gate stops at
// the 5-sigma depth, well inside the table, so this branch only serves diagnostics).
double errRomanM(const roman& ro, double mag)
{
    if (mag <= ro.mag[0]) return ro.err[0];
    int i = 1;
    while (i < NaRoman - 1 and ro.mag[i] < mag) ++i;
    const double f = (mag - ro.mag[i - 1]) / (ro.mag[i] - ro.mag[i - 1]);
    return std::exp(std::log(ro.err[i - 1]) + f * (std::log(ro.err[i]) - std::log(ro.err[i - 1])));
}
