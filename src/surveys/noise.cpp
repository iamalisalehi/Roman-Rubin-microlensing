// Photometric and astrometric error models of the two surveys (surveys/noise.h).
#include "surveys/noise.h"

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

// Per-epoch astrometric error, per coordinate [mas], from the photometric error of the SAME epoch:
// the photon-noise part of that error gives the SNR, and the position error of a PSF fit is
// sigma = kappa * FWHM / SNR (Lindegren 1978), plus the instrument's floor in quadrature. One noise
// model, two outputs. Constants and sources in config/parameters.h.
static double astromFromPhot(double errPhotMag, double photFloorMag, double kappaFwhmMas, double floorMas)
{
    // Photon-noise part of the magnitude error, and the flux SNR it implies (sigma_m = 1.0857 / SNR).
    const double rand2 = errPhotMag * errPhotMag - photFloorMag * photFloorMag;
    const double stat  = (rand2 > 0.0) ? kappaFwhmMas * std::sqrt(rand2) / 1.0857 : 0.0;
    return std::sqrt(stat * stat + floorMas * floorMas);
}

// Rubin, one visit: errPhotMag = errlsstM of this visit (its band, its depth); fwhmArcsec = this
// visit's geometric PSF FWHM.
double errlsstA(double errPhotMag, double fwhmArcsec)
{
    CHECK(errPhotMag > 0.0);
    CHECK(fwhmArcsec > 0.0);
    const double error = astromFromPhot(errPhotMag, delta2, LSST_AST_KAPPA * fwhmArcsec * ARCSEC_TO_MAS,
                                        LSST_AST_FLOOR);
    CHECK(error >= LSST_AST_FLOOR);
    CHECK(std::isfinite(error));
    return error;
}

// Roman, one F146 exposure: errPhotMag = errRomanM of this exposure.
double errRomanA(double errPhotMag)
{
    CHECK(errPhotMag > 0.0);
    const double error = astromFromPhot(errPhotMag, ROMAN_PHOT_FLOOR, ROMAN_AST_K, ROMAN_AST_FLOOR);
    CHECK(error >= ROMAN_AST_FLOOR);
    CHECK(std::isfinite(error));
    return error;
}

// Roman's per-exposure F146 photometric error [mag] at AB magnitude `mag` the
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
