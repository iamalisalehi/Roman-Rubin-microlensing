#ifndef PHYSICAL_CONSTANTS_H
#define PHYSICAL_CONSTANTS_H

// Physical constants and unit conversions. These never change with the model, the survey or the
// data -- do not edit them to explore a scenario (see config/parameters.h for that).

#include <cmath>

constexpr double RAa = 180.0 / M_PI;
constexpr double pi  = M_PI;
constexpr double velocity = 299792458.0;//velosity of light
constexpr double Msun = 1.98892 * std::pow(10., 30); //in [kg].
constexpr double Rsun = 6.957 * std::pow(10.0, 8.0); ///solar radius [meter]
constexpr double KP = 3.08568025 * std::pow(10., 19); // in meter.
constexpr double G = 6.67384 * std::pow(10., -11.0);// in [m^3/s^2*kg].
constexpr double AU = 1.4960 * std::pow(10.0, 11.0);
constexpr double year = 365.2425;//days
constexpr double eps = double(0.000000000000005463263454624313654);
constexpr double ARCSEC_TO_MAS    = 1000.0;
constexpr double AU_KM        = 1.496e8;
constexpr double omegae = double(2.0 * M_PI / year); //radian per day
constexpr double vearth = omegae;                    //radian per day

#endif // PHYSICAL_CONSTANTS_H
