// The seeded random-number generator (one mt19937_64, re-seeded per sightline) and the scalar draws built on it.
#ifndef ROMAN_UTIL_RANDOM_H
#define ROMAN_UTIL_RANDOM_H

#include "common.h"

////
#include <random>
inline std::mt19937_64 rng{seed};

// Deviation 77: the generator is RE-SEEDED at the start of every sightline from (base seed,
// sightline index), so a sightline's draws do not depend on which sightlines ran before it. A scan
// split into chunks (--start-index / --end-index), or resumed, reproduces the unsplit run exactly.
// SplitMix64 (Steele, Lea & Flood 2014) mixes the pair into well-separated 64-bit seeds.
inline std::uint64_t splitmix64(std::uint64_t x)
{
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}
inline std::uint64_t sightlineSeed(std::uint64_t base, long index)
{
    return splitmix64(splitmix64(base) ^ static_cast<std::uint64_t>(index));
}
//inline std::mt19937_64 rng{std::random_device{}()};

double RandN(double , double);
double RandR(double , double);
int    RandPois(double);

#endif // ROMAN_UTIL_RANDOM_H
