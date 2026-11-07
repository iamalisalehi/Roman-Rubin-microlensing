// The seeded random-number generator (one mt19937_64, re-seeded per sightline) and the scalar draws built on it.
#ifndef ROMAN_UTIL_RANDOM_H
#define ROMAN_UTIL_RANDOM_H

#include "common.h"

#include <random>
inline std::mt19937_64 rng{seed};

// The generator is re-seeded at the start of every sightline from (base seed, sightline index), so
// a sightline's draws do not depend on which ran before it; a scan split into chunks
// (--start-index / --end-index), or resumed, reproduces the unsplit run exactly.
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

// A second generator for the positions of unresolved neighbours (func_source), re-seeded per sightline
// like `rng` but from a salted base, so drawing them leaves every draw of the main stream unchanged.
inline std::mt19937_64 rngBlend{seed};
constexpr std::uint64_t BLEND_STREAM_SALT = 0xB7E151628AED2A6BULL;

double RandN(double , double);
double RandR(double , double);
int    RandPois(double);
double RandBlendUnit();   // uniform on [0, 1), from rngBlend

#endif // ROMAN_UTIL_RANDOM_H
