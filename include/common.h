// What nearly every file needs: system and GSL includes, the project configuration headers,
// the Nx/Ny constants, the CHECK/MIN macros, the selected lens population, and the structural constants.
#ifndef ROMAN_COMMON_H
#define ROMAN_COMMON_H

#include <stdlib.h>
#include <stdio.h>
#include <dirent.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <sys/timeb.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <string>
#include <array>
#include <vector>
#include <memory>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <utility>
#include <map>
#include <sstream>
#include <random>
#include <stdexcept>
#include <cstdint>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

// Project configuration, in dependency order: physical constants, parameter types, hand-edited
// parameters (config/parameters.h), and data-file descriptions (config/data_products.h).
#include "physical_constants.h"
#include "parameter_types.h"
#include "parameters.h"
#include "data_products.h"

namespace fs = std::filesystem;

using std::cout;
using std::endl;
using std::cin;

// Photometric Fisher parameters, in index order:
//   0 u0   1 tE   2 fb0   3 piE   4 xi   5 t0   6 mbs0   7 fb1   8 mbs1
// fb0/mbs0 are Rubin's source-flux fraction and baseline magnitude; fb1/mbs1 are Roman's.
// fb is the fraction of aperture flux coming from the SOURCE, despite the name (see
// src/events/source.cpp): F_src = fb * 10^(-0.4 mbs), F_bl = (1-fb) * 10^(-0.4 mbs).
// Indices 0-4 keep the meanings hard-coded throughout co.resu[].
constexpr int Nx = 9;
constexpr int Ny = 4;
#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define CHECK(cond) \
    do { \
        if (!(cond)) \
            throw std::runtime_error("Check failed: " #cond); \
    } while (0)


// Set once, from --population, before any `lens` is constructed -- the mass-efficiency grid
// is built in that constructor from the bounds below.
inline const LensPopulation* gPop = &POPULATIONS[0];

inline double mlMin() { return gPop->mlMin; }
inline double mlMax() { return gPop->mlMax; }

constexpr int nq = 15;     //resu

// `coun` bounds the length of one event's light-curve buffers (lens::timn/magn/errm/
// soux/souy/erra/tele). The running count of accepted epochs (Rubin + Roman combined) cannot
// exceed the total number of visit rows of either instrument, so Nl + NlRoman is a safe bound.
// `lens` is allocated once, so the cost is fixed: (Nl+NlRoman) * 7 buffers * 8 bytes ~= 17.7 MB.
constexpr int    coun  = Nl + NlRoman;

constexpr double tetp   = double(M_PI / 3.0);        //parallax

// Finite-difference stencils used by FisherM. For each parameter it evaluates the model at
// theta + Delta*s[h] for h = 0,1, forms (model - stored)/(Delta*s[h]), and averages the two.
//
// sig is a CENTRAL difference: the two terms (f(x+D) - f(x))/D and (f(x-D) - f(x))/(-D) average
// to (f(x+D) - f(x-D))/(2D), error O(D^2).
//
// sig2 averages two FORWARD differences, at D/2 and D: f'(x) + (3/8) f''(x) D + O(D^2), i.e. first
// order and biased (~22x the error of sig at equal step on a smooth test function). It is used
// only for the astrometric tetE and piE.
constexpr std::array<double, 2> sig  = {+1.0 ,-1.0};
constexpr std::array<double, 2> sig2 = {+0.5 ,+1.0};

enum class GalacticComponent : int {
    THIN_DISK    = 0,
    BULGE        = 1,
    THICK_DISK   = 2,
    HALO         = 3,
};

#endif // ROMAN_COMMON_H
