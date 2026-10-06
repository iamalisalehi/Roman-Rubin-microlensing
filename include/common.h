// What nearly every file needs: system and GSL includes, the project configuration headers,
// the CHECK/MIN/Nx/Ny macros, the selected lens population, and the structural constants.
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

/// for reading the extinction maps
namespace fs = std::filesystem;
///

using std::cout;
using std::endl;
using std::cin;

// Photometric Fisher parameters, in index order:
//   0 u0   1 tE   2 fb0   3 piE   4 xi   5 t0   6 mbs0   7 fb1   8 mbs1
// fb0/mbs0 are Rubin's source-flux fraction and baseline magnitude; fb1/mbs1 are Roman's.
// (fb is the fraction of aperture flux coming from the SOURCE, despite the name -- see
// Lensing.cpp. Together with the baseline magnitude it is a bijective reparametrization of
// the source-flux / blend-flux pair: F_src = fb * 10^(-0.4 mbs), F_bl = (1-fb) * 10^(-0.4 mbs).)
// t0, mbs0, fb1 and mbs1 were appended rather than inserted so that indices 0-4 keep the
// meanings hard-coded throughout co.resu[]. See DEVIATIONS.md.
#define Nx 9
#define Ny 4
#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define CHECK(cond) \
    do { \
        if (!(cond)) \
            throw std::runtime_error("Check failed: " #cond); \
    } while (0)


// Lens-population table (POPULATIONS), MassFunction and LensPopulation: see config/parameters.h,
// section 6, and include/parameter_types.h.

// Set once, from --population, before any `lens` is constructed -- the mass-efficiency grid
// is built in that constructor from the bounds below.
inline const LensPopulation* gPop = &POPULATIONS[0];

inline double mlMin() { return gPop->mlMin; }
inline double mlMax() { return gPop->mlMax; }

// Physical constants: include/physical_constants.h. Model/survey parameters: config/parameters.h.
// Data-file row counts and catalogue mean masses: config/data_products.h.


///=============================LSST constant===============================///

//constexpr double cade2 = 10.0;//ELT [days]

//constexpr int YZ = 3578;   //No.yzma.txt rows
//constexpr int met = 70;   //No rows metal.txt
//constexpr int Nel = 7;//rows in "sigma_ELT.txt" (https://academic.oup.com/mnras/article/494/3/4413/5813442)
//constexpr int nfiles = 2518;
//constexpr int nlines = 3686;
//constexpr int nex = 2518 * 3686; //number of ext files * lines in each file
// NFILES/NROWS (the 2,518 x 3,686 per-pointing extinction tables) are gone: the single table file
// carries its own sizes (Deviation 70; struct extin).
constexpr int nq = 15;     //resu


// `coun` bounds the length of one event's light-curve buffers (lens::timn/magn/errm/
// soux/souy/erra/tele). `ndw`, the running count of accepted epochs (Rubin + Roman
// combined) for a single star, can never exceed the total number of visit rows that
// exist for either instrument — Nl + NlRoman is therefore an unconditionally safe
// upper bound at any sky position, in any season. One-time memory cost since `lens`
// is allocated once, not per event: (Nl+NlRoman) * 7 buffers * 8 bytes ~= 17.7 MB.
constexpr int    coun  = Nl + NlRoman;

constexpr double tetp   = double(M_PI / 3.0);        //parallax


constexpr double ROMAN_AST_SLOPE_BKG = 0.4;


// Finite-difference stencils used by FisherM. For each parameter it evaluates the model at
// theta + Delta*s[h] for h = 0,1, forms (model - stored)/(Delta*s[h]), and averages the two.
//
// sig gives a CENTRAL difference: the h=0 and h=1 terms are
//     (f(x+D) - f(x))/D   and   (f(x-D) - f(x))/(-D),
// whose average is (f(x+D) - f(x-D))/(2D). The O(D) error terms cancel, leaving O(D^2).
//
// sig2 averages two FORWARD differences, at D/2 and D. Nothing cancels: the result is
// f'(x) + (3/8) f''(x) D + O(D^2) -- first order, and biased. On a smooth test function it
// carries ~22x the error of sig at equal step and gains only one decade of accuracy per decade
// of step reduction, against sig's two (Step C3).
//
// sig2 was applied to exactly the strictly-positive parameters -- tE, piE (photometric) and
// tetE, piE (astrometric) -- which suggests the legacy intent was to avoid ever stepping them
// through zero. That is not a real constraint here: the steps are a small fraction of the
// parameter, so theta - Delta stays positive, and after Step C3 they are smaller still.
// The photometric tE and piE now use sig. The two astrometric uses are unchanged pending a
// step-size sweep of Delta2 -- see OPEN_ITEMS.md.
constexpr std::array<double, 2> sig  = {+1.0 ,-1.0};
constexpr std::array<double, 2> sig2 = {+0.5 ,+1.0};


// Lens mass range: now a property of the selected population (see POPULATIONS above), not a
// constant, because a bulge population and a black-hole population share no sensible bounds.
// The same numbers set the Mls grid the mass-efficiency histogram is binned on, so they must
// bracket the masses actually drawn or that output collapses into one bin.
//
// Kept as names because the rest of the code reads them as names; they are now functions.
#define Ml_min (mlMin())
#define Ml_max (mlMax())

// Bulge distance grid (Num, MaxD, step, dd): see config/parameters.h, section 2.
////=================================== Bulge ====================================
//const double RaLMC  =  80.89375;
//const double DecLMC = -68.2438888888889;
//const double DLMC =  49.97;///KPC


///============================================================================
enum class GalacticComponent : int {
    THIN_DISK    = 0,
    BULGE        = 1,
    THICK_DISK   = 2,
    HALO         = 3,
};

#endif // ROMAN_COMMON_H
