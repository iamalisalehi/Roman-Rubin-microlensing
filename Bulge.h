#ifndef LMC_H
#define LMC_H

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
////
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


struct ScaRect { double l0, l1, b0, b1; };   // offsets from the field centre [deg], l0<l1, b0<b1

struct GbtdsLayout {
    std::array<std::vector<ScaRect>, GBTDS_NLAYOUT> sca;
    // Bounding rectangle of each layout's detectors, for a cheap early reject and for the fine
    // stratum of the sightline grid.
    std::array<double, GBTDS_NLAYOUT> dlMin{}, dlMax{}, dbMin{}, dbMax{};
    double rField  = 0.0;   // largest centre-to-detector-corner distance, either layout [deg]
    double scaSide = 0.0;   // smallest detector side [deg]
};

// A field placement: one (centre, roll) that RomanBaseline.dat visits. 6 fields x 2 rolls = 12.
struct FieldPlacement { double l, b; int layout; };

constexpr double tetp   = double(M_PI / 3.0);        //parallax


constexpr double ROMAN_AST_SLOPE_BKG = 0.4;


// The two images of one source at impact parameter u, with their own (unblended) magnitudes.
// `magBase`/`blendFrac` are the BLENDED baseline magnitude and the source's flux share in the
// filter being observed, i.e. s.magb[i] and s.blend[i] -- per FILTER, not per telescope.
struct ImagePair {
    double sep;            //angular separation of the two images [mas]; -1 if undefined
    double magPlus;        //apparent magnitude of the major image [mag]
    double magMinus;       //apparent magnitude of the minor image [mag]
    bool   bothDetectable; //both within [saturation, single-visit depth] in this filter
};

inline ImagePair imagePair(double u, double tetE, double magBase, double blendFrac,
                           double thrMag, double satMag)
{
    ImagePair ip{-1.0, 99.0, 99.0, false};
    if (!(u > 0.0) or !(tetE > 0.0) or !(blendFrac > 0.0)) return ip;

    const double root = std::sqrt(u * u + 4.0);
    ip.sep = tetE * root;                                   //[mas]

    const double Aplus  = (u * u + 2.0) / (2.0 * u * root) + 0.5;
    const double Aminus = (u * u + 2.0) / (2.0 * u * root) - 0.5;
    if (!(Aminus > 0.0) or !(Aplus > 0.0)) return ip;       //minor image formally extinguished

    // The SOURCE's own magnitude, recovered from the blended baseline: blendFrac is the
    // source's share of the aperture flux, so m_source = m_base - 2.5 log10(blendFrac) and is
    // always the FAINTER of the two (blendFrac <= 1). Each image then carries its own
    // magnification. The blend light is deliberately NOT added back: an image that is resolved
    // from its twin is resolved from the neighbours too, and re-adding the full blend would
    // make a faint minor image look detectable on light that is not its own.
    const double magSource = magBase - 2.5 * std::log10(blendFrac);
    ip.magPlus  = magSource - 2.5 * std::log10(Aplus);
    ip.magMinus = magSource - 2.5 * std::log10(Aminus);

    ip.bothDetectable = (ip.magPlus  <= thrMag and ip.magPlus  >= satMag and
                         ip.magMinus <= thrMag and ip.magMinus >= satMag);
    return ip;
}

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
struct GSLMatrixDeleter {
        void operator()(gsl_matrix* m) const {
        gsl_matrix_free(m);
    }
};

using gsl_matrix_uptr = std::unique_ptr<gsl_matrix,GSLMatrixDeleter>;

///============================================================================
enum class GalacticComponent : int {
    THIN_DISK    = 0,
    BULGE        = 1,
    THICK_DISK   = 2,
    HALO         = 3,
};

///============================================================================
static std::vector<double> make_grid(double min, double max)
{
    std::vector<double> v(GG + 1);

    for (int i = 0; i <= GG; ++i)
    {
        v[i] = min + (max - min) * static_cast<double>(i) / GG;
    }

    return v;
}

// Geometric spacing, for a quantity whose range spans decades. A linear grid over 3-1000
// Msun would put every black hole below 13 Msun -- most of a log-uniform population -- into
// the first bin and report one number for the entire low-mass end, which is where the
// interesting behaviour is. Only used where the population asks for it.
static std::vector<double> make_grid_log(double min, double max)
{
    std::vector<double> v(GG + 1);
    const double lo = std::log(min), hi = std::log(max);

    for (int i = 0; i <= GG; ++i)
    {
        v[i] = std::exp(lo + (hi - lo) * static_cast<double>(i) / GG);
    }

    return v;
}

///============================================================================
struct source {
    int nums, cl;

    double mass, logT, typ, age, ros;
    double mus1, mus2, mus;
    double xv, yv, zv;
    double Av; //, Avv;
    double SV_n1, LV_n1, VSun_n1;
    double SV_n2, LV_n2, VSun_n2;
    double pos1b, pos2b, pos1c, pos2c, def1a, def2a, def1c, def2c;
    double errM, errA, FWHM;
    double ut, ut0, Astar, xi, ux, uy;
    double Ds, TET, FI, lat, lon, vs;
    double Nstart, Rostart, Romaxs, Romins, nstart, nstarti;
    double od_thin, od_thick, od_bulge, od_halo, opt;

    std::array<double, 2> fb, mbs; // small fixed arrays
    // Share of each telescope's baseline flux that is the LENS's own light (Deviation 74): 0 for
    // dark lenses; for a luminous lens it is already inside the blend (fb counts the source only).
    std::array<double, 2> fLens{};

    std::vector<double> nssim; // size Num
    std::vector<double> nsdet; // size Num
    std::vector<double> rho_thin;
    std::vector<double> rho_thick;
    std::vector<double> rho_halo;
    std::vector<double> rho_stars;
    std::vector<double> rho_bulge;
    std::vector<double> Rostar0;
    std::vector<double> Rostari;
    std::vector<double> Nstari;

    std::vector<double> nsbl;  // size M+1
    std::vector<double> blend; // size M+1
    std::vector<double> Fluxb; // size M+1
    std::vector<double> magb;  // size M+1
    std::vector<double> Ai;    // size M+1
    std::vector<double> Mab;   // size M+1
    std::vector<double> Map;   // size M+1

    GalacticComponent struc;
    // Constructor
    source()
        : nssim(Num), nsdet(Num),
          rho_thin(Num), rho_thick(Num), rho_halo(Num), rho_stars(Num), rho_bulge(Num),
          Rostar0(Num), Rostari(Num), Nstari(Num),
          nsbl(M), blend(M), Fluxb(M), magb(M), Ai(M), Mab(M), Map(M)
    {}
};

struct lens {
    int numl;
    bool luminous = false;          // a living (main-sequence) star, whose light is blended (Dev. 74)

    double Ml, Dl, vl, Vt, xls, u0, A0, mi1, mi2;
    double rhomaxl, tE, RE, t0, murel, DeltaT;
    double piE, pirel, tetE, pos1, pos2, mul1, mul2, mul;
    double betal, betas, deltal, deltas, deltao;

    std::array<double, 2> Nhalo, Nself; // small fixed arrays

    std::vector<double> nstE; // size GG+1
    std::vector<double> ndtE; // size GG+1
    std::vector<double> NstE; // size GG+1
    std::vector<double> NdtE; // size GG+1
    std::vector<double> NsMl; // size GG+1
    std::vector<double> NdMl; // size GG+1
    std::vector<double> Nspi; // size GG+1
    std::vector<double> Ndpi; // size GG+1
    std::vector<double> Nsu0; // size GG+1
    std::vector<double> Ndu0; // size GG+1
    std::vector<double> Nsmb; // size GG+1
    std::vector<double> Ndmb; // size GG+1
    std::vector<double> Nsfb; // size GG+1
    std::vector<double> Ndfb; // size GG+1
    std::vector<double> Nsmu; // size GG+1
    std::vector<double> Ndmu; // size GG+1
    std::vector<double> timn; // size coun
    std::vector<double> magn; // size coun
    std::vector<double> soux; // size coun
    std::vector<double> souy; // size coun
    std::vector<double> errm; // size coun
    std::vector<double> erra; // size coun
    std::vector<int> tele;    // size coun
    std::vector<int> rseas;   // size coun: Roman season index of the epoch (-1 for Rubin), Dev. 71
    std::vector<int> rroll;   // size coun: Roman roll of the epoch, 0 spring / 1 autumn (-1 Rubin)
    
    std::vector<double> tEs;  // size GG+1
    std::vector<double> Mls;  // size GG+1
    std::vector<double> pis;  // size GG+1
    std::vector<double> u0s;  // size GG+1
    std::vector<double> mbs;  // size GG+1
    std::vector<double> fbs;  // size GG+1
    std::vector<double> mus;  // size GG+1

    GalacticComponent struc;
    // Constructor
    lens()
        : nstE(GG+1), ndtE(GG+1), NstE(GG+1), NdtE(GG+1),
          NsMl(GG+1), NdMl(GG+1),  Nspi(GG+1), Ndpi(GG+1),
          Nsu0(GG+1), Ndu0(GG+1),
          Nsmb(GG+1), Ndmb(GG+1),
          Nsfb(GG+1), Ndfb(GG+1),
          Nsmu(GG+1), Ndmu(GG+1),
          timn(coun), magn(coun), soux(coun), souy(coun), errm(coun), erra(coun),
          tele(coun), rseas(coun, -1), rroll(coun, -1),

          tEs(make_grid(tE_min, tE_max)),
          Mls(gPop->logGrid ? make_grid_log(Ml_min, Ml_max) : make_grid(Ml_min, Ml_max)),
          pis(make_grid(pi_min, pi_max)),
          u0s(make_grid(u0_min, u0_max)),
          mbs(make_grid(mb_min, mb_max)),
          fbs(make_grid(fb_min, fb_max)),
          mus(make_grid(mu_min, mu_max))
    {}
};

struct astromet{
   double Ve_n1, Ve_n2;
   double ue_n1, ue_n2;

   // Step H1/H3. Multiplies L2_OFFSET_AU in lightcurve(), so 1 = Roman at L2 (the physical
   // configuration) and 0 = Roman at the centre of the Earth (what the code did before H1).
   // It exists so the satellite-parallax experiment of PHASE_H_PLAN.md Step H3 -- two runs
   // identical but for the spatial baseline -- is a command-line flag rather than a rebuild,
   // and so the H1 regression can prove the term is a clean no-op when switched off.
   double satScale = 1.0;
};

struct CMD {
    // Thin disk
    std::vector<double> logT_thin;
    std::vector<double> mass_thin;
    std::vector<std::array<double, M>> Mab_thin; // M × N1
    std::vector<double> typ_thin;
    std::vector<double> cl_thin;
    std::vector<double> age_thin;

    // Bulge
    std::vector<double> logT_bulge;
    std::vector<double> mass_bulge;
    std::vector<std::array<double, M>> Mab_bulge; // M × N2
    std::vector<double> typ_bulge; 
    std::vector<double> cl_bulge;
    std::vector<double> age_bulge;

    // Thick disk
    std::vector<double> logT_thick;
    std::vector<double> mass_thick;
    std::vector<std::array<double, M>> Mab_thick; // M × N3
    std::vector<double> typ_thick;
    std::vector<double> cl_thick;
    std::vector<double> age_thick;

    // Halo
    std::vector<double> logT_halo;
    std::vector<double> mass_halo;
    std::vector<std::array<double, M>> Mab_halo; // M × N4
    std::vector<double> typ_halo;
    std::vector<double> cl_halo;
    std::vector<double> age_halo;

    // Constructor
    CMD()
        : logT_thin(N1),  mass_thin(N1),  Mab_thin(N1),  typ_thin(N1),  cl_thin(N1),  age_thin(N1),
          logT_bulge(N2), mass_bulge(N2), Mab_bulge(N2), typ_bulge(N2), cl_bulge(N2), age_bulge(N2),
          logT_thick(N3), mass_thick(N3), Mab_thick(N3), typ_thick(N3), cl_thick(N3), age_thick(N3),
          logT_halo(N4),  mass_halo(N4),  Mab_halo(N4),  typ_halo(N4),  cl_halo(N4),  age_halo(N4)
    {}
};

// The extinction tables (Deviation 70): files/ext/ext_tables.dat, built by maps.py from DECaPS and
// Marshall through dustmaps. One shared distance grid; per table position its (l, b) and A_V on that
// grid, row-major in `ext`. Sizes come from the file's header at run time -- no NFILES/NROWS to keep
// in step -- and A_V is stored as float (the maps' own errors are tenths of a magnitude), so the
// ~16,000 positions x 399 distances take ~25 MB, against 148 MB for the 2,518 per-pointing tables
// with their own distance columns before.
struct extin
{
    std::vector<double> dist;    // n_dist, kpc, strictly increasing
    std::vector<double> l, b;    // n_tables
    std::vector<float>  ext;     // n_tables * n_dist, A_V [mag]
    int nTables = 0, nDist = 0;
    std::string built, k;        // provenance, from the header ("# ext_tables ...", "# k ...")
};

struct lsst {
    std::vector<double> mag;   // Na
    std::vector<double> err;   // Na
    std::vector<int> filter;    // Nl
    
    std::vector<int> ct;        // Nl -- one slot per Rubin visit; see matchVisibleEpochs
    std::vector<double> RA;     // Nl
    std::vector<double> DEC;    // Nl
    std::vector<double> l;      // Nl
    std::vector<double> b;      // Nl
    std::vector<double> tim;    // Nl
    std::vector<double> sig5;   // Nl
    std::vector<double> dist;   // Nl
    std::vector<double> rot;    // Nl -- OpSim rotSkyPos [deg] (Deviation 80)

    //ID  RA  Dec  l  b  start  filter  airmass  seeing  skyBrightness visittime sigma5 targetname distance
    lsst()
        : mag(Na), err(Na),
          filter(Nl),
          ct(Nl),
          RA(Nl), DEC(Nl), l(Nl), b(Nl), tim(Nl), sig5(Nl), dist(Nl), rot(Nl)
    {}
};

struct roman {
    std::vector<double> mag;   // NaRoman: mag-vs-error lookup (sigma_roman.txt)
    std::vector<double> err;   // NaRoman

    // --- New: per-visit epoch bookkeeping, mirrors lsst's fields ---
    std::vector<int> ct;        // NlRoman -- visible-epoch indices for the current sightline.
                                // MUST be the full visit count, not a round number: a
                                // sightline inside a GBTDS field matches ~51,500 visits, and
                                // truncating keeps only the earliest, silently ending Roman's
                                // mission 8 days in. (1.24 MB, allocated once.)
    std::vector<double> RA;     // NlRoman
    std::vector<double> DEC;    // NlRoman
    std::vector<double> l;      // NlRoman
    std::vector<double> b;      // NlRoman
    std::vector<double> tim;    // NlRoman
    std::vector<double> sig5;   // NlRoman — only needed if the Roman photometric error
                                // model varies per-visit; otherwise mag/err alone may suffice.
    std::vector<int> field;     // NlRoman -- GBTDS field index, 0-4 contiguous block, 5 GC
    std::vector<int> layout;    // NlRoman -- roll of this visit, 0 spring / 1 autumn; selects
                                // which detector layout is placed at (l, b)

    // NOTE: no `filter` array — currently only F146 (constant filter index 6) is modeled
    // for Roman. If F087/F213 are added later, give roman a `filter` array like lsst's.

    roman()
        : mag(NaRoman), err(NaRoman),
          ct(NlRoman),
          RA(NlRoman), DEC(NlRoman), l(NlRoman), b(NlRoman), tim(NlRoman), sig5(NlRoman),
          field(NlRoman), layout(NlRoman)
    {}
};

// ---------------------------------------------------------------------------------------------
// Which subset of a single event's light curve a Fisher matrix was accumulated from.
//
// The model, and therefore every derivative, is identical across all three -- the same event,
// the same physics. Only the set of epochs summed over differs. Because each epoch contributes
// an independent, positive semi-definite term to the information sum, and the epochs partition
// cleanly by observatory, the matrices satisfy exactly
//
//     F[SJOINT] == F[SRUBIN] + F[SROMAN]
//
// element by element. tests/fisher_fixture.cpp asserts this; if it ever fails, the partitioning
// is wrong. It also follows that sigma from the joint matrix can never exceed sigma from either
// single-survey matrix: adding information can only sharpen a forecast.
enum SurveyIdx { SJOINT = 0, SRUBIN = 1, SROMAN = 2, NSURV = 3 };

// ---------------------------------------------------------------------------
// Roman observing-season geometry (Step D1).
//
// Roman can only look at the bulge when the Sun angle permits, so its ten-year visit
// list is a comb of ~70-day observing seasons separated by ~110-day gaps. Where an
// event's peak falls relative to those seasons is the independent variable of the
// gap-filling result, so it has to be computed at simulation time and stored.
//
// The windows are DERIVED from the epoch times actually present in RomanBaseline.dat,
// never restated as constants here. The schedule already lives in
// Baseline/generateRomanBaseline.py; a second copy would drift silently, and it is
// expected to change (OPEN_ITEMS.md: the GBTDS footprint and cadence are still being
// reconciled against STScI's current pages). Deriving means the C++ can never disagree
// with the visit list it is actually integrating.
//
// Clustering rule and the SEASON_GAP_MIN_DAYS constant: see config/parameters.h, section 5.

// Where t0 sits relative to Roman's mission. Three states, not two.
//
// An event peaking before Roman launches, or after it ends, is Rubin-only BY
// CONSTRUCTION: there is no Roman data anywhere near it and nothing for a joint fit to
// rescue. An event peaking in a MID-MISSION gap is a completely different object --
// Roman brackets it, with dense photometry on both sides, so a long-tE event's wings
// are still measured even though its peak was missed. That second case is the one the
// joint-fit science claim is about. Collapsing the two into a single "Roman had no data
// at t0" boolean would dilute the headline result with events that were never
// candidates, which is the easiest available way to wash the effect out.
enum T0Zone {
    T0_IN_SEASON   = 0, //Roman was observing at t0
    T0_IN_GAP      = 1, //between two Roman seasons -- the gap-filling regime
    T0_OFF_MISSION = 2, //before Roman's first epoch or after its last
};

struct RomanSchedule {
    std::vector<std::pair<double,double>> seasons; //[start, end] in simulation days
    double missionStart       = 0.0;
    double missionEnd         = 0.0;
    double maxInSeasonSpacing = 0.0; //largest spacing kept INSIDE a season
    double minSeasonGap       = 0.0; //smallest spacing treated as a gap
    double minSeasonLength    = 0.0; //shortest season found, end - start [d]. Zero means some
                                     //"season" holds a single epoch, which is not a season at
                                     //all -- see the guard in main(). A schedule sampled more
                                     //coarsely than SEASON_GAP_MIN_DAYS degenerates that way
                                     //while leaving the other two margins looking healthy.

    // Signed days from t0 to the nearest season boundary.
    //   negative -> t0 is INSIDE a season; |value| is how deep into it the peak sits
    //   positive -> t0 is outside every season; value is the distance to the nearest edge
    // Signed this way so the gap-filling plot reads left to right: x < 0 is "Roman was
    // watching", x > 0 is "Roman was not", and the joint-over-Roman precision gain is
    // expected to grow with x. Use zone() to tell a mid-mission gap from off-mission;
    // both give a positive dtToSeasonEdge and they must not be pooled.
    double dtToSeasonEdge(double t0) const;
    int    zone(double t0) const;
    int    seasonOf(double t) const;   // index into seasons containing t, or -1 (Deviation 71)
};

// Cluster ro.tim into seasons. One pass over a sorted, de-duplicated copy of the epoch
// times; called once per run, not per event.
RomanSchedule buildRomanSchedule(const roman& ro);

// Column names of the per-event table, in write order. Kept next to the writer's data so
// the two cannot drift: a header that disagrees with its columns is worse than none.
const char* eventTableHeader();
// ---------------------------------------------------------------------------


// Maps a per-epoch telescope tag (lens::tele[i]: 0 = Rubin, 1 = Roman) to its survey index.
inline int surveyOfTele(int tele) { return (tele == 0) ? SRUBIN : SROMAN; }

// Which photometric parameters a given survey partition can actually constrain.
//
// Rubin's epochs carry no information whatsoever about Roman's flux parameters (fb1, mbs1) and
// vice versa: perturbing them leaves the model magnitude of the other telescope's epochs exactly
// unchanged, so those rows and columns of the single-survey information matrices are identically
// zero. Inverting the full Nx x Nx matrix for a single survey would therefore be singular by
// construction, not by accident. Each partition instead inverts only the submatrix it can
// constrain; parameters outside the subset report sigma = -1.
//
// The joint matrix keeps all Nx parameters, which is the point -- it is the only one that sees
// both telescopes' flux scales at once, and the chromatic difference between them is part of what
// breaks the u0-tE-fb degeneracy.
// A parameter is only active for a partition if that partition's data can constrain it.
//
// Rubin's epochs carry no information about Roman's flux parameters (fb1, mbs1) and vice versa:
// perturbing them leaves the other telescope's model magnitudes exactly unchanged, so those rows
// and columns are identically zero. Inverting the full Nx x Nx matrix for such a partition would
// be singular by construction, not by accident.
//
// The joint set additionally depends on which surveys actually contributed epochs for THIS event.
// A short event peaking in a Roman gap has no Roman data at all, so the joint fit cannot solve for
// Roman's flux scale either and must fall back to Rubin's parameter set. Without this the joint
// matrix would go singular on exactly the gap-peaking events the project is about.
// kMinTeleEpochs (the minimum epochs a telescope needs to enter the photometric matrices,
// Deviation 76): see config/parameters.h, section 7.

inline std::vector<int> activePhotParams(int surv, int nRubinEpochs, int nRomanEpochs)
{
    const std::vector<int> shared = {0, 1, 3, 4, 5};  //u0, tE, piE, xi, t0
    const std::vector<int> rubinFlux = {2, 6};        //fb0, mbs0
    const std::vector<int> romanFlux = {7, 8};        //fb1, mbs1

    std::vector<int> out = shared;
    const bool wantRubin = (surv == SRUBIN) || (surv == SJOINT && nRubinEpochs > 0);
    const bool wantRoman = (surv == SROMAN) || (surv == SJOINT && nRomanEpochs > 0);
    if (wantRubin) out.insert(out.end(), rubinFlux.begin(), rubinFlux.end());
    if (wantRoman) out.insert(out.end(), romanFlux.begin(), romanFlux.end());
    std::sort(out.begin(), out.end());
    return out;
}

struct covarian {
    int    sign;
    int    flagi;
    
    double deter;
    double sigmul1, sigmul2, f1, f2;
    double magw;
    double derm1f, derm2f, dera1f, derb1f, dera2f, derb2f, diff;

    std::array<double, nq> resu;
    std::array<double, 2> derm1, derm2, dera1, derb1, dera2, derb2, bb;
    std::vector<double> Delta1, Delta2; //size Nx, Ny

    // One photometric and one astrometric Fisher matrix per survey subset (see SurveyIdx).
    // Index with SJOINT / SRUBIN / SROMAN.
    std::array<gsl_matrix_uptr, NSURV> inputA; //size Nx each
    std::array<gsl_matrix_uptr, NSURV> inverA; //size Nx each
    std::array<gsl_matrix_uptr, NSURV> inputB; //size Ny each
    std::array<gsl_matrix_uptr, NSURV> inverB; //size Ny each

    std::array<std::vector<double>, NSURV> Era; //size Nx each: 1-sigma photometric
    std::array<std::vector<double>, NSURV> Erb; //size Ny each: 1-sigma astrometric

    // Step 3c (Deviation 71): the astrometric forecast under the three noise variants (AV_W, AV_N,
    // AV_P; see AST_SIGC in this file). Index [variant][survey][param]. Variant W is identical to
    // Erb/okB/condB/relMl; the Rubin partition is identical in all three (only Roman's noise and
    // offsets differ). -1 sentinels as everywhere else.
    std::array<std::array<std::array<double, Ny>, NSURV>, NAVAR> ErbV{};
    std::array<std::array<int, NSURV>, NAVAR>    okBV{};
    std::array<std::array<double, NSURV>, NAVAR> condBV{};
    std::array<std::array<double, NSURV>, NAVAR> relMlV{};
    std::array<std::array<std::array<double, Ny * Ny>, NSURV>, NAVAR> FBV{};  // the matrices

    // Per-survey epoch counts and characterizability flags. A partition with no epochs at all
    // (a short event peaking in a Roman gap genuinely has no Roman data), or with fewer epochs
    // than free parameters, is rank-deficient by construction: it must be reported as
    // not-characterizable, NOT inverted into a meaningless ~1e10 sigma.
    std::array<int, NSURV> nepochA; //epochs contributing to each photometric matrix
    std::array<int, NSURV> okA;     //1 = inverted and usable, 0 = not characterizable
    std::array<int, NSURV> okB;     //same for the astrometric matrix

    // Condition number (lambda_max/lambda_min) of the NORMALIZED information matrix, per survey.
    // Normalized, not raw: dividing row and column i by sqrt(F_ii) removes the spread that comes
    // merely from the parameters being expressed in different units (tE in days ~30, u0
    // dimensionless ~0.3, xi in radians), leaving only genuine parameter degeneracy. A large
    // value after normalization is a physical statement -- some combination of parameters is
    // unconstrained by this data, classically the u0-tE-fb degeneracy -- not a units artifact.
    // Stored rather than only thresholded so the cut can be chosen during analysis.
    // -1.0 means "not computed" (partition rejected before it got this far).
    std::array<double, NSURV> condA;
    std::array<double, NSURV> condB;

    // Fractional 1-sigma on the lens mass, per survey (Step D1). Ml is not a fitted
    // parameter: it follows from Ml = tetE / (kappa * piE), a pure ratio, so the two
    // fractional errors add in quadrature. Stored per survey because its two ingredients
    // come from different instruments' different strengths -- tetE from the ASTROMETRIC
    // matrix (sub-mas centroid motion: Roman) and piE from the PHOTOMETRIC one over a
    // long time baseline (annual parallax distortion: Rubin). A mass the joint fit
    // measures and neither survey measures alone is the black-hole result.
    // -1.0 means "not measurable in this partition", never a measured value.
    std::array<double, NSURV> relMl;

    // Multiplicative override on each photometric parameter's finite-difference step (Step C3).
    // Default 1.0 is an exact no-op, so production behaviour is byte-identical; the step-size
    // convergence sweep drives these at runtime. A runtime knob rather than a compile flag
    // because a sweep needs many step values within a single process -- a compile flag would
    // mean one rebuild, and for the live binary one multi-minute CMD reload, per sweep point.
    std::array<double, Nx> deltaScale;
    // Same knob for the astrometric steps Delta2[] (Deviation 75, step M3): 1.0 = production.
    std::array<double, Ny> deltaScaleB{1.0, 1.0, 1.0, 1.0};

    gsl_matrix_uptr summA; //size Nx (diagnostic only, joint)
    gsl_matrix_uptr summB; //size Ny (diagnostic only, joint)

    // Constructor
    covarian()
        : Delta1(Nx), Delta2(Ny),
          summA(gsl_matrix_alloc(Nx, Nx)),
          summB(gsl_matrix_alloc(Ny, Ny))
    {
        for (int q = 0; q < NSURV; ++q) {
            inputA[q].reset(gsl_matrix_alloc(Nx, Nx));
            inverA[q].reset(gsl_matrix_alloc(Nx, Nx));
            inputB[q].reset(gsl_matrix_alloc(Ny, Ny));
            inverB[q].reset(gsl_matrix_alloc(Ny, Ny));
            if (!inputA[q] || !inverA[q] || !inputB[q] || !inverB[q]) {
                throw std::runtime_error("GSL matrix allocation failed");
            }
            Era[q].assign(Nx, 0.0);
            Erb[q].assign(Ny, 0.0);
            nepochA[q] = 0;
            okA[q] = 0;
            okB[q] = 0;
            condA[q] = -1.0;
            condB[q] = -1.0;
            relMl[q] = -1.0;
        }
        for (int k = 0; k < Nx; ++k) {
            deltaScale[k] = 1.0;
        }
        if (!summA || !summB) {
            throw std::runtime_error("GSL matrix allocation failed");
        }
    }

//        deter = 0; sign = 0; flagi = 0;
//        sigmul1 = sigmul2 = f1 = f2 = 0;
//        magw = derm1f = derm2f = dera1f = derb1f = dera2f = derb2f = diff = 0;

//        for(int i=0; i<2; ++i) {
//            derm1[i]=derm2[i]=dera1[i]=derb1[i]=dera2[i]=derb2[i]=bb[i]=0;
//        }
//        for(int i=0; i<nq; ++i) resu[i] = 0;
//    }

    // Disable copy, allow move
    covarian(const covarian&) = delete;
    covarian& operator=(const covarian&) = delete;
    covarian(covarian&&) = default;
    covarian& operator=(covarian&&) = default;
};
/*
struct yfilter{
    std::vector<double> Age; //Size YZ
    std::vector<double> B;   //Size YZ
    std::vector<double> M;   //Size YZ
    std::vector<double> mm;  //Size YZ

    std::vector<int> number;   //Size met
    std::vector<int> count;    //Size met
    std::vector<double> Metal; //Size met

        // Constructor
    yfilter()
        : Age(YZ), B(YZ), M(YZ), mm(YZ),
          number(met), count(met), Metal(met)
    {}
};
*/
//==========================================//
// How the three Fisher partitions came out for one event.
//
// This exists so that events where a survey contributes real information but cannot characterize
// the event on its own are LABELLED rather than silently lost. They are the strongest evidence
// for the joint fit: data that is insufficient alone still sharpens the combined result. A naive
// analysis computing sigma_joint / sigma_roman would hit a not-characterizable sentinel on
// exactly these events and, if it dropped the row, would discard the best synergy cases and bias
// the reported gain downward -- the mirror image of the selection bias the plan warns about at
// Step C5. Never drop a row on the basis of a missing single-survey sigma; classify it.
// Detection taxonomy. A microlensing event is only meaningfully "detected" if the joint fit
// detects it: the joint stream contains strictly more data than either survey alone, so a
// single-telescope detection that the joint test misses is not a real category but a symptom
// of an inconsistent threshold (see DET_ANOMALY below).
//
// That leaves four ways an event can be detected, distinguished by which surveys ALSO detect
// it on their own -- which is the quantity the joint-fit science case is about:
enum DetClass {
    DET_NONE         = 0, //nothing detected it
    DET_JOINT_ONLY   = 1, //only the combined stream -- neither telescope alone would have found it
    DET_RUBIN_JOINT  = 2, //Rubin alone, and the joint fit
    DET_ROMAN_JOINT  = 3, //Roman alone, and the joint fit
    DET_BOTH_JOINT   = 4, //both telescopes alone, and the joint fit
    DET_ANOMALY      = 5, //a telescope detected it but the joint test did NOT -- see below
    NDETCLASS        = 6,
};

// DET_ANOMALY must stay empty. Adding data cannot destroy signal, so if either survey alone
// clears its detection bar the combined stream must clear its own. That this counter is NOT
// currently zero is a real finding, not a bookkeeping artifact: the detection test compares
// dchi against 2*ndw, i.e. it thresholds the MEAN per-epoch chi-squared improvement. Pooling a
// survey with many low-signal epochs therefore raises the joint bar without adding signal, and
// can veto a detection the other survey made alone. Counted explicitly rather than folded into
// a neighbouring class, so the inconsistency stays visible instead of being silently absorbed.
inline int detClass(int detL, int detR, int detJ)
{
    if (!detJ) return (detL or detR) ? DET_ANOMALY : DET_NONE;
    if (detL and detR) return DET_BOTH_JOINT;
    if (detL)          return DET_RUBIN_JOINT;
    if (detR)          return DET_ROMAN_JOINT;
    return DET_JOINT_ONLY;
}

inline const char* detClassName(int c)
{
    switch (c) {
        case DET_NONE:        return "none";
        case DET_JOINT_ONLY:  return "joint-only";
        case DET_RUBIN_JOINT: return "Rubin+joint";
        case DET_ROMAN_JOINT: return "Roman+joint";
        case DET_BOTH_JOINT:  return "both+joint";
        case DET_ANOMALY:     return "ANOMALY(single-not-joint)";
        default:              return "?";
    }
}

enum SynergyClass {
    SYN_NONE       = 0, //joint not characterizable either -- no usable Fisher information
    SYN_BOTH_ALONE = 1, //both surveys characterize alone; joint/single ratio defined both ways
    SYN_RUBIN_ONLY = 2, //only Rubin alone; Roman's epochs still sharpen the joint fit
    SYN_ROMAN_ONLY = 3, //only Roman alone; Rubin's epochs still sharpen the joint fit
    SYN_JOINT_ONLY = 4, //NEITHER survey alone, but the joint fit works -- pure joint-fit rescue
};

// Classify from the photometric characterizability flags. SYN_JOINT_ONLY and the two
// *_ONLY classes are the scientifically interesting ones; see the note above.
inline int synergyClass(const covarian& co)
{
    if (!co.okA[SJOINT])                    return SYN_NONE;
    if ( co.okA[SRUBIN] &&  co.okA[SROMAN]) return SYN_BOTH_ALONE;
    if ( co.okA[SRUBIN] && !co.okA[SROMAN]) return SYN_RUBIN_ONLY;
    if (!co.okA[SRUBIN] &&  co.okA[SROMAN]) return SYN_ROMAN_ONLY;
    return SYN_JOINT_ONLY;
}

struct EventRecord {
    int    counter, flagL;
    double tE, RE, piE, tetE, Vt, u0, Ml;
    double opt, Dl, Ds, vl, vs;
    double mbase, fblend;
    int    gg, struc;
    double FWHM, vsave, DeltaT, murel;
    double resu0, resu1, resu2, resu3, resu5, resu9, resu10, resu13, resu14;
    double Map2, nsbl2;
    int    flagi;
    double Ai2;

    // ---- per-survey bookkeeping (Step C5) ----
    // Without these, nothing downstream can tell a Rubin-only detection from a Roman-only one,
    // nor apply the "both surveys have data" cut that the joint-gain ratio requires. Appended
    // rather than interleaved so the positional aggregate initialiser above stays valid.
    int    ndwL, ndwR;          //epochs actually contributed by each survey
    int    detL, detR, detJ;    //independent per-survey and joint detection outcomes
    int    okJoint, okRubin, okRoman;   //1 = that Fisher partition is characterizable
    double sigtE_J,   sigtE_L,   sigtE_R;    //absolute 1-sigma on tE   [days]
    double sigpiE_J,  sigpiE_L,  sigpiE_R;   //absolute 1-sigma on piE  []
    double sigtetE_J, sigtetE_L, sigtetE_R;  //absolute 1-sigma on tetE [mas]
    int    detCls;              //DetClass -- which combination of surveys detected it
    int    synClass;            //SynergyClass -- see the note above the enum. Do NOT drop rows
                                //with a not-characterizable single-survey sigma; use this.
    double condA_J, condA_L, condA_R;   //photometric condition numbers (normalized matrix);
                                //-1 where that partition was not characterizable. Stored so the
                                //ill-conditioning cut can be chosen in analysis, not baked in.

    // ---- the rest of the per-event row (Step D1) ----
    // Appended, never interleaved: the positional aggregate initialiser at the push_back
    // site depends on this order, and so does every column index downstream.
    double t0;                  //time of closest approach [d] on the simulation clock
                                //(day 0 = 2026-04-11, the first Rubin bulge visit)
    double xi;                  //source-trajectory angle [rad]; sets how the annual
                                //parallax ellipse projects onto the source track
    double lon, lat;            //Galactic sightline [deg]. Present so "per field" is a cut
                                //on this table rather than 1706 separate files.
    double mbs1, fb1;           //Roman F146 baseline magnitude [mag] and source-flux
                                //fraction []. mbs0/fb0 above are Rubin's. Today these
                                //duplicate magb[6]/blend[6] and magb[2]/blend[2]
                                //respectively, because RUBIN_BANDS is {r} -- they are the
                                //quantities the Fisher matrix actually fits (params 6-8),
                                //and they stop duplicating the moment a band is added.
    std::array<double, M> magb; //blended baseline magnitude per FILTER [mag]
                                //(0-5 = LSST ugrizy, 6 = Roman F146) -- all stars in the
                                //seeing disc, not just the source
    std::array<double, M> blend;//fraction of aperture flux from the SOURCE, per filter []
    double relMl_J, relMl_L, relMl_R;  //fractional sigma(Ml); -1 = not measurable
    int    okB_J, okB_L, okB_R; //astrometric matrix inverted? NOT implied by okA: a row can
                                //have okRubin=1 while sigtetE_L is -1, and without this
                                //nothing in the row explains why.
    double condB_J, condB_L, condB_R;  //normalized astrometric condition numbers; -1 as above
    double dtEdge;              //signed days from t0 to the nearest Roman season edge;
                                //NEGATIVE means t0 fell inside a season. See RomanSchedule.
    int    t0zone;              //T0Zone: 0 in-season, 1 mid-mission gap, 2 off-mission

    // ---- Step R1: resolving the two images, counted per RECORDED EPOCH ----
    // Number of that survey's epochs at which both images were within the filter's
    // [saturation, single-visit depth] AND separated by more than the stated bar. Counts, not
    // fractions: the paper's criterion is "at least three such epochs", which only a count can
    // answer. nres5/nres20 use D*sigma_a at the paper's two anchors; nresPSF uses the filter's
    // PSF FWHM. See the Step R1 block above for why all three are kept.
    long   nres5_L, nres20_L, nresPSF_L;
    long   nres5_R, nres20_R, nresPSF_R;
    // Largest separation reached at an epoch where both images were detectable [mas].
    // -1 means that never happened -- a SENTINEL, never a measurement. It must not enter a
    // mean or a histogram, exactly like every other -1 in this row.
    double dsepMax_L, dsepMax_R;
};
///===================== FUNCTION ===========================================//
int    Funcu0(lens & l);
int    FunctE(lens & l);
int    FuncMl(lens & l);
int    FuncPi(lens & l);
int    FuncMu(lens & l);
int    FuncMb(lens & l, double);
int    FuncFb(lens & l, double);
int    nearestSightline(const extin& ex, double lon, double lat);

void   read_cmd(CMD & cm);
void   readExtinction(extin& ex, const std::string& path); // files/ext/ext_tables.dat; exits on bad input
void   optical_depth(source & s);
void   func_source(source & s, CMD & cm, const extin& ex, int sightlineIdx);
void   func_lens( lens & l, source & s, const extin & ex, int sightlineIdx);
void   vrel(source & s, lens & l);
void   Disk_model(source & s, int);
void   ErrorCal(covarian & co, lens &l, source &s);
void   lightcurve(source & s, lens & l, astromet & as, double, int tele);
void   FisherM(source & s, lens & l, astromet & as,  covarian & co, int);
//double ylsst(yfilter & yf, double , double , double);
double interpExtinctionAlongSightline(const extin& ex, int k, double dist);
double errlsstM(double,int,double);
double errlsstA(lsst & ls,  double);
double errRomanA(double magF146);   //Roman WFI per-exposure astrometric error [mas] (Step H4)
//double errELT(lsst & ls,double,int);
// TODO(Ali): wire this to whatever sigma_roman.txt actually represents (a fixed
// mag-vs-error lookup, or a per-visit-depth-dependent formula like errlsstM). Signature
// below assumes the simpler case (no per-visit depth); adjust if you need `sig5`.
double errRomanM(const roman & ro, double mag);

// matchVisibleEpochs is a template (the coverage test is a predicate) and is defined in
// Bulge_LSST.cpp, its only caller.

// LSSTCam's footprint (Deviation 80): the active-silicon map OpSim/MAF use (rubin_scheduler
// LsstCameraFootprint, fov_map.npz), exported to Baseline/lsstcam_fov/fov_map.txt. A sky point is on
// a Rubin visit's silicon if its gnomonic projection about the boresight, rotated by rotSkyPos,
// falls on an active pixel -- rubin_scheduler's own algorithm, mirrored. Replaces the 1.75-deg circle.
// RUBIN_MAX_RADIUS (1.94 deg, rubin_scheduler's max_radius): see config/parameters.h, section 4.
struct LsstCamMap { int n = 0; double x0 = 0.0, step = 0.0; std::vector<unsigned char> on; };
inline LsstCamMap gLsstCam;
void   readLsstCamMap(const std::string& path);
bool   onLsstCam(double ra, double dec, double ra0, double dec0, double rotSkyPos);   // all deg
void   galToIcrs(double l, double b, double& ra, double& dec);                       // deg (J2000)

// GBTDS detector layout (helper.cpp). readGbtdsLayout exits with the file named on any
// malformed input; inDetector tests a sky offset (dl, db) from a field centre against one
// layout's 18 detector rectangles.
GbtdsLayout readGbtdsLayout();
bool inDetector(const GbtdsLayout& g, int layout, double dl, double db);

double CCM89_a(double lambda_um);
double CCM89_b(double lambda_um);
double AlAv(double lambda_um, double Rv);
//void   getCofactorA(double input[Nx][Nx], double temp[Nx][Nx], int , int , int );
//void   getCofactorB(double input[Ny][Ny], double temp[Ny][Ny], int , int , int );
//double determinantA(double input[Nx][Nx], int);
//double determinantB(double input[Ny][Ny], int);
//void   inverse(covarian & co, int );
// Inverts one of the Fisher matrices in place. `flag` selects photometric (0, Nx) or
// astrometric (1, Ny); `surv` selects which SurveyIdx partition. Returns 1 if the matrix was
// non-singular and the inverse is usable, 0 if it was singular -- in which case the caller must
// treat that partition as not-characterizable rather than reading numbers out of it.
int    invert_matrix(covarian & co, int flag, int surv);
int    invertNormalized(const gsl_matrix* in, gsl_matrix* out, const std::vector<int>& act,
                        double& cond, double* deter);   // invert_matrix's core (Deviation 71)
void   print_mat_contents(gsl_matrix *matrix, int);

double RandN(double , double);
double RandR(double , double);
int    RandPois(double);
double drawKroupaInitialMass();
double remnantMass(double initialMass);
double drawLogUniformMass(double lo, double hi);
double drawNeutronStarMass();
double drawPowerLawMass(double lo, double hi, double alpha);
double drawLensMass(bool* luminous = nullptr);   // *luminous: a living star (Deviation 74)

// Luminous lenses (Deviation 74): main-sequence mass -> absolute magnitude (AB, ugrizy + F146) per
// Galactic component, from CMD/components/lens_ml.dat (CMD/lens_ml_table.py). Bins of 0.02 Msun
// over 0.08-1.00 Msun; lighter lenses are brown dwarfs and treated as dark.
struct LensML {
    std::array<std::vector<double>, 4> mmid;
    std::array<std::vector<std::array<double, 7>>, 4> mab;
};
inline LensML gLensML;
void   readLensML(const std::string& path);
bool   lensAbsMag(int comp, double mass, std::array<double, 7>& mab);   // false if dark

#endif // LMC_H
