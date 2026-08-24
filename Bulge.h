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

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

////
#include <random>
constexpr int seed = 42;
inline std::mt19937_64 rng{seed};
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


constexpr int    IMnum = 2;
constexpr double u0m   = 3.0;


constexpr double RAa = 180.0 / M_PI;
constexpr double pi  = M_PI;
constexpr double binary_fraction = double(2.0 / 3.0);
constexpr double velocity = 299792458.0;//velosity of light
constexpr double Msun = 1.98892 * std::pow(10., 30); //in [kg].
constexpr double Rsun = 6.957 * std::pow(10.0, 8.0); ///solar radius [meter]
constexpr double KP = 3.08568025 * std::pow(10., 19); // in meter.
constexpr double G = 6.67384 * std::pow(10., -11.0);// in [m^3/s^2*kg].
constexpr double AU = 1.4960 * std::pow(10.0, 11.0);
constexpr double vro_sun = 226.0;
constexpr double VSunR = 11.1;
constexpr double VSunT = vro_sun*(1.00762 + 0.00712) + 12.24;
constexpr double VSunZ = 7.25;
constexpr double year = 365.2425;//days
constexpr double eps = double(0.000000000000005463263454624313654);

///============================ Besancon constant ==========================///
constexpr double Dsun = 8.0;
constexpr std::array<double, 8> rho0 = {4.0, 7.9, 6.2, 4.0, 5.8, 4.9, 6.6, 3.96}; //considering WD
constexpr std::array<double, 8> d0   = {0.073117, 0.0216524, 0.0217405, 0.0217901, 0.0218061, 0.0218118, 0.0218121, 0.0218121};
constexpr std::array<double, 8> epci = {0.014, 0.0268, 0.0375, 0.0551, 0.0696, 0.0785, 0.0791, 0.0791};
constexpr std::array<double, 8> corr = {1.0, 7.9/4.48419, 6.2/3.52112, 4.0/2.27237, 5.8/3.29525, 4.9/2.78402, 6.6/3.74991, 3.96/2.24994};
constexpr std::array<double, 4> Rv   = {3.1, 2.5 ,3.1 ,3.1}; //Disk, Bulge, Thick, Halo

///=============================LSST constant===============================///
constexpr int    M = 6 + 1;    //No. of filter  ugrizy  of LSST + Roman's F146 filter
constexpr double Tobs = 10.0 * year;///LSST observational time 10 years
constexpr double delta2 = 0.005;///systematic errors


// TODO: ADD 7th value for F146 filter
constexpr std::array<double, M> gama = {0.037,0.038,0.039,0.039,0.040,0.040}; //a function of sky brightness and airmass in wavelengths.
constexpr std::array<double, M> seeing = {0.77,  0.73,  0.70,  0.67,  0.65,  0.63}; //seeing
constexpr std::array<double, M> msky   = {22.9,  22.3,  21.2,  20.5,  19.6,  18.6};
constexpr std::array<double, M> Cm     = {22.92, 24.29, 24.33, 24.20, 24.07, 23.69};
constexpr std::array<double, M> Dci    = {0.67,  0.21,  0.11,  0.08,  0.05,  0.04};
constexpr std::array<double, M> km     = {0.451, 0.163, 0.087, 0.065, 0.043, 0.138}; //sky extinction


constexpr std::array<double, M> sigma = {0.022, 0.02, 0.017, 0.017, 0.027, 0.027, 0.04}; // PLACEHOLDER: K-band value, not F146
constexpr std::array<double, M> thre  = {23.4, 24.6, 24.3, 23.6, 22.9, 21.7, 29.0}; //depth of single visit in ugrizy + F146 Filter (value needs to change)
constexpr std::array<double, M> satu  = {15.2, 16.3, 16.0, 15.3, 14.6, 13.4, 12.0}; //saturation limit of single visit in ugrizy + F146 Filter (value needs to change)
constexpr std::array<double, M> FWHM  = {1.22087, 1.10136, 0.993103, 0.967076, 0.951766, 0.936578, 0.105}; //LSST [arcsec] ugrizy + F146 Filter
//constexpr std::array<double, M> a0    = {0.9429, 1.0138, 0.94027, 0.8139, 0.6641, 0.5703, 0.1615}; //for calculating the extinction + F146 Filter (value needs to change)
//constexpr std::array<double, M> b0    = {1.9788, 0.5575, -0.2197, -0.4982, -0.6097, -0.5236, -0.1483}; // PLACEHOLDER: K-band value, not F146
constexpr std::array<double, M> lambda_um = {0.367, 0.482, 0.622, 0.755, 0.869, 0.971, 1.464};

// Which LSST filter(s) (indices 0-5 = u,g,r,i,z,y) form the single "Rubin representative
// band" standing in for the whole LSST light curve in the Fisher matrix, the recorded
// per-epoch Rubin model magnitude, and the LSST astrometric-error evaluation magnitude.
// {2} = r-band only, matching the pre-existing hardcoded behavior. Listing more than one
// index combines them by *summing* their baseline fluxes and source fluxes separately
// (see Lensing.cpp) -- an equal-weighted flux sum, not throughput-weighted (no per-filter
// throughput curve exists in this codebase yet). Does NOT change which real filter's noise
// model (errlsstM) applies to a given epoch -- that always reflects the epoch's own actual
// filter, regardless of this setting.
inline const std::vector<int> RUBIN_REF_BANDS = {2};

constexpr double cade1 = 3.0 ;//LSST[days]
//constexpr double cade2 = 10.0;//ELT [days]

//constexpr int YZ = 3578;   //No.yzma.txt rows
//constexpr int met = 70;   //No rows metal.txt
//constexpr int Nel = 7;//rows in "sigma_ELT.txt" (https://academic.oup.com/mnras/article/494/3/4413/5813442)
//constexpr int nfiles = 2518;
//constexpr int nlines = 3686;
//constexpr int nex = 2518 * 3686; //number of ext files * lines in each file
constexpr int NFILES = 2518; //number of extinction files
constexpr int NROWS = 3686; //number of rows in extintion files
constexpr int nrd = 10000; //rows in "convert_coordinate_2.dat"
constexpr int Na = 96;     //rows in "sigmaA_LSST.txt"
constexpr int NaRoman = 123;  // rows in sigma_roman.txt
constexpr int nq = 15;     //resu
constexpr int N1 = 396593, N2 = 3568010, N3 = 646090, N4 = 3171; //CMD_BESANCON: ThinDisk, Bulge, ThickDisk, Halo
// Data rows in BulgeBaseline.dat, EXCLUDING the header. Regenerated 2026-08-22 from
// baseline_v5.1.0_10yrs.db; readbaselineBulge.py prints the value to use here. The
// previous 7373 counted a doubled file (append-mode bug) and read 3687 phantom rows.
constexpr int Nl = 3686;

// Data rows in RomanBaseline.dat, EXCLUDING the header. generateRomanBaseline.py prints
// the value to use here; it must be updated whenever the season pattern, cadence or
// mission start day changes, or the read guard in Bulge_LSST.cpp will fire.
// Current: 6 high-cadence seasons (F146 every 12.1 min) + 4 low-cadence (every 5 days),
// on STScI's real alternating spring/fall visibility windows.
constexpr int NlRoman = 302406;

// `coun` bounds the length of one event's light-curve buffers (lens::timn/magn/errm/
// soux/souy/erra/tele). `ndw`, the running count of accepted epochs (Rubin + Roman
// combined) for a single star, can never exceed the total number of visit rows that
// exist for either instrument — Nl + NlRoman is therefore an unconditionally safe
// upper bound at any sky position, in any season. One-time memory cost since `lens`
// is allocated once, not per event: (Nl+NlRoman) * 7 buffers * 8 bytes ~= 17.7 MB.
constexpr int    coun  = Nl + NlRoman;

// Roman WFI field of view is much smaller than Rubin's and covers a small number of
// discrete GBTDS fields, not a rolling-cadence footprint. Do not reuse `FoV` for Roman.
// Matching RADIUS in degrees, NOT an area: matchVisibleEpochs tests
// sqrt(dl^2 + db^2) <= FoVRoman. STScI's GBTDS design page quotes 1.7 deg^2 monitored over
// six WFI fields = 0.2833 deg^2 per field, so the equal-area radius is
// sqrt(0.2833/pi) = 0.3003 deg. The previous 0.28 took a deg^2 figure as if it were already
// a radius, giving pi*0.28^2 = 0.2463 deg^2 -- about 13% too small. This value decides which
// sightlines see Roman at all, so it matters more than its size suggests.
constexpr double FoVRoman = 0.3003;

constexpr double tetp   = double(M_PI / 3.0);        //parallax
constexpr double omegae = double(2.0 * M_PI / year); //radian per day
constexpr double vearth = omegae;                    //radian per day
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


constexpr int GG = 100;
constexpr double tE_min  = 0.0;///days
constexpr double tE_max  = 50.0*year;//days
constexpr double Ml_min  = 3.0;
constexpr double Ml_max  = 5000.0;
constexpr double pi_min  = -0.45;
constexpr double pi_max  = 0.85;
constexpr double u0_min  = 0.0;
constexpr double u0_max  = u0m;
constexpr double mb_min  = 15.0;
constexpr double mb_max  = 26.0;
constexpr double fb_min  = 0.0;
constexpr double fb_max  = 1.0;
constexpr double mu_min  = 0.0; //mu_relative
constexpr double mu_max  = 100.0;

////=================================== Bulge ====================================
constexpr int    Num  = 9500;
constexpr double MaxD = 12.0; //kpc
constexpr double step = double(MaxD / Num / 1.0); //step in kpc
//const double RaLMC  =  80.89375;
//const double DecLMC = -68.2438888888889;
//const double DLMC =  49.97;///KPC
constexpr double DBulge = 8; //Kpc
constexpr double l1  = -0.219 - 0.2 - 3.5 / 2;
constexpr double l2  = 1.4134 + 0.2 + 3.5 / 2;
constexpr double lx  = 1.0053 - 0.2 - 3.5 / 2;
constexpr double b1  = -1.64  - 0.2 - 3.5 / 2;
constexpr double b2  = -0.85  + 0.2 + 3.5 / 2;
constexpr double bx  = -1.64  + 0.2 + 3.5 / 2;
constexpr double wid = 3.5 / 2;
constexpr double dd  = 0.02;
constexpr double FoV = double(3.5 / 2.0);  //the radius of teh Rubin Field of View

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
          tele(coun),

          tEs(make_grid(tE_min, tE_max)),
          Mls(make_grid(Ml_min, Ml_max)),
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

// One extinction curve for one line of sight
struct ExtinctionProfile
{
    std::array<double, NROWS> dist;
    std::array<double, NROWS> ext;
};

// One Galactic line of sight
struct Sightline
{
    double l;
    double b;

    ExtinctionProfile profile;
};

// Whole Bayestar dataset
struct extin
{
    std::array<Sightline,NFILES> sightlines;
};

/*
struct extin {
    std::vector<double> l;
    std::vector<double> b;
    std::vector<double> dist;
    std::vector<double> ext;

    extin()
        : l(nex), b(nex), dist(nex), ext(nex)
    {}
};
*/
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

    //ID  RA  Dec  l  b  start  filter  airmass  seeing  skyBrightness visittime sigma5 targetname distance
    lsst()
        : mag(Na), err(Na),
          filter(Nl),
          ct(Nl),
          RA(Nl), DEC(Nl), l(Nl), b(Nl), tim(Nl), sig5(Nl), dist(Nl)
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

    // NOTE: no `filter` array — currently only F146 (constant filter index 6) is modeled
    // for Roman. If F087/F213 are added later, give roman a `filter` array like lsst's.

    roman()
        : mag(NaRoman), err(NaRoman),
          ct(NlRoman),
          RA(NlRoman), DEC(NlRoman), l(NlRoman), b(NlRoman), tim(NlRoman), sig5(NlRoman)
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
// Clustering rule: consecutive distinct epoch times more than SEASON_GAP_MIN_DAYS apart
// begin a new season. This is safe by a wide margin on the current schedule -- the
// largest spacing INSIDE a season is 5.0 d (the low-cadence seasons' five-day sampling)
// and the smallest gap BETWEEN seasons is 108.2 d -- but the margin is checked at
// runtime rather than assumed; see the guard in main().
constexpr double SEASON_GAP_MIN_DAYS = 20.0;

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
void   readBayestar(extin& ex, const std::string& folder); // Read files
void   optical_depth(source & s);
void   func_source(source & s, CMD & cm, const extin& ex, int sightlineIdx);
void   func_lens( lens & l, source & s);
void   vrel(source & s, lens & l);
void   Disk_model(source & s, int);
void   ErrorCal(covarian & co, lens &l, source &s);
void   lightcurve(source & s, lens & l, astromet & as, double);
void   FisherM(source & s, lens & l, astromet & as,  covarian & co, int);
//double ylsst(yfilter & yf, double , double , double);
double interpExtinctionAlongSightline(const extin& ex, int k, double dist);
double errlsstM(double,int,double);
double errlsstA(lsst & ls,  double);
//double errELT(lsst & ls,double,int);
// TODO(Ali): wire this to whatever sigma_roman.txt actually represents (a fixed
// mag-vs-error lookup, or a per-visit-depth-dependent formula like errlsstM). Signature
// below assumes the simpler case (no per-visit depth); adjust if you need `sig5`.
double errRomanM(const roman & ro, double mag);

// Finds all epochs of a given instrument's baseline within `fov` of (lon,lat),
// filling `ct` with indices into the instrument's own tim/l/b arrays (sorted by time,
// since the baseline files are pre-sorted). Returns the number of visible epochs (ndd)
// and reports the smallest gap between consecutive visible epochs via minCadence.
// `label` ("LSST"/"Roman") is used only in diagnostic messages if tied/out-of-order
// timestamps are found and skipped.
int matchVisibleEpochs(const char* label, double lon, double lat, double fov,
                       const std::vector<double>& l_arr,
                       const std::vector<double>& b_arr,
                       const std::vector<double>& tim_arr,
                       int nEpochs,
                       std::vector<int>& ct,
                       double& minCadence);

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
void   print_mat_contents(gsl_matrix *matrix, int);

double RandN(double , double);
double RandR(double , double);

#endif // LMC_H
