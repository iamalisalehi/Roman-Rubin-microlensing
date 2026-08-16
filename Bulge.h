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

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

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

#define Nx 6
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
constexpr int Nl = 7373; //BulgeBaseline.dat 18/07/2026

// TODO(Ali): set to the actual row count of RomanBaseline.dat once it is generated
// from the ROTAC 2025 overguide season/cadence design (6 high-cadence + 4 low-cadence
// seasons; F146 ~12 min in high-cadence seasons, 3-day in low-cadence seasons).
constexpr int NlRoman = 309084;

// `coun` bounds the length of one event's light-curve buffers (lens::timn/magn/errm/
// soux/souy/erra/tele). `ndw`, the running count of accepted epochs (Rubin + Roman
// combined) for a single star, can never exceed the total number of visit rows that
// exist for either instrument — Nl + NlRoman is therefore an unconditionally safe
// upper bound at any sky position, in any season. One-time memory cost since `lens`
// is allocated once, not per event: (Nl+NlRoman) * 7 buffers * 8 bytes ~= 17.7 MB.
constexpr int    coun  = Nl + NlRoman;

// Roman WFI field of view is much smaller than Rubin's and covers a small number of
// discrete GBTDS fields, not a rolling-cadence footprint. Do not reuse `FoV` for Roman.
// TODO(Ali): replace with the actual per-field WFI FoV radius (deg) used in RomanBaseline.dat
constexpr double FoVRoman = 0.28; // PLACEHOLDER — Roman WFI FoV is ~0.28 deg^2 total

constexpr double tetp   = double(M_PI / 3.0);        //parallax
constexpr double omegae = double(2.0 * M_PI / year); //radian per day
constexpr double vearth = omegae;                    //radian per day
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
    
    std::vector<int> ct;        // 1000
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
          ct(1000),
          RA(Nl), DEC(Nl), l(Nl), b(Nl), tim(Nl), sig5(Nl), dist(Nl)
    {}
};

struct roman {
    std::vector<double> mag;   // NaRoman: mag-vs-error lookup (sigma_roman.txt)
    std::vector<double> err;   // NaRoman

    // --- New: per-visit epoch bookkeeping, mirrors lsst's fields ---
    std::vector<int> ct;        // 1000 — visible-epoch indices for the current sightline
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
          ct(1000),
          RA(NlRoman), DEC(NlRoman), l(NlRoman), b(NlRoman), tim(NlRoman), sig5(NlRoman)
    {}
};

struct covarian {
    int    sign;
    int    flagi;
    
    double deter;
    double sigmul1, sigmul2, f1, f2;
    double magw;
    double derm1f, derm2f, dera1f, derb1f, dera2f, derb2f, diff;

    std::array<double, nq> resu;
    std::array<double, 2> derm1, derm2, dera1, derb1, dera2, derb2, bb;
    std::vector<double> Era, Erb; //size Nx, Ny
    std::vector<double> Delta1, Delta2; //size Nx, Ny
    
    gsl_matrix_uptr inputA; //size Nx
    gsl_matrix_uptr inverA; //size Nx

    gsl_matrix_uptr inputB; //size Ny
    gsl_matrix_uptr inverB; //size Ny
   
    gsl_matrix_uptr summA; //size Nx
    gsl_matrix_uptr summB; //size Ny

    // Constructor
    covarian()
        : Era(Nx), Erb(Ny),
          Delta1(Nx), Delta2(Ny),
          inputA(gsl_matrix_alloc(Nx, Nx)),
          inverA(gsl_matrix_alloc(Nx, Nx)),
          inputB(gsl_matrix_alloc(Ny, Ny)),
          inverB(gsl_matrix_alloc(Ny, Ny)),
          summA(gsl_matrix_alloc(Nx, Nx)),
          summB(gsl_matrix_alloc(Ny, Ny))
    {
        if (!inputA || !inverA || !inputB || !inverB || !summA || !summB) {
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
void   invert_matrix(covarian & co, int);
void   print_mat_contents(gsl_matrix *matrix, int);

double RandN(double , double);
double RandR(double , double);

#endif // LMC_H
