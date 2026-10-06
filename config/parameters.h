#ifndef PARAMETERS_H
#define PARAMETERS_H

// config/parameters.h -- THE ONE PLACE TO CHANGE THE MODEL / SURVEY SETUP.
//
// Everything here is a hand-edited choice: the Galactic model, filters, the Rubin and Roman
// surveys and footprints, the lens populations and mass functions, the detection threshold,
// and the output grids. Values are COMPILE-TIME: re-run `make` after editing.
//
// Numbers that describe the data files themselves (row counts, catalogue mean masses) live in
// config/data_products.h instead; physical constants live in include/physical_constants.h.

#include <array>
#include <vector>
#include <cmath>

#include "physical_constants.h"
#include "parameter_types.h"
#include "data_products.h"   // LSST_AST_TABLE_FLOOR (generated from files/sigmaA_LSST.txt)

// ==========================================================================================
// (1) RUN DEFAULTS
// ==========================================================================================

constexpr int seed = 42;   // default base seed; --seed overrides it at run time (Deviation 77)
// Per-sightline Monte Carlo budgets (defaults of --events / --lenses / --nerr / --maxdraws; the do/while over
// stars stops once events, lenses and nerr are ALL met, or at the draw cap; see Bulge_LSST.cpp RunConfig).
constexpr int    DEFAULT_EVENTS_TARGET = 850;    //detected events per sightline (icon)
constexpr int    DEFAULT_LENSES_TARGET = 150;    //Fisher-characterised events per sightline (nlens)
constexpr double DEFAULT_NERR_TARGET   = 2.0;    //accumulated Fisher-error weight
constexpr double DEFAULT_MAXDRAWS      = 5.0e4;  //hard cap on stars drawn at one sightline
// Dense model-curve sampling of the s2 sample light curves (defaults of the sample spec keys step_te /
// span_te / dt_coarse, Bulge_LSST.cpp SampleSpec): fine steps of DEFAULT_S2_STEP_TE * tE within
// +-DEFAULT_S2_SPAN_TE * tE of the peak, steps of DEFAULT_S2_DT_COARSE days over the rest of the mission.
constexpr double DEFAULT_S2_STEP_TE   = 0.01;   //peak-window step, in units of tE
constexpr double DEFAULT_S2_SPAN_TE   = 3.0;    //peak-window half-width, in units of tE
constexpr double DEFAULT_S2_DT_COARSE = 2.0;    //step over the rest of the mission [days]

// ==========================================================================================
// (2) GALACTIC MODEL
// ==========================================================================================

constexpr double binary_fraction = double(2.0 / 3.0);
constexpr double vro_sun = 226.0;
constexpr double VSunR = 11.1;
// 1.00762 and 0.00712 are ROT_A and ROT_B (section 2b) evaluated at R = Dsun; kept as literals here so
// the rounding of this expression is exactly what it always was.
constexpr double VSunT = vro_sun*(1.00762 + 0.00712) + 12.24;
constexpr double VSunZ = 7.25;
///============================ Besancon constant ==========================///
constexpr double Dsun = 8.0;
constexpr std::array<double, 8> rho0 = {4.0, 7.9, 6.2, 4.0, 5.8, 4.9, 6.6, 3.96}; //considering WD
constexpr std::array<double, 8> d0   = {0.073117, 0.0216524, 0.0217405, 0.0217901, 0.0218061, 0.0218118, 0.0218121, 0.0218121};
constexpr std::array<double, 8> epci = {0.014, 0.0268, 0.0375, 0.0551, 0.0696, 0.0785, 0.0791, 0.0791};
constexpr std::array<double, 8> corr = {1.0, 7.9/4.48419, 6.2/3.52112, 4.0/2.27237, 5.8/3.29525, 4.9/2.78402, 6.6/3.74991, 3.96/2.24994};
constexpr std::array<double, 4> Rv   = {3.1, 2.5 ,3.1 ,3.1}; //Disk, Bulge, Thick, Halo
constexpr int    Num  = 9500;
constexpr double MaxD = 12.0; //kpc
constexpr double step = double(MaxD / Num / 1.0); //step in kpc
// Source distance draw: the distance-grid index is drawn uniformly in [SRC_IDX_MIN, Num - SRC_IDX_END_MARGIN]
// (Lensing.cpp, func_source), i.e. it skips the first 5 grid cells (~0 kpc) and the last 2.
constexpr double SRC_IDX_MIN = 5.0;
constexpr double SRC_IDX_END_MARGIN = 2.0;
constexpr double dd  = 0.02;   // native sightline grid step [deg]; the scan steps by stride*dd

// ---- (2a) Galactic model: density laws (Disk_model, Bulge_LSST.cpp) ----
// Mass fractions of each component are set by the *_NORM factors below; the per-population number
// densities then follow from Disk_model's rho / <m> with the MEANMASS_* of config/data_products.h.

// Star-count completeness factors per component.
constexpr double DENS_FD = 1.0;  //see the program mass_averaged.cpp.  we do not apply any limitation
constexpr double DENS_FB = 1.0;  //just stars brighter than V=11.5, but we change to consider all stars
constexpr double DENS_FH = 1.0;  //No limitation

// Thin disc (8 age bins; rho0, d0, epci, corr above are its per-bin parameters).
constexpr double THIN_RDD = 2.17;       //scale length of the old bins (ii > 0) [kpc]; 2.53; ///2.17;
constexpr double THIN_RHH = 1.33;       //hole scale length of the old bins (ii > 0) [kpc]; 1.32; //1.33;
constexpr double THIN_YOUNG_L1 = 25.0;  //young (ii == 0) bin: outer scale, exp(-rdi/25)
constexpr double THIN_YOUNG_L2 = 9.0;   //young (ii == 0) bin: inner scale, exp(-rdi/9)
constexpr double THIN_CORE = 0.25;      //softening inside sqrt(0.25 + rdi/R^2) of the old bins
constexpr double THIN_NORM = 1.2;       // totalmass= 4.25e10

// Thick disc.
constexpr double THICK_RHO00 = 1.34 * 0.001 + 3.04 * 0.0001;  //local density, Msun/pc^3
constexpr double THICK_RHO_DIV = 0.999719;                    //normalisation divisor of rho00
constexpr double THICK_SCALE_LEN = 2.5;                       //radial scale length [kpc]
constexpr double THICK_H1 = 0.4;                              //height of the parabolic core [kpc]
constexpr double THICK_H2 = 0.8;                              //exponential scale height [kpc]
constexpr double THICK_NNF = 0.4 / 0.8;                       //THICK_H1 / THICK_H2
constexpr double THICK_NORM = 2.67;                           //total_mass=0.8e10

// Stellar halo.
constexpr double HALO_FLATTEN = 0.76;                         //axis ratio
constexpr double HALO_CORE = 0.5;                             //core radius [kpc]
constexpr double HALO_RHO0 = (0.932 * 0.00001 / 867.067);     //local density, Msun/pc^3
constexpr double HALO_SLOPE = -2.44;                          //power-law index
constexpr double HALO_NORM = 5281.0;                          //Total_mass=1.2e9

// Bulge / bar: two triaxial components, S (boxy, 1/cosh^2 profile) and E (exponential).
constexpr double BAR_ANGLE_DEG = 12.89;                       //bar angle [deg]
constexpr double BAR_MASS_RESCALE = 0.24529; // calibrated so bulge column density toward
                                             // Baade's Window (l=1, b=-3.9) matches the
                                             // Han & Gould (2003) HST benchmark: 2086 Msun/pc^2
constexpr double BAR_CUTOFF_K = 4.0;                          //Gaussian cut-off exp(-K (r2-Rc)^2) beyond Rc
constexpr double BAR_NORM = 0.45;                             //total mass= 1.7e10
constexpr double BAR_S_RX0 = 1.46;
constexpr double BAR_S_RY0 = 0.49;
constexpr double BAR_S_RZ0 = 0.39;
constexpr double BAR_S_RC  = 3.43;
constexpr double BAR_S_CP  = 3.007;
constexpr double BAR_S_CN  = 3.329;
constexpr double BAR_S_MASS_NUM = 35.45;                      //mBarre = NUM / DEN * BAR_MASS_RESCALE
constexpr double BAR_S_MASS_DEN = 3.84723;
constexpr double BAR_E_RX0 = 4.44;
constexpr double BAR_E_RY0 = 1.31;
constexpr double BAR_E_RZ0 = 0.80;
constexpr double BAR_E_RC  = 6.83;
constexpr double BAR_E_CP  = 2.786;
constexpr double BAR_E_CN  = 3.917;
constexpr double BAR_E_MASS_NUM = 2.27;                       //mBarre = NUM / DEN * BAR_MASS_RESCALE
constexpr double BAR_E_MASS_DEN = 87.0;                       // normalized

// ---- (2b) Galactic model: kinematics (vrel, Lensing.cpp) ----
// Velocity dispersions [km/s] in the Galactic (R, T, Z) frame, drawn as truncated Gaussians.
constexpr double VEL_NSIGMA_TRUNC = 3.5;   //truncation of the Gaussian velocity draws [sigma]
// Thin disc: one row per age bin, chosen by the bin's share of Rho[0..3] = rho0*corr/d0.
constexpr std::array<double, 4> THIN_AGE     = {0.075, 0.575, 1.5, 2.5};   //Gyr
constexpr std::array<double, 4> THIN_SIGMA_R = {16.7, 19.8, 27.2, 30.2};
constexpr std::array<double, 4> THIN_SIGMA_T = {10.8, 12.8, 17.6, 19.5};
constexpr std::array<double, 4> THIN_SIGMA_Z = {6.0, 8.0, 10.0, 13.2};
constexpr double THICK_SIGMA_R = 67.0,  THICK_SIGMA_T = 51.0,  THICK_SIGMA_Z = 42.0;
constexpr double HALO_SIGMA_R  = 131.0, HALO_SIGMA_T  = 106.0, HALO_SIGMA_Z  = 85.0;
constexpr double BULGE_SIGMA_R = 113.0, BULGE_SIGMA_T = 115.0, BULGE_SIGMA_Z = 100.0;
// Rotation of the disc: V_T += vro_sun * (ROT_A * (R/Dsun)^ROT_SLOPE + ROT_B), thin and thick disc only.
constexpr double ROT_A = 1.00762;
constexpr double ROT_SLOPE = 0.0394;
constexpr double ROT_B = 0.00712;

// ==========================================================================================
// (3) FILTERS: LSST ugrizy + Roman F146
// ==========================================================================================

// NOT freely editable: indices 0-5 = LSST ugrizy and 6 = Roman F146 are hard-coded throughout the
// code, and every per-filter array below must have exactly M entries.
constexpr int    M = 6 + 1;    //No. of filter  ugrizy  of LSST + Roman's F146 filter
constexpr double delta2 = 0.005;///systematic errors


// gamma of the LSST photometric error model, sigma_rand^2 = (0.04 - gamma) x + gamma x^2 (Ivezic
// et al. 2019 eq. 5, Table 2; Deviation 73 -- the values were an older version, 0.037-0.040). Rubin
// only (index 6, F146, is unused). delta2 above is the 5 mmag bright-end repeatability (Ivezic
// et al. 2019, requirement 3). The unused seeing/msky/Cm/Dci/km constants (an older Table 2) and
// cade1 were deleted in Deviation 73.
constexpr std::array<double, M> gama = {0.038, 0.039, 0.039, 0.039, 0.039, 0.039, 0.0};
// The 0.04 of the same equation: sigma_rand^2 = (LSST_ERR_C04 - gamma) x + gamma x^2 (Ivezic et al. 2019
// eq. 5), used by errlsstM in helper.cpp.
constexpr double LSST_ERR_C04 = 0.04;

// Truncation of the per-epoch measurement-noise draws, RandN(err, n): a draw is rejected if it lies
// beyond n sigma (Bulge_LSST.cpp, photometric and astrometric noise of every Rubin and Roman epoch).
constexpr double NOISE_TRUNC_NSIGMA = 3.0;
// Truncation of the per-band extinction scatter A_i = Av*A_i/A_V + RandN(sigma[i], n) (Lensing.cpp).
constexpr double EXT_SCATTER_TRUNC_NSIGMA = 1.0;
// Per-epoch outlier flag and run test (Bulge_LSST.cpp, per instrument): an epoch is flagged when its noisy
// magnitude departs from the baseline by more than OUTLIER_FLAG_NSIGMA times its error; a run is declared
// (flag_det_*) when the sum of the last three epochs' flags exceeds OUTLIER_RUN_THRESHOLD, i.e. all three
// consecutive epochs are flagged (and more than 2 epochs have been taken).
constexpr double OUTLIER_FLAG_NSIGMA = 3.0;
constexpr double OUTLIER_RUN_THRESHOLD = 2.0;
// Light-curve time-loop step [days] while the epoch is outside the Rubin visit window (the loop then
// advances by min(this, Roman's cadence); inside the window it follows the local visit spacing).
constexpr double TIME_STEP_OUTSIDE_LSST_DAYS = 3.0;

// Rubin's saturation, relative to each visit's own 5-sigma depth (Deviation 73): saturation =
// fiveSigmaDepth - RUBIN_SATU_BELOW_M5. Ivezic et al. 2019 give "the LSST saturation limit at r ~ 16"
// for a 24.35 design depth (8.35 mag); the previous fixed constants were exactly depth - 8.3.
constexpr double RUBIN_SATU_BELOW_M5 = 8.3;


constexpr std::array<double, M> sigma = {0.022, 0.02, 0.017, 0.017, 0.027, 0.027, 0.04}; // PLACEHOLDER: K-band value, not F146
// thre / satu: single-visit 5-sigma depth and saturation, AB. ugrizy are the SRD minimum depths
// (Ivezic et al. 2019 Table 1, "min.") and depth - 8.3 (M1b replaces their use in the Rubin gate
// by each visit's own depth). F146 (Deviation 72): depth 25.45 = STScI's 5-sigma point-source
// sensitivity in 57 s, 25.37 AB (roman-technical-information, AB_mag_limiting_sensitivity.ecsv,
// 2x minimum zodi), scaled to the GBTDS's 66-s exposure, +1.25 log10(66/57); saturation 14.8 =
// Penny et al. 2019 Table 3 ("W149 saturation ~14.8", brightest pixel 1e5 e- before the first read,
// usable up-the-ramp). Both were placeholders (29.0, 12.0) until Deviation 72.
constexpr double ROMAN_DEPTH5_AB = 25.37 + 1.25 * 0.06368;   // log10(66/57) = 0.06368 -> 25.4496
constexpr double ROMAN_SATU_AB   = 14.8;
constexpr std::array<double, M> thre  = {23.4, 24.6, 24.3, 23.6, 22.9, 21.7, ROMAN_DEPTH5_AB};
constexpr std::array<double, M> satu  = {15.2, 16.3, 16.0, 15.3, 14.6, 13.4, ROMAN_SATU_AB};
// 5-sigma depth [AB] assumed for the r band at an epoch-less draw (errlsstM in the no-light-curve fallback of
// Bulge_LSST.cpp), where no visit supplies its own depth.
constexpr double RUBIN_R_DEPTH5_FALLBACK = 24.43;
// PSF FWHM [arcsec]: the image-resolution bar (Step R1) and the blending disc (Lensing.cpp).
// ugrizy (Deviation 73): the median GEOMETRIC seeing, OpSim seeingFwhmGeom (= 0.822 seeingFwhmEff +
// 0.052, verified exactly in the database), of the 12,308 bulge visits in Baseline/BulgeBaseline.dat
// (baseline_v5.1.0); per-visit 16-84% spans ~0.77-1.4". The previous values (1.221 ... 0.937) were an
// older OpSim's, 1-10% wider. F146: 0.105", STScI SummaryPSFstats (centre and corner).
constexpr std::array<double, M> FWHM  = {1.1140, 1.0420, 0.9819, 0.9487, 0.9320, 0.8992, 0.105};
// The blending disc has radius FWHM * BLEND_RADIUS_FWHM_FRAC (= the HWHM), Lensing.cpp func_source: the
// expected number of stars in it is lambda, and the source's own blend is 1 + Poisson(lambda).
constexpr double BLEND_RADIUS_FWHM_FRAC = 0.5;
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

// ==========================================================================================
// (4) RUBIN / LSST SURVEY AND FOOTPRINT
// ==========================================================================================

constexpr double Tobs = 10.0 * year;///LSST observational time 10 years
// ---------------------------------------------------------------------------------------
// Rubin/LSST astrometric error: renormalising files/sigmaA_LSST.txt from mission-averaged
// to PER VISIT.
//
// Sources:
//   [3] Ivezic et al. 2019, "LSST: From Science Drivers to Reference Design and Anticipated
//       Data Products", ApJ 873, 111 (arXiv:0805.2366): the survey's proper-motion and
//       parallax requirements are derived from "an assumed astrometric accuracy of 10 mas
//       per observation per coordinate", with ~1000 observations needed to reach the ~0.6-1
//       mas mission parallax accuracy.
//   [4] SITCOMTN-159 (Rubin commissioning, Operations Rehearsal 3): single-visit source
//       positions carry a systematic uncertainty of 3-7 mas to be added in quadrature to
//       the pipeline uncertainty.
//
// THE DEFECT. files/sigmaA_LSST.txt is a MISSION-AVERAGED curve, not a per-visit one, and
// errlsstA() was feeding it straight into l.erra[], which FisherM divides by per epoch and
// then sums over ~2,300 epochs -- applying the sqrt(N) averaging a second time. Two
// independent checks identify the factor:
//
//   * the table's bright-star floor is 0.3739576 mas; 10.0/0.3739576 = 26.74, i.e. exactly
//     [3]'s 10 mas per visit divided by sqrt(715), and ~715 visits is the right order for a
//     10-year all-band LSST count;
//   * scaled by that same 26.74 the faint end reads 132 mas at r = 24.44, against an
//     independent seeing-limited estimate FWHM/SNR ~ 700 mas / 5 ~ 140 mas.
//
// The shape of the table is therefore right and only its normalisation is wrong, so the fix
// renormalises rather than replaces: multiply by LSST_AST_FLOOR / LSST_AST_TABLE_FLOOR so
// the bright end sits on [3]'s published per-visit figure and the magnitude dependence
// already encoded in the file is preserved untouched.
//
// Uncorrected, this made Rubin's per-epoch astrometry 26.7x better than reality and its
// astrometric Fisher information ~715x too large -- which is why a ground-based telescope
// was out-centroiding Roman (0.374 mas against Roman's sourced 1.1 mas) before this.
// 10 mas is the conservative choice of the two sources: [4] suggests the delivered
// single-visit floor may be nearer 3-7 mas, which would make Rubin better than assumed here.
// ---------------------------------------------------------------------------------------
constexpr double LSST_AST_FLOOR       = 10.0;      //mas per visit per coordinate [3]
// LSST_AST_TABLE_FLOOR (0.3739576 mas, the bright-star floor of the table as shipped) is a property of the
// table, so it lives in the generated config/data_products.h.
constexpr double LSST_AST_RENORM      = LSST_AST_FLOOR / LSST_AST_TABLE_FLOOR; //26.74
constexpr double FoV = double(3.5 / 2.0);  //the radius of teh Rubin Field of View
// The scan region (Deviation 69; replaces the l1/l2/b1/b2 box and its lx/bx corner cut, which
// were built from an older field layout and cut the wrong corner). A sky point is scanned if a
// Rubin pointing that ALSO images a Roman field could image it: such a pointing is centred within
// FoV + rField of a Roman field centre, and images points within FoV of its own centre, so the
// region is every point within SCAN_RUBIN_REACH + rField of any field centre, spring or autumn.
// rField comes from the detector layout at run time (GbtdsLayout::rField, ~0.48 deg).
constexpr double SCAN_RUBIN_REACH = 2.0 * FoV;
constexpr double RUBIN_MAX_RADIUS = 1.94;   //deg: rubin_scheduler's max_radius; no active pixel beyond

// ==========================================================================================
// (5) ROMAN SURVEY AND FOOTPRINT
// ==========================================================================================

// Roman's footprint is the ADOPTED GBTDS layout (Step 1 of the 2026-10 pre-production fixes,
// Deviation 69): six fields, each a mosaic of 18 rectangular detectors (SCAs) with gaps between
// them, placed differently in spring and autumn because the telescope rolls by 180 deg between
// the two seasons. The field centres travel with every visit in RomanBaseline.dat (columns l, b,
// layout); the detector rectangles, as (l, b) offsets from the centre, are read at start-up from
// the vendored files below (Baseline/gbtds_layout/README.md: source, commit, checks). A sightline
// sees a Roman visit only if it falls ON a detector of that visit's layout -- not, as until
// Deviation 69, if it lies within an equal-area circle of 0.3003 deg about the field centre.
constexpr int    GBTDS_NLAYOUT = 2;    // 0 = spring roll, 1 = autumn roll
constexpr int    GBTDS_NSCA    = 18;   // detectors per field
inline const char* const GBTDS_SCA_FILES[GBTDS_NLAYOUT] = {
    "./Baseline/gbtds_layout/sca_layout_spring.txt",
    "./Baseline/gbtds_layout/sca_layout_fall.txt"};
// Clustering rule: consecutive distinct epoch times more than SEASON_GAP_MIN_DAYS apart
// begin a new season. This is safe by a wide margin on the current schedule -- the
// largest spacing INSIDE a season is 5.0 d (the low-cadence seasons' five-day sampling)
// and the smallest gap BETWEEN seasons is 108.2 d -- but the margin is checked at
// runtime rather than assumed; see the guard in main().
constexpr double SEASON_GAP_MIN_DAYS = 20.0;
// Sun-Earth L2, where Roman flies, as a fraction of an AU (Step H1).
//
// L2 is on the Sun-Earth line, ~1.5e6 km beyond the Earth, so to leading order Roman's
// heliocentric position is Earth's scaled by (1 + L2_OFFSET_AU). This is what makes the two
// observatories different places and gives the joint fit a *spatial* baseline to go with its
// temporal one: the difference in impact parameter the two see is
//     delta_u ~ L2_OFFSET_AU * piE  ~  1e-3  for a typical bulge event.
// Small, and concentrated in high-magnification, short-tE events -- PHASE_H_PLAN.md 0.3.
//
// Roman's halo orbit about L2 (amplitude ~1e5-1e6 km) is NOT modelled; this is the mean
// offset only. OPEN_ITEMS.md.
constexpr double L2_KM        = 1.5e6;
constexpr double L2_OFFSET_AU = L2_KM / AU_KM;   // ~0.01003
// Roman's photometric error (Deviation 72). files/sigma_roman.txt is Penny et al. 2019 Fig. 4: single-
// epoch precision vs W149 (= F146) AB magnitude for a 46.8-s Cycle-7 exposure, with a 1 mmag floor.
// At load time its photon-noise part, sqrt(err^2 - floor^2), is shifted in magnitude so that the
// 5-sigma point (err = 1.0857/5) falls at ROMAN_DEPTH5_AB, the current STScI depth for the GBTDS's
// 66-s exposure, and the floor is re-added; errRomanM then interpolates log(err) linearly in
// magnitude (was nearest neighbour). Penny's curve rescaled by exposure time alone would sit at
// ~25.71 (0.26 mag deeper than STScI's current figure), because it was made for an older design;
// anchoring to STScI's number moves it by only ~-0.07 mag.
constexpr double ROMAN_PHOT_FLOOR = 0.001;   //mag, Penny et al. 2019 Table 2 "Error floor 1.0 mmag"
// ---------------------------------------------------------------------------------------
// Roman WFI per-exposure astrometric precision, F146 (a.k.a. W149), in milliarcseconds.
// Step H4; used by errRomanA() in helper.cpp.
//
// Sources:
//   [1] Sanderson et al. 2019, "Astrometry with the Wide-Field Infrared Survey Telescope",
//       arXiv:1712.05420 sec 1.1: "single-exposure precision for well-exposed point sources
//       is 0.01 pixel, or about 1.1 mas", improving by ~10x when ~100 exposures are stacked.
//   [2] "Black hole astrometric binaries in the Roman Galactic Bulge Time Domain Survey",
//       arXiv:2608.24998, Fig. 5 and surrounding text: 1% centroiding => "a floor of 1.1 mas
//       for Roman"; the floor "impacts bright sources F146_Vega < 20.62 mag"; sources become
//       background dominated near "F146_Vega < 23.5 mag, which corresponds to sigma_ast ~ 10
//       mas"; pixels are "0.11 arcsec"; each GBTDS exposure is "66 seconds" at a "12.1
//       minute" cadence. Their curve derives from Pandeia and the Roman astrometry
//       simulation tool of Bellini et al. 2024.
//
// PER EXPOSURE, and this is a factor of ten waiting to be got wrong. The 0.1 mas figure that
// appears in both sources is the DAILY-BINNED precision -- ~100 exposures stacked. Our
// l.erra[] is a per-epoch error and one row of RomanBaseline.dat is one 12.1-minute exposure
// (measured: median inter-epoch gap 0.008403 d = 12.1 min, 50,401 epochs per field x 6
// fields = 302,406 = NlRoman). So 1.1 mas is the right floor here. Using 0.1 would overstate
// Roman's astrometry tenfold and flatter every tetE and lens-mass forecast in the project.
// ---------------------------------------------------------------------------------------
// MAGNITUDE SYSTEM (Deviation 72). [2]'s anchors (20.62, 23.5) are F146 VEGA magnitudes; the
// simulator's magnitudes are AB (the MIST bolometric-correction tables in CMD/ are "Roman (AB)").
// errRomanA therefore converts first: m_Vega = m_AB - F146_AB_MINUS_VEGA. The offset is synphot's
// AB magnitude of Vega (CALSPEC alpha_lyr_stis_011) through STScI's F146 effective area
// (roman-technical-information Roman_effarea_v8_SCA01_20240301): 1.0324 mag (EXOZIPPy #313 gets
// 1.037; it lies between STScI-000825's 2MASS J 0.913 and H 1.391, as F146 spans both). Before
// Deviation 72 the AB magnitude was used as if Vega: every Roman astrometric error ~2.2-2.6x too big.
constexpr double F146_AB_MINUS_VEGA = 1.0324;
constexpr double ROMAN_PIX_MAS   = 110.0;  //0.11 arcsec pixels [2]
constexpr double ROMAN_AST_FLOOR = 0.01 * ROMAN_PIX_MAS;  //1.1 mas: 1% centroiding [1][2]
constexpr double ROMAN_AST_MFLR  = 20.62;  //mag below which the floor dominates [2]
constexpr double ROMAN_AST_MBKG  = 23.5;   //mag where the background starts to dominate [2]
constexpr double ROMAN_AST_SBKG  = 10.0;   //mas, sigma_ast at ROMAN_AST_MBKG [2]
// Slope between the two anchors above. NOT a physical constant and NOT a free choice:
// log10(10.0/1.1)/(23.5-20.62) = 0.3329 per mag. Source-dominated photon noise
// (SNR ~ sqrt(counts)) would give 0.2/mag and pure background domination (SNR ~ counts)
// gives 0.4/mag; 0.333 sits between them because the transition is already under way across
// this range. A later step could replace this interpolation with a Pandeia-derived table,
// exactly as errlsstA reads files/sigmaA_LSST.txt.
constexpr double ROMAN_AST_SLOPE_SRC = 0.33285;
// ---------------------------------------------------------------------------------------
// Step 3c (Deviation 71). The astrometric noise model, three ways.
//
// PER COORDINATE. errRomanA / errlsstA give the 1D (x or y) per-exposure precision, as every
// source defines it (Lam et al. 2026 fn. 14; McKinnon & van der Marel 2026; Ivezic et al. "per
// observation per coordinate"). Each coordinate's variance is erra^2 -- not 2 erra^2, which the
// code used until Deviation 71 (every astrometric sigma was sqrt(2) too large).
//
// THE REFERENCE POSITION IS FREE. Real astrometric fits solve for the source's position offset;
// the model's -u0 tetE sin(xi) term otherwise lets every exposure measure tetE from the source's
// absolute position (fixture: sigma(tetE) x1.4-10 optimistic). An offset per "frame group" is
// marginalised in closed form: F_g = sum_blocks F_k - (sum b_k)(sum b_k)^T / sum c_k, per coordinate.
//
// WHETHER THE 1.1 mas FLOOR AVERAGES DOWN is unknown: the literature adds it as white noise
// (justified by the GBTDS's sub-pixel dithers) and quantifies no correlated part (Deviation 71).
// So every event carries three forecasts:
//   W  white (the literature's assumption): one free offset per telescope; Roman's errors white.
//   N  nominal: one free offset per Roman ROLL (crowding biases flip with the PSF orientation) and
//      a per-coordinate error AST_SIGC_N shared by all Roman exposures of the same day (distortion
//      residuals at "a few x 0.1% of a pixel", Bellini 2024 via Lam et al. 2026).
//   P  pessimistic: one free offset per Roman SEASON (each season its own frame) and the WHOLE
//      floor, AST_SIGC_P = 1.1 mas, shared within each day (it averages only across days).
// Rubin's errors are white with one offset in all three. The day blocks enter by Sherman-Morrison:
// for a block with weights w_i = 1/erra_i^2 and derivatives d_i, F_k = S_wdd - s^2 S_wd S_wd^T /
// (1 + s^2 S_w), b_k = S_wd / (1 + s^2 S_w), c_k = S_w / (1 + s^2 S_w), s = sigma_c.
// The main table columns (sigtetE_*, relMl_*, okB_*, condB_*) are W.
// ---------------------------------------------------------------------------------------
constexpr double AST_SIGC_N = 0.3;               //mas per coordinate per Roman day
constexpr double AST_SIGC_P = ROMAN_AST_FLOOR;   //1.1 mas
constexpr std::array<double, NAVAR> AST_SIGC = {0.0, AST_SIGC_N, AST_SIGC_P};
constexpr int    AST_MAX_SEASONS = 16;           //Roman seasons a P-variant offset can be keyed on

// ==========================================================================================
// (6) LENS POPULATIONS AND MASS FUNCTIONS
// ==========================================================================================

// The bulge entry is first and is the default: it reproduces the pre-2026-09-17 behaviour
// exactly, tag included.
inline constexpr LensPopulation POPULATIONS[] = {
    {"bulge", "5", MassFunction::KROUPA_REMNANTS, 0.01, 30.0, false, 5,
     "Kroupa IMF + remnants; the present-day bulge population"},
    {"bh",    "bh", MassFunction::LOG_UNIFORM,     3.0, 1000.0, true, 0,
     "black holes, flat in log M over 3-1000 Msun"},
    {"ns",    "ns", MassFunction::NEUTRON_STAR,    1.0,    2.5, false, 0,
     "neutron stars, Gaussian about 1.35 Msun (Ozel & Freire 2016)"},
    {"macho-uniform", "1", MassFunction::UNIFORM,      3.0, 5000.0, true, 1, "legacy MACHO search"},
    {"macho-m05",     "2", MassFunction::POWER_LAW_05, 3.0, 5000.0, true, 2, "legacy MACHO search"},
    {"macho-m1",      "3", MassFunction::POWER_LAW_10, 3.0, 5000.0, true, 3, "legacy MACHO search"},
    {"macho-m2",      "4", MassFunction::POWER_LAW_20, 3.0, 5000.0, true, 4, "legacy MACHO search"},
};
constexpr double u0m   = 3.0;
// Event-parameter draws (Lensing.cpp, func_lens): u0 is uniform in [U0_MIN_DRAW, u0m]; the peak time t0 is
// uniform in [T0_MARGIN_DAYS, Tobs - T0_MARGIN_DAYS] days, i.e. at least 2 days inside the survey ends.
constexpr double U0_MIN_DRAW = 0.001;
constexpr double T0_MARGIN_DAYS = 2.0;
// ---- Kroupa (2001) initial mass function, dN/dM ~ M^-alpha, with the standard breaks ----
// Coefficients enforcing continuity at the breaks are derived in drawKroupaInitialMass();
// only the breaks and slopes are named here.
constexpr double KROUPA_MI_MIN = 0.01;   //below the hydrogen-burning limit: brown dwarfs
constexpr double KROUPA_MI_MAX = 120.0;  //initial mass; nothing this heavy survives to today
constexpr double KROUPA_BREAK1 = 0.08;   //hydrogen-burning limit
constexpr double KROUPA_BREAK2 = 0.50;
constexpr double KROUPA_ALPHA1 = 0.3;    //0.01 - 0.08
constexpr double KROUPA_ALPHA2 = 1.3;    //0.08 - 0.50
constexpr double KROUPA_ALPHA3 = 2.3;    //0.50 - 120   (Salpeter-like)

// ---- Initial-to-final mass, for the remnants ----
// A bulge population is ~10 Gyr old, so everything born above the turnoff is already dead.
// Which remnant it left is set by its INITIAL mass, and that is what makes the long-tE tail
// this project cares about: a black hole lens is heavy, so tE ~ sqrt(Ml) is long, and a long
// event is exactly the one that spans Roman's season gaps.
// ---- Neutron star masses, for the NEUTRON_STAR population ----
// The observed distribution is narrow and well measured. Ozel & Freire (2016) review the
// measured sample: double neutron stars cluster at 1.33 +/- 0.09 Msun, slow pulsars sit
// slightly higher and recycled ones spread wider, and the population as a whole is often
// summarised as a Gaussian near 1.35 Msun with a ~0.15 Msun spread. The truncation is
// physics, not tidiness: below ~1.1 Msun no supernova is known to leave a neutron star, and
// above ~2.2 Msun the equation of state gives a black hole instead.
//
// UNCERTAINTY, STATED: the real distribution is arguably bimodal (a recycled population
// above the canonical peak), and a single Gaussian will understate the high-mass tail. The
// alternative is one line in drawNeutronStarMass(); this is the simpler model, chosen
// deliberately and recorded rather than presented as settled.
constexpr double NS_MEAN_MASS = 1.35;  //Msun
constexpr double NS_MASS_SIG  = 0.15;  //Msun
constexpr double NS_MASS_TRUNC_NSIGMA = 4.0;   //truncation of the Gaussian NS-mass draw [sigma]
constexpr double NS_MASS_LO   = 1.10;  //Msun, below which no NS is expected to form
constexpr double NS_MASS_HI   = 2.20;  //Msun, near the maximum the equation of state allows

constexpr double MS_TURNOFF   = 1.0;   //Msun; below this the star is still on the main sequence
constexpr double WD_MI_MAX    = 8.0;   //Mi 1-8   -> white dwarf
constexpr double NS_MI_MAX    = 20.0;  //Mi 8-20  -> neutron star; above -> black hole
constexpr double NS_MASS      = 1.4;   //Msun, the canonical value
constexpr double BH_MI_SLOPE  = 0.24;  //Ml = 0.24*Mi, giving ~4.8-28.8 Msun over Mi 20-120

// ==========================================================================================
// (7) DETECTION AND CHARACTERISATION
// ==========================================================================================

// ---------------------------------------------------------------------------------------
// Step R1. Resolving the two lensing-induced images.
//
// A point lens makes TWO images of the source, at
//     theta_pm = 0.5 * (u +- sqrt(u^2+4)) * theta_E,
// so their angular separation is
//     Delta_theta(u) = theta_E * sqrt(u^2 + 4).
// Ordinarily we never see them apart: for a bulge event with a stellar lens theta_E ~ 0.3 mas
// and 2*theta_E is far below any of our resolutions, so what reaches the detector is the sum
// of the two fluxes -- the magnification A -- and their flux-weighted centroid, which is the
// astrometric shift of Step H5. A BLACK-HOLE lens is the interesting case, because
// theta_E scales as sqrt(Ml): at 1000 Msun it is ~30x the stellar value.
//
// TWO THINGS FIGHT EACH OTHER, which is the whole reason this needs counting per epoch rather
// than once per event. The separation is SMALLEST at closest approach, sqrt(u0^2+4)*theta_E,
// and grows without bound as the source moves away. But the minor image's magnification,
//     A_minus = (u^2+2) / (2 u sqrt(u^2+4)) - 1/2,
// collapses toward zero on the same motion. So the images are least separated exactly when
// both are bright, and are well separated only once one of them has faded. A criterion applied
// at peak would be far too optimistic, and one applied at maximum separation would be
// meaningless. Sajadian & Makler (arXiv:2608.16448, their sec. 3) resolve this by counting
// DATA POINTS that satisfy both conditions at once, and calling the images resolvable if at
// least three do. That count cannot be reconstructed from a per-event summary row, which is
// why it is accumulated here inside the epoch loop.
//
// THE RESOLUTION THRESHOLD. That paper takes the Rubin criterion as Delta_theta >= sigma_r
// with sigma_r = D * sigma_a, where sigma_a is the per-visit astrometric precision and D runs
// over [5, 100] with the source's signal-to-noise: ~5 at the faint limit (SNR 5), ~20 at
// SNR 100 for two stars of similar brightness, and higher still -- by 40% for a 2-3 mag
// difference, 3x beyond that -- when one image is much fainter than the other (Ivezic, priv.
// comm. quoted therein). We record the count at BOTH anchors rather than picking one, because
// the answer depends strongly on D and a single number would hide that.
//
// We add a third, independent criterion the paper does not need: Delta_theta >= the PSF FWHM.
// Their target is Rubin alone, where the empirical D*sigma_a captures seeing-limited reality.
// Roman is diffraction limited at 0.105 arcsec in F146, and for a space telescope the PSF
// width is the honest physical bar -- D*sigma_a on Roman's 1.1 mas floor would claim a
// resolution of a few mas, which no 2.4 m telescope delivers. Quoting all three makes the
// assumption visible instead of buried.
// ---------------------------------------------------------------------------------------
constexpr double RESOLVE_D_FAINT  = 5.0;   //D at the faint detection limit, SNR ~ 5
constexpr double RESOLVE_D_BRIGHT = 20.0;  //D at SNR ~ 100, images of similar brightness
// ---------------------------------------------------------------------------------------
// Step H7. The detection threshold.
//
// A microlensing detection is declared when the lensing model beats a flat-baseline model by
// enough chi-squared. The statistic is
//
//     dchi = chi2(flat baseline) - chi2(lensing model)
//
// which is a chi-squared DIFFERENCE between nested models. Under the null hypothesis of no
// lensing it is distributed as chi-squared with p degrees of freedom, where p is the number of
// extra parameters the lensing model carries -- NOT the number of epochs. Its expectation is p
// and its variance 2p, both fixed and small. The bar therefore does not scale with the epoch
// count, and a detection threshold that does scale with it is thresholding the MEAN per-epoch
// improvement rather than the total significance.
//
// This is the whole content of Step H7. The previous form, dchi > 2*ndw, made the bar grow
// with the number of epochs, so pooling a survey with many low-signal epochs raised the joint
// bar without adding signal: in Roman's footprint Roman's own bar was 2*50401 = 100,802 while
// the joint bar was 105,530, and Rubin's ~2,300 near-flat epochs lifted it by ~4,728. An event
// clearing Roman's bar by less than that failed the JOINT test -- measured at 21.3% of all
// detections inside the footprint on the 2026-09-05 v2 run.
//
// A fixed bar also restores monotonicity by construction, which is the property that actually
// matters: chi1 and chi3 are accumulated over both instruments' epochs, so chi1 = chi1_L +
// chi1_R and chi3 = chi3_L + chi3_R exactly, hence dchi = dchi_L + dchi_R. With the same
// threshold on all three tests, either survey clearing the bar alone forces the joint sum over
// it too, so detL or detR implies detJ and DET_ANOMALY cannot occur. (This requires the SIGNED
// difference; see the note on fabs at the test site in Bulge_LSST.cpp.)
//
// The value 500 is Penny et al. 2019 (ApJS 241, 3), the reference Roman/WFIRST microlensing
// yield forecast, which adopts dchi2 > 500 against a flat baseline. Matching it keeps this
// project's yields comparable with the number the Roman community already quotes. It is
// deliberately conservative -- a nominal 3-sigma bar on a few parameters would be nearer 20 --
// because the real false-alarm population is systematics, variable stars and blending, not
// Gaussian noise, and a high bar is the standard defence. Overridable with --dchi-det for
// sensitivity tests, and recorded in run_provenance.txt because every yield in this project
// is conditioned on it.
// ---------------------------------------------------------------------------------------
constexpr double DCHI_DET_DEFAULT = 500.0; //delta-chi2 against a flat baseline [Penny+2019]
// Step-size multipliers on the astrometric Delta2[] (tetE, mus1, mus2, piE), chosen by the M3 sweep
// (Deviation 75; ./fishertest --sweep-astro). 1.0 = the legacy 25%-of-value steps. tetE, mus1, mus2:
// the modelled centroid is LINEAR in them (blending and lens light included), so every finite
// difference is exact -- the sweep is flat to all digits over 1e-8..2. piE (through the parallax-
// bent trajectory) is not: plateau 1e-8..1e-2, 0.4% off at the legacy step, 2.5% at 2x; 1e-2 chosen.
constexpr std::array<double, 4> kFDStepScaleB = {1.0, 1.0, 1.0, 1.0e-2};

// Photometric finite-difference steps of FisherM (Delta1[], Bulge_LSST.cpp). u0, tE, piE, xi, t0 came from
// the legacy LMC code and are ~25% of the parameter (u0 is perturbed by 0.15, tE and t0 by 0.25*tE). Step
// C3 (DEVIATIONS.md 10; ./fishertest --sweep, tests/c3_step_sweep.py) showed that sat far up the
// truncation-error branch (sigma(u0) off ~57%, sigma(t0) ~71%), and scaled them all by kFDStepScale into
// the convergence plateau (flat to <0.2% over 1e-6..1e-3 of the legacy steps; 1e-4 sits two decades clear
// of the round-off wall). The legacy values are kept and the scale factored out so the change stays
// auditable and the sweep, defined in these units, stays comparable. Re-run the sweep after changing any.
constexpr double FD_STEP_U0       = 0.1507586576;      //u0 (absolute)
constexpr double FD_STEP_TE_FRAC  = 0.254674;          //tE and t0: step = tE * this [days]
constexpr double FD_STEP_PIE_FRAC = 0.2509463534656;   //piE: step = this * piE
constexpr double FD_STEP_XI_DEG   = 3.0;               //xi [deg]; converted to radians at the use
constexpr double FD_STEP_MBS      = 0.05;              //mbs0 / mbs1 [mag]; the model is linear in mbs, so exact
constexpr double kFDStepScale     = 1.0e-4;            //Step C3 plateau scale applied to every Delta1[] and fb step
// fb0 / fb1 steps (blend fraction, bounded to [0,1]). The step pair is chosen from the epoch's own
// telescope's fb so that fb + step never leaves the range: fb < FB_BIN_LO -> {+SMALL, +LARGE};
// fb < FB_BIN_HI -> {-SMALL, +SMALL}; else {-SMALL, -LARGE}. Production behaviour (then * kFDStepScale).
constexpr double FB_BIN_LO     = 0.15;
constexpr double FB_BIN_HI     = 0.85;
constexpr double FB_STEP_SMALL = 0.07;
constexpr double FB_STEP_LARGE = 0.15;
// A telescope contributing fewer than kMinTeleEpochs epochs to an event is left out of the
// photometric matrices altogether (Deviation 76): its flux pair (fb, mbs) cannot be constrained,
// and with one epoch the joint matrix went singular (condition 1.5e16 on the bulge event whose joint
// mass error was 8.3x Roman's -- one Rubin epoch). Its epochs carry next to no information on the
// event anyway. The counts passed to activePhotParams are therefore >= kMinTeleEpochs or zero.
constexpr int kMinTeleEpochs = 3;

// ==========================================================================================
// (8) OUTPUT HISTOGRAM GRIDS
// ==========================================================================================

constexpr int GG = 100;
constexpr double tE_min  = 0.0;///days
constexpr double tE_max  = 50.0*year;//days
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

// ==========================================================================================
// (9) DATA-FILE CONVENTIONS
// ==========================================================================================

// Magnitude of a "dark" list entry (brown dwarf, white dwarf: no MIST track) -- zero light.
constexpr double DARK_MAG = 99.0;

// ==========================================================================================
// (10) DATA FILES
// ==========================================================================================

// Paths are relative to the repository root: ./roman must be run from there.
inline constexpr const char* PATH_BULGE_BASELINE = "./Baseline/BulgeBaseline.dat";
inline constexpr const char* PATH_ROMAN_BASELINE = "./Baseline/RomanBaseline.dat";
inline constexpr const char* PATH_LSSTCAM_FOV    = "./Baseline/lsstcam_fov/fov_map.txt";
inline constexpr const char* PATH_SIGMA_A_LSST   = "./files/sigmaA_LSST.txt";
inline constexpr const char* PATH_SIGMA_ROMAN    = "./files/sigma_roman.txt";
inline constexpr const char* PATH_EXT_TABLES     = "./files/ext/ext_tables.dat";
inline constexpr const char* PATH_LENS_ML        = "./CMD/components/lens_ml.dat";
inline constexpr const char* PATH_CMD_THIN       = "./CMD/components/thin_disk.dat";
inline constexpr const char* PATH_CMD_BULGE      = "./CMD/components/bulge.dat";
inline constexpr const char* PATH_CMD_THICK      = "./CMD/components/thick_disk.dat";
inline constexpr const char* PATH_CMD_HALO       = "./CMD/components/halo.dat";
// Output directory (LpLMC / EfLMC / MapLMC / run_provenance) and the Disk_model debug dumps;
// both end in a slash.
inline constexpr const char* PATH_OUT_DIR        = "./files/MONTLMC/files/";
inline constexpr const char* PATH_DENSITY_DIR    = "./files/density/";

#endif // PARAMETERS_H
