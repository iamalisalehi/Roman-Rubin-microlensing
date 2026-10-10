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
#include "data_products.h"   // row counts and catalogue means (generated from the data files)

// ==========================================================================================
// (1) RUN DEFAULTS
// ==========================================================================================

constexpr int seed = 42;   // default base seed; --seed overrides it at run time
// Per-sightline Monte Carlo budgets (defaults of --events / --lenses / --nerr / --maxdraws). The star-draw
// loop stops once events, lenses and nerr are ALL met, or at the draw cap.
constexpr int    DEFAULT_EVENTS_TARGET = 850;    //detected events per sightline (icon)
constexpr int    DEFAULT_LENSES_TARGET = 150;    //Fisher-characterised events per sightline (nlens)
constexpr double DEFAULT_NERR_TARGET   = 2.0;    //accumulated Fisher-error weight
constexpr double DEFAULT_MAXDRAWS      = 5.0e4;  //hard cap on stars drawn at one sightline
// Dense model-curve sampling of the sample light curves (defaults of the sample spec keys step_te /
// span_te / dt_coarse): fine steps of DEFAULT_S2_STEP_TE * tE within +-DEFAULT_S2_SPAN_TE * tE of the
// peak, steps of DEFAULT_S2_DT_COARSE days over the rest of the mission.
constexpr double DEFAULT_S2_STEP_TE   = 0.01;   //peak-window step, in units of tE
constexpr double DEFAULT_S2_SPAN_TE   = 3.0;    //peak-window half-width, in units of tE
constexpr double DEFAULT_S2_DT_COARSE = 2.0;    //step over the rest of the mission [days]

// ==========================================================================================
// (2) GALACTIC MODEL
// ==========================================================================================

constexpr double binary_fraction = double(2.0 / 3.0);
constexpr double vro_sun = 226.0;
constexpr double VSunR = 11.1;
// 1.00762 and 0.00712 are ROT_A and ROT_B (section 2b) evaluated at R = Dsun, kept as literals.
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
// Source distance draw: the distance-grid index is uniform in [SRC_IDX_MIN, Num - SRC_IDX_END_MARGIN]
// (func_source), skipping the first 5 grid cells (~0 kpc) and the last 2.
constexpr double SRC_IDX_MIN = 5.0;
constexpr double SRC_IDX_END_MARGIN = 2.0;
constexpr double dd  = 0.02;   // native sightline grid step [deg]; the scan steps by stride*dd

// ---- (2a) Galactic model: density laws (Disk_model, src/galaxy/density.cpp) ----
// Mass fractions of each component are set by the *_NORM factors below; number densities follow
// from Disk_model's rho / <m> with the MEANMASS_* of config/data_products.h.

// Star-count completeness factors per component.
constexpr double DENS_FD = 1.0;  //disc: all stars counted (no V < 11.5 limit)
constexpr double DENS_FB = 1.0;  //bulge: all stars counted
constexpr double DENS_FH = 1.0;  //halo: no limitation

// Thin disc (8 age bins; rho0, d0, epci, corr above are its per-bin parameters).
constexpr double THIN_RDD = 2.17;       //scale length of the old bins (ii > 0) [kpc]
constexpr double THIN_RHH = 1.33;       //hole scale length of the old bins (ii > 0) [kpc]
constexpr double THIN_YOUNG_L1 = 25.0;  //young (ii == 0) bin: outer scale, exp(-rdi/25)
constexpr double THIN_YOUNG_L2 = 9.0;   //young (ii == 0) bin: inner scale, exp(-rdi/9)
constexpr double THIN_CORE = 0.25;      //softening inside sqrt(0.25 + rdi/R^2) of the old bins
constexpr double THIN_NORM = 1.2;       //total mass 4.25e10 Msun

// Thick disc.
constexpr double THICK_RHO00 = 1.34 * 0.001 + 3.04 * 0.0001;  //local density, Msun/pc^3
constexpr double THICK_RHO_DIV = 0.999719;                    //normalisation divisor of rho00
constexpr double THICK_SCALE_LEN = 2.5;                       //radial scale length [kpc]
constexpr double THICK_H1 = 0.4;                              //height of the parabolic core [kpc]
constexpr double THICK_H2 = 0.8;                              //exponential scale height [kpc]
constexpr double THICK_NNF = 0.4 / 0.8;                       //THICK_H1 / THICK_H2
constexpr double THICK_NORM = 2.67;                           //total mass 0.8e10 Msun

// Stellar halo.
constexpr double HALO_FLATTEN = 0.76;                         //axis ratio
constexpr double HALO_CORE = 0.5;                             //core radius [kpc]
constexpr double HALO_RHO0 = (0.932 * 0.00001 / 867.067);     //local density, Msun/pc^3
constexpr double HALO_SLOPE = -2.44;                          //power-law index
constexpr double HALO_NORM = 5281.0;                          //total mass 1.2e9 Msun

// Bulge / bar: two triaxial components, S (boxy, 1/cosh^2 profile) and E (exponential).
constexpr double BAR_ANGLE_DEG = 12.89;                       //bar angle [deg]
constexpr double BAR_MASS_RESCALE = 0.24529; // calibrated so bulge column density toward
                                             // Baade's Window (l=1, b=-3.9) matches the
                                             // Han & Gould (2003) HST benchmark: 2086 Msun/pc^2
constexpr double BAR_CUTOFF_K = 4.0;                          //Gaussian cut-off exp(-K (r2-Rc)^2) beyond Rc
constexpr double BAR_NORM = 0.45;                             //total mass 1.7e10 Msun
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

// ---- (2b) Galactic model: kinematics (vrel, src/galaxy/kinematics.cpp) ----
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
constexpr int    M = 6 + 1;    //No. of filters: LSST ugrizy + Roman's F146
constexpr double delta2 = 0.005;//systematic error: 5 mmag bright-end repeatability (Ivezic et al. 2019, requirement 3)

// gamma of the LSST photometric error model, sigma_rand^2 = (0.04 - gamma) x + gamma x^2 (Ivezic
// et al. 2019 eq. 5, Table 2). Rubin only (index 6, F146, is unused).
constexpr std::array<double, M> gama = {0.038, 0.039, 0.039, 0.039, 0.039, 0.039, 0.0};
// The 0.04 of the same equation (used by errlsstM).
constexpr double LSST_ERR_C04 = 0.04;

// Truncation of the per-epoch measurement-noise draws (photometric and astrometric, Rubin and Roman):
// a draw beyond n sigma is rejected.
constexpr double NOISE_TRUNC_NSIGMA = 3.0;
// Truncation of the per-band extinction scatter A_i = Av*A_i/A_V + RandN(sigma[i], n).
constexpr double EXT_SCATTER_TRUNC_NSIGMA = 1.0;
// Per-epoch outlier flag and run test, per instrument: an epoch is flagged when its noisy magnitude
// departs from the baseline by more than OUTLIER_FLAG_NSIGMA times its error; a run is declared
// (flag_det_*) when the last three epochs' flags sum above OUTLIER_RUN_THRESHOLD, i.e. all three
// consecutive epochs are flagged (and more than 2 epochs have been taken).
constexpr double OUTLIER_FLAG_NSIGMA = 3.0;
constexpr double OUTLIER_RUN_THRESHOLD = 2.0;
// Light-curve time-loop step [days] outside the Rubin visit window (the loop advances by
// min(this, Roman's cadence); inside the window it follows the local visit spacing).
constexpr double TIME_STEP_OUTSIDE_LSST_DAYS = 3.0;

// Rubin's saturation relative to each visit's own 5-sigma depth: saturation = fiveSigmaDepth -
// RUBIN_SATU_BELOW_M5. Ivezic et al. 2019 give a saturation limit at r ~ 16 for a 24.35 design depth.
constexpr double RUBIN_SATU_BELOW_M5 = 8.3;


constexpr std::array<double, M> sigma = {0.022, 0.02, 0.017, 0.017, 0.027, 0.027, 0.04}; // PLACEHOLDER: K-band value, not F146
// thre / satu: single-visit 5-sigma depth and saturation, AB. ugrizy are the SRD minimum depths
// (Ivezic et al. 2019 Table 1, "min.") and depth - 8.3; the Rubin detection gate uses each visit's
// own depth instead. F146: depth 25.45 = STScI's 5-sigma point-source sensitivity in 57 s, 25.37 AB
// (roman-technical-information, AB_mag_limiting_sensitivity.ecsv, 2x minimum zodi), scaled to the
// GBTDS's 66-s exposure, +1.25 log10(66/57); saturation 14.8 = Penny et al. 2019 Table 3
// ("W149 saturation ~14.8", brightest pixel 1e5 e- before the first read, usable up-the-ramp).
constexpr double ROMAN_DEPTH5_AB = 25.37 + 1.25 * 0.06368;   // log10(66/57) = 0.06368 -> 25.4496
constexpr double ROMAN_SATU_AB   = 14.8;
constexpr std::array<double, M> thre  = {23.4, 24.6, 24.3, 23.6, 22.9, 21.7, ROMAN_DEPTH5_AB};
constexpr std::array<double, M> satu  = {15.2, 16.3, 16.0, 15.3, 14.6, 13.4, ROMAN_SATU_AB};
// 5-sigma depth [AB] assumed for the r band at an epoch-less draw, where no visit supplies its own depth.
constexpr double RUBIN_R_DEPTH5_FALLBACK = 24.43;
// PSF FWHM [arcsec]: the image-resolution bar and the blending disc. ugrizy: median geometric seeing,
// OpSim seeingFwhmGeom (= 0.822 seeingFwhmEff + 0.052), of the 12,308 bulge visits in
// Baseline/BulgeBaseline.dat (baseline_v5.1.0); per-visit 16-84% spans ~0.77-1.4". F146: 0.105",
// STScI SummaryPSFstats (centre and corner).
constexpr std::array<double, M> FWHM  = {1.1140, 1.0420, 0.9819, 0.9487, 0.9320, 0.8992, 0.105};
// The blending disc has radius FWHM * BLEND_RADIUS_FWHM_FRAC (= the HWHM) in func_source: the expected
// number of stars in it is lambda, and the source's own blend is 1 + Poisson(lambda).
constexpr double BLEND_RADIUS_FWHM_FRAC = 0.5;
constexpr std::array<double, M> lambda_um = {0.367, 0.482, 0.622, 0.755, 0.869, 0.971, 1.464};

// Which LSST filter(s) (indices 0-5 = u,g,r,i,z,y) form the single "Rubin representative band"
// standing in for the whole LSST light curve in the Fisher matrix, the recorded per-epoch Rubin model
// magnitude, and the LSST astrometric-error evaluation magnitude. {2} = r only. Listing several sums
// their baseline and source fluxes (equal-weighted, not throughput-weighted). It does not change
// which filter's noise model (errlsstM) applies to an epoch; that is always the epoch's own filter.
inline const std::vector<int> RUBIN_REF_BANDS = {2};

// ==========================================================================================
// (4) RUBIN / LSST SURVEY AND FOOTPRINT
// ==========================================================================================

constexpr double Tobs = 10.0 * year;//LSST observational time: 10 years
// ---------------------------------------------------------------------------------------
// Rubin/LSST per-visit astrometric error, per coordinate (errlsstA), in milliarcseconds:
//     sigma = sqrt( (LSST_AST_KAPPA * FWHM_visit / SNR_visit)^2 + LSST_AST_FLOOR^2 )
// The first term is the statistical limit of a PSF fit, sigma = kappa FWHM / SNR [5]. FWHM_visit is the
// visit's geometric PSF FWHM (OpSim seeingFwhmGeom = LSST_FWHM_GEOM_SLOPE * seeingFwhmEff +
// LSST_FWHM_GEOM_OFFSET [6], from the visit list's seeing column); SNR_visit = 1.0857 / sigma_rand, where
// sigma_rand is the photon-noise part of the same visit's photometric error (errlsstM, with the visit's
// own band and depth; the 5 mmag calibration term delta2 excluded). kappa depends on the PSF shape and on
// whether the source or the sky dominates the noise: 0.425 / 0.60 for a Gaussian, 0.51 / 0.63 for a
// long-exposure Kolmogorov seeing profile (source / background limited; [5]'s eqs. 11 and 14 evaluated
// on his eq. 43 profile). Bulge sources are background limited wherever this term is not buried under
// the floor, so 0.63.
// The floor is the per-visit systematic: atmospheric image motion that no analysis of the image can
// remove [5 sec. 6.4], white between visits. 10 mas is the conservative choice [3]; [4] suggests 3-7.
//
// Sources:
//   [3] Ivezic et al. 2019, ApJ 873, 111 (arXiv:0805.2366): requirements derived from "an assumed
//       astrometric accuracy of 10 mas per observation per coordinate".
//   [4] SITCOMTN-159 (Rubin commissioning, Operations Rehearsal 3): single-visit positions carry a
//       3-7 mas systematic, added in quadrature to the pipeline uncertainty.
//   [5] Lindegren 1978, IAU Coll. 48, p. 197 (Zenodo 10493823): eqs. 6-14, Table 1; Fritz et al. 2010,
//       MNRAS 401, 1177, eq. 2.
//   [6] the OpSim relation, also used for FWHM[] in section 3.
// ---------------------------------------------------------------------------------------
constexpr double LSST_AST_FLOOR        = 10.0;   //mas per visit per coordinate [3]
constexpr double LSST_AST_KAPPA        = 0.63;   //sigma * SNR / FWHM, Kolmogorov PSF, background limited [5]
constexpr double LSST_FWHM_GEOM_SLOPE  = 0.822;  //seeingFwhmGeom = SLOPE * seeingFwhmEff + OFFSET [arcsec] [6]
constexpr double LSST_FWHM_GEOM_OFFSET = 0.052;
constexpr double FoV = double(3.5 / 2.0);  //radius of the Rubin field of view [deg]
// The scan region. A sky point is scanned if a Rubin pointing that ALSO images a Roman field could
// image it: such a pointing is centred within FoV + rField of a Roman field centre and images points
// within FoV of its own centre, so the region is every point within SCAN_RUBIN_REACH + rField of any
// field centre, spring or autumn. rField comes from the detector layout at run time
// (GbtdsLayout::rField, ~0.48 deg).
constexpr double SCAN_RUBIN_REACH = 2.0 * FoV;
constexpr double RUBIN_MAX_RADIUS = 1.94;   //deg: rubin_scheduler's max_radius; no active pixel beyond

// ==========================================================================================
// (5) ROMAN SURVEY AND FOOTPRINT
// ==========================================================================================

// Roman's footprint is the adopted GBTDS layout: six fields, each a mosaic of 18 rectangular
// detectors (SCAs) with gaps between them, placed differently in spring and autumn because the
// telescope rolls by 180 deg between the two seasons. The field centres travel with every visit in
// RomanBaseline.dat (columns l, b, layout); the detector rectangles, as (l, b) offsets from the
// centre, are read at start-up from the vendored files below (Baseline/gbtds_layout/README.md gives
// source and checks). A sightline sees a Roman visit only if it falls ON a detector of that layout.
constexpr int    GBTDS_NLAYOUT = 2;    // 0 = spring roll, 1 = autumn roll
constexpr int    GBTDS_NSCA    = 18;   // detectors per field
inline const char* const GBTDS_SCA_FILES[GBTDS_NLAYOUT] = {
    "./Baseline/gbtds_layout/sca_layout_spring.txt",
    "./Baseline/gbtds_layout/sca_layout_fall.txt"};
// Clustering rule: consecutive distinct epoch times more than SEASON_GAP_MIN_DAYS apart begin a new
// season. On the current schedule the largest spacing inside a season is 5.0 d and the smallest gap
// between seasons is 108.2 d; the margin is checked at run time.
constexpr double SEASON_GAP_MIN_DAYS = 20.0;
// Sun-Earth L2, where Roman flies, as a fraction of an AU.
//
// L2 is on the Sun-Earth line, ~1.5e6 km beyond the Earth, so to leading order Roman's heliocentric
// position is Earth's scaled by (1 + L2_OFFSET_AU). This gives the joint fit a spatial baseline: the
// difference in impact parameter the two observatories see is
//     delta_u ~ L2_OFFSET_AU * piE  ~  1e-3  for a typical bulge event,
// concentrated in high-magnification, short-tE events. Roman's halo orbit about L2 (amplitude
// ~1e5-1e6 km) is not modelled; this is the mean offset only.
constexpr double L2_KM        = 1.5e6;
constexpr double L2_OFFSET_AU = L2_KM / AU_KM;   // ~0.01003
// Roman's photometric error. files/sigma_roman.txt is Penny et al. 2019 Fig. 4: single-epoch precision
// vs W149 (= F146) AB magnitude for a 46.8-s Cycle-7 exposure, with a 1 mmag floor. At load time its
// photon-noise part, sqrt(err^2 - floor^2), is shifted in magnitude so that the 5-sigma point
// (err = 1.0857/5) falls at ROMAN_DEPTH5_AB, the STScI depth for the GBTDS's 66-s exposure, and the
// floor is re-added; errRomanM then interpolates log(err) linearly in magnitude.
constexpr double ROMAN_PHOT_FLOOR = 0.001;   //mag, Penny et al. 2019 Table 2 "Error floor 1.0 mmag"
// ---------------------------------------------------------------------------------------
// Roman WFI per-exposure astrometric precision, F146, per coordinate (errRomanA), in milliarcseconds:
//     sigma = sqrt( (ROMAN_AST_K / SNR)^2 + ROMAN_AST_FLOOR^2 )
// SNR = 1.0857 / sigma_rand, sigma_rand the photon-noise part of the same exposure's photometric error
// (errRomanM, the 1 mmag floor removed in quadrature). The relation sigma = k / SNR is the statistical
// limit of a PSF fit [5] (eqs. 7-8 of [7]); k = 0.792 pixel was measured by [7] on the stpsf F146 PSF
// sampled on 0.11" pixels (FWHM 1.044 pixel; kappa = k / FWHM = 0.758, larger than a Gaussian's 0.425
// because the PSF is undersampled; their zero sub-pixel-offset value, the pessimistic one). With our
// SNR this reproduces [7]'s tabulated 66-s GBTDS curve (medium background) to 2-3% over F146 AB 19-26,
// so Roman's photometric and astrometric errors now share one SNR. Neither includes crowding or
// geometric distortion [2][7].
// The floor: "single-exposure precision for well-exposed point sources is 0.01 pixel, or about 1.1 mas"
// [1]; adopted by [2] and [7] (added in quadrature by [7]). A systematic, so it does not improve for
// brighter stars.
//
// Sources:
//   [1] Sanderson et al. 2019, arXiv:1712.05420 sec 1.1, improving ~10x when ~100 exposures are stacked.
//   [2] Lam et al. 2026, arXiv:2608.24998, sec. 4.2.2 and Fig. 5 (same tool as [7], at the 90th-percentile
//       background); pixels "0.11 arcsec"; each GBTDS exposure "66 seconds" at a "12.1 minute" cadence.
//   [7] McKinnon & van der Marel 2026, arXiv:2602.00310, Table 1 (F146: k = 0.792 px) and the tool's
//       data/roman_IM_66_6_gbtds_mid_5stripe_medium_F146_pos_errs.csv (Zenodo 10.5281/zenodo.18407082).
//
// PER EXPOSURE: the 0.1 mas figure in [1] is the daily-binned precision (~100 exposures stacked).
// l.erra[] is a per-epoch error and one row of RomanBaseline.dat is one 12.1-minute exposure, so the
// per-exposure numbers are the right ones here.
// ---------------------------------------------------------------------------------------
constexpr double ROMAN_PIX_MAS   = 110.0;                 //0.11 arcsec pixels [2]
constexpr double ROMAN_AST_FLOOR = 0.01 * ROMAN_PIX_MAS;  //1.1 mas: 1% centroiding [1][2]
constexpr double ROMAN_AST_K     = 0.792 * ROMAN_PIX_MAS; //87.12 mas: sigma = k / SNR, F146 [7]
// ---------------------------------------------------------------------------------------
// The astrometric noise model, three ways.
//
// PER COORDINATE. errRomanA / errlsstA give the 1D (x or y) per-exposure precision, as every source
// defines it (Lam et al. 2026 fn. 14; McKinnon & van der Marel 2026; Ivezic et al. "per observation
// per coordinate"). Each coordinate's variance is erra^2.
//
// THE REFERENCE POSITION IS FREE. Real astrometric fits solve for the source's position offset;
// otherwise the model's -u0 tetE sin(xi) term lets every exposure measure tetE from the source's
// absolute position and sigma(tetE) comes out optimistic. An offset per "frame group" is marginalised
// in closed form: F_g = sum_blocks F_k - (sum b_k)(sum b_k)^T / sum c_k, per coordinate.
//
// WHETHER ROMAN'S ERRORS AVERAGE DOWN is unknown. The 1.1 mas floor is the exposure-to-exposure
// scatter of dithered HST data, and the GBTDS dithers by several pixels and sub-pixel steps, so the
// floor itself is white; no source quantifies an ADDITIONAL error shared by many exposures. An
// error constant within a frame group costs nothing (its offset is free): what matters is correlation
// on timescales shorter than an event. Each event carries three forecasts; the per-exposure error
// (erra, floor included) is the same in all three, and N and P ADD a day-shared term to it:
//   W  white (the literature's assumption): one free offset per telescope; nothing added.
//   N  nominal: one free offset per Roman ROLL (static crowding and distortion biases change with the
//      PSF orientation) plus AST_SIGC_N per coordinate shared by all Roman exposures of the same day
//      (time-varying distortion residuals at "a few x 0.1% of a pixel", Bellini 2024 via Lam et al. 2026).
//   P  pessimistic: one free offset per Roman SEASON plus AST_SIGC_P = 1.1 mas shared within each day,
//      a day-shared error as large as the floor (a day of ~120 exposures then measures no better than
//      ~1.1 mas). No identified source is this large; it is an extreme bound.
// Rubin's errors are white with one offset in all three. The day blocks enter by Sherman-Morrison:
// for a block with weights w_i = 1/erra_i^2 and derivatives d_i, F_k = S_wdd - s^2 S_wd S_wd^T /
// (1 + s^2 S_w), b_k = S_wd / (1 + s^2 S_w), c_k = S_w / (1 + s^2 S_w), s = sigma_c.
// The main table columns (sigtetE_*, relMl_*, okB_*, condB_*) are W.
// ---------------------------------------------------------------------------------------
constexpr double AST_SIGC_N = 0.3;               //mas per coordinate per Roman day
constexpr double AST_SIGC_P = ROMAN_AST_FLOOR;   //1.1 mas
constexpr std::array<double, NAVAR> AST_SIGC = {0.0, AST_SIGC_N, AST_SIGC_P};
constexpr int    AST_MAX_SEASONS = 16;           //Roman seasons a P-variant offset can be keyed on
// THE BLEND OFFSET IS FREE. The unresolved neighbours sit at a light centroid offset from the source
// (s.blendOff, drawn in func_source), and as the source brightens the measured centroid slides from
// them toward it: a magnification-locked shift that real fits must model with two extra parameters per
// telescope (the offset's two coordinates), since the neighbours' positions are not known in advance.
// true: those parameters are fitted and marginalised in every variant; false: they are treated as known
// (the forecast then shows only what the blending costs in light, not in position).
constexpr bool   AST_BLEND_OFFSET_FREE = true;

// ==========================================================================================
// (6) LENS POPULATIONS AND MASS FUNCTIONS
// ==========================================================================================

// The bulge entry is first and is the default.
inline constexpr LensPopulation POPULATIONS[] = {
    {"bulge", "5", MassFunction::KROUPA_REMNANTS, 0.01, 30.0, false, 5,
     "Kroupa IMF + remnants; the present-day bulge population"},
    {"bh",    "bh", MassFunction::LOG_UNIFORM,     3.0, 1000.0, true, 0,
     "black holes, flat in log M over 3-1000 Msun"},
    {"ns",    "ns", MassFunction::NEUTRON_STAR,    1.0,    2.5, false, 0,
     "neutron stars, Gaussian about 1.35 Msun (Ozel & Freire 2016)"},
    {"besancon", "bes", MassFunction::BESANCON_CATALOGUE, 0.01, 30.0, false, 0,
     "ordinary stars; the lens is a random member of its component's Besancon list (mass and light)"},
    {"macho-uniform", "1", MassFunction::UNIFORM,      3.0, 5000.0, true, 1, "legacy MACHO search"},
    {"macho-m05",     "2", MassFunction::POWER_LAW_05, 3.0, 5000.0, true, 2, "legacy MACHO search"},
    {"macho-m1",      "3", MassFunction::POWER_LAW_10, 3.0, 5000.0, true, 3, "legacy MACHO search"},
    {"macho-m2",      "4", MassFunction::POWER_LAW_20, 3.0, 5000.0, true, 4, "legacy MACHO search"},
};
constexpr double u0m   = 3.0;
// Event-parameter draws (func_lens): u0 is uniform in [U0_MIN_DRAW, u0m]; the peak time t0 is uniform
// in [T0_MARGIN_DAYS, Tobs - T0_MARGIN_DAYS] days, i.e. at least 2 days inside the survey ends.
constexpr double U0_MIN_DRAW = 0.001;
constexpr double T0_MARGIN_DAYS = 2.0;
// ---- Kroupa (2001) initial mass function, dN/dM ~ M^-alpha, with the standard breaks ----
// Continuity coefficients are derived in drawKroupaInitialMass(); only breaks and slopes are named here.
constexpr double KROUPA_MI_MIN = 0.01;   //below the hydrogen-burning limit: brown dwarfs
constexpr double KROUPA_MI_MAX = 120.0;  //initial mass; nothing this heavy survives to today
constexpr double KROUPA_BREAK1 = 0.08;   //hydrogen-burning limit
constexpr double KROUPA_BREAK2 = 0.50;
constexpr double KROUPA_ALPHA1 = 0.3;    //0.01 - 0.08
constexpr double KROUPA_ALPHA2 = 1.3;    //0.08 - 0.50
constexpr double KROUPA_ALPHA3 = 2.3;    //0.50 - 120   (Salpeter-like)

// ---- Initial-to-final mass, for the remnants ----
// A bulge population is ~10 Gyr old, so everything born above the turnoff is already dead. The remnant
// type is set by the initial mass. Black-hole lenses are heavy, so tE ~ sqrt(Ml) is long and the event
// can span Roman's season gaps.
// ---- Neutron star masses, for the NEUTRON_STAR population ----
// Ozel & Freire (2016): double neutron stars cluster at 1.33 +/- 0.09 Msun, and the population is
// often summarised as a Gaussian near 1.35 Msun with a ~0.15 Msun spread. The truncation is physical:
// below ~1.1 Msun no supernova is known to leave a neutron star, and above ~2.2 Msun the equation of
// state gives a black hole. The real distribution is arguably bimodal (a recycled population above
// the canonical peak), so a single Gaussian understates the high-mass tail.
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
// Resolving the two lensing-induced images.
//
// A point lens makes two images of the source at theta_pm = 0.5 * (u +- sqrt(u^2+4)) * theta_E, so
//     Delta_theta(u) = theta_E * sqrt(u^2 + 4).
// For a stellar lens theta_E ~ 0.3 mas, far below any resolution, and the detector sees the summed flux
// and flux-weighted centroid. A black-hole lens is different because theta_E scales as sqrt(Ml).
//
// The separation is smallest at closest approach and grows as the source moves away, but the minor
// image's magnification A_minus = (u^2+2) / (2 u sqrt(u^2+4)) - 1/2 collapses on the same motion.
// A criterion at peak is too optimistic and one at maximum separation meaningless. Sajadian & Makler
// (arXiv:2608.16448, sec. 3) count data points that satisfy both conditions at once and call the images
// resolvable if at least three do; that count is accumulated per epoch inside the epoch loop.
//
// RESOLUTION THRESHOLD. That paper takes the Rubin criterion as Delta_theta >= sigma_r = D * sigma_a,
// with sigma_a the per-visit astrometric precision and D in [5, 100] depending on the source's S/N:
// ~5 at the faint limit (SNR 5), ~20 at SNR 100 for two similar stars, higher still when one image is
// much fainter (Ivezic, priv. comm. quoted therein). The count is recorded at both anchors
// (RESOLVE_D_FAINT, RESOLVE_D_BRIGHT). A third criterion, Delta_theta >= the PSF FWHM, is added for
// Roman, which is diffraction limited at 0.105 arcsec in F146; D*sigma_a on Roman's 1.1 mas floor would
// claim a resolution of a few mas, which no 2.4 m telescope delivers.
// ---------------------------------------------------------------------------------------
constexpr double RESOLVE_D_FAINT  = 5.0;   //D at the faint detection limit, SNR ~ 5
constexpr double RESOLVE_D_BRIGHT = 20.0;  //D at SNR ~ 100, images of similar brightness
// ---------------------------------------------------------------------------------------
// The detection threshold.
//
// A detection is declared when the lensing model beats a flat-baseline model by enough chi-squared:
//     dchi = chi2(flat baseline) - chi2(lensing model),
// a difference between nested models. Under the null it is chi-squared distributed with p degrees of
// freedom, p being the number of extra parameters of the lensing model (not the number of epochs), so
// the bar must not scale with the epoch count: a bar that does thresholds the mean per-epoch
// improvement rather than the total significance, and pooling many low-signal epochs would raise the
// joint bar without adding signal.
//
// A fixed bar also makes the tests monotone: chi1 and chi3 are accumulated over both instruments, so
// dchi = dchi_L + dchi_R exactly. With the same threshold on all three tests, either survey clearing
// the bar alone forces the joint sum over it too, so DET_ANOMALY cannot occur. This requires the
// SIGNED difference (see the note on fabs at the test site in src/sim/detect.cpp).
//
// The value 500 is Penny et al. 2019 (ApJS 241, 3), which adopts dchi2 > 500 against a flat baseline,
// keeping yields comparable with the Roman community's. It is deliberately conservative (a nominal
// 3-sigma bar on a few parameters would be nearer 20) because the real false-alarm population is
// systematics, variable stars and blending. Overridable with --dchi-det and recorded in
// run_provenance.txt, since every yield is conditioned on it.
// ---------------------------------------------------------------------------------------
constexpr double DCHI_DET_DEFAULT = 500.0; //delta-chi2 against a flat baseline [Penny+2019]
// Step-size multipliers on the astrometric Delta2[] (tetE, mus1, mus2, piE), chosen with
// ./fishertest --sweep-astro; 1.0 = 25%-of-value steps. The modelled centroid is linear in tetE, mus1
// and mus2 (blending and lens light included), so their finite differences are exact. piE acts through
// the parallax-bent trajectory and is not: the plateau is 1e-8..1e-2 and 1e-2 is used.
constexpr std::array<double, 4> kFDStepScaleB = {1.0, 1.0, 1.0, 1.0e-2};

// Photometric finite-difference steps of FisherM (Delta1[]). The base steps are ~25% of the parameter
// (u0 by 0.15, tE and t0 by 0.25*tE), far up the truncation-error branch; kFDStepScale shrinks them
// into the convergence plateau (flat to <0.2% over 1e-6..1e-3 of the base steps; 1e-4 sits two
// decades clear of the round-off wall). Re-run ./fishertest --sweep (tests/c3_step_sweep.py) after
// changing any of them.
constexpr double FD_STEP_U0       = 0.1507586576;      //u0 (absolute)
constexpr double FD_STEP_TE_FRAC  = 0.254674;          //tE and t0: step = tE * this [days]
constexpr double FD_STEP_PIE_FRAC = 0.2509463534656;   //piE: step = this * piE
constexpr double FD_STEP_XI_DEG   = 3.0;               //xi [deg]; converted to radians at the use
constexpr double FD_STEP_MBS      = 0.05;              //mbs0 / mbs1 [mag]; the model is linear in mbs, so exact
constexpr double kFDStepScale     = 1.0e-4;            //plateau scale applied to every Delta1[] and fb step
// fb0 / fb1 steps (blend fraction, bounded to [0,1]). The step pair is chosen from the epoch's own
// telescope's fb so that fb + step never leaves the range: fb < FB_BIN_LO -> {+SMALL, +LARGE};
// fb < FB_BIN_HI -> {-SMALL, +SMALL}; else {-SMALL, -LARGE}. (Then * kFDStepScale.)
constexpr double FB_BIN_LO     = 0.15;
constexpr double FB_BIN_HI     = 0.85;
constexpr double FB_STEP_SMALL = 0.07;
constexpr double FB_STEP_LARGE = 0.15;
// A telescope contributing fewer than kMinTeleEpochs epochs to an event is left out of the photometric
// matrices altogether: its flux pair (fb, mbs) cannot be constrained and with one epoch the joint
// matrix goes singular. The counts passed to activePhotParams are therefore >= kMinTeleEpochs or zero.
constexpr int kMinTeleEpochs = 3;

// ==========================================================================================
// (8) OUTPUT HISTOGRAM GRIDS
// ==========================================================================================

constexpr int GG = 100;
constexpr double tE_min  = 0.0;//days
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
inline constexpr const char* PATH_SIGMA_ROMAN    = "./files/sigma_roman.txt";
inline constexpr const char* PATH_EXT_TABLES     = "./files/ext/ext_tables.dat";
inline constexpr const char* PATH_LENS_ML        = "./CMD/components/lens_ml.dat";
inline constexpr const char* PATH_CMD_THIN       = "./CMD/components/thin_disk.dat";
inline constexpr const char* PATH_CMD_BULGE      = "./CMD/components/bulge.dat";
inline constexpr const char* PATH_CMD_THICK      = "./CMD/components/thick_disk.dat";
inline constexpr const char* PATH_CMD_HALO       = "./CMD/components/halo.dat";
// Output directory (LpLMC / EfLMC / MapLMC / run_provenance) and the Disk_model debug dumps; both end in a slash.
inline constexpr const char* PATH_OUT_DIR        = "./files/MONTLMC/files/";
inline constexpr const char* PATH_DENSITY_DIR    = "./files/density/";

#endif // PARAMETERS_H
