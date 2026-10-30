// The Fisher-matrix forecast: per-survey partitions, active parameter sets, the covarian storage,
// and the detection/synergy classifications of one event.
#ifndef ROMAN_FISHER_FISHER_H
#define ROMAN_FISHER_FISHER_H

#include "common.h"
#include "types.h"
#include "fisher/linalg.h"
#include "events/vbm_model.h"

// ---------------------------------------------------------------------------------------------
// Which subset of a single event's light curve a Fisher matrix was accumulated from.
//
// The model and its derivatives are identical across the three; only the set of epochs summed
// over differs. Epochs contribute independent, positive semi-definite terms and partition by
// observatory, so F[SJOINT] == F[SRUBIN] + F[SROMAN] element by element (asserted in
// tests/fisher_fixture.cpp), and sigma_joint can never exceed either single-survey sigma.
enum SurveyIdx { SJOINT = 0, SRUBIN = 1, SROMAN = 2, NSURV = 3 };

// Maps a per-epoch telescope tag (lens::tele[i]: 0 = Rubin, 1 = Roman) to its survey index.
inline int surveyOfTele(int tele) { return (tele == 0) ? SRUBIN : SROMAN; }

// Which photometric parameters a survey partition can constrain.
//
// Rubin's epochs carry no information about Roman's flux parameters (fb1, mbs1) and vice versa,
// so those rows and columns of a single-survey matrix are identically zero and the full Nx x Nx
// matrix would be singular by construction. Each partition therefore inverts only its active
// submatrix; parameters outside it report sigma = -1.
//
// The joint set depends on which surveys contributed epochs to THIS event: a short event peaking
// in a Roman gap has no Roman data, so the joint fit falls back to Rubin's parameter set rather
// than going singular. kMinTeleEpochs (minimum epochs for a telescope to enter the photometric
// matrices) is in config/parameters.h.
//
// ln rho (9) is in the set when the partition has data and the source is a disc (rhoFree; RHO_PRIOR_SIGLN in
// config/parameters.h), which FisherM records per partition in covarian::rhoFit. Whether the data MEASURE it is
// covarian::rhoMeas (SIGLNRHO_MAX).

inline std::vector<int> activePhotParams(int surv, int nRubinEpochs, int nRomanEpochs, bool rhoFree)
{
    const std::vector<int> shared = {0, 1, 3, 4, 5};  //u0, tE, piE, xi, t0
    const std::vector<int> rubinFlux = {2, 6};        //fb0, mbs0
    const std::vector<int> romanFlux = {7, 8};        //fb1, mbs1

    std::vector<int> out = shared;
    const bool wantRubin = (surv == SRUBIN) || (surv == SJOINT && nRubinEpochs > 0);
    const bool wantRoman = (surv == SROMAN) || (surv == SJOINT && nRomanEpochs > 0);
    if (wantRubin) out.insert(out.end(), rubinFlux.begin(), rubinFlux.end());
    if (wantRoman) out.insert(out.end(), romanFlux.begin(), romanFlux.end());
    if (rhoFree)   out.push_back(9);                  //ln rho
    std::sort(out.begin(), out.end());
    return out;
}

struct covarian {
    int    sign;
    int    flagi;
    
    double deter;
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

    // Astrometric forecast under the three noise variants (AV_W, AV_N, AV_P; see AST_SIGC).
    // Index [variant][survey][param]. Variant W equals Erb/okB/condB/relMl; the Rubin partition
    // is the same in all three (only Roman's noise and offsets differ). -1 sentinels as elsewhere.
    std::array<std::array<std::array<double, Ny>, NSURV>, NAVAR> ErbV{};
    std::array<std::array<int, NSURV>, NAVAR>    okBV{};
    std::array<std::array<double, NSURV>, NAVAR> condBV{};
    std::array<std::array<double, NSURV>, NAVAR> relMlV{};
    // cov(thetaE, piE) of each variant's astrometric inverse [mas, per unit piE], the covariance the mass and
    // distance errors need when the astrometric piE is the route (ErrorCal, lnCovTetPi). 0 where okBV = 0.
    std::array<std::array<double, NSURV>, NAVAR> covTetPiV{};
    std::array<std::array<std::array<double, Ny * Ny>, NSURV>, NAVAR> FBV{};  // the matrices

    // Per-survey epoch counts and characterizability flags. A partition with no epochs, or fewer
    // epochs than free parameters, is rank-deficient and is reported as not-characterizable
    // rather than inverted into a meaningless ~1e10 sigma.
    std::array<int, NSURV> nepochA; //epochs contributing to each photometric matrix
    std::array<int, NSURV> okA;     //1 = inverted and usable, 0 = not characterizable
    std::array<int, NSURV> okB;     //same for the astrometric matrix

    // Condition number (lambda_max/lambda_min) of the NORMALIZED information matrix, per survey.
    // Rows and columns are divided by sqrt(F_ii), so the value reflects genuine parameter
    // degeneracy (classically u0-tE-fb) rather than the parameters' differing units. Stored so
    // the cut can be chosen in analysis. -1.0 = not computed (partition rejected earlier).
    std::array<double, NSURV> condA;
    std::array<double, NSURV> condB;

    // Fractional 1-sigma on the lens mass, per survey. Ml = tetE / (kappa * piE) is not fitted, so
    // the fractional errors of tetE (astrometric matrix: sub-mas centroid motion, Roman; and, when the
    // finite source is seen, photometric: theta* / rho) and piE (photometric matrix: parallax over a long
    // baseline, Rubin) add in quadrature, less twice their covariance c = cov(ln tetE, ln piE): -wp cov(ln rho,
    // ln piE) with the photometric piE, (1 - wp) cov(tetE, piE) / (tetE piE) with the astrometric one (wp the
    // photometric weight in the combined tetE; ErrorCal, lnCovTetPi). relMlV[v] is the same for noise variant v.
    // -1.0 = not measurable in this partition, never a measured value.
    std::array<double, NSURV> relMl;

    // The lens's motion on the sky, per partition (ErrorCal; the model is lensMotion in events/vbm_model.h, the
    // propagation lensMotionCovariance / lensMotionPropagate below):
    //   relMuL   sigma(|mu_L|) / |mu_L|, mu_L the lens's heliocentric proper motion;
    //   sigVL    sqrt(lambda_max(covVL)) [km/s]: the 1-sigma error of the lens's velocity on the sky, v_L = D_L mu_L,
    //            along its worst-measured direction (the direction xi is the weak parameter, so the error ellipse is
    //            elongated). Independent of the reference frame: the Sun's velocity added to go to the Galactic rest
    //            frame shifts v_L, not its covariance. A per-axis bound, not the error of |v_L|;
    //   covMuL, covVL   the 2x2 covariances (row major) of mu_L [mas/day] and v_L [km/s] on (e_l, e_b).
    // -1 = not measured: the partition has no photometric or no astrometric matrix, or no piE or thetaE by any
    // route (relMl < 0); the covariances are then 0. lmCheck is the largest relative departure of the closed-form
    // lensMotion at the true parameters from the draw's own (mul1, mul2) and v_L, set by ErrorCal; -1 before it.
    std::array<double, NSURV> relMuL;
    std::array<double, NSURV> sigVL;
    std::array<std::array<double, 4>, NSURV> covMuL{}, covVL{};
    double lmCheck = -1.0;

    // The finite source (--source-model espl; see RHO_PRIOR_SIGLN in config/parameters.h). Per partition:
    //   dchiFS      the chi-squared by which that partition's epochs separate the uniform disc from a point
    //               source at the true parameters, a diagnostic of how strongly the disc is seen; -1 = not computed
    //               (point source, or no epochs).
    //   rhoFit      1 = ln rho is a parameter of the partition's photometric matrix (a disc, and the partition has
    //               data, unless rhoFixed); Era[q][9] is its posterior sigma, the prior included. 0 = rho held at
    //               its true value (point source, no data, rhoFixed, or the 10-parameter inversion failed:
    //               rhoFallback), Era[q][9] = -1.
    //   sigLnRhoData  the sigma(ln rho) the partition's DATA alone give, 1 / sqrt(1 / Era[q][9]^2 - 1 / prior^2) (the
    //               prior touches only F[9][9], so this is exact); -1 = not fitted, or no information in the data.
    //   rhoMeas     1 = the data measure rho: rhoFit and sigLnRhoData <= sigLnRhoMax. Only then the photometric
    //               thetaE = theta* / rho is used.
    //   relTetEph   fractional 1-sigma on thetaE = theta* / rho, sqrt(Era[q][9]^2 + THETASTAR_SIGLN^2); -1 where rho
    //               is not measured.
    //   relTetEc    the fractional 1-sigma on thetaE from the astrometric and photometric routes together,
    //               inverse-variance combined; the one route that exists if only one does; -1 = neither.
    //   corrRhoPiE  correlation of ln rho with the photometric piE in the partition's inverse matrix: how
    //               far the combination of the two thetaE and piE errors as uncorrelated (ErrorCal) is off.
    //               NaN where rho is not measured. A diagnostic, not an output.
    // rhoFixed, rhoPriorSigLn and sigLnRhoMax are the production settings unless a test overrides them: rhoFixed
    // holds rho at its true value everywhere (ln rho never a parameter); rhoPriorSigLn infinite drops the prior;
    // sigLnRhoMax infinite makes rho measured wherever it is fitted and the data give it any information.
    std::array<double, NSURV> dchiFS;
    std::array<int, NSURV>    rhoFit;
    std::array<int, NSURV>    rhoMeas;
    std::array<int, NSURV>    rhoFallback;   //1 = ln rho was a parameter but its matrix did not invert (see FisherM)
    std::array<double, NSURV> sigLnRhoData;
    std::array<double, NSURV> relTetEph;
    std::array<double, NSURV> relTetEc;
    std::array<double, NSURV> corrRhoPiE;
    bool   rhoFixed = false;
    double rhoPriorSigLn = RHO_PRIOR_SIGLN;
    double sigLnRhoMax = SIGLNRHO_MAX;

    // Multiplicative override on each photometric parameter's finite-difference step; 1.0 is the
    // production value. A runtime knob so the step-size sweep can run many values in one process.
    std::array<double, Nx> deltaScale;
    // Same for the astrometric steps Delta2[].
    std::array<double, Ny> deltaScaleB{1.0, 1.0, 1.0, 1.0};
    // Day-shared astrometric error of each noise variant [mas per coordinate]: AST_SIGC unless
    // --ast-sigc overrides N and P (main.cpp).
    std::array<double, NAVAR> astSigc = AST_SIGC;

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
            relMuL[q] = -1.0;
            sigVL[q] = -1.0;
            dchiFS[q] = -1.0;
            rhoFit[q] = 0;
            rhoMeas[q] = 0;
            rhoFallback[q] = 0;
            sigLnRhoData[q] = -1.0;
            relTetEph[q] = -1.0;
            relTetEc[q] = -1.0;
            corrRhoPiE[q] = std::numeric_limits<double>::quiet_NaN();
        }
        for (int k = 0; k < Nx; ++k) {
            deltaScale[k] = 1.0;
        }
        if (!summA || !summB) {
            throw std::runtime_error("GSL matrix allocation failed");
        }
    }

    // Disable copy, allow move
    covarian(const covarian&) = delete;
    covarian& operator=(const covarian&) = delete;
    covarian(covarian&&) = default;
    covarian& operator=(covarian&&) = default;
};
// Detection taxonomy. An event counts as detected only if the joint fit detects it: the joint
// stream contains strictly more data than either survey alone. The classes below say which
// surveys ALSO detect it on their own, which is what the joint-fit science case is about.
// Events that one survey cannot characterize alone are labelled, not dropped: they are the
// strongest evidence for the joint fit, and dropping rows with a missing single-survey sigma
// would bias the reported gain downward.
enum DetClass {
    DET_NONE         = 0, //nothing detected it
    DET_JOINT_ONLY   = 1, //only the combined stream -- neither telescope alone would have found it
    DET_RUBIN_JOINT  = 2, //Rubin alone, and the joint fit
    DET_ROMAN_JOINT  = 3, //Roman alone, and the joint fit
    DET_BOTH_JOINT   = 4, //both telescopes alone, and the joint fit
    DET_ANOMALY      = 5, //a telescope detected it but the joint test did NOT -- see below
    NDETCLASS        = 6,
};

// DET_ANOMALY should stay empty: adding data cannot destroy signal. It is not currently zero
// because the detection test thresholds the MEAN per-epoch chi-squared improvement (dchi vs
// 2*ndw), so pooling a survey with many low-signal epochs raises the joint bar and can veto a
// detection the other survey made alone. It is counted separately so the inconsistency stays
// visible.
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

// Classify from the photometric characterizability flags. SYN_JOINT_ONLY and the two *_ONLY
// classes are the scientifically interesting ones.
inline int synergyClass(const covarian& co)
{
    if (!co.okA[SJOINT])                    return SYN_NONE;
    if ( co.okA[SRUBIN] &&  co.okA[SROMAN]) return SYN_BOTH_ALONE;
    if ( co.okA[SRUBIN] && !co.okA[SROMAN]) return SYN_RUBIN_ONLY;
    if (!co.okA[SRUBIN] &&  co.okA[SROMAN]) return SYN_ROMAN_ONLY;
    return SYN_JOINT_ONLY;
}

// The active photometric parameters of partition `surv` of an event whose matrices are in `co`.
inline std::vector<int> activePhotParams(const covarian& co, int surv)
{
    return activePhotParams(surv, co.nepochA[SRUBIN], co.nepochA[SROMAN], co.rhoFit[surv] != 0);
}

// The parameters of the lens-motion model lensMotion, in the order its covariance is given.
enum LensMotionPar { LM_TE = 0, LM_XI, LM_PIE, LM_TETE, LM_MUS1, LM_MUS2, NLM };

// The covariance of those six parameters in partition q, from the photometric matrix inverA[q] (tE, xi, piE,
// ln rho) and the astrometric one inverB[q] (thetaE, mus1, mus2, piE), which are independent: the astrometric
// matrix holds the photometric parameters fixed and vice versa, so the blocks do not talk. piE and thetaE are
// taken by the routes relMl takes. False (S = 0) where relMl[q] < 0 or a matrix is missing.
bool   lensMotionCovariance(const covarian& co, const lens& l, int q, double S[NLM][NLM]);
// Linear propagation of S through lensMotion at the point p (central differences): the 2x2 covariances of mu_L
// [mas/day] and v_L [km/s].
void   lensMotionPropagate(const MotionContext& c, const double p[NLM], const double S[NLM][NLM],
                           double Smu[2][2], double Sv[2][2]);

// How often, over the run so far, ln rho was a parameter of a partition whose matrix failed to invert with it (the
// partition is then forecast with rho held at its true value) and how often the nine-parameter matrix failed too.
long   rhoFallbackCount();
long   rhoFallbackLostCount();

void   ErrorCal(covarian & co, lens &l, source &s);
void   FisherM(source & s, lens & l, astromet & as,  covarian & co, int);

#endif // ROMAN_FISHER_FISHER_H
