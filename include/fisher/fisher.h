// The Fisher-matrix forecast: per-survey partitions, active parameter sets, the covarian storage,
// and the detection/synergy classifications of one event.
#ifndef ROMAN_FISHER_FISHER_H
#define ROMAN_FISHER_FISHER_H

#include "common.h"
#include "types.h"
#include "fisher/linalg.h"

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

void   ErrorCal(covarian & co, lens &l, source &s);
void   FisherM(source & s, lens & l, astromet & as,  covarian & co, int);

#endif // ROMAN_FISHER_FISHER_H
