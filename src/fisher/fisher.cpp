// FisherM (information matrices by finite differences) and ErrorCal (sigmas from the inverses).
#include "fisher/fisher.h"
#include "events/lightcurve.h"
#include "surveys/noise.h"


void FisherM(source & s, lens & l, astromet & as,  covarian & co, int ndw)
{
    co.flagi =+ 1;
    int tt; // tev;

// (K + H + J)/3 = F146
// Photometry
    co.Delta1[0] = FD_STEP_U0;//u0
    co.Delta1[1] = l.tE * FD_STEP_TE_FRAC;//tE [days]
    co.Delta1[3] = FD_STEP_PIE_FRAC * l.piE;//piE
    co.Delta1[4] = FD_STEP_XI_DEG * M_PI / 180.0; //for xi [radian]
    co.Delta1[5] = l.tE * FD_STEP_TE_FRAC; // t0 [days] -- same fractional-of-tE step as tE
    co.Delta1[6] = FD_STEP_MBS; // mbs0 [mag] -- Rubin baseline magnitude
    co.Delta1[8] = FD_STEP_MBS; // mbs1 [mag] -- Roman baseline magnitude

    // The legacy step constants are ~25% of the parameter value, far up the truncation-error
    // branch of the step-size sweep (./fishertest --sweep, tests/c3_step_sweep.py). kFDStepScale
    // (see config/parameters.h) scales them into the plateau, where sigma is flat to <0.2%; below
    // ~1e-6 of the legacy value round-off dominates. Oversized steps manufacture information in
    // near-degenerate directions, so the correct sigmas for short events are much larger.
    for (int q = 0; q < Nx; ++q) co.Delta1[q] *= kFDStepScale * co.deltaScale[q];
    // mbs enters with unit slope, so its finite difference is exact for any step. fb0 (index 2)
    // and fb1 (index 7) use the telescope-keyed co.bb[] steps set in the data loop below.

    // Three matrices, zeroed together: joint, Rubin-only, Roman-only (see SurveyIdx).
    for (int q = 0; q < NSURV; ++q) {
        co.nepochA[q] = 0;
        co.okA[q] = 0;
        co.okB[q] = 0;
        // Reset here: a partition rejected for too few epochs never reaches invert_matrix.
        co.condA[q] = -1.0;
        co.condB[q] = -1.0;
        for (int j = 0; j < Nx; ++j) {
            for (int k = 0; k < Nx; ++k) {
                gsl_matrix_set(co.inputA[q].get(), j, k, 0.0);
                gsl_matrix_set(co.inverA[q].get(), j, k, 0.0);
            }
        }
    }



    // A telescope with fewer than kMinTeleEpochs epochs is left out (see config/parameters.h).
    std::array<int, 2> nTele{0, 0};
    for (int i = 0; i < ndw; ++i) nTele[int(l.tele[i]) == 1 ? 1 : 0] += 1;

    for (int i = 0; i < ndw; ++i) {//data
        tt = int(l.tele[i]);//telescope[0,1] LSST, ELT
        if (nTele[tt == 1 ? 1 : 0] < kMinTeleEpochs) continue;
        const int surv = surveyOfTele(tt); // which single-survey matrix this epoch also feeds
        co.nepochA[SJOINT] += 1;
        co.nepochA[surv]   += 1;

        // The fb step is binned on THIS epoch's telescope blend fraction s.fb[tt] (not always
        // s.fb[0]): bins and steps are chosen so fb[tt] + bb never leaves [0,1], which holds only
        // if bb is computed from the fb it perturbs. Roman's F146 blend fraction is routinely
        // higher than Rubin's r-band one.
        if (s.fb[tt] < FB_BIN_LO)      {co.bb[0] =+ FB_STEP_SMALL; co.bb[1] =+ FB_STEP_LARGE;}
        else if (s.fb[tt] < FB_BIN_HI) {co.bb[0] =- FB_STEP_SMALL; co.bb[1] =+ FB_STEP_SMALL;}
        else                      {co.bb[0] =- FB_STEP_SMALL; co.bb[1] =- FB_STEP_LARGE;}

        // Apply the plateau scaling and clamp so the step cannot leave the physical range [0,1]
        // (the bin edges guarantee this only at unscaled sizes). Note the outer fb bins are two
        // same-signed (one-sided) steps, so fb carries an O(h) bias there; kFDStepScale makes it
        // negligible.
        {
            const double fscale = kFDStepScale * co.deltaScale[(tt == 0) ? 2 : 7];
            for (int b = 0; b < 2; ++b) {
                double step = co.bb[b] * fscale;
                const double lo = 1.0e-6 - s.fb[tt];        //keeps fb strictly above 0
                const double hi = 1.0 - 1.0e-6 - s.fb[tt];  //keeps fb strictly below 1
                if (step < lo) step = lo;
                if (step > hi) step = hi;
                co.bb[b] = step;
            }
        }

        // Unperturbed model magnitude for THIS epoch, recomputed under the current observer
        // configuration rather than taken from the cached l.magn[i]. Each derivative is
        // (model(theta + Delta) - reference) / Delta, so the reference must share the observing
        // geometry (as.satScale). When the two differ, their difference leaks into every
        // derivative as a constant dm/Delta; it cancels in symmetric stencils but not in the
        // cross-telescope rows, fb's outer bins, or any one-sided stencil.
        lightcurve(s, l, as, l.timn[i], tt);
        s.Astar = (s.ut * s.ut + 2.0) / std::sqrt(s.ut * s.ut * (s.ut * s.ut + 4.0));
        const double magn0 = s.mbs[tt] - 2.5 * std::log10(s.Astar * s.fb[tt] + 1.0 - s.fb[tt]);

        for (int j = 0; j < Nx; ++j) {
            for (int h = 0; h < 2; ++h) {

                if (j == 0) {co.diff = double(+co.Delta1[j] * sig[h]) ;      l.u0 += co.diff;}
                if (j == 1) {co.diff = double(+co.Delta1[j] * sig[h]) ;      l.tE += co.diff;}
                // fb0 is Rubin's blend fraction: perturb s.fb[0] on Rubin epochs only (perturbing
                // s.fb[tt] would make parameters 2 and 7 identical on Roman epochs).
                if (j == 2) {co.diff = (tt == 0) ? double(co.bb[h]) : 1.0;
                             if (tt == 0) s.fb[0] += co.diff;}
                // Per-telescope flux parameters affect only their own telescope's epochs; the
                // other telescope's derivative is zero. co.diff is a nonzero dummy there.
                if (j == 6) {co.diff = (tt == 0) ? double(+co.Delta1[j] * sig[h]) : 1.0;
                                 if (tt == 0) s.mbs[0] += co.diff;}
                if (j == 7) {co.diff = (tt == 1) ? double(co.bb[h]) : 1.0;
                                 if (tt == 1) s.fb[1]  += co.diff;}
                if (j == 8) {co.diff = (tt == 1) ? double(+co.Delta1[j] * sig[h]) : 1.0;
                                 if (tt == 1) s.mbs[1] += co.diff;}
                if (j == 3) {co.diff = double(+co.Delta1[j] * sig[h]) ;     l.piE += co.diff;}
                if (j == 4) {co.diff = double(+co.Delta1[j] * sig[h]) ;      s.xi += co.diff;}
                if (j == 5) {co.diff = double(+co.Delta1[j] * sig[h]) ;      l.t0 += co.diff;}

                // Perturbed model evaluated with the observer that produced THIS datum.
                lightcurve(s, l, as, l.timn[i], int(l.tele[i]));
                s.Astar = (s.ut * s.ut + 2.0) / std::sqrt(s.ut * s.ut * (s.ut * s.ut + 4.0));
                co.magw = s.mbs[tt] - 2.5 * std::log10(s.Astar * s.fb[tt] + 1.0 - s.fb[tt]);
                co.derm1[h] = double(co.magw - magn0) / co.diff;

                CHECK(l.tE > 0.0);
                CHECK(s.fb[tt] > 0.0);
                CHECK(s.fb[tt] <= 1.0);
                CHECK(l.piE > 0.0);
                CHECK(l.u0 != 0.0);
                CHECK(co.diff != 0.0);
                CHECK(l.errm[i] > 0.0);

                if (j == 0)      l.u0 -= co.diff;
                if (j == 1)      l.tE -= co.diff;
                if (j == 2 and tt == 0)  s.fb[0] -= co.diff;
                if (j == 3)     l.piE -= co.diff;
                if (j == 4)      s.xi -= co.diff;
                if (j == 5)      l.t0 -= co.diff;
                if (j == 6 and tt == 0) s.mbs[0] -= co.diff;
                if (j == 7 and tt == 1)  s.fb[1] -= co.diff;
                if (j == 8 and tt == 1) s.mbs[1] -= co.diff;
            }

            co.derm1f = double(co.derm1[0] + co.derm1[1]) * 0.5;

            for (int k = 0; k <= j; ++k) {
                for (int h = 0; h < 2;  ++h) {

                    if (k == 0) {co.diff = double(+co.Delta1[k] * sig[h]) ;     l.u0 += co.diff;}
                    if (k == 1) {co.diff = double(+co.Delta1[k] * sig[h]) ;     l.tE += co.diff;}
                    if (k == 2) {co.diff = (tt == 0) ? double(co.bb[h]) : 1.0;
                                 if (tt == 0) s.fb[0] += co.diff;}
                    // Per-telescope flux parameters affect only their own telescope's epochs; the
                    // other telescope's derivative is zero. co.diff is a nonzero dummy there.
                    if (k == 6) {co.diff = (tt == 0) ? double(+co.Delta1[k] * sig[h]) : 1.0;
                                     if (tt == 0) s.mbs[0] += co.diff;}
                    if (k == 7) {co.diff = (tt == 1) ? double(co.bb[h]) : 1.0;
                                     if (tt == 1) s.fb[1]  += co.diff;}
                    if (k == 8) {co.diff = (tt == 1) ? double(+co.Delta1[k] * sig[h]) : 1.0;
                                     if (tt == 1) s.mbs[1] += co.diff;}
                    if (k == 3) {co.diff = double(+co.Delta1[k] * sig[h]) ;    l.piE += co.diff;}
                    if (k == 4) {co.diff = double(+co.Delta1[k] * sig[h]) ;     s.xi += co.diff;}
                    if (k == 5) {co.diff = double(+co.Delta1[k] * sig[h]) ;     l.t0 += co.diff;}

                    lightcurve(s, l, as, l.timn[i], int(l.tele[i]));
                    s.Astar = (s.ut * s.ut + 2.0) / std::sqrt(s.ut * s.ut * (s.ut * s.ut + 4.0));
                    co.magw = s.mbs[tt] - 2.5 * std::log10(s.Astar * s.fb[tt] + 1.0 - s.fb[tt]);
                    co.derm2[h] = double(co.magw - magn0) / co.diff;

                    CHECK(l.tE > 0.0);
                    CHECK(s.fb[tt] >= 0.0);
                    CHECK(s.fb[tt] <= 1.0);
                    CHECK(l.piE > 0.0);
                    CHECK(l.u0 != 0.0);
                    CHECK(co.diff != 0.0);
                    CHECK(l.errm[i] > 0.0);

                    if (k == 0)      l.u0 -= co.diff;
                    if (k == 1)      l.tE -= co.diff;
                    if (k == 2 and tt == 0)  s.fb[0] -= co.diff;
                    if (k == 3)     l.piE -= co.diff;
                    if (k == 4)      s.xi -= co.diff;
                    if (k == 5)      l.t0 -= co.diff;
                    if (k == 6 and tt == 0) s.mbs[0] -= co.diff;
                    if (k == 7 and tt == 1)  s.fb[1] -= co.diff;
                    if (k == 8 and tt == 1) s.mbs[1] -= co.diff;
                }

                co.derm2f = double(co.derm2[0] + co.derm2[1]) * 0.5;

                // Accumulate the information sum F_jk = sum_i (dm_i/dtheta_j)(dm_i/dtheta_k)/sigma_i^2.
                // Each epoch feeds the joint matrix and exactly one single-survey matrix, so
                // F[SJOINT] == F[SRUBIN] + F[SROMAN] holds element by element.
                {
                    const double contrib = co.derm1f * co.derm2f / (l.errm[i] * l.errm[i]);
                    gsl_matrix_set(co.inputA[SJOINT].get(), j, k,
                                   gsl_matrix_get(co.inputA[SJOINT].get(), j, k) + contrib);
                    gsl_matrix_set(co.inputA[surv].get(), j, k,
                                   gsl_matrix_get(co.inputA[surv].get(), j, k) + contrib);
                }
            }
        }//end of for J
    }//end of data for

    for (int q = 0; q < NSURV; ++q) {
        for (int j = 0; j < Nx; ++j) {
            for (int k = (j + 1); k < Nx; ++k) {
                double element = gsl_matrix_get(co.inputA[q].get(), k, j);
                gsl_matrix_set(co.inputA[q].get(), j, k, element);
            }
        }
        // A partition with fewer epochs than its ACTIVE parameters is rank-deficient (e.g. a short
        // event peaking in a Roman gap has no Roman data): report it as not characterizable
        // rather than inverting it into a meaningless ~1e10 sigma.
        co.okA[q] = (co.nepochA[q] >= static_cast<int>(
                         activePhotParams(q, co.nepochA[SRUBIN], co.nepochA[SROMAN]).size()))
                        ? invert_matrix(co, 0, q) : 0;
    }

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///                   Astrometry                                            ///
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    co.Delta2[0] = l.tetE * 0.2500467345692345   * kFDStepScaleB[0] * co.deltaScaleB[0];
    co.Delta2[1] = s.mus1 * 0.2500465923465493   * kFDStepScaleB[1] * co.deltaScaleB[1];
    co.Delta2[2] = s.mus2 * 0.25000465936443     * kFDStepScaleB[2] * co.deltaScaleB[2];
    co.Delta2[3] = 0.25000375423543264 * l.piE   * kFDStepScaleB[3] * co.deltaScaleB[3];
    for (int q = 0; q < NSURV; ++q)
    for(int j = 0; j < Ny; ++j){
        for(int k = 0; k < Ny; ++k){
            gsl_matrix_set(co.inputB[q].get(), j, k, 0.0);
            gsl_matrix_set(co.inverB[q].get(), j, k, 0.0);
        }
    }

    // ---- Per-coordinate weights, free reference positions, and the three noise variants
    // (AST_SIGC in config/parameters.h), from ONE pass over the epochs. ----
    // Per epoch the derivatives of each sky coordinate c (0 = x, 1 = y) are computed once.
    // Sums kept per coordinate: S_w = sum w, S_wd[c][j] = sum w d_cj, S_wdd[c][jk]. Rubin's epochs
    // go into one white group; Roman's go into day blocks (the sigma_c correlation unit), tagged
    // with season and roll so the variants can group them differently.
    struct AstSums {
        double Sw = 0.0;
        std::array<std::array<double, Ny>, 2>      Swd{};
        std::array<std::array<double, Ny * Ny>, 2> Swdd{};
        int season = -1, roll = -1;
        void add(double w, const std::array<double, Ny>& dx, const std::array<double, Ny>& dy) {
            Sw += w;
            for (int c = 0; c < 2; ++c) {
                const auto& d = c ? dy : dx;
                for (int j = 0; j < Ny; ++j) {
                    Swd[c][j] += w * d[j];
                    for (int k = 0; k < Ny; ++k) Swdd[c][j * Ny + k] += w * d[j] * d[k];
                }
            }
        }
    };
    AstSums rubin;
    std::map<long, AstSums> romanDays;
    std::array<double, Ny> dx{}, dy{};

    for (int i = 0; i < ndw; ++i) {
        // Reference and perturbations use the observer that produced THIS datum.
        lightcurve(s, l, as, l.timn[i], int(l.tele[i]));
        const double soux0 = s.pos1c, souy0 = s.pos2c;

        for (int j = 0; j < Ny; ++j) {
            for (int h = 0; h < 2; ++h) {
                // Central stencil for all four astrometric parameters.
                co.diff = double(co.Delta2[j] * sig[h]);
                if (j == 0) l.tetE += co.diff;
                if (j == 1) s.mus1 += co.diff;
                if (j == 2) s.mus2 += co.diff;
                if (j == 3) l.piE  += co.diff;

                lightcurve(s, l, as, l.timn[i], int(l.tele[i]));
                co.dera1[h] = double(s.pos1c - soux0) / co.diff;
                co.derb1[h] = double(s.pos2c - souy0) / co.diff;

                CHECK(l.tetE > 0.0);
                CHECK(l.piE > 0.0);
                CHECK(co.diff != 0.0);
                // A null derivative (epoch at t0 or far from the peak) is legitimate and contributes nothing.

                if (j==0) l.tetE -= co.diff;
                if (j==1) s.mus1 -= co.diff;
                if (j==2) s.mus2 -= co.diff;
                if (j==3)  l.piE -= co.diff;
            }
            dx[j] = (co.dera1[0] + co.dera1[1]) * 0.5;
            dy[j] = (co.derb1[0] + co.derb1[1]) * 0.5;
        }
        CHECK(l.erra[i] > 0.0);
        const double w = 1.0 / (l.erra[i] * l.erra[i]);       // per coordinate
        if (int(l.tele[i]) == 0) {
            rubin.add(w, dx, dy);
        } else {
            CHECK(l.rseas[i] >= 0 and l.rseas[i] < AST_MAX_SEASONS);
            CHECK(l.rroll[i] == 0 or l.rroll[i] == 1);
            AstSums& d = romanDays[long(std::floor(l.timn[i]))];
            if (d.season < 0) { d.season = l.rseas[i]; d.roll = l.rroll[i]; }
            CHECK(d.season == l.rseas[i] and d.roll == l.rroll[i]);
            d.add(w, dx, dy);
        }
    }

    // A frame group: sum of block matrices, the block offset-vectors b per coordinate, and c.
    struct Group {
        std::array<double, Ny * Ny> F{};
        std::array<std::array<double, Ny>, 2> b{};
        double c = 0.0;
        void fold(const AstSums& k, double sc) {          // one block, correlated at sigma_c = sc
            const double s2 = sc * sc, den = 1.0 + s2 * k.Sw;
            for (int cc = 0; cc < 2; ++cc)
                for (int j = 0; j < Ny; ++j) {
                    b[cc][j] += k.Swd[cc][j] / den;
                    for (int m = 0; m < Ny; ++m)
                        F[j * Ny + m] += k.Swdd[cc][j * Ny + m] - s2 * k.Swd[cc][j] * k.Swd[cc][m] / den;
                }
            c += k.Sw / den;
        }
        void marginaliseInto(std::array<double, Ny * Ny>& out) const {   // free offset
            for (int j = 0; j < Ny; ++j)
                for (int m = 0; m < Ny; ++m)
                    out[j * Ny + m] += F[j * Ny + m]
                        - (c > 0.0 ? (b[0][j] * b[0][m] + b[1][j] * b[1][m]) / c : 0.0);
        }
    };

    std::array<double, Ny * Ny> FL{};
    { Group g; g.fold(rubin, 0.0); g.marginaliseInto(FL); }
    for (int v = 0; v < NAVAR; ++v) {
        std::vector<Group> groups(v == AV_W ? 1 : (v == AV_N ? 2 : AST_MAX_SEASONS));
        for (const auto& [day, k] : romanDays)
            groups[v == AV_W ? 0 : (v == AV_N ? k.roll : k.season)].fold(k, AST_SIGC[v]);
        std::array<double, Ny * Ny> FR{};
        for (const auto& g : groups) g.marginaliseInto(FR);
        for (int jm = 0; jm < Ny * Ny; ++jm) {
            co.FBV[v][SRUBIN][jm] = FL[jm];
            co.FBV[v][SROMAN][jm] = FR[jm];
            co.FBV[v][SJOINT][jm] = FL[jm] + FR[jm];
        }
    }
    // The main matrices are variant W.
    for (int q = 0; q < NSURV; ++q)
        for (int j = 0; j < Ny; ++j)
            for (int m = 0; m < Ny; ++m)
                gsl_matrix_set(co.inputB[q].get(), j, m, co.FBV[AV_W][q][j * Ny + m]);

    for (int q = 0; q < NSURV; ++q) {
        co.okB[q] = (co.nepochA[q] >= Ny) ? invert_matrix(co, 1, q) : 0;
    }

    // Variants N and P, through the same normalised inversion.
    {
        static const std::vector<int> kAllAst = {0, 1, 2, 3};
        gsl_matrix* Fin  = gsl_matrix_alloc(Ny, Ny);
        gsl_matrix* Fout = gsl_matrix_alloc(Ny, Ny);
        for (int v = 0; v < NAVAR; ++v)
            for (int q = 0; q < NSURV; ++q) {
                co.okBV[v][q] = 0; co.condBV[v][q] = -1.0;
                for (int k = 0; k < Ny; ++k) co.ErbV[v][q][k] = -1.0;
                if (v == AV_W or co.nepochA[q] < Ny) continue;      // W is filled in ErrorCal
                for (int j = 0; j < Ny; ++j)
                    for (int m = 0; m < Ny; ++m) gsl_matrix_set(Fin, j, m, co.FBV[v][q][j * Ny + m]);
                co.okBV[v][q] = invertNormalized(Fin, Fout, kAllAst, co.condBV[v][q], nullptr);
                if (co.okBV[v][q])
                    for (int k = 0; k < Ny; ++k)
                        co.ErbV[v][q][k] = std::sqrt(std::fabs(gsl_matrix_get(Fout, k, k)));
            }
        gsl_matrix_free(Fin);
        gsl_matrix_free(Fout);
    }


    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, co.inverB[SJOINT].get(), co.inputB[SJOINT].get(), 0.0, co.summB.get());
}

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Error parameters Calculation               //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void ErrorCal(covarian & co, lens & l , source & s){
  double corr1;

  // Per-survey 1-sigma forecasts. -1.0 marks a partition that is not characterizable (too few
  // epochs or singular matrix) and, for the photometric sigmas, a parameter outside the
  // partition's active subset. Downstream code must test okA/okB and skip negative sigmas.
  for (int q = 0; q < NSURV; ++q) {
      for (int k = 0; k < Nx; ++k) co.Era[q][k] = -1.0;
      if (co.okA[q]) {
          for (int k : activePhotParams(q, co.nepochA[SRUBIN], co.nepochA[SROMAN])) {
              co.Era[q][k] = std::sqrt(std::fabs(gsl_matrix_get(co.inverA[q].get(), k, k)));
          }
      }
      for (int k = 0; k < Ny; ++k) {
          co.Erb[q][k] = co.okB[q] ? std::sqrt(std::fabs(gsl_matrix_get(co.inverB[q].get(), k, k)))
                                   : -1.0;
      }
  }

  // Fractional 1-sigma on the lens mass, per survey: Ml = tetE / (kappa * piE), kappa = 8.144
  // mas/Msun, so the fractional errors add in quadrature. piE has a photometric (Era[3]) and an
  // astrometric (Erb[3]) route; take the better available one, and report no mass if either
  // ingredient is missing (-1.0 means "not measured", never a measurement).
  for (int q = 0; q < NSURV; ++q) {
      const double rp1 = (co.Era[q][3] >= 0.0) ? co.Era[q][3] / (std::fabs(l.piE)  + eps) : -1.0;
      const double rp2 = (co.Erb[q][3] >= 0.0) ? co.Erb[q][3] / (std::fabs(l.piE)  + eps) : -1.0;
      const double rp  = (rp1 >= 0.0 and rp2 >= 0.0) ? MIN(rp1, rp2) : std::max(rp1, rp2);
      const double rt  = (co.Erb[q][0] >= 0.0) ? co.Erb[q][0] / (std::fabs(l.tetE) + eps) : -1.0;
      co.relMl[q] = (rp >= 0.0 and rt >= 0.0) ? std::sqrt(rp * rp + rt * rt) : -1.0;
  }

  // Astrometric noise variants: W is the main result; N and P were inverted in FisherM. Lens-mass
  // errors follow the same rule, with each variant's own tetE and astrometric piE.
  for (int q = 0; q < NSURV; ++q) {
      co.okBV[AV_W][q] = co.okB[q];
      co.condBV[AV_W][q] = co.condB[q];
      for (int k = 0; k < Ny; ++k) co.ErbV[AV_W][q][k] = co.Erb[q][k];
      for (int v = 0; v < NAVAR; ++v) {
          const auto& E = co.ErbV[v][q];
          const double rp1 = (co.Era[q][3] >= 0.0) ? co.Era[q][3] / (std::fabs(l.piE) + eps) : -1.0;
          const double rp2 = (E[3] >= 0.0) ? E[3] / (std::fabs(l.piE) + eps) : -1.0;
          const double rp  = (rp1 >= 0.0 and rp2 >= 0.0) ? MIN(rp1, rp2) : std::max(rp1, rp2);
          const double rt  = (E[0] >= 0.0) ? E[0] / (std::fabs(l.tetE) + eps) : -1.0;
          co.relMlV[v][q] = (rp >= 0.0 and rt >= 0.0) ? std::sqrt(rp * rp + rt * rt) : -1.0;
      }
  }

  // Joint event summary: resu[] keeps its index layout for the output columns. Per-survey
  // numbers live in Era[]/Erb[].
  const auto& Era = co.Era[SJOINT];
  const auto& Erb = co.Erb[SJOINT];

  corr1 = double(gsl_matrix_get(co.inverA[SJOINT].get(), 1, 4) / std::fabs(Era[1] * Era[4])); //correlation tE  && xi

  co.resu[0]  = double(Era[0] / (std::fabs(l.u0)    + eps));
  co.resu[1]  = double(Era[1] / (std::fabs(l.tE)    + eps));
  co.resu[2]  = double(Era[2] / (std::fabs(s.fb[0]) + eps));
  co.resu[3]  = double(Era[3] / (std::fabs(l.piE)   + eps));
  co.resu[4]  = double(Era[4] / (std::fabs(s.xi)    + eps));
  co.resu[5]  = double(Erb[0] / (std::fabs(l.tetE)  + eps));
  co.resu[6]  = double(Erb[1] / (std::fabs(s.mus1)  + eps));
  co.resu[7]  = double(Erb[2] / (std::fabs(s.mus2)  + eps));
  co.resu[8]  = double(Erb[3] / (std::fabs(l.piE)   + eps));

  // piE is measured photometrically (resu[3]) and astrometrically (resu[8]). Take the better of
  // the routes that produced a measurement; an unguarded MIN would pick the -1.0 sentinel over a
  // real sigma. Same rule as relMl[].
  if      (co.resu[3] < 0.0 and co.resu[8] < 0.0) co.resu[3] = -1.0;
  else if (co.resu[3] < 0.0)                      co.resu[3] = co.resu[8];
  else if (co.resu[8] >= 0.0)                     co.resu[3] = MIN(co.resu[3], co.resu[8]);

  //Lens Mass: Ml = tetE/(kappa*piE); fractional errors add in quadrature. -1 if either is missing.
  co.resu[9]  = (co.resu[3] >= 0.0 and co.resu[5] >= 0.0)
              ? std::sqrt(co.resu[3] * co.resu[3] + co.resu[5] * co.resu[5]) : -1.0;
  co.resu[10] = std::fabs(co.resu[9] * (s.Ds - l.Dl) / s.Ds); //Dl

  co.f1 = co.resu[5] * co.resu[5] + co.resu[1] * co.resu[1] + (Era[4] * std::tan(s.xi)) * (Era[4] * std::tan(s.xi))
    - 2.0 * co.resu[1] * Era[4] * std::tan(s.xi) * corr1;
  co.f2 = co.resu[5] * co.resu[5] + co.resu[1] * co.resu[1] + (Era[4] / std::tan(s.xi)) * (Era[4] / std::tan(s.xi))
    - 2.0 * co.resu[1] * Era[4] / std::tan(s.xi) * corr1;

  co.sigmul1 = std::sqrt(Erb[1] * Erb[1] + l.murel * l.murel * std::cos(s.xi) * std::cos(s.xi) * std::fabs(co.f1));
  co.sigmul2 = std::sqrt(Erb[2] * Erb[2] + l.murel * l.murel * std::sin(s.xi) * std::sin(s.xi) * std::fabs(co.f2));

  co.resu[11] = double(co.sigmul1 / (std::fabs(l.mul1) + eps)); //mu_lens_1
  co.resu[12] = double(co.sigmul2 / (std::fabs(l.mul2) + eps)); //mu_lens_2
  co.resu[13] = std::sqrt(co.sigmul1 * co.sigmul1 + co.sigmul2 * co.sigmul2) / (std::fabs(l.mul) + eps); //mu_lens
  co.resu[14] = std::sqrt(Erb[1]  * Erb[1]  +  Erb[2] * Erb[2])  / (std::fabs(s.mus) + eps); // mu_source
}
