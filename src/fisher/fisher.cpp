// FisherM (information matrices by finite differences) and ErrorCal (sigmas from the inverses).
#include "fisher/fisher.h"
#include "events/lightcurve.h"
#include "surveys/noise.h"
#include <cstring>


void FisherM(source & s, lens & l, astromet & as,  covarian & co, int ndw)
{
    co.flagi =+ 1;

    // The parameter point. Every perturbation below is undone by restoring the saved value (not by
    // subtracting the step, which can leave it an ulp off), so FisherM returns it bit for bit.
    const auto point = [&] {
        return std::array<double, 12>{l.u0, l.tE, l.piE, s.xi, l.t0, s.fb[0], s.fb[1],
                                      s.mbs[0], s.mbs[1], l.tetE, s.mus1, s.mus2};
    };
    const auto point0 = point();

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
    // and fb1 (index 7) use the telescope-keyed steps bb[] set per telescope below.

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
    // The model is evaluated one observer at a time over all of its epochs, so the epochs are split
    // by telescope here (stored order kept) and the results put back at their stored index.
    std::array<int, 2> nTele{0, 0};
    std::array<std::vector<int>, 2>    epo;  //stored index of each of the telescope's epochs
    std::array<std::vector<double>, 2> tim;  //and its time
    for (int i = 0; i < ndw; ++i) {
        const int tt = int(l.tele[i]);//telescope[0,1] LSST, ELT
        CHECK(tt == 0 or tt == 1);
        nTele[tt] += 1;
        epo[tt].push_back(i);
        tim[tt].push_back(l.timn[i]);
    }

    // dm[i][j] = d(model magnitude of epoch i)/d(theta_j), from the central stencil sig. Each
    // derivative is needed once per epoch, however many F_jk it enters, so the reference model and
    // each perturbed model are evaluated once over a telescope's epochs (perturbation outer, epoch
    // inner). Of the nine parameters only u0, tE, piE, xi and t0 move the source relative to the
    // lens; fb and mbs enter m = mbs - 2.5 log10(fb A + 1 - fb) through the flux alone, so their
    // perturbed magnitudes come from the reference A. That is 1 + 2*5 = 11 model evaluations per
    // epoch; re-deriving the second factor of each of the 45 F_jk would take 1 + 2*9 + 2*45 = 109.
    std::vector<std::array<double, Nx>> dm(ndw);
    std::vector<ModelPoint> ref, pert;
    std::vector<double> lgf;  //2.5 log10(fb A + 1 - fb) of the reference model, per epoch
    const std::array<std::pair<int, double*>, 5> lcPar{{
        {0, &l.u0}, {1, &l.tE}, {3, &l.piE}, {4, &s.xi}, {5, &l.t0}}};

    for (int tt = 0; tt < 2; ++tt) {
        const int n = nTele[tt];
        if (n < kMinTeleEpochs) continue;
        const int*    ix = epo[tt].data();
        const double* t  = tim[tt].data();
        ref.resize(n);  pert.resize(n);  lgf.resize(n);
        const double fb = s.fb[tt], mbs = s.mbs[tt];  //this telescope's flux pair; never written here
        CHECK(fb > 0.0);
        CHECK(fb <= 1.0);

        // Unperturbed model magnitudes, recomputed under the current observer configuration rather
        // than taken from the cached l.magn[]. Each derivative is (model(theta + Delta) - reference)
        // / Delta, so the reference must share the observing geometry (as.satScale). When the two
        // differ, their difference leaks into every derivative as a constant dm/Delta; it cancels in
        // symmetric stencils but not in the cross-telescope rows, fb's outer bins, or any one-sided
        // stencil. Every model here is evaluated with the observer that produced these data.
        evaluateModel(s, l, as, tt, t, n, ref.data());
        for (int m = 0; m < n; ++m) lgf[m] = 2.5 * std::log10(ref[m].A * fb + 1.0 - fb);

        // u0, tE, piE, xi, t0: one new model per perturbation.
        for (const auto& [j, par] : lcPar) {
            const double keep = *par;
            for (int h = 0; h < 2; ++h) {
                const double diff = co.Delta1[j] * sig[h];
                *par = keep + diff;
                CHECK(l.tE > 0.0);
                CHECK(l.piE > 0.0);
                CHECK(l.u0 != 0.0);
                CHECK(diff != 0.0);
                evaluateModel(s, l, as, tt, t, n, pert.data());
                *par = keep;
                for (int m = 0; m < n; ++m) {
                    const double magw = mbs - 2.5 * std::log10(pert[m].A * fb + 1.0 - fb);
                    const double der  = (magw - (mbs - lgf[m])) / diff;
                    double& d = dm[ix[m]][j];
                    d = (h == 0) ? der : (d + der) * 0.5;
                }
            }
        }

        // This telescope's flux pair, fb (index 2 or 7) and mbs (6 or 8), from the reference A.
        // The fb step is binned on THIS telescope's blend fraction s.fb[tt] (not always s.fb[0]):
        // bins and steps are chosen so fb[tt] + bb never leaves [0,1], which holds only if bb is
        // computed from the fb it perturbs. Roman's F146 blend fraction is routinely higher than
        // Rubin's r-band one.
        const int jfb = (tt == 0) ? 2 : 7, jmbs = (tt == 0) ? 6 : 8;
        std::array<double, 2> bb;
        if (fb < FB_BIN_LO)      bb = {+FB_STEP_SMALL, +FB_STEP_LARGE};
        else if (fb < FB_BIN_HI) bb = {-FB_STEP_SMALL, +FB_STEP_SMALL};
        else                     bb = {-FB_STEP_SMALL, -FB_STEP_LARGE};

        // Apply the plateau scaling and clamp so the step cannot leave the physical range [0,1]
        // (the bin edges guarantee this only at unscaled sizes). Note the outer fb bins are two
        // same-signed (one-sided) steps, so fb carries an O(h) bias there; kFDStepScale makes it
        // negligible.
        {
            const double fscale = kFDStepScale * co.deltaScale[jfb];
            for (int b = 0; b < 2; ++b) {
                double step = bb[b] * fscale;
                const double lo = 1.0e-6 - fb;        //keeps fb strictly above 0
                const double hi = 1.0 - 1.0e-6 - fb;  //keeps fb strictly below 1
                if (step < lo) step = lo;
                if (step > hi) step = hi;
                bb[b] = step;
            }
        }

        for (int h = 0; h < 2; ++h) {
            const double dfb  = bb[h];
            const double dmbs = co.Delta1[jmbs] * sig[h];
            const double fbp  = fb + dfb, mbsp = mbs + dmbs;  //the perturbed values
            CHECK(dfb != 0.0);
            CHECK(fbp > 0.0);
            CHECK(fbp <= 1.0);
            CHECK(dmbs != 0.0);
            for (int m = 0; m < n; ++m) {
                const double mag0  = mbs - lgf[m];
                const double magfb = mbs - 2.5 * std::log10(ref[m].A * fbp + 1.0 - fbp);
                const double magmb = mbsp - lgf[m];
                const double derfb = (magfb - mag0) / dfb;
                const double dermb = (magmb - mag0) / dmbs;
                double& dfbm = dm[ix[m]][jfb];
                double& dmbm = dm[ix[m]][jmbs];
                dfbm = (h == 0) ? derfb : (dfbm + derfb) * 0.5;
                dmbm = (h == 0) ? dermb : (dmbm + dermb) * 0.5;
            }
        }

        // The other telescope's flux pair does not enter these epochs: its derivative is exactly 0.
        for (int m = 0; m < n; ++m) dm[ix[m]][(tt == 0) ? 7 : 2] = dm[ix[m]][(tt == 0) ? 8 : 6] = 0.0;
    }

    // Accumulate the information sum F_jk = sum_i (dm_i/dtheta_j)(dm_i/dtheta_k)/sigma_i^2 (lower
    // triangle, epochs in stored order). Each epoch feeds the joint matrix and exactly one
    // single-survey matrix, so F[SJOINT] == F[SRUBIN] + F[SROMAN] holds element by element.
    std::array<std::array<double, Nx * Nx>, NSURV> FA{};
    for (int i = 0; i < ndw; ++i) {//data
        const int tt = int(l.tele[i]);
        if (nTele[tt] < kMinTeleEpochs) continue;
        const int surv = surveyOfTele(tt); // which single-survey matrix this epoch also feeds
        co.nepochA[SJOINT] += 1;
        co.nepochA[surv]   += 1;
        CHECK(l.errm[i] > 0.0);
        const auto& d = dm[i];
        for (int j = 0; j < Nx; ++j) {
            for (int k = 0; k <= j; ++k) {
                const double contrib = d[j] * d[k] / (l.errm[i] * l.errm[i]);
                FA[SJOINT][j * Nx + k] += contrib;
                FA[surv][j * Nx + k]   += contrib;
            }
        }
    }//end of data for
    for (int q = 0; q < NSURV; ++q)
        for (int j = 0; j < Nx; ++j)
            for (int k = 0; k <= j; ++k) gsl_matrix_set(co.inputA[q].get(), j, k, FA[q][j * Nx + k]);

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
    // (AST_SIGC in config/parameters.h, co.astSigc), from ONE pass over the epochs. ----
    // Per epoch the derivatives of each sky coordinate c (0 = x, 1 = y) are computed once.
    // Sums kept per coordinate: S_w = sum w, S_wd[c][j] = sum w d_cj, S_wdd[c][jk]. Rubin's epochs
    // go into one white group; Roman's go into day blocks (the sigma_c correlation unit), tagged
    // with season and roll so the variants can group them differently.
    // Each coordinate carries one extra parameter, slot Ny: the light centroid of the telescope's
    // unresolved neighbours in that coordinate (s.blendOff), unknown to the analyst. Its derivative is
    // dcentroid/dblendOff = fn / (fs A + 1 - fs), largest at baseline and falling as the source
    // brightens; it is marginalised below, per telescope and per coordinate (it appears in no other
    // telescope's or coordinate's data), when AST_BLEND_OFFSET_FREE.
    constexpr int D = Ny + 1;
    struct AstSums {
        double Sw = 0.0;
        std::array<std::array<double, D>, 2>     Swd{};
        std::array<std::array<double, D * D>, 2> Swdd{};
        int season = -1, roll = -1;
        void add(double w, const std::array<double, D>& dx, const std::array<double, D>& dy) {
            Sw += w;
            for (int c = 0; c < 2; ++c) {
                const auto& d = c ? dy : dx;
                for (int j = 0; j < D; ++j) {
                    Swd[c][j] += w * d[j];
                    for (int k = 0; k < D; ++k) Swdd[c][j * D + k] += w * d[j] * d[k];
                }
            }
        }
    };
    AstSums rubin;
    std::map<long, AstSums> romanDays;

    // Per epoch i, ddx[i][j] and ddy[i][j]: the derivatives of the centroid's x and y, by the
    // central stencil sig for all four astrometric parameters, and in slot Ny the neighbours' share
    // h. As for the photometry, the reference and each of the 2*Ny perturbed models are evaluated
    // once over a telescope's epochs: 1 + 2*4 = 9 model evaluations per epoch. Every epoch enters,
    // whatever its telescope's photometric count.
    std::vector<std::array<double, D>> ddx(ndw), ddy(ndw);
    std::vector<double> sx0, sy0;  //reference centroid without the neighbours' share, per epoch
    const std::array<double*, Ny> astPar{&l.tetE, &s.mus1, &s.mus2, &l.piE};

    for (int tt = 0; tt < 2; ++tt) {
        const int n = nTele[tt];
        if (n == 0) continue;
        const int*    ix = epo[tt].data();
        const double* t  = tim[tt].data();
        ref.resize(n);  pert.resize(n);  sx0.resize(n);  sy0.resize(n);

        // The neighbours' share of the centroid, fn blendOff / (fs A + 1 - fs), and its derivative with
        // respect to blendOff, at impact parameter u. The parameter derivatives below are taken on the
        // centroid WITHOUT this share: its time profile is the magnification curve, which belongs to the
        // photometric matrix, and read here (with u0, tE, t0 held fixed) it would hand piE information
        // the split cannot marginalise. The share enters only through its unknown amplitude, slot Ny.
        const double fn = std::max(0.0, 1.0 - s.fb[tt] - s.fLens[tt]);
        auto nbrShare = [&](double u, double& bx, double& by) {
            const double u2 = u * u;
            const double h  = fn / (s.fb[tt] * (u2 + 2.0) / std::sqrt(u2 * (u2 + 4.0)) + 1.0 - s.fb[tt]);
            bx = h * s.blendOff[tt][0];  by = h * s.blendOff[tt][1];
            return h;
        };

        // Reference and perturbations use the observer that produced these data.
        evaluateModel(s, l, as, tt, t, n, ref.data());
        for (int m = 0; m < n; ++m) {
            double nbx0, nby0;
            ddx[ix[m]][Ny] = ddy[ix[m]][Ny] = nbrShare(ref[m].u, nbx0, nby0);
            sx0[m] = ref[m].pos1c - nbx0;
            sy0[m] = ref[m].pos2c - nby0;
        }

        for (int j = 0; j < Ny; ++j) {
            const double keep = *astPar[j];
            for (int h = 0; h < 2; ++h) {
                const double diff = co.Delta2[j] * sig[h];
                *astPar[j] = keep + diff;
                CHECK(l.tetE > 0.0);
                CHECK(l.piE > 0.0);
                CHECK(diff != 0.0);
                evaluateModel(s, l, as, tt, t, n, pert.data());
                *astPar[j] = keep;
                // A null derivative (epoch at t0 or far from the peak) is legitimate and contributes nothing.
                for (int m = 0; m < n; ++m) {
                    double nbx, nby;
                    nbrShare(pert[m].u, nbx, nby);
                    const double derx = (pert[m].pos1c - nbx - sx0[m]) / diff;
                    const double dery = (pert[m].pos2c - nby - sy0[m]) / diff;
                    double& dx = ddx[ix[m]][j];
                    double& dy = ddy[ix[m]][j];
                    dx = (h == 0) ? derx : (dx + derx) * 0.5;
                    dy = (h == 0) ? dery : (dy + dery) * 0.5;
                }
            }
        }
    }

    for (int i = 0; i < ndw; ++i) {
        CHECK(l.erra[i] > 0.0);
        const double w = 1.0 / (l.erra[i] * l.erra[i]);       // per coordinate
        if (int(l.tele[i]) == 0) {
            rubin.add(w, ddx[i], ddy[i]);
        } else {
            CHECK(l.rseas[i] >= 0 and l.rseas[i] < AST_MAX_SEASONS);
            CHECK(l.rroll[i] == 0 or l.rroll[i] == 1);
            AstSums& d = romanDays[long(std::floor(l.timn[i]))];
            if (d.season < 0) { d.season = l.rseas[i]; d.roll = l.rroll[i]; }
            CHECK(d.season == l.rseas[i] and d.roll == l.rroll[i]);
            d.add(w, ddx[i], ddy[i]);
        }
    }

    // A frame group: per coordinate, the sum of block matrices and offset-vectors b, and c.
    using PerCoord = std::array<std::array<double, D * D>, 2>;
    struct Group {
        PerCoord F{};
        std::array<std::array<double, D>, 2> b{};
        double c = 0.0;
        void fold(const AstSums& k, double sc) {          // one block, correlated at sigma_c = sc
            const double s2 = sc * sc, den = 1.0 + s2 * k.Sw;
            for (int cc = 0; cc < 2; ++cc)
                for (int j = 0; j < D; ++j) {
                    b[cc][j] += k.Swd[cc][j] / den;
                    for (int m = 0; m < D; ++m)
                        F[cc][j * D + m] += k.Swdd[cc][j * D + m] - s2 * k.Swd[cc][j] * k.Swd[cc][m] / den;
                }
            c += k.Sw / den;
        }
        void marginaliseInto(PerCoord& out) const {       // free offset
            for (int cc = 0; cc < 2; ++cc)
                for (int jm = 0; jm < D * D; ++jm)
                    out[cc][jm] += F[cc][jm] - (c > 0.0 ? b[cc][jm / D] * b[cc][jm % D] / c : 0.0);
        }
    };
    // One telescope's per-coordinate matrices -> its Ny x Ny information on the physical parameters,
    // with the blend offset marginalised (Schur complement of slot Ny). raw = that telescope's
    // sum w h^2 per coordinate before any marginalisation: when what is left after the offsets is a
    // vanishing fraction of it, h was constant (no magnification while observed), the blend offset is
    // indistinguishable from the reference position and carries no cross-information, so nothing is
    // subtracted.
    auto physical = [](const PerCoord& G, const std::array<double, 2>& raw, std::array<double, Ny * Ny>& out) {
        for (int cc = 0; cc < 2; ++cc) {
            const double hh = G[cc][Ny * D + Ny];
            const bool   fit = AST_BLEND_OFFSET_FREE and raw[cc] > 0.0 and hh > 1.0e-9 * raw[cc];
            for (int j = 0; j < Ny; ++j)
                for (int m = 0; m < Ny; ++m)
                    out[j * Ny + m] += G[cc][j * D + m] - (fit ? G[cc][j * D + Ny] * G[cc][m * D + Ny] / hh : 0.0);
        }
    };

    std::array<double, Ny * Ny> FL{};
    {
        Group g; g.fold(rubin, 0.0);
        PerCoord G{}; g.marginaliseInto(G);
        physical(G, {rubin.Swdd[0][Ny * D + Ny], rubin.Swdd[1][Ny * D + Ny]}, FL);
    }
    std::array<double, 2> rawR{};
    for (const auto& [day, k] : romanDays)
        for (int cc = 0; cc < 2; ++cc) rawR[cc] += k.Swdd[cc][Ny * D + Ny];
    for (int v = 0; v < NAVAR; ++v) {
        std::vector<Group> groups(v == AV_W ? 1 : 2);          // W: one frame; N, P: one per roll
        for (const auto& [day, k] : romanDays)
            groups[v == AV_W ? 0 : k.roll].fold(k, co.astSigc[v]);
        PerCoord G{};
        for (const auto& g : groups) g.marginaliseInto(G);
        std::array<double, Ny * Ny> FR{};
        physical(G, rawR, FR);
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

    const auto point1 = point();
    CHECK(std::memcmp(point1.data(), point0.data(), sizeof(point1)) == 0);  //bit for bit, NaN included
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
