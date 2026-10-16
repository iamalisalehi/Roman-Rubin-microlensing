// FisherM (information matrices by finite differences) and ErrorCal (sigmas from the inverses).
#include "fisher/fisher.h"
#include "events/lightcurve.h"
#include "events/vbm_model.h"
#include "surveys/noise.h"
#include <cstring>


void FisherM(source & s, lens & l, astromet & as,  covarian & co, int ndw)
{
    co.flagi =+ 1;

    // The parameter point. Every perturbation below is undone by restoring the saved value (not by
    // subtracting the step, which can leave it an ulp off), so FisherM returns it bit for bit.
    const auto point = [&] {
        return std::array<double, 13>{l.u0, l.tE, l.piE, s.xi, l.t0, s.fb[0], s.fb[1],
                                      s.mbs[0], s.mbs[1], l.tetE, s.mus1, s.mus2, s.rho};
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
    co.Delta1[9] = FD_STEP_LNRHO; // ln rho -- the source radius in Einstein radii, perturbed as rho * exp(Delta)

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
        co.dchiFS[q] = -1.0;
        co.rhoOn[q] = 0;
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
    // inner). Of the nine parameters of a point source only u0, tE, piE, xi and t0 move the source
    // relative to the lens; fb and mbs enter m = mbs - 2.5 log10(fb A + 1 - fb) through the flux alone,
    // so their perturbed magnitudes come from the reference A. That is 1 + 2*5 = 11 model evaluations
    // per epoch; re-deriving the second factor of each of the 45 F_jk would take 1 + 2*9 + 2*45 = 109.
    // The tenth, ln rho, changes the shape of the magnification near the lens and costs two more
    // evaluations per epoch, but only for a partition that sees the finite source (below).
    std::vector<std::array<double, Nx>> dm(ndw);
    std::vector<ModelPoint> ref, pert;
    std::array<std::vector<double>, 2> lgfT;  //2.5 log10(fb A + 1 - fb) of the reference model, per epoch
    const std::array<std::pair<int, double*>, 5> lcPar{{
        {0, &l.u0}, {1, &l.tE}, {3, &l.piE}, {4, &s.xi}, {5, &l.t0}}};

    // Whether the model has a source size to fit: a uniform disc with a radius. Under a point source
    // (--source-model pspl) ln rho is not a parameter and dchiFS is not computed.
    const bool hasDisc = (vbmSourceModel() == SourceModel::ESPL) and s.rho > 0.0;
    std::array<double, 2> dchiT{0.0, 0.0};  //each telescope's share of dchiFS

    for (int tt = 0; tt < 2; ++tt) {
        const int n = nTele[tt];
        if (n < kMinTeleEpochs) continue;
        const int*    ix = epo[tt].data();
        const double* t  = tim[tt].data();
        ref.resize(n);  pert.resize(n);
        std::vector<double>& lgf = lgfT[tt];
        lgf.resize(n);
        const double fb = s.fb[tt], mbs = s.mbs[tt];  //this telescope's flux pair; never written here
        CHECK(fb > 0.0);
        CHECK(fb <= 1.0);

        // Unperturbed model magnitudes, recomputed under the current observer configuration rather
        // than taken from the cached l.magn[]. Each derivative is (model(theta + Delta) - reference)
        // / Delta, so the reference must share the observing geometry (as.satScale). When the two
        // differ, their difference leaks into every derivative as a constant dm/Delta; it cancels in
        // symmetric stencils but not in the cross-telescope rows, fb's outer bins, or any one-sided
        // stencil. Every model here is evaluated with the observer that produced these data.
        evaluateModel(s, l, as, tt, t, n, false, ref.data());
        for (int m = 0; m < n; ++m) lgf[m] = 2.5 * std::log10(ref[m].A * fb + 1.0 - fb);

        // The finite-source delta-chi2 of this telescope's epochs: the chi-squared by which they tell the disc
        // from a point source at the true parameters, sum ((m_disc - m_point) / sigma_m)^2. The two models
        // share the track (the source radius does not enter it), so the point source's magnification at an
        // epoch is magnifOf of the disc model's own u, bit for bit what the library gives for rho = 0
        // (tests/fisher_fixture.cpp --finite-source asserts it), and no second light curve is needed.
        if (hasDisc) {
            double d = 0.0;
            for (int m = 0; m < n; ++m) {
                CHECK(l.errm[ix[m]] > 0.0);
                const double dmag = 2.5 * std::log10(magnifOf(ref[m].u) * fb + 1.0 - fb) - lgf[m];
                const double r = dmag / l.errm[ix[m]];
                d += r * r;
            }
            dchiT[tt] = d;
        }

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
                evaluateModel(s, l, as, tt, t, n, false, pert.data());
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

    // Is the finite source seen? dchiFS of a partition is the sum of its telescopes' shares (the joint is
    // Rubin + Roman; a telescope left out above has none). Partitions at or above the bar DCHI_FS try ln rho as a
    // parameter; the joint is at or above both single surveys, so its derivative is taken whenever either tries.
    // A partition below the bar holds rho at its true value. Whether a trying partition keeps ln rho is decided
    // after its matrix is inverted, below (SIGLNRHO_MAX).
    if (hasDisc) {
        const bool useL = nTele[0] >= kMinTeleEpochs, useR = nTele[1] >= kMinTeleEpochs;
        if (useL)         co.dchiFS[SRUBIN] = dchiT[0];
        if (useR)         co.dchiFS[SROMAN] = dchiT[1];
        if (useL or useR) co.dchiFS[SJOINT] = (useL ? dchiT[0] : 0.0) + (useR ? dchiT[1] : 0.0);
        for (int q = 0; q < NSURV; ++q) co.rhoOn[q] = (co.dchiFS[q] >= co.dchiFSMin) ? 1 : 0;
    }

    // ln rho: one perturbed model per step and telescope, as for u0 and tE above; none at all (for the great
    // majority of events) unless some partition tries it. The source radius is perturbed multiplicatively,
    // rho -> rho exp(+-Delta), so the step is in ln rho. The uniform disc's table is bilinear in ln rho (rows
    // 0.092 apart), so for a step inside a cell the derivative is the slope of the cell holding rho.
    if (co.rhoOn[SJOINT]) {
        const double keep = s.rho;
        for (int tt = 0; tt < 2; ++tt) {
            const int n = nTele[tt];
            if (n < kMinTeleEpochs) continue;
            const int*    ix = epo[tt].data();
            const double* t  = tim[tt].data();
            const std::vector<double>& lgf = lgfT[tt];
            pert.resize(n);
            const double fb = s.fb[tt], mbs = s.mbs[tt];
            for (int h = 0; h < 2; ++h) {
                const double diff = co.Delta1[9] * sig[h];
                s.rho = keep * std::exp(diff);
                CHECK(std::isfinite(s.rho) and s.rho > 0.0);
                CHECK(diff != 0.0);
                evaluateModel(s, l, as, tt, t, n, false, pert.data());
                s.rho = keep;
                for (int m = 0; m < n; ++m) {
                    const double magw = mbs - 2.5 * std::log10(pert[m].A * fb + 1.0 - fb);
                    const double der  = (magw - (mbs - lgf[m])) / diff;
                    double& d = dm[ix[m]][9];
                    d = (h == 0) ? der : (d + der) * 0.5;
                }
            }
        }
    }

    // Accumulate the information sum F_jk = sum_i (dm_i/dtheta_j)(dm_i/dtheta_k)/sigma_i^2 (lower
    // triangle, epochs in stored order). Each epoch feeds the joint matrix and exactly one
    // single-survey matrix, so F[SJOINT] == F[SRUBIN] + F[SROMAN] holds element by element.
    // ln rho's derivative is zero unless some partition tries it (rhoOn[SJOINT] is the union), and its rows then
    // add nothing: leave them out of the sum.
    const int nj = co.rhoOn[SJOINT] ? Nx : Nx - 1;
    std::array<std::array<double, Nx * Nx>, NSURV> FA{};
    for (int i = 0; i < ndw; ++i) {//data
        const int tt = int(l.tele[i]);
        if (nTele[tt] < kMinTeleEpochs) continue;
        const int surv = surveyOfTele(tt); // which single-survey matrix this epoch also feeds
        co.nepochA[SJOINT] += 1;
        co.nepochA[surv]   += 1;
        CHECK(l.errm[i] > 0.0);
        const auto& d = dm[i];
        for (int j = 0; j < nj; ++j) {
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
        co.okA[q] = (co.nepochA[q] >= static_cast<int>(activePhotParams(co, q).size()))
                        ? invert_matrix(co, 0, q) : 0;

        // A partition that passed the dchiFS pre-filter keeps ln rho only if the data MEASURE the source's size:
        // the matrix inverted with it, and sigma(ln rho) is at most SIGLNRHO_MAX (config/parameters.h). Otherwise
        // the finite-source term is degenerate with u0 and no linear forecast of it holds, and the partition is
        // forecast as a point-source event would be fitted, with rho held at its true value: the nine-parameter
        // inversion, which is what an event below the pre-filter gets, bit for bit (the matrix entries of the
        // other parameters are the same ones).
        if (co.rhoOn[q]) {
            const bool measured = co.okA[q]
                and std::sqrt(std::fabs(gsl_matrix_get(co.inverA[q].get(), 9, 9))) <= co.sigLnRhoMax;
            if (!measured) {
                co.rhoOn[q] = 0;
                co.okA[q] = (co.nepochA[q] >= static_cast<int>(activePhotParams(co, q).size()))
                                ? invert_matrix(co, 0, q) : 0;
            }
        }
    }

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///                   Astrometry                                            ///
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    // The astrometric matrix holds every photometric parameter fixed at its true value, and so the source
    // radius s.rho with them, although rho = theta* / thetaE would move with thetaE: the centroid's
    // derivatives with respect to tetE, mus1, mus2 and piE are taken at fixed (u0, tE, t0, xi, fb, mbs, rho).
    // Whatever the light curve says about rho is the photometric matrix's to use (ln rho, index 9 above).
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
        // respect to blendOff, at the model's magnification A of the epoch. The parameter derivatives below
        // are taken on the centroid WITHOUT this share: its time profile is the magnification curve, which
        // belongs to the photometric matrix, and read here (with u0, tE, t0 held fixed) it would hand piE
        // information the split cannot marginalise. The share enters only through its unknown amplitude,
        // slot Ny.
        // For a point source A is magnifOf(u) bit for bit, and the share keeps the form it has always had
        // (fs (u^2 + 2) / sqrt(u^2 (u^2 + 4)) rounds differently from fs A), so that point-source results do
        // not move in the last digit; for a finite source it is the model's A, which u alone no longer gives.
        const double fn = std::max(0.0, 1.0 - s.fb[tt] - s.fLens[tt]);
        const bool   finite = (vbmSourceModel() == SourceModel::ESPL);
        auto nbrShare = [&](const ModelPoint& p, double& bx, double& by) {
            const double u2   = p.u * p.u;
            const double flux = finite ? s.fb[tt] * p.A + 1.0 - s.fb[tt]
                                       : s.fb[tt] * (u2 + 2.0) / std::sqrt(u2 * (u2 + 4.0)) + 1.0 - s.fb[tt];
            const double h = fn / flux;
            bx = h * s.blendOff[tt][0];  by = h * s.blendOff[tt][1];
            return h;
        };

        // Reference and perturbations use the observer that produced these data.
        evaluateModel(s, l, as, tt, t, n, true, ref.data());
        for (int m = 0; m < n; ++m) {
            double nbx0, nby0;
            ddx[ix[m]][Ny] = ddy[ix[m]][Ny] = nbrShare(ref[m], nbx0, nby0);
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
                evaluateModel(s, l, as, tt, t, n, true, pert.data());
                *astPar[j] = keep;
                // A null derivative (epoch at t0 or far from the peak) is legitimate and contributes nothing.
                for (int m = 0; m < n; ++m) {
                    double nbx, nby;
                    nbrShare(pert[m], nbx, nby);
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
                co.covTetPiV[v][q] = 0.0;
                if (v == AV_W or co.nepochA[q] < Ny) continue;      // W is filled in ErrorCal
                for (int j = 0; j < Ny; ++j)
                    for (int m = 0; m < Ny; ++m) gsl_matrix_set(Fin, j, m, co.FBV[v][q][j * Ny + m]);
                co.okBV[v][q] = invertNormalized(Fin, Fout, kAllAst, co.condBV[v][q], nullptr);
                if (co.okBV[v][q]) {
                    for (int k = 0; k < Ny; ++k)
                        co.ErbV[v][q][k] = std::sqrt(std::fabs(gsl_matrix_get(Fout, k, k)));
                    co.covTetPiV[v][q] = gsl_matrix_get(Fout, 0, 3);   //cov(thetaE, piE), for the mass
                }
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
// The fractional 1-sigma of one quantity measured two independent ways, from the two: the inverse-variance
// combination 1 / sqrt(1/a^2 + 1/b^2). -1 = not measured; with one route absent the other is returned
// untouched, and with both absent -1.
static double combineRel(double a, double b)
{
  if (a >= 0.0 and b >= 0.0) {
      const double h = std::hypot(a, b);
      return (h > 0.0) ? a * b / h : 0.0;
  }
  return std::max(a, b);
}

// The lens mass Ml = thetaE / (kappa piE) and the relative parallax pi_rel = thetaE piE are products of the same
// two measured quantities: ln Ml = ln thetaE - ln piE + const, ln pi_rel = ln thetaE + ln piE, so with rt, rpi
// the fractional errors of the thetaE and the piE used and c = cov(ln thetaE, ln piE),
//     var(ln Ml) = rt^2 + rpi^2 - 2c,     var(ln pi_rel) = rt^2 + rpi^2 + 2c.
// Both enter the output, with the same c of opposite sign: a positive c shrinks the mass error and grows the distance
// error, a negative one the reverse. On the astrometric piE route c is mostly negative (sightline 795: 74% of the
// joint events on that route), so there the mass error is larger and the distance error smaller than in quadrature.
// The thetaE used is the inverse-variance combination of the astrometric and the photometric (theta* / rho) one,
// ln thetaE_c = (1 - wp) ln thetaE_astr + wp ln thetaE_phot, wp = rt^2 / rph^2 (rph = relTetEph, the photometric
// route's own fractional error; wp = 0 where ln rho is not a parameter of the partition). The piE used is the
// photometric or the astrometric one, whichever is the better. c then has one source per route:
//   - photometric piE: it shares a matrix with ln rho, ln thetaE_phot = ln theta* - ln rho (theta* independent), so
//     c = -wp cov(ln rho, ln piE) = -wp corr s rpi (corr = corrRhoPiE, s = sigma(ln rho)). The astrometric thetaE
//     is in another matrix: no covariance with it;
//   - astrometric piE: it shares the astrometric matrix with the astrometric thetaE, so c = (1 - wp)
//     cov(thetaE, piE) / (thetaE piE) from that matrix's inverse (covB; piE > 0 always); the photometric thetaE
//     is in the other matrix: no covariance with this piE.
// |c| <= rt rpi (Cauchy-Schwarz) so both variances are >= 0; a negative one within rounding is set to 0, any other
// is a bug (CHECK). With c = 0 (no ln rho on the photometric route) the expression is sqrt(rpi^2 + rt^2), as it
// always was. This function is the single source of c: relMl, relMlV, resu[9] and resu[10] all take it.
static double lnCovTetPi(bool photPi, double wp, double rpi, double corr, double s, double covB, double tetE, double piE)
{
  if (photPi) return std::isfinite(corr) ? -(wp * corr * s * rpi) : 0.0;     //corr is NaN only for rpi = 0
  return (1.0 - wp) * covB / ((std::fabs(tetE) + eps) * (std::fabs(piE) + eps));
}

// The variance of ln Ml (c2 = -2c) or of ln pi_rel (c2 = +2c), guarded as above.
static double lnVar(double rt, double rpi, double c2)
{
  const double base = rt * rt + rpi * rpi;
  const double var  = base + c2;
  if (var < 0.0) {
      CHECK(var >= -1.0e-12 * base);   //rounding only
      return 0.0;
  }
  return var;
}

// Fractional 1-sigma of the lens mass Ml = thetaE / (kappa piE): sqrt(rt^2 + rpi^2 - 2c), c from lnCovTetPi.
static double relMass(double rt, double rpi, double c)
{
  return std::sqrt(lnVar(rt, rpi, -2.0 * c));
}

// ---------------------------------------------------------------------------------------------------------
// The lens's motion on the sky and its error (lensMotion, events/vbm_model.h). The six parameters of the model,
// p = (tE, xi, piE, thetaE, mus1, mus2), are not fitted together: the photometric matrix gives tE, xi, piE (and ln
// rho), the astrometric one thetaE, mus1, mus2, piE, and each holds the other's parameters fixed, so the two
// matrices are independent blocks. p is a linear map of independent variables,
//     a = (tE, xi, piE_phot, ln rho)  from inverA[q],   b = (thetaE_ast, mus1, mus2, piE_ast)  from inverB[q],
//     e = the error of ln theta*  (variance THETASTAR_SIGLN^2),
// and its covariance is L Sigma L^T. The pieces of L:
//   - piE is the photometric one (a) or the astrometric one (b), the route relMl takes (the better one);
//   - thetaE is the inverse-variance combination of the astrometric thetaE and the photometric theta* / rho
//     that relTetEc is: delta thetaE = (1 - wp) delta thetaE_ast + wp thetaE (-delta ln rho + e), with the photometric
//     weight wp = relTetEc^2 / relTetEph^2 (the weight of ln thetaE_phot in ln thetaE_c, as lnCovTetPi has it). So
//     thetaE carries the covariance of thetaE_ast with mus (and with piE_ast) from the astrometric matrix, and
//     that of ln rho with tE, xi and piE_phot from the photometric one;
//   - mus1, mus2 are the astrometric ones.
// Sigma(thetaE)/thetaE^2 comes out as relTetEc^2 and, on either route, the variance of ln thetaE - ln piE as relMl^2
// and that of ln thetaE + ln piE as (rel_Dl / (1 - D_L/D_S))^2 (relMass, lnCovTetPi), which the fixture asserts.
bool lensMotionCovariance(const covarian& co, const lens& l, int q, double S[NLM][NLM])
{
  for (int i = 0; i < NLM; ++i) for (int j = 0; j < NLM; ++j) S[i][j] = 0.0;
  if (!(co.okA[q] and co.okB[q] and co.relMl[q] >= 0.0)) return false;
  const auto& Era = co.Era[q];
  const auto& Erb = co.Erb[q];
  const double rp1 = (Era[3] >= 0.0) ? Era[3] / (std::fabs(l.piE) + eps) : -1.0;
  const double rp2 = (Erb[3] >= 0.0) ? Erb[3] / (std::fabs(l.piE) + eps) : -1.0;
  const bool photPi = (rp1 >= 0.0 and (rp2 < 0.0 or rp1 < rp2));      //the route relMl takes
  const bool rho    = (Era[9] >= 0.0);                                 //ln rho is a parameter of this partition
  if (Era[1] < 0.0 or Era[4] < 0.0 or Erb[0] < 0.0 or Erb[1] < 0.0 or Erb[2] < 0.0 or (rp1 < 0.0 and rp2 < 0.0))
      return false;
  const double wp = (rho and co.relTetEph[q] > 0.0 and co.relTetEc[q] >= 0.0)
                  ? co.relTetEc[q] * co.relTetEc[q] / (co.relTetEph[q] * co.relTetEph[q]) : 0.0;

  // Sigma of (a[0..3], b[4..7], e[8]).
  static const int ia[4] = {1, 4, 3, 9};
  constexpr int NB = 9;
  double Sb[NB][NB] = {};
  const int na = rho ? 4 : 3;
  for (int i = 0; i < na; ++i)
      for (int j = 0; j < na; ++j) Sb[i][j] = gsl_matrix_get(co.inverA[q].get(), ia[i], ia[j]);
  for (int i = 0; i < Ny; ++i)
      for (int j = 0; j < Ny; ++j) Sb[4 + i][4 + j] = gsl_matrix_get(co.inverB[q].get(), i, j);
  Sb[8][8] = THETASTAR_SIGLN * THETASTAR_SIGLN;

  double L[NLM][NB] = {};
  L[LM_TE][0] = 1.0;
  L[LM_XI][1] = 1.0;
  L[LM_PIE][photPi ? 2 : 7] = 1.0;
  L[LM_TETE][4] = 1.0 - wp;
  if (rho) { L[LM_TETE][3] = -wp * l.tetE; L[LM_TETE][8] = wp * l.tetE; }
  L[LM_MUS1][5] = 1.0;
  L[LM_MUS2][6] = 1.0;

  double LS[NLM][NB] = {};
  for (int i = 0; i < NLM; ++i)
      for (int k = 0; k < NB; ++k)
          for (int m = 0; m < NB; ++m) LS[i][k] += L[i][m] * Sb[m][k];
  for (int i = 0; i < NLM; ++i)
      for (int j = 0; j < NLM; ++j)
          for (int k = 0; k < NB; ++k) S[i][j] += LS[i][k] * L[j][k];
  return true;
}

void lensMotionPropagate(const MotionContext& c, const double p[NLM], const double S[NLM][NLM],
                         double Smu[2][2], double Sv[2][2])
{
  // Central differences of the closed-form function: smooth, so a relative step of 1e-5 is good to ~1e-10 (the
  // function is linear in mus). xi and the proper motions have a floor to the step, a zero or tiny value must not
  // make it vanish.
  static const double floorStep[NLM] = {0.0, 1.0, 0.0, 0.0, 1.0e-3, 1.0e-3};
  double Jmu[2][NLM], Jv[2][NLM];
  for (int k = 0; k < NLM; ++k) {
      const double h = 1.0e-5 * std::max(std::fabs(p[k]), floorStep[k]);
      double pp[NLM], pm[NLM];
      for (int i = 0; i < NLM; ++i) pp[i] = pm[i] = p[i];
      pp[k] += h; pm[k] -= h;
      const LensMotion a = lensMotion(c, pp[LM_TE], pp[LM_XI], pp[LM_PIE], pp[LM_TETE], pp[LM_MUS1], pp[LM_MUS2]);
      const LensMotion b = lensMotion(c, pm[LM_TE], pm[LM_XI], pm[LM_PIE], pm[LM_TETE], pm[LM_MUS1], pm[LM_MUS2]);
      for (int i = 0; i < 2; ++i) {
          Jmu[i][k] = (a.mu[i] - b.mu[i]) / (2.0 * h);
          Jv[i][k]  = (a.v[i]  - b.v[i])  / (2.0 * h);
      }
  }
  for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j) {
          Smu[i][j] = Sv[i][j] = 0.0;
          for (int k = 0; k < NLM; ++k)
              for (int m = 0; m < NLM; ++m) {
                  Smu[i][j] += Jmu[i][k] * S[k][m] * Jmu[j][m];
                  Sv[i][j]  += Jv[i][k]  * S[k][m] * Jv[j][m];
              }
      }
}

void ErrorCal(covarian & co, lens & l , source & s){

  // Per-survey 1-sigma forecasts. -1.0 marks a partition that is not characterizable (too few
  // epochs or singular matrix) and, for the photometric sigmas, a parameter outside the
  // partition's active subset. Downstream code must test okA/okB and skip negative sigmas.
  for (int q = 0; q < NSURV; ++q) {
      for (int k = 0; k < Nx; ++k) co.Era[q][k] = -1.0;
      if (co.okA[q]) {
          for (int k : activePhotParams(co, q)) {
              co.Era[q][k] = std::sqrt(std::fabs(gsl_matrix_get(co.inverA[q].get(), k, k)));
          }
      }
      for (int k = 0; k < Ny; ++k) {
          co.Erb[q][k] = co.okB[q] ? std::sqrt(std::fabs(gsl_matrix_get(co.inverB[q].get(), k, k)))
                                   : -1.0;
      }
  }

  // The Einstein radius has two routes. Astrometry gives it directly (Erb[0]). When the finite source is seen
  // (ln rho a parameter of the partition, index 9) the light curve gives it as thetaE = theta* / rho, with
  // the fractional error of ln rho and of the source's angular radius theta* (THETASTAR_SIGLN) added in
  // quadrature. The two are combined by inverse variance, as independent measurements (they come from different
  // matrices, and theta* is independent of the fit). The lens mass below takes the combined thetaE; its error
  // carries the covariance of ln rho with the photometric piE, which share a matrix (lnCovTetPi), co.corrRhoPiE
  // being its correlation coefficient.
  for (int q = 0; q < NSURV; ++q) {
      const double rastr = (co.Erb[q][0] >= 0.0) ? co.Erb[q][0] / (std::fabs(l.tetE) + eps) : -1.0;
      co.relTetEph[q] = (co.Era[q][9] >= 0.0) ? std::sqrt(co.Era[q][9] * co.Era[q][9] + THETASTAR_SIGLN * THETASTAR_SIGLN) : -1.0;
      co.relTetEc[q]  = combineRel(rastr, co.relTetEph[q]);
      co.corrRhoPiE[q] = std::numeric_limits<double>::quiet_NaN();
      if (co.Era[q][9] >= 0.0 and co.Era[q][3] > 0.0)
          co.corrRhoPiE[q] = gsl_matrix_get(co.inverA[q].get(), 9, 3) / (co.Era[q][9] * co.Era[q][3]);
  }

  // Fractional 1-sigma on the lens mass, per survey: Ml = tetE / (kappa * piE), kappa = 8.144 mas/Msun, so the
  // fractional errors add in quadrature minus twice their covariance (relMass, lnCovTetPi: ln rho with the
  // photometric piE, or thetaE with the astrometric piE, by the route). piE has a photometric (Era[3]) and an
  // astrometric (Erb[3]) route; take the better available one, and report no mass if either
  // ingredient is missing (-1.0 means "not measured", never a measurement).
  for (int q = 0; q < NSURV; ++q) {
      const double rp1 = (co.Era[q][3] >= 0.0) ? co.Era[q][3] / (std::fabs(l.piE)  + eps) : -1.0;
      const double rp2 = (co.Erb[q][3] >= 0.0) ? co.Erb[q][3] / (std::fabs(l.piE)  + eps) : -1.0;
      const double rp  = (rp1 >= 0.0 and rp2 >= 0.0) ? MIN(rp1, rp2) : std::max(rp1, rp2);
      const double rt  = co.relTetEc[q];
      const bool   photPi = (rp1 >= 0.0 and (rp2 < 0.0 or rp1 < rp2));   //the route MIN picked
      if (!(rp >= 0.0 and rt >= 0.0)) { co.relMl[q] = -1.0; continue; }
      const double wp   = (co.relTetEph[q] >= 0.0) ? rt * rt / (co.relTetEph[q] * co.relTetEph[q]) : 0.0;
      const double covB = co.okB[q] ? gsl_matrix_get(co.inverB[q].get(), 0, 3) : 0.0;
      co.relMl[q] = relMass(rt, rp, lnCovTetPi(photPi, wp, rp, co.corrRhoPiE[q], co.Era[q][9], covB, l.tetE, l.piE));
  }

  // Astrometric noise variants: W is the main result; N and P were inverted in FisherM. Lens-mass
  // errors follow the same rule, with each variant's own tetE, astrometric piE and their covariance.
  for (int q = 0; q < NSURV; ++q) {
      co.okBV[AV_W][q] = co.okB[q];
      co.condBV[AV_W][q] = co.condB[q];
      for (int k = 0; k < Ny; ++k) co.ErbV[AV_W][q][k] = co.Erb[q][k];
      co.covTetPiV[AV_W][q] = co.okB[q] ? gsl_matrix_get(co.inverB[q].get(), 0, 3) : 0.0;
      for (int v = 0; v < NAVAR; ++v) {
          const auto& E = co.ErbV[v][q];
          const double rp1 = (co.Era[q][3] >= 0.0) ? co.Era[q][3] / (std::fabs(l.piE) + eps) : -1.0;
          const double rp2 = (E[3] >= 0.0) ? E[3] / (std::fabs(l.piE) + eps) : -1.0;
          const double rp  = (rp1 >= 0.0 and rp2 >= 0.0) ? MIN(rp1, rp2) : std::max(rp1, rp2);
          const double rt  = combineRel((E[0] >= 0.0) ? E[0] / (std::fabs(l.tetE) + eps) : -1.0, co.relTetEph[q]);
          const bool   photPi = (rp1 >= 0.0 and (rp2 < 0.0 or rp1 < rp2));
          if (!(rp >= 0.0 and rt >= 0.0)) { co.relMlV[v][q] = -1.0; continue; }
          const double wp = (co.relTetEph[q] >= 0.0) ? rt * rt / (co.relTetEph[q] * co.relTetEph[q]) : 0.0;
          co.relMlV[v][q] = relMass(rt, rp, lnCovTetPi(photPi, wp, rp, co.corrRhoPiE[q], co.Era[q][9],
                                                       co.covTetPiV[v][q], l.tetE, l.piE));
      }
  }

  // Joint event summary: resu[] keeps its index layout for the output columns. Per-survey
  // numbers live in Era[]/Erb[].
  const auto& Era = co.Era[SJOINT];
  const auto& Erb = co.Erb[SJOINT];

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
  const bool photPiJ = (co.resu[3] >= 0.0 and (co.resu[8] < 0.0 or co.resu[3] < co.resu[8]));   //the route MIN picks
  if      (co.resu[3] < 0.0 and co.resu[8] < 0.0) co.resu[3] = -1.0;
  else if (co.resu[3] < 0.0)                      co.resu[3] = co.resu[8];
  else if (co.resu[8] >= 0.0)                     co.resu[3] = MIN(co.resu[3], co.resu[8]);

  // Lens mass: Ml = tetE/(kappa*piE), with the combined thetaE (resu[5] stays the astrometric one: the column and
  // the proper-motion errors below read it) and the piE route above; var(ln Ml) = rt^2 + rpi^2 - 2c (relMass).
  // Lens distance: Dl = 1 / (pi_rel + pi_S) with pi_S = 1/Ds known and pi_rel = tetE piE, so d ln Dl = -(1 - Dl/Ds)
  // d ln pi_rel and rel_Dl = (1 - Dl/Ds) sqrt(var(ln pi_rel)), var(ln pi_rel) = rt^2 + rpi^2 + 2c, the same c
  // of opposite sign in the two. -1 where the mass is.
  {
      const double rt = co.relTetEc[SJOINT], rpi = co.resu[3];
      co.resu[9] = co.resu[10] = -1.0;
      if (rpi >= 0.0 and rt >= 0.0) {
          const double wp   = (co.relTetEph[SJOINT] >= 0.0) ? rt * rt / (co.relTetEph[SJOINT] * co.relTetEph[SJOINT]) : 0.0;
          const double covB = co.okB[SJOINT] ? gsl_matrix_get(co.inverB[SJOINT].get(), 0, 3) : 0.0;
          const double c    = lnCovTetPi(photPiJ, wp, rpi, co.corrRhoPiE[SJOINT], Era[9], covB, l.tetE, l.piE);
          co.resu[9]  = relMass(rt, rpi, c);
          co.resu[10] = std::fabs(std::sqrt(lnVar(rt, rpi, 2.0 * c)) * (s.Ds - l.Dl) / s.Ds);
      }
  }

  // The lens's proper motion and velocity on the sky, per partition (lensMotion; the covariance of its parameters
  // is lensMotionCovariance). Only the default astrometric matrix (variant W) enters, not the noise variants N, P.
  // mu_rel = thetaE / tE is the geocentric one (tE, xi are geocentric); the lens's heliocentric motion adds the
  // source's (mus, astrometry) and the Earth's velocity times pi_rel = thetaE piE. First the identity that pins the
  // model: at the true parameters it must return the draw's mul1, mul2 and v_L.
  const MotionContext mc = motionContext(l.t0, 1.0 / s.Ds);
  const double pTrue[NLM] = {l.tE, s.xi, l.piE, l.tetE, s.mus1, s.mus2};
  const LensMotion m0 = lensMotion(mc, pTrue[LM_TE], pTrue[LM_XI], pTrue[LM_PIE], pTrue[LM_TETE], pTrue[LM_MUS1], pTrue[LM_MUS2]);
  {
      const double vconv = l.Dl * AU / (1000.0 * 24.0 * 3600.0);   //[km/s per mas/day]
      const double dmu = std::hypot(m0.mu[0] - l.mul1, m0.mu[1] - l.mul2) / (std::hypot(l.mul1, l.mul2) + eps);
      const double dv  = std::hypot(m0.v[0] - l.mul1 * vconv, m0.v[1] - l.mul2 * vconv) / (std::hypot(l.mul1, l.mul2) * vconv + eps);
      co.lmCheck = std::max(dmu, dv);
  }
  const double muAbs = std::hypot(m0.mu[0], m0.mu[1]);
  for (int q = 0; q < NSURV; ++q) {
      co.relMuL[q] = -1.0;
      co.sigVL[q]  = -1.0;
      co.covMuL[q].fill(0.0);
      co.covVL[q].fill(0.0);
      double Sp[NLM][NLM], Smu[2][2], Sv[2][2];
      if (!lensMotionCovariance(co, l, q, Sp) or !(muAbs > 0.0)) continue;
      lensMotionPropagate(mc, pTrue, Sp, Smu, Sv);
      // sigma(|x|) of a vector x with covariance C: the error along x, x_hat^T C x_hat.
      const double ux = m0.mu[0] / muAbs, uy = m0.mu[1] / muAbs;
      co.relMuL[q] = std::sqrt(std::fabs(ux * ux * Smu[0][0] + 2.0 * ux * uy * Smu[0][1] + uy * uy * Smu[1][1])) / muAbs;
      // sigVL: the largest eigenvalue of the symmetric 2x2 Sigma_v, (a+d)/2 + sqrt(((a-d)/2)^2 + b^2), whatever
      // the frame the velocity is later referred to.
      co.sigVL[q]  = std::sqrt(0.5 * (Sv[0][0] + Sv[1][1]) + std::hypot(0.5 * (Sv[0][0] - Sv[1][1]), Sv[0][1]));
      co.covMuL[q] = {Smu[0][0], Smu[0][1], Smu[1][0], Smu[1][1]};
      co.covVL[q]  = {Sv[0][0], Sv[0][1], Sv[1][0], Sv[1][1]};
  }
  const auto& cm = co.covMuL[SJOINT];
  const bool lmJ = (co.relMuL[SJOINT] >= 0.0);
  co.resu[11] = lmJ ? double(std::sqrt(cm[0]) / (std::fabs(l.mul1) + eps)) : -1.0; //mu_lens_1
  co.resu[12] = lmJ ? double(std::sqrt(cm[3]) / (std::fabs(l.mul2) + eps)) : -1.0; //mu_lens_2
  co.resu[13] = co.relMuL[SJOINT];                                                 //mu_lens
  co.resu[14] = std::sqrt(Erb[1]  * Erb[1]  +  Erb[2] * Erb[2])  / (std::fabs(s.mus) + eps); // mu_source
}
