// FisherM (information matrices by finite differences) and ErrorCal (sigmas from the inverses).
#include "fisher/fisher.h"
#include "events/lightcurve.h"
#include "events/vbm_model.h"
#include "surveys/noise.h"
#include <cstring>
#include <limits>

namespace {
long gRhoFallback = 0, gRhoFallbackLost = 0;
}
long rhoFallbackCount()     { return gRhoFallback; }
long rhoFallbackLostCount() { return gRhoFallbackLost; }

namespace {
// The slots of an epoch's compact derivative row: the six parameters every epoch is sensitive to (kShared, in this
// order), then the epoch's own flux pair. The other flux parameters have derivative exactly 0 for it, so they are never
// stored and the accumulation touches only these eight entries (an Nx-wide row, Nx = 20, would be mostly zeros).
constexpr int kShared[6] = {0, 1, 3, 4, 5, 9};   //u0, tE, piE, xi, t0, ln rho
constexpr int kSlotRho = 5, kSlotFb = 6, kSlotMbs = 7;
using EpochDeriv = std::array<double, 8>;

// The derivatives of one flux pair (fb, mbs) over the epochs `pos` (positions in the telescope's arrays), from the
// reference magnification. m = mbs - 2.5 log10(fb A + 1 - fb) depends on fb and mbs through the flux alone and A is the
// same in every band, so no new model is evaluated. The fb step is binned on THIS pair's blend fraction: bins and steps
// are chosen so fb + bb never leaves [0,1], which holds only if bb is computed from the fb it perturbs. Roman's F146
// blend fraction is routinely higher than Rubin's, and each Rubin band has its own. One function serves every pair, so
// r and F146 are computed exactly as they were before the other bands existed.
void fluxPairDerivs(double fb, double mbs, double fbScale, double mbsStep, const std::vector<int>& pos,
                    const ModelPoint* ref, const double* lgf, const int* ix, std::vector<EpochDeriv>& dm)
{
    std::array<double, 2> bb;
    if (fb < FB_BIN_LO)      bb = {+FB_STEP_SMALL, +FB_STEP_LARGE};
    else if (fb < FB_BIN_HI) bb = {-FB_STEP_SMALL, +FB_STEP_SMALL};
    else                     bb = {-FB_STEP_SMALL, -FB_STEP_LARGE};

    // Apply the plateau scaling and clamp so the step cannot leave the physical range [0,1]
    // (the bin edges guarantee this only at unscaled sizes). Note the outer fb bins are two
    // same-signed (one-sided) steps, so fb carries an O(h) bias there; kFDStepScale makes it
    // negligible.
    for (int b = 0; b < 2; ++b) {
        double step = bb[b] * fbScale;
        const double lo = 1.0e-6 - fb;        //keeps fb strictly above 0
        const double hi = 1.0 - 1.0e-6 - fb;  //keeps fb strictly below 1
        if (step < lo) step = lo;
        if (step > hi) step = hi;
        bb[b] = step;
    }

    for (int h = 0; h < 2; ++h) {
        const double dfb  = bb[h];
        const double dmbs = mbsStep * sig[h];
        const double fbp  = fb + dfb, mbsp = mbs + dmbs;  //the perturbed values
        CHECK(dfb != 0.0);
        CHECK(fbp > 0.0);
        CHECK(fbp <= 1.0);
        CHECK(dmbs != 0.0);
        for (const int m : pos) {
            const double mag0  = mbs - lgf[m];
            const double magfb = mbs - 2.5 * std::log10(ref[m].A * fbp + 1.0 - fbp);
            const double magmb = mbsp - lgf[m];
            const double derfb = (magfb - mag0) / dfb;
            const double dermb = (magmb - mag0) / dmbs;
            double& dfbm = dm[ix[m]][kSlotFb];
            double& dmbm = dm[ix[m]][kSlotMbs];
            dfbm = (h == 0) ? derfb : (dfbm + derfb) * 0.5;
            dmbm = (h == 0) ? dermb : (dmbm + dermb) * 0.5;
        }
    }
}
} // namespace


// Rubin's information from several bands: per coordinate, each band's own neighbours' offset is marginalised from that
// band's sums (Schur complement of slot Ny); what is left of the bands is added and the ONE shared reference position
// is taken out of the total. A band with no offset (f_N = 0, Swhh = 0, or AST_BLEND_OFFSET_FREE off) goes in as it is.
// When every band sits at baseline each band's offset has already absorbed its own mean and c vanishes: the reference
// is degenerate and nothing more is subtracted (c against the summed weight, as physical() compares against raw).
void fisher_detail::rubinBandInformation(const std::array<AstSums, kRubinBands>& bands, std::array<double, Ny * Ny>& FL)
{
    constexpr int D = kAstD;
    for (int cc = 0; cc < 2; ++cc) {
        std::array<double, Ny * Ny> P{};
        std::array<double, Ny>      r{};
        double c = 0.0, SwTot = 0.0;
        for (int b = 0; b < kRubinBands; ++b) {
            const AstSums& S = bands[b];
            if (not (S.Sw > 0.0)) continue;
            const double Swh = S.Swd[cc][Ny], Swhh = S.Swdd[cc][Ny * D + Ny];
            SwTot += S.Sw;
            const bool off = AST_BLEND_OFFSET_FREE and Swhh > 0.0;
            for (int j = 0; j < Ny; ++j) {
                const double Swhd_j = S.Swdd[cc][j * D + Ny];
                for (int m = 0; m < Ny; ++m)
                    P[j * Ny + m] += S.Swdd[cc][j * D + m] - (off ? Swhd_j * S.Swdd[cc][m * D + Ny] / Swhh : 0.0);
                r[j] += S.Swd[cc][j] - (off ? Swhd_j * Swh / Swhh : 0.0);
            }
            c += S.Sw - (off ? Swh * Swh / Swhh : 0.0);
        }
        const bool ref = c > 1.0e-9 * SwTot;
        for (int j = 0; j < Ny; ++j)
            for (int m = 0; m < Ny; ++m) FL[j * Ny + m] += P[j * Ny + m] - (ref ? r[j] * r[m] / c : 0.0);
    }
}

void fisher_detail::astPhysical(const AstPerCoord& G, const std::array<double, 2>& raw, std::array<double, Ny * Ny>& out)
{
    constexpr int D = kAstD;
    for (int cc = 0; cc < 2; ++cc) {
        const double hh = G[cc][Ny * D + Ny];
        const bool   fit = AST_BLEND_OFFSET_FREE and raw[cc] > 0.0 and hh > 1.0e-9 * raw[cc];
        for (int j = 0; j < Ny; ++j)
            for (int m = 0; m < Ny; ++m)
                out[j * Ny + m] += G[cc][j * D + m] - (fit ? G[cc][j * D + Ny] * G[cc][m * D + Ny] / hh : 0.0);
    }
}

void fisher_detail::romanInformation(const std::map<long, AstSums>& days, bool oneFrame, double sigc, double flInfo,
                                     std::array<double, Ny * Ny>& out)
{
    constexpr int D = kAstD;
    const bool lensFree = flInfo >= 0.0;
    std::array<double, 2> raw{};
    for (const auto& [day, k] : days)
        for (int cc = 0; cc < 2; ++cc) raw[cc] += k.Swdd[cc][Ny * D + Ny];
    std::vector<AstGroup> groups(oneFrame ? 1 : 2);        // W: one frame; N, P: one per roll
    for (auto& g : groups) g.nact = lensFree ? D : D - 1;
    for (const auto& [day, k] : days)
        groups[oneFrame ? 0 : k.roll].fold(k, sigc);
    AstPerCoord G{};
    for (const auto& g : groups) g.marginaliseInto(G);
    if (not lensFree) { astPhysical(G, raw, out);  return; }

    // The offset out as in astPhysical, keeping (physical, f_L); the two coordinates added, the prior on f_L added, and f_L
    // out (a scalar Schur complement).
    constexpr int E = Ny + 1;
    std::array<double, E * E> T{};
    int idx[E];
    for (int j = 0; j < Ny; ++j) idx[j] = j;
    idx[Ny] = Ny + 1;                                          //the sums' slot of f_L
    for (int cc = 0; cc < 2; ++cc) {
        const double hh = G[cc][Ny * D + Ny];
        const bool   fit = AST_BLEND_OFFSET_FREE and raw[cc] > 0.0 and hh > 1.0e-9 * raw[cc];
        for (int a = 0; a < E; ++a)
            for (int m = 0; m < E; ++m)
                T[a * E + m] += G[cc][idx[a] * D + idx[m]] - (fit ? G[cc][idx[a] * D + Ny] * G[cc][idx[m] * D + Ny] / hh : 0.0);
    }
    T[Ny * E + Ny] += flInfo;
    const double tff = T[Ny * E + Ny];
    const bool good = tff > 0.0 and std::isfinite(tff);
    for (int j = 0; j < Ny; ++j)
        for (int m = 0; m < Ny; ++m)
            out[j * Ny + m] += good ? T[j * E + m] - T[j * E + Ny] * T[Ny * E + m] / tff : std::numeric_limits<double>::quiet_NaN();
}

// The same with the lens's light fraction of every band a parameter beside the four physical ones (slot Ny + 1 of the
// band's sums; each band has sums of its own, so that slot is always that band's fraction). Per coordinate: the band's
// offset is taken out of the band's (physical, fraction) block exactly as above, the blocks are embedded in the kept space
// (physical, then one fraction per band that has one), added over the bands, and the shared reference position taken out
// of the total. The coordinates are added (the fraction is the same number in x and y), the prior precision flInfo[b] put
// on the fractions' diagonal, and the fractions marginalised (Cholesky of their block, which the prior makes positive
// definite). The l column of a band without a fraction is dropped, its physical part kept.
bool fisher_detail::rubinBandInformationLensFree(const std::array<AstSums, kRubinBands>& bands,
                                                 const std::array<double, kRubinBands>& flInfo,
                                                 std::array<double, Ny * Ny>& FL)
{
    constexpr int D = kAstD, SL = Ny + 1;                       //SL: the sums' slot of l
    constexpr int KMAX = Ny + kRubinBands;
    std::array<int, kRubinBands> slot{};                          //band's index in the kept space, -1: none
    int nf = 0;
    for (int b = 0; b < kRubinBands; ++b) slot[b] = (bands[b].Sw > 0.0 and flInfo[b] >= 0.0) ? Ny + nf++ : -1;
    const int K = Ny + nf;
    std::array<double, KMAX * KMAX> T{};
    for (int cc = 0; cc < 2; ++cc) {
        std::array<double, KMAX * KMAX> P{};
        std::array<double, KMAX>        r{};
        double c = 0.0, SwTot = 0.0;
        for (int b = 0; b < kRubinBands; ++b) {
            const AstSums& S = bands[b];
            if (not (S.Sw > 0.0)) continue;
            const double Swh = S.Swd[cc][Ny], Swhh = S.Swdd[cc][Ny * D + Ny];
            SwTot += S.Sw;
            const bool off = AST_BLEND_OFFSET_FREE and Swhh > 0.0;
            int from[Ny + 1], to[Ny + 1];                         //the band's components: slot in the sums, in the kept space
            for (int j = 0; j < Ny; ++j) { from[j] = j;  to[j] = j; }
            from[Ny] = SL;  to[Ny] = slot[b];
            for (int a = 0; a <= Ny; ++a) {
                if (to[a] < 0) continue;
                const double Swha = S.Swdd[cc][from[a] * D + Ny];
                for (int m = 0; m <= Ny; ++m) {
                    if (to[m] < 0) continue;
                    P[to[a] * KMAX + to[m]] += S.Swdd[cc][from[a] * D + from[m]] - (off ? Swha * S.Swdd[cc][from[m] * D + Ny] / Swhh : 0.0);
                }
                r[to[a]] += S.Swd[cc][from[a]] - (off ? Swha * Swh / Swhh : 0.0);
            }
            c += S.Sw - (off ? Swh * Swh / Swhh : 0.0);
        }
        const bool ref = c > 1.0e-9 * SwTot;
        for (int a = 0; a < K; ++a)
            for (int m = 0; m < K; ++m) T[a * KMAX + m] += P[a * KMAX + m] - (ref ? r[a] * r[m] / c : 0.0);
    }
    for (int b = 0; b < kRubinBands; ++b) if (slot[b] >= 0) T[slot[b] * KMAX + slot[b]] += flInfo[b];

    // Out with the fractions: Cholesky T_ff = Lc Lc^T, then X = T_ff^-1 T_fp, FL += T_pp - T_pf X.
    std::array<double, KMAX * KMAX> Lc{};
    for (int a = 0; a < nf; ++a)
        for (int m = 0; m <= a; ++m) {
            double v = T[(Ny + a) * KMAX + Ny + m];
            for (int p = 0; p < m; ++p) v -= Lc[a * KMAX + p] * Lc[m * KMAX + p];
            if (a == m) {
                if (not (v > 0.0) or not std::isfinite(v)) return false;
                Lc[a * KMAX + a] = std::sqrt(v);
            } else {
                Lc[a * KMAX + m] = v / Lc[m * KMAX + m];
            }
        }
    std::array<double, KMAX * Ny> X{};                             //row f, column m of T_ff^-1 T_fp
    for (int m = 0; m < Ny; ++m) {
        for (int a = 0; a < nf; ++a) {                              //forward: Lc y = T_fp[:, m]
            double v = T[(Ny + a) * KMAX + m];
            for (int p = 0; p < a; ++p) v -= Lc[a * KMAX + p] * X[p * Ny + m];
            X[a * Ny + m] = v / Lc[a * KMAX + a];
        }
        for (int a = nf - 1; a >= 0; --a) {                         //back: Lc^T x = y
            double v = X[a * Ny + m];
            for (int p = a + 1; p < nf; ++p) v -= Lc[p * KMAX + a] * X[p * Ny + m];
            X[a * Ny + m] = v / Lc[a * KMAX + a];
        }
    }
    for (int j = 0; j < Ny; ++j)
        for (int m = 0; m < Ny; ++m) {
            double v = T[j * KMAX + m];
            for (int a = 0; a < nf; ++a) v -= T[j * KMAX + Ny + a] * X[a * Ny + m];
            FL[j * Ny + m] += v;
        }
    return true;
}

void FisherM(source & s, lens & l, astromet & as,  covarian & co, int ndw)
{
    co.flagi =+ 1;
    // r's pair is kept twice, per telescope (s.fb[0], s.mbs[0]: the astrometric part, resu[2]) and per filter (s.blend[2],
    // s.magb[2]: the photometric fit). They are the same numbers because RUBIN_REF_BANDS = {2}; a change that parts them
    // must change both readers.
    CHECK(s.fb[0] == s.blend[2] and s.mbs[0] == s.magb[2]);
    CHECK(s.fLens[0] == s.fLensF[2] and s.fLens[1] == s.fLensF[6]);   //the lens's share, likewise per telescope and per filter

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
    co.Delta1[8] = FD_STEP_MBS; // mbs1 [mag] -- Roman baseline magnitude
    for (int b = 0; b < kRubinBands; ++b)
        co.Delta1[mbsIndexOfFilter(b)] = FD_STEP_MBS; // mbs of Rubin's u g r i z y (11, 13, 6 = mbs0, 15, 17, 19) [mag]
    co.Delta1[9] = FD_STEP_LNRHO; // ln rho -- the source radius in Einstein radii, perturbed as rho * exp(Delta)

    // The legacy step constants are ~25% of the parameter value, far up the truncation-error
    // branch of the step-size sweep (./fishertest --sweep, tests/c3_step_sweep.py). kFDStepScale
    // (see config/parameters.h) scales them into the plateau, where sigma is flat to <0.2%; below
    // ~1e-6 of the legacy value round-off dominates. Oversized steps manufacture information in
    // near-degenerate directions, so the correct sigmas for short events are much larger.
    for (int q = 0; q < Nx; ++q) co.Delta1[q] *= kFDStepScale * co.deltaScale[q];
    // mbs enters with unit slope, so its finite difference is exact for any step. The blend fractions
    // (fb0 at 2, fb1 at 7, the Rubin bands' at 10, 12, 14, 16, 18) use the steps bb[] set per flux pair below.

    // Three matrices, zeroed together: joint, Rubin-only, Roman-only (see SurveyIdx).
    for (int q = 0; q < NSURV; ++q) {
        co.nepochA[q] = 0;
        co.okA[q] = 0;
        co.okB[q] = 0;
        // Reset here: a partition rejected for too few epochs never reaches invert_matrix.
        co.condA[q] = -1.0;
        co.condB[q] = -1.0;
        co.dchiFS[q] = -1.0;
        co.rhoFit[q] = 0;
        co.rhoMeas[q] = 0;
        co.rhoFallback[q] = 0;
        co.sigLnRhoData[q] = -1.0;
        for (int j = 0; j < Nx; ++j) {
            for (int k = 0; k < Nx; ++k) {
                gsl_matrix_set(co.inputA[q].get(), j, k, 0.0);
                gsl_matrix_set(co.inverA[q].get(), j, k, 0.0);
            }
        }
    }



    // A telescope with fewer than kMinTeleEpochs epochs is left out (see config/parameters.h), and so is each of
    // Rubin's bands with fewer than that: a band's own blend fraction and baseline cannot be fixed by one or two
    // visits (co.nRubinBand, 0 for a band left out; activePhotParams reads it). The model is evaluated one observer
    // at a time over all of its epochs, so the epochs are split by telescope here (stored order kept) and the results
    // put back at their stored index. epo/tim/nTele hold EVERY epoch of a telescope. pix/pt/nPh hold the epochs of the
    // photometric matrices; the astrometric ones take Roman's from epo/tim (every epoch, whatever its photometric count)
    // and Rubin's from pix/pt, so that a Rubin band left out of the photometry leaves the astrometry too.
    std::array<int, 2> nTele{0, 0};
    std::array<std::vector<int>, 2>    epo;  //stored index of each of the telescope's epochs
    std::array<std::vector<double>, 2> tim;  //and its time
    std::array<int, kRubinBands> nBand{};    //Rubin's recorded visits per band, u g r i z y
    for (int i = 0; i < ndw; ++i) {
        const int tt = int(l.tele[i]);//telescope[0,1] LSST, ELT
        CHECK(tt == 0 or tt == 1);
        nTele[tt] += 1;
        epo[tt].push_back(i);
        tim[tt].push_back(l.timn[i]);
        if (tt == 0) {
            CHECK(l.filt[i] >= 0 and l.filt[i] < kRubinBands);   //a Rubin epoch is in one of ugrizy
            nBand[l.filt[i]] += 1;
        } else {
            CHECK(l.filt[i] == 6);                               //a Roman epoch is in F146
        }
    }
    for (int b = 0; b < kRubinBands; ++b) co.nRubinBand[b] = (nBand[b] >= kMinTeleEpochs) ? nBand[b] : 0;
    const bool useL = nRubinBandsFitted(co) > 0, useR = nTele[1] >= kMinTeleEpochs;
    std::array<int, 2> nPh{0, 0};
    std::array<std::vector<int>, 2>    pix;  //stored index of each epoch of the photometric matrices
    std::array<std::vector<double>, 2> pt;   //and its time
    std::vector<char> inFit(ndw, 0);
    for (int i = 0; i < ndw; ++i) {
        const int tt = int(l.tele[i]);
        if (tt == 0 ? co.nRubinBand[l.filt[i]] == 0 : not useR) continue;
        inFit[i] = 1;
        nPh[tt] += 1;
        pix[tt].push_back(i);
        pt[tt].push_back(l.timn[i]);
    }

    // dm[i][slot] = d(model magnitude of epoch i)/d(theta_j), from the central stencil sig, the epoch's compact row
    // (kShared, then its own fb and mbs). Each derivative is needed once per epoch, however many F_jk it enters, so
    // the reference model and each perturbed model are evaluated once over a telescope's epochs (perturbation outer,
    // epoch inner). Of the nine parameters of a point source only u0, tE, piE, xi and t0 move the source relative to
    // the lens; the flux pair of the epoch's band, fb and mbs, enters m = mbs - 2.5 log10(fb A + 1 - fb) through the
    // flux alone, so its perturbed magnitudes come from the reference A. That is 1 + 2*5 = 11 model evaluations per
    // epoch whatever the number of bands: the magnification is achromatic, the bands differ only in the (fb, mbs) that
    // turn A into a magnitude. Re-deriving the second factor of each F_jk would take 1 + 2*9 + 2*45 = 109.
    // The tenth, ln rho, changes the shape of the magnification near the lens. Its derivative is the analytic
    // quadrupole term far from the lens and a finite difference of the library's disc table near it, so the
    // model is evaluated twice more only for the epochs near the lens, which most events have none of (below).
    std::vector<EpochDeriv> dm(ndw);
    std::vector<ModelPoint> ref, pert;
    const std::array<std::pair<int, double*>, 5> lcPar{{
        {0, &l.u0}, {1, &l.tE}, {3, &l.piE}, {4, &s.xi}, {5, &l.t0}}};   //slots 0-4 of kShared

    // Whether the model has a source size to fit: a uniform disc with a radius. Under a point source
    // (--source-model pspl) ln rho is not a parameter and dchiFS is not computed.
    const bool hasDisc = (vbmSourceModel() == SourceModel::ESPL) and s.rho > 0.0;
    const bool fitRho  = hasDisc and not co.rhoFixed;   //ln rho is a parameter of every partition with data
    std::array<double, 2> dchiT{0.0, 0.0};  //each telescope's share of dchiFS

    for (int tt = 0; tt < 2; ++tt) {
        const int n = nPh[tt];
        if (n == 0) continue;   //a telescope (or Rubin with no band) left out has no epochs here
        const int*    ix = pix[tt].data();
        const double* t  = pt[tt].data();
        ref.resize(n);  pert.resize(n);

        // Each epoch's own flux pair: its band's blend fraction and baseline magnitude. Rubin's r band and Roman's
        // F146 are s.blend/s.magb[2] and [6], bit for bit the s.fb/s.mbs[0] and [1] they were read from before.
        std::vector<double> fbE(n), mbsE(n), lgf(n);
        for (int m = 0; m < n; ++m) {
            fbE[m]  = s.blend[l.filt[ix[m]]];
            mbsE[m] = s.magb[l.filt[ix[m]]];
            CHECK(fbE[m] > 0.0);
            CHECK(fbE[m] <= 1.0);
        }

        // Unperturbed model magnitudes, recomputed under the current observer configuration rather
        // than taken from the cached l.magn[]. Each derivative is (model(theta + Delta) - reference)
        // / Delta, so the reference must share the observing geometry (as.satScale). When the two
        // differ, their difference leaks into every derivative as a constant dm/Delta; it cancels in
        // symmetric stencils but not in the cross-telescope rows, fb's outer bins, or any one-sided
        // stencil. Every model here is evaluated with the observer that produced these data.
        evaluateModel(s, l, as, tt, t, n, false, ref.data());
        for (int m = 0; m < n; ++m) lgf[m] = 2.5 * std::log10(ref[m].A * fbE[m] + 1.0 - fbE[m]);

        // The finite-source delta-chi2 of this telescope's epochs: the chi-squared by which they tell the disc
        // from a point source at the true parameters, sum ((m_disc - m_point) / sigma_m)^2. The two models
        // share the track (the source radius does not enter it), so the point source's magnification at an
        // epoch is magnifOf of the disc model's own u, bit for bit what the library gives for rho = 0
        // (tests/fisher_fixture.cpp --finite-source asserts it), and no second light curve is needed.
        if (hasDisc) {
            double d = 0.0;
            for (int m = 0; m < n; ++m) {
                CHECK(l.errm[ix[m]] > 0.0);
                const double dmag = 2.5 * std::log10(magnifOf(ref[m].u) * fbE[m] + 1.0 - fbE[m]) - lgf[m];
                const double r = dmag / l.errm[ix[m]];
                d += r * r;
            }
            dchiT[tt] = d;
        }

        // ln rho. Far from the lens (u >= rhoDerivAnalyticU(rho), RHO_DERIV_QUAD_U rho) the disc's magnification is the point source's plus a
        // quadrupole term exactly proportional to rho^2, so dA/d ln rho = 2 quad(u, rho) and no model is evaluated.
        // Nearer, the derivative is the finite difference of the library's disc magnification in ln rho (rho ->
        // rho exp(+-Delta); its table is bilinear in ln rho, rows 0.092 apart, so inside a cell it is the slope of the
        // cell holding rho), computed for those epochs only. The library leaves its table for the quadrupole at
        // quad ~ Tol / 10; rhoDerivAnalyticU keeps every difference on the table side of that switch, where the two
        // branches differ by ~Tol and a step of 1e-7 would turn the jump into a derivative a million times the true one.
        // The magnitude is the epoch's own band's: dm = conv dA / (fb A + 1 - fb) with that band's fb.
        if (fitRho) {
            const double uAn = rhoDerivAnalyticU(s.rho);
            std::vector<int> near;                            //this telescope's epochs that need the table
            for (int m = 0; m < n; ++m) {
                const double u = ref[m].u;
                if (u >= uAn) {
                    const double conv = -2.5 / std::log(10.0) * fbE[m];   //dm = conv dA / (fb A + 1 - fb)
                    const double dA = 2.0 * discQuadrupole(u, s.rho);
                    dm[ix[m]][kSlotRho] = conv * dA / (ref[m].A * fbE[m] + 1.0 - fbE[m]);
                } else near.push_back(m);
            }
            if (!near.empty()) {
                const int nn = int(near.size());
                std::vector<double> tn(nn);
                for (int k = 0; k < nn; ++k) tn[k] = t[near[k]];
                std::vector<ModelPoint> pn(nn);
                const double keep = s.rho;
                for (int h = 0; h < 2; ++h) {
                    const double diff = co.Delta1[9] * sig[h];
                    s.rho = keep * std::exp(diff);
                    CHECK(std::isfinite(s.rho) and s.rho > 0.0);
                    CHECK(diff != 0.0);
                    evaluateModel(s, l, as, tt, tn.data(), nn, false, pn.data());
                    s.rho = keep;
                    for (int k = 0; k < nn; ++k) {
                        const int m = near[k];
                        const double magw = mbsE[m] - 2.5 * std::log10(pn[k].A * fbE[m] + 1.0 - fbE[m]);
                        const double der  = (magw - (mbsE[m] - lgf[m])) / diff;
                        double& d = dm[ix[m]][kSlotRho];
                        d = (h == 0) ? der : (d + der) * 0.5;
                    }
                }
            }
        }

        // u0, tE, piE, xi, t0: one new model per perturbation.
        for (int slot = 0; slot < 5; ++slot) {
            const int j = lcPar[slot].first;
            double* par = lcPar[slot].second;
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
                    const double magw = mbsE[m] - 2.5 * std::log10(pert[m].A * fbE[m] + 1.0 - fbE[m]);
                    const double der  = (magw - (mbsE[m] - lgf[m])) / diff;
                    double& d = dm[ix[m]][slot];
                    d = (h == 0) ? der : (d + der) * 0.5;
                }
            }
        }

        // The flux pairs: one per Rubin band present among these epochs, or Roman's F146. A band's epochs share its
        // fb and mbs, so the pair's step and derivatives are computed once per band, over that band's epochs. A band
        // whose epochs all sit at baseline (A = 1) has dm/dfb proportional to (A - 1) = 0: an empty row, which the
        // prior on fb (below) holds invertible.
        for (int b = 0; b <= kRubinBands; ++b) {
            if ((b == kRubinBands) != (tt == 1)) continue;   //Rubin's bands for tt = 0, F146 (filter 6) for tt = 1
            std::vector<int> pos;
            for (int m = 0; m < n; ++m) if (l.filt[ix[m]] == b) pos.push_back(m);
            if (pos.empty()) continue;
            const int jfb = fbIndexOfFilter(b), jmbs = mbsIndexOfFilter(b);
            fluxPairDerivs(s.blend[b], s.magb[b], kFDStepScale * co.deltaScale[jfb], co.Delta1[jmbs], pos,
                           ref.data(), lgf.data(), ix, dm);
        }
    }

    // dchiFS of a partition is the sum of its telescopes' shares (the joint is Rubin + Roman; a telescope left out
    // above has none). A diagnostic of how strongly the partition's data see the disc; it gates nothing.
    if (hasDisc) {
        if (useL)         co.dchiFS[SRUBIN] = dchiT[0];
        if (useR)         co.dchiFS[SROMAN] = dchiT[1];
        if (useL or useR) co.dchiFS[SJOINT] = (useL ? dchiT[0] : 0.0) + (useR ? dchiT[1] : 0.0);
    }

    // Accumulate the information sum F_jk = sum_i (dm_i/dtheta_j)(dm_i/dtheta_k)/sigma_i^2 (lower
    // triangle, epochs in stored order). Each epoch feeds the joint matrix and exactly one
    // single-survey matrix, so F[SJOINT] == F[SRUBIN] + F[SROMAN] holds element by element.
    // Only an epoch's non-zero derivatives enter: its other flux pairs and, under a point source, ln rho are
    // exactly 0, and a product with an exact 0 adds nothing, so the sum is unchanged to the last bit. The entries are
    // taken in ascending parameter order and each product is d[j] d[k] / sigma^2, as it was over the full row.
    std::array<std::array<double, Nx * Nx>, NSURV> FA{};
    for (int i = 0; i < ndw; ++i) {//data
        if (not inFit[i]) continue;
        const int tt = int(l.tele[i]);
        const int surv = surveyOfTele(tt); // which single-survey matrix this epoch also feeds
        co.nepochA[SJOINT] += 1;
        co.nepochA[surv]   += 1;
        CHECK(l.errm[i] > 0.0);
        const auto& d = dm[i];
        std::array<int, 8>    idx;
        std::array<double, 8> val;
        int nz = 0;
        for (int sl = 0; sl < 8; ++sl) {
            const int j = (sl < 6) ? kShared[sl] : (sl == kSlotFb ? fbIndexOfFilter(l.filt[i]) : mbsIndexOfFilter(l.filt[i]));
            if (j == 9 and not fitRho) continue;   //ln rho is a parameter only under a disc
            if (d[sl] == 0.0) continue;            //structurally zero: adds nothing
            int a = nz++;
            while (a > 0 and idx[a - 1] > j) { idx[a] = idx[a - 1]; val[a] = val[a - 1]; --a; }   //keep ascending
            idx[a] = j;  val[a] = d[sl];
        }
        const double w2 = l.errm[i] * l.errm[i];
        for (int a = 0; a < nz; ++a) {
            for (int b = 0; b <= a; ++b) {
                const double contrib = val[a] * val[b] / w2;
                FA[SJOINT][idx[a] * Nx + idx[b]] += contrib;
                FA[surv][idx[a] * Nx + idx[b]]   += contrib;
            }
        }
    }//end of data for
    // ln rho is a parameter of every partition that has data (under a disc), with the Gaussian prior RHO_PRIOR_SIGLN
    // added once to each matrix after the data: the joint matrix carries one prior, not Rubin's plus Roman's. The
    // data part of F[9][9] is kept for a partition that has to give ln rho up (below).
    const double priorInfo = std::isfinite(co.rhoPriorSigLn) ? 1.0 / (co.rhoPriorSigLn * co.rhoPriorSigLn) : 0.0;
    std::array<double, NSURV> f99data{};
    for (int q = 0; q < NSURV; ++q) {
        co.rhoFit[q] = (fitRho and co.nepochA[q] > 0) ? 1 : 0;
        f99data[q] = FA[q][9 * Nx + 9];
        if (co.rhoFit[q]) FA[q][9 * Nx + 9] += priorInfo;
    }
    // The same for the blend fraction of every Rubin band in the fit, r included: sigma(f_S) = fbPriorSig, added once
    // to each matrix that holds the band (the joint and the Rubin-only one; Roman's F146 has none). It leaves a band
    // that never sees the bump (an empty fb row) invertible and moves the physical parameters not at all, since their
    // derivatives vanish at baseline as well; where the band does see the bump it is a weak 1/sigma^2 against the data.
    const double fbPriorInfo = std::isfinite(co.fbPriorSig) ? 1.0 / (co.fbPriorSig * co.fbPriorSig) : 0.0;
    for (int b = 0; b < kRubinBands; ++b) {
        if (co.nRubinBand[b] == 0) continue;
        const int jfb = fbIndexOfFilter(b);
        FA[SJOINT][jfb * Nx + jfb] += fbPriorInfo;
        FA[SRUBIN][jfb * Nx + jfb] += fbPriorInfo;
    }
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
        // A partition with fewer epochs than its ACTIVE data parameters (ln rho, held by its prior, not counted) is
        // rank-deficient (e.g. a short event peaking in a Roman gap has no Roman data): report it as not
        // characterizable rather than inverting it into a meaningless ~1e10 sigma.
        const bool enough = co.nepochA[q] >= static_cast<int>(activePhotParams(co, q).size()) - co.rhoFit[q];
        co.okA[q] = enough ? invert_matrix(co, 0, q) : 0;

        // If the matrix with ln rho does not invert, the partition gives it up and is forecast with rho held at its
        // true value: the nine-parameter inversion, whose matrix entries are the same ones. Counted, since it should
        // not happen: the prior keeps ln rho's own row well conditioned.
        if (co.rhoFit[q] and not co.okA[q]) {
            co.rhoFit[q] = 0;
            if (enough) {
                co.rhoFallback[q] = 1;
                ++gRhoFallback;
                gsl_matrix_set(co.inputA[q].get(), 9, 9, f99data[q]);
                co.okA[q] = invert_matrix(co, 0, q);
                if (not co.okA[q]) ++gRhoFallbackLost;
            }
        }

        // What the data alone say about ln rho. The prior adds 1/sp^2 to F[9][9] and nothing else, so the marginal
        // precision of ln rho is the data's plus the prior's, exactly: 1/sigma_post^2 = 1/sigma_data^2 + 1/sp^2.
        if (co.rhoFit[q]) {
            const double iPost = 1.0 / gsl_matrix_get(co.inverA[q].get(), 9, 9);
            const double iData = iPost - priorInfo;
            co.sigLnRhoData[q] = (iData > 1.0e-12 * iPost) ? 1.0 / std::sqrt(iData) : -1.0;
            co.rhoMeas[q] = (co.sigLnRhoData[q] > 0.0 and co.sigLnRhoData[q] <= co.sigLnRhoMax) ? 1 : 0;
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
    // go into one white group per band (only the bands in the photometric fit: a band with fewer than
    // kMinTeleEpochs visits is left out here too); Roman's go into day blocks (the sigma_c correlation
    // unit), tagged with season and roll so the variants can group them differently.
    // Each visit's centroid is diluted with the fractions of ITS band (f_S, f_L, f_N = 1 - f_S - f_L of the
    // baseline in that band). Each coordinate carries one extra parameter per telescope, slot Ny: the light
    // centroid of the unresolved neighbours in that coordinate (s.blendOff), unknown to the analyst; Rubin has
    // one per band (the neighbours' light centroid differs from band to band, and a shared offset would be
    // pinned down by colour differences a real analysis does not have), all sharing one reference position.
    // Its derivative is dcentroid/dblendOff = fn / (fs A + 1 - fs), largest at baseline and falling as the
    // source brightens; it is marginalised below, per telescope or Rubin band and per coordinate (it appears
    // in no other telescope's, band's or coordinate's data), when AST_BLEND_OFFSET_FREE.
    // The lens's own light (fraction f_L of the band's baseline, at the lens, -thetaE u from the source) is held at its
    // true value (co.lensLightFree false, the default), so its drift with the lens-source motion measures thetaE; or it is
    // a parameter of the fit, slot Ny + 1: one per Rubin band with astrometric visits (each band's sums are its own) and
    // one for Roman's F146, the same number in x and y. Its derivative is l = dcentroid/df_L = -thetaE u / (fs A + 1 - fs)
    // (ModelPoint::lensDx, lensDy at the reference; the neighbours' part -blendOff/den lies along the offset's derivative
    // and is absorbed by the free offset). A Gaussian prior of sigma = lensLightPriorFrac (1 - f_S) holds each fraction
    // in [0, 1 - f_S] (frac 1/sqrt(12): uniform on that range; 1: the whole blend; --lens-light-prior). A band with
    // no room (1 - f_S < 1e-6) has no such parameter. The photometric matrices and the detection are not touched.
    constexpr int D = fisher_detail::kAstD;
    using fisher_detail::AstSums;
    std::array<AstSums, kRubinBands> rubinBand;
    const int nact = co.lensLightFree ? D : D - 1;               //the lens-light slot is accumulated in free mode only
    std::map<long, AstSums> romanDays;

    // Per epoch i, ddx[i][j] and ddy[i][j]: the derivatives of the centroid's x and y, by the
    // central stencil sig for all four astrometric parameters, in slot Ny the neighbours' share
    // h and in slot Ny + 1 the lens's light fraction's derivative l (the reference's, free mode only). As for the photometry, the reference and each of the 2*Ny perturbed models are evaluated
    // once over a telescope's epochs: 1 + 2*4 = 9 model evaluations per epoch. Roman's every epoch enters,
    // Rubin's those of the bands of the photometric fit (inFit).
    std::vector<std::array<double, D>> ddx(ndw), ddy(ndw);
    std::vector<double> sx0, sy0;  //reference centroid without the neighbours' share, per epoch
    const std::array<double*, Ny> astPar{&l.tetE, &s.mus1, &s.mus2, &l.piE};

    for (int tt = 0; tt < 2; ++tt) {
        // Roman: every epoch. Rubin: the epochs of the bands in the photometric fit, each with its own filter.
        const int n = (tt == 0) ? nPh[0] : nTele[1];
        if (n == 0) continue;
        const int*    ix = (tt == 0) ? pix[0].data() : epo[1].data();
        const double* t  = (tt == 0) ? pt[0].data()  : tim[1].data();
        std::vector<int> filt;                                   //Rubin's epoch bands
        if (tt == 0) for (int m = 0; m < n; ++m) filt.push_back(l.filt[ix[m]]);
        const int* fp = (tt == 0) ? filt.data() : nullptr;
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
        // fs, fL of epoch m: the telescope's, or for Rubin the epoch's band's (the numbers measuredCentroid used).
        const bool   finite = (vbmSourceModel() == SourceModel::ESPL);
        auto nbrShare = [&](int m, const ModelPoint& p, double& bx, double& by) {
            const double fs = fp ? s.blend[fp[m]]  : s.fb[tt];
            const double fL = fp ? s.fLensF[fp[m]] : s.fLens[tt];
            const double fn = std::max(0.0, 1.0 - fs - fL);
            const double u2   = p.u * p.u;
            const double flux = finite ? fs * p.A + 1.0 - fs
                                       : fs * (u2 + 2.0) / std::sqrt(u2 * (u2 + 4.0)) + 1.0 - fs;
            const double h = fn / flux;
            bx = h * s.blendOff[tt][0];  by = h * s.blendOff[tt][1];
            return h;
        };

        // Reference and perturbations use the observer that produced these data.
        evaluateModel(s, l, as, tt, t, n, true, ref.data(), fp, co.lensLightFree);
        for (int m = 0; m < n; ++m) {
            double nbx0, nby0;
            ddx[ix[m]][Ny] = ddy[ix[m]][Ny] = nbrShare(m, ref[m], nbx0, nby0);
            ddx[ix[m]][Ny + 1] = ref[m].lensDx;                  //the lens's light fraction, read in free mode only
            ddy[ix[m]][Ny + 1] = ref[m].lensDy;
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
                evaluateModel(s, l, as, tt, t, n, true, pert.data(), fp);
                *astPar[j] = keep;
                // A null derivative (epoch at t0 or far from the peak) is legitimate and contributes nothing.
                for (int m = 0; m < n; ++m) {
                    double nbx, nby;
                    nbrShare(m, pert[m], nbx, nby);
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
            if (inFit[i]) rubinBand[l.filt[i]].add(w, ddx[i], ddy[i], nact);   //a band left out of the photometry is out here too
        } else {
            CHECK(l.rseas[i] >= 0 and l.rseas[i] < AST_MAX_SEASONS);
            CHECK(l.rroll[i] == 0 or l.rroll[i] == 1);
            AstSums& d = romanDays[long(std::floor(l.timn[i]))];
            if (d.season < 0) { d.season = l.rseas[i]; d.roll = l.rroll[i]; }
            CHECK(d.season == l.rseas[i] and d.roll == l.rroll[i]);
            d.add(w, ddx[i], ddy[i], nact);
        }
    }

    // The precision of the Gaussian prior on a band's lens-light fraction, 1 / (frac (1 - f_S))^2, or -1 where the band has
    // no room for lens light (no such parameter). Read in free mode only.
    const auto lensLightInfo = [&](double fS) {
        const double room = 1.0 - fS;
        if (co.lensLightFloor > 0.0) {                              //a floor on sigma: every band carries the parameter
            const double sig = std::max(co.lensLightPriorFrac * std::max(room, 0.0), co.lensLightFloor);
            return 1.0 / (sig * sig);
        }
        return room >= 1.0e-6 ? 1.0 / (co.lensLightPriorFrac * room * co.lensLightPriorFrac * room) : -1.0;
    };
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    std::array<double, Ny * Ny> FL{};
    {
        int kAst = 0, only = -1;                                   //Rubin's bands with astrometric epochs
        for (int b = 0; b < kRubinBands; ++b) if (rubinBand[b].Sw > 0.0) { ++kAst; only = b; }
        if (co.lensLightFree) {
            // Free fractions: every band's, whatever their number. A block that does not invert leaves NaN, so that the
            // partitions holding it are not characterised (and not silently better or worse).
            std::array<double, kRubinBands> flInfo;
            for (int b = 0; b < kRubinBands; ++b) flInfo[b] = lensLightInfo(s.blend[b]);
            if (kAst >= 1 and not fisher_detail::rubinBandInformationLensFree(rubinBand, flInfo, FL)) FL.fill(kNaN);
        } else if (kAst >= 2) {
            fisher_detail::rubinBandInformation(rubinBand, FL);
        } else {                                                    //no band, or one: its offset and the reference
            const AstSums none;
            const AstSums& rubin = (kAst == 1) ? rubinBand[only] : none;
            fisher_detail::AstGroup g;  g.nact = nact;  g.fold(rubin, 0.0);
            fisher_detail::AstPerCoord G{}; g.marginaliseInto(G);
            fisher_detail::astPhysical(G, {rubin.Swdd[0][Ny * D + Ny], rubin.Swdd[1][Ny * D + Ny]}, FL);
        }
    }
    const double romanFlInfo = co.lensLightFree ? lensLightInfo(s.fb[1]) : -1.0;   //-1: no f_L parameter in Roman's fit
    for (int v = 0; v < NAVAR; ++v) {
        std::array<double, Ny * Ny> FR{};
        fisher_detail::romanInformation(romanDays, v == AV_W, co.astSigc[v], romanFlInfo, FR);
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
// route's own fractional error; wp = 0 where rho is not measured in the partition). The piE used is the
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
  const bool rho    = (co.rhoMeas[q] != 0 and Era[9] >= 0.0);          //the photometric thetaE is used
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

  // The Einstein radius has two routes. Astrometry gives it directly (Erb[0]). When the data measure the source's
  // size (rhoMeas: the partition's data-only sigma(ln rho) is at most SIGLNRHO_MAX) the light curve gives it as
  // thetaE = theta* / rho, with the fractional error of ln rho (the posterior one) and of the source's angular
  // radius theta* (THETASTAR_SIGLN) added in quadrature. The two are combined by inverse variance, as independent measurements (they come from different
  // matrices, and theta* is independent of the fit). The lens mass below takes the combined thetaE; its error
  // carries the covariance of ln rho with the photometric piE, which share a matrix (lnCovTetPi), co.corrRhoPiE
  // being its correlation coefficient.
  for (int q = 0; q < NSURV; ++q) {
      const double rastr = (co.Erb[q][0] >= 0.0) ? co.Erb[q][0] / (std::fabs(l.tetE) + eps) : -1.0;
      const bool rhoUsed = co.rhoMeas[q] and co.okA[q] and co.Era[q][9] >= 0.0;
      co.relTetEph[q] = rhoUsed ? std::sqrt(co.Era[q][9] * co.Era[q][9] + THETASTAR_SIGLN * THETASTAR_SIGLN) : -1.0;
      co.relTetEc[q]  = combineRel(rastr, co.relTetEph[q]);
      co.corrRhoPiE[q] = std::numeric_limits<double>::quiet_NaN();
      if (rhoUsed and co.Era[q][3] > 0.0)
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
  // Era[2] is the -1 sentinel when Rubin's r band is not in the fit (no r visits, or fewer than kMinTeleEpochs): it must
  // not be divided into a negative number that reads as a measurement.
  co.resu[2]  = (Era[2] >= 0.0) ? double(Era[2] / (std::fabs(s.fb[0]) + eps)) : -1.0;
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
