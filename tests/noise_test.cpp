// Unit test for the per-epoch astrometric error models (src/surveys/noise.cpp errRomanA, errlsstA),
// run against the vendored Roman photometric table files/sigma_roman.txt. Exit 0 = all held.
//
// Pinned: (1) Roman, sigma = sqrt((k/SNR)^2 + floor^2) with the SNR of errRomanM, against McKinnon &
// van der Marel 2026's tabulated 66-s GBTDS F146 curve (medium background, statistical part; their
// tool's data/roman_IM_66_6_gbtds_mid_5stripe_medium_F146_pos_errs.csv), floor added in quadrature;
// (2) Roman reaches its 1.1 mas floor for bright stars and never decreases with magnitude;
// (3) Rubin, sigma = sqrt((kappa FWHM / SNR)^2 + 10^2) with the SNR of errlsstM, at the 5-sigma
// depth (an independent evaluation of the formula, which pins units: FWHM in arcsec, sigma in mas),
// at the bright end (the floor), and its scaling with the visit's seeing and depth;
// (4) Rubin photometry, errlsstM = hypot(LSST_PHOT_FLOOR, sigma_rand) with sigma_rand written out from
// Ivezic et al. 2019 eq. 5, the bright end sitting on the floor, and errlsstA independent of the floor.
#include "common.h"
#include "surveys/noise.h"
#include "surveys/visits.h"
#include "run/inputs.h"
#include <cstdio>
#include <cmath>
#include <cstring>
#include <random>

static int nfail = 0;
static void expect(bool ok, const char* what, double got, double want) {
    std::printf("%-58s got %9.4f  want %9.4f  %s\n", what, got, want, ok ? "ok" : "FAIL");
    if (!ok) ++nfail;
}

// The oracle: errRomanM exactly as it was before the table logs were stored (Plan B, Task 2). Kept
// verbatim on purpose; the optimised errRomanM must return the same bits.
static double errRomanM_ref(const roman& ro, double mag)
{
    if (mag <= ro.mag[0]) return ro.err[0];
    int i = 1;
    while (i < NaRoman - 1 and ro.mag[i] < mag) ++i;
    const double f = (mag - ro.mag[i - 1]) / (ro.mag[i] - ro.mag[i - 1]);
    return std::exp(std::log(ro.err[i - 1]) + f * (std::log(ro.err[i]) - std::log(ro.err[i - 1])));
}

int main() {
    roman ro;
    if (readRomanErrorTable(ro) != 0) { std::printf("cannot read files/sigma_roman.txt\n"); return 1; }

    // (1) McKinnon & van der Marel 2026, F146 AB -> statistical position error [mas], medium background.
    const double mAB[]  = {20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0};
    const double mvdm[] = {0.421, 0.705, 1.254, 2.436, 5.191, 12.049, 29.065};
    for (int i = 0; i < 7; ++i) {
        const double got = errRomanA(errRomanM(ro, mAB[i]));
        const double want = std::hypot(mvdm[i], ROMAN_AST_FLOOR);
        char what[80]; std::snprintf(what, sizeof what, "Roman F146 AB %.0f vs McKinnon & vdM (4%%)", mAB[i]);
        expect(std::fabs(got / want - 1.0) < 0.04, what, got, want);
    }

    // (2) bright end and monotonicity
    expect(std::fabs(errRomanA(errRomanM(ro, 16.0)) / ROMAN_AST_FLOOR - 1.0) < 0.01,
           "Roman AB 16 sits on the 1.1 mas floor", errRomanA(errRomanM(ro, 16.0)), ROMAN_AST_FLOOR);
    bool mono = true;
    for (double m = 15.0; m < 26.0; m += 0.05)
        if (errRomanA(errRomanM(ro, m + 0.05)) < errRomanA(errRomanM(ro, m)) - 1e-12) mono = false;
    expect(mono, "Roman error never decreases with magnitude (15-26)", mono, 1.0);

    // (3) Rubin: a median r visit, depth 23.93, geometric FWHM 0.98".
    const double m5 = 23.93, fw = 0.98;
    // At m = m5 the photon-noise magnitude error of errlsstM is exactly 0.2 mag, so SNR = 1.0857/0.2.
    const double want5 = std::hypot(LSST_AST_KAPPA * fw * 1000.0 * 0.2 / 1.0857, LSST_AST_FLOOR);
    expect(std::fabs(errlsstA(errlsstM(m5, 2, m5), fw) / want5 - 1.0) < 1e-6,
           "Rubin r at the 5-sigma depth: kappa FWHM/SNR (+) 10 mas", errlsstA(errlsstM(m5, 2, m5), fw), want5);
    expect(std::fabs(errlsstA(errlsstM(16.0, 2, m5), fw) / LSST_AST_FLOOR - 1.0) < 0.005,
           "Rubin r = 16 sits on the 10 mas floor", errlsstA(errlsstM(16.0, 2, m5), fw), LSST_AST_FLOOR);
    const double s1 = errlsstA(errlsstM(23.0, 2, m5), fw), s2 = errlsstA(errlsstM(23.0, 2, m5), 2.0 * fw);
    const double stat1 = std::sqrt(s1 * s1 - LSST_AST_FLOOR * LSST_AST_FLOOR);
    const double stat2 = std::sqrt(s2 * s2 - LSST_AST_FLOOR * LSST_AST_FLOOR);
    expect(std::fabs(stat2 / stat1 - 2.0) < 1e-9, "Rubin statistical term doubles with the seeing FWHM", stat2 / stat1, 2.0);
    expect(errlsstA(errlsstM(23.0, 2, m5 - 0.5), fw) > s1, "Rubin error grows in a shallower visit",
           errlsstA(errlsstM(23.0, 2, m5 - 0.5), fw), s1);

    // (4) Rubin photometric error: sigma^2 = sigma_sys^2 + sigma_rand^2 (Ivezic et al. 2019 eq. 4-5),
    // sigma_rand^2 = (0.04 - gamma) x + gamma x^2, x = 10^(0.4 (m - m5)); r band (index 2), gamma 0.039.
    const double gam = 0.039;
    const double ms[] = {m5 - 8.3, m5 - 4.0, m5};
    const char* lab[] = {"m5 - 8.3", "m5 - 4", "m5"};
    for (int i = 0; i < 3; ++i) {
        const double x = std::pow(10.0, 0.4 * (ms[i] - m5));
        const double srand_ = std::sqrt((0.04 - gam) * x + gam * x * x);
        const double want = std::hypot(LSST_PHOT_FLOOR, srand_);
        const double got = errlsstM(ms[i], 2, m5);
        char what[80]; std::snprintf(what, sizeof what, "Rubin r errlsstM = hypot(floor, sigma_rand) at %s", lab[i]);
        expect(std::fabs(got / want - 1.0) < 1e-12, what, got, want);
    }
    const double eSat = errlsstM(m5 - RUBIN_SATU_BELOW_M5, 2, m5);
    expect(std::fabs(eSat / LSST_PHOT_FLOOR - 1.0) < 0.03, "Rubin r at saturation sits on the photometric floor (3%)",
           eSat, LSST_PHOT_FLOOR);
    // errlsstA takes the photon-noise part of errlsstM, so the floor must drop out: compare with the
    // astrometric error computed from sigma_rand directly.
    for (int i = 0; i < 3; ++i) {
        const double x = std::pow(10.0, 0.4 * (ms[i] - m5));
        const double srand_ = std::sqrt((0.04 - gam) * x + gam * x * x);
        const double stat = LSST_AST_KAPPA * fw * 1000.0 * srand_ / 1.0857;
        const double want = std::hypot(stat, LSST_AST_FLOOR);
        const double got = errlsstA(errlsstM(ms[i], 2, m5), fw);
        char what[80]; std::snprintf(what, sizeof what, "Rubin r errlsstA independent of the photometric floor at %s", lab[i]);
        expect(std::fabs(got / want - 1.0) < 1e-9, what, got, want);
    }

    // (5) errRomanM with stored table logs is bit-identical to the pre-change body.
    {
        bool inc = true, logsOk = (ro.logErr.size() == ro.err.size());
        for (int i = 1; i < NaRoman; ++i) if (!(ro.mag[i] > ro.mag[i - 1])) inc = false;
        for (int i = 0; logsOk and i < NaRoman; ++i) {
            const double want = std::log(ro.err[i]);
            if (std::memcmp(&ro.logErr[i], &want, sizeof want) != 0) logsOk = false;
        }
        expect(inc, "Roman table magnitudes strictly increasing", inc, 1.0);
        expect(logsOk, "Roman logErr[i] == log(err[i]) bitwise", logsOk, 1.0);
        long n = 0, bad = 0;
        auto cmp = [&](double m) {
            const double a = errRomanM(ro, m), b = errRomanM_ref(ro, m);
            ++n;
            if (std::memcmp(&a, &b, sizeof a) != 0) ++bad;
        };
        for (int i = 0; i < NaRoman; ++i) {
            cmp(ro.mag[i]);
            if (i + 1 < NaRoman) cmp(0.5 * (ro.mag[i] + ro.mag[i + 1]));
        }
        std::mt19937_64 gen(20261113);
        std::uniform_real_distribution<double> U(ro.mag[0] - 2.0, ro.mag[NaRoman - 1] + 2.0);
        for (int k = 0; k < 1000000; ++k) cmp(U(gen));
        std::printf("errRomanM vs reference: %ld comparisons, %ld differ\n", n, bad);
        expect(bad == 0, "errRomanM bit-identical to the pre-change body", (double)bad, 0.0);
    }

    std::printf("%s (%d failed)\n", nfail ? "FAILED" : "all held", nfail);
    return nfail ? 1 : 0;
}
