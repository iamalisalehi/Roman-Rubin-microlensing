// Fisher characterisation of a detected event, and the no-satellite forecast.
#include "sim/characterize.h"
#include "sim/detect.h"
#include "fisher/fisher.h"
#include "surveys/noise.h"

Characterization characterizeEvent(SimContext& ctx, SightlineState& st, const EfficiencyBins& bins,
                                   const LightCurveStats& lc) {
    const RunConfig& cfg = ctx.cfg;
    source& s = ctx.s;
    lens& l = ctx.l;
    astromet& as = ctx.as;
    covarian& co = ctx.co;
    covarian& coNS = ctx.coNS;
    const std::string& fnLDt = ctx.outs.fnLDt;

    double errg;
    Characterization ch;

    for (int i = 0; i < nq; ++i) co.resu[i] = -1.0;
    // Per-survey results start at the not-measured sentinel so an event that never reaches
    // FisherM cannot inherit the previous event's sigmas.
    for (int q = 0; q < NSURV; ++q) {
        co.okA[q] = 0; co.okB[q] = 0; co.nepochA[q] = 0;
        co.flagi = 0;
        co.condA[q] = -1.0; co.condB[q] = -1.0;
        co.relMl[q] = -1.0;
        for (int k = 0; k < Nx; ++k) co.Era[q][k] = -1.0;
        for (int k = 0; k < Ny; ++k) co.Erb[q][k] = -1.0;
        for (int v = 0; v < NAVAR; ++v) {
            co.okBV[v][q] = 0; co.condBV[v][q] = -1.0; co.relMlV[v][q] = -1.0;
            for (int k = 0; k < Ny; ++k) co.ErbV[v][q][k] = -1.0;
        }
    }

    if (lc.flagf == 0 or lc.ndw <= 2) {
        errg    = errlsstM(s.magb[2], 2, double(RUBIN_R_DEPTH5_FALLBACK)); //r-band
        s.errA = errlsstA(errg, FWHM[2]); //r-band, a median r visit
        s.errM = std::fabs(std::pow(10.0, - 0.4 * errg) - 1.0); //r-band
        ch.vMean = s.mus;   // the delta-chi-squared statistics keep their zero defaults
    }

    // The no-satellite forecast of this event is ch.pair; its defaults are the -1.0 not-measured
    // sentinel and it is filled below only if the Fisher step runs.

    if (lc.flagf > 0 and lc.ndw > 2) { //if star is visible
        cout << "************** DETECTABLE!!!!!! ********" << endl;
        st.icon +=1;
        ch.vMean = double(lc.vsave / (lc.ndw + 0.000065645));
        s.errM   = double(s.errM   / (lc.ndw + 0.000065645));
        s.errA   = double(s.errA   / (lc.ndw + 0.000065645));
        ch.det = detectEvent(ctx, lc);

        if (ch.det.detL or ch.det.detR or ch.det.detJ) { //lensing — detected by Rubin, Roman, or the joint test
            ch.fisher = 1; //Lensing
            st.nlens += ch.fisher;
            FisherM(s, l, as, co, lc.ndw);

            if (co.flagi > 0) {
                // flagi is always +1 (FisherM's F*F^-1 checks are disabled); the effective
                // conditioning test is okA.
                st.nerr += co.okA[SJOINT] ? 1.0 : 0.0;
                ErrorCal(co, l, s);

                std::ofstream fil0_append(fnLDt, std::ios::app);
                fil0_append << std::fixed << std::setprecision(5)
                            << l.Ml       << " " << l.Dl       << " " << l.mul       << " " << s.fb[0]     << " " << s.mbs[0]   << " "
<<std::setprecision(7)           << l.tE       << " " << l.murel    << " " << l.u0        << " " << s.lon       << " " << s.lat      << " "
                            << l.piE      << " " << l.tetE     << " " << co.resu[1]  << " " << co.resu[2]  << " " << co.resu[3] << " "
                            << co.resu[5] << " " << co.resu[9] << " " << co.resu[10] << " " << co.resu[13] << "\n";
                fil0_append.close();
            }

            // No-satellite forecast: the same event (same draw, epochs and photometry) re-evaluated with
            // the observer moved to Earth. A controlled comparison of this kind is needed because two
            // runs differing only in L2_OFFSET_AU fork the RNG stream and select different detected
            // populations (shifting the median sigma_tE by ~15%). satScale carries L2_OFFSET_AU into
            // lightcurve(); the Fisher matrix uses model derivatives at the true parameters, not
            // realised noise, so re-evaluating is legitimate.
            if (cfg.pairSat) {
                const double keepScale = as.satScale;
                as.satScale = 0.0;
                FisherM(s, l, as, coNS, lc.ndw);
                if (coNS.flagi > 0) ErrorCal(coNS, l, s);
                as.satScale = keepScale;

                // Both halves are recorded. Photometric (okA, Era): the two observers see different |u|
                // at the same instant, a simultaneous parallax baseline that lands on piE. Astrometric
                // (okB, Erb): the deflection theta_E*u/(u^2+2) is a vector, so moving the observer changes
                // its direction too; Erb[0] is sigma(theta_E), Erb[3] the astrometric piE. relMl
                // (Ml = theta_E/(kappa piE)) combines them. SROMAN is recorded as well as SJOINT since
                // the L2 offset is Roman's geometry. Condition numbers of both matrices are kept.
                ch.pair.okNS  = coNS.okA[SJOINT];
                ch.pair.okNSb = coNS.okB[SJOINT];
                ch.pair.sigtE_ns   = ch.pair.okNS  ? coNS.Era[SJOINT][1] : -1.0;
                ch.pair.sigpiE_ns  = ch.pair.okNS  ? coNS.Era[SJOINT][3] : -1.0;
                ch.pair.sigpiER_ns = coNS.okA[SROMAN] ? coNS.Era[SROMAN][3] : -1.0;
                ch.pair.sigtetE_ns = ch.pair.okNSb ? coNS.Erb[SJOINT][0] : -1.0;
                ch.pair.sigpiEb_ns = ch.pair.okNSb ? coNS.Erb[SJOINT][3] : -1.0;
                ch.pair.relMl_ns   = coNS.relMl[SJOINT];
                ch.pair.condA_ns   = coNS.condA[SJOINT];
                ch.pair.condB_ns   = coNS.condB[SJOINT];
            }
        }

        tallyDetection(ctx, st, bins, lc, ch.det);
    }
    return ch;
}

void writeSatellitePair(SimContext& ctx, const SightlineState& st, const Characterization& ch,
                        const PeakCoverage& pk) {
    const RunConfig& cfg = ctx.cfg;
    source& s = ctx.s;
    lens& l = ctx.l;
    covarian& co = ctx.co;
    const std::string& fnPair = ctx.outs.fnPair;

    // One row per detected event with both forecasts. Written here rather than beside the Fisher
    // call because duSat and the peak-coverage counts (the comparison axes) are computed first.
    if (cfg.pairSat and (ch.det.detL or ch.det.detR or ch.det.detJ)) {
        std::ofstream fpair(fnPair, std::ios::app);
        const int okAs = co.okA[SJOINT], okBs = co.okB[SJOINT];
        fpair << std::setprecision(7)
              << s.lon << " " << s.lat << " "
              << l.tE  << " " << l.u0  << " " << l.piE << " " << l.tetE << " "
              << pk.duSat << " "
              << okAs << " " << ch.pair.okNS << " " << okBs << " " << ch.pair.okNSb << " "
              << (okAs ? co.Era[SJOINT][1] : -1.0) << " " << ch.pair.sigtE_ns   << " "
              << (okAs ? co.Era[SJOINT][3] : -1.0) << " " << ch.pair.sigpiE_ns  << " "
              << (co.okA[SROMAN] ? co.Era[SROMAN][3] : -1.0) << " " << ch.pair.sigpiER_ns << " "
              << (okBs ? co.Erb[SJOINT][0] : -1.0) << " " << ch.pair.sigtetE_ns << " "
              << (okBs ? co.Erb[SJOINT][3] : -1.0) << " " << ch.pair.sigpiEb_ns << " "
              << co.relMl[SJOINT] << " " << ch.pair.relMl_ns << " "
              << co.condA[SJOINT] << " " << ch.pair.condA_ns << " "
              << co.condB[SJOINT] << " " << ch.pair.condB_ns << " "
              << pk.nepLpk << " " << pk.nepRpk << " " << st.wArea << " "
              << l.Ml << " " << l.Dl << " " << s.Ds << " " << l.Vt << "\n";
        fpair.close();
    }
}

