// Fisher characterisation of a detected event, and the Step H3 no-satellite forecast.
#include "sim/characterize.h"
#include "sim/detect.h"
#include "fisher/fisher.h"
#include "surveys/noise.h"

void characterizeEvent(SimContext& ctx, SightlineState& st, EventState& ev) {
    const RunConfig& cfg = ctx.cfg;
    source& s = ctx.s;
    lens& l = ctx.l;
    astromet& as = ctx.as;
    lsst& ls = ctx.ls;
    covarian& co = ctx.co;
    covarian& coNS = ctx.coNS;
    const std::string& fnLDt = ctx.outs.fnLDt;

    double errg;

    for (int i = 0; i < nq; ++i) co.resu[i] = -1.0;
    // Same sentinel discipline for the per-survey results: an event that never
    // reaches FisherM must not inherit the previous event's sigmas.
    for (int q = 0; q < NSURV; ++q) {
        co.okA[q] = 0; co.okB[q] = 0; co.nepochA[q] = 0;
        co.flagi = 0;   // Deviation 78: was stale on uncharacterised events
        co.condA[q] = -1.0; co.condB[q] = -1.0;
        co.relMl[q] = -1.0;
        for (int k = 0; k < Nx; ++k) co.Era[q][k] = -1.0;
        for (int k = 0; k < Ny; ++k) co.Erb[q][k] = -1.0;
        for (int v = 0; v < NAVAR; ++v) {
            co.okBV[v][q] = 0; co.condBV[v][q] = -1.0; co.relMlV[v][q] = -1.0;
            for (int k = 0; k < Ny; ++k) co.ErbV[v][q][k] = -1.0;
        }
    }

    ev.FFG[0] = 0;
    ev.FFG[1] = 0;
    ev.FFG[2] = 0;
    ev.detL = 0; ev.detR = 0; ev.detJ = 0;

    //cout << "flagf: " << flagf << "ndw: " << ndw << endl;

    if (ev.flagf == 0 or ev.ndw <= 2) {
        errg    = errlsstM(s.magb[2], 2, double(RUBIN_R_DEPTH5_FALLBACK)); //r-band
        s.errA = errlsstA(ls, s.magb[2]); //r-band
        s.errM = std::fabs(std::pow(10.0, - 0.4 * errg) - 1.0); //r-band
        ev.dchiL = 0.0;
        ev.dchiP = 0.0;
        ev.dchiA = 0.0;
        ev.dchiL_L = 0.0; ev.dchiP_L = 0.0; ev.dchiA_L = 0.0;
        ev.dchiL_R = 0.0; ev.dchiP_R = 0.0; ev.dchiA_R = 0.0;
        ev.vsave = s.mus;
    }

    // Step H3's second forecast for this event. Declared one scope out from the
    // detection block below, because the row is written further down, where the
    // satellite observable and the coverage counts have been computed.
    // Step H3's no-satellite forecast for this same event. Declared here, before
    // the detection block, so the row write further down can see them whether or
    // not the Fisher step ran. -1.0 is the not-measured sentinel used everywhere
    // else in this code: a sigma that never inverted is not a large sigma.
    ev.sigtE_ns   = -1.0, ev.sigpiE_ns  = -1.0, ev.sigpiER_ns = -1.0;
    ev.sigtetE_ns = -1.0, ev.sigpiEb_ns = -1.0, ev.relMl_ns   = -1.0;
    ev.condA_ns   = -1.0, ev.condB_ns   = -1.0;
    ev.okNS       = 0,    ev.okNSb      = 0;

    if (ev.flagf > 0 and ev.ndw > 2) { //if star is visible
        cout << "************** DETECTABLE!!!!!! ********" << endl;
        st.icon +=1;
        ev.vsave = double(ev.vsave / (ev.ndw + 0.000065645));
        s.errM   = double(s.errM   / (ev.ndw + 0.000065645));
        s.errA   = double(s.errA   / (ev.ndw + 0.000065645));
        detectEvent(ctx, ev);

        if (ev.detL or ev.detR or ev.detJ) { //lensing — detected by Rubin, Roman, or the joint test
            ev.FFG[0] = 1; //Lensing
            st.nlens += ev.FFG[0];
            FisherM(s, l, as, co, ev.ndw);

            if (co.flagi > 0) {
                // flagi is always +1: FisherM's F*F^-1 checks are commented out
                // (Deviation 40). The conditioning test that works is okA.
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

            // -----------------------------------------------------------------
            // Step H3. The satellite-parallax experiment as a CONTROLLED
            // comparison on one event, instead of a difference between two runs.
            //
            // Why it must be done here. The obvious design -- two runs differing
            // only in L2_OFFSET_AU, matched by row -- cannot work. The RNG is one
            // stream and the per-event path draws conditionally on `acceptRubin
            // or acceptRoman`, so moving the observer changes what is detectable,
            // the streams fork at the first footprint sightline, and the two
            // runs' DETECTED populations then differ. Measured on the v3 pair,
            // that selection difference shifts the median sigma_tE -- which the
            // satellite cannot physically touch -- by 15%, several times the
            // effect being looked for. DEVIATIONS.md 35.
            //
            // Here the event is fixed: same draw, same epochs, same photometry,
            // only the observer moves. satScale is the single knob carrying
            // L2_OFFSET_AU into lightcurve(), and FisherM rebuilds the light
            // curve from it, so flipping it and recomputing gives the forecast
            // THIS event would have had with Roman at Earth.
            //
            // Re-evaluating is legitimate for an event whose data were generated
            // with the offset on: the Fisher matrix is built from model
            // derivatives at the true parameters, not from realised noise. The
            // question is what each observing geometry can constrain.
            if (cfg.pairSat) {
                const double keepScale = as.satScale;
                as.satScale = 0.0;
                FisherM(s, l, as, coNS, ev.ndw);
                if (coNS.flagi > 0) ErrorCal(coNS, l, s);
                as.satScale = keepScale;

                // Both halves of the forecast are recorded, because moving the
                // observer acts on both and they carry different physics.
                //
                //   PHOTOMETRIC (okA, Era): the magnification depends on |u|, and
                //   the two observers see different |u| at the same instant. That
                //   difference is the simultaneous parallax baseline, and it lands
                //   on piE.
                //
                //   ASTROMETRIC (okB, Erb): the centroid deflection is a vector,
                //   theta_E * u/(u^2+2), so moving the observer changes its
                //   direction as well as its size. Erb[0] is sigma(theta_E) and
                //   Erb[3] the astrometric route to piE, which is independent of
                //   the photometric one.
                //
                // relMl combines them: Ml = theta_E/(kappa piE) is the quantity a
                // microlensing survey actually wants, and it needs one observable
                // from each matrix, so it is the only place where a change in
                // either shows up as a change in the science.
                //
                // SROMAN as well as SJOINT: the L2 offset is Roman's geometry, so
                // Roman's own matrix is where any effect must appear first and
                // undiluted by Rubin's epochs.
                //
                // Condition numbers for both matrices and both observers, so the
                // conditioning question can be answered from the data instead of
                // hypothesised. This is the diagnostic OPEN_ITEMS asked for.
                ev.okNS  = coNS.okA[SJOINT];
                ev.okNSb = coNS.okB[SJOINT];
                ev.sigtE_ns   = ev.okNS  ? coNS.Era[SJOINT][1] : -1.0;
                ev.sigpiE_ns  = ev.okNS  ? coNS.Era[SJOINT][3] : -1.0;
                ev.sigpiER_ns = coNS.okA[SROMAN] ? coNS.Era[SROMAN][3] : -1.0;
                ev.sigtetE_ns = ev.okNSb ? coNS.Erb[SJOINT][0] : -1.0;
                ev.sigpiEb_ns = ev.okNSb ? coNS.Erb[SJOINT][3] : -1.0;
                ev.relMl_ns   = coNS.relMl[SJOINT];
                ev.condA_ns   = coNS.condA[SJOINT];
                ev.condB_ns   = coNS.condB[SJOINT];
            }
        }

        tallyDetection(ctx, st, ev);
    }
}

void writeSatellitePair(SimContext& ctx, const SightlineState& st, const EventState& ev) {
    const RunConfig& cfg = ctx.cfg;
    source& s = ctx.s;
    lens& l = ctx.l;
    covarian& co = ctx.co;
    const std::string& fnPair = ctx.outs.fnPair;

    // Step H3. One row per DETECTED event, carrying both forecasts for that same
    // event. Written here rather than beside the Fisher call because duSat and the
    // contemporaneous-coverage counts are computed above, and they are the axes
    // every H3 figure uses.
    if (cfg.pairSat and (ev.detL or ev.detR or ev.detJ)) {
        std::ofstream fpair(fnPair, std::ios::app);
        const int okAs = co.okA[SJOINT], okBs = co.okB[SJOINT];
        fpair << std::setprecision(7)
              << s.lon << " " << s.lat << " "
              << l.tE  << " " << l.u0  << " " << l.piE << " " << l.tetE << " "
              << ev.duSat << " "
              << okAs << " " << ev.okNS << " " << okBs << " " << ev.okNSb << " "
              << (okAs ? co.Era[SJOINT][1] : -1.0) << " " << ev.sigtE_ns   << " "
              << (okAs ? co.Era[SJOINT][3] : -1.0) << " " << ev.sigpiE_ns  << " "
              << (co.okA[SROMAN] ? co.Era[SROMAN][3] : -1.0) << " " << ev.sigpiER_ns << " "
              << (okBs ? co.Erb[SJOINT][0] : -1.0) << " " << ev.sigtetE_ns << " "
              << (okBs ? co.Erb[SJOINT][3] : -1.0) << " " << ev.sigpiEb_ns << " "
              << co.relMl[SJOINT] << " " << ev.relMl_ns << " "
              << co.condA[SJOINT] << " " << ev.condA_ns << " "
              << co.condB[SJOINT] << " " << ev.condB_ns << " "
              << ev.nepLpk << " " << ev.nepRpk << " " << st.wArea << " "
              << l.Ml << " " << l.Dl << " " << s.Ds << " " << l.Vt << "\n";
        fpair.close();
    }
}

