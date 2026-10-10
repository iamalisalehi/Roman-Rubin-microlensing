// The per-event output row, the satellite-parallax quantities and the sample light-curve dump.
#include "sim/record.h"
#include "sim/characterize.h"
#include "events/lightcurve.h"
#include "run/sample_dump.h"
#include "fisher/fisher.h"

void recordEvent(SimContext& ctx, SightlineState& st, const EfficiencyBins& bins, const LightCurveStats& lc,
                 const Characterization& ch) {
    const RomanSchedule& sched = ctx.sched;
    source& s = ctx.s;
    lens& l = ctx.l;
    astromet& as = ctx.as;
    covarian& co = ctx.co;
    std::ofstream& filg_in = ctx.outs.filg_in;
    const std::string& testf = ctx.outs.testf;
    PeakCoverage pk;

    // The observed peak, from which the gap geometry is measured.
    std::tie(pk.t0obs, pk.uminObs) = observedPeak(s, l, as);

    st.records.push_back(EventRecord{
        st.icon, static_cast<int>(ch.fisher),
        l.tE, l.RE/AU, l.piE, l.tetE, l.Vt, l.u0, l.Ml,
        s.opt*1.0e6, l.Dl, s.Ds, l.vl, s.vs,
        s.mbs[0], s.fb[0],
        bins.gg, static_cast<int>(l.struc),
        s.FWHM/year, ch.vMean/s.mus, l.DeltaT/s.errA, l.murel*year,
        co.resu[0], co.resu[1], co.resu[2], co.resu[3], co.resu[5],
        co.resu[9], co.resu[10], co.resu[13], co.resu[14],
        s.Map[2], s.nsbl[2],
        co.flagi, s.Ai[2],
        lc.ndw_L, lc.ndw_R,
        ch.det.detL, ch.det.detR, ch.det.detJ,
        co.okA[SJOINT], co.okA[SRUBIN], co.okA[SROMAN],
        co.Era[SJOINT][1], co.Era[SRUBIN][1], co.Era[SROMAN][1],   //sigma(tE)
        co.Era[SJOINT][3], co.Era[SRUBIN][3], co.Era[SROMAN][3],   //sigma(piE)
        co.Erb[SJOINT][0], co.Erb[SRUBIN][0], co.Erb[SROMAN][0],   //sigma(tetE)
        ch.det.dclsEvent,
        synergyClass(co),
        co.condA[SJOINT], co.condA[SRUBIN], co.condA[SROMAN],
        // ---- columns appended after the original layout ----
        l.t0, s.xi, s.lon, s.lat,
        s.mbs[1], s.fb[1],
        {s.magb[0], s.magb[1], s.magb[2], s.magb[3], s.magb[4], s.magb[5], s.magb[6]},
        {s.blend[0], s.blend[1], s.blend[2], s.blend[3], s.blend[4], s.blend[5], s.blend[6]},
        co.relMl[SJOINT], co.relMl[SRUBIN], co.relMl[SROMAN],
        co.okB[SJOINT], co.okB[SRUBIN], co.okB[SROMAN],
        co.condB[SJOINT], co.condB[SRUBIN], co.condB[SROMAN],
        sched.dtToSeasonEdge(pk.t0obs), sched.zone(pk.t0obs),
        lc.nres5_L, lc.nres20_L, lc.nresPSF_L, lc.nres5_R, lc.nres20_R, lc.nresPSF_R,
        lc.dsepMax_L, lc.dsepMax_R
    });

    // du_sat: observer separation in Einstein radii at t0, piE * D_perp / AU. D_perp is the separation
    // projected perpendicular to the line of sight (~0.87-0.99 of the L2 offset over the year), so
    // L2_OFFSET_AU * piE is only a ceiling; lightcurve() is asked for both observers so the projection
    // is not re-derived here. nepL_pk / nepR_pk: epochs of each survey within +-2 tE of t0; satellite
    // parallax needs contemporaneous coverage. Calling lightcurve() here is safe: FisherM has run, and
    // the state it touches is not read by the row below and is recomputed by the next draw.
    pk.duSat = 0.0;
    {
        lightcurve(s, l, as, l.t0, 0);
        const double r1 = as.ue_n1, r2 = as.ue_n2;
        lightcurve(s, l, as, l.t0, 1);
        const double dn1 = as.ue_n1 - r1, dn2 = as.ue_n2 - r2;
        pk.duSat = l.piE * std::sqrt(dn1 * dn1 + dn2 * dn2);
    }
    pk.nepLpk = 0, pk.nepRpk = 0;
    {
        const double win = 2.0 * l.tE;
        for (int i = 0; i < lc.ndw; ++i) {
            if (std::fabs(l.timn[i] - pk.t0obs) > win) continue;
            if (int(l.tele[i]) == 1) pk.nepRpk += 1;
            else                      pk.nepLpk += 1;
        }
    }

    filg_in.open(testf, std::ios::app);
    filg_in << st.icon           << " " << ch.fisher         << " " << l.tE                      << " "
            << l.RE / AU     << " " << l.piE         << " " << l.tetE                    << " "
            << l.Vt          << " " << l.u0          << " " << l.Ml                      << " " 
            << s.opt * 1.0e6 << " " << l.Dl          << " " << s.Ds                      << " "
            << l.vl          << " " << s.vs          << " " << s.mbs[0]                  << " "
            << s.fb[0]       << " " << bins.gg             << " " << static_cast<int>(l.struc) << " "
            << s.FWHM / year << " " << ch.vMean / s.mus << " " << l.DeltaT / s.errA        << " " << l.murel * year << " "
            << co.resu[0]    << " " << co.resu[1]    << " " << co.resu[2]                << " " 
            << co.resu[3]    << " " << co.resu[5]    << " "
            << co.resu[9]    << " " << co.resu[10]   << " " << co.resu[13]               << " " << co.resu[14]    << " "
            << s.Map[2]      << " " << s.nsbl[2]     << " " << co.flagi                  << " " << s.Ai[2]        << " "
            // per-survey bookkeeping, appended so existing column indices hold
            << lc.ndw_L << " " << lc.ndw_R << " "
            << ch.det.detL  << " " << ch.det.detR  << " " << ch.det.detJ << " "
            << co.okA[SJOINT] << " " << co.okA[SRUBIN] << " " << co.okA[SROMAN] << " "
            << co.Era[SJOINT][1] << " " << co.Era[SRUBIN][1] << " " << co.Era[SROMAN][1] << " "
            << co.Era[SJOINT][3] << " " << co.Era[SRUBIN][3] << " " << co.Era[SROMAN][3] << " "
            << co.Erb[SJOINT][0] << " " << co.Erb[SRUBIN][0] << " " << co.Erb[SROMAN][0] << " "
            << ch.det.dclsEvent << " " << synergyClass(co) << " "
            << co.condA[SJOINT] << " " << co.condA[SRUBIN] << " " << co.condA[SROMAN] << " "
            // appended after column 57 so earlier column indices hold
            << l.t0 << " " << s.xi << " " << s.lon << " " << s.lat << " "
            << s.mbs[1] << " " << s.fb[1] << " ";
    for (int i = 0; i < M; ++i) filg_in << s.magb[i]  << " ";
    for (int i = 0; i < M; ++i) filg_in << s.blend[i] << " ";
    filg_in << co.relMl[SJOINT] << " " << co.relMl[SRUBIN] << " " << co.relMl[SROMAN] << " "
            << co.okB[SJOINT]   << " " << co.okB[SRUBIN]   << " " << co.okB[SROMAN]   << " "
            << co.condB[SJOINT] << " " << co.condB[SRUBIN] << " " << co.condB[SROMAN] << " "
            // Gap geometry: dt_edge is negative when t0 fell inside a Roman season; t0zone separates a
            // mid-mission gap (1) from before-launch/after-end (2), and only the former is gap-filling.
            << sched.dtToSeasonEdge(pk.t0obs) << " " << sched.zone(pk.t0obs) << " "
            // Sky area this sightline stands for [deg^2]; not constant once --stride-roman is used, so
            // statistics pooled over sightlines must weight by it.
            << st.wArea << " "
            // satellite-parallax observable and contemporaneous coverage
            << pk.duSat << " " << pk.nepLpk << " " << pk.nepRpk << " "
            // Image resolution: counts of qualifying epochs per survey, then the largest separation
            // while both images were detectable (-1 = never). The three bars differ in what counts as resolved.
            << lc.nres5_L << " " << lc.nres20_L << " " << lc.nresPSF_L << " " << lc.dsepMax_L << " "
            << lc.nres5_R << " " << lc.nres20_R << " " << lc.nresPSF_R << " " << lc.dsepMax_R << " "
            // astrometric noise variants N and P (joint, Roman)
            << co.ErbV[AV_N][SJOINT][0] << " " << co.ErbV[AV_N][SROMAN][0] << " "
            << co.ErbV[AV_P][SJOINT][0] << " " << co.ErbV[AV_P][SROMAN][0] << " "
            << co.relMlV[AV_N][SJOINT] << " " << co.relMlV[AV_N][SROMAN] << " "
            << co.relMlV[AV_P][SJOINT] << " " << co.relMlV[AV_P][SROMAN] << " "
            << co.okBV[AV_N][SJOINT] << " " << co.okBV[AV_N][SROMAN] << " "
            << co.okBV[AV_P][SJOINT] << " " << co.okBV[AV_P][SROMAN] << " "
            << int(l.luminous) << " " << s.fLens[0] << " " << s.fLens[1] << " "
            << pk.t0obs << " " << pk.uminObs << " "
            << std::hypot(s.blendOff[0][0], s.blendOff[0][1]) << " "
            << std::hypot(s.blendOff[1][0], s.blendOff[1][1]) << " "
            // finite-source size: radius [Rsun], angular radius [mas], and rho = theta* / thetaE
            << s.Rstar << " " << s.thetaStar << " " << s.rho << "\n";
    filg_in.close();

    commitSampleDump(ctx, lc, ch, pk);
    writeSatellitePair(ctx, st, ch, pk);
}

void commitSampleDump(SimContext& ctx, const LightCurveStats& lc, const Characterization& ch,
                      const PeakCoverage& pk) {
    const RomanSchedule& sched = ctx.sched;
    source& s = ctx.s;
    lens& l = ctx.l;
    astromet& as = ctx.as;
    covarian& co = ctx.co;
    SampleSpec& dumpSpec = ctx.outs.dumpSpec;
    std::vector<DumpEpoch>& dumpBuf = ctx.outs.dumpBuf;
    long& dumpSeq = ctx.outs.dumpSeq;

    // Commit this draw's buffered light curve if it fills a requested sample class. Done here, not
    // in the time loop, because the classes are defined by detL/detR/detJ and the Fisher sigmas.
    // First match wins: an event is written once, under the first class in the spec file whose quota
    // is still open.
    if (dumpSpec.on and not dumpBuf.empty()) {
        double maxShift = 0.0;
        for (const DumpEpoch& e : dumpBuf)
            maxShift = std::max(maxShift,
                                std::sqrt(e.def1c * e.def1c + e.def2c * e.def2c));

        const SampleFacts facts{ch.det.detL, ch.det.detR, ch.det.detJ, l.tE, sched.zone(pk.t0obs),
                                pk.nepLpk, pk.nepRpk, lc.ndw_L, lc.ndw_R,
                                co.okB[SROMAN], maxShift};

        for (SampleClass& c : dumpSpec.classes) {
            if (c.kept >= c.quota)        continue;
            if (!sampleMatch(c, facts))   continue;

            std::ostringstream pr;
            pr << std::setprecision(10)
               << "# Sample event for class '" << c.name << "'.\n"
               << "# One key per line. Angles in deg, times in days, masses in Msun,\n"
               << "# distances in kpc, angular scales in mas, Rstar in Rsun. A sigma of -1 means the\n"
               << "# parameter was NOT measured by that survey (inactive or no epochs),\n"
               << "# never that it was measured to be -1.\n"
               << "class "        << c.name        << "\n"
               << "population "   << gPop->name    << "\n"
               << "lon "          << s.lon        << "\n"
               << "lat "          << s.lat        << "\n"
               << "tE "           << l.tE         << "\n"
               << "t0 "           << l.t0         << "\n"
               << "u0 "           << l.u0         << "\n"
               << "xi "           << s.xi         << "\n"
               << "piE "          << l.piE        << "\n"
               << "tetE "         << l.tetE       << "\n"
               << "Ml "           << l.Ml         << "\n"
               << "Dl "           << l.Dl         << "\n"
               << "Ds "           << s.Ds         << "\n"
               << "Rstar "        << s.Rstar      << "\n"
               << "thetaStar "    << s.thetaStar  << "\n"
               << "rho "          << s.rho        << "\n"
               << "Vt "           << l.Vt         << "\n"
               << "murel_yr "     << l.murel * year << "\n"
               << "mus1 "         << s.mus1       << "\n"
               << "mus2 "         << s.mus2       << "\n"
               << "mul1 "         << l.mul1       << "\n"
               << "mul2 "         << l.mul2       << "\n"
               << "lens_struc "   << int(l.struc) << "\n"
               << "mbs0 "         << s.mbs[0]     << "\n"
               << "fb0 "          << s.fb[0]      << "\n"
               << "mbs1 "         << s.mbs[1]     << "\n"
               << "fb1 "          << s.fb[1]      << "\n";
            pr << "magb";  for (int i = 0; i < M; ++i) pr << " " << s.magb[i];
            pr << "\nblend"; for (int i = 0; i < M; ++i) pr << " " << s.blend[i];
            pr << "\n"
               << "ndw_L "        << lc.ndw_L         << "\n"
               << "ndw_R "        << lc.ndw_R         << "\n"
               << "nep_pk_L "     << pk.nepLpk        << "\n"
               << "nep_pk_R "     << pk.nepRpk        << "\n"
               << "detL "         << ch.det.detL          << "\n"
               << "detR "         << ch.det.detR          << "\n"
               << "detJ "         << ch.det.detJ          << "\n"
               << "dt_edge "      << sched.dtToSeasonEdge(pk.t0obs) << "\n"
               << "t0zone "       << sched.zone(pk.t0obs)           << "\n"
               << "du_sat "       << pk.duSat         << "\n"
               << "max_shift "    << maxShift      << "\n"
               << "dsep_max_L "   << lc.dsepMax_L     << "\n"
               << "dsep_max_R "   << lc.dsepMax_R     << "\n"
               << "okA_J "  << co.okA[SJOINT] << " okA_L " << co.okA[SRUBIN]
               << " okA_R " << co.okA[SROMAN] << "\n"
               << "okB_J "  << co.okB[SJOINT] << " okB_L " << co.okB[SRUBIN]
               << " okB_R " << co.okB[SROMAN] << "\n"
               << "sigtE_J "   << co.Era[SJOINT][1] << " sigtE_L " << co.Era[SRUBIN][1]
               << " sigtE_R "  << co.Era[SROMAN][1] << "\n"
               << "sigpiE_J "  << co.Era[SJOINT][3] << " sigpiE_L " << co.Era[SRUBIN][3]
               << " sigpiE_R " << co.Era[SROMAN][3] << "\n"
               << "sigtetE_J "  << co.Erb[SJOINT][0] << " sigtetE_L " << co.Erb[SRUBIN][0]
               << " sigtetE_R " << co.Erb[SROMAN][0] << "\n"
               << "relMl_J "  << co.relMl[SJOINT] << " relMl_L " << co.relMl[SRUBIN]
               << " relMl_R " << co.relMl[SROMAN] << "\n";

            std::ostringstream idss;
            idss << std::setw(3) << std::setfill('0') << (++dumpSeq);
            writeSampleEvent(dumpSpec, c.name, idss.str(), dumpBuf, pr.str(),
                             s, l, as);
            c.kept += 1;
            std::cout << "  [sample] " << c.name << " " << idss.str()
                      << "  tE=" << l.tE << "d  Ml=" << l.Ml << "Msun"
                      << "  ndw_L=" << lc.ndw_L << " ndw_R=" << lc.ndw_R
                      << "  (" << c.kept << "/" << c.quota << ")" << std::endl;
            break;
        }
    }
}

