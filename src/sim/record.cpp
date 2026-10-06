// The per-event output row, the Step H2 quantities and the Step S1 sample dump.
#include "sim/record.h"
#include "sim/characterize.h"
#include "events/lightcurve.h"
#include "run/sample_dump.h"
#include "fisher/fisher.h"

void recordEvent(SimContext& ctx, SightlineState& st, EventState& ev) {
    const RomanSchedule& sched = ctx.sched;
    source& s = ctx.s;
    lens& l = ctx.l;
    astromet& as = ctx.as;
    covarian& co = ctx.co;
    std::ofstream& filg_in = ctx.outs.filg_in;
    const std::string& testf = ctx.outs.testf;

    // Deviation 76: the observed peak, and the gap geometry measured from it.
    std::tie(ev.t0obs, ev.uminObs) = observedPeak(s, l, as);

    st.records.push_back(EventRecord{
        st.icon, static_cast<int>(ev.FFG[0]),
        l.tE, l.RE/AU, l.piE, l.tetE, l.Vt, l.u0, l.Ml,
        s.opt*1.0e6, l.Dl, s.Ds, l.vl, s.vs,
        s.mbs[0], s.fb[0],
        ev.gg, static_cast<int>(l.struc),
        s.FWHM/year, ev.vsave/s.mus, l.DeltaT/s.errA, l.murel*year,
        co.resu[0], co.resu[1], co.resu[2], co.resu[3], co.resu[5],
        co.resu[9], co.resu[10], co.resu[13], co.resu[14],
        s.Map[2], s.nsbl[2],
        co.flagi, s.Ai[2],
        ev.ndw_L, ev.ndw_R,
        ev.detL, ev.detR, ev.detJ,
        co.okA[SJOINT], co.okA[SRUBIN], co.okA[SROMAN],
        co.Era[SJOINT][1], co.Era[SRUBIN][1], co.Era[SROMAN][1],   //sigma(tE)
        co.Era[SJOINT][3], co.Era[SRUBIN][3], co.Era[SROMAN][3],   //sigma(piE)
        co.Erb[SJOINT][0], co.Erb[SRUBIN][0], co.Erb[SROMAN][0],   //sigma(tetE)
        ev.dclsEvent,
        synergyClass(co),
        co.condA[SJOINT], co.condA[SRUBIN], co.condA[SROMAN],
        // ---- the rest of the row (Step D1) ----
        l.t0, s.xi, s.lon, s.lat,
        s.mbs[1], s.fb[1],
        {s.magb[0], s.magb[1], s.magb[2], s.magb[3], s.magb[4], s.magb[5], s.magb[6]},
        {s.blend[0], s.blend[1], s.blend[2], s.blend[3], s.blend[4], s.blend[5], s.blend[6]},
        co.relMl[SJOINT], co.relMl[SRUBIN], co.relMl[SROMAN],
        co.okB[SJOINT], co.okB[SRUBIN], co.okB[SROMAN],
        co.condB[SJOINT], co.condB[SRUBIN], co.condB[SROMAN],
        sched.dtToSeasonEdge(ev.t0obs), sched.zone(ev.t0obs),
        ev.nres5_L, ev.nres20_L, ev.nresPSF_L, ev.nres5_R, ev.nres20_R, ev.nresPSF_R,
        ev.dsepMax_L, ev.dsepMax_R
    });

    // ------------------------------------------------------------------------------
    // Step H2. Two quantities that make the satellite-parallax effect visible in the
    // table instead of only implicit in the Fisher matrix.
    //
    // du_sat  the observer separation in Einstein radii at t0, piE * D_perp / AU.
    //         This is the amplitude of the satellite effect for this event and the
    //         natural x-axis of every Step H3 figure. It is NOT L2_OFFSET_AU * piE:
    //         D_perp is the separation projected perpendicular to the line of sight
    //         and runs ~0.87-0.99 of the full L2 offset around the year, so the
    //         simple product is a ceiling (DEVIATIONS.md 28.2). Computed by asking
    //         lightcurve() for both observers rather than re-deriving the projection
    //         here, so the two can never drift apart.
    //
    // nepL_pk, nepR_pk  epochs from each survey within +-2 tE of t0, i.e. while the
    //         event is actually magnified. Satellite parallax needs CONTEMPORANEOUS
    //         coverage: an event Roman saw in season 3 and Rubin saw in season 7 has
    //         none of it, however large ndw_L and ndw_R are. The plan asked for a
    //         single `nep_both` flag; two counts are the same cost and strictly more
    //         informative, and the flag is just (nepL_pk > 0 and nepR_pk > 0).
    //
    // Safe to call lightcurve() here: FisherM has already run, and the only state it
    // touches (s->ux/uy, s->def*, s->pos*, as->ue_n*) is not read by the row written
    // below and is recomputed from scratch by the next draw.
    // ------------------------------------------------------------------------------
    ev.duSat = 0.0;
    {
        lightcurve(s, l, as, l.t0, 0);
        const double r1 = as.ue_n1, r2 = as.ue_n2;
        lightcurve(s, l, as, l.t0, 1);
        const double dn1 = as.ue_n1 - r1, dn2 = as.ue_n2 - r2;
        ev.duSat = l.piE * std::sqrt(dn1 * dn1 + dn2 * dn2);
    }
    ev.nepLpk = 0, ev.nepRpk = 0;
    {
        const double win = 2.0 * l.tE;
        for (int i = 0; i < ev.ndw; ++i) {
            if (std::fabs(l.timn[i] - ev.t0obs) > win) continue;   // Deviation 76
            if (int(l.tele[i]) == 1) ev.nepRpk += 1;
            else                      ev.nepLpk += 1;
        }
    }

    filg_in.open(testf, std::ios::app);
    filg_in << st.icon           << " " << ev.FFG[0]         << " " << l.tE                      << " "
            << l.RE / AU     << " " << l.piE         << " " << l.tetE                    << " "
            << l.Vt          << " " << l.u0          << " " << l.Ml                      << " " 
            << s.opt * 1.0e6 << " " << l.Dl          << " " << s.Ds                      << " "
            << l.vl          << " " << s.vs          << " " << s.mbs[0]                  << " "
            << s.fb[0]       << " " << ev.gg             << " " << static_cast<int>(l.struc) << " "
            << s.FWHM / year << " " << ev.vsave / s.mus << " " << l.DeltaT / s.errA        << " " << l.murel * year << " "
            << co.resu[0]    << " " << co.resu[1]    << " " << co.resu[2]                << " " 
            << co.resu[3]    << " " << co.resu[5]    << " "
            << co.resu[9]    << " " << co.resu[10]   << " " << co.resu[13]               << " " << co.resu[14]    << " "
            << s.Map[2]      << " " << s.nsbl[2]     << " " << co.flagi                  << " " << s.Ai[2]        << " "
            // per-survey bookkeeping (Step C5), appended so existing column indices hold
            << ev.ndw_L << " " << ev.ndw_R << " "
            << ev.detL  << " " << ev.detR  << " " << ev.detJ << " "
            << co.okA[SJOINT] << " " << co.okA[SRUBIN] << " " << co.okA[SROMAN] << " "
            << co.Era[SJOINT][1] << " " << co.Era[SRUBIN][1] << " " << co.Era[SROMAN][1] << " "
            << co.Era[SJOINT][3] << " " << co.Era[SRUBIN][3] << " " << co.Era[SROMAN][3] << " "
            << co.Erb[SJOINT][0] << " " << co.Erb[SRUBIN][0] << " " << co.Erb[SROMAN][0] << " "
            << ev.dclsEvent << " " << synergyClass(co) << " "
            << co.condA[SJOINT] << " " << co.condA[SRUBIN] << " " << co.condA[SROMAN] << " "
            // the rest of the row (Step D1) -- appended, so columns 1-57 keep their indices
            << l.t0 << " " << s.xi << " " << s.lon << " " << s.lat << " "
            << s.mbs[1] << " " << s.fb[1] << " ";
    for (int i = 0; i < M; ++i) filg_in << s.magb[i]  << " ";
    for (int i = 0; i < M; ++i) filg_in << s.blend[i] << " ";
    filg_in << co.relMl[SJOINT] << " " << co.relMl[SRUBIN] << " " << co.relMl[SROMAN] << " "
            << co.okB[SJOINT]   << " " << co.okB[SRUBIN]   << " " << co.okB[SROMAN]   << " "
            << co.condB[SJOINT] << " " << co.condB[SRUBIN] << " " << co.condB[SROMAN] << " "
            // Gap geometry. dt_edge is NEGATIVE when t0 fell inside a Roman season;
            // t0zone distinguishes a mid-mission gap (1) from before-launch/after-end
            // (2), which must never be pooled -- only the former is gap-filling.
            << sched.dtToSeasonEdge(ev.t0obs) << " " << sched.zone(ev.t0obs) << " "
            // Sky area this event's sightline stands for, deg^2 (Step E1). Constant
            // across an unstratified run; NOT constant once --stride-roman is used,
            // and then any statistic pooled over sightlines must weight by it.
            << st.wArea << " "
            // Step H2: the satellite-parallax observable and contemporaneous coverage.
            << ev.duSat << " " << ev.nepLpk << " " << ev.nepRpk << " "
            // Step R1: resolving the two images. Counts of qualifying epochs per
            // survey, then the largest separation reached while both were detectable
            // (-1 = never). The three bars differ only in what counts as "resolved".
            << ev.nres5_L << " " << ev.nres20_L << " " << ev.nresPSF_L << " " << ev.dsepMax_L << " "
            << ev.nres5_R << " " << ev.nres20_R << " " << ev.nresPSF_R << " " << ev.dsepMax_R << " "
            // Deviation 71: astrometric noise variants N and P (joint, Roman).
            << co.ErbV[AV_N][SJOINT][0] << " " << co.ErbV[AV_N][SROMAN][0] << " "
            << co.ErbV[AV_P][SJOINT][0] << " " << co.ErbV[AV_P][SROMAN][0] << " "
            << co.relMlV[AV_N][SJOINT] << " " << co.relMlV[AV_N][SROMAN] << " "
            << co.relMlV[AV_P][SJOINT] << " " << co.relMlV[AV_P][SROMAN] << " "
            << co.okBV[AV_N][SJOINT] << " " << co.okBV[AV_N][SROMAN] << " "
            << co.okBV[AV_P][SJOINT] << " " << co.okBV[AV_P][SROMAN] << " "
            << int(l.luminous) << " " << s.fLens[0] << " " << s.fLens[1] << " "
            << ev.t0obs << " " << ev.uminObs << "\n";
    filg_in.close();

    commitSampleDump(ctx, ev);
    writeSatellitePair(ctx, st, ev);
}

void commitSampleDump(SimContext& ctx, const EventState& ev) {
    const RomanSchedule& sched = ctx.sched;
    source& s = ctx.s;
    lens& l = ctx.l;
    astromet& as = ctx.as;
    covarian& co = ctx.co;
    SampleSpec& dumpSpec = ctx.outs.dumpSpec;
    std::vector<DumpEpoch>& dumpBuf = ctx.outs.dumpBuf;
    long& dumpSeq = ctx.outs.dumpSeq;

    // ------------------------------------------------------------------------------
    // Step S1. Commit this draw's buffered light curve if it fills a requested
    // sample class. Here, and not in the time loop, because detL/detR/detJ and the
    // Fisher sigmas -- which is what the classes are defined in terms of -- do not
    // exist until now.
    //
    // First match wins: an event is written once, under the first class in the spec
    // file whose quota is still open. Otherwise a `both` event would also land in
    // `any` and the same light curve would be drawn twice in one figure.
    // ------------------------------------------------------------------------------
    if (dumpSpec.on and not dumpBuf.empty()) {
        double maxShift = 0.0;
        for (const DumpEpoch& e : dumpBuf)
            maxShift = std::max(maxShift,
                                std::sqrt(e.def1c * e.def1c + e.def2c * e.def2c));

        const SampleFacts facts{ev.detL, ev.detR, ev.detJ, l.tE, sched.zone(ev.t0obs),
                                ev.nepLpk, ev.nepRpk, ev.ndw_L, ev.ndw_R,
                                co.okB[SROMAN], maxShift};

        for (SampleClass& c : dumpSpec.classes) {
            if (c.kept >= c.quota)        continue;
            if (!sampleMatch(c, facts))   continue;

            std::ostringstream pr;
            pr << std::setprecision(10)
               << "# Sample event for class '" << c.name << "'.\n"
               << "# One key per line. Angles in deg, times in days, masses in Msun,\n"
               << "# distances in kpc, angular scales in mas. A sigma of -1 means the\n"
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
               << "ndw_L "        << ev.ndw_L         << "\n"
               << "ndw_R "        << ev.ndw_R         << "\n"
               << "nep_pk_L "     << ev.nepLpk        << "\n"
               << "nep_pk_R "     << ev.nepRpk        << "\n"
               << "detL "         << ev.detL          << "\n"
               << "detR "         << ev.detR          << "\n"
               << "detJ "         << ev.detJ          << "\n"
               << "dt_edge "      << sched.dtToSeasonEdge(ev.t0obs) << "\n"
               << "t0zone "       << sched.zone(ev.t0obs)           << "\n"
               << "du_sat "       << ev.duSat         << "\n"
               << "max_shift "    << maxShift      << "\n"
               << "dsep_max_L "   << ev.dsepMax_L     << "\n"
               << "dsep_max_R "   << ev.dsepMax_R     << "\n"
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
                      << "  ndw_L=" << ev.ndw_L << " ndw_R=" << ev.ndw_R
                      << "  (" << c.kept << "/" << c.quota << ")" << std::endl;
            break;
        }
    }
}

