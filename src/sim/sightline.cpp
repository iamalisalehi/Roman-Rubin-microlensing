// Per-sightline stages of the Monte Carlo: setup before the draws, aggregation after them.
#include "sim/sightline.h"
#include "util/random.h"
#include "galaxy/density.h"
#include "galaxy/extinction.h"
#include "surveys/schedule.h"
#include "surveys/footprints.h"
#include "events/lightcurve.h"
#include "fisher/fisher.h"

SightlineStart setupSightline(SimContext& ctx, SightlineState& st, const Sightline& sightline) {
    const RunConfig& cfg = ctx.cfg;
    const GbtdsLayout& gl = ctx.gl;
    source& s = ctx.s;
    lens& l = ctx.l;
    extin& ex = ctx.ex;
    lsst& ls = ctx.ls;
    roman& ro = ctx.ro;
    RunTotals& run = ctx.run;

    st.iScan += 1;
    s.lon = sightline.lon;
    s.lat = sightline.lat;
    // The sky area THIS sightline stands for. Written into every event row it produces;
    // any statistic pooled across sightlines has to weight by it (see the list build).
    st.wArea = sightline.area;
    // nri/nde keep their old meaning -- longitude column index, and index within that
    // column -- which is what the map file and the event table record. Assigned rather
    // than incremented: the old loop bumped nri once per iLon whether or not any
    // sightline in that column survived the corner cut, so an incrementing counter here
    // would silently renumber the columns the moment a fully-cut column existed. Under
    // stratification the column index is a fine-grid one, so it steps by kSub between
    // coarse columns.
    st.nri = sightline.col;
    if (sightline.col != st.lastCol) { st.nde = -1; st.lastCol = sightline.col; }
    st.nde += 1;

    // --start-index. Resuming an interrupted run. The nri/nde bookkeeping above runs for
    // skipped sightlines too, deliberately: nde counts position within a longitude column
    // and is written into the map file, so skipping it would renumber the first partial
    // column and a resumed run would not concatenate onto the interrupted one. Only the
    // simulation is skipped, not the numbering.
    //
    // The index to resume AT is the number of sightlines entered, which is the number of
    // "NEW STEP" lines in the interrupted run's log -- this print sits above the
    // no-coverage and barren skips, so it counts every sightline the scan reached. The map
    // file counts only sightlines that produced something to aggregate and is therefore a
    // smaller number; resuming at it would redo the difference and duplicate rows.
    //
    // This is safe precisely because `scan` is built deterministically from --stride,
    // --stride-roman and --stub and does not depend on the RNG. Note that the random
    // sequence IS advanced by the sightlines a full run would have simulated, so a
    // resumed run is not bit-identical to an uninterrupted one -- it is a valid
    // continuation with a different draw sequence, not a replay.
    if (st.iScan < cfg.startIndex) return SightlineStart::Skip;
    if (cfg.endIndex >= 0 and st.iScan >= cfg.endIndex) return SightlineStart::Stop;   // Deviation 77
    // Deviation 77: this sightline's own random stream (see sightlineSeed in include/util/random.h). The note
    // above, that a resumed run is "not a replay", no longer applies: it IS one.
    rng.seed(sightlineSeed(cfg.seedBase, st.iScan));
    cout << ">>>>>>>>>>> NEW STEP " << st.nde << " <<<<<<<<\t nri:  " << st.nri << endl;
    cout << "longtitude: " << s.lon << "\t latitude: " << s.lat << endl;


    // Deviation 80: on LSSTCam's active silicon for this visit's pointing and rotation
    // (was: within a 1.75-deg circle). A cheap (l, b) distance cut first.
    double slRA, slDec;
    galToIcrs(s.lon, s.lat, slRA, slDec);
    auto rubinCovers = [&](int i) {
        const double dl = s.lon - ls.l[i], db = s.lat - ls.b[i];
        if (dl * dl + db * db > (RUBIN_MAX_RADIUS + 0.1) * (RUBIN_MAX_RADIUS + 0.1)) return false;
        return onLsstCam(slRA, slDec, ls.RA[i], ls.DEC[i], ls.rot[i]);
    };
    auto romanCovers = [&](int i) {
        return inDetector(gl, ro.layout[i], s.lon - ro.l[i], s.lat - ro.b[i]);
    };
    st.ndd  = matchVisibleEpochs("LSST",  rubinCovers, ls.tim, Nl,      ls.ct, st.minc);
    st.nddR = matchVisibleEpochs("Roman", romanCovers, ro.tim, NlRoman, ro.ct, st.mincR);

    // Rubin's depth per band at this sightline (Deviation 73): the median of its matched
    // visits' own 5-sigma depths, for the pre-selection below. A band with no visit here
    // gets -inf, so it cannot count toward "detectable in >= 2 bands".
    {
        std::array<std::vector<double>, 6> d5;
        for (int k = 0; k < st.ndd; ++k) {
            const int v = int(ls.ct[k]);
            d5[ls.filter[v]].push_back(ls.sig5[v]);
        }
        for (int b = 0; b < 6; ++b) {
            if (d5[b].empty()) { st.rubinDepthMed[b] = -std::numeric_limits<double>::infinity(); continue; }
            std::nth_element(d5[b].begin(), d5[b].begin() + d5[b].size() / 2, d5[b].end());
            st.rubinDepthMed[b] = d5[b][d5[b].size() / 2];
        }
    }
    cout << "ndd (LSST): "  << st.ndd  << "\t minc (LSST): "  << st.minc  << endl;
    cout << "ndd (Roman): " << st.nddR << "\t minc (Roman): " << st.mincR << endl;

    // Empty sky. Neither survey visits this sightline, so no light curve can ever
    // have a datum, no event can be detected, and nerr can never advance -- the
    // do/while below would spin forever. Skipping is not an approximation: a
    // sightline with no epochs contributes exactly zero events to the yield.
    //
    // It has to be a `continue` rather than a run that finds nothing, because the
    // aggregation block at the end of this loop asserts CHECK(numd[1] != 0.0) and
    // CHECK(nerr != 0.0) -- reaching it with an empty field aborts the whole run.
    //
    // The skipped area is NOT written to the map file, so anything converting the
    // per-sightline Neven density into a total count must use the aggregated
    // sightline count reported at the end of the run, not the grid size.
    if (st.ndd == 0 and st.nddR == 0) {
        run.nSkipNoCoverage += 1;
        run.areaNoCoverage  += st.wArea;
        cout << "  SKIP: no Rubin and no Roman epochs at this sightline" << endl;
        return SightlineStart::Skip;
    }

    st.icon  = 0;
    st.nlens = 0;
    st.nDetClass.fill(0);
    st.nsim  = 0.0;
    st.nerr  = 0.0;
    for (int i = 0; i < Num; ++i) { s.nssim[i] = 0.0;  s.nsdet[i] = 0.0; }
    for (int i = 0; i <= GG; ++i) { l.nstE[i]  = 0.0;  l.ndtE[i]  = 0.0; }

    s.TET = (360.0 - s.lon) / RAa;///radian s.lon/RA;//
    s.FI  = s.lat / RAa;

    Disk_model(s, 1);
    st.sightlineIdx = nearestSightline(ex, s.lon, s.lat);

    st.records.clear();
    return SightlineStart::Simulate;
}

void finishSightline(SimContext& ctx, SightlineState& st) {
    const RunConfig& cfg = ctx.cfg;
    source& s = ctx.s;
    lens& l = ctx.l;
    covarian& co = ctx.co;
    std::ofstream& fil2 = ctx.outs.fil2;
    std::ofstream& fil2b = ctx.outs.fil2b;
    std::ofstream& fil3 = ctx.outs.fil3;
    RunTotals& run = ctx.run;

    // Per-sightline aggregates, filled from the replayed records below.
    std::array<double, 2> tE, RE, piE, tetE, Vt, u0, Ml, opd, Dl, Ds, vl, vs;
    std::array<double, 2> numd, Struc, murel, vsn, DelT, mbs, fb, fwhm, Map, nbl, Ext;
    double test, mbase, fblend, vsave, Gamma, Neven, EFF, EffiD, EffiL;
    double ErtE, ErpiE, ErtetE, Erml, Erdl, Ermul, Ermus, Eru0, Erfb, nErAvg;
    int    flagL, gg = -1;

    // Did the sightline actually meet its budget, or did the cap stop it? The
    // distinction matters: a capped sightline's Poisson precision is whatever it
    // reached, not what was asked for, and averaging it in as though it were a
    // full sample would understate the error bars.
    const bool budgetMet = (st.icon  >= cfg.iconTarget and
                            st.nlens >= cfg.nlensTarget and
                            st.nerr  >= cfg.nerrTarget);
    if (!budgetMet) {
        run.nCapped += 1;
        cout << "  CAP: stopped at nsim = " << st.nsim << " with icon = " << st.icon
             << "/" << cfg.iconTarget << ", nlens = " << st.nlens << "/" << cfg.nlensTarget
             << ", nerr = " << st.nerr << "/" << cfg.nerrTarget << endl;
    }

    // Some epochs existed but nothing survived to be aggregated. Same reasoning as
    // the no-coverage skip above: the CHECK block below requires at least one
    // detected AND one characterised event, so this must not fall through.
    if (st.nlens < 1 or st.nerr <= 0.0) {
        run.nSkipBarren += 1;
        run.areaBarren  += st.wArea;
        cout << "  BARREN: " << st.nlens << " detected, " << st.nerr
             << " characterised -- nothing to aggregate at this sightline" << endl;
        return;
    }
    run.nAggregated    += 1;
    run.areaAggregated += st.wArea;

    for (int i = 0; i <= GG; ++i) {
        l.NstE[i] += l.nstE[i];
        l.NdtE[i] += l.ndtE[i];
        l.ndtE[i]  = double(l.ndtE[i] / (l.nstE[i] + eps)); // [0.0, 1.0]

        fil2 << std::fixed << std::setprecision(4)
             << l.tEs[i]  << " "
             << double(l.NdtE[i]  * 100.0 / (l.NstE[i] + eps)) << " "
             << l.Mls[i]  << " "
             << double(l.NdMl[i]  * 100.0 / (l.NsMl[i] + eps)) << " "
             << l.pis[i]  << " "
             << double(l.Ndpi[i]  * 100.0 / (l.Nspi[i] + eps)) << " "
             << l.u0s[i]  << " "
             << double(l.Ndu0[i]  * 100.0 / (l.Nsu0[i] + eps)) << " "
             << l.mbs[i]  << " "
             << double(l.Ndmb[i]  * 100.0 / (l.Nsmb[i] + eps)) << " "
             << l.fbs[i]  << " "
             << double(l.Ndfb[i]  * 100.0 / (l.Nsfb[i] + eps)) << " "
             << l.mus[i]  << " "
             << double(l.Ndmu[i]  * 100.0 / (l.Nsmu[i] + eps)) << " "
             // The + eps every other column has, and these two lacked: before the
             // counters above existed this was 0/0 and printed "-nan" on every row.
             << double(l.Nhalo[1] * 100.0 / (l.Nhalo[0] + eps)) << " "
             << double(l.Nself[1] * 100.0 / (l.Nself[0] + eps)) << "\n";

        fil2b << std::fixed  << std::setprecision(1)
              << l.NdtE[i]  << " " << l.NstE[i]  << " "
              << l.NdMl[i]  << " " << l.NsMl[i]  << " "
              << l.Ndpi[i]  << " " << l.Nspi[i]  << " "
              << l.Ndu0[i]  << " " << l.Nsu0[i]  << " "
              << l.Ndmb[i]  << " " << l.Nsmb[i]  << " "
              << l.Ndfb[i]  << " " << l.Nsfb[i]  << " "
              << l.Ndmu[i]  << " " << l.Nsmu[i]  << " "
              << l.Nhalo[1] << " " << l.Nhalo[0] << " "
              << l.Nself[1] << " " << l.Nself[0] << "\n";
    }

    s.nstart = 0.0;
    for (int i = 0; i < Num; ++i) {
        test = double(s.nsdet[i]   / (s.nssim[i] + eps));
        s.nstart  += s.Rostari[i] * (s.Nstart   / s.Rostart) * test;
    } //Number/deg^{2}

    for (int i = 0; i < 2; ++i){
        tE[i]   = 0.0, RE[i]  = 0.0, piE[i]  = 0.0, tetE[i]  = 0.0; Vt[i]    = 0.0;
        u0[i]   = 0.0, Ml[i]  = 0.0, opd[i]  = 0.0; Dl[i]    = 0.0, Ds[i]    = 0.0;
        vl[i]   = 0.0, vs[i]  = 0.0, mbs[i]  = 0.0, fb[i]    = 0.0; numd[i]  = 0.0;
        fwhm[i] = 0.0; vsn[i] = 0.0; DelT[i] = 0.0; Struc[i] = 0.0; murel[i] = 0.0;
        Map[i]  = 0.0; nbl[i] = 0.0; Ext[i]  = 0.0;
    }

    EffiL = 0.0; EFF   = 0.0; Gamma  = 0.0; Neven = 0.0; EffiD = 0.0;
    ErtE  = 0.0; ErpiE = 0.0; ErtetE = 0.0; Erml  = 0.0; Erfb  = 0.0; nErAvg = 0.0;
    Erdl  = 0.0; Ermul = 0.0; Ermus  = 0.0; Eru0  = 0.0;

    int tempStruc;
    for (const auto& r : st.records) {
//        counter = r.counter;
        flagL = r.flagL;
        l.tE = r.tE;   l.RE = r.RE;   l.piE = r.piE;  l.tetE = r.tetE;
        l.Vt = r.Vt;   l.u0 = r.u0;   l.Ml  = r.Ml;
        s.opt = r.opt; l.Dl = r.Dl;   s.Ds  = r.Ds;   l.vl = r.vl; s.vs = r.vs;
        mbase = r.mbase; fblend = r.fblend; gg = r.gg; tempStruc = r.struc;
        s.FWHM = r.FWHM; vsave = r.vsave; l.DeltaT = r.DeltaT; l.murel = r.murel;
        co.resu[0]=r.resu0;  co.resu[1]=r.resu1;   co.resu[2]=r.resu2;
        co.resu[3]=r.resu3;  co.resu[5]=r.resu5;   co.resu[9]=r.resu9;
        co.resu[10]=r.resu10; co.resu[13]=r.resu13; co.resu[14]=r.resu14;
        s.Map[2]=r.Map2; s.nsbl[2]=r.nsbl2; co.flagi=r.flagi; s.Ai[2]=r.Ai2;
        // Must be replayed too: this loop runs after the whole field, so co->okA
        // otherwise holds whatever the LAST FisherM call left, not this event's.
        co.okA[SJOINT] = r.okJoint;

        l.struc = static_cast<GalacticComponent>(tempStruc);

        // ------------------ Update cumulative sums ------------------
        Struc[0] += 1.0;
        tE[0]    += l.tE;    RE[0]  += l.RE;     piE[0]  += l.piE;
        tetE[0]  += l.tetE;  Vt[0]  += l.Vt;     u0[0]   += l.u0;
        Ml[0]    += l.Ml;    opd[0] += s.opt;    Dl[0]   += l.Dl;
        Ds[0]    += s.Ds;    vl[0]  += l.vl;     vs[0]   += s.vs;
        mbs[0]   += mbase;    fb[0]  += fblend;    numd[0] += 1.0;
        fwhm[0]  += s.FWHM;  vsn[0] += vsave;     DelT[0] += l.DeltaT;
        murel[0] += l.murel; Map[0] += s.Map[2]; nbl[0]  += s.nsbl[2];
        Ext[0]   += s.Ai[2];

        // ------------------ Flagged case ---------------------------
        if (flagL > 0) {
            Struc[1] += 1.0;
            tE[1]    += l.tE;    RE[1]  += l.RE;     piE[1]  += l.piE;
            tetE[1]  += l.tetE;  Vt[1]  += l.Vt;     u0[1]   += l.u0;
            Ml[1]    += l.Ml;    opd[1] += s.opt;    Dl[1]   += l.Dl;
            Ds[1]    += s.Ds;    vl[1]  += l.vl;     vs[1]   += s.vs;
            mbs[1]   += mbase;    fb[1]  += fblend;    numd[1] += 1.0;
            fwhm[1]  += s.FWHM;  vsn[1] += vsave;     DelT[1] += l.DeltaT;
            murel[1] += l.murel; Map[1] += s.Map[2]; nbl[1]  += s.nsbl[2];
            Ext[1]   += s.Ai[2];

            EFF += static_cast<double>(l.ndtE[gg] / (l.tE / year)); // 1/years

            // Average only over events the joint fit could actually characterize. Step C4
            // reports sigma = -1 as an explicit "not characterizable" sentinel, and ErrorCal
            // divides it by the parameter value, so an unguarded sum pulls in large negative
            // numbers: one event in this field contributed resu[3] = -697, dragging the mean
            // fractional piE error negative and tripping CHECK(ErpiE > 0.0).
            //
            // Skipping those events is not the selection bias DEVIATIONS entry 8 warns about.
            // That rule is about never dropping an event from the joint-vs-single RATIO
            // statistics, where a missing single-survey sigma is itself the result. Here we are
            // forming a mean precision, and an event with no measurement has no precision to
            // average -- including it would be averaging a sentinel. The count of events the
            // mean is actually over is tracked separately so the denominator is honest.
            //
            // okA[SJOINT] is NOT sufficient on its own. It says the joint photometric
            // matrix inverted -- not that every parameter was in the fit. Since the joint
            // refactor gave each survey partition its own active parameter subset
            // (activePhotParams in include/fisher/fisher.h), fb0 and mbs0 enter the joint fit only when the
            // event has Rubin epochs, and fb1/mbs1 only when it has Roman ones. An event
            // detected by Roman with no Rubin data therefore has a perfectly valid joint
            // fit in which Era[2] is still the -1.0 sentinel, and ErrorCal divides that by
            // fb0 regardless: one such event contributed resu[2] = -8499 at l=0.881,
            // b=-0.94 and dragged a 102-event mean to -83, tripping CHECK(Erfb > 0.0).
            //
            // So test the values themselves. An event is averaged only if all nine are real
            // measurements, which keeps every mean over the same event set and keeps
            // nErAvg meaningful as a single denominator. The cost is small and measured:
            // across 48,959 characterised events in the 2026-08-29 run exactly ONE was
            // excluded by this, and no column other than resu[2] was ever negative.
            const bool allMeasured = (co.flagi > 0 and co.okA[SJOINT]
                                      and co.resu[0]  >= 0.0 and co.resu[1]  >= 0.0
                                      and co.resu[2]  >= 0.0 and co.resu[3]  >= 0.0
                                      and co.resu[5]  >= 0.0 and co.resu[9]  >= 0.0
                                      and co.resu[10] >= 0.0 and co.resu[13] >= 0.0
                                      and co.resu[14] >= 0.0);
            if (allMeasured) {
                Eru0  += co.resu[0];  ErtE   += co.resu[1];  Erfb  += co.resu[2];
                ErpiE += co.resu[3];  ErtetE += co.resu[5];  Erml  += co.resu[9];
                Erdl  += co.resu[10]; Ermul  += co.resu[13]; Ermus += co.resu[14];
                nErAvg += 1.0;
            }
        }
    }

    for (int i = 0; i < 2; ++i) { // What is the purpose of this block?
        tE[i]    = double(tE[i]            / (numd[i] + eps)) / year; //[years]  
        RE[i]    = double(RE[i]            / (numd[i] + eps));//[AU]  
        piE[i]   = double(piE[i]           / (numd[i] + eps));//[]     
        tetE[i]  = double(tetE[i]          / (numd[i] + eps));//[mas]  
        Vt[i]    = double(Vt[i]            / (numd[i] + eps));//[km/s]  
        u0[i]    = double(u0[i]            / (numd[i] + eps));//[]  
        Ml[i]    = double(Ml[i]            / (numd[i] + eps));//[Msun]  
        opd[i]   = double(opd[i]           / (numd[i] + eps)) * u0m * u0m;//[] x10^{6}
        Dl[i]    = double(Dl[i]            / (numd[i] + eps));//[kpc]  
        Ds[i]    = double(Ds[i]            / (numd[i] + eps));//[kpc]  
        vl[i]    = double(vl[i]            / (numd[i] + eps));//[km/s]  
        vs[i]    = double(vs[i]            / (numd[i] + eps));//[km/s]  
        mbs[i]   = double(mbs[i]           / (numd[i] + eps));//[mag]  
        fb[i]    = double(fb[i]            / (numd[i] + eps));//[]
        fwhm[i]  = double(fwhm[i]          / (numd[i] + eps));//[years]
        vsn[i]   = double(vsn[i]           / (numd[i] + eps));//[]
        DelT[i]  = double(DelT[i]          / (numd[i] + eps));//[]
        murel[i] = double(murel[i]         / (numd[i] + eps));//[mas/year]
        Map[i]   = double(Map[i]           / (numd[i] + eps));//[mag]
        nbl[i]   = double(nbl[i]           / (numd[i] + eps));//[]
        Ext[i]   = double(Ext[i]           / (numd[i] + eps));//[mag]
        Struc[i] = double(Struc[i] * 100.0 / (numd[i] + eps));
    }

    EFF = double(EFF / (numd[1] + eps));//1/[years]
    Gamma = double(2.0 / M_PI * opd[0] * 1.0e-6 * EFF) / u0m; // 1/[star*year]  
    EffiL = double(numd[1] * 100.0 / (numd[0] + eps)); // probability of detecting lensing  
    EffiD = double(st.icon * 100.0 / (st.nsim    + eps)); // % of drawn stars that are visible (Deviation 78; was numd[0], i.e. 100% by construction)
    Neven = double(s.nstart * Gamma * 10.0);//deg^{-2}
   
    Eru0   = double(Eru0   / (nErAvg + eps));  
    Erfb   = double(Erfb   / (nErAvg + eps));  
    ErtE   = double(ErtE   / (nErAvg + eps));  
    ErpiE  = double(ErpiE  / (nErAvg + eps));  
    ErtetE = double(ErtetE / (nErAvg + eps));  
    Erml   = double(Erml   / (nErAvg + eps));  
    Erdl   = double(Erdl   / (nErAvg + eps));  
    Ermul  = double(Ermul  / (nErAvg + eps));  
    Ermus  = double(Ermus  / (nErAvg + eps));  
   
    //Maps
    fil3 << std::fixed << std::setprecision(6)
         << tE[0]    << " " << tE[1]    << " "
         << RE[0]    << " " << RE[1]    << " "
         << piE[0]   << " " << piE[1]   << " "
         << tetE[0]  << " " << tetE[1]  << " "
         << Vt[0]    << " " << Vt[1]    << " "
         << u0[0]    << " " << u0[1]    << " "
         << Ml[0]    << " " << Ml[1]    << " "
         << opd[0]   << " " << opd[1]   << " "
         << Dl[0]    << " " << Dl[1]    << " "
         << Ds[0]    << " " << Ds[1]    << " "
         << vl[0]    << " " << vl[1]    << " "
         << vs[0]    << " " << vs[1]    << " "
         << mbs[0]   << " " << mbs[1]   << " "
         << fb[0]    << " " << fb[1]    << " "
         << fwhm[0]  << " " << fwhm[1]  << " "
         << vsn[0]   << " " << vsn[1]   << " "
         << DelT[0]  << " " << DelT[1]  << " "
         << Struc[0] << " " << Struc[1] << " "
         << murel[0] << " " << murel[1] << " "
         << Map[0]   << " " << Map[1]   << " "
         << nbl[0]   << " " << nbl[1]   << " "
         << Ext[0]   << " " << Ext[1]   << " ";

    fil3 << std::setprecision(8)
         << EffiD      << " " << EffiL        << " "
         << std::log10(EFF) << " " << std::log10(Gamma) << " " << std::log10(Neven) << " "
         << Eru0       << " " << ErtE         << " " << Erfb         << " " << ErpiE << " " << ErtetE << " "
         << Erml       << " " << Erdl         << " " << Ermul        << " " << Ermus << " ";

    fil3 << std::setprecision(1)
         << st.nsim              << " " << numd[0]          << " " << numd[1] << " "
         << st.nerr              << " " << st.nri              << " " << st.nde     << " "
         << std::log10(s.Rostart) << " " << std::log10(s.Nstart) << " " << std::log10(s.nstart)
         // Step E1. Three columns appended, in this order:
         //   w_area  deg^2 of sky this sightline stands for -- no longer a run-wide constant
         //   lon,lat where it is. The map file had NO position column at all, so a row in it
         //           could not be tied to the events it produced, and the draw count `nsim`
         //           it records -- the denominator any pooled yield needs -- was unreachable
         //           from the event table. With these, an event joins its sightline on
         //           (lon, lat) and the correct pooled weight, w_area/nsim, is computable.
         << " " << std::setprecision(8) << st.wArea
         << " " << std::setprecision(6) << s.lon << " " << s.lat
         << "\n";

    // Flush the per-sightline outputs now rather than when the stream is destroyed. Every
    // production pause so far has been a kill, and a kill discards whatever is still
    // buffered: the 2026-09-06 chunk-1 stop lost the map rows of six completed sightlines
    // and left a half-written seventh, onto which the resuming run's first row was then
    // appended -- one 122-field line that made the whole file unreadable until the Python
    // reader learned to skip it (Deviation 41). The same stop cost EfLMC5/EfLMC5B the same
    // six blocks, which is why fil2/fil2b are flushed here too. A sightline costs minutes
    // of CPU, so three flushes per sightline are free, and what reaches disk is then what
    // the log says was finished.
    fil3.flush();
    fil2.flush();
    fil2b.flush();

    cout << "nsim:  "  << st.nsim    << "\t Ndetected:  " << st.icon    << "\t Nlensing:  " << st.nlens << "\t NError:  " << st.nerr << endl;
    cout << "Detection classes:";
    for (int c = 0; c < NDETCLASS; ++c)
        cout << "  " << detClassName(c) << ": " << st.nDetClass[c];
    cout << endl;
    if (st.nDetClass[DET_ANOMALY] > 0) {
        cout << "  WARNING: " << st.nDetClass[DET_ANOMALY] << " event(s) detected by one telescope "
             << "but NOT by the joint test. Adding data cannot destroy signal, so this is a "
             << "threshold inconsistency -- see DetClass in include/fisher/fisher.h." << endl;
    }
    cout << "numd0:  " << numd[0] << "\t numd1:  "     << numd[1] << endl;
    cout << "EFF:  "   << EFF     << "\t Gamma:  "     << Gamma   << "\t Neven:  "    << Neven << endl;
    cout << "l.tE:  "  << l.tE   << "\t l.Ml:  "      << l.Ml   << endl;
   
    CHECK(EFF > 0.0);
    CHECK(Gamma > 0.0);
    CHECK(EffiL > 0.0);
    CHECK(EffiD > 0.0);
    CHECK(Neven > 0.0);
    
    CHECK(DelT[0] > 0.0);
    CHECK(DelT[1] > 0.0);
    CHECK(vsn[0] > 0.0);
    CHECK(vsn[1] > 0.0);
    CHECK(fwhm[0] > 0.0);
    CHECK(fwhm[1] > 0.0);
    
    CHECK(fb[0] > 0.0);
    CHECK(fb[0] <= 1.0);
    CHECK(fb[1] > 0.0);
    CHECK(fb[1] <= 1.0);
    
    CHECK(mbs[0] > 0.0);
    CHECK(mbs[1] > 0.0);
    
    CHECK(vs[0] > 0.0);
    CHECK(vs[1] > 0.0);
    CHECK(vl[0] > 0.0);
    CHECK(vl[1] > 0.0);
    
    CHECK(Ds[0] > 0.0);
    CHECK(Ds[1] > 0.0);
    CHECK(Dl[0] > 0.0);
    CHECK(Dl[1] >= 0.0);
    
    CHECK(piE[1] >= 0.0);
    CHECK(tetE[1] > 0.0);
    
    CHECK(murel[0] > 0.0);
    CHECK(murel[1] > 0.0);
    
    if (nErAvg > 0.0) CHECK(ErtE > 0.0);
    if (nErAvg > 0.0) CHECK(ErpiE > 0.0);
    if (nErAvg > 0.0) CHECK(ErtetE > 0.0);
    if (nErAvg > 0.0) CHECK(Erml > 0.0);
    if (nErAvg > 0.0) CHECK(Erdl > 0.0);
    if (nErAvg > 0.0) CHECK(Ermul > 0.0);
    if (nErAvg > 0.0) CHECK(Ermus > 0.0);
    if (nErAvg > 0.0) CHECK(Eru0 > 0.0);
    if (nErAvg > 0.0) CHECK(Erfb > 0.0);
    
    CHECK(st.nerr != 0.0);
    CHECK(numd[0] != 0.0);
    CHECK(numd[1] != 0.0);
    CHECK(st.nsim != 0.0);
    
    // NOT an equality. `icon` counts stars that were OBSERVABLE (flagf > 0 and ndw > 2);
    // `numd[0]` counts every record pushed, and the push site sits OUTSIDE that gate, so a
    // star that was drawn but never observable still gets a record. The two are equal only
    // where every draw is observable, which is true in the dense stub patch that every run
    // before commit 81a6b04 used and false as soon as the scan reaches sparse sky -- at
    // l=-3.499, b=-1.98 it is 5 observable out of 15 drawn.
    //
    // Loosening this assertion changes no computed value. It does expose a real question
    // about what the per-sightline denominators mean -- EffiD = numd[0]/nsim is 100% by
    // construction, and the [0] means average over drawn rather than observed stars. That
    // is a science decision, recorded in OPEN_ITEMS.md, not something to change here.
    CHECK(st.icon <= numd[0]);
    CHECK(numd[1] == st.nlens);
     
    cout << "==============================================================" << endl;  
}

