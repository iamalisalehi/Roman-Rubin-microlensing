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
    // Sky area this sightline stands for; written into every event row, and statistics pooled across
    // sightlines must weight by it.
    st.wArea = sightline.area;
    // nri/nde: longitude column index, and index within that column (recorded in the map file and
    // event table). nri is assigned, not incremented, so a fully cut column cannot renumber the rest;
    // under stratification it is a fine-grid index stepping by kSub between coarse columns.
    st.nri = sightline.col;
    if (sightline.col != st.lastCol) { st.nde = -1; st.lastCol = sightline.col; }
    st.nde += 1;

    // --start-index resumes an interrupted run. The nri/nde bookkeeping above also runs for skipped
    // sightlines, so a resumed run concatenates onto the interrupted one. The index to resume at is
    // the number of "NEW STEP" lines in the interrupted log (it counts every sightline the scan
    // reached, including no-coverage and barren ones); the map file counts fewer and would duplicate rows.
    if (st.iScan < cfg.startIndex) return SightlineStart::Skip;
    if (cfg.endIndex >= 0 and st.iScan >= cfg.endIndex) return SightlineStart::Stop;
    // Each sightline has its own random stream (sightlineSeed in include/util/random.h), so a
    // resumed run reproduces an uninterrupted one.
    rng.seed(sightlineSeed(cfg.seedBase, st.iScan));
    cout << ">>>>>>>>>>> NEW STEP " << st.nde << " <<<<<<<<\t nri:  " << st.nri << endl;
    cout << "longtitude: " << s.lon << "\t latitude: " << s.lat << endl;


    // Coverage is on LSSTCam's active silicon for this visit's pointing and rotation, after a cheap
    // (l, b) distance cut.
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

    // Rubin's depth per band at this sightline: the median of its matched visits' own 5-sigma depths,
    // used by the pre-selection. A band with no visit gets -inf, so it cannot count toward
    // "detectable in >= 2 bands".
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

    // Empty sky: with no epochs no event can be detected and nerr never advances, so the draw loop
    // would not terminate. Skipping is exact (zero events), and the aggregation CHECKs would abort on
    // an empty field. The skipped area is not written to the map file, so converting Neven to a total
    // count must use the aggregated sightline count, not the grid size.
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

    s.TET = (360.0 - s.lon) / RAa;///radian
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

    // Whether the sightline met its budget or the cap stopped it: a capped sightline's Poisson
    // precision is what it reached, not what was asked for.
    const bool budgetMet = (st.icon  >= cfg.iconTarget and
                            st.nlens >= cfg.nlensTarget and
                            st.nerr  >= cfg.nerrTarget);
    if (!budgetMet) {
        run.nCapped += 1;
        cout << "  CAP: stopped at nsim = " << st.nsim << " with icon = " << st.icon
             << "/" << cfg.iconTarget << ", nlens = " << st.nlens << "/" << cfg.nlensTarget
             << ", nerr = " << st.nerr << "/" << cfg.nerrTarget << endl;
    }

    // Epochs existed but nothing survived to be aggregated: skip, since the CHECKs below need at
    // least one detected and one characterised event.
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
        // Replayed too: this loop runs after the whole field, so co.okA would otherwise hold the
        // last FisherM call's value.
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

            // Average only over events whose Fisher results are all real measurements; the -1.0
            // sentinel (and ErrorCal dividing it by the parameter value) would otherwise pull large
            // negative numbers into the means. This is a mean precision over measured events, not a
            // joint-vs-single ratio, so skipping unmeasured events does not bias it; nErAvg counts
            // the events actually averaged. okA[SJOINT] alone is not sufficient: it says the joint
            // matrix inverted, not that every parameter was fitted (activePhotParams in
            // include/fisher/fisher.h), so e.g. a Roman-only detection has Era[2] = -1.0 in a valid
            // joint fit. Testing the values keeps every mean over the same event set.
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

    for (int i = 0; i < 2; ++i) { // turn the sums into means ([0] all records, [1] flagged)
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
    EffiD = double(st.icon * 100.0 / (st.nsim    + eps)); // % of drawn stars that are visible
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
         // Appended columns: w_area (deg^2 this sightline stands for), lon, lat, so an event joins its
         // sightline on (lon, lat) and the pooled weight w_area/nsim is computable.
         << " " << std::setprecision(8) << st.wArea
         << " " << std::setprecision(6) << s.lon << " " << s.lat
         // Then six columns: this sightline's median 5-sigma depth in each of ugrizy (rubinDepthMed,
         // -inf for a band with no visit). preselectEvent compares peaks with these depths, so
         // romanlib.acceptance_probability needs them to rebuild which draws were kept.
         << std::fixed << std::setprecision(6);
    for (int b = 0; b < 6; ++b) fil3 << " " << st.rubinDepthMed[b];
    fil3 << "\n";

    // Flush per sightline: a killed run otherwise loses buffered rows and leaves a half-written one
    // that the resumed run would append to. A sightline costs minutes, so the flushes are free.
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
    
    // Not an equality: icon counts observable stars (flagf > 0 and ndw > 2), while numd[0] counts
    // every record pushed, including draws that were never observable. Consequently EffiD =
    // numd[0]/nsim is 100% by construction and the [0] means average over drawn, not observed,
    // stars; which of the two the per-sightline denominators should be is an open question.
    CHECK(st.icon <= numd[0]);
    CHECK(numd[1] == st.nlens);
     
    cout << "==============================================================" << endl;  
}

