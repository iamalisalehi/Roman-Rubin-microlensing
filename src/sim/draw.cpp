// Drawing an event: source, lens, efficiency bins, resets, per-survey pre-selection.
#include "sim/draw.h"
#include "util/random.h"
#include "events/source.h"
#include "events/lens.h"
#include "galaxy/density.h"
#include "run/histograms.h"
#include "surveys/noise.h"

EfficiencyBins drawEvent(SimContext& ctx, SightlineState& st, int prevNdw) {
    source& s = ctx.s;
    lens& l = ctx.l;
    CMD& cm = ctx.cm;
    extin& ex = ctx.ex;
    std::vector<DumpEpoch>& dumpBuf = ctx.outs.dumpBuf;
    RunTotals& run = ctx.run;
    EfficiencyBins bins;

    st.nsim += 1.0;
    func_source(s, cm, ex, st.sightlineIdx);
    func_lens(l, s, cm, ex, st.sightlineIdx);
//                std::cerr << "nsim=" << nsim << "  Ds=" << s->Ds << "  mass=" << s->mass
//                          << "  nums=" << s->nums << "  Ml=" << l->Ml << "  u0=" << l->u0 << "\n";
    optical_depth(s);

    // tE-histogram bin for THIS event, computed for every draw rather than only
    // for detected ones: nstE is the denominator of the detection efficiency, so
    // it has to count everything simulated. Previously gg was only evaluated
    // inside the detection branch and neither counter was ever incremented, so
    // ndtE stayed identically zero and EFF, Gamma and Neven with it (the run then
    // aborted on CHECK(EFF > 0.0) as soon as a field managed to complete).
    bins.gg = FunctE(l);
    l.nstE[bins.gg] += 1.0;
    run.nSimTot += 1;

    // The other six efficiency axes, on the same principle as tE: the DENOMINATOR
    // has to count everything simulated, so the bin is computed here, before the
    // detection test, and the numerator is incremented in the detection branch
    // using these same indices. Every input is already set for this draw --
    // func_source filled s->Map and s->blend, func_lens filled u0, pirel and
    // murel -- so this is six array lookups and no new physics.
    //
    // Unlike the tE pair there is no per-sightline lowercase pair for these; the
    // N* arrays accumulate over the whole run, which is what the EfLMC writer
    // reports and why nothing resets them per sightline.
    bins.ss = FuncMl(l);                    // lens mass
    bins.qq = FuncPi(l);                    // log10 relative parallax
    bins.ww = Funcu0(l);                    // impact parameter
    bins.vv = FuncMu(l);                    // relative proper motion
    bins.zz = FuncMb(l, s.Map[2]);         // source baseline magnitude, r band
    bins.pp = FuncFb(l, s.blend[2]);       // blend fraction, r band
    l.NsMl[bins.ss] += 1.0;
    l.Nspi[bins.qq] += 1.0;
    l.Nsu0[bins.ww] += 1.0;
    l.Nsmu[bins.vv] += 1.0;
    l.Nsmb[bins.zz] += 1.0;
    l.Nsfb[bins.pp] += 1.0;

    s.nssim[s.nums] += 1.0;
    dumpBuf.clear(); //Step S1: this draw's epoch buffer. See the note on `ndw`.
    // (Step B2: the old single `test = RandR(0.0,1.0)` draw consumed here by
    // `test <= s->blend[2]` is gone — testL/testR are now drawn fresh right
    // before the per-survey pre-selection check, below.)

    // Clear only the prefix the PREVIOUS event dirtied -- `prevNdw`, the
    // ndw of its LightCurveStats, which the caller carries from draw to draw. Clearing all
    // `coun` slots (Nl + NlRoman = 306,092, times seven arrays = 2.1M writes)
    // to reset the ~2,000 an event actually uses was ~150x of pure waste per
    // draw, and got 16x more expensive when NlRoman became the real visit count.
    //
    // Safe because every slot in [0, ndw) is fully written before it is read --
    // both fill branches write all seven arrays at index ndw before incrementing
    // it -- and FisherM reads only [0, ndw). Slots past the previous ndw are
    // therefore untouched since construction, i.e. already zero.
    //
    // Kept as a prefix clear rather than deleted outright so the clean-slate
    // invariant survives: if a future edit ever advances ndw without filling
    // every array, that shows up as a zero instead of as the previous event's
    // photometry silently entering this event's Fisher matrix.
    //
    // NOTE: the untouched tail of tele[] is 0 (its constructed value), not the
    // -1 the old full clear wrote. Unobservable today -- tele is only read at
    // [0, ndw) -- but relevant if anything ever scans the whole array.
    for (int i = 0; i < prevNdw; ++i) {
        l.timn[i] = 0.0;  l.magn[i] = 0.0; l.soux[i] = 0.0;  l.souy[i] = 0.0;
        l.errm[i] = 0.0;  l.erra[i] = 0.0; l.tele[i] = -1;
    }

    // The per-event accumulators (LightCurveStats, Detection) are fresh value structs, zeroed by their
    // defaults; what lives in the shared objects is reset here.
    s.def1c = 0.0;
    s.def2c = 0.0;
    s.errM = 0.0; s.errA  = 0.0;
    return bins;
}

bool preselectEvent(SimContext& ctx, const SightlineState& st) {
    source& s = ctx.s;
    lens& l = ctx.l;

    double fdetRubin, testL, testR, Mpeak;
    bool   rubinDetectable, romanDetectable, acceptRubin, acceptRoman;

    fdetRubin = 0.0;
    romanDetectable = false;

    for (int i = 0; i < M; ++i) {
        Mpeak = s.magb[i] - 2.5 * std::log10(l.A0 * s.blend[i] + 1.0 - s.blend[i]);
//                        cout << "i=" << i << "  Mab=" << s->Mab[i] << "  Map=" << s->Map[i]
//                             << "  blend=" << s->blend[i] << "  Mpeak=" << Mpeak << endl;
        if (i < 6) { // LSST ugrizy
            if (Mpeak <= st.rubinDepthMed[i] and s.magb[i] > st.rubinDepthMed[i] - RUBIN_SATU_BELOW_M5)
                fdetRubin += 1.0;
        } else {     // i == 6, Roman F146 — single band, no ">=2 filters" bar applies
            if (Mpeak <= thre[i] and s.magb[i] > satu[i])    romanDetectable = true;
        }
    }
    rubinDetectable = (fdetRubin > 1.0); // at least 2 of the 6 LSST bands

    testL = RandR(0.0, 1.0);
    testR = RandR(0.0, 1.0);
    // Independent draws per survey — reusing one draw for both would correlate
    // the Rubin-accept and Roman-accept decisions for no physical reason. Each
    // draw is weighted by that survey's OWN blend fraction (Step B2), replacing
    // the old single test <= s->blend[2] (LSST r-band only, for both surveys).
    acceptRubin = rubinDetectable and (testL <= s.blend[2]);
    acceptRoman = romanDetectable and (testR <= s.blend[6]);


    return acceptRubin or acceptRoman;
}

