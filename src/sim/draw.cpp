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
    optical_depth(s);

    // Efficiency bins are computed for every draw, before the detection test: the denominators
    // (nstE, Ns*) count everything simulated and the detection branch increments the numerators
    // with the same indices.
    bins.gg = FunctE(l);
    l.nstE[bins.gg] += 1.0;
    run.nSimTot += 1;

    // The other six axes. Their N* arrays accumulate over the whole run (the EfLMC writer reports
    // them), so nothing resets them per sightline.
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
    dumpBuf.clear(); // this draw's sample-dump epoch buffer
    // Clear only the prefix the previous event dirtied (prevNdw): clearing all `coun` slots
    // (Nl + NlRoman, times seven arrays) every draw is far more expensive than the ~2,000 slots an
    // event uses. Safe because every slot in [0, ndw) is fully written before it is read and FisherM
    // reads only [0, ndw); slots past it are still zero from construction. The prefix clear is kept
    // so that an edit advancing ndw without filling every array shows up as a zero rather than as
    // the previous event's photometry. The untouched tail of tele[] is 0, not -1.
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
        if (i < 6) { // LSST ugrizy
            if (Mpeak <= st.rubinDepthMed[i] and s.magb[i] > st.rubinDepthMed[i] - RUBIN_SATU_BELOW_M5)
                fdetRubin += 1.0;
        } else {     // i == 6, Roman F146, single band, no ">=2 filters" bar applies
            if (Mpeak <= thre[i] and s.magb[i] > satu[i])    romanDetectable = true;
        }
    }
    rubinDetectable = (fdetRubin > 1.0); // at least 2 of the 6 LSST bands

    testL = RandR(0.0, 1.0);
    testR = RandR(0.0, 1.0);
    // Independent accept draws per survey (reusing one would correlate the two decisions), each
    // weighted by that survey's own blend fraction.
    acceptRubin = rubinDetectable and (testL <= s.blend[2]);
    acceptRoman = romanDetectable and (testR <= s.blend[6]);

    return acceptRubin or acceptRoman;
}

