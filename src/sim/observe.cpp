// The light curve as the two surveys record it: the time loop over epochs, with Rubin's ugrizy block,
// Roman's F146 block and the adaptive step size. The hot path -- ~10^5 epochs per draw in the
// Roman footprint -- so the accounting loop's working variables are locals of this function (set to the
// LightCurveStats defaults on entry, copied into the returned struct on exit) rather than members
// reached through a reference.
//
// The loop runs in three phases, because the model is meant to come from a library that computes whole
// light curves and the time grid does not depend on the model:
//   1. the time grid: every epoch, with the visit each survey would take there. The step length comes
//      only from the two visit cadences, and the visit pointers advance whether or not a visit passes
//      the magnitude gate, so the grid is the same for any model and is laid out first (buildTimeGrid);
//   2. the model: the grid handed over in blocks, as arrays of epochs, one observer per call --
//      Rubin's frame at every epoch, Roman's at the epochs where Roman takes a visit (evaluateTrack);
//   3. the accounting: epoch by epoch in order, the magnitude gates, the noise draws, the chi-squared
//      sums, the recorded epochs, the resolution test and the sample dump, all reading the model
//      values phase 2 produced.
// The random numbers are drawn only in phase 3, in the order of the epochs.
#include "sim/observe.h"
#include "util/random.h"
#include "events/lightcurve.h"
#include "surveys/noise.h"
#include "run/sample_dump.h"

namespace {

// One epoch of the time grid, and what the two surveys do there.
struct TimeStep {
    double tim;         // [days]
    double dt;          // the step length in force at this epoch (chosen at the previous one) [days]
    int    sq, sqR;     // the next visit due in Rubin's / Roman's list; sqR = -1 if Roman has none left
    bool   rubin;       // Rubin's branch takes the visit sq here
    bool   roman;       // Roman's branch takes the visit sqR here
};

// Epochs handed to the model per call. Large enough that a library's per-call overhead is negligible,
// small enough that the model values (~150 bytes an epoch) stay in cache until the accounting reads them.
constexpr std::size_t MODEL_BLOCK = 1024;

// The per-event working arrays, kept between events so the hot path does not allocate.
struct TimeLoopBuffers {
    std::vector<TimeStep>   grid;                 // the whole event's time grid
    std::vector<double>     tRubin, tRoman;       // one block's epochs, per observer
    std::vector<TrackPoint> atRubin, atRoman;     // the model at those epochs
};

// Phase 1. Walks the event's clock from 100 days before the survey to 100 days after, a step at a time,
// the step being the sooner of the two surveys' next-visit gaps (60 days outside the survey). Records
// each epoch with the visit each survey takes there.
void buildTimeGrid(const SightlineState& st, const lsst& ls, const roman& ro, std::vector<TimeStep>& grid)
{
    const int    ndd = st.ndd, nddR = st.nddR;
    const double minc = st.minc, mincR = st.mincR;
    const double initial = 0.0;   // legacy padding of the time window [days]; always 0

    double dt = 60.0;///days
    double cade = 0.0, cadeR = 0.0;
    int    gi = 0, giR = 0;       // the next visit due in each survey's matched list

    grid.clear();
    for (double tim = float(0.0 * year - 100.0 - initial);  tim < float(10.0 * year + 100.0 + initial); tim = tim + dt) {
        TimeStep step;
        step.tim = tim;
        step.dt  = dt;
        step.sq  = int(ls.ct[gi]);
        step.sqR = (nddR > 0 and giR < nddR) ? int(ro.ct[giR]) : -1;
        const int sq = step.sq, sqR = step.sqR;

        step.rubin = (tim >= 0.0 and tim <= Tobs and tim >= ls.tim[int(ls.ct[0])] and tim <= ls.tim[int(ls.ct[ndd - 1])] and
                      gi < ndd and sq >= 0 and sq <= static_cast<int>(Nl) and tim >= ls.tim[sq]);
        step.roman = (nddR > 0 and tim >= 0.0 and tim <= Tobs and
                      tim >= ro.tim[int(ro.ct[0])] and tim <= ro.tim[int(ro.ct[nddR - 1])] and
                      giR < nddR and sqR >= 0 and sqR <= static_cast<int>(NlRoman) and tim >= ro.tim[sqR]);
        // A survey's pointer moves past its visit when it takes it, whether or not the visit then passes
        // the magnitude gate.
        if (step.rubin) { gi  += 1; CHECK(gi  <= ndd);  }
        if (step.roman) { giR += 1; CHECK(giR <= nddR); }
        grid.push_back(step);

        // ---------------- Adaptive step size ----------------
        // dt must respect whichever instrument's next epoch is sooner; LSST's cadence alone
        // would step over Roman's much denser epochs.
        if ( tim < -50.0 or tim > (10.0 * year + 50.0)) {
            dt = 60.0; //days
        }

        else {
            bool lsstInWindow  = (tim >= ls.tim[int(ls.ct[0])] and tim <= ls.tim[int(ls.ct[ndd - 1])]);
            bool romanInWindow = (nddR > 0 and tim >= ro.tim[int(ro.ct[0])] and tim <= ro.tim[int(ro.ct[nddR - 1])]);

            if (lsstInWindow) {
                if (gi > 0 and gi < ndd) cade = float(ls.tim[int(ls.ct[gi])] - ls.tim[int(ls.ct[gi - 1])]); //days
                else cade = minc;
            } else {
                cade = TIME_STEP_OUTSIDE_LSST_DAYS; //days
            }

            if (romanInWindow) {
                if (giR > 0 and giR < nddR) cadeR = float(ro.tim[int(ro.ct[giR])] - ro.tim[int(ro.ct[giR - 1])]); //days
                else cadeR = mincR;
            } else {
                cadeR = cade; // Roman not active right now — don't let it constrain dt
            }

            dt = double(std::min(cade, cadeR));
        }
    }
}

} // namespace

LightCurveStats simulateLightCurve(SimContext& ctx, const SightlineState& st) {
    source& s = ctx.s;
    lens& l = ctx.l;
    astromet& as = ctx.as;
    lsst& ls = ctx.ls;
    roman& ro = ctx.ro;
    const RomanSchedule& sched = ctx.sched;
    SampleSpec& dumpSpec = ctx.outs.dumpSpec;
    std::vector<DumpEpoch>& dumpBuf = ctx.outs.dumpBuf;

    static thread_local TimeLoopBuffers buf;
    std::vector<TimeStep>&   grid    = buf.grid;
    std::vector<double>&     tRubin  = buf.tRubin;
    std::vector<double>&     tRoman  = buf.tRoman;
    std::vector<TrackPoint>& atRubin = buf.atRubin;
    std::vector<TrackPoint>& atRoman = buf.atRoman;

    // Hot working variables: this event's accumulators, starting from the LightCurveStats defaults and
    // copied into the returned struct at the end.
    int    ndw = 0, ndw_L = 0, ndw_R = 0;
    int    flag_det = 0, flag_det_L = 0, flag_det_R = 0;   // the loop sets the joint flag, but detectEvent
                                                            // recomputes it, so it is not returned
    long   nres5_L = 0, nres20_L = 0, nresPSF_L = 0;
    long   nres5_R = 0, nres20_R = 0, nresPSF_R = 0;
    double dsepMax_L = -1.0, dsepMax_R = -1.0;
    double chi1 = 0, chi2 = 0, chi3 = 0, chi1a = 0, chi2a = 0, chi3a = 0;
    double chi1_L = 0, chi2_L = 0, chi3_L = 0;
    double chi1a_L = 0, chi2a_L = 0, chi3a_L = 0;
    double chi1_R = 0, chi2_R = 0, chi3_R = 0;
    double chi1a_R = 0, chi2a_R = 0, chi3a_R = 0;
    double vsave = 0.0;

    // Loop-private state (reset to these values at the start of every event, never read outside).
    double flag0 = 0.0, flag1 = 0.0, flag2 = 0.0;
    double flag0_L = 0.0, flag1_L = 0.0, flag2_L = 0.0;
    double flag0_R = 0.0, flag1_R = 0.0, flag2_R = 0.0;
    double def1p = 0.0, def2p = 0.0;
    int    sq, sqR, fi;
    double errs, errg, errsR, errgR, magnio, magnioR, deltaA;
    double trajm, trajp, sil, sil2, silR, sil2R;
    std::array<double, M> magni, magni0;

    cout << "************** DETECTABLE!!!!!! ********" << endl;
    s.nsdet[s.nums] += 1.0;

    // Phase 1: the time grid.
    buildTimeGrid(st, ls, ro, grid);

    const std::size_t nSteps = grid.size();
    for (std::size_t first = 0; first < nSteps; first += MODEL_BLOCK) {
        const std::size_t last = std::min(nSteps, first + MODEL_BLOCK);   // this block is [first, last)

        // Phase 2: the model over the block. Rubin's frame at every epoch; Roman's only where Roman's
        // branch takes a visit (the accounting reads the next Roman value each time it does).
        tRubin.clear();
        tRoman.clear();
        for (std::size_t j = first; j < last; ++j) {
            tRubin.push_back(grid[j].tim);
            if (grid[j].roman) tRoman.push_back(grid[j].tim);
        }
        atRubin.resize(tRubin.size());
        atRoman.resize(tRoman.size());
        evaluateTrack(s, l, as, 0, tRubin.data(), int(tRubin.size()), atRubin.data());
        if (not tRoman.empty())
            evaluateTrack(s, l, as, 1, tRoman.data(), int(tRoman.size()), atRoman.data());

        // Phase 3: the accounting over the block.
        std::size_t jRoman = 0;   // which of atRoman is the next Roman visit's
        for (std::size_t j = first; j < last; ++j) {
            const TimeStep& step = grid[j];
            const double tim = step.tim;
            const double dt  = step.dt;
            // Rubin's geometry (observer on Earth): r. The Roman branch uses its own, q. Only the deflection
            // is carried through every epoch -- the source's velocity is a difference between consecutive
            // epochs' deflections, whichever survey observes them; everything else is computed in the
            // branch that uses it.
            const TrackPoint& r = atRubin[j - first];
            const double def1prev = def1p, def2prev = def2p;
            def1p = r.def1c; //pervious
            def2p = r.def2c;
            sq  = step.sq;
            sqR = step.sqR;

            // ---------------- LSST (ugrizy) ----------------

            if (step.rubin) {

                fi = int(ls.filter[sq]);

                // The model magnitude in this visit's filter, with and without the microlensing parallax.
                magni0[fi] = s.magb[fi] - 2.5 * std::log10(r.ANoPlx * s.blend[fi] + 1.0 - s.blend[fi]);
                magni[fi]  = s.magb[fi] - 2.5 * std::log10(r.A * s.blend[fi] + 1.0 - s.blend[fi]);

                // This visit's own depth and saturation (not the SRD minimum), the same depth that
                // sets errg below.
                const double m5v   = double(ls.sig5[sq]);
                const double satuv = m5v - RUBIN_SATU_BELOW_M5;
                if (magni[fi] >= satuv and magni[fi] <= m5v) {
                    errg = errlsstM(magni[fi], int(fi), m5v); //[mag]
                    errs = errlsstA(errg, double(ls.fwhm[sq])); //[mas], this visit's band, depth and seeing

                    // Rubin's representative-band model magnitude (RUBIN_REF_BANDS, config/parameters.h);
                    // equals magni[2] (r-band) for the default {2}, since s.mbs[0]/s.fb[0] are built from
                    // the same combination in func_source.
                    const double magniRubinRef = s.mbs[0] - 2.5 * std::log10(r.A * s.fb[0] + 1.0 - s.fb[0]);
                    const double vs1 = double(s.mus1 + (r.def1c - def1prev) / dt); //[mas/days]
                    const double vs2 = double(s.mus2 + (r.def2c - def2prev) / dt); //[mas/days]
                    // Both instruments' astrometric chi-square terms need these.
                    trajm = std::sqrt(r.pos1b * r.pos1b + r.pos2b * r.pos2b); //stright + parallax
                    trajp = std::sqrt(r.pos1c * r.pos1c + r.pos2c * r.pos2c); //stright + parallax+lensing

                    // Could Rubin have told the two images apart at this epoch? Inside the magnitude
                    // gate on purpose: the criterion counts recorded data points.
                    {
                        const ImagePair ip = imagePair(r.u, l.tetE, s.magb[fi],
                                                       s.blend[fi], m5v, satuv);
                        if (ip.bothDetectable) {
                            if (ip.sep >= RESOLVE_D_FAINT  * errs)          nres5_L   += 1;
                            if (ip.sep >= RESOLVE_D_BRIGHT * errs)          nres20_L  += 1;
                            if (ip.sep >= FWHM[fi] * ARCSEC_TO_MAS)         nresPSF_L += 1;
                            if (ip.sep >  dsepMax_L)                        dsepMax_L  = ip.sep;
                        }
                    }

                    deltaA = std::fabs(std::pow(10.0, -0.4 * errg) - 1.0) * (s.blend[fi] * r.A + 1.0 - s.blend[fi]);
                    magnio = magni[fi] + RandN(errg, NOISE_TRUNC_NSIGMA);

                    chi1 += std::fabs((magnio -   magni[fi]) * (magnio -   magni[fi]) / (errg * errg)); //real
                    chi2 += std::fabs((magnio -  magni0[fi]) * (magnio -  magni0[fi]) / (errg * errg)); //real no parallax
                    chi3 += std::fabs((magnio - s.magb[fi]) * (magnio - s.magb[fi]) / (errg * errg)); //baseline
                    chi1_L += std::fabs((magnio -   magni[fi]) * (magnio -   magni[fi]) / (errg * errg));
                    chi2_L += std::fabs((magnio -  magni0[fi]) * (magnio -  magni0[fi]) / (errg * errg));
                    chi3_L += std::fabs((magnio - s.magb[fi]) * (magnio - s.magb[fi]) / (errg * errg));

                    // erra is the per-coordinate astrometric sigma.
                    sil  = RandN(errs, NOISE_TRUNC_NSIGMA);
                    sil2 = RandN(errs, NOISE_TRUNC_NSIGMA);

                    chi1a += std::fabs((trajp + sil - trajp) * (trajp + sil - trajp) / (errs * errs)); //real trajectory
                    chi2a += std::fabs((trajp + sil - trajm) * (trajp + sil - trajm) / (errs * errs)); //real without lensing
                    chi3a += std::fabs( sil2  * sil2 / (errs * errs));
                    chi1a_L += std::fabs((trajp + sil - trajp) * (trajp + sil - trajp) / (errs * errs));
                    chi2a_L += std::fabs((trajp + sil - trajm) * (trajp + sil - trajm) / (errs * errs));
                    chi3a_L += std::fabs( sil2  * sil2 / (errs * errs));

                    CHECK(ndw < coun); // array-bounds guard — see note on `coun` sizing
                    l.timn[ndw] = tim;
                    l.magn[ndw] = magniRubinRef; // Rubin representative-band model value (RUBIN_REF_BANDS)
                    l.errm[ndw] = errg;
                    l.soux[ndw] = r.pos1c;
                    l.souy[ndw] = r.pos2c;
                    l.erra[ndw] = errs;
                    l.tele[ndw] = 0; // 0 = LSST
                    l.rseas[ndw] = -1; l.rroll[ndw] = -1;

                    // Sample dump: reuses the values computed above (magnio is the noisy datum the
                    // chi-squared sums just used); nothing new is drawn.
                    if (dumpSpec.on)
                        dumpBuf.push_back(DumpEpoch{
                            tim, 0, int(fi),
                            magnio, magni[fi], magni0[fi], errg,
                            r.u, r.uNoPlx, r.A, r.ANoPlx,
                            r.def1c, r.def2c, r.def1a, r.def2a,
                            r.pos1b, r.pos2b, r.pos1c, r.pos2c,
                            r.lens1, r.lens2, errs});

                    flag2 = 0.0;
                    if (std::fabs(magnio - s.magb[fi]) > std::fabs(OUTLIER_FLAG_NSIGMA * errg))    flag2 = 1.0;
                    if (ndw_L > 2 and float(flag0 + flag1 + flag2) > OUTLIER_RUN_THRESHOLD)   flag_det = 1;
                    flag0 = flag1;
                    flag1 = flag2;
                    flag2_L = flag2; // identical test; kept as an explicit Rubin-labeled copy
                    if (ndw_L > 2 and float(flag0_L + flag1_L + flag2_L) > OUTLIER_RUN_THRESHOLD)   flag_det_L = 1;
                    flag0_L = flag1_L;
                    flag1_L = flag2_L;


                    CHECK(sq >= 0);
                    CHECK(sq <= int(Nl - 1));
                    CHECK(fi >= 0);
                    CHECK(fi <= 5);
                    CHECK(errs >= 0.0);
                    CHECK(errg >= 0.0);
                    CHECK(deltaA >= 0.0);

                    s.errM += deltaA;
                    s.errA += errs;
                    vsave += std::sqrt(vs1 * vs1 + vs2 * vs2);

                    ndw += 1;
                    ndw_L += 1;
                }//magnitude limit
            }
            // ---------------- Roman (F146) ----------------
            // Feeds chi1_R/chi2_R/chi3_R (and the astrometric _R sums) and the joint chi1/chi2/chi3.
            // s.errM/s.errA stay Rubin-only: they are LSST running-error accumulators unrelated to detection.
            if (step.roman) {

                constexpr int fiR = 6; // F146, the only Roman band modeled

                // Switch to Roman's observer position. The Rubin-frame values above are off by the L2
                // offset, so the magnification, F146 magnitude and astrometric positions are rebuilt
                // here from the Roman-frame model. magni/magni0/trajm/trajp are per-timestep scratch,
                // rewritten in place.
                const TrackPoint& q = atRoman[jRoman];
                {
                    magni0[fiR] = s.magb[fiR] - 2.5 * std::log10(q.ANoPlx * s.blend[fiR] + 1.0 - s.blend[fiR]);
                    magni[fiR]  = s.magb[fiR] - 2.5 * std::log10(q.A * s.blend[fiR] + 1.0 - s.blend[fiR]);
                    trajm = std::sqrt(q.pos1b * q.pos1b + q.pos2b * q.pos2b);
                    trajp = std::sqrt(q.pos1c * q.pos1c + q.pos2c * q.pos2c);
                }

                if (magni[fiR] >= satu[fiR] and magni[fiR] <= thre[fiR]) {
                    errgR = errRomanM(ro, magni[fiR]); //[mag]

                    magnioR = magni[fiR] + RandN(errgR, NOISE_TRUNC_NSIGMA);
                    chi1 += std::fabs((magnioR -   magni[fiR]) * (magnioR -   magni[fiR]) / (errgR * errgR));
                    chi2 += std::fabs((magnioR -  magni0[fiR]) * (magnioR -  magni0[fiR]) / (errgR * errgR));
                    chi3 += std::fabs((magnioR - s.magb[fiR]) * (magnioR - s.magb[fiR]) / (errgR * errgR));
                    chi1_R += std::fabs((magnioR -   magni[fiR]) * (magnioR -   magni[fiR]) / (errgR * errgR));
                    chi2_R += std::fabs((magnioR -  magni0[fiR]) * (magnioR -  magni0[fiR]) / (errgR * errgR));
                    chi3_R += std::fabs((magnioR - s.magb[fiR]) * (magnioR - s.magb[fiR]) / (errgR * errgR));

                    // Roman's per-exposure astrometric error (one row of RomanBaseline.dat is one
                    // exposure), from this exposure's photometric SNR; config/parameters.h.
                    errsR = errRomanA(errgR); //[mas]

                    // Resolution test, Roman side. q.u is Roman's own impact parameter here
                    // (trajectory in the L2 frame).
                    {
                        const ImagePair ip = imagePair(q.u, l.tetE, s.magb[fiR],
                                                       s.blend[fiR], thre[fiR], satu[fiR]);
                        if (ip.bothDetectable) {
                            if (ip.sep >= RESOLVE_D_FAINT  * errsR)         nres5_R   += 1;
                            if (ip.sep >= RESOLVE_D_BRIGHT * errsR)         nres20_R  += 1;
                            if (ip.sep >= FWHM[fiR] * ARCSEC_TO_MAS)        nresPSF_R += 1;
                            if (ip.sep >  dsepMax_R)                        dsepMax_R  = ip.sep;
                        }
                    }

                    silR  = RandN(errsR, NOISE_TRUNC_NSIGMA);
                    sil2R = RandN(errsR, NOISE_TRUNC_NSIGMA);
                    chi1a += std::fabs((trajp + silR - trajp) * (trajp + silR - trajp) / (errsR * errsR));
                    chi2a += std::fabs((trajp + silR - trajm) * (trajp + silR - trajm) / (errsR * errsR));
                    chi3a += std::fabs( sil2R * sil2R / (errsR * errsR));
                    chi1a_R += std::fabs((trajp + silR - trajp) * (trajp + silR - trajp) / (errsR * errsR));
                    chi2a_R += std::fabs((trajp + silR - trajm) * (trajp + silR - trajm) / (errsR * errsR));
                    chi3a_R += std::fabs( sil2R * sil2R / (errsR * errsR));

                    flag2_R = 0.0;
                    if (std::fabs(magnioR - s.magb[fiR]) > std::fabs(OUTLIER_FLAG_NSIGMA * errgR))    flag2_R = 1.0;
                    if (ndw_R > 2 and float(flag0_R + flag1_R + flag2_R) > OUTLIER_RUN_THRESHOLD)   flag_det_R = 1;
                    flag0_R = flag1_R;
                    flag1_R = flag2_R;

                    CHECK(ndw < coun); // array-bounds guard — see note on `coun` sizing
                    l.timn[ndw] = tim;
                    l.magn[ndw] = magni[fiR]; // F146 magnitude, not magni[2]: FisherM's tt==1 branch
                                                // compares against the F146-based s.mbs[1]/s.fb[1]
                    l.errm[ndw] = errgR;
                    l.soux[ndw] = q.pos1c;
                    l.souy[ndw] = q.pos2c;
                    l.erra[ndw] = errsR;
                    l.tele[ndw] = 1; // 1 = Roman/F146
                    // Season and roll of this exposure, for the astrometric noise variants' day blocks
                    // and frame groups.
                    l.rseas[ndw] = sched.seasonOf(ro.tim[sqR]);
                    l.rroll[ndw] = ro.layout[sqR];
                    CHECK(l.rseas[ndw] >= 0);

                    // Sample dump, Roman side: u, def*, pos* and the magnifications are the L2-frame values.
                    if (dumpSpec.on)
                        dumpBuf.push_back(DumpEpoch{
                            tim, 1, fiR,
                            magnioR, magni[fiR], magni0[fiR], errgR,
                            q.u, q.uNoPlx, q.A, q.ANoPlx,
                            q.def1c, q.def2c, q.def1a, q.def2a,
                            q.pos1b, q.pos2b, q.pos1c, q.pos2c,
                            q.lens1, q.lens2, errsR});

                    CHECK(sqR >= 0);
                    CHECK(sqR <= int(NlRoman - 1));
                    CHECK(errgR >= 0.0);

                    ndw += 1;
                    ndw_R += 1;
                }//magnitude limit
                jRoman += 1;
            }

        }//end of block's epochs
    }//end of blocks

    // Hand the accumulators back.
    (void)flag_det;
    LightCurveStats lc;
    lc.flagf = 1;
    lc.ndw = ndw; lc.ndw_L = ndw_L; lc.ndw_R = ndw_R;
    lc.flag_det_L = flag_det_L; lc.flag_det_R = flag_det_R;
    lc.nres5_L = nres5_L; lc.nres20_L = nres20_L; lc.nresPSF_L = nresPSF_L;
    lc.nres5_R = nres5_R; lc.nres20_R = nres20_R; lc.nresPSF_R = nresPSF_R;
    lc.dsepMax_L = dsepMax_L; lc.dsepMax_R = dsepMax_R;
    lc.chi1 = chi1; lc.chi2 = chi2; lc.chi3 = chi3; lc.chi1a = chi1a; lc.chi2a = chi2a; lc.chi3a = chi3a;
    lc.chi1_L = chi1_L; lc.chi2_L = chi2_L; lc.chi3_L = chi3_L;
    lc.chi1a_L = chi1a_L; lc.chi2a_L = chi2a_L; lc.chi3a_L = chi3a_L;
    lc.chi1_R = chi1_R; lc.chi2_R = chi2_R; lc.chi3_R = chi3_R;
    lc.chi1a_R = chi1a_R; lc.chi2a_R = chi2a_R; lc.chi3a_R = chi3a_R;
    lc.vsave = vsave;
    return lc;
}
