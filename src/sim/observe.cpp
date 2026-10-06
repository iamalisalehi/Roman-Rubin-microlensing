// The light curve as the two surveys record it: the time loop over epochs, with Rubin's ugrizy block,
// Roman's F146 block and the adaptive step size. The hot path -- ~10^5 iterations per draw in the
// Roman footprint -- so the loop's working variables are locals of this function (loaded from the
// EventState on entry, stored back on exit) rather than members reached through a reference.
#include "sim/observe.h"
#include "util/random.h"
#include "events/lightcurve.h"
#include "surveys/noise.h"
#include "run/sample_dump.h"

void simulateLightCurve(SimContext& ctx, const SightlineState& st, EventState& ev) {
    source& s = ctx.s;
    lens& l = ctx.l;
    astromet& as = ctx.as;
    lsst& ls = ctx.ls;
    roman& ro = ctx.ro;
    const RomanSchedule& sched = ctx.sched;
    SampleSpec& dumpSpec = ctx.outs.dumpSpec;
    std::vector<DumpEpoch>& dumpBuf = ctx.outs.dumpBuf;

    // Per-sightline inputs, loop-invariant.
    const int    ndd = st.ndd, nddR = st.nddR;
    const double minc = st.minc, mincR = st.mincR;

    // Hot working variables: this event's accumulators, loaded here and stored back at the end.
    int    ndw = ev.ndw, ndw_L = ev.ndw_L, ndw_R = ev.ndw_R;
    int    flag_det = ev.flag_det, flag_det_L = ev.flag_det_L, flag_det_R = ev.flag_det_R;
    long   nres5_L = ev.nres5_L, nres20_L = ev.nres20_L, nresPSF_L = ev.nresPSF_L;
    long   nres5_R = ev.nres5_R, nres20_R = ev.nres20_R, nresPSF_R = ev.nresPSF_R;
    double dsepMax_L = ev.dsepMax_L, dsepMax_R = ev.dsepMax_R;
    double chi1 = ev.chi1, chi2 = ev.chi2, chi3 = ev.chi3, chi1a = ev.chi1a, chi2a = ev.chi2a, chi3a = ev.chi3a;
    double chi1_L = ev.chi1_L, chi2_L = ev.chi2_L, chi3_L = ev.chi3_L;
    double chi1a_L = ev.chi1a_L, chi2a_L = ev.chi2a_L, chi3a_L = ev.chi3a_L;
    double chi1_R = ev.chi1_R, chi2_R = ev.chi2_R, chi3_R = ev.chi3_R;
    double chi1a_R = ev.chi1a_R, chi2a_R = ev.chi2a_R, chi3a_R = ev.chi3a_R;
    double vsave = ev.vsave;
    const double initial = ev.initial;

    // Loop-private state (reset to these values at the start of every event, never read outside).
    double flag0 = 0.0, flag1 = 0.0, flag2 = 0.0;
    double flag0_L = 0.0, flag1_L = 0.0, flag2_L = 0.0;
    double flag0_R = 0.0, flag1_R = 0.0, flag2_R = 0.0;
    double def1p = 0.0, def2p = 0.0;
    double dt = 60.0;///days
    double cade = 0.0, cadeR = 0.0;
    int    gi, giR, sq, sqR, fi;
    double errs, errg, errsR, errgR, magnio, magnioR, deltaA;
    double Astar0, As1, As0, trajm, trajp, vs1, vs2, sil, sil2, silR, sil2R;
    std::array<double, M> magni, magni0;

    cout << "************** DETECTABLE!!!!!! ********" << endl;
    s.nsdet[s.nums] += 1.0;
    ev.flagf = 1;
    // (The legacy magC/datC/BHLSSTMONTS demo dump that was gated here could never
    // fire -- save < 0 with save = 0 -- and was removed in Deviation 78, with the
    // random draw that fed it.)

    gi = 0;
    gi = 0; giR = 0;
    for (double tim = float(0.0 * year - 100.0 - initial);  tim < float(10.0 * year + 100.0 + initial); tim = tim + dt) {
        // Rubin's geometry (observer on Earth). The quantities derived below and
        // shared across both branches -- Astar0/As0/As1, vs1/vs2, def1p/def2p,
        // trajm/trajp, magni[] -- are all in this frame; the Roman branch
        // recomputes the ones it needs in its own frame (Step H1).
        lightcurve(s, l, as, tim, 0);
        Astar0   = double(s.ut0 * s.ut0 + 2.0) / std::sqrt(s.ut0 * s.ut0 * (s.ut0 * s.ut0 + 4.0)); //MAgnification equation
        s.Astar = double(s.ut  * s.ut  + 2.0) / std::sqrt(s.ut  * s.ut  * (s.ut  * s.ut  + 4.0)); //MAgnification equation
        As0      = double(Astar0   * s.blend[2] + 1.0 - s.blend[2]);
        As1      = double(s.Astar * s.blend[2] + 1.0 - s.blend[2]); //LSST r-band
        vs1      = double(s.mus1 + (s.def1c - def1p) / dt); //[mas/days]
        vs2      = double(s.mus2 + (s.def2c - def2p) / dt); //[mas/days]
        def1p    = s.def1c; //pervious
        def2p    = s.def2c;

        // Computed once per timestep (not just inside the LSST branch) since both
        // instruments' astrometric chi-square terms need it, and it only depends on
        // s->pos1b/pos2b/pos1c/pos2c, which lightcurve() already refreshed above.
        trajm = std::sqrt(s.pos1b * s.pos1b + s.pos2b * s.pos2b); //stright + parallax
        trajp = std::sqrt(s.pos1c * s.pos1c + s.pos2c * s.pos2c); //stright + parallax+lensing


        for (int i = 0; i < M; ++i) {
            magni0[i] = s.magb[i] - 2.5 * std::log10(Astar0   * s.blend[i] + 1.0 - s.blend[i]);
            magni[i]  = s.magb[i] - 2.5 * std::log10(s.Astar * s.blend[i] + 1.0 - s.blend[i]);
        }
        // Rubin's representative-band model magnitude (RUBIN_REF_BANDS, config/parameters.h) --
        // replaces the old hardcoded magni[2] (r-band) at the two use sites below.
        // Reduces to exactly magni[2] when RUBIN_REF_BANDS = {2} (the default), since
        // s->mbs[0]/s->fb[0] were built from the same combination in func_source.
        double magniRubinRef = s.mbs[0] - 2.5 * std::log10(s.Astar * s.fb[0] + 1.0 - s.fb[0]);
        sq = int(ls.ct[gi]);
        sq  = int(ls.ct[gi]);
        sqR = (nddR > 0 and giR < nddR) ? int(ro.ct[giR]) : -1;

        // ---------------- LSST (ugrizy) ----------------

        if (tim >= 0.0 and tim <= Tobs and tim >= ls.tim[int(ls.ct[0])] and tim <= ls.tim[int(ls.ct[ndd - 1])] and
            gi < ndd and sq >= 0 and sq <= static_cast<int>(Nl) and tim >= ls.tim[sq]) {

            fi = int(ls.filter[sq]);

            // Deviation 73: this visit's own depth and saturation, not the SRD
            // minimum -- the same depth that sets errg below.
            const double m5v   = double(ls.sig5[sq]);
            const double satuv = m5v - RUBIN_SATU_BELOW_M5;
            if (magni[fi] >= satuv and magni[fi] <= m5v) {
                errg = errlsstM(magni[fi], int(fi), m5v); //[mag]
                errs = errlsstA(ls, magniRubinRef); ///[mas]

                // Step R1. Could Rubin have told the two images apart at THIS
                // epoch? Inside the magnitude gate on purpose: the paper's
                // criterion counts recorded data points, and an epoch the
                // survey never recorded is not one.
                {
                    const ImagePair ip = imagePair(s.ut, l.tetE, s.magb[fi],
                                                   s.blend[fi], m5v, satuv);
                    if (ip.bothDetectable) {
                        if (ip.sep >= RESOLVE_D_FAINT  * errs)          nres5_L   += 1;
                        if (ip.sep >= RESOLVE_D_BRIGHT * errs)          nres20_L  += 1;
                        if (ip.sep >= FWHM[fi] * ARCSEC_TO_MAS)         nresPSF_L += 1;
                        if (ip.sep >  dsepMax_L)                        dsepMax_L  = ip.sep;
                    }
                }

                deltaA = std::fabs(std::pow(10.0, -0.4 * errg) - 1.0) * (s.blend[fi] * s.Astar + 1.0 - s.blend[fi]);
                magnio = magni[fi] + RandN(errg, NOISE_TRUNC_NSIGMA);

                chi1 += std::fabs((magnio -   magni[fi]) * (magnio -   magni[fi]) / (errg * errg)); //real
                chi2 += std::fabs((magnio -  magni0[fi]) * (magnio -  magni0[fi]) / (errg * errg)); //real no parallax
                chi3 += std::fabs((magnio - s.magb[fi]) * (magnio - s.magb[fi]) / (errg * errg)); //baseline
                chi1_L += std::fabs((magnio -   magni[fi]) * (magnio -   magni[fi]) / (errg * errg));
                chi2_L += std::fabs((magnio -  magni0[fi]) * (magnio -  magni0[fi]) / (errg * errg));
                chi3_L += std::fabs((magnio - s.magb[fi]) * (magnio - s.magb[fi]) / (errg * errg));

                // erra is the per-coordinate sigma (Deviation 71; was drawn and
                // divided as if each coordinate had variance 2 errs^2).
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
                l.soux[ndw] = s.pos1c;
                l.souy[ndw] = s.pos2c;
                l.erra[ndw] = errs;
                l.tele[ndw] = 0; // 0 = LSST
                l.rseas[ndw] = -1; l.rroll[ndw] = -1;

                // Step S1. Everything here was computed above for the
                // detection test; nothing new is drawn. magnio in
                // particular is the noisy datum chi1/chi2/chi3 just used.
                if (dumpSpec.on)
                    dumpBuf.push_back(DumpEpoch{
                        tim, 0, int(fi),
                        magnio, magni[fi], magni0[fi], errg,
                        s.ut, s.ut0, s.Astar, Astar0,
                        s.def1c, s.def2c, s.def1a, s.def2a,
                        s.pos1b, s.pos2b, s.pos1c, s.pos2c,
                        l.pos1,  l.pos2,  errs});

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
                CHECK(gi <= ndd);

                s.errM += deltaA;
                s.errA += errs;
                vsave += std::sqrt(vs1 * vs1 + vs2 * vs2);

                ndw += 1;
                ndw_L += 1;
            }//magnitude limit
            gi += 1;
        }
        // ---------------- Roman (F146) — new ----------------
        // Roman now feeds its own chi1_R/chi2_R/chi3_R (+ astrometric _R,
        // placeholder error — see errsR below) and the joint chi1/chi2/chi3,
        // alongside the Rubin-only _L versions built in the LSST branch above.
        // See Step B1 of JOINT_FIT_REFACTOR_PLAN.md. s->errM/s->errA remain
        // Rubin-only by design (ORIENTATION.md: they're specifically the
        // LSST-only running-error accumulators, unrelated to detection).
        if (nddR > 0 and tim >= 0.0 and tim <= Tobs and
            tim >= ro.tim[int(ro.ct[0])] and tim <= ro.tim[int(ro.ct[nddR - 1])] and
            giR < nddR and sqR >= 0 and sqR <= static_cast<int>(NlRoman) and tim >= ro.tim[sqR]) {

            constexpr int fiR = 6; // F146 — the only Roman band modeled so far

            // ---- Step H1: switch to Roman's observer position ----
            // Everything above was computed with the observer on Earth, which is
            // right for Rubin and wrong for Roman by the L2 offset. Roman sees a
            // slightly different trajectory, so the magnification, the F146
            // magnitude and the astrometric positions all have to be rebuilt here
            // before any of them is recorded. `magni`/`magni0`/`trajm`/`trajp` are
            // per-timestep scratch that the Rubin branch above has already
            // consumed and that the next timestep overwrites, so they are
            // rewritten in place rather than shadowed.
            //
            // No RNG is consumed by this call, so the random stream is untouched
            // and L2_OFFSET_AU = 0 reproduces the pre-H1 run exactly.
            lightcurve(s, l, as, tim, 1);
            {
                const double Astar0R = double(s.ut0 * s.ut0 + 2.0)
                                     / std::sqrt(s.ut0 * s.ut0 * (s.ut0 * s.ut0 + 4.0));
                s.Astar = double(s.ut * s.ut + 2.0)
                         / std::sqrt(s.ut * s.ut * (s.ut * s.ut + 4.0));
                magni0[fiR] = s.magb[fiR] - 2.5 * std::log10(Astar0R  * s.blend[fiR] + 1.0 - s.blend[fiR]);
                magni[fiR]  = s.magb[fiR] - 2.5 * std::log10(s.Astar * s.blend[fiR] + 1.0 - s.blend[fiR]);
                trajm = std::sqrt(s.pos1b * s.pos1b + s.pos2b * s.pos2b);
                trajp = std::sqrt(s.pos1c * s.pos1c + s.pos2c * s.pos2c);
            }

            if (magni[fiR] >= satu[fiR] and magni[fiR] <= thre[fiR]) {
                errgR = errRomanM(ro, magni[fiR]); //[mag] (Deviation 72)

                magnioR = magni[fiR] + RandN(errgR, NOISE_TRUNC_NSIGMA);
                chi1 += std::fabs((magnioR -   magni[fiR]) * (magnioR -   magni[fiR]) / (errgR * errgR));
                chi2 += std::fabs((magnioR -  magni0[fiR]) * (magnioR -  magni0[fiR]) / (errgR * errgR));
                chi3 += std::fabs((magnioR - s.magb[fiR]) * (magnioR - s.magb[fiR]) / (errgR * errgR));
                chi1_R += std::fabs((magnioR -   magni[fiR]) * (magnioR -   magni[fiR]) / (errgR * errgR));
                chi2_R += std::fabs((magnioR -  magni0[fiR]) * (magnioR -  magni0[fiR]) / (errgR * errgR));
                chi3_R += std::fabs((magnioR - s.magb[fiR]) * (magnioR - s.magb[fiR]) / (errgR * errgR));

                // Step H4: Roman's own per-exposure astrometric error, replacing the
                // errlsstA() placeholder (Rubin's curve at Roman's magnitude, which had
                // no reason to be right). Constants and sources in config/parameters.h; the model
                // is per EXPOSURE, which is what one row of RomanBaseline.dat is.
                errsR = errRomanA(magni[fiR]); //[mas]

                // Step R1, Roman side. s->ut is Roman's OWN impact parameter
                // here: lightcurve(..., 1) rebuilt the trajectory in the L2
                // frame above, so this is not the Rubin value reused.
                {
                    const ImagePair ip = imagePair(s.ut, l.tetE, s.magb[fiR],
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
                l.magn[ndw] = magni[fiR]; // F146 magnitude — NOT magni[2]; FisherM's
                                            // tt==1 branch compares against s.mbs[1]/s.fb[1],
                                            // which are F146-based, so the reference point
                                            // recorded here must be F146 too.
                l.errm[ndw] = errgR;
                l.soux[ndw] = s.pos1c;
                l.souy[ndw] = s.pos2c;
                // Step H4. This used to store `errs` -- the RUBIN astrometric error,
                // left over from whichever Rubin epoch last set it, and in general from
                // a different timestep entirely. So the astrometric Fisher matrix was
                // being weighted by a stale value from the other telescope, not even by
                // the errlsstA(magni[fiR]) the comment above it described: errsR was
                // computed for the chi-squared terms and then thrown away. Now Roman's
                // own per-exposure error is both used and stored.
                l.erra[ndw] = errsR;
                l.tele[ndw] = 1; // 1 = Roman/F146
                // Season and roll of this exposure, for the astrometric noise
                // variants' day blocks and frame groups (Deviation 71).
                l.rseas[ndw] = sched.seasonOf(ro.tim[sqR]);
                l.rroll[ndw] = ro.layout[sqR];
                CHECK(l.rseas[ndw] >= 0);

                // Step S1, Roman side. s->ut, s->def* and s->pos* are the
                // L2-frame values lightcurve(..., 1) rebuilt above, not the
                // Rubin ones from earlier in this timestep. Astar0 is
                // recomputed from s->ut0 because the Roman branch's own
                // Astar0R went out of scope before the magnitude gate.
                if (dumpSpec.on)
                    dumpBuf.push_back(DumpEpoch{
                        tim, 1, fiR,
                        magnioR, magni[fiR], magni0[fiR], errgR,
                        s.ut, s.ut0, s.Astar, magnifOf(s.ut0),
                        s.def1c, s.def2c, s.def1a, s.def2a,
                        s.pos1b, s.pos2b, s.pos1c, s.pos2c,
                        l.pos1,  l.pos2,  errsR});

                CHECK(sqR >= 0);
                CHECK(sqR <= int(NlRoman - 1));
                CHECK(errgR >= 0.0);
                CHECK(giR <= nddR);

                ndw += 1;
                ndw_R += 1;
            }//magnitude limit
            giR += 1;
        }

        // ---------------- Adaptive step size ----------------
        // dt must respect whichever instrument's *next* epoch is sooner.
        // Sizing dt to LSST's cadence alone (the original behavior) would
        // silently step right over Roman's much denser epochs.
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
                cade = TIME_STEP_OUTSIDE_LSST_DAYS; //days — matches the original "outside LSST window" fallback
            }

            if (romanInWindow) {
                if (giR > 0 and giR < nddR) cadeR = float(ro.tim[int(ro.ct[giR])] - ro.tim[int(ro.ct[giR - 1])]); //days
                else cadeR = mincR;
            } else {
                cadeR = cade; // Roman not active right now — don't let it constrain dt
            }

            dt = double(std::min(cade, cadeR));
        }

    }//end of loop time

    // Hand the accumulators back.
    ev.ndw = ndw; ev.ndw_L = ndw_L; ev.ndw_R = ndw_R;
    ev.flag_det = flag_det; ev.flag_det_L = flag_det_L; ev.flag_det_R = flag_det_R;
    ev.nres5_L = nres5_L; ev.nres20_L = nres20_L; ev.nresPSF_L = nresPSF_L;
    ev.nres5_R = nres5_R; ev.nres20_R = nres20_R; ev.nresPSF_R = nresPSF_R;
    ev.dsepMax_L = dsepMax_L; ev.dsepMax_R = dsepMax_R;
    ev.chi1 = chi1; ev.chi2 = chi2; ev.chi3 = chi3; ev.chi1a = chi1a; ev.chi2a = chi2a; ev.chi3a = chi3a;
    ev.chi1_L = chi1_L; ev.chi2_L = chi2_L; ev.chi3_L = chi3_L;
    ev.chi1a_L = chi1a_L; ev.chi2a_L = chi2a_L; ev.chi3a_L = chi3a_L;
    ev.chi1_R = chi1_R; ev.chi2_R = chi2_R; ev.chi3_R = chi3_R;
    ev.chi1a_R = chi1a_R; ev.chi2a_R = chi2a_R; ev.chi3a_R = chi3a_R;
    ev.vsave = vsave;
}
