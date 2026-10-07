// The detection test (Step H7) and the detection taxonomy / efficiency tallies.
#include "sim/detect.h"
#include "fisher/fisher.h"

void detectEvent(SimContext& ctx, EventState& ev) {
    const RunConfig& cfg = ctx.cfg;
    source& s = ctx.s;

    // Step H7. The lensing statistic is SIGNED: chi3 is the flat-baseline
    // residual and chi1 the lensing-model residual, so chi3 - chi1 > 0 means
    // the lensing model fits better and is the only direction that can count
    // as evidence of lensing. The previous fabs() meant a lensing model that
    // fitted WORSE than the baseline by more than the bar would have been
    // reported as a detection. That never bit in practice, because the data
    // are lensed and chi3 > chi1 for any event with real signal, but it also
    // breaks the monotonicity H7 exists to establish: with fabs, an instrument
    // whose epochs all fall outside the event contributes a small NEGATIVE
    // difference from noise, and |dchi_L + dchi_R| can then fall below
    // max(|dchi_L|, |dchi_R|). Signed, the sum is exact and monotone.
    //
    // dchiP and dchiA are not detection tests (nothing thresholds them); since
    // Deviation 78 they are signed like dchiL, so a negative value (the simpler
    // model fitting better) is visible rather than folded into a "size".
    ev.dchiL   = ev.chi3  - ev.chi1;             //lensing_effect (signed)
    // Deviation 78: signed like dchiL (positive = the full model fits better).
    ev.dchiP   = ev.chi2  - ev.chi1;             //parallax_effect (signed)
    ev.dchiA   = ev.chi2a - ev.chi1a;            //deflection_effect (signed)
    ev.dchiL_L = ev.chi3_L - ev.chi1_L;
    ev.dchiP_L = (ev.chi2_L  - ev.chi1_L);   // signed (Deviation 78)
    ev.dchiA_L = (ev.chi2a_L - ev.chi1a_L);   // signed (Deviation 78)
    ev.dchiL_R = ev.chi3_R - ev.chi1_R;
    ev.dchiP_R = (ev.chi2_R  - ev.chi1_R);   // signed (Deviation 78)
    ev.dchiA_R = (ev.chi2a_R - ev.chi1a_R);   // signed (Deviation 78)
    // flag_det is the JOINT run-test flag the table reports (Deviation 78): it was set in
    // the Rubin branch only, so a Roman-only persistent signal left it at 0.
    ev.flag_det = (ev.flag_det_L > 0 or ev.flag_det_R > 0) ? 1 : 0;

    // Three independent detection tests. detL/detR use each instrument's own
    // epoch count (ndw_L/ndw_R) and run-test result — mixing in the joint ndw
    // here would let Roman's dense epochs silently raise the bar for a purely
    // Rubin-driven signal (and vice versa). detJ requires a persistent run in
    // EITHER instrument's own cadence (flag_det_L or flag_det_R) rather than an
    // interleaved run across two different cadences, which wouldn't mean anything
    // (see Step B1 teaching brief). FFG[0] — the boolean that gates FisherM below,
    // preserving the existing detection-then-Fisher ordering (plan §0.4) — is the
    // union of all three, so an event Rubin alone clearly detects is never dropped
    // from characterization just because the joint-scaled threshold happens to miss.
    if (s.FWHM < Tobs and ev.dchiL_L > cfg.dchiDet and ev.flag_det_L > 0 and ev.ndw_L > 10) ev.detL = 1;
    if (s.FWHM < Tobs and ev.dchiL_R > cfg.dchiDet and ev.flag_det_R > 0 and ev.ndw_R > 10) ev.detR = 1;
    if (s.FWHM < Tobs and ev.dchiL   > cfg.dchiDet and (ev.flag_det_L > 0 or ev.flag_det_R > 0) and ev.ndw > 10) ev.detJ = 1;

    // Monotonicity. Adding data to an analysis cannot destroy information, so
    // if either survey alone clears its detection bar, the combined stream --
    // which contains that survey's data in full, plus more -- must clear its
    // own.
    //
    // Since Step H7 this holds BY CONSTRUCTION rather than by patch. chi1 and
    // chi3 are accumulated over both instruments' epochs in the same loops that
    // fill chi1_L/chi1_R and chi3_L/chi3_R, so
    //
    //     chi1 = chi1_L + chi1_R   and   chi3 = chi3_L + chi3_R
    //
    // exactly, hence dchiL = dchiL_L + dchiL_R with the signed statistic above.
    // All three tests now use the SAME fixed bar, so if either survey's own
    // dchi clears it, the sum -- which is that dchi plus a quantity that is
    // positive whenever the other survey sees any signal at all -- clears it
    // too. The auxiliary gates agree: flag_det_L > 0 implies (flag_det_L or
    // flag_det_R), and ndw_L > 10 implies ndw = ndw_L + ndw_R > 10.
    //
    // The old bar was 2*ndw, a threshold on the MEAN per-epoch chi-squared
    // improvement rather than on total significance, which made pooling a
    // survey with many low-signal epochs RAISE the joint bar without adding
    // signal. In Roman's footprint that lifted the joint bar ~4,728 above
    // Roman's own and vetoed 21.3% of all detections there.
    //
    // detJ_raw and the monotone patch are retained deliberately (H7 acceptance
    // criterion 4) so the anomaly rate stays MEASURABLE and can be shown to
    // have gone to zero rather than been hidden. If DET_ANOMALY is ever
    // non-zero again, the construction above has been broken and the counter
    // is how that gets noticed.
    ev.detJ_raw = ev.detJ;
    if (ev.detL or ev.detR) ev.detJ = 1;
}

void tallyDetection(SimContext& ctx, SightlineState& st, EventState& ev) {
    const RunConfig& cfg = ctx.cfg;
    source& s = ctx.s;
    lens& l = ctx.l;
    RunTotals& run = ctx.run;

    // Detection taxonomy (DetClass, include/fisher/fisher.h). The joint fit is what makes a
    // detection meaningful -- it sees strictly more data than either survey
    // alone -- so the classes are distinguished by which telescopes ALSO
    // detect the event unaided. DET_ANOMALY catches the case that should be
    // impossible, a single-telescope detection the joint test misses.
    ev.dclsEvent = detClass(ev.detL, ev.detR, ev.detJ);
    const int dcls = ev.dclsEvent;
    st.nDetClass[dcls]    += 1;
    run.NDetClassTot[dcls] += 1;
    run.NDetClassTE[ev.gg][dcls] += 1;
    // Diagnostic only: how often the raw joint test contradicted a
    // single-survey detection, before monotonicity was imposed.
    // Both counters, and the run total is the one that was missing: it is
    // read at the end of the run to decide whether to explain DET_ANOMALY,
    // but was never incremented, so every run has reported exactly zero
    // anomalies since the counter was introduced -- including the 2026-08-30
    // production run, whose "0 ANOMALY" is this bug and not a measurement.
    // The per-sightline counter was always right, which is why the anomaly
    // was visible in the log all along and invisible in every summary.
    // Deviation 79: the joint delta-chi2 is accumulated separately from the two
    // per-survey ones; it must equal their sum on every event.
    if (std::fabs((ev.dchiL_L + ev.dchiL_R) - ev.dchiL) > 1e-6 * std::max(1.0, std::fabs(ev.dchiL))) {
        if (run.nDchiMismatch < 10)
            std::cerr << "DCHI_MISMATCH sightline " << st.iScan << " joint " << ev.dchiL
                      << " L+R " << ev.dchiL_L + ev.dchiL_R << "\n";
        run.nDchiMismatch += 1;
    }
    if (!ev.detJ_raw and (ev.detL or ev.detR)) {
        st.nDetClass[DET_ANOMALY]    += 1;
        run.NDetClassTot[DET_ANOMALY] += 1;
        // Deviation 79: log everything needed to explain it. The joint delta-chi2
        // is the SUM of the two surveys' (each accumulated over its own epochs), and
        // either term can be negative through noise -- a survey whose data happen to
        // fit a flat baseline slightly better than the true model -- so a single-
        // survey detection just over the bar can leave the sum just under it. That
        // is physics, not a sign bug; "sum - joint" must be ~0 if the bookkeeping is
        // right, and a non-zero value there WOULD be a bug.
        std::cerr << std::setprecision(10)
                  << "DET_ANOMALY_DETAIL sightline " << st.iScan << " lon " << s.lon
                  << " lat " << s.lat << " | dchiL_L " << ev.dchiL_L << " dchiL_R "
                  << ev.dchiL_R << " sum " << ev.dchiL_L + ev.dchiL_R << " joint " << ev.dchiL
                  << " (sum - joint " << (ev.dchiL_L + ev.dchiL_R) - ev.dchiL << ") bar "
                  << cfg.dchiDet << " | ndw_L " << ev.ndw_L << " ndw_R " << ev.ndw_R
                  << " ndw " << ev.ndw << " | flag_det_L " << ev.flag_det_L
                  << " flag_det_R " << ev.flag_det_R << " | detL " << ev.detL << " detR "
                  << ev.detR << " | tE " << l.tE << " u0 " << l.u0 << " t0 " << l.t0
                  << "\n";
    }

    // Detected-event count for this tE bin. The joint detection is the one
    // that defines "detected" here, per the taxonomy above; the per-class
    // breakdown lives in nDetClass. gg was computed for every simulated event
    // before the detection test, so nstE (the denominator) counts all draws
    // and this counts the numerator -- which is what makes EFF an efficiency.
    if (ev.detJ) l.ndtE[ev.gg] += 1.0;

    // The numerators for the other six axes, using the bin indices computed
    // for this draw before the detection test. Same definition of "detected"
    // as the tE curve -- the joint test -- so all seven efficiencies describe
    // the same thing and can be read side by side.
    if (ev.detJ) {
        l.NdMl[ev.ss] += 1.0;
        l.Ndpi[ev.qq] += 1.0;
        l.Ndu0[ev.ww] += 1.0;
        l.Ndmu[ev.vv] += 1.0;
        l.Ndmb[ev.zz] += 1.0;
        l.Ndfb[ev.pp] += 1.0;

        // Nhalo and Nself were declared, written out, and never counted, so
        // the EfLMC writer divided 0 by 0 and every row of every EfLMC file
        // ever produced ends in two "-nan" columns. Counted here, with the
        // definitions stated rather than guessed:
        //   Nhalo[1]/Nhalo[0]  fraction of detections whose LENS is a halo
        //                      star -- the population a MACHO search cares
        //                      about, and negligible toward the bulge.
        //   Nself[1]/Nself[0]  fraction that are bulge self-lensing: bulge
        //                      lens AND bulge source, the dominant channel
        //                      here and the one whose kinematics set tE.
        l.Nhalo[0] += 1.0;
        l.Nself[0] += 1.0;
        if (l.struc == GalacticComponent::HALO) l.Nhalo[1] += 1.0;
        if (l.struc == GalacticComponent::BULGE and
            s.struc == GalacticComponent::BULGE) l.Nself[1] += 1.0;
    }
}

