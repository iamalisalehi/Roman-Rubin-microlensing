// The detection test and the detection taxonomy / efficiency tallies.
#include "sim/detect.h"
#include "fisher/fisher.h"

Detection detectEvent(SimContext& ctx, const LightCurveStats& lc) {
    const RunConfig& cfg = ctx.cfg;
    source& s = ctx.s;
    Detection d;

    // The lensing statistic is signed: chi3 - chi1 > 0 means the lensing model fits better, the only
    // direction that counts as evidence of lensing. Signing also makes the joint statistic exactly the
    // sum of the per-survey ones (dchiL = dchiL_L + dchiL_R), which makes detection monotone.
    // dchiP and dchiA are not thresholded; they are signed so that a negative value (the simpler
    // model fitting better) stays visible.
    d.dchiL   = lc.chi3  - lc.chi1;             //lensing_effect (signed)
    d.dchiP   = lc.chi2  - lc.chi1;             //parallax_effect (signed)
    d.dchiA   = lc.chi2a - lc.chi1a;            //deflection_effect (signed)
    d.dchiL_L = lc.chi3_L - lc.chi1_L;
    d.dchiP_L = (lc.chi2_L  - lc.chi1_L);
    d.dchiA_L = (lc.chi2a_L - lc.chi1a_L);
    d.dchiL_R = lc.chi3_R - lc.chi1_R;
    d.dchiP_R = (lc.chi2_R  - lc.chi1_R);
    d.dchiA_R = (lc.chi2a_R - lc.chi1a_R);

    // Three independent tests. detL/detR use each instrument's own epoch count and run-test result, so
    // Roman's dense epochs cannot raise the bar for a Rubin-driven signal (and vice versa). detJ needs
    // a persistent run in either instrument's own cadence; an interleaved run across two cadences has
    // no meaning. The union of the three gates FisherM, so an event Rubin alone detects is never
    // dropped from characterisation.
    if (s.FWHM < Tobs and d.dchiL_L > cfg.dchiDet and lc.flag_det_L > 0 and lc.ndw_L > 10) d.detL = 1;
    if (s.FWHM < Tobs and d.dchiL_R > cfg.dchiDet and lc.flag_det_R > 0 and lc.ndw_R > 10) d.detR = 1;
    if (s.FWHM < Tobs and d.dchiL   > cfg.dchiDet and (lc.flag_det_L > 0 or lc.flag_det_R > 0) and lc.ndw > 10) d.detJ = 1;

    // Monotonicity: adding data cannot destroy information, so a single-survey detection implies a
    // joint one. With the signed statistic and one fixed bar this holds by construction (chi1 and chi3
    // are accumulated over both instruments, so they equal the per-survey sums, and the auxiliary gates
    // agree). detJ_raw and the patch below are kept so the violation rate stays measurable:
    // a non-zero DET_ANOMALY count means the construction has been broken.
    d.detJ_raw = d.detJ;
    if (d.detL or d.detR) d.detJ = 1;
    return d;
}

void tallyDetection(SimContext& ctx, SightlineState& st, const EfficiencyBins& bins, const LightCurveStats& lc,
                    Detection& det) {
    const RunConfig& cfg = ctx.cfg;
    source& s = ctx.s;
    lens& l = ctx.l;
    RunTotals& run = ctx.run;

    // Detection taxonomy (DetClass, include/fisher/fisher.h): classes record which telescopes also
    // detect the event unaided. DET_ANOMALY is the case that should be impossible, a single-telescope
    // detection the joint test misses.
    det.dclsEvent = detClass(det.detL, det.detR, det.detJ);
    const int dcls = det.dclsEvent;
    st.nDetClass[dcls]    += 1;
    run.NDetClassTot[dcls] += 1;
    run.NDetClassTE[bins.gg][dcls] += 1;
    // The joint delta-chi2 is accumulated separately from the two per-survey ones; it must equal
    // their sum on every event.
    if (std::fabs((det.dchiL_L + det.dchiL_R) - det.dchiL) > 1e-6 * std::max(1.0, std::fabs(det.dchiL))) {
        if (run.nDchiMismatch < 10)
            std::cerr << "DCHI_MISMATCH sightline " << st.iScan << " joint " << det.dchiL
                      << " L+R " << det.dchiL_L + det.dchiL_R << "\n";
        run.nDchiMismatch += 1;
    }
    if (!det.detJ_raw and (det.detL or det.detR)) {
        st.nDetClass[DET_ANOMALY]    += 1;
        run.NDetClassTot[DET_ANOMALY] += 1;
        // Log what is needed to explain it. Either per-survey term can be negative through noise, so
        // a single-survey detection just over the bar can leave the sum just under it; "sum - joint"
        // must be ~0 if the bookkeeping is right.
        std::cerr << std::setprecision(10)
                  << "DET_ANOMALY_DETAIL sightline " << st.iScan << " lon " << s.lon
                  << " lat " << s.lat << " | dchiL_L " << det.dchiL_L << " dchiL_R "
                  << det.dchiL_R << " sum " << det.dchiL_L + det.dchiL_R << " joint " << det.dchiL
                  << " (sum - joint " << (det.dchiL_L + det.dchiL_R) - det.dchiL << ") bar "
                  << cfg.dchiDet << " | ndw_L " << lc.ndw_L << " ndw_R " << lc.ndw_R
                  << " ndw " << lc.ndw << " | flag_det_L " << lc.flag_det_L
                  << " flag_det_R " << lc.flag_det_R << " | detL " << det.detL << " detR "
                  << det.detR << " | tE " << l.tE << " u0 " << l.u0 << " t0 " << l.t0
                  << "\n";
    }

    // Numerator of the tE efficiency: the joint detection defines "detected" here. bins.gg was
    // computed for every draw, so nstE (the denominator) counts all of them.
    if (det.detJ) l.ndtE[bins.gg] += 1.0;

    // Numerators for the other six axes, with the same definition of "detected".
    if (det.detJ) {
        l.NdMl[bins.ss] += 1.0;
        l.Ndpi[bins.qq] += 1.0;
        l.Ndu0[bins.ww] += 1.0;
        l.Ndmu[bins.vv] += 1.0;
        l.Ndmb[bins.zz] += 1.0;
        l.Ndfb[bins.pp] += 1.0;

        // Nhalo[1]/Nhalo[0]: fraction of detections whose lens is a halo star.
        // Nself[1]/Nself[0]: fraction that are bulge self-lensing (bulge lens and bulge source).
        l.Nhalo[0] += 1.0;
        l.Nself[0] += 1.0;
        if (l.struc == GalacticComponent::HALO) l.Nhalo[1] += 1.0;
        if (l.struc == GalacticComponent::BULGE and
            s.struc == GalacticComponent::BULGE) l.Nself[1] += 1.0;
    }
}

