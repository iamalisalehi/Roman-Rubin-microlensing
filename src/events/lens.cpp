// Lens draw (func_lens) and the lens mass functions: Kroupa IMF, remnants, log-uniform, neutron stars, power laws.
#include "events/lens.h"
#include "util/random.h"
#include "galaxy/catalogue.h"
#include "galaxy/kinematics.h"

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Func lens  calculations                    //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void func_lens(lens & l, source & s, const extin & ex, int sightlineIdx){

    double test, tt, Am, DD;
    double mmin = Ml_min;
    double mmax = Ml_max;
    l.rhomaxl = 0.0;

    for (int k = 1; k < int(s.nums - 1) ;++k) {
        l.Dl = double(k * step);
        tt = std::sqrt((s.Ds - l.Dl) * l.Dl / s.Ds) * s.Rostar0[k];
        if (tt > l.rhomaxl)  l.rhomaxl = tt;
    }

    do {
        l.numl = int(RandR(1.0, s.nums - 1.0));
        test   =     RandR(0.0, l.rhomaxl);
        l.Dl   = double(l.numl * step);
        tt     = std::sqrt((s.Ds - l.Dl) * l.Dl / s.Ds) * s.Rostar0[l.numl];

        CHECK(tt <= l.rhomaxl);
        CHECK(l.Dl > 0.0);
        CHECK(s.Ds > 0.0);
        CHECK(l.Dl < s.Ds);

    } while(test > tt);
    l.Dl = double(l.numl * step);

    double   randflag = RandR(0.0, s.Rostar0[l.numl]);
    if      (randflag <= (s.rho_thin[l.numl])) {
        l.struc = GalacticComponent::THIN_DISK;
    }
    else if (randflag <= (s.rho_thin[l.numl] + s.rho_bulge[l.numl])) {
        l.struc = GalacticComponent::BULGE;
    }
    else if (randflag <= (s.rho_thin[l.numl] + s.rho_bulge[l.numl] + s.rho_thick[l.numl])) {
        l.struc = GalacticComponent::THICK_DISK;
    }
    else if (randflag <= (s.rho_thin[l.numl] + s.rho_bulge[l.numl] + s.rho_thick[l.numl] + s.rho_halo[l.numl])) {
        l.struc = GalacticComponent::HALO;
    }
    else {
        throw std::runtime_error("Invalid randflag: " + std::to_string(randflag) + ", rho_star0: " + std::to_string(s.Rostar0[l.numl]));
    }
     //cout<<"Dl:  "<<l.Dl<<"\t struc_lens :  "<<l.struc<<endl;

    // Which mass function this draws from is a property of the population selected by
    // --population (POPULATIONS in config/parameters.h), not of the build. Every branch that used to
    // live here -- uniform, three power laws, Kroupa+remnants -- now lives in drawLensMass()
    // beside the two new ones, so a population is one table entry rather than an `if` here
    // plus a constant there plus a filename suffix somewhere else.
    l.Ml = drawLensMass(&l.luminous);

    // The bounds bracket the masses the population can produce and also set the Mls
    // efficiency grid, so a draw outside them would land outside every efficiency bin.
    CHECK(l.Ml >= mmin * 0.999);
    CHECK(l.Ml <= mmax * 1.001);

    l.xls    = l.Dl / s.Ds;
    l.RE     = std::sqrt(4.0 * G * l.Ml * Msun * s.Ds * KP) / velocity;
    l.RE     = l.RE * std::sqrt(l.xls * (1.0 - l.xls)); //meter
    l.tetE   = l.RE / AU / l.Dl; //[mas]
    s.ros    = 1.0 * Rsun * l.xls / l.RE;
    l.pirel  = 1.0 / l.Dl - 1.0 / s.Ds; //[mas]
    l.piE    = l.pirel / l.tetE; //[]
    // ---- Deviation 74: a luminous lens's own light joins the blend ----
    // Apparent magnitudes from its main-sequence absolute magnitudes (CMD/components/lens_ml.dat),
    // distance modulus, and the dust in front of the LENS (the same tables, at Dl). No random
    // scatter is drawn, so the RNG stream is unchanged. The source fraction fb and baseline mbs of
    // both telescopes are rebuilt from the enlarged blend, and fLens records the lens's share,
    // which pulls the astrometric centroid toward the lens (lightcurve()).
    s.fLens = {0.0, 0.0};
    if (l.luminous) {
        std::array<double, 7> mab;
        if (lensAbsMag(static_cast<int>(l.struc), l.Ml, mab)) {
            const double AvL = interpExtinctionAlongSightline(ex, sightlineIdx, l.Dl);
            std::array<double, M> fluxL{};
            for (int i = 0; i < M; ++i) {
                const double AiL = std::max(0.0, AvL * AlAv(lambda_um[i], Rv[static_cast<int>(l.struc)]));
                fluxL[i] = std::pow(10.0, -0.4 * (mab[i] + 5.0 * std::log10(l.Dl * 100.0) + AiL));
                s.Fluxb[i] += fluxL[i];
                s.magb[i]   = -2.5 * std::log10(s.Fluxb[i]);
                s.blend[i]  = std::pow(10.0, -0.4 * s.Map[i]) / s.Fluxb[i];
            }
            double fluxTotRubin = 0.0, fluxSrcRubin = 0.0, fluxLensRubin = 0.0;
            for (int band : RUBIN_REF_BANDS) {
                fluxTotRubin  += s.Fluxb[band];
                fluxSrcRubin  += std::pow(10.0, -0.4 * s.Map[band]);
                fluxLensRubin += fluxL[band];
            }
            s.mbs[0]   = -2.5 * std::log10(fluxTotRubin);
            s.fb[0]    = fluxSrcRubin / fluxTotRubin;
            s.fb[1]    = s.blend[6];
            s.mbs[1]   = s.magb[6];
            s.fLens[0] = fluxLensRubin / fluxTotRubin;
            s.fLens[1] = fluxL[6] / s.Fluxb[6];
            CHECK(s.fb[0] + s.fLens[0] <= 1.000001 and s.fb[1] + s.fLens[1] <= 1.000001);
        } else {
            l.luminous = false;
        }
    }

    l.u0     = RandR(U0_MIN_DRAW, u0m);
    l.t0     = RandR(T0_MARGIN_DAYS, Tobs - T0_MARGIN_DAYS);
    l.DeltaT = std::sqrt(4.0 + l.u0 * l.u0) * l.tetE; //[mas]

    vrel(s,l);
    l.tE  = std::fabs(l.RE / (l.Vt * 1000.0 * 3600.0 * 24.0));//[days]
    l.A0  = (l.u0 * l.u0 + 2.0) / std::sqrt(l.u0 * l.u0 * (l.u0 * l.u0 + 4.0));
    l.mi1 = s.magb[2] - 2.5 * std::log10(std::fabs(l.A0 + 1.0) * 0.5 * s.blend[2] + 1.0 - s.blend[2]);
    l.mi2 = s.magb[2] - 2.5 * std::log10(std::fabs(l.A0 - 1.0) * 0.5 * s.blend[2] + 1.0 - s.blend[2]);
    Am = l.A0 + 1.0;
    DD = double(4.0 - Am * Am + Am * std::sqrt(Am * Am - 4.0)) / (Am * Am * 0.5 - 2.0);
    s.FWHM = 2.0 * l.tE * std::sqrt(std::fabs(DD - l.u0 * l.u0));//days
    // if(DD<double(l.u0*l.u0))  s.FWHM=l.tE;

    //cout<<"Ml:  "<<l.Ml<<"\t RE:  "<<l.RE/AU<<"\t tE:  "<<l.tE<<endl;
    //cout<<"ros:  "<<s.ros<<"\t murel:  "<<l.murel<<"\t tetE:  "<<l.tetE<<endl;
    //cout<<"u0:  "<<l.u0<<"\t pirel:  "<<l.pirel<<"\t piE:  "<<l.piE<<endl;
    //cout<<"mus1:  "<<s.mus1<<"\t mus2:  "<<s.mus2<<"\t mul1:  "<<l.mul1<<"\t mul2:  "<<l.mul2<<endl;
    //cout<<"DD:  "<<DD<<"\t FWHM:  "<<s.FWHM<<endl;

    CHECK(l.tE > 0.0);
    CHECK(l.Dl <= s.Ds);
    CHECK(l.Vt > 0.0);
    CHECK(l.t0 <= Tobs);
    CHECK(l.t0 >= 0.0);
    CHECK(l.Ml <= mmax);
    CHECK(l.Ml >= mmin);
    CHECK(l.pirel > 0.0);
    CHECK(l.RE >= 0.0);
    CHECK(l.xls < 1.0);
    CHECK(l.Dl >= 0.0);
    CHECK(s.Ds >= 0.0);
    CHECK(s.Ds <= MaxD);
    CHECK(l.u0 >= 0.0);
    CHECK(DD >= double(l.u0 * l.u0));
    CHECK(s.FWHM >= 0.0);
    CHECK(Am >= 2.0);

//    cout<<">>>>>>>>>>>> End of func_lens <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
}

///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//            Kroupa initial mass function + stellar remnants         //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//

// Draw an INITIAL stellar mass from the Kroupa (2001) broken power law,
//
//     dN/dM  =  k_i * M^-alpha_i
//
// with alpha = 0.3, 1.3, 2.3 on [0.01, 0.08], [0.08, 0.5], [0.5, 120] Msun.
//
// Sampled by inverse CDF, not rejection. Rejection sampling is what the legacy IMnum 2-4
// branches use, and each of them carries hardcoded acceptance bounds (CHECK(f >= 1.4) and
// friends) that are only correct for the 3-5000 Msun range they were written against --
// change the range and the CHECKs fire. Inverse CDF has no such constants: the segment
// weights and the inverse are both computed from the breaks themselves, so the sampler
// stays correct if the breaks are ever retuned.
//
// The k_i are fixed by requiring the IMF be CONTINUOUS at each break, not chosen freely:
//     k1 * 0.08^-0.3 = k2 * 0.08^-1.3  =>  k2 = k1 * 0.08^(1.3-0.3)
//     k2 * 0.50^-1.3 = k3 * 0.50^-2.3  =>  k3 = k2 * 0.50^(2.3-1.3)
// A discontinuous IMF would put a step in the mass distribution at 0.08 and 0.5 Msun, and
// therefore a step in the tE distribution, which is the observable being compared to data.
double drawKroupaInitialMass()
{
    const double lo[3] = {KROUPA_MI_MIN, KROUPA_BREAK1, KROUPA_BREAK2};
    const double hi[3] = {KROUPA_BREAK1, KROUPA_BREAK2, KROUPA_MI_MAX};
    const double al[3] = {KROUPA_ALPHA1, KROUPA_ALPHA2, KROUPA_ALPHA3};

    // Continuity coefficients, k1 fixed to 1 (overall normalisation is irrelevant here --
    // we only ever sample from this, never evaluate an absolute number density).
    double k[3];
    k[0] = 1.0;
    k[1] = k[0] * std::pow(KROUPA_BREAK1, KROUPA_ALPHA2 - KROUPA_ALPHA1);
    k[2] = k[1] * std::pow(KROUPA_BREAK2, KROUPA_ALPHA3 - KROUPA_ALPHA2);

    // Number of stars per segment: k_i * integral of M^-alpha over the segment.
    // None of the Kroupa slopes equals 1, so the (1-alpha) form never divides by zero;
    // the CHECK keeps that assumption honest if the slopes are ever changed.
    double w[3], total = 0.0;
    for (int i = 0; i < 3; ++i) {
        CHECK(std::fabs(1.0 - al[i]) > 1e-9);
        const double p = 1.0 - al[i];
        w[i] = k[i] * (std::pow(hi[i], p) - std::pow(lo[i], p)) / p;
        CHECK(w[i] > 0.0);
        total += w[i];
    }

    // Pick a segment in proportion to its star count, then invert that segment's CDF.
    double u = RandR(0.0, total);
    int seg = 2;
    for (int i = 0; i < 3; ++i) {
        if (u <= w[i]) { seg = i; break; }
        u -= w[i];
    }

    const double p = 1.0 - al[seg];
    const double v = RandR(0.0, 1.0);
    const double lop = std::pow(lo[seg], p);
    const double Mi = std::pow(lop + v * (std::pow(hi[seg], p) - lop), 1.0 / p);

    CHECK(Mi >= KROUPA_MI_MIN * 0.999);
    CHECK(Mi <= KROUPA_MI_MAX * 1.001);
    return Mi;
}


// Map an initial mass to what is still there to act as a lens ~10 Gyr later.
//
// This is where the long-tE tail comes from, and it is not a detail. A star born at 25 Msun
// is long gone, but it left a ~6 Msun black hole -- and since the Einstein radius and the
// event timescale both scale as sqrt(Ml), that black hole lenses for roughly five times as
// long as the 0.3 Msun dwarf next to it. Long events are the ones that span Roman's season
// gaps, so the remnant prescription is what populates the regime this whole project is
// about. Dropping remnants would leave the short-tE yield science intact and quietly
// remove the long-tE precision science.
double remnantMass(double initialMass)
{
    const double Mi = initialMass;

    // Still burning hydrogen: a ~10 Gyr population has a turnoff near 1 Msun, so anything
    // lighter is unevolved and lenses at its birth mass.
    if (Mi < MS_TURNOFF) return Mi;

    // White dwarf. Kalirai et al. (2008) semi-empirical initial-final mass relation,
    // Mf = 0.109*Mi + 0.394, calibrated on open-cluster white dwarfs.
    if (Mi < WD_MI_MAX) return 0.109 * Mi + 0.394;

    // Neutron star. The observed mass distribution is narrow, so a single canonical value
    // is a better model than a spread invented to look sophisticated.
    if (Mi < NS_MI_MAX) return NS_MASS;

    // Black hole. Rough proportional fallback: Ml = 0.24*Mi gives ~4.8 Msun at Mi = 20 and
    // ~28.8 at Mi = 120, spanning the observed stellar-mass black hole range. This is the
    // crudest step in the chain -- black hole remnant masses depend on metallicity and
    // mass loss in ways no single slope captures -- and is flagged in OPEN_ITEMS.md.
    return BH_MI_SLOPE * Mi;
}


// Flat in log M over [lo, hi]. The standard choice for a black-hole lens population, and it
// is a statement about ignorance rather than about stars: no mass function is measured over
// 3-1000 Msun, so weighting every decade equally is the prior that does not invent structure.
//
// Note what this does to a sample: half the draws land above sqrt(lo*hi) = 55 Msun, where
// tE ~ sqrt(Ml) makes events long. That is the regime Roman's season gaps bite hardest, and
// it is why this population is worth simulating separately rather than as the tail of a
// bulge run, where such lenses are a fraction of a per cent of the draws.
double drawLogUniformMass(double lo, double hi)
{
    CHECK(lo > 0.0);
    CHECK(hi > lo);
    const double x = RandR(std::log(lo), std::log(hi));
    const double M = std::exp(x);
    CHECK(M >= lo * 0.999);
    CHECK(M <= hi * 1.001);
    return M;
}


// A measured neutron-star mass distribution: Gaussian, truncated to the range in which
// neutron stars actually exist (Ozel & Freire 2016; constants in config/parameters.h).
//
// Redraw rather than clamp. Clamping to an edge piles probability onto 1.10 and 2.20 exactly
// -- a spike at the boundary that no physical population has, and one that would show up in
// every tE histogram as two spurious lines, since tE ~ sqrt(Ml).
double drawNeutronStarMass()
{
    for (int guard = 0; guard < 1000; ++guard) {
        const double M = NS_MEAN_MASS + RandN(NS_MASS_SIG, NS_MASS_TRUNC_NSIGMA);
        if (M >= NS_MASS_LO and M <= NS_MASS_HI) return M;
    }
    // Unreachable in practice: the truncation is +/-1.7 sigma at the tighter end, so a
    // single draw succeeds ~91% of the time and 1000 failures has probability ~1e-1000.
    // Kept so a future edit to the constants cannot produce a silent infinite loop.
    throw std::runtime_error("drawNeutronStarMass: no draw inside [NS_MASS_LO, NS_MASS_HI]");
}


// The one entry point func_lens uses. Which mass function runs is a property of the
// population selected by --population, not of the build.
double drawLensMass(bool* luminous)
{
    if (luminous) *luminous = false;
    switch (gPop->mf) {
    case MassFunction::KROUPA_REMNANTS: {
        // Draw what the star was BORN as, then ask what is left of it. The order matters:
        // the mass function describes formation, the lens is whatever survived, and
        // collapsing the two loses the black holes that make the long-tE regime exist.
        const double Mi = drawKroupaInitialMass();
        // Deviation 74: below the turnoff it is still a main-sequence star, and shines (white
        // dwarfs, neutron stars, black holes and brown dwarfs are treated as dark).
        if (luminous) *luminous = (Mi < MS_TURNOFF and Mi >= KROUPA_BREAK1);
        return remnantMass(Mi);
    }
    case MassFunction::LOG_UNIFORM:
        return drawLogUniformMass(mlMin(), mlMax());
    case MassFunction::NEUTRON_STAR:
        return drawNeutronStarMass();
    case MassFunction::UNIFORM:
        return RandR(mlMin(), mlMax());
    case MassFunction::POWER_LAW_05:
        return drawPowerLawMass(mlMin(), mlMax(), 0.5);
    case MassFunction::POWER_LAW_10:
        return drawPowerLawMass(mlMin(), mlMax(), 1.0);
    case MassFunction::POWER_LAW_20:
        return drawPowerLawMass(mlMin(), mlMax(), 2.0);
    }
    throw std::runtime_error("drawLensMass: unhandled mass function");
}


// dN/dM ~ M^-alpha on [lo, hi], by inverse CDF.
//
// The legacy code sampled these by rejection against a bounding box, which is correct but
// wasteful and, more to the point, consumes a variable number of RNG draws per event -- so
// two runs differing only in mass function diverge in their random streams for reasons that
// have nothing to do with the physics. Inverse CDF takes exactly one draw, always.
double drawPowerLawMass(double lo, double hi, double alpha)
{
    CHECK(lo > 0.0);
    CHECK(hi > lo);
    const double v = RandR(0.0, 1.0);

    // alpha == 1 is not a corner case to guard against, it is one of the legacy populations
    // (macho-m1, dN/dM ~ 1/M). There the (1-alpha) exponent form divides by zero, and the
    // CDF is logarithmic instead: M = lo * (hi/lo)^v, which is exactly log-uniform.
    if (std::fabs(1.0 - alpha) < 1e-9)
        return lo * std::pow(hi / lo, v);

    const double p = 1.0 - alpha;
    const double lop = std::pow(lo, p);
    return std::pow(lop + v * (std::pow(hi, p) - lop), 1.0 / p);
}
