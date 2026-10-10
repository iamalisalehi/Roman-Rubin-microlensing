// Lens draw (func_lens) and the lens mass functions: Kroupa IMF, remnants, log-uniform, neutron stars, power laws.
#include "events/lens.h"
#include "util/random.h"
#include "galaxy/catalogue.h"
#include "galaxy/kinematics.h"

// Draws a lens along the sightline: distance, component, mass, kinematics, and (for a luminous
// lens) its own light added to the blend.
void func_lens(lens & l, source & s, const CMD & cm, const extin & ex, int sightlineIdx){

    double test, tt, Am, DD;
    double mmin = mlMin();
    double mmax = mlMax();
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

    // The mass function is set by --population (POPULATIONS in config/parameters.h); see drawLensMass().
    std::array<double, 7> mab{};
    bool mabFromCatalogue = false;
    if (gPop->mf == MassFunction::BESANCON_CATALOGUE) {
        l.Ml = drawCatalogueLens(cm, l.struc, &l.luminous, mab);
        mabFromCatalogue = true;
    } else {
        l.Ml = drawLensMass(&l.luminous);
    }

    // The bounds also set the Mls efficiency grid; a draw outside them would fall outside every bin.
    CHECK(l.Ml >= mmin * 0.999);
    CHECK(l.Ml <= mmax * 1.001);

    l.xls    = l.Dl / s.Ds;
    l.RE     = std::sqrt(4.0 * G * l.Ml * Msun * s.Ds * KP) / velocity;
    l.RE     = l.RE * std::sqrt(l.xls * (1.0 - l.xls)); //meter
    l.tetE   = l.RE / AU / l.Dl; //[mas]
    // Finite-source size. theta* = Rstar / Ds in tetE's units (1 AU at 1 kpc = 1 mas), so rho =
    // theta* / thetaE = Rstar * Rsun * xls / RE: the source radius projected onto the lens plane, in
    // Einstein radii. The two forms are the same number.
    s.thetaStar = s.Rstar * Rsun / AU / s.Ds; //[mas]
    s.rho       = s.thetaStar / l.tetE;       //[]
    l.pirel  = 1.0 / l.Dl - 1.0 / s.Ds; //[mas]
    l.piE    = l.pirel / l.tetE; //[]
    // A luminous lens's own light joins the blend. Apparent magnitudes come from its main-sequence
    // absolute magnitudes (CMD/components/lens_ml.dat), the distance modulus, and the dust in front
    // of the LENS (same tables, at Dl). No random scatter is drawn. fb and mbs of both telescopes
    // are rebuilt from the enlarged blend; fLens records the lens's share, which pulls the
    // astrometric centroid toward the lens (lightcurve()).
    s.fLens = {0.0, 0.0};
    if (l.luminous) {
        if (mabFromCatalogue or lensAbsMag(static_cast<int>(l.struc), l.Ml, mab)) {
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
}

// Kroupa initial mass function and stellar remnants.

// Draw an INITIAL stellar mass from the Kroupa (2001) broken power law dN/dM = k_i M^-alpha_i,
// alpha = 0.3, 1.3, 2.3 on [0.01, 0.08], [0.08, 0.5], [0.5, 120] Msun, by inverse CDF (the
// segment weights and the inverse follow from the breaks, so retuning them stays correct).
//
// The k_i are fixed by continuity at each break:
//     k2 = k1 * 0.08^(1.3-0.3),   k3 = k2 * 0.50^(2.3-1.3)
double drawKroupaInitialMass()
{
    const double lo[3] = {KROUPA_MI_MIN, KROUPA_BREAK1, KROUPA_BREAK2};
    const double hi[3] = {KROUPA_BREAK1, KROUPA_BREAK2, KROUPA_MI_MAX};
    const double al[3] = {KROUPA_ALPHA1, KROUPA_ALPHA2, KROUPA_ALPHA3};

    // Continuity coefficients, k1 = 1 (only sampled, never evaluated as a number density).
    double k[3];
    k[0] = 1.0;
    k[1] = k[0] * std::pow(KROUPA_BREAK1, KROUPA_ALPHA2 - KROUPA_ALPHA1);
    k[2] = k[1] * std::pow(KROUPA_BREAK2, KROUPA_ALPHA3 - KROUPA_ALPHA2);

    // Number of stars per segment: k_i * integral of M^-alpha. The CHECK guards the (1-alpha)
    // division should a slope ever equal 1.
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
// This sets the long-tE tail: tE scales as sqrt(Ml), so a ~6 Msun black hole left by a 25 Msun star
// lenses ~5x longer than a 0.3 Msun dwarf, and long events are the ones that span Roman's season gaps.
double remnantMass(double initialMass)
{
    const double Mi = initialMass;

    // Still burning hydrogen: a ~10 Gyr population has a turnoff near 1 Msun, so anything
    // lighter is unevolved and lenses at its birth mass.
    if (Mi < MS_TURNOFF) return Mi;

    // White dwarf. Kalirai et al. (2008) semi-empirical initial-final mass relation,
    // Mf = 0.109*Mi + 0.394, calibrated on open-cluster white dwarfs.
    if (Mi < WD_MI_MAX) return 0.109 * Mi + 0.394;

    // Neutron star: the observed mass distribution is narrow, so a single canonical value is used.
    if (Mi < NS_MI_MAX) return NS_MASS;

    // Black hole: rough proportional fallback, Ml = 0.24*Mi (~4.8 Msun at Mi = 20, ~28.8 at 120).
    // This is the crudest step in the chain; remnant masses depend on metallicity and mass loss.
    return BH_MI_SLOPE * Mi;
}


// Flat in log M over [lo, hi]: an uninformative prior for a black-hole population, since no mass
// function is measured over 3-1000 Msun. Half the draws land above sqrt(lo*hi) = 55 Msun, where
// tE is long and Roman's season gaps matter most.
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
// Redraw rather than clamp: clamping would pile probability onto the two edges and show up as
// spurious lines in tE histograms.
double drawNeutronStarMass()
{
    for (int guard = 0; guard < 1000; ++guard) {
        const double M = NS_MEAN_MASS + RandN(NS_MASS_SIG, NS_MASS_TRUNC_NSIGMA);
        if (M >= NS_MASS_LO and M <= NS_MASS_HI) return M;
    }
    // Unreachable in practice (a single draw succeeds ~91% of the time); guards against an infinite
    // loop if the constants are edited.
    throw std::runtime_error("drawNeutronStarMass: no draw inside [NS_MASS_LO, NS_MASS_HI]");
}


// A lens that is a random member of its Galactic component's Besancon list (population "besancon").
//
// The four lists in CMD/components/ are the complete present-day Besancon population of each
// component, the same population Nstart = rho/<m> counts (<m> = MEANMASS_*), so the lens mass
// function is consistent with the density normalisation. The lens takes the row's mass and its
// absolute magnitudes (ugrizy, F146); a row with DARK_MAG in every band (a white dwarf) is dark.
// Unlike Kroupa+remnants there are no brown dwarfs (lists start at 0.073 Msun in the thin disc,
// ~0.155 in the others), no neutron stars or black holes, and no bulge white dwarfs; the bh and
// ns populations cover those.
//
// One RNG call. (func_source's own pick, RandR(0, N-1), can never reach the last row; this one can.)
double drawCatalogueLens(const CMD& cm, GalacticComponent comp, bool* luminous, std::array<double, 7>& mab)
{
    const std::vector<double>* mass = nullptr;
    const std::vector<std::array<double, M>>* Mab = nullptr;
    switch (comp) {
    case GalacticComponent::THIN_DISK:  mass = &cm.mass_thin;  Mab = &cm.Mab_thin;  break;
    case GalacticComponent::BULGE:      mass = &cm.mass_bulge; Mab = &cm.Mab_bulge; break;
    case GalacticComponent::THICK_DISK: mass = &cm.mass_thick; Mab = &cm.Mab_thick; break;
    case GalacticComponent::HALO:       mass = &cm.mass_halo;  Mab = &cm.Mab_halo;  break;
    }
    CHECK(mass != nullptr);
    const int n = static_cast<int>(mass->size());
    int j = int(RandR(0.0, double(n)));
    if (j >= n) j = n - 1;
    for (int i = 0; i < 7; ++i) mab[i] = (*Mab)[j][i];
    if (luminous) *luminous = ((*Mab)[j][2] != DARK_MAG);
    return (*mass)[j];
}


// The entry point func_lens uses; the mass function is set by --population.
double drawLensMass(bool* luminous)
{
    if (luminous) *luminous = false;
    switch (gPop->mf) {
    case MassFunction::KROUPA_REMNANTS: {
        // Draw the birth mass, then ask what is left of it: the mass function describes formation.
        const double Mi = drawKroupaInitialMass();
        // Below the turnoff it is still a main-sequence star and shines (white dwarfs, neutron
        // stars, black holes and brown dwarfs are treated as dark).
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
    case MassFunction::BESANCON_CATALOGUE:
        throw std::runtime_error("drawLensMass: the Besancon-catalogue population is drawn in func_lens from the catalogue");
    }
    throw std::runtime_error("drawLensMass: unhandled mass function");
}


// dN/dM ~ M^-alpha on [lo, hi], by inverse CDF.
//
// Inverse CDF takes exactly one RNG draw, so random streams stay aligned between runs that differ
// only in mass function.
double drawPowerLawMass(double lo, double hi, double alpha)
{
    CHECK(lo > 0.0);
    CHECK(hi > lo);
    const double v = RandR(0.0, 1.0);

    // alpha == 1 (dN/dM ~ 1/M) makes the (1-alpha) form divide by zero; the CDF is then
    // logarithmic, M = lo * (hi/lo)^v, i.e. log-uniform.
    if (std::fabs(1.0 - alpha) < 1e-9)
        return lo * std::pow(hi / lo, v);

    const double p = 1.0 - alpha;
    const double lop = std::pow(lo, p);
    return std::pow(lop + v * (std::pow(hi, p) - lop), 1.0 / p);
}
