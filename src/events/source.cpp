// Source-star draw (func_source).
#include "events/source.h"
#include "util/random.h"

void func_source(source& s, CMD& cm, const extin& ex, int sightlineIdx) {
    int nums, num;
    double rho, rf;
    double Ds, Av = -1, Alv;
    double maxnb = 0.0;
    std::array<double, M> Map{}, Mab{}, Ai{}; 
    GalacticComponent struc;

    for (int i = 0; i < M; ++i) {
        // Reset the blended-flux accumulator for this star: Fluxb is a running sum over the source
        // (k == 1) and its unresolved neighbours, built below with `+=`, and must not carry flux
        // over from earlier draws.
        s.Fluxb[i] = 0.0;

        // Number of stars sharing this filter's seeing disc with the source.
        //
        // lambda is the mean number of field stars in the disc: a surface density times the disc
        // area. Nstart is a total number density (mass density over mean mass), complete down the
        // whole IMF. The count is 1 + Poisson(lambda): the disc is conditioned to contain the
        // source, and for a Poisson field conditioning on a point leaves the rest Poisson with the
        // same rate (Slivnyak's theorem). max(1, Poisson) would absorb the first neighbour into the
        // source and leave Roman (lambda ~ 0.14 in the 0.105" F146 disc) unblended ~99% of the
        // time. A Gaussian approximation (RandN) is also inadequate at such small lambda.
        const double lambda = std::fabs(s.Nstart * std::pow(FWHM[i] * BLEND_RADIUS_FWHM_FRAC, 2) * M_PI
                                        / (3600.0 * 3600.0));
        s.nsbl[i] = 1.0 + double(RandPois(lambda));
        if (s.nsbl[i] > maxnb) {
            maxnb = s.nsbl[i];
        }
    }

    // Positions of the unresolved neighbours (k >= 2; the source sits at the origin), drawn from their own
    // stream (rngBlend) so the main draws are untouched. Neighbour k is in band i's disc when
    // nsbl[i] >= k; it is placed uniformly in area inside the smallest disc that holds it and outside the
    // largest smaller disc that does not, so the bands' neighbour sets stay nested. Sx/Sy/NF accumulate
    // each band's neighbour light and its first moment, giving the neighbours' light centroid.
    std::array<double, M> discR{}, Sx{}, Sy{}, NF{};
    for (int i = 0; i < M; ++i) discR[i] = FWHM[i] * BLEND_RADIUS_FWHM_FRAC * ARCSEC_TO_MAS;   //[mas]

    for (int k = 1; k <= int(maxnb + 0.000000034756346); ++k) {
        double nx = 0.0, ny = 0.0;
        if (k >= 2) {
            double rHi = 1.0e30, rLo = 0.0;
            for (int i = 0; i < M; ++i) if (s.nsbl[i] >= k and discR[i] < rHi) rHi = discR[i];
            for (int i = 0; i < M; ++i) if (s.nsbl[i] <  k and discR[i] < rHi and discR[i] > rLo) rLo = discR[i];
            const double r   = std::sqrt(rLo * rLo + RandBlendUnit() * (rHi * rHi - rLo * rLo));
            const double phi = 2.0 * M_PI * RandBlendUnit();
            nx = r * std::cos(phi);
            ny = r * std::sin(phi);
        }
        do {
            nums = int(RandR(SRC_IDX_MIN, Num - SRC_IDX_END_MARGIN));
            rho  = RandR(s.Romins, s.Romaxs);
            Ds   = double(nums * step);
        } while (rho > s.Rostari[nums] or Ds < 0.0 or Ds > MaxD); //distance beyond MaxD

        rf = RandR(0.0, s.Rostar0[nums]);
        if      (rf <=  s.rho_thin[nums]) {
            struc = GalacticComponent::THIN_DISK;
        }
        else if (rf <= (s.rho_thin[nums] + s.rho_bulge[nums])) {
            struc = GalacticComponent::BULGE;
        }
        else if (rf <= (s.rho_thin[nums] + s.rho_bulge[nums] + s.rho_thick[nums])) {
            struc = GalacticComponent::THICK_DISK;
        }
        else if (rf <= (s.rho_thin[nums] + s.rho_bulge[nums] + s.rho_thick[nums] + s.rho_halo[nums])) {
            struc = GalacticComponent::HALO;
        }
        else {
           std::cerr << "Selected Galactic Component can not be initialized with rf =" << rf << ".\n";
           std::exit(EXIT_FAILURE);
       }

        if (struc == GalacticComponent::THIN_DISK) { //thin disk
            num = int(RandR(0.0, N1 - 1.0));
            for (int i = 0; i < M; ++i) {
                Mab[i] = cm.Mab_thin[num][i];
            }
            if (k == 1) {
                s.mass = cm.mass_thin[num];
                s.age  = cm.age_thin[num];
                s.logT = cm.logT_thin[num];
                s.Rstar = cm.Rstar_thin[num];
                s.cl   = cm.cl_thin[num];
                s.typ  = cm.typ_thin[num];
            }
        }

        if (struc == GalacticComponent::BULGE) {// bulge
            num = int(RandR(0.0, N2 - 1.0));
            for (int i = 0; i < M; ++i) {
                Mab[i] = cm.Mab_bulge[num][i];
            }
            if (k == 1) {
                s.mass = cm.mass_bulge[num];
                s.age  = cm.age_bulge[num];
                s.logT = cm.logT_bulge[num];
                s.Rstar = cm.Rstar_bulge[num];
                s.cl   = cm.cl_bulge[num];
                s.typ  = cm.typ_bulge[num];
            }
        }

        if (struc == GalacticComponent::THICK_DISK) { //thick disk
            num = int(RandR(0.0, N3 - 1.0));
            for (int i = 0; i < M; ++i) {
                Mab[i] = cm.Mab_thick[num][i];
            }
            if (k == 1) {
                s.mass = cm.mass_thick[num];
                s.age  = cm.age_thick[num];
                s.logT = cm.logT_thick[num];
                s.Rstar = cm.Rstar_thick[num];
                s.cl   = cm.cl_thick[num];
                s.typ  = cm.typ_thick[num];
            }
        }

        if (struc == GalacticComponent::HALO) {// stellar halo
            num = int(RandR(0.0, N4 - 1.0));
            for (int i = 0; i < M; ++i) {
                Mab[i] = cm.Mab_halo[num][i];
            }
            if (k == 1) {
                s.mass = cm.mass_halo[num];
                s.age  = cm.age_halo[num];
                s.logT = cm.logT_halo[num];
                s.Rstar = cm.Rstar_halo[num];
                s.cl   = cm.cl_halo[num];
                s.typ  = cm.typ_halo[num];
            }
        }

        Av = interpExtinctionAlongSightline(ex, sightlineIdx, Ds);
        for (int i = 0; i < M; ++i) {
            Alv = AlAv(lambda_um[i], Rv[static_cast<int>(struc)]); // A_lambda / A_V
            Ai[i] = Av * Alv + RandN(sigma[i], EXT_SCATTER_TRUNC_NSIGMA); //extinction in other bands

            if (Ai[i] < 0.0) {
                Ai[i] = 0.0;
            }
            Map[i] = Mab[i] + 5.0 * std::log10(Ds * 100.0) + Ai[i];

            if(s.nsbl[i] >= k) {
                const double fk = std::pow(10.0, -0.4 * Map[i]);
                s.Fluxb[i] += fk;
                if (k >= 2) { Sx[i] += fk * nx;  Sy[i] += fk * ny;  NF[i] += fk; }
            }
        }

        if (k == 1) {
            s.struc = struc;
            s.Ds = Ds;
            s.nums = nums;
            s.Av = Av;

            for (int i = 0; i < (M); ++i) {
                s.Ai[i]  = Ai[i];
                s.Map[i] = Map[i];
                s.Mab[i] = Mab[i];
            }
        }
    } //loop over the stars

    // The neighbours' light centroid per telescope: Rubin's reference band(s), Roman's F146.
    {
        double sx = 0.0, sy = 0.0, nf = 0.0;
        for (int band : RUBIN_REF_BANDS) { sx += Sx[band]; sy += Sy[band]; nf += NF[band]; }
        s.blendOff[0] = (nf > 0.0) ? std::array<double, 2>{sx / nf, sy / nf} : std::array<double, 2>{0.0, 0.0};
        s.blendOff[1] = (NF[6] > 0.0) ? std::array<double, 2>{Sx[6] / NF[6], Sy[6] / NF[6]}
                                      : std::array<double, 2>{0.0, 0.0};
    }

    for (int i = 0; i < M; ++i) {
        s.magb[i]  = -2.5 * std::log10(s.Fluxb[i]);
        s.blend[i] = std::pow(10.0, -0.4 * s.Map[i]) / s.Fluxb[i];

        CHECK(s.Fluxb[i] > 0.0);
        CHECK(s.nsbl[i]  >= 1.0);
        CHECK(s.blend[i] <= 1.00001);
        CHECK(s.blend[i] > 0.0);
        CHECK(std::fabs(s.Map[i] - s.Mab[i] - s.Ai[i] - 5.0 * std::log10(s.Ds * 100.0)) <= 0.1);
        CHECK(s.Ds > 0.0);
        CHECK(Av >= 0.0);
    }

    // Fisher-matrix inputs. Rubin's representative band is configurable (RUBIN_REF_BANDS,
    // config/parameters.h): the listed filters' fluxes are combined into one synthetic baseline flux
    // and one synthetic source-only flux.
    double fluxTotRubin = 0.0, fluxSrcRubin = 0.0;
    for (int band : RUBIN_REF_BANDS) {
        fluxTotRubin += s.Fluxb[band];
        fluxSrcRubin += std::pow(10.0, -0.4 * s.Map[band]);
    }
    s.mbs[0] = -2.5 * std::log10(fluxTotRubin);
    s.fb[0]  = fluxSrcRubin / fluxTotRubin;
    s.fb[1]  = s.blend[6]; // F146
    s.mbs[1] = s.magb[6];  // F146
}
