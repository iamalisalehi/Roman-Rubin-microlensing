// Source-star draw (func_source).
#include "events/source.h"
#include "util/random.h"


///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Func source calculations                   //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void func_source(source& s, CMD& cm, const extin& ex, int sightlineIdx) {
    int nums, num;
    double rho, rf;
    double Ds, Av = -1, Alv;
    double maxnb = 0.0;
    std::array<double, M> Map{}, Mab{}, Ai{}; 
    GalacticComponent struc;

    for (int i = 0; i < M; ++i) {
        // Reset the blended-flux accumulator for THIS star. Fluxb is a running sum over the
        // source (k == 1) and its unresolved neighbours (k = 2..nsbl), built by the k-loop
        // below with `+=`. Without this line it is never cleared, so every star inherits the
        // blend flux of every star drawn before it: magb brightens without limit as
        // magb(N) ~ magb(1) - 2.5*log10(N), and blend = F_source/Fluxb decays as 1/N until no
        // event can pass the acceptance test and the Monte Carlo silently stops detecting
        // anything at all. Measured on a 300k-draw run before this fix: the blend fraction
        // fell from 0.130 on the first star to 0.0017 on the second, and the blended baseline
        // magnitude drifted from 17.6 to 3.5 -- brighter than Sirius, and past the saturation
        // guard. Present since the initial commit, inherited from the legacy LMC code.
        s.Fluxb[i] = 0.0;

        // Number of stars sharing this filter's seeing disc with the source.
        //
        // lambda is the MEAN number of field stars in the disc: a surface density times the
        // disc area. Nstart is a total number density -- mass density divided by the mean
        // mass of the relevant population -- so it is complete down the whole IMF and does
        // NOT stop at any survey's detection limit. The geometry is right too. What was
        // wrong was the statistics.
        //
        // The old code added RandN(sqrt(lambda), 2.0) and clamped the result up to 1. For
        // Rubin, where lambda ~ 12.7 in r band, that is a passable Gaussian approximation to
        // a Poisson count. For Roman it is qualitatively wrong: lambda ~ 0.14 in the 0.105"
        // F146 disc, the truncated Gaussian can add at most 2*sqrt(0.14) = 0.75, and the
        // clamp then rounds every draw to exactly 1 -- the source alone. blend[6] came out
        // as exactly 1.0 for all 1,950 events of the 2026-08-30 run: not a narrow
        // distribution but a constant, and Roman could not be blended at all, by construction.
        //
        // A Poisson draw is the correct count and fixes it without special-casing anything:
        // at lambda = 0.14, P(at least one neighbour) = 1 - exp(-0.14) = 13.3%, so Roman is
        // mostly unblended -- which is physically right for a 0.105" PSF -- but not always,
        // which is the part that was missing.
        //
        // Why 1 + Poisson(lambda) rather than max(1, Poisson(lambda)): we are looking AT a
        // source, so the disc is conditioned to contain it. For a Poisson field, conditioning
        // on a point at a given location leaves the rest of the field Poisson with the same
        // rate (Slivnyak's theorem), so the total is the source plus Poisson(lambda)
        // neighbours. max(1, .) would instead absorb the first neighbour into the source and
        // keep Roman unblended 99% of the time -- the same bug in a new place.
        //
        // This moves Rubin slightly as well: its mean count goes from lambda to 1 + lambda,
        // about +8% in r band. That is the same Palm-conditioning correction, previously
        // missing, and it is not separable from the Roman fix -- they are one line.
        const double lambda = std::fabs(s.Nstart * std::pow(FWHM[i] * BLEND_RADIUS_FWHM_FRAC, 2) * M_PI
                                        / (3600.0 * 3600.0));
        s.nsbl[i] = 1.0 + double(RandPois(lambda));
        if (s.nsbl[i] > maxnb) {
            maxnb = s.nsbl[i];
        }
        //cout << "nsbl[i]: " << s.nsbl[i] << "\t maxnb: " << maxnb << endl;
    }

    for (int k = 1; k <= int(maxnb + 0.000000034756346); ++k) {
        do {
            nums = int(RandR(SRC_IDX_MIN, Num - SRC_IDX_END_MARGIN));
            rho  = RandR(s.Romins, s.Romaxs);
            Ds   = double(nums * step);
        } while (rho > s.Rostari[nums] or Ds < 0.0 or Ds > MaxD); //distance larger than 20.0

//        cout<<"k:  "<<k<<"\t Ds:  "<<Ds<<"\t nums:  "<<nums<<endl;

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
       }  // cout<<"struc :  "<<struc<<endl;

        if (struc == GalacticComponent::THIN_DISK) { //thin disk
            num = int(RandR(0.0, N1 - 1.0));
            for (int i = 0; i < M; ++i) {
                Mab[i] = cm.Mab_thin[num][i];
            }
            if (k == 1) {
                s.mass = cm.mass_thin[num];
                s.age  = cm.age_thin[num];
                s.logT = cm.logT_thin[num];
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
                s.cl   = cm.cl_halo[num];
                s.typ  = cm.typ_halo[num];
            }
        }

        Av = interpExtinctionAlongSightline(ex, sightlineIdx, Ds);
//        cout << "Av: " << Av << endl;
        for (int i = 0; i < M; ++i) {
            Alv = AlAv(lambda_um[i], Rv[static_cast<int>(struc)]); // A_lambda / A_V
            Ai[i] = Av * Alv + RandN(sigma[i], EXT_SCATTER_TRUNC_NSIGMA); //extinction in other bands

            if (Ai[i] < 0.0) {
                Ai[i] = 0.0;
            }
//        cout << "Alv: " << Alv << "\t filter: " << i << "\tRv[struc]: " << Rv[static_cast<int>(struc)]
//             << "\tstruc: " << static_cast<int>(struc) << endl;
//        cout << "Av: " << Av << "\t Ai[i]: " << Ai[i] << endl;
            Map[i] = Mab[i] + 5.0 * std::log10(Ds * 100.0) + Ai[i];

            if(s.nsbl[i] >= k) {
                s.Fluxb[i] += std::pow(10.0, -0.4 * Map[i]);
            }
//        cout << "filter:  " << i << "\tAlv: " << Alv << "\tExt: " << Ai[i]
//             << "\tMab: " << Mab[i] << "\tMap: " << Map[i] << endl;
        }
//        cout << "*****************************" << endl;

        if (k == 1) {
//            cout << "mass:  " << s.mass << "\t age:  " << s.age << "\t type:  " << s.typ  << endl;
//            cout << "cl:  "   << s.cl   << "\t age:  " << s.age << "\t tef:  "  << s.logT << endl;
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
        //cout<<"=========================="<<endl;
    } //loop over the stars

    for (int i = 0; i < M; ++i) {
        s.magb[i]  = -2.5 * std::log10(s.Fluxb[i]);
        s.blend[i] = std::pow(10.0, -0.4 * s.Map[i]) / s.Fluxb[i];

        CHECK(s.Fluxb[i] > 0.0);
        CHECK(s.nsbl[i]  >= 1.0);
        CHECK(s.blend[i] <= 1.00001);
        CHECK(s.blend[i] > 0.0);
//        cout << "i = " << i << "\tnsbl[i]: " << s.nsbl[i] << "\tblend[i]: " << s.blend[i] << endl;
//        if (s.nsbl[i] == 1) {
//            CHECK(s.blend[i] >= 1.0);
//        }
        CHECK(std::fabs(s.Map[i] - s.Mab[i] - s.Ai[i] - 5.0 * std::log10(s.Ds * 100.0)) <= 0.1);
        CHECK(s.Ds > 0.0);
        CHECK(Av >= 0.0);
    }

    // For purposes of Fisher Matrix calculations.
    // Rubin's representative band is configurable (RUBIN_REF_BANDS, config/parameters.h) rather than
    // hardcoded to r: combine the listed filters' fluxes into one synthetic baseline flux
    // and one synthetic source-only flux, the same way a single filter already worked.
    // RUBIN_REF_BANDS = {2} reduces to exactly the old r-only behavior, bit-for-bit.
    double fluxTotRubin = 0.0, fluxSrcRubin = 0.0;
    for (int band : RUBIN_REF_BANDS) {
        fluxTotRubin += s.Fluxb[band];
        fluxSrcRubin += std::pow(10.0, -0.4 * s.Map[band]);
    }
    s.mbs[0] = -2.5 * std::log10(fluxTotRubin);
    s.fb[0]  = fluxSrcRubin / fluxTotRubin;
    s.fb[1]  = s.blend[6]; // F146
    s.mbs[1] = s.magb[6];  // F146
//    cout << "Ds:   " << s.Ds << "\t nums: " << s.nums << endl;
//    debug note: enabling the line above showed that func_source runs many many times,
//    but no star is detectable by criteria stated in main function.
//cout<<">>>>>>>>>>>>  End of Func_source <<<<<<<<<<<<<<<<<<<<<"<<endl;
}
