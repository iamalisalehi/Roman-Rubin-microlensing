// Relative kinematics of a lens and its source (vrel).
#include "galaxy/kinematics.h"
#include "galaxy/density.h"
#include "util/random.h"

// Row of the age-velocity table for a catalogue age [Gyr]; ages of 10 Gyr and more take the last row.
static int thinRowOfAge(double age)
{
    for (int r = 0; r < THIN_NBIN - 1; ++r)
        if (age < THIN_AGE_EDGES[r + 1]) return r;
    return THIN_NBIN - 1;
}

// Galactocentric radius Rb and height zb [kpc] of the point at distance D [kpc] along the sightline, as in Disk_model.
static void sightlinePoint(const source & s, double D, double & Rb, double & zb)
{
    zb = std::sin(s.FI) * D + Z_SUN;   // height above the Galactic plane
    double yb = std::cos(s.FI) * std::sin(s.TET) * D;
    double xb = D * std::cos(s.FI) * std::cos(s.TET) - Dsun;
    Rb = std::sqrt(xb * xb + yb * yb);
}

// Row for a thin-disc star with no catalogue age at distance D along the sightline: draw its density bin
// with probability proportional to thinBinDensity there (as Disk_model builds it), take that bin's row.
static int thinRowAtPosition(const source & s, double D)
{
    double Rb, zb;
    sightlinePoint(s, D, Rb, zb);
    double w[THIN_NBIN], tot = 0.0;
    for (int ii = 0; ii < THIN_NBIN; ++ii) {
        w[ii] = thinBinDensity(ii, Rb, zb);
        tot += w[ii];
    }
    double test = RandR(0.0, tot), acc = 0.0;
    for (int ii = 0; ii < THIN_NBIN - 1; ++ii) {
        acc += w[ii];
        if (test <= acc) return ii;
    }
    return THIN_NBIN - 1;
}

// Thick-disc episode (0 young, 1 old) of a star with no catalogue age at distance D along the sightline: drawn
// with probability proportional to the two episodes' densities there.
static int thickEpisodeAtPosition(const source & s, double D)
{
    double Rb, zb;
    sightlinePoint(s, D, Rb, zb);
    const double w0 = thickEpisodeDensity(0, Rb, zb), w1 = thickEpisodeDensity(1, Rb, zb);
    return (RandR(0.0, w0 + w1) <= w0) ? 0 : 1;
}

// Draws source and lens velocities from the component dispersions (plus disk rotation), projects
// them and the Sun's motion onto the lens plane, and sets the relative speed l.Vt [km/s], the
// proper motions [mas/day] and the source-trajectory angle s.xi [rad].
void vrel(source & s, lens & l){
    if (l.Dl == 0.0) l.Dl = 0.00034735;
    double Rlc = std::sqrt(l.Dl * l.Dl * std::cos(s.FI) * std::cos(s.FI) + Dsun * Dsun - 2. * Dsun * l.Dl * std::cos(s.TET) * std::cos(s.FI));
    double Rsc = std::sqrt(s.Ds * s.Ds * std::cos(s.FI) * std::cos(s.FI) + Dsun * Dsun - 2. * Dsun * s.Ds * std::cos(s.TET) * std::cos(s.FI));
    if (Rlc == 0.0) Rlc = 0.0000000000034346123;
    if (Rsc == 0.0) Rsc = 0.000000000004762654134;

    double LVx, SVx, vt2;
    double SVT, SVR, SVZ, LVT, LVR, LVZ;
    double test;
    double  VSunx, vls2, vls1;
    double  tetd ;

    double NN = VEL_NSIGMA_TRUNC;
    double sigma_R_DiskL = 0.0,   sigma_T_DiskL = 0.0,   sigma_Z_DiskL = 0.0;   // used only for thin-disc objects
    double sigma_R_DiskS = 0.0,   sigma_T_DiskS = 0.0,   sigma_Z_DiskS = 0.0;
    double sigma_R_TDiskL = 0.0,  sigma_T_TDiskL = 0.0,  sigma_Z_TDiskL = 0.0;  // used only for thick-disc objects
    double sigma_R_TDiskS = 0.0,  sigma_T_TDiskS = 0.0,  sigma_Z_TDiskS = 0.0;
    double sigma_R_halo  = HALO_SIGMA_R, sigma_T_halo  = HALO_SIGMA_T, sigma_Z_halo  = HALO_SIGMA_Z;
    double sigma_R_Bulge = BULGE_SIGMA_R, sigma_T_Bulge = BULGE_SIGMA_T, sigma_Z_Bulge = BULGE_SIGMA_Z;

    // Thin-disc dispersions: a source takes the row of its own catalogue age; a lens (no age), or a source
    // without one, draws a density bin at its position.
    if (s.struc == GalacticComponent::THIN_DISK) {
        int row = (s.age >= 0.0) ? thinRowOfAge(s.age) : thinRowAtPosition(s, s.Ds);
        sigma_R_DiskS = THIN_SIGMA_R[row]; sigma_T_DiskS = THIN_SIGMA_T[row]; sigma_Z_DiskS = THIN_SIGMA_Z[row];
    }
    if (l.struc == GalacticComponent::THIN_DISK) {
        int row = thinRowAtPosition(s, l.Dl);
        sigma_R_DiskL = THIN_SIGMA_R[row]; sigma_T_DiskL = THIN_SIGMA_T[row]; sigma_Z_DiskL = THIN_SIGMA_Z[row];
    }
    // Thick-disc dispersions, the same way: a source by its catalogue age (10 or 12 Gyr), a lens by position.
    if (s.struc == GalacticComponent::THICK_DISK) {
        int k = (s.age >= 0.0) ? (s.age < THICK_AGE_SPLIT ? 0 : 1) : thickEpisodeAtPosition(s, s.Ds);
        sigma_R_TDiskS = THICK_SIGMA_R[k]; sigma_T_TDiskS = THICK_SIGMA_T[k]; sigma_Z_TDiskS = THICK_SIGMA_Z[k];
    }
    if (l.struc == GalacticComponent::THICK_DISK) {
        int k = thickEpisodeAtPosition(s, l.Dl);
        sigma_R_TDiskL = THICK_SIGMA_R[k]; sigma_T_TDiskL = THICK_SIGMA_T[k]; sigma_Z_TDiskL = THICK_SIGMA_Z[k];
    }

    // Source velocity.
    if (s.struc == GalacticComponent::THIN_DISK) {///Galactic disk
        SVR = RandN(sigma_R_DiskS, NN);
        SVT = RandN(sigma_T_DiskS, NN);
        SVZ = RandN(sigma_Z_DiskS, NN);
    }

    else if (s.struc == GalacticComponent::BULGE) {///Galactic bulge
        SVR = RandN(sigma_R_Bulge, NN);
        SVT = RandN(sigma_T_Bulge, NN);
        SVZ = RandN(sigma_Z_Bulge, NN);
    }

    else if (s.struc == GalacticComponent::THICK_DISK) {///thick disk
        SVR = RandN(sigma_R_TDiskS, NN);
        SVT = RandN(sigma_T_TDiskS, NN);
        SVZ = RandN(sigma_Z_TDiskS, NN);
    }

    else if (s.struc == GalacticComponent::HALO){///stellar halo
        SVR = RandN(sigma_R_halo, NN);
        SVT = RandN(sigma_T_halo, NN);
        SVZ = RandN(sigma_Z_halo, NN);
    }

    else {
        std::cerr << "Selected Galactic Component " << static_cast<int>(s.struc) << " does not exist.\n";
        std::exit(EXIT_FAILURE);
    }

    if (s.struc == GalacticComponent::THIN_DISK or s.struc == GalacticComponent::THICK_DISK) {
        SVT = SVT + vro_sun * (ROT_A * std::pow(Rsc / Dsun, ROT_SLOPE) + ROT_B);
    }

    s.vs = std::sqrt(SVR * SVR + SVT * SVT + SVZ * SVZ);

    // Lens velocity.
    if (l.struc == GalacticComponent::THIN_DISK) {///Galactic disk
        LVR = RandN(sigma_R_DiskL, NN);
        LVT = RandN(sigma_T_DiskL, NN);
        LVZ = RandN(sigma_Z_DiskL, NN);
    }

    else if (l.struc == GalacticComponent::BULGE) {///Galactic bulge
        LVR = RandN(sigma_R_Bulge, NN);
        LVT = RandN(sigma_T_Bulge, NN);
        LVZ = RandN(sigma_Z_Bulge, NN);
    }

    else if (l.struc == GalacticComponent::THICK_DISK){///thick disk
        LVR = RandN(sigma_R_TDiskL, NN);
        LVT = RandN(sigma_T_TDiskL, NN);
        LVZ = RandN(sigma_Z_TDiskL, NN);
    }

    else if (l.struc == GalacticComponent::HALO) {///stellar halo
        LVR = RandN(sigma_R_halo, NN);
        LVT = RandN(sigma_T_halo, NN);
        LVZ = RandN(sigma_Z_halo, NN);
    }

    else {
        std::cerr << "Selected Galactic Component" << static_cast<int>(s.struc) << "does not exist.\n";
        std::exit(EXIT_FAILURE);
    }

    if (l.struc == GalacticComponent::THIN_DISK or l.struc == GalacticComponent::THICK_DISK) {
        LVT = LVT + vro_sun * (ROT_A * std::pow(Rlc / Dsun, ROT_SLOPE) + ROT_B);
    }

    l.vl = std::sqrt(LVT * LVT + LVZ * LVZ + LVR * LVR);

    // Angles beta of the lens and source about the Galactic centre.

    l.betal = 0.0; l.betas = 0.0;
    tetd = s.TET;
    test = double(l.Dl * std::cos(s.FI) * std::sin(tetd) / Rlc);

    double tol = 1e-11;

    if (std::fabs(test - 1.0) < tol)      l.betal = pi/2.0;
    else if (std::fabs(test + 1.0) < tol)  l.betal = -pi/2.0;
    else                          l.betal = std::asin(test);

    test = double(s.Ds * std::cos(s.FI) * std::sin(tetd) / Rsc);

    if (std::fabs(test - 1.0) < 0.01)      l.betas =  pi/2.0;
    else if (std::fabs(test + 1.0) < 0.01) l.betas = -pi/2.0;
    else                              l.betas = std::asin(test);

    if (Dsun < std::fabs(l.Dl * std::cos(s.FI) * std::cos(tetd)) and l.betal >= 0.0) l.betal =   pi - l.betal;
    if (Dsun < std::fabs(s.Ds * std::cos(s.FI) * std::cos(tetd)) and l.betas >= 0.0) l.betas =   pi - l.betas;
    if (Dsun < std::fabs(l.Dl * std::cos(s.FI) * std::cos(tetd)) and l.betal < 0.0)  l.betal = -(pi - std::fabs(l.betal));
    if (Dsun < std::fabs(s.Ds * std::cos(s.FI) * std::cos(tetd)) and l.betas < 0.0)  l.betas = -(pi - std::fabs(l.betas));


    if (!(std::fabs(l.Dl * std::cos(s.FI) * std::sin(tetd)) <= Rlc)) {
        std::cerr << "Check failed: |Dl*cos(FI)*sin(tetd)| <= Rlc\n"
                  << "Dl: " << l.Dl << "\t Ds: " << s.Ds << '\n'
                  << "FI: " << s.FI << "\t TET: " << tetd << "\t betal: " << l.betal << '\n'
                  << "Rlc: " << Rlc << "\t Rsc: " << Rsc << "\t betas: " << l.betas << '\n'
                  << "value: " << std::fabs(l.Dl * std::cos(s.FI) * std::sin(tetd)) << '\n';
        throw std::runtime_error("Lens position exceeds Rlc");
    }
    
    if (!(std::fabs(test) <= 1.0)) {
        std::cerr << "Check failed: |test| <= 1.0\n"
                  << "test: " << test << '\n'
                  << "Dl: " << l.Dl << "\t Ds: " << s.Ds << '\n'
                  << "FI: " << s.FI << "\t TET: " << tetd << "\t betas: " << l.betas << '\n';
        throw std::runtime_error("Invalid test value");
    }
    // Angles delta between each line of sight and the local velocity frame.

    if (s.TET > pi)  tetd = s.TET - 2.0 * pi;

    l.deltal = pi - std::fabs(tetd) - std::fabs(l.betal);
    l.deltas = pi - std::fabs(tetd) - std::fabs(l.betas);

    if (l.betal < 0.0) l.deltal = -1.0 * l.deltal;
    if (l.betas < 0.0) l.deltas = -1.0 * l.deltas;

    l.deltao = pi - std::fabs(tetd);
    if (tetd < 0.0) l.deltao = -1.0 * l.deltao;

    // Project onto the lens plane.
    s.SV_n1 =+ SVR * std::sin(l.deltas) - SVT * std::cos(l.deltas);
    s.LV_n1 =+ LVR * std::sin(l.deltal) - LVT * std::cos(l.deltal);
    s.VSun_n1 =+ VSunR * std::sin(l.deltao) - VSunT * std::cos(l.deltao);

    SVx = -SVR * std::cos(l.deltas) - SVT * std::sin(l.deltas);
    LVx = -LVR * std::cos(l.deltal) - LVT * std::sin(l.deltal);
    VSunx= -VSunR * std::cos(l.deltao) - VSunT * std::sin(l.deltao);

    s.SV_n2 = -std::sin(s.FI) * (SVx) + std::cos(s.FI) * SVZ;
    s.LV_n2 = -std::sin(s.FI) * (LVx) + std::cos(s.FI) * LVZ;
    s.VSun_n2 = -std::sin(s.FI) * (VSunx) + std::cos(s.FI) * (VSunZ);


    vls1 = l.xls * s.SV_n1 - s.LV_n1 + (1.0 - l.xls) * s.VSun_n1;  ///Source - lens
    vls2 = l.xls * s.SV_n2 - s.LV_n2 + (1.0 - l.xls) * s.VSun_n2;  /// Source -lens
    l.Vt = std::sqrt(std::fabs(vls1 * vls1 + vls2 * vls2 ));


    s.mus1 = double(s.SV_n1 - s.VSun_n1) * 1000.0 * 3600.0 * 24.0 / (s.Ds * AU);//[mas/days]
    s.mus2 = double(s.SV_n2 - s.VSun_n2) * 1000.0 * 3600.0 * 24.0 / (s.Ds * AU);//[mas/days]
    l.mul1 = double(s.LV_n1 - s.VSun_n1) * 1000.0 * 3600.0 * 24.0 / (l.Dl * AU);//[mas/days]
    l.mul2 = double(s.LV_n2 - s.VSun_n2) * 1000.0 * 3600.0 * 24.0 / (l.Dl * AU);//[mas/days]
    s.mus  = std::sqrt(s.mus1 * s.mus1 + s.mus2 * s.mus2);//constant
    l.mul  = std::sqrt(l.mul1 * l.mul1 + l.mul2 * l.mul2);

    l.murel = std::sqrt((s.mus1 - l.mul1) * (s.mus1 - l.mul1) + (s.mus2 - l.mul2) * (s.mus2 - l.mul2));//[mas/days]
    vt2 = l.murel * l.Dl * AU / (1000.0 * 24.0 * 3600.0);//[km/s]



    if (vls1 > 0.0 and vls2 >= 0.0)          s.xi  = std::atan(std::fabs(vls2) / std::fabs(vls1));//in radian
    else if (vls1 <= 0.0 and vls2 >  0.0)    s.xi  = std::atan(std::fabs(vls1) / std::fabs(vls2)) + M_PI / 2.0;
    else if (vls1 <  0.0 and vls2 <= 0.0)    s.xi  = std::atan(std::fabs(vls2) / std::fabs(vls1)) + M_PI;
    else if (vls1 >= 0.0 and vls2 <  0.0)    s.xi  = std::atan(std::fabs(vls1) / std::fabs(vls2)) + 3.0 * M_PI / 2.0;
    else if (s.xi >  2.0 * M_PI)             s.xi -= 2.0 * M_PI;

    CHECK(!(vls1 == 0.0 && vls2 == 0.0 && std::fabs(vt2 - l.Vt) > 0.1));

    CHECK(l.Vt >= 0.0);
    CHECK(l.Vt <= 1.0e6);
}
