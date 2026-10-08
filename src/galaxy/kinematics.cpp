// Relative kinematics of a lens and its source (vrel).
#include "galaxy/kinematics.h"
#include "util/random.h"

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
    double test, age = -1;
    double  VSunx, vls2, vls1;
    double  tetd ;

    double NN = VEL_NSIGMA_TRUNC;
    double sigma_R_Disk,          sigma_T_Disk,          sigma_Z_Disk;
    double sigma_R_DiskL,         sigma_T_DiskL,         sigma_Z_DiskL;
    double sigma_R_DiskS,         sigma_T_DiskS,         sigma_Z_DiskS;
    double sigma_R_TDisk = THICK_SIGMA_R,  sigma_T_TDisk = THICK_SIGMA_T,  sigma_Z_TDisk = THICK_SIGMA_Z;
    double sigma_R_halo  = HALO_SIGMA_R, sigma_T_halo  = HALO_SIGMA_T, sigma_Z_halo  = HALO_SIGMA_Z;
    double sigma_R_Bulge = BULGE_SIGMA_R, sigma_T_Bulge = BULGE_SIGMA_T, sigma_Z_Bulge = BULGE_SIGMA_Z;

    double Rho[4] = {0.0};
    double maxr = 0.0;
    for (int i = 0; i < 4; ++i) {
        Rho[i] = rho0[i] * corr[i] / d0[i];
        maxr += Rho[i];
    }


    for (int i = 0; i < 2; ++i) {
        test = RandR(0.0, maxr); ///total ages

        if      (test <=  Rho[0])                       {sigma_R_Disk = THIN_SIGMA_R[0]; sigma_T_Disk = THIN_SIGMA_T[0]; sigma_Z_Disk = THIN_SIGMA_Z[0]; age = THIN_AGE[0];}
        else if (test <= (Rho[0]+Rho[1]))               {sigma_R_Disk = THIN_SIGMA_R[1]; sigma_T_Disk = THIN_SIGMA_T[1]; sigma_Z_Disk = THIN_SIGMA_Z[1]; age = THIN_AGE[1];}
        else if (test <= (Rho[0]+Rho[1]+Rho[2]))        {sigma_R_Disk = THIN_SIGMA_R[2]; sigma_T_Disk = THIN_SIGMA_T[2]; sigma_Z_Disk = THIN_SIGMA_Z[2]; age = THIN_AGE[2];}
        else if (test <= (Rho[0]+Rho[1]+Rho[2]+Rho[3])) {sigma_R_Disk = THIN_SIGMA_R[3]; sigma_T_Disk = THIN_SIGMA_T[3]; sigma_Z_Disk = THIN_SIGMA_Z[3]; age = THIN_AGE[3];}

        if (i == 0) {
            sigma_R_DiskS = sigma_R_Disk;
            sigma_T_DiskS = sigma_T_Disk;
            sigma_Z_DiskS = sigma_Z_Disk;
        }
        if (i == 1) {
            sigma_R_DiskL = sigma_R_Disk;
            sigma_T_DiskL = sigma_T_Disk;
            sigma_Z_DiskL = sigma_Z_Disk;
        }
    }
    CHECK(age >= 0);


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
        SVR = RandN(sigma_R_TDisk, NN);
        SVT = RandN(sigma_T_TDisk, NN);
        SVZ = RandN(sigma_Z_TDisk, NN);
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
        LVR = RandN(sigma_R_TDisk, NN);
        LVT = RandN(sigma_T_TDisk, NN);
        LVZ = RandN(sigma_Z_TDisk, NN);
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
