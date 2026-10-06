// Reading the CMD component files and the lens mass-luminosity table.
#include "galaxy/catalogue.h"

///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Read CMD                                   //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void read_cmd(CMD & cm)
{
// mass  logT  Mbol  Age  Pop  Roman_F146  LSST_u  LSST_g  LSST_r  LSST_i  LSST_z  LSST_y  CL  Type //20/07/2026
// 0     1     2     3    4    5           6       7       8       9       10      11      12  13
    double Mbol, Pop;
    double dummy;
    std::string header;
//    double logL, gravity, metal, B, V, R, I, J, H;

    // ================================ THIN DISK =============================
    std::ifstream fp2(PATH_CMD_THIN);
    if (!fp2.is_open()) {
        throw std::runtime_error(std::string("cannot read ") + PATH_CMD_THIN);
    }

    std::getline(fp2, header);   // Skip the header line
    for (size_t j = 0; j < N1; ++j) {
          if (!(fp2 >> cm.mass_thin[j]
                    >> cm.logT_thin[j]
                    >> Mbol
                    >> cm.age_thin[j]
                    >> Pop
                    >> cm.Mab_thin[j][6]
                    >> cm.Mab_thin[j][0]
                    >> cm.Mab_thin[j][1]
                    >> cm.Mab_thin[j][2]
                    >> cm.Mab_thin[j][3]
                    >> cm.Mab_thin[j][4]
                    >> cm.Mab_thin[j][5]
                    >> cm.cl_thin[j]
                    >> cm.typ_thin[j])) {
              throw std::runtime_error("Unexpected end of thin_disk.dat");
          }
        // The CHECKs follow the read: placed before it they tested the zero-initialised slot, never the value read.
        // M_r <= 30: the thin disc's luminous brown-dwarf-limit stars reach M_r = 28.86 (bos10, Deviation 88);
        // Typ <= 9.2: white dwarfs are Typ 9.0-9.2 and sit in the list as dark entries.
        CHECK(cm.mass_thin[j]   >= 0.0);
        CHECK(cm.logT_thin[j]   >= 0.0);
        CHECK(cm.Mab_thin[j][2] <= 30.0 or cm.Mab_thin[j][2] == DARK_MAG);
        CHECK(cm.age_thin[j]    <= 10);
        CHECK(cm.cl_thin[j]     <= 7);
        CHECK(cm.typ_thin[j]    <= 9.2);
      }

    // Make sure there's no extra data.
    CHECK(!(fp2 >> dummy));
    fp2.close();

    // ================================ BULGE ==================================
    fp2.open(PATH_CMD_BULGE);
    if (!fp2.is_open()) {
        throw std::runtime_error(std::string("cannot read ") + PATH_CMD_BULGE);
    }

    std::getline(fp2, header);   // Skip the header line
    for (size_t j = 0; j < N2; ++j) {
          if (!(fp2 >> cm.mass_bulge[j]
                    >> cm.logT_bulge[j]
                    >> Mbol
                    >> cm.age_bulge[j]
                    >> Pop
                    >> cm.Mab_bulge[j][6]
                    >> cm.Mab_bulge[j][0]
                    >> cm.Mab_bulge[j][1]
                    >> cm.Mab_bulge[j][2]
                    >> cm.Mab_bulge[j][3]
                    >> cm.Mab_bulge[j][4]
                    >> cm.Mab_bulge[j][5]
                    >> cm.cl_bulge[j]
                    >> cm.typ_bulge[j])) {
              throw std::runtime_error("Unexpected end of bulge.dat");
          }
        // CHECKs follow the read (see the thin disk); same bounds.
        CHECK(cm.mass_bulge[j]   >= 0.0);
        CHECK(cm.logT_bulge[j]   >= 0.0);
        CHECK(cm.Mab_bulge[j][2] <= 30.0 or cm.Mab_bulge[j][2] == DARK_MAG);
        CHECK(cm.age_bulge[j]    <= 10);
        CHECK(cm.cl_bulge[j]     <= 7);
        CHECK(cm.typ_bulge[j]    <= 9.2);
    }

    CHECK(!(fp2 >> dummy));
    fp2.close();

    // ================================ THICK DISK =============================
    fp2.open(PATH_CMD_THICK);
    if (!fp2.is_open()) {
        throw std::runtime_error(std::string("cannot read ") + PATH_CMD_THICK);
    }

    std::getline(fp2, header);   // Skip the header line
    for (size_t j = 0; j < N3; ++j) {
          if (!(fp2 >> cm.mass_thick[j]
                    >> cm.logT_thick[j]
                    >> Mbol
                    >> cm.age_thick[j]
                    >> Pop
                    >> cm.Mab_thick[j][6]
                    >> cm.Mab_thick[j][0]
                    >> cm.Mab_thick[j][1]
                    >> cm.Mab_thick[j][2]
                    >> cm.Mab_thick[j][3]
                    >> cm.Mab_thick[j][4]
                    >> cm.Mab_thick[j][5]
                    >> cm.cl_thick[j]
                    >> cm.typ_thick[j])) {
              throw std::runtime_error("Unexpected end of thick_disk.dat");
          }
        // CHECKs follow the read (see the thin disk); same bounds.
        CHECK(cm.mass_thick[j]   >= 0.0);
        CHECK(cm.logT_thick[j]   >= 0.0);
        CHECK(cm.Mab_thick[j][2] <= 30.0 or cm.Mab_thick[j][2] == DARK_MAG);
//        CHECK(cm.age_thick[j]    <= 8);
        CHECK(cm.age_thick[j]    <= 13);
        CHECK(cm.cl_thick[j]     <= 7);
        CHECK(cm.typ_thick[j]    <= 9.2);
    }

    CHECK(!(fp2 >> dummy));
    fp2.close();

    // ================================ STELLAR HALO ===========================
    fp2.open(PATH_CMD_HALO);
    if (!fp2.is_open()) {
        throw std::runtime_error(std::string("cannot read ") + PATH_CMD_HALO);
    }

    std::getline(fp2, header);   // Skip the header line
    for (size_t j = 0; j < N4; ++j) {
          if (!(fp2 >> cm.mass_halo[j]
                    >> cm.logT_halo[j]
                    >> Mbol
                    >> cm.age_halo[j]
                    >> Pop
                    >> cm.Mab_halo[j][6]
                    >> cm.Mab_halo[j][0]
                    >> cm.Mab_halo[j][1]
                    >> cm.Mab_halo[j][2]
                    >> cm.Mab_halo[j][3]
                    >> cm.Mab_halo[j][4]
                    >> cm.Mab_halo[j][5]
                    >> cm.cl_halo[j]
                    >> cm.typ_halo[j])) {
              throw std::runtime_error("Unexpected end of halo.dat");
          }
        // CHECKs follow the read (see the thin disk); same bounds.
        CHECK(cm.mass_halo[j]   >= 0.0);
        CHECK(cm.logT_halo[j]   >= 0.0);
        CHECK(cm.Mab_halo[j][2] <= 30.0 or cm.Mab_halo[j][2] == DARK_MAG);
//        CHECK(cm.age_halo[j]    <= 9);
        CHECK(cm.age_halo[j]    <= 14);
        CHECK(cm.cl_halo[j]     <= 7);
        CHECK(cm.typ_halo[j]    <= 9.2);
    }

    CHECK(!(fp2 >> dummy));
    fp2.close();

    std::cout << ">>>>>>>>>>>>>>>>>> END OF CMD READING <<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
}

void readLensML(const std::string& path)
{
    std::ifstream fin(path);
    if (!fin) {
        std::cerr << "ERROR: cannot read " << path << " (run CMD/lens_ml_table.py)\n";
        std::exit(EXIT_FAILURE);
    }
    std::string line;
    int n = 0;
    while (std::getline(fin, line)) {
        if (line.empty() or line[0] == '#') continue;
        std::istringstream ss(line);
        int comp, cnt; double lo, hi; std::array<double, 7> m{};
        if (!(ss >> comp >> lo >> hi >> cnt)) continue;
        for (auto& v : m) ss >> v;
        if (!ss or comp < 0 or comp > 3 or !(hi > lo)) {
            std::cerr << "ERROR: malformed line in " << path << ": '" << line << "'\n";
            std::exit(EXIT_FAILURE);
        }
        // File order is u g r i z y F146, the simulator's filter order.
        gLensML.mmid[comp].push_back(0.5 * (lo + hi));
        gLensML.mab[comp].push_back(m);
        ++n;
    }
    for (int c = 0; c < 4; ++c)
        if (gLensML.mmid[c].size() < 2) {
            std::cerr << "ERROR: " << path << " has fewer than 2 bins for component " << c << "\n";
            std::exit(EXIT_FAILURE);
        }
    std::cout << "Loaded " << n << " luminous-lens mass-magnitude rows from " << path << "\n";
}

bool lensAbsMag(int comp, double mass, std::array<double, 7>& mab)
{
    const auto& x = gLensML.mmid[comp];
    const auto& y = gLensML.mab[comp];
    if (mass < KROUPA_BREAK1) return false;
    if (mass <= x.front()) { mab = y.front(); return true; }
    if (mass >= x.back())  { mab = y.back();  return true; }
    size_t i = 1;
    while (x[i] < mass) ++i;
    const double f = (mass - x[i - 1]) / (x[i] - x[i - 1]);
    for (int b = 0; b < 7; ++b) mab[b] = y[i - 1][b] + f * (y[i][b] - y[i - 1][b]);
    return true;
}
