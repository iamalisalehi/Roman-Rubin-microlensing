// Unit test for the buffered number reader (src/util/parse.cpp) and the read_cmd that uses it.
// Exit 0 = all held. Needs CMD/components/*.dat (a few GB-seconds of reading; run from the repo root).
//
// (a) the four component files, token by token, against today's `std::ifstream >> double`;
// (b) the whole read_cmd against read_cmd_ref (a verbatim copy of the istream version, the oracle),
//     every array of the CMD struct hashed byte by byte; the two CMD structs are never alive together;
// (c) synthetic token forms (nan, inf, signs, exponents, subnormals, split tokens, CRLF, tabs, tiny
//     buffers that force refills) against istream, values bit for bit and accept/reject;
// (d) is the mutation check, done by hand: make the fast path round differently and (a), (b) must fail.
#include "common.h"
#include "galaxy/catalogue.h"
#include "util/parse.h"
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <memory>
#include <string>

static int nfail = 0;
static void expect(bool ok, const std::string& what)
{
    std::printf("%-70s %s\n", what.c_str(), ok ? "ok" : "FAIL");
    std::fflush(stdout);
    if (!ok) ++nfail;
}

static bool sameBits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

// ---- the oracle: read_cmd as it was before TextReader (copied verbatim, names suffixed) ----
// Stellar radius [Rsun] from the bolometric magnitude and log10 Teff: L/Lsun = 10^(-0.4 (Mbol - MbolSun))
// and L = 4 pi R^2 sigma Teff^4 give R/Rsun = 10^(-0.2 (Mbol - MbolSun)) (TeffSun / Teff)^2. The same
// zero points built the star lists, so this agrees with their own logg (g ~ M / R^2).
static double stellarRadius_ref(double Mbol, double logT)
{
    const double tRatio = TeffSun / std::pow(10.0, logT);
    return std::pow(10.0, -0.2 * (Mbol - MbolSun)) * tRatio * tRatio;
}

void read_cmd_ref(CMD & cm)
{
// mass  logT  Mbol  Age  Pop  Roman_F146  LSST_u  LSST_g  LSST_r  LSST_i  LSST_z  LSST_y  CL  Type
// 0     1     2     3    4    5           6       7       8       9       10      11      12  13
    double Mbol, Pop;
    double dummy;
    std::string header;

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
        cm.Rstar_thin[j] = stellarRadius_ref(Mbol, cm.logT_thin[j]);
        // M_r <= 30: the thin disc's luminous brown-dwarf-limit stars reach M_r = 28.86;
        // Typ <= 9.2: white dwarfs are Typ 9.0-9.2 and sit in the list as dark entries.
        CHECK(cm.mass_thin[j]   >= 0.0);
        CHECK(cm.logT_thin[j]   >= 0.0);
        CHECK(cm.Mab_thin[j][2] <= 30.0 or cm.Mab_thin[j][2] == DARK_MAG);
        CHECK(cm.age_thin[j]    <= 10);
        CHECK(cm.cl_thin[j]     <= 7);
        CHECK(cm.typ_thin[j]    <= 9.2);
        CHECK(std::isfinite(cm.Rstar_thin[j]) and cm.Rstar_thin[j] > 0.0);
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
        cm.Rstar_bulge[j] = stellarRadius_ref(Mbol, cm.logT_bulge[j]);
        // CHECKs follow the read (see the thin disk); same bounds.
        CHECK(cm.mass_bulge[j]   >= 0.0);
        CHECK(cm.logT_bulge[j]   >= 0.0);
        CHECK(cm.Mab_bulge[j][2] <= 30.0 or cm.Mab_bulge[j][2] == DARK_MAG);
        CHECK(cm.age_bulge[j]    <= 10);
        CHECK(cm.cl_bulge[j]     <= 7);
        CHECK(cm.typ_bulge[j]    <= 9.2);
        CHECK(std::isfinite(cm.Rstar_bulge[j]) and cm.Rstar_bulge[j] > 0.0);
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
        cm.Rstar_thick[j] = stellarRadius_ref(Mbol, cm.logT_thick[j]);
        // CHECKs follow the read (see the thin disk); same bounds.
        CHECK(cm.mass_thick[j]   >= 0.0);
        CHECK(cm.logT_thick[j]   >= 0.0);
        CHECK(cm.Mab_thick[j][2] <= 30.0 or cm.Mab_thick[j][2] == DARK_MAG);
        CHECK(cm.age_thick[j]    <= 13);
        CHECK(cm.cl_thick[j]     <= 7);
        CHECK(cm.typ_thick[j]    <= 9.2);
        CHECK(std::isfinite(cm.Rstar_thick[j]) and cm.Rstar_thick[j] > 0.0);
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
        cm.Rstar_halo[j] = stellarRadius_ref(Mbol, cm.logT_halo[j]);
        // CHECKs follow the read (see the thin disk); same bounds.
        CHECK(cm.mass_halo[j]   >= 0.0);
        CHECK(cm.logT_halo[j]   >= 0.0);
        CHECK(cm.Mab_halo[j][2] <= 30.0 or cm.Mab_halo[j][2] == DARK_MAG);
        CHECK(cm.age_halo[j]    <= 14);
        CHECK(cm.cl_halo[j]     <= 7);
        CHECK(cm.typ_halo[j]    <= 9.2);
        CHECK(std::isfinite(cm.Rstar_halo[j]) and cm.Rstar_halo[j] > 0.0);
    }

    CHECK(!(fp2 >> dummy));
    fp2.close();

    std::cout << ">>>>>>>>>>>>>>>>>> END OF CMD READING <<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
}

// ---- (a) real files, token by token ----
static void testFileTokens(const char* path, size_t rows)
{
    std::ifstream fp(path);
    TextReader rd(path);
    std::string header;
    std::getline(fp, header);
    rd.skipLine();
    size_t n = 0, bad = 0;
    bool flagsDiffer = false;
    for (;;) {
        double a = 0, b = 0;
        const bool oa = static_cast<bool>(fp >> a);
        const bool ob = rd.next(b);
        if (oa != ob) { flagsDiffer = true; break; }
        if (!oa) break;
        if (!sameBits(a, b)) ++bad;
        ++n;
    }
    expect(!flagsDiffer, std::string("(a) ") + path + ": EOF at the same token");
    expect(bad == 0, std::string("(a) ") + path + ": " + std::to_string(n) + " tokens, " + std::to_string(bad) + " differ");
    expect(n == 14 * rows, std::string("(a) ") + path + ": 14 x " + std::to_string(rows) + " tokens");
}

// ---- (b) hash of the whole struct ----
static void fnv(std::uint64_t& h, const double* p, size_t n)
{
    const unsigned char* b = reinterpret_cast<const unsigned char*>(p);
    for (size_t i = 0; i < n * sizeof(double); ++i) { h ^= b[i]; h *= 1099511628211ull; }
}
static void fnvVec(std::uint64_t& h, const std::vector<double>& v) { fnv(h, v.data(), v.size()); }
static void fnvMab(std::uint64_t& h, const std::vector<std::array<double, M>>& v)
{
    for (const auto& r : v) fnv(h, r.data(), r.size());
}
static std::uint64_t hashCMD(const CMD& c)
{
    std::uint64_t h = 14695981039346656037ull;
    fnvVec(h, c.logT_thin);  fnvVec(h, c.Rstar_thin);  fnvVec(h, c.mass_thin);  fnvMab(h, c.Mab_thin);  fnvVec(h, c.typ_thin);  fnvVec(h, c.cl_thin);  fnvVec(h, c.age_thin);
    fnvVec(h, c.logT_bulge); fnvVec(h, c.Rstar_bulge); fnvVec(h, c.mass_bulge); fnvMab(h, c.Mab_bulge); fnvVec(h, c.typ_bulge); fnvVec(h, c.cl_bulge); fnvVec(h, c.age_bulge);
    fnvVec(h, c.logT_thick); fnvVec(h, c.Rstar_thick); fnvVec(h, c.mass_thick); fnvMab(h, c.Mab_thick); fnvVec(h, c.typ_thick); fnvVec(h, c.cl_thick); fnvVec(h, c.age_thick);
    fnvVec(h, c.logT_halo);  fnvVec(h, c.Rstar_halo);  fnvVec(h, c.mass_halo);  fnvMab(h, c.Mab_halo);  fnvVec(h, c.typ_halo);  fnvVec(h, c.cl_halo);  fnvVec(h, c.age_halo);
    return h;
}

// ---- (c) synthetic tokens ----
static std::string tmpDir()
{
    const char* d = std::getenv("PARSETEST_TMP");
    return d ? d : "/tmp";
}

// Both readers over the same file, read in lockstep until the first reject or the end.
// Returns false if the accept flags or any value differ; `log` gets the per-read values.
static bool lockstep(const std::string& file, size_t bufBytes, std::string* log)
{
    std::ifstream fp(file);
    TextReader rd(file.c_str(), bufBytes);
    std::string h;
    std::getline(fp, h);
    rd.skipLine();
    for (int k = 0; k < 8; ++k) {
        double a = 0, b = 0;
        const bool oa = static_cast<bool>(fp >> a);
        const bool ob = rd.next(b);
        if (log) {
            char s[160];
            std::snprintf(s, sizeof s, "[%s %.17g | %s %.17g] ", oa ? "acc" : "rej", a, ob ? "acc" : "rej", b);
            *log += s;
        }
        if (oa != ob) return false;
        if (!oa) return true;
        if (!sameBits(a, b)) return false;
    }
    return true;
}

static void testSynthetic()
{
    const char* toks[] = {"1", "-0", "+1.5", "1e-05", "1E+3", ".5", "5.", "00012", "1e-310", "1e-400", "1e400",
                          "nan", "inf", "-inf", "0x1p3", "1e", "1.5e+", "-", "+", "1-2",
                          "0.1", "123456789.123456789", "-7.25e-3", "1.7976931348623157e308", "4.9e-324"};
    // each token alone on a line (with CRLF, tabs, trailing spaces around it), then more numbers after it
    const char* pre[]  = {"", "  ", "\t", "\r\n", " \t \v\f"};
    const char* post[] = {"\n", "\r\n", "   \n", "\t\n", " "};
    const size_t bufs[] = {4, 5, 16, size_t(1) << 20};
    const std::string path = tmpDir() + "/parse_test_synthetic.dat";
    int mismatches = 0;
    std::printf("(c) token | istream, first reads | TextReader, first reads (default buffer, no padding)\n");
    for (const char* t : toks) {
        std::string shown;
        for (size_t pi = 0; pi < 5; ++pi) {
            for (size_t qi = 0; qi < 5; ++qi) {
                std::ofstream f(path, std::ios::binary);
                f << "# header line, skipped\r\n" << pre[pi] << t << post[qi] << "42 -0 7.5e1\n";
                f.close();
                for (size_t bsz : bufs) {
                    const bool ok = lockstep(path, bsz, pi == 0 && qi == 0 && bsz == bufs[3] ? &shown : nullptr);
                    if (!ok) ++mismatches;
                }
            }
        }
        std::printf("   %-26s %s\n", t, shown.c_str());
    }
    expect(mismatches == 0, "(c) synthetic tokens x padding x buffer sizes: istream and TextReader agree (" + std::to_string(mismatches) + " mismatches)");

    // A file with no trailing newline, and an empty file after the header.
    {
        std::ofstream f(path, std::ios::binary);
        f << "hdr\n1.5 2.5";
        f.close();
        expect(lockstep(path, 4, nullptr), "(c) no newline after the last token");
        std::ofstream g(path, std::ios::binary);
        g << "hdr only";
        g.close();
        expect(lockstep(path, 4, nullptr), "(c) header only, then EOF");
    }
    // A long token (longer than the buffer) is read whole.
    {
        std::string longTok = "0." + std::string(200, '1');
        std::ofstream f(path, std::ios::binary);
        f << "hdr\n" << longTok << " 3\n";
        f.close();
        expect(lockstep(path, 4, nullptr), "(c) a 202-character token through a 4-byte buffer");
    }
}

int main()
{
    testSynthetic();

    testFileTokens(PATH_CMD_HALO,  N4);
    testFileTokens(PATH_CMD_THIN,  N1);
    testFileTokens(PATH_CMD_THICK, N3);
    testFileTokens(PATH_CMD_BULGE, N2);

    std::uint64_t hNew, hRef;
    {
        auto cm = std::make_unique<CMD>();
        read_cmd(*cm);
        hNew = hashCMD(*cm);
    }
    {
        auto cm = std::make_unique<CMD>();
        read_cmd_ref(*cm);
        hRef = hashCMD(*cm);
    }
    std::printf("(b) hash read_cmd %016llx  read_cmd_ref %016llx\n", (unsigned long long)hNew, (unsigned long long)hRef);
    expect(hNew == hRef, "(b) read_cmd and read_cmd_ref fill every CMD array identically");

    std::printf("%s (%d failure%s)\n", nfail ? "FAIL" : "all held", nfail, nfail == 1 ? "" : "s");
    return nfail ? 1 : 0;
}
