// Readers and coverage tests for the GBTDS detector layout and the LSSTCam footprint map.
#include "surveys/footprints.h"

// ---------------------------------------------------------------------------------------
// GBTDS detector layout (Deviation 69). Each file lists 18 detector outlines as 5 vertices
// `sca dl db` (the fifth closes the rectangle), offsets in deg from the field centre, in the
// convention of the upstream tool: the sky outline is the offset ADDED to the field's (l, b).
// Every detector must be an axis-aligned rectangle in (l, b); anything else is refused rather
// than approximated, because the coverage test below assumes it.
// ---------------------------------------------------------------------------------------
GbtdsLayout readGbtdsLayout() {
    GbtdsLayout g;
    g.scaSide = std::numeric_limits<double>::max();
    for (int k = 0; k < GBTDS_NLAYOUT; ++k) {
        const char* path = GBTDS_SCA_FILES[k];
        std::ifstream fin(path);
        if (!fin) {
            std::cerr << "ERROR: cannot read the GBTDS detector layout " << path << "\n";
            std::exit(EXIT_FAILURE);
        }
        std::map<int, std::vector<std::pair<double,double>>> vert;
        std::string line;
        while (std::getline(fin, line)) {
            std::istringstream ss(line);
            int id; double dl, db;
            if (!(ss >> id)) continue;                         // blank separator line
            if (!(ss >> dl >> db) or !std::isfinite(dl) or !std::isfinite(db)) {
                std::cerr << "ERROR: malformed line in " << path << ": '" << line << "'\n";
                std::exit(EXIT_FAILURE);
            }
            vert[id].emplace_back(dl, db);
        }
        if (int(vert.size()) != GBTDS_NSCA) {
            std::cerr << "ERROR: " << path << " has " << vert.size() << " detectors, expected "
                      << GBTDS_NSCA << "\n";
            std::exit(EXIT_FAILURE);
        }
        g.dlMin[k] = g.dbMin[k] =  std::numeric_limits<double>::max();
        g.dlMax[k] = g.dbMax[k] = -std::numeric_limits<double>::max();
        for (const auto& [id, v] : vert) {
            double l0 = v[0].first, l1 = l0, b0 = v[0].second, b1 = b0;
            for (const auto& p : v) {
                l0 = std::min(l0, p.first);  l1 = std::max(l1, p.first);
                b0 = std::min(b0, p.second); b1 = std::max(b1, p.second);
            }
            // Axis-aligned: every vertex sits on one of the four edges' coordinates.
            const double tol = 1e-6;
            for (const auto& p : v) {
                const bool onL = std::fabs(p.first - l0) < tol or std::fabs(p.first - l1) < tol;
                const bool onB = std::fabs(p.second - b0) < tol or std::fabs(p.second - b1) < tol;
                if (!(onL and onB)) {
                    std::cerr << "ERROR: detector " << id << " in " << path
                              << " is not an axis-aligned rectangle in (l, b)\n";
                    std::exit(EXIT_FAILURE);
                }
            }
            g.sca[k].push_back({l0, l1, b0, b1});
            g.dlMin[k] = std::min(g.dlMin[k], l0); g.dlMax[k] = std::max(g.dlMax[k], l1);
            g.dbMin[k] = std::min(g.dbMin[k], b0); g.dbMax[k] = std::max(g.dbMax[k], b1);
            g.scaSide  = std::min({g.scaSide, l1 - l0, b1 - b0});
            for (double cl : {l0, l1})
                for (double cb : {b0, b1})
                    g.rField = std::max(g.rField, std::hypot(cl, cb));
        }
    }
    return g;
}

// ---------------------------------------------------------------------------------------
// LSSTCam footprint (Deviation 80). See include/surveys/footprints.h.
// ---------------------------------------------------------------------------------------
void readLsstCamMap(const std::string& path)
{
    std::ifstream fin(path);
    if (!fin) { std::cerr << "ERROR: cannot read " << path << " (run Baseline/lsstcam_fov/export_fov_map.py)\n"; std::exit(EXIT_FAILURE); }
    std::string line;
    std::getline(fin, line);                       // provenance
    std::getline(fin, line);                       // '# n x0 step max_radius'
    { std::istringstream ss(line.substr(1)); double rmax; ss >> gLsstCam.n >> gLsstCam.x0 >> gLsstCam.step >> rmax;
      if (!ss or gLsstCam.n <= 0 or !(gLsstCam.step > 0)) { std::cerr << "ERROR: bad header in " << path << "\n"; std::exit(EXIT_FAILURE); } }
    gLsstCam.on.assign(size_t(gLsstCam.n) * gLsstCam.n, 0);
    for (int ix = 0; ix < gLsstCam.n; ++ix) {
        if (!std::getline(fin, line) or int(line.size()) < gLsstCam.n) {
            std::cerr << "ERROR: " << path << " row " << ix << " short or missing\n"; std::exit(EXIT_FAILURE);
        }
        for (int iy = 0; iy < gLsstCam.n; ++iy) gLsstCam.on[size_t(ix) * gLsstCam.n + iy] = (line[iy] == '1');
    }
    size_t nOn = 0; for (auto v : gLsstCam.on) nOn += v;
    std::cout << "Loaded LSSTCam footprint " << path << ": " << nOn * gLsstCam.step * gLsstCam.step
              << " deg^2 active\n";
}

// rubin_scheduler.utils.LsstCameraFootprint.__call__, one point: gnomonic projection about the
// boresight (projections.gnomonic_project_toxy), rotate by rotSkyPos (its rotate()), nearest pixel,
// image[ix][iy].
bool onLsstCam(double ra, double dec, double ra0, double dec0, double rotSkyPos)
{
    const double d2r = M_PI / 180.0;
    const double a = ra * d2r, d = dec * d2r, a0 = ra0 * d2r, d0 = dec0 * d2r;
    const double cosc = std::sin(d0) * std::sin(d) + std::cos(d0) * std::cos(d) * std::cos(a - a0);
    if (cosc <= 0.0) return false;
    double x = std::cos(d) * std::sin(a - a0) / cosc;
    double y = (std::cos(d0) * std::sin(d) - std::sin(d0) * std::cos(d) * std::cos(a - a0)) / cosc;
    const double r = rotSkyPos * d2r, c = std::cos(r), s = std::sin(r);
    const double xr = c * x - s * y, yr = s * x + c * y;
    const double stepRad = gLsstCam.step * d2r, x0Rad = gLsstCam.x0 * d2r;
    const long ix = std::lround((xr - x0Rad) / stepRad), iy = std::lround((yr - x0Rad) / stepRad);
    if (ix < 0 or iy < 0 or ix >= gLsstCam.n or iy >= gLsstCam.n) return false;
    return gLsstCam.on[size_t(ix) * gLsstCam.n + size_t(iy)] != 0;
}

// Galactic -> ICRS (J2000), via the transpose of the Hipparcos ICRS->Galactic rotation matrix
// (ESA 1997, vol. 1, sec. 1.5.3); checked at start-up against every visit's OpSim RA/Dec vs l/b.
void galToIcrs(double l, double b, double& ra, double& dec)
{
    static const double A[3][3] = {{-0.0548755604162154, -0.8734370902348850, -0.4838350155487132},
                                   {+0.4941094278755837, -0.4448296299600112, +0.7469822444972189},
                                   {-0.8676661490190047, -0.1980763734312015, +0.4559837761750669}};
    const double d2r = M_PI / 180.0;
    const double g[3] = {std::cos(b * d2r) * std::cos(l * d2r), std::cos(b * d2r) * std::sin(l * d2r), std::sin(b * d2r)};
    double e[3];
    for (int i = 0; i < 3; ++i) e[i] = A[0][i] * g[0] + A[1][i] * g[1] + A[2][i] * g[2];   // A^T g
    ra  = std::atan2(e[1], e[0]) / d2r; if (ra < 0.0) ra += 360.0;
    dec = std::asin(std::max(-1.0, std::min(1.0, e[2]))) / d2r;
}

bool inDetector(const GbtdsLayout& g, int layout, double dl, double db) {
    if (dl < g.dlMin[layout] or dl > g.dlMax[layout] or db < g.dbMin[layout] or db > g.dbMax[layout])
        return false;
    for (const ScaRect& r : g.sca[layout])
        if (dl >= r.l0 and dl <= r.l1 and db >= r.b0 and db <= r.b1) return true;
    return false;
}
