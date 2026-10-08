// Extinction tables (read, nearest sightline, interpolation along the line of sight) and CCM89.
#include "galaxy/extinction.h"

// Finds the extinction table (one of `nTables` unique l,b pointings) closest to
// (lon, lat). Call once per field pointing, not once per star.
int nearestSightline(const extin& ex, double lon, double lat) {
    int best = -1;

    double bestD2 = std::numeric_limits<double>::max();

    for (int k = 0; k < ex.nTables; k++) {
        double dl = ex.l[k] - lon;
        double db = ex.b[k] - lat;

        double d2 = dl * dl + db * db;

        if(d2 < bestD2){
            bestD2 = d2;
            best = k;
        }
    }

    CHECK(best >= 0);

    return best;
}

// Linearly interpolates A_V at distance `dist` (kpc) along table k; held flat beyond the grid.
double interpExtinctionAlongSightline(const extin& ex, int k, double dist) {
    const float* ext = ex.ext.data() + size_t(k) * size_t(ex.nDist);
    const int n = ex.nDist;

    if(dist <= ex.dist[0])
        return ext[0];

    if(dist >= ex.dist[n-1])
        return ext[n-1];

    int lo = 0;
    int hi = n-1;

    while(hi-lo > 1){
        int mid = lo + (hi-lo)/2;

        if(ex.dist[mid] <= dist)
            lo = mid;
        else
            hi = mid;
    }

    double t = (dist - ex.dist[lo]) / (ex.dist[hi] - ex.dist[lo]);

    return double(ext[lo]) + t*(double(ext[hi])-double(ext[lo]));
}

// Reads files/ext/ext_tables.dat line by line straight into `ex`. Format, written by maps.py:
//     # ext_tables v1 -- built by maps.py <date> ...      (provenance)
//     # k <A_Ks/A_V>
//     # n_tables <N>
//     # n_dist <M>
//     # dist d_1 ... d_M                                   (kpc, strictly increasing)
//     l b A_V(d_1) ... A_V(d_M)                            (N lines)
// It exits naming the file and line on a missing header field, a row with too few or too many
// numbers, a non-finite or negative A_V, a profile that decreases with distance, or a row count
// that differs from the header.
void readExtinction(extin& ex, const std::string& path) {
    auto fail = [&](long line, const std::string& why) {
        std::cerr << "ERROR: " << path << (line > 0 ? ":" + std::to_string(line) : std::string())
                  << ": " << why << ". Rebuild with maps.py.\n";
        std::exit(EXIT_FAILURE);
    };
    std::ifstream fin(path);
    if (!fin) fail(0, "cannot open the extinction tables");

    std::string line;
    long lineNo = 0, rows = 0;
    ex = extin();
    while (std::getline(fin, line)) {
        ++lineNo;
        if (line.empty()) continue;
        if (line[0] == '#') {
            std::istringstream ss(line.substr(1));
            std::string key;
            ss >> key;
            if (key == "ext_tables") ex.built = line.substr(1);
            else if (key == "k") ss >> ex.k;
            else if (key == "n_tables") ss >> ex.nTables;
            else if (key == "n_dist") ss >> ex.nDist;
            else if (key == "dist") {
                double d;
                while (ss >> d) ex.dist.push_back(d);
            }
            continue;
        }
        if (rows == 0) {                                  // header complete: check, then allocate
            if (ex.nTables <= 0 or ex.nDist < 2) fail(lineNo, "header lacks n_tables / n_dist");
            if (int(ex.dist.size()) != ex.nDist) fail(lineNo, "'# dist' does not have n_dist values");
            for (int i = 1; i < ex.nDist; ++i)
                if (!(ex.dist[i] > ex.dist[i - 1])) fail(lineNo, "distance grid not increasing");
            ex.l.resize(ex.nTables);
            ex.b.resize(ex.nTables);
            ex.ext.resize(size_t(ex.nTables) * size_t(ex.nDist));
        }
        if (rows >= ex.nTables) fail(lineNo, "more table rows than n_tables");
        const char* p = line.c_str();
        char* end = nullptr;
        ex.l[rows] = std::strtod(p, &end);
        if (end == p) fail(lineNo, "cannot parse l");
        p = end;
        ex.b[rows] = std::strtod(p, &end);
        if (end == p) fail(lineNo, "cannot parse b");
        p = end;
        float* row = ex.ext.data() + size_t(rows) * size_t(ex.nDist);
        for (int i = 0; i < ex.nDist; ++i) {
            const double v = std::strtod(p, &end);
            if (end == p) fail(lineNo, "row has fewer than n_dist values");
            if (!std::isfinite(v) or v < 0.0) fail(lineNo, "non-finite or negative A_V");
            if (i > 0 and float(v) < row[i - 1]) fail(lineNo, "A_V decreases with distance");
            row[i] = float(v);
            p = end;
        }
        while (*p == ' ' or *p == '\t' or *p == '\r') ++p;
        if (*p != '\0') fail(lineNo, "row has more than n_dist values");
        ++rows;
    }
    if (rows != ex.nTables)
        fail(0, "read " + std::to_string(rows) + " table rows, header says " + std::to_string(ex.nTables));

    std::cout << "Loaded " << ex.nTables << " extinction tables x " << ex.nDist
              << " distances (k = " << ex.k << ") from " << path << "\n";
}

double CCM89_a(double lambda_um)
{
    double x = 1.0 / lambda_um;

    if (x < 1.1) {
        return 0.574 * std::pow(x, 1.61);
    }

    if (x < 3.3) {
        double y = x - 1.82;

        return 1.0
             + 0.17699 * y
             - 0.50447 * std::pow(y, 2)
             - 0.02427 * std::pow(y, 3)
             + 0.72085 * std::pow(y, 4)
             + 0.01979 * std::pow(y, 5)
             - 0.77530 * std::pow(y, 6)
             + 0.32999 * std::pow(y, 7);
    }

    throw std::runtime_error("CCM89 wavelength out of supported range");
}

double CCM89_b(double lambda_um)
{
    double x = 1.0 / lambda_um;

    if (x < 1.1) {
        return -0.527 * std::pow(x, 1.61);
    }

    if (x < 3.3) {
        double y = x - 1.82;

        return 1.41338 * y
             + 2.28305 * std::pow(y, 2)
             + 1.07233 * std::pow(y, 3)
             - 5.38434 * std::pow(y, 4)
             - 0.62251 * std::pow(y, 5)
             + 5.30260 * std::pow(y, 6)
             - 2.09002 * std::pow(y, 7);
    }

    throw std::runtime_error("CCM89 wavelength out of supported range");
}

// A_lambda / A_V from Cardelli, Clayton & Mathis (1989). CCM89_a/b take a WAVELENGTH in
// micron and form x = 1/lambda themselves, so lambda is passed straight through.
// tests/extinction_test.cpp pins the values.
double AlAv(double lambda_um, double Rv)
{
    return CCM89_a(lambda_um) + CCM89_b(lambda_um) / Rv;
}
