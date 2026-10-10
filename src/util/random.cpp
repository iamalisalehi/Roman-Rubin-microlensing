// Scalar random draws on the global generator `rng` (util/random.h).
#include "util/random.h"

// Gaussian truncated to +-nnd sigma, by redrawing.
double RandN(double sigm, double nnd) {
    std::normal_distribution<double> normal(0.0, 1.0);

    double x;
    do {
        x = sigm * normal(rng);
    } while (std::abs(x) > sigm * nnd); //[-N sigma:N sigma]

    return x;
}

// Poisson deviate, used for the number of unresolved neighbours in a seeing disc
// (src/events/source.cpp). A Gaussian of width sqrt(mean) is adequate only for mean >> 1.
int RandPois(double mean) {
    CHECK(mean >= 0.0);
    CHECK(std::isfinite(mean));
    std::poisson_distribution<int> pois(mean);
    return pois(rng);
}

// Uniform on [0, 1) from the neighbour-position stream.
double RandBlendUnit() {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rngBlend);
}

// Uniform on [down, up).
double RandR(double down, double up) {
    std::uniform_real_distribution<double> dist(down, up);

    return dist(rng);
}
