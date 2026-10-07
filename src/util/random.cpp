// Scalar random draws on the global generator `rng` (util/random.h).
#include "util/random.h"

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double RandN(double sigm, double nnd) {
    std::normal_distribution<double> normal(0.0, 1.0);

    double x;
    do {
        x = sigm * normal(rng);
    } while (std::abs(x) > sigm * nnd); //[-N sigma:N sigma]

    return x;
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
// Poisson deviate, for counting things that are counted rather than measured -- currently
// the number of unresolved neighbours in a seeing disc (src/events/source.cpp). A Gaussian of width
// sqrt(mean) is the right approximation only for mean >> 1; at mean < 1 it is qualitatively
// wrong, because it cannot produce the integer 1 with the right probability and, truncated,
// cannot produce it at all. See the blending block in src/events/source.cpp.
int RandPois(double mean) {
    CHECK(mean >= 0.0);
    CHECK(std::isfinite(mean));
    std::poisson_distribution<int> pois(mean);
    return pois(rng);
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double RandR(double down, double up) {
    std::uniform_real_distribution<double> dist(down, up);

    return dist(rng);
}
