// Which sky points each survey images: the GBTDS detector layout, the LSSTCam active-silicon map,
// the coordinate transform, and the instrument-agnostic visible-epoch matcher.
#ifndef ROMAN_SURVEYS_FOOTPRINTS_H
#define ROMAN_SURVEYS_FOOTPRINTS_H

#include "common.h"

struct ScaRect { double l0, l1, b0, b1; };   // offsets from the field centre [deg], l0<l1, b0<b1

struct GbtdsLayout {
    std::array<std::vector<ScaRect>, GBTDS_NLAYOUT> sca;
    // Bounding rectangle of each layout's detectors, for a cheap early reject and for the fine
    // stratum of the sightline grid.
    std::array<double, GBTDS_NLAYOUT> dlMin{}, dlMax{}, dbMin{}, dbMax{};
    double rField  = 0.0;   // largest centre-to-detector-corner distance, either layout [deg]
    double scaSide = 0.0;   // smallest detector side [deg]
};

// A field placement: one (centre, roll) that RomanBaseline.dat visits. 6 fields x 2 rolls = 12.
struct FieldPlacement { double l, b; int layout; };

// LSSTCam's footprint (Deviation 80): the active-silicon map OpSim/MAF use (rubin_scheduler
// LsstCameraFootprint, fov_map.npz), exported to Baseline/lsstcam_fov/fov_map.txt. A sky point is on
// a Rubin visit's silicon if its gnomonic projection about the boresight, rotated by rotSkyPos,
// falls on an active pixel -- rubin_scheduler's own algorithm, mirrored. Replaces the 1.75-deg circle.
// RUBIN_MAX_RADIUS (1.94 deg, rubin_scheduler's max_radius): see config/parameters.h, section 4.
struct LsstCamMap { int n = 0; double x0 = 0.0, step = 0.0; std::vector<unsigned char> on; };
inline LsstCamMap gLsstCam;
void   readLsstCamMap(const std::string& path);
bool   onLsstCam(double ra, double dec, double ra0, double dec0, double rotSkyPos);   // all deg
void   galToIcrs(double l, double b, double& ra, double& dec);                       // deg (J2000)

// GBTDS detector layout (helper.cpp). readGbtdsLayout exits with the file named on any
// malformed input; inDetector tests a sky offset (dl, db) from a field centre against one
// layout's 18 detector rectangles.
GbtdsLayout readGbtdsLayout();
bool inDetector(const GbtdsLayout& g, int layout, double dl, double db);

///==============================================================//
///     Instrument-agnostic epoch matching (LSST or Roman)       //
///                                                              //
///==============================================================//
// Extracted from the sky-position loop in main() so it can be called once per
// instrument instead of being hardwired to `ls`. Behavior is unchanged for LSST;
// calling it a second time with Roman's own visits and coverage test is what gives Roman its
// own epoch list instead of inheriting LSST's cadence.
//
// `covers(i)` says whether visit i images the current sky position. Rubin's is the circle
// sqrt(dl^2 + db^2) <= FoV about the pointing centre (unchanged arithmetic); Roman's is "on one
// of the 18 detectors of visit i's layout, placed at the field centre" (Deviation 69).
//
// `label` is purely diagnostic ("LSST" / "Roman") — it identifies which call
// printed a given warning, since both calls share this one function.
//
// Two visits can legitimately share an identical recorded timestamp for a given
// sky position — most commonly two different fields/pointings whose FoV circles
// overlap and which happen to share the same observing schedule (this is exactly
// what can happen between adjacent Roman fields, whose detector mosaics interlock and can
// overlap at the edges). There is no meaningful
// cadence between two simultaneous visits, so rather than treat this as fatally
// corrupt data, we keep the first and skip the duplicate — but we log it, because
// if this fires constantly (not just occasionally near field boundaries) that's a
// sign of a genuine sorting/data problem in the baseline file that needs fixing at
// the source, not papering over here.
template <typename Covers>
int matchVisibleEpochs(const char* label, Covers covers,
                       const std::vector<double>& tim_arr,
                       int nEpochs,
                       std::vector<int>& ct,
                       double& minCadence) {
    int ndd = 0;
    int nTies = 0;
    minCadence = 100000.0;

    const int ctCap = static_cast<int>(ct.size());
    for (int i = 0; i < ctCap; ++i) ct[i] = -1;

    for (int i = 0; i < nEpochs; ++i) {
        if (covers(i)) {

            if (ndd > 0) {
                double cade = tim_arr[i] - tim_arr[ct[ndd - 1]];

                if (cade <= 0.0) {
                    nTies++;
                    if (nTies <= 5) {
                        std::cerr << "[matchVisibleEpochs:" << label << "] skipping tied/out-of-order "
                                  << "epoch at index " << i << " (tim=" << tim_arr[i]
                                  << "), previous kept epoch tim=" << tim_arr[ct[ndd - 1]] << "\n";
                    }
                    continue; // keep the earlier one, don't record this as a new epoch
                }

                if (minCadence > cade) minCadence = cade;
                // Unreachable when ct is sized to the instrument's full visit count, which
                // is the point: ndd <= nEpochs by construction. Kept as a guard, and made
                // loud -- a SILENT truncation here previously ended Roman's mission on day
                // 8.4 of 1715 and Rubin's on day 1573 of 3459, with no diagnostic at all.
                if (ndd >= ctCap - 1) {
                    std::cerr << "[matchVisibleEpochs:" << label << "] FATAL: ct capacity "
                              << ctCap << " exhausted at epoch " << i << " of " << nEpochs
                              << " -- the visit stream would be silently truncated at tim="
                              << tim_arr[i] << " d. Size ct to the full visit count.\n";
                    std::exit(1);
                }
                CHECK(tim_arr[i] >= 0.0);
                CHECK(tim_arr[i] <= Tobs);
                CHECK(minCadence > 0.0);
            }
            ct[ndd] = i;
            ndd += 1;
        }
    }

    if (nTies > 5) {
        std::cerr << "[matchVisibleEpochs:" << label << "] ... " << (nTies - 5)
                  << " more skipped (total " << nTies << " tied/out-of-order epochs for this sky position)\n";
    } else if (nTies > 0) {
        std::cerr << "[matchVisibleEpochs:" << label << "] " << nTies << " tied/out-of-order epoch(s) skipped\n";
    }

    return ndd;
}

#endif // ROMAN_SURVEYS_FOOTPRINTS_H
