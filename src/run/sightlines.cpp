// Building the sightline list. Moved out of main() verbatim.
#include "run/sightlines.h"

int buildSightlines(const RunConfig& cfg, const GbtdsLayout& gl, const GridSteps& steps,
                    const lsst& ls, const roman& ro, SightlineGrid& grid) {
    const double gridStep = steps.gridStep;
    const int    kSub     = steps.kSub;
    const double fineStep = steps.fineStep;

    // ----------------------------------------------------------------------
    // Roman field placements and the coverage guard.
    //
    // RomanBaseline.dat repeats a small set of distinct pointings once per visit: the six
    // GBTDS fields at each of the two rolls, i.e. 12 (centre, layout) placements. Collect them,
    // then verify the sightline grid puts at least one sightline ON A DETECTOR of each. A grid
    // that misses a placement produces a run with no Roman epochs there, which does not crash
    // and does not warn: it quietly reports Rubin-only results in joint-labelled columns.
    // ----------------------------------------------------------------------
    std::vector<FieldPlacement> romanFields;
    for (int i = 0; i < NlRoman; ++i) {
        bool seen = false;
        for (const auto& f : romanFields)
            if (std::fabs(f.l - ro.l[i]) < 1e-6 and std::fabs(f.b - ro.b[i]) < 1e-6
                and f.layout == ro.layout[i]) {
                seen = true; break;
            }
        if (!seen) romanFields.push_back({ro.l[i], ro.b[i], ro.layout[i]});
    }

    // Scan bounds (Deviation 69). The region is every point within scanReach of a field centre
    // (see SCAN_RUBIN_REACH in config/parameters.h); the grid's bounding box is that, with its origin on a
    // multiple of the coarse step so sightlines sit at round coordinates. The stub is a 0.1 x
    // 0.1 deg test patch inside field 3, which both rolls image (spring centre l 0.500, autumn
    // 0.350; it was l 0.5-0.6, b -1.0..-0.9 before, which the adopted fields do not reach).
    const double scanReach = SCAN_RUBIN_REACH + gl.rField;
    double fLmin = 1e9, fLmax = -1e9, fBmin = 1e9, fBmax = -1e9;
    for (const auto& f : romanFields) {
        fLmin = std::min(fLmin, f.l); fLmax = std::max(fLmax, f.l);
        fBmin = std::min(fBmin, f.b); fBmax = std::max(fBmax, f.b);
    }
    const double lonMin = cfg.stubPatch ?  0.40  : gridStep * std::floor((fLmin - scanReach) / gridStep);
    const double lonMax = cfg.stubPatch ?  0.50  : fLmax + scanReach;
    const double latMin = cfg.stubPatch ? -1.45  : gridStep * std::floor((fBmin - scanReach) / gridStep);
    const double latMax = cfg.stubPatch ? -1.35  : fBmax + scanReach;

    // Inside the scan region: within scanReach of any field centre (the stub keeps its square).
    auto inScan = [&](double lon, double lat) {
        if (cfg.stubPatch) return true;
        for (const auto& f : romanFields)
            if (std::hypot(lon - f.l, lat - f.b) <= scanReach) return true;
        return false;
    };

    // Every Rubin pointing in the visit list must be one that can image the scan region, i.e.
    // centred within scanReach + FoV of a field centre. A pointing farther out means the list
    // was extracted for a different region (readbaselineBulge.py uses the same rule).
    if (not cfg.stubPatch) {
        int nFar = 0;
        for (int i = 0; i < Nl; ++i) {
            bool near = false;
            for (const auto& f : romanFields)
                if (std::hypot(ls.l[i] - f.l, ls.b[i] - f.b) <= scanReach + RUBIN_MAX_RADIUS + 1e-3) { near = true; break; }
            if (!near) nFar += 1;
        }
        if (nFar > 0) {
            std::cerr << "ERROR: " << nFar << " of " << Nl << " Rubin visits in BulgeBaseline.dat "
                      << "are centred farther than " << scanReach + RUBIN_MAX_RADIUS << " deg from every "
                      << "Roman field, so the visit list was built for another region. "
                      << "Regenerate it with Baseline/readbaselineBulge.py.\n";
            return 1;
        }
    }

    // ----------------------------------------------------------------------
    // The sightline list, with the sky area each sightline stands for (Step E1).
    //
    // WHY THIS IS NOT JUST A NESTED LOOP ANY MORE. The scan region is 68 deg^2; Roman's six
    // GBTDS fields cover ~1.7 deg^2 of it, about 2.6%. A uniform grid therefore spends 97% of
    // its wall clock on sightlines Roman never visits, where the joint Fisher matrix IS
    // Rubin's matrix and nothing whatever can be learned about the combination of the two
    // surveys. The 2026-08-30 production run bought 74,812 joint detections and only 1,950 of
    // them -- 2.6%, exactly the area fraction -- inside the footprint. Every result that is
    // currently sample-limited (F3 panel (a), the F4 precision fractions, the long-tE piE
    // null) is limited by that same 1,950.
    //
    // So the grid is stratified. Two strata:
    //   R  sightlines within FoVRoman of a GBTDS field centre, on the FINE grid (fineStep)
    //   O  everything else, on the COARSE grid (gridStep), one sightline per coarse cell
    // and each carries `area`, the deg^2 of sky it represents. An R sightline carries one
    // fine cell; an O sightline carries however many fine cells of its coarse block fall
    // outside the footprint -- NOT the whole coarse cell, or a block straddling the footprint
    // edge would count its overlapping part twice, once in each stratum.
    //
    // THE ONE INVARIANT: sum(area) over the list equals the scanned area. It is asserted
    // below rather than trusted, because every absolute yield in deg^-2 downstream is that
    // sum in disguise, and an area bookkeeping error does not look like an error -- it looks
    // like a survey that found more events than it did.
    //
    // WHAT STRATIFICATION DOES AND DOES NOT BIAS. Nothing computed *at* a sightline changes:
    // detection efficiency, per-event precision, and every ratio conditional on the sample
    // are untouched, because which sightlines were visited is not an input to any of them.
    // What changes is any quantity POOLED ACROSS sightlines -- a survey-wide yield, a
    // histogram over all events, the "all joint detections" panel of F4. Those must weight
    // each event by its sightline's `w_area` or they will describe a sky that is 2.6% Roman
    // by area and (say) 40% Roman by sample. The weight is written into every row of the
    // per-event table for exactly this reason.
    //
    // With --stride-roman absent, kSub == 1: the fine grid IS the coarse grid, every block is
    // a single cell that is its own representative, every area is gridStep^2, and the list is
    // the same points in the same order as the old nested loop -- so the RNG stream, and
    // therefore the run, is bit-identical to before this step.
    // ----------------------------------------------------------------------
    // The fine stratum (Deviation 69): every fine cell that overlaps the detector outline's
    // bounding rectangle of any placement. A grid point stands for the cell extending from it
    // in +l and +b, so a cell overlaps when its point lies within one fine step below the
    // rectangle; the same margin is kept on the other side. Chip-gap cells are therefore in the
    // fine stratum too: which sightlines actually see Roman is decided by the detector test in
    // matchVisibleEpochs, not here.
    auto insideFootprint = [&](double lon, double lat) {
        for (const auto& f : romanFields) {
            const double dl = lon - f.l, db = lat - f.b;
            const int k = f.layout;
            if (dl >= gl.dlMin[k] - fineStep and dl <= gl.dlMax[k] + fineStep
                and db >= gl.dbMin[k] - fineStep and db <= gl.dbMax[k] + fineStep) return true;
        }
        return false;
    };

    const int nLonGrid = int(std::floor((lonMax - lonMin) / gridStep + 1e-9)) + 1;
    const int nLatGrid = int(std::floor((latMax - latMin) / gridStep + 1e-9)) + 1;
    // Each coarse cell is subdivided into exactly kSub x kSub fine cells, so the fine grid has
    // nLonGrid*kSub columns -- not (nLonGrid-1)*kSub+1. The difference is the edge convention:
    // a grid POINT stands for the cell extending from it, which is what makes the areas tile
    // exactly and the sum below come out on the nose.
    const int    nLonFine = nLonGrid * kSub;
    const int    nLatFine = nLatGrid * kSub;
    const double cellArea = (gridStep * gridStep) / double(kSub * kSub); // deg^2, one fine cell

    // Roman coverage class of a sky point (Deviation 69): bit 0 = on a detector in the spring
    // roll, bit 1 = in the autumn roll. 0 none, 1 spring only, 2 autumn only, 3 both.
    auto romanClassAt = [&](double lon, double lat) {
        int c = 0;
        for (const auto& p : romanFields)
            if (inDetector(gl, p.layout, lon - p.l, lat - p.b)) c |= (1 << p.layout);
        return c;
    };

    // Pass 1: for each coarse block, how many of its fine cells survive the corner cut and lie
    // OUTSIDE the footprint, and which of them represents that area.
    const size_t nBlocks = size_t(nLonGrid) * size_t(nLatGrid);
    std::vector<int>  blockOutCount(nBlocks, 0);
    std::vector<long> blockRep(nBlocks, -1);      // fine cell index (iF*nLatFine + jF)
    long nFineKept = 0;
    for (int iF = 0; iF < nLonFine; ++iF) {
        const double lon = lonMin + iF * fineStep;
        for (int jF = 0; jF < nLatFine; ++jF) {
            const double lat = latMin + jF * fineStep;
            if (!inScan(lon, lat)) continue;       // same region test as the scan below
            nFineKept += 1;
            if (insideFootprint(lon, lat)) continue;
            const size_t b = size_t(iF / kSub) * size_t(nLatGrid) + size_t(jF / kSub);
            if (blockRep[b] < 0) blockRep[b] = long(iF) * nLatFine + jF;
            blockOutCount[b] += 1;
        }
    }

    // Pass 2: the list itself, in the same iLon-major order the old loop walked.
    std::vector<Sightline> scan;
    scan.reserve(size_t(nFineKept));
    std::vector<int> fieldHits(romanFields.size(), 0);
    long nSightlines = 0, nSightlinesRoman = 0;
    double areaFootprint = 0.0, areaOutside = 0.0;
    // Point-sampled sky area on a detector, per roll: the grid's estimate of what Roman images
    // in a spring / an autumn season, to compare with the exact detector area.
    std::array<double, GBTDS_NLAYOUT> areaOnDetector{};
    for (int iF = 0; iF < nLonFine; ++iF) {
        const double lon = lonMin + iF * fineStep;
        for (int jF = 0; jF < nLatFine; ++jF) {
            const double lat = latMin + jF * fineStep;
            if (!inScan(lon, lat)) continue;
            const bool anyField = insideFootprint(lon, lat);
            if (anyField) {
                std::array<bool, GBTDS_NLAYOUT> onDet{};
                for (size_t f = 0; f < romanFields.size(); ++f) {
                    const auto& p = romanFields[f];
                    if (inDetector(gl, p.layout, lon - p.l, lat - p.b)) {
                        fieldHits[f] += 1; onDet[p.layout] = true;
                    }
                }
                for (int k = 0; k < GBTDS_NLAYOUT; ++k) if (onDet[k]) areaOnDetector[k] += cellArea;
                scan.push_back({lon, lat, cellArea, true, iF, int(onDet[0]) | (int(onDet[1]) << 1)});
                areaFootprint    += cellArea;
                nSightlinesRoman += 1;
            } else {
                const size_t b = size_t(iF / kSub) * size_t(nLatGrid) + size_t(jF / kSub);
                if (blockRep[b] != long(iF) * nLatFine + jF) continue;  // not this block's rep
                const double a = blockOutCount[b] * cellArea;
                scan.push_back({lon, lat, a, false, iF, 0});
                areaOutside += a;
            }
            nSightlines += 1;
        }
    }

    // ----------------------------------------------------------------------
    // Post-stratifying the footprint's area weights (Deviation 69).
    //
    // A footprint sightline is a POINT: it either lands on a detector in a given roll or in a
    // chip gap. The detectors are 0.125 deg across with gaps of 0.008-0.026 deg, so a grid of
    // 0.1 deg aliases against them: at --stride-roman 5 the grid put 1.46 deg^2 on a detector per
    // roll where the exact figure is 1.68 -- every Roman yield per deg^2 13% low -- and even a
    // 0.02 deg grid lands 2% high. Rather than chase that with CPU, each footprint sightline's
    // area is rescaled so that each coverage CLASS (none / spring only / autumn only / both)
    // carries its exact sky area inside the footprint stratum:
    //     area_i <- cellArea * exact(c_i) / grid(c_i)
    // exact(c) from sub-sampling every fine cell on a ~0.002 deg raster. The classes partition
    // the stratum, so the stratum's total area -- and the invariant below -- is unchanged. The
    // sightlines of a class still sample its sky at the grid points; only how much sky each
    // stands for changes.
    // ----------------------------------------------------------------------
    std::array<double, 4> classGrid{}, classExact{};
    std::array<long, 4>   classN{};
    for (const auto& sl : scan)
        if (sl.inFootprint) { classGrid[sl.romanClass] += cellArea; classN[sl.romanClass] += 1; }
    {
        const int nSub = std::max(1, int(std::ceil(fineStep / 0.002 - 1e-9)));
        const double h = fineStep / nSub;
        for (const auto& sl : scan) {
            if (!sl.inFootprint) continue;
            for (int a = 0; a < nSub; ++a)
                for (int c = 0; c < nSub; ++c)
                    classExact[romanClassAt(sl.lon + (a + 0.5) * h, sl.lat + (c + 0.5) * h)] += h * h;
        }
    }
    for (int c = 0; c < 4; ++c) {
        if (classExact[c] > 0.0 and classN[c] == 0) {
            std::cerr << "ERROR: Roman coverage class " << c << " has " << classExact[c]
                      << " deg^2 of sky in the footprint but no sightline samples it. Use a finer "
                      << "--stride-roman.\n";
            return 1;
        }
    }
    areaFootprint = 0.0;
    for (auto& sl : scan)
        if (sl.inFootprint) {
            sl.area = cellArea * classExact[sl.romanClass] / classGrid[sl.romanClass];
            areaFootprint += sl.area;
        }
    std::cout << "Roman coverage classes in the footprint stratum (none / spring / autumn / both):"
              << "\n  sightlines  " << classN[0] << " / " << classN[1] << " / " << classN[2]
              << " / " << classN[3]
              << "\n  grid area   " << classGrid[0] << " / " << classGrid[1] << " / "
              << classGrid[2] << " / " << classGrid[3] << " deg^2"
              << "\n  exact area  " << classExact[0] << " / " << classExact[1] << " / "
              << classExact[2] << " / " << classExact[3] << " deg^2 (weights rescaled to these)"
              << std::endl;

    // The invariant. Every fine cell inside the scan region is represented exactly once,
    // either by itself (footprint) or by its block's representative (outside).
    const double areaScanned = areaFootprint + areaOutside;
    {
        const double areaExpect = nFineKept * cellArea;
        if (std::fabs(areaScanned - areaExpect) > 1e-9 * std::max(1.0, areaExpect)) {
            std::cerr << "ERROR: sightline area weights sum to " << areaScanned
                      << " deg^2 but the scanned grid is " << areaExpect << " deg^2. "
                      << "Every absolute yield downstream is this sum; refusing to run.\n";
            return 1;
        }
        std::cout << "Sightline grid: " << scan.size() << " sightlines ("
                  << nSightlinesRoman << " inside Roman's footprint at " << fineStep
                  << " deg, " << (scan.size() - size_t(nSightlinesRoman)) << " outside at "
                  << gridStep << " deg), covering " << areaScanned << " deg^2 ("
                  << areaFootprint << " footprint + " << areaOutside << " outside)." << std::endl;
        // The grid samples the detector mosaic at points, so its on-detector area is an
        // estimate; the exact value is 6 fields x 18 detectors. Their difference is the
        // sampling error of every per-deg^2 Roman yield from this grid.
        for (int k = 0; k < GBTDS_NLAYOUT; ++k) {
            double exact = 0.0;
            for (const ScaRect& r : gl.sca[k]) exact += (r.l1 - r.l0) * (r.b1 - r.b0);
            int nPlaced = 0;
            for (const auto& f : romanFields) nPlaced += (f.layout == k);
            exact *= nPlaced;
            std::cout << "Roman detector area, " << (k == 0 ? "spring" : "autumn") << " roll: "
                      << areaOnDetector[k] << " deg^2 on the grid vs " << exact
                      << " deg^2 exact (" << nPlaced << " fields)." << std::endl;
        }
    }

    const size_t nFieldsCovered = std::count_if(fieldHits.begin(), fieldHits.end(),
                                                [](int h){ return h > 0; });
    // Fatal for a full-region scan, advisory for --stub: the stub patch is 0.1x0.1 deg
    // by design and cannot possibly reach all twelve placements.
    if (nFieldsCovered < romanFields.size() and not cfg.stubPatch) {
        std::cerr << "ERROR: the sightline grid (stride " << cfg.stride << ", step "
                  << gridStep << " deg) misses " << (romanFields.size() - nFieldsCovered)
                  << " of " << romanFields.size() << " Roman field placements:\n";
        for (size_t f = 0; f < romanFields.size(); ++f)
            if (fieldHits[f] == 0)
                std::cerr << "    field at (l,b) = (" << romanFields[f].l << ", "
                          << romanFields[f].b << "), " << (romanFields[f].layout ? "autumn" : "spring")
                          << " roll, has no sightline on a detector\n";
        std::cerr << "Those fields would contribute no Roman epochs, and the run would "
                     "report joint columns built from Rubin data alone. Use a smaller "
                     "--stride.\n";
        return 1;
    }
    if (nFieldsCovered < romanFields.size() and cfg.stubPatch) {
        std::cout << "NOTE: --stub reaches " << nFieldsCovered << " of "
                  << romanFields.size() << " Roman field placements. Expected for a patch this "
                  << "small; joint statistics from it describe those fields only.\n";
    }


    grid.scan = std::move(scan);
    grid.romanFields = std::move(romanFields);
    grid.scanReach = scanReach;
    grid.lonMin = lonMin; grid.lonMax = lonMax; grid.latMin = latMin; grid.latMax = latMax;
    grid.cellArea = cellArea;
    grid.nSightlines = nSightlines; grid.nSightlinesRoman = nSightlinesRoman;
    grid.areaFootprint = areaFootprint; grid.areaOutside = areaOutside; grid.areaScanned = areaScanned;
    grid.areaOnDetector = areaOnDetector;
    grid.classGrid = classGrid; grid.classExact = classExact;
    grid.nFieldsCovered = nFieldsCovered;
    return 0;
}

void printDryRunStrata(const SightlineGrid& grid, const GridSteps& steps) {
    const auto& scan = grid.scan;
    const long   nSightlinesRoman = grid.nSightlinesRoman;
    const double areaFootprint = grid.areaFootprint, areaOutside = grid.areaOutside;
    const double areaScanned = grid.areaScanned, cellArea = grid.cellArea;
    const double fineStep = steps.fineStep, gridStep = steps.gridStep;
        // Per-stratum sightline counts, and what they cost. The footprint sightlines are the
        // ones worth buying and the ones that are expensive; printing both together is what
        // makes --stride-roman a decision rather than a guess.
        const long nOutside = long(scan.size()) - nSightlinesRoman;
        std::cout << "\n---- dry run: sightline strata ----\n"
                  << "  footprint  " << nSightlinesRoman << " sightlines at " << fineStep
                  << " deg, " << areaFootprint << " deg^2 ("
                  << (areaFootprint > 0.0 ? cellArea : 0.0) << " deg^2 each)\n"
                  << "  outside    " << nOutside << " sightlines at " << gridStep
                  << " deg, " << areaOutside << " deg^2\n"
                  << "  total      " << scan.size() << " sightlines, " << areaScanned
                  << " deg^2\n"
                  << "  footprint share: " << (100.0 * double(nSightlinesRoman) / double(scan.size()))
                  << "% of sightlines, " << (100.0 * areaFootprint / areaScanned)
                  << "% of area\n"
                  << "No stars drawn. Remove --dry-run to run the simulation.\n";
}
