import numpy as np
from astropy.coordinates import SkyCoord
from dustmaps.bayestar import BayestarQuery
from dustmaps.decaps import DECaPSQueryLite

bay    = BayestarQuery(max_samples=1)
decaps = DECaPSQueryLite(mean_only=True)

RV_BAYESTAR = 3.1    # standard total-to-selective extinction ratio
RV_DECAPS   = 3.32   # dustmaps docs specify this value for the DECaPS E(B-V) calibration

l_raw, b_raw, _ = np.loadtxt('./Baseline/BulgeBaseline.dat', usecols=[3, 4, 13], unpack=True)

# Deduplicate: BulgeBaseline.dat is a per-visit file (the same field is revisited
# constantly over the 10-year baseline), so without this we'd query the dust map
# once per visit instead of once per unique sightline.
coords_rounded = np.round(np.column_stack([l_raw, b_raw]), 2)
unique_coords  = np.unique(coords_rounded, axis=0)
l_list, b_list = unique_coords[:, 0], unique_coords[:, 1]

nfiles = len(l_list)
print(f"{nfiles} unique sightlines from {len(l_raw)} visits -- update NFILES in Bulge.h to match")

d_max     = 20.0
nlines    = 3686
dist_grid = np.linspace(0.05, d_max, nlines)

for i in range(nfiles):
    l = np.full(nlines, l_list[i])
    b = np.full(nlines, b_list[i])
    coords = SkyCoord(l, b, distance=dist_grid, frame='galactic', unit=('deg', 'deg', 'kpc'))
    dec = coords.icrs.dec.deg[0]

    if dec > -30.0:
        ebv, flags = bay(coords, mode='mean', return_flags=True)
        ebv = ebv * 0.884
        rv = RV_BAYESTAR
        # Bayestar's own reliable-distance flag: samples beyond what this
        # sightline's photometry actually constrains come back NaN already,
        # but the flag lets you tell "no data at all" apart from "extrapolated."
        reliable_dist = flags['reliable_dist'] if 'reliable_dist' in flags.dtype.names else None
    else:
        ebv, flags = decaps(coords, mode='mean', return_flags=True)
        rv = RV_DECAPS
        reliable_dist = None

    av = rv * ebv

    # Repair NaN samples: interior gaps get filled by interpolating along this
    # sightline's own valid distance samples (extinction is monotonically
    # non-decreasing with distance, so this is physically well-justified).
    # A sightline that's entirely NaN is handled in a second pass below.
    nan_mask = np.isnan(av)
    is_total_dropout = False
    if nan_mask.any():
        valid = ~nan_mask
        if valid.sum() >= 2:
            av[nan_mask] = np.interp(dist_grid[nan_mask], dist_grid[valid], av[valid])
        elif valid.sum() == 1:
            av[nan_mask] = av[valid][0]
        else:
            is_total_dropout = True

    out = np.stack((l, b, dist_grid, av), axis=-1)
    fname = f'./files/ext/bayestar_{l_list[i]:.2f}_{b_list[i]:.2f}.txt'
    np.savetxt(fname, out)

    tag = " [TOTAL DROPOUT -- needs neighbor fallback]" if is_total_dropout else ""
    print(f'{fname} written{tag}')
