import numpy as np
from astropy.coordinates import SkyCoord
from dustmaps.bayestar import BayestarQuery
from dustmaps.decaps import DECaPSQueryLite

bay    = BayestarQuery(max_samples=1)
decaps = DECaPSQueryLite(mean_only=True)   # E(B-V) directly, lighter memory footprint

RV = 3.1  # standard total-to-selective extinction ratio

l_list, b_list, _ = np.loadtxt('./BulgeBaseline.dat', usecols=[3, 4, 13], unpack=True)
nfiles = len(l_list)              # must match `nfiles` in Bulge.h
nlines = 3686                     # must match `nlines` in Bulge.h
d_max  = 20.0                     # kpc — matches MaxD in Bulge.h
dist_grid = np.linspace(0.05, d_max, nlines)  # already sorted ascending, no extra sort needed

coords_check = SkyCoord(l_list, b_list, frame='galactic', unit='deg')
frac_south = (coords_check.icrs.dec.deg < -30.0).mean()
print(f"{frac_south*100:.1f}% of sightlines fall south of dec=-30 (need DECaPS)")

for i in range(nfiles):
    l = np.full(nlines, l_list[i])
    b = np.full(nlines, b_list[i])
    coords = SkyCoord(l, b, distance=dist_grid, frame='galactic', unit=('deg', 'deg', 'kpc'))

    dec = coords.icrs.dec.deg[0]  # same for every point on one sightline
    if dec > -30.0:
        ebv = bay(coords, mode='mean') * 0.884   # Bayestar19 -> E(B-V)
    else:
        ebv = decaps(coords, mode='mean')        # already E(B-V)

    av = RV * ebv
    out = np.stack((l, b, dist_grid, av), axis=-1)
    fname = f'./files/ext/bayestar_{l_list[i]:.2f}_{b_list[i]:.2f}.txt'
    np.savetxt(fname, out, fmt="%.3e\t%.3e\t%.3e\t%.3e")
    print(f'{fname} written')
