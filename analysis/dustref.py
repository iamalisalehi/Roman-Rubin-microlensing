"""The reference dust: A_V(d) along a sightline from DECaPS where it can see, Marshall where it cannot.

The simulator's extinction tables were too thin within 1 deg of the plane when built from Bayestar
(optical, unreliable at large distance) and DECaPS (saturates toward the Galactic-centre field). The
hybrid profile below was validated against the independent VVV reddening map. maps.py builds the
simulator's tables from this module, so the tables and their checks share one definition.

Profile ("hybrid", per sightline, on DGRID):
  1. DECaPS (Zucker et al. 2025; dustmaps DECaPSQueryLite, mean), A_V = 3.32 E(B-V) (the dustmaps
     convention), out to the largest distance its `reliable_dist` flag accepts;
  2. beyond it, DECaPS's last reliable value plus Marshall's further increase, A_Ks / k;
  3. from the first distance at which Marshall's A_Ks / k reaches DECAPS_AV_MAX (DECaPS's stated
     sensitivity limit) on, Marshall's A_Ks / k itself;
  4. forced non-decreasing (a running maximum), since step 3 can step down where steps 1-2 had
     already passed Marshall's value and a column of dust cannot shrink with distance.
  Variants "decaps" (DECaPS everywhere, flag ignored) and "marshall" (A_Ks / k everywhere) give the
  systematic error. k = A_Ks/A_V calibrates the near-infrared map onto DECaPS's scale: the median
  A_Ks(Marshall)/A_V(DECaPS) at 8 kpc over Roman's five-field block where DECaPS is reliable
  (`calibrate_k`).

Both maps are queried through the dustmaps library (data in ./dustmaps, the repo-local data_dir).
"""

import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DGRID = np.arange(0.05, 20.0, 0.05)            # kpc; 399 distances (config/parameters.h MaxD = 12 kpc)
RV_DECAPS = 3.32
DECAPS_AV_MAX = 12.0
K_NOMINAL_LEGACY = 0.0805                      # calibrate_k on the notional layout (block at b -1.2)


def _coords(l, b, d):
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    return SkyCoord(l=np.asarray(l) * u.deg, b=np.asarray(b) * u.deg,
                    distance=np.asarray(d) * u.kpc, frame="galactic")


class ReferenceDust:
    """Queries the two maps and assembles the reference profile.

    Single-sightline methods cache by (l, b) rounded to 1e-3 deg; `batch` queries
    many sightlines at once, which is how maps.py builds the tables (DECaPS's cost is per distinct
    sky position, ~0.13 s each, so the batch is ordered spatially by the caller).
    """

    def __init__(self):
        from dustmaps.config import config
        config["data_dir"] = os.path.join(ROOT, "dustmaps")
        from dustmaps.decaps import DECaPSQueryLite
        from dustmaps.marshall import MarshallQuery
        self.mq = MarshallQuery()
        self.dq = DECaPSQueryLite(mean_only=True)
        self._aks, self._dec = {}, {}
        self.dsat = {}

    # ---- raw maps, many sightlines at once: arrays of shape (n, len(DGRID)) ----
    def batch(self, l, b):
        l, b = np.asarray(l, float), np.asarray(b, float)
        n, m = l.size, DGRID.size
        c = _coords(np.repeat(l, m), np.repeat(b, m), np.tile(DGRID, n))
        v, fl = self.dq(c, mode="mean", return_flags=True)
        dav = RV_DECAPS * np.asarray(v, float).reshape(n, m)
        rel = np.asarray(fl["reliable_dist"], bool).reshape(n, m)
        aks = np.asarray(self.mq(c), float).reshape(n, m)
        return dav, rel, np.array([_hold_aks(a) for a in aks])

    # ---- single sightline, cached ----
    def aks_profile(self, l, b):
        k = (round(l, 3), round(b, 3))
        if k not in self._aks:
            self._aks[k] = _hold_aks(np.asarray(self.mq(_coords(np.full(DGRID.size, l),
                                                                  np.full(DGRID.size, b), DGRID)), float))
        return self._aks[k]

    def decaps_profile(self, l, b):
        k = (round(l, 3), round(b, 3))
        if k not in self._dec:
            v, fl = self.dq(_coords(np.full(DGRID.size, l), np.full(DGRID.size, b), DGRID),
                            mode="mean", return_flags=True)
            self._dec[k] = (RV_DECAPS * np.asarray(v, float), np.asarray(fl["reliable_dist"], bool))
        return self._dec[k]

    def reference_profile(self, l, b, variant, aks_av, monotone=True):
        dav, rel = self.decaps_profile(l, b)
        prof, dsat = combine(dav, rel, self.aks_profile(l, b), variant, aks_av, monotone)
        if variant == "hybrid":
            self.dsat[(round(l, 3), round(b, 3), aks_av)] = dsat
        return prof


def _hold_aks(a):
    """Marshall A_Ks on DGRID; NaN beyond its coverage is held at the last valid value."""
    ok = np.isfinite(a)
    if ok.sum() >= 2:
        return np.interp(DGRID, DGRID[ok], a[ok])
    if ok.sum() == 1:
        return np.full(DGRID.size, a[ok][0])
    return np.full(DGRID.size, np.nan)


def combine(dav, rel, aks, variant, aks_av, monotone=True):
    """The reference profile from one sightline's raw maps. Returns (A_V on DGRID, d_switch),
    d_switch the distance from which Marshall is used outright (inf if never; nan if n/a)."""
    nir = aks / aks_av
    if variant == "marshall":
        prof, dsat = nir.copy(), np.nan
    elif not np.isfinite(dav).any():
        prof, dsat = nir.copy(), np.nan                  # outside DECaPS (not the case in this scan)
    else:
        dav = np.interp(DGRID, DGRID[np.isfinite(dav)], dav[np.isfinite(dav)])
        if variant == "decaps":
            prof, dsat = dav, np.nan
        elif np.all(np.isnan(aks)):
            prof, dsat = dav, np.nan
        elif not rel.any():
            prof, dsat = nir.copy(), np.nan
        else:
            imax = int(np.flatnonzero(rel).max())
            prof = dav.copy()
            beyond = np.arange(DGRID.size) > imax
            prof[beyond] = dav[imax] + np.maximum(aks[beyond] - aks[imax], 0.0) / aks_av
            hit = np.flatnonzero(nir >= DECAPS_AV_MAX)
            dsat = DGRID[hit[0]] if hit.size else np.inf
            if hit.size:
                prof[hit[0]:] = nir[hit[0]:]
    if monotone:
        prof = np.maximum.accumulate(prof)
    return prof, dsat


def calibrate_k(dust, l, b, d_kpc=8.0):
    """k = median A_Ks(Marshall)/A_V(DECaPS) at d_kpc over the given sightlines, using only those
    where DECaPS's flag calls d_kpc reliable. Returns (k, 16th, 84th percentile, n used)."""
    i = int(np.argmin(np.abs(DGRID - d_kpc)))
    dav, rel, aks = dust.batch(l, b)
    ok = rel[:, i] & np.isfinite(dav[:, i]) & np.isfinite(aks[:, i]) & (dav[:, i] > 0)
    r = aks[ok, i] / dav[ok, i]
    return float(np.median(r)), float(np.percentile(r, 16)), float(np.percentile(r, 84)), int(ok.sum())
