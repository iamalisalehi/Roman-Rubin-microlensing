# Sample-event dump spec -- black-hole lens population.
#   ./roman --population bh --dump-samples samples/bh.spec ...
#
# CLASS ORDER MATTERS. An event is written once, under the FIRST class whose quota is
# still open, so the rare classes must come first. `bh_long`/`bh_short` are just
# "the joint fit detected it" with a tE window, and would swallow everything if they
# were listed above the selective ones -- the v3 bh run detected 110,144 events, of
# which only 5.2% were seen by both telescopes and 5.3% by Roman alone.

dir        samples/bh
step_te    0.01      # peak-window step as a fraction of tE: ~600 points across the peak
span_te    3.0       # peak window is t0 +- 3 tE
dt_coarse  2.0       # 2-day steps over the rest of the decade, to carry the parallax wobble

# name          quota  cuts
class astrometric   2  shift_min=1.1   # centroid shift clears Roman's 1.1 mas floor AND
                                       #   Roman's astrometric matrix actually inverted
class roman_only    2                  # Roman alone, with Rubin epochs on the peak
class both          3                  # where the joint fit earns its keep
class gap_filler    2                  # peak in a MID-MISSION Roman gap; Rubin alone caught it
class rubin_only    2                  # Rubin alone, with Roman epochs on the peak
class bh_long       2  te_min=500      # upper tail (median tE for bh is 332 d)
class bh_short      2  te_max=100      # lower tail
