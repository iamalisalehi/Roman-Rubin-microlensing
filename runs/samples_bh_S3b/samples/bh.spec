# Sample-event dump spec -- bh, BATCH b (2026-09-25). Batch a ran the whole scan and left BOTH
# gap classes empty (gap_filler 0/2, rubin_only 0/2), and its two bh_short events are weak
# (u0 > 1, sigma_tE 100% and 160%). Black-hole events are long (median tE ~330 d), so one rarely
# fits inside a ~110-d Roman gap, and one that overlaps a season Roman detects from the wings.
# Both gap classes are therefore restricted to tE <= 150 d, the only events that CAN fill them.
# rubin_only is listed first so gap_filler cannot take its events (bulge batch a lost them that
# way). 100 detections per sightline (batch a: 30). If a gap class is still empty after this
# pass, that is the finding. Batch a's events are kept.

dir        samples/bh
step_te    0.01
span_te    3.0
dt_coarse  2.0

class rubin_only    2  te_max=150
class gap_filler    2  te_max=150
class bh_short      2  te_max=100
