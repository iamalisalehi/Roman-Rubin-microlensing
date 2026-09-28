# Sample-event dump spec -- ns, BATCH b (2026-09-25): the rubin_only class, which batch a left
# EMPTY through sightline 1664. Same reason and same remedy as bulge batch c
# (runs/samples_bulge_S3c): with the extinction fixed, when Roman sees the peak Roman detects the
# event, so rubin_only (Roman epochs within +-2 tE, no Roman detection) comes only from gap events
# whose wings reach a Roman season edge, and batch a's gap_filler -- listed first -- took those.
# 100 detections per sightline (batch a: 30). Batch a's events are kept.

dir        samples/ns
step_te    0.01
span_te    3.0
dt_coarse  2.0

class rubin_only    2
