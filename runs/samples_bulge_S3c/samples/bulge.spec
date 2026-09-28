# Sample-event dump spec -- bulge, BATCH c (2026-09-25): the rubin_only class, which batch a
# left EMPTY over the whole scan.
#
# Why it is rare, measured: among footprint events Rubin detects and Roman does not, almost none
# peak inside a Roman season (production bulge table: 1 of 917; batch a: 0 of 119) -- with the
# extinction fixed, when Roman sees the peak Roman detects it. rubin_only needs Roman epochs
# within +-2 tE of the peak AND no Roman detection, so it can only come from gap events whose
# wings reach a season edge, and in batch a the earlier-listed gap_filler class took those. This
# batch asks for rubin_only alone, with 100 detections per sightline (batch a: 30) to see ~3x as
# many candidates. An empty result here is itself the finding. Batch a's events are kept.

dir        samples/bulge
step_te    0.01
span_te    3.0
dt_coarse  2.0

class rubin_only    2
