# Sample-event dump spec -- neutron-star lens population.
#   ./roman --population ns --dump-samples samples/ns.spec ...
#
# The ns astrometric signal is small on purpose: the v3 run's median centroid shift is
# 0.265 mas against Roman's 1.1 mas floor, and only 0.57% of events clear it. The
# `astrometric` class here therefore has NO shift_min cut -- it asks only that Roman's
# astrometric matrix inverted. If it comes back empty, that is the physics, not a bug.

dir        samples/ns
step_te    0.01
span_te    3.0
dt_coarse  2.0

class astrometric   2
class roman_only    2
class both          3
class gap_filler    2
class rubin_only    2
class ns_typical    3                  # the reference case (median tE for ns is 64 d)
