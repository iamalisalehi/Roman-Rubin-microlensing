# Sample-event dump spec -- the default bulge population (Kroupa IMF + remnants).
#   ./roman --population bulge --dump-samples samples/bulge.spec ...
#
# These are the ordinary events, median tE ~15 d -- the regime where Roman's 12.1-minute
# cadence dominates inside a season and Rubin's decade of coverage dominates outside one.
# The short-tE gap-filling contrast is sharpest here, so gap_filler gets the largest quota.

dir        samples/bulge
step_te    0.01
span_te    3.0
dt_coarse  2.0

class astrometric   2
class roman_only    2
class gap_filler    3                  # thesis novelty claim 2, made drawable
class both          3
class rubin_only    2
class any           2                  # an unremarkable event, for scale
