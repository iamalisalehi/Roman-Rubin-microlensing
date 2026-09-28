# Sample-event dump spec -- bulge, BATCH b (2026-09-24): redraw the astrometric class only.
#
# Batch a's two astrometric events both peaked at ~0.27 mas of centroid shift, barely above the
# ordinary events, and one was a 1,030-day, u0 = 2.1 event peaking after Roman's mission; the
# best astrometric example batch a produced (0.84 mas, theta_E to 0.5%) landed in `both` because
# the astrometric quota was already full. This batch asks for what the class is meant to show:
# a centroid shift that stands out (>= 0.5 mas; bulge median ~0.12) on an event short enough to
# sit inside Roman's coverage. Batch a's events are kept (user, 2026-09-24).

dir        samples/bulge
step_te    0.01
span_te    3.0
dt_coarse  2.0

class astrometric   3  shift_min=0.5 te_max=300
