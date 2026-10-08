"""Export LSSTCam's active-silicon map for src/surveys/footprints.cpp (readLsstCamMap).

Source: the map OpSim/MAF itself use, rubin_scheduler.utils.LsstCameraFootprint's default
`fov_map.npz`, from the rubin_sim_data bundle utils_2023_11_02.tgz
(https://s3df.slac.stanford.edu/data/rubin/sim-data/rubin_sim_data/). A 1000 x 1000 boolean image
of active science silicon on the focal plane, x = y = -1.75..1.75 deg in the gnomonic tangent plane,
indexed image[ix][iy]; rubin_scheduler masks pixels beyond max_radius = 1.94 deg, and so does this script.

Output fov_map.txt: header lines '# n x0 step max_radius' then n rows of n characters '0'/'1', row
ix, column iy. 9.12 deg^2 active.

    .roman/bin/python Baseline/lsstcam_fov/export_fov_map.py
"""
import io
import os
import tarfile
import urllib.request

import numpy as np

URL = "https://s3df.slac.stanford.edu/data/rubin/sim-data/rubin_sim_data/utils_2023_11_02.tgz"
MAX_RADIUS = 1.94
here = os.path.dirname(os.path.abspath(__file__))

raw = urllib.request.urlopen(URL, timeout=120).read()
with tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz") as t:
    z = np.load(io.BytesIO(t.extractfile("utils/fov_map.npz").read()))
    image, x = z["image"].copy(), z["x"].copy()
n = len(x)
step = (x[-1] - x[0]) / (n - 1)
X = np.ones(image.shape) * x                 # as in LsstCameraFootprint: full_x varies along axis 1
Y = np.ones(image.shape) * x[np.newaxis].T
image[np.hypot(X, Y) > MAX_RADIUS] = False
with open(os.path.join(here, "fov_map.txt"), "w") as f:
    f.write(f"# LSSTCam active silicon, rubin_sim_data utils_2023_11_02 fov_map.npz, cropped at {MAX_RADIUS} deg\n")
    f.write(f"# {n} {x[0]:.12f} {step:.15f} {MAX_RADIUS}\n")
    for row in image:
        f.write("".join("1" if v else "0" for v in row) + "\n")
print(f"wrote fov_map.txt: {n}x{n}, x0 {x[0]}, step {step:.6f} deg, active {image.sum() * step**2:.3f} deg^2")
