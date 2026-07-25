import numpy as np
import glob

total_dropout, partial_gap = [], []
for fname in glob.glob('../files/ext/*.txt'):
    ext = np.loadtxt(fname)[:, 3]
    n_nan = np.isnan(ext).sum()
    if n_nan == len(ext):
        total_dropout.append(fname)
    elif n_nan > 0:
        partial_gap.append((fname, n_nan, len(ext)))

print(f"{len(total_dropout)} sightlines entirely NaN")
print(f"{len(partial_gap)} sightlines with isolated NaN gaps")
for f, n, tot in partial_gap[:10]:
    print(f"  {f}: {n}/{tot} NaN")
