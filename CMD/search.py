import numpy as np
import pandas as pd

def read_input(path):
    with open(path, "r") as f:
        for line in f:
            if line.startswith("#"):
                header_input = line.lstrip("#").split()
            else:
                break

    input_data = pd.DataFrame(np.loadtxt(path), columns=header_input)

    return input_data

df = read_input('./Besancon/bos6.dat')

for p in sorted(df["Pop"].unique()):
    s = df[df["Pop"] == p]
    print(
        p,
        s["Age"].mean(),
        s["[M/H]"].mean(),
        len(s)
    )
