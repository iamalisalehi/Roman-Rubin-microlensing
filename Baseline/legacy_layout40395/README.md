# Visit lists of the runs made before Deviation 69 (archived 2026-10-01)

The two baselines every production run up to and including `runs/prod_{bulge,bh,ns}_20260924`
read: Roman on the notional `layout_40395` (six 0.3003-deg circles, no rolls, mission start sim
day 730 = 2028-04-10; 7 columns, 302,406 rows) and Rubin's pointings cut by the old l/b box
(3,686 rows). Kept so that analyses of those runs (e.g. `analysis/u2_resolution_depth.py`,
`analysis/f1_results_table.py`) can still be reproduced: pass these files explicitly.
`analysis/gbtds_geometry.read_roman_visits()` recognises the 7-column format and applies the
old circle test. The live files one directory up describe the adopted layout.
