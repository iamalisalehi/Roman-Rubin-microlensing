# V1: rebuilt extinction tables vs VVV (analysis/v1_ext_vvv.py)

Tables: `files/ext/ext_tables.dat`; 2013 scan sightlines (production grid); E(J-Ks)/A_V = 0.1347 on the five-field block.

(a) A_V(8 kpc) tables / VVV per |b| bin (target ~0.9-1.0; old tables 0.17 at |b|<0.5)

```
                   n  vvv_av  tables_av  ratio  ratio_p16  ratio_p84
bbin                                                                
(-0.001, 0.5]  278.0  23.951     23.904  0.988      0.871      1.170
(0.5, 1.0]     310.0  13.198     13.564  0.979      0.857      1.137
(1.0, 1.5]     276.0   7.111      6.912  0.985      0.907      1.081
(1.5, 2.5]     464.0   5.385      4.891  0.905      0.829      1.003
(2.5, 6.0]     685.0   2.457      2.257  0.883      0.829      0.981
```

(b) GC field / five-field contrast: tables 4.74, VVV 5.52 (target ~4.5-4.9; old tables 0.57)
