# Roman-LSST Microlensing Simulation

This project contains tools for simulating and studying microlensing events detectable by the Nancy Grace Roman Space Telescope and the Vera C. Rubin Observatory Legacy Survey of Space and Time (LSST).

The code focuses on Galactic bulge microlensing simulations, including:
- generation and handling of stellar populations and color-magnitude diagrams (CMDs),
- extinction and dust-map handling,
- stellar photometry conversions,
- Roman and LSST observing simulations,
- microlensing event modeling and analysis.

## Structure

```

.
├── Baseline/      # Survey baseline data and configurations
├── CMD/           # Color-magnitude diagram related files
├── dustmaps/      # Dust/extinction map utilities
├── files/         # Auxiliary files
├── pics/          # Generated figures
├── *.cpp, *.h     # C++ simulation code
├── *.py           # Python analysis and plotting scripts
└── Makefile       # Build instructions

````

## Requirements

The project uses:

- C++ compiler with C++11 or newer support
- Python 3
- Scientific Python packages (NumPy, SciPy, Astropy, Matplotlib, etc.)

Additional dependencies may be required depending on the specific simulation modules.

## Building

Compile the C++ code using:

```bash
make
````

## Acknowledgements

This project builds upon previous microlensing simulation work developed by Prof. Sedighe Sajadian. I thank her for her guidance and for the foundation provided by the original codebase.

## Notes

Large simulation datasets, catalogs, and generated files are intentionally not included in the repository. They should be obtained separately or generated using the provided tools.
