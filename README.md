# inverse_msmd

Inverse analysis with mixed-solvent molecular dynamics (MSMD)

## Overview

`inverse_msmd` is a Python package for inverse analysis of mixed-solvent molecular dynamics (MSMD) simulations. It provides functionality for protein structure superimposition and matching score calculations based on interaction profiles.

## Features

- **Structure Superimposition**: BioPython-based structure superimposition with scikit-learn compatible interface
- **PDB Manipulation**: Read, get/set attributes, and save PDB files
- **Spatial Calculations**: Utilities for spatial calculations such as volume estimation of spheres
- **Path Handling**: Path processing with environment variable and tilde expansion support

## Installation

### Dependencies

- Python >= 3.8
- numpy >= 1.20.0
- biopython >= 1.79
- scikit-learn >= 0.24.0
- scipy >= 1.7.0
- gridData >= 0.6.0

### Installation Methods

#### Development Mode (Recommended)

```bash
# Clone the repository
git clone https://github.com/akiyamalab/inverse_msmd.git
cd inverse_msmd

# Install in development mode
pip install -e .
```

#### Normal Installation

```bash
pip install -r requirements.txt
pip install .
```

## Usage

### Basic Example

```python
from inverse_msmd import SuperImposer, PDB
import numpy as np

# Read PDB files
protein = PDB.get_structure("protein.pdb")
probe = PDB.get_structure("probe.pdb")

# Get coordinates
protein_coords = PDB.get_attr(protein, "coord")
probe_coords = PDB.get_attr(probe, "coord")

# Select coordinates based on atom pairs
atom_pairs = np.loadtxt("atom_matching.txt", int)
probe_coords_target = probe_coords[atom_pairs[0]]
protein_coords_target = protein_coords[atom_pairs[1]]

# Superimpose structures
si = SuperImposer()
si.fit(protein_coords_target, probe_coords_target)

# Set transformed coordinates
transformed_coords = si.transform(protein_coords)
PDB.set_attr(protein, "coord", transformed_coords)

# Save results
PDB.save(protein, "aligned_protein.pdb")
```

### Example Scripts

The package includes example scripts in the `scripts/` directory:

- `scripts/superimposition.py`: Example of structure superimposition
- `scripts/calculate_matching.py`: Example of matching score calculation

## Directory Structure

```
inverse_msmd/
├── pyproject.toml          # Package configuration
├── README.md               # This file
├── LICENSE                 # License file
├── requirements.txt        # Dependency list
├── inverse_msmd/           # Main package
│   ├── __init__.py        # Public API
│   └── utils/             # Utility modules
│       ├── __init__.py
│       ├── bio_utils.py   # BioPython utilities
│       ├── spatial_utils.py  # Spatial calculation utilities
│       └── path_utils.py  # Path handling utilities
├── scripts/                # Example scripts
│   ├── superimposition.py
│   └── calculate_matching.py
└── profile_superimposition/  # Original implementation (reference)
```

## API Reference

### SuperImposer

A class for structure superimposition (scikit-learn compatible interface)

**Methods:**
- `fit(coords, reference_coords)`: Calculate superimposition parameters
- `transform(coords)`: Transform coordinates
- `inverse_transform(coords)`: Perform inverse transformation

### PDB

Utility functions for PDB structure manipulation

**Main Functions:**
- `get_structure(filepath)`: Read a PDB file
- `get_attr(model, attr, sele=None)`: Get attributes
- `set_attr(model, attr, lst, sele=None)`: Set attributes
- `save(structs, path)`: Save structure(s) to a PDB file

### Other Utilities

- `estimate_volume(points, radii, granularity=10)`: Estimate the volume of a set of spheres
- `expandpath(path)`: Return a path with environment variables and tilde expanded

## License

MIT License

## Author

Keisuke Yanagisawa (yanagisawa@comp.isct.ac.jp)

## Citation

If you use this package, please provide appropriate citation.

## Contributing

Bug reports and feature requests are welcome via GitHub Issues.
