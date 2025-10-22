"""
Structure superimposition script.

This script demonstrates how to use inverse_msmd package to superimpose
protein structures based on atom matching data.
"""
from inverse_msmd.utils.bio_utils import SuperImposer, PDB
import numpy as np

protein_file = "../atom_auto_matching/4hw3_A.pdb"
matching_ids = "A08 A08_1 A08_2 A17 A57 A58 A58_1 B23 B39 E21 E21_1 E21_2 E22 E23 E24 E24_1 E24_2 E25 E25_1 E26 E27 E28 E28_1 E41 E42 E43 E44 E45 E46 E47 E48 E49 E50 E51 E52".split(" ")
probe_pdb = "/home/7/16D30270/workspace/9998_share/msmd/probe/{probe}/{probe}.pdb"

# print(np.loadtxt(f"atom_matching_{probe_id}", int)) # probe - refmol

protein = PDB.get_structure(f"{protein_file}")

for matching in matching_ids:
    cid = matching[:3]
    probe = PDB.get_structure(probe_pdb.format(probe=cid))
    import os
    os.system(f"cp {probe_pdb.format(probe=cid)} {cid}.pdb")


    protein_coords = PDB.get_attr(protein, "coord")
    probe_coords = PDB.get_attr(probe, "coord")

    atom_pairs = np.loadtxt(f"../atom_auto_matching/atom_matching_{matching}", int)
    atom_pairs += [[0],[1193]] # reference sdf to reference pdb

    probe_coords_target   = probe_coords[atom_pairs[0]]
    protein_coords_target = protein_coords[atom_pairs[1]]


    si = SuperImposer()
    si.fit(protein_coords_target, probe_coords_target)

    PDB.set_attr(protein, "coord", si.transform(protein_coords))
    PDB.save(protein, f"4hw3_aligned_to_{matching}.pdb")