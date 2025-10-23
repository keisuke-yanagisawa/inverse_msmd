#!/usr/bin/env python
"""リガンド置換関数のデバッグ"""

from rdkit import Chem
import numpy as np
from inverse_msmd.substructure_replacement import (
    find_substructure_in_ligand,
    match_substructures
)
from inverse_msmd.utils.mol_utils import read_mol_from_pdb_smi

# テストデータ読み込み
ligand_mol = next(Chem.SDMolSupplier("data/atom_matching/4hw3_A_lig.sdf"))
e23_mol = read_mol_from_pdb_smi("data/sample_probes/E23.pdb", "data/sample_probes/E23.smi")
e24_mol = read_mol_from_pdb_smi("data/sample_probes/E24.pdb", "data/sample_probes/E24.smi")

# マッチング情報取得
matches = find_substructure_in_ligand(ligand_mol, e23_mol)
atom_pair_patterns = match_substructures(e23_mol, e24_mol)

match = matches[0]
atom_pairs = atom_pair_patterns[0]

print("=== 置換前の座標確認 ===")
ligand_no_h = Chem.RemoveHs(ligand_mol)
e24_no_h = Chem.RemoveHs(e24_mol)

ligand_coords = ligand_no_h.GetConformer().GetPositions()
e24_coords = e24_no_h.GetConformer().GetPositions()

print(f"リガンド座標（最初の3原子）:")
for i in range(3):
    print(f"  {i}: {ligand_coords[i]}")

print(f"\nE24座標（最初の3原子）:")
for i in range(3):
    print(f"  {i}: {e24_coords[i]}")

print(f"\nE24座標（全て）:")
for i in range(len(e24_coords)):
    print(f"  {i}: {e24_coords[i]}")

# 置換関数を呼び出す
from inverse_msmd.substructure_replacement import replace_ligand_substructure
print("\n=== replace_ligand_substructure実行 ===")
replaced = replace_ligand_substructure(ligand_mol, match, e24_mol, atom_pairs)

print(f"\n置換後の分子:")
print(f"  原子数: {replaced.GetNumAtoms()}")
print(f"  コンフォーマー数: {replaced.GetNumConformers()}")

if replaced.GetNumConformers() > 0:
    replaced_coords = replaced.GetConformer().GetPositions()
    print(f"  座標形状: {replaced_coords.shape}")
    print(f"\n置換後の座標（最後の12原子=E24部分）:")
    start_idx = replaced.GetNumAtoms() - 12
    for i in range(start_idx, replaced.GetNumAtoms()):
        print(f"  {i}: {replaced_coords[i]}")
else:
    print("  警告: コンフォーマーが存在しません！")