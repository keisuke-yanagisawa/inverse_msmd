#!/usr/bin/env python3
"""Atom Matching機能のテスト"""

from rdkit import Chem
import numpy as np
from inverse_msmd.substructure_replacement import (
    find_substructure_in_ligand,
    match_substructures
)
from inverse_msmd.utils.mol_utils import read_mol_from_pdb_smi

# テストデータの読み込み
ligand_mol = next(Chem.SDMolSupplier("data/atom_matching/4hw3_A_lig.sdf"))
e23_mol = read_mol_from_pdb_smi(
    "data/sample_probes/E23.pdb",
    "data/sample_probes/E23.smi"
)
e24_mol = read_mol_from_pdb_smi(
    "data/sample_probes/E24.pdb",
    "data/sample_probes/E24.smi"
)

print(f"E23原子数: {e23_mol.GetNumAtoms()}")
print(f"E24原子数: {e24_mol.GetNumAtoms()}")

# 水素を除去した後の原子数も確認
e23_no_h = Chem.RemoveHs(e23_mol)
e24_no_h = Chem.RemoveHs(e24_mol)
print(f"E23原子数（水素なし）: {e23_no_h.GetNumAtoms()}")
print(f"E24原子数（水素なし）: {e24_no_h.GetNumAtoms()}")

# リガンド中のE23を探索
ligand_e23_matches = find_substructure_in_ligand(ligand_mol, e23_mol)
print(f"\nリガンド中のE23マッチ数: {len(ligand_e23_matches)}")

# 最初のマッチを使用してテスト
first_match = ligand_e23_matches[0]
print(f"使用するマッチ: {first_match}")

# E23とE24のatom matchingを実行
print("\nAtom Matchingを実行中...")
atom_pair_patterns = match_substructures(e23_mol, e24_mol)

print(f"\n✓ Atom Matchingパターン数: {len(atom_pair_patterns)}")
for i, pairs in enumerate(atom_pair_patterns):
    print(f"\nパターン {i}:")
    print(f"  形状: {pairs.shape}")
    print(f"  E23側インデックス: {pairs[0]}")
    print(f"  E24側インデックス: {pairs[1]}")
    
    # 検証
    assert pairs.shape[0] == 2, "atom pairsは2行である必要があります"
    assert pairs.shape[1] > 0, "少なくとも1つの原子ペアが必要です"
    assert np.all(pairs[0] < e23_no_h.GetNumAtoms()), "E23のインデックスが範囲外"
    assert np.all(pairs[1] < e24_no_h.GetNumAtoms()), "E24のインデックスが範囲外"

print("\n✓ 全ての検証に合格")