#!/usr/bin/env python
"""座標保持のデバッグ用スクリプト"""

from rdkit import Chem
import numpy as np
from inverse_msmd.utils.mol_utils import read_mol_from_pdb_smi

# E24を読み込み
e24_mol = read_mol_from_pdb_smi(
    "data/sample_probes/E24.pdb",
    "data/sample_probes/E24.smi"
)

print("=== E24分子（水素あり）===")
print(f"原子数: {e24_mol.GetNumAtoms()}")
coords_with_h = e24_mol.GetConformer().GetPositions()
print(f"座標形状: {coords_with_h.shape}")
print("\n最初の5原子の座標:")
for i in range(min(5, len(coords_with_h))):
    atom = e24_mol.GetAtomWithIdx(i)
    print(f"  原子 {i} ({atom.GetSymbol()}): {coords_with_h[i]}")

# 水素を除去
e24_no_h = Chem.RemoveHs(e24_mol)

print("\n=== E24分子（水素なし）===")
print(f"原子数: {e24_no_h.GetNumAtoms()}")

# コンフォーマーがあるか確認
if e24_no_h.GetNumConformers() > 0:
    coords_no_h = e24_no_h.GetConformer().GetPositions()
    print(f"座標形状: {coords_no_h.shape}")
    print("\n最初の5原子の座標:")
    for i in range(min(5, len(coords_no_h))):
        atom = e24_no_h.GetAtomWithIdx(i)
        print(f"  原子 {i} ({atom.GetSymbol()}): {coords_no_h[i]}")
else:
    print("警告: 水素除去後にコンフォーマーが失われました！")
    
# 手動で座標をマッピング
print("\n=== 手動マッピング ===")
h_to_no_h = {}
no_h_idx = 0
for atom in e24_mol.GetAtoms():
    if atom.GetAtomicNum() != 1:  # 水素でない
        h_to_no_h[atom.GetIdx()] = no_h_idx
        no_h_idx += 1

print(f"マッピング数: {len(h_to_no_h)}")
print("水素なし原子の座標:")
for h_idx in sorted(h_to_no_h.keys(), key=lambda x: h_to_no_h[x])[:5]:
    no_h_idx = h_to_no_h[h_idx]
    atom = e24_mol.GetAtomWithIdx(h_idx)
    print(f"  元のIdx={h_idx} -> 新Idx={no_h_idx} ({atom.GetSymbol()}): {coords_with_h[h_idx]}")