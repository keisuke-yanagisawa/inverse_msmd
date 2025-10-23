#!/usr/bin/env python3
"""部分構造探索機能のテスト"""

from rdkit import Chem
from inverse_msmd.substructure_replacement import find_substructure_in_ligand
from inverse_msmd.utils.mol_utils import read_mol_from_pdb_smi

# テストデータの読み込み
print("テストデータを読み込み中...")
ligand_mol = next(Chem.SDMolSupplier("data/atom_matching/4hw3_A_lig.sdf"))
substructure_mol = read_mol_from_pdb_smi(
    "data/sample_probes/E23.pdb",
    "data/sample_probes/E23.smi"
)

print(f"リガンド原子数: {ligand_mol.GetNumAtoms()}")
print(f"部分構造原子数（水素あり）: {substructure_mol.GetNumAtoms()}")

# 水素を除去した後の原子数も確認
ligand_no_h = Chem.RemoveHs(ligand_mol)
substructure_no_h = Chem.RemoveHs(substructure_mol)
print(f"リガンド原子数（水素なし）: {ligand_no_h.GetNumAtoms()}")
print(f"部分構造原子数（水素なし）: {substructure_no_h.GetNumAtoms()}")

# 部分構造探索の実行
print("\n部分構造探索を実行中...")
matches = find_substructure_in_ligand(ligand_mol, substructure_mol)

# 結果の確認
print(f"\n✓ マッチ数: {len(matches)}")
for i, match in enumerate(matches):
    print(f"  マッチ {i}: 原子インデックス = {match}")
    print(f"    原子数 = {len(match)}")

# 検証（水素を除いた原子数で比較）
assert len(matches) > 0, "少なくとも1つのマッチが見つかる必要があります"
for match in matches:
    # マッチした原子数は水素を除いた部分構造の原子数と一致する必要がある
    assert len(match) == substructure_no_h.GetNumAtoms(), \
        f"マッチの原子数（{len(match)}）が水素なし部分構造の原子数（{substructure_no_h.GetNumAtoms()}）と一致しません"
    
print("\n✓ 全ての検証に合格")