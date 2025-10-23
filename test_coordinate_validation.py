#!/usr/bin/env python
"""座標保持の検証テスト"""

from rdkit import Chem
import numpy as np
from inverse_msmd.utils.mol_utils import read_mol_from_pdb_smi

print("=" * 70)
print("座標保持の検証テスト")
print("=" * 70)

# 元のリガンドを読み込み
print("\n1. 元のリガンドを読み込み...")
original_ligand = next(Chem.SDMolSupplier("data/atom_matching/4hw3_A_lig.sdf"))
original_coords = original_ligand.GetConformer().GetPositions()
print(f"   原子数: {original_ligand.GetNumAtoms()}")
print(f"   座標形状: {original_coords.shape}")

# 置換後のリガンドを読み込み（パターン0を使用）
print("\n2. 置換後のリガンドを読み込み...")
replaced_ligand = next(Chem.SDMolSupplier("test_output/integrated_simple/pattern_0_ligand_replaced.sdf"))
replaced_coords = replaced_ligand.GetConformer().GetPositions()
print(f"   原子数: {replaced_ligand.GetNumAtoms()}")
print(f"   座標形状: {replaced_coords.shape}")

# E23とE24を読み込み
print("\n3. 部分構造を読み込み...")
e23_mol = read_mol_from_pdb_smi("data/sample_probes/E23.pdb", "data/sample_probes/E23.smi")
e24_mol = read_mol_from_pdb_smi("data/sample_probes/E24.pdb", "data/sample_probes/E24.smi")

# 水素を除去
original_ligand_no_h = Chem.RemoveHs(original_ligand)
replaced_ligand_no_h = Chem.RemoveHs(replaced_ligand)
e23_no_h = Chem.RemoveHs(e23_mol)

print(f"   E23原子数（水素なし）: {e23_no_h.GetNumAtoms()}")
print(f"   E24原子数（水素なし）: {Chem.RemoveHs(e24_mol).GetNumAtoms()}")

# マッチ部分を特定
print("\n4. マッチ部分を特定...")
matches = original_ligand_no_h.GetSubstructMatches(e23_no_h)
print(f"   マッチ数: {len(matches)}")
if len(matches) == 0:
    print("   エラー: マッチが見つかりません")
    exit(1)

match = matches[0]
print(f"   使用するマッチ: {match}")

# 置換されていない原子のインデックスを特定
match_set = set(match)
non_replaced_indices = [i for i in range(original_ligand_no_h.GetNumAtoms()) if i not in match_set]
print(f"\n5. 置換されていない原子: {len(non_replaced_indices)}個")
print(f"   インデックス: {non_replaced_indices}")

# 置換されていない原子の座標を比較
# 注意: 元のリガンドと置換後のリガンドでは変換が適用されているため、
# 直接の座標比較ではなく、変換後の座標を計算する必要があります

print("\n6. 座標の比較...")
print("   注意: 置換されていない原子は座標変換が適用されているため、")
print("   元の座標とは異なることが期待されます。")
print("   しかし、すべての非置換原子は同じ変換を受けているはずです。")

# 元のリガンドの座標（水素なし）
original_coords_no_h = original_ligand_no_h.GetConformer().GetPositions()

# 最初の非置換原子で確認
if len(non_replaced_indices) > 0:
    idx = non_replaced_indices[0]
    print(f"\n   非置換原子 {idx} の座標:")
    print(f"     元の座標: {original_coords_no_h[idx]}")
    
    # 置換後のリガンドで対応する原子を見つける
    # 原子の削除と追加があるため、インデックスが変わっている可能性がある
    # SMARTSパターンでマッチングして対応を見つける
    print("\n7. 置換後のリガンドでの原子の対応を確認...")
    
    # 単純な比較: 原子数の違いを確認
    atom_count_diff = replaced_ligand_no_h.GetNumAtoms() - original_ligand_no_h.GetNumAtoms()
    print(f"   原子数の差: {atom_count_diff}")
    print(f"   期待される差: E24原子数 - E23原子数 = {Chem.RemoveHs(e24_mol).GetNumAtoms() - e23_no_h.GetNumAtoms()}")

print("\n" + "=" * 70)
print("検証完了")
print("=" * 70)
print("\n重要: このテストでは、置換されていない原子が同じ座標変換を")
print("受けていることを確認する必要があります。元の座標と直接比較")
print("するのではなく、相対的な位置関係が保持されているかを確認")
print("する必要があります。")