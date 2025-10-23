#!/usr/bin/env python
"""座標保持の詳細検証"""

from rdkit import Chem
import numpy as np

print("=" * 70)
print("座標保持の詳細検証")
print("=" * 70)

# 元のリガンドを読み込み
print("\n1. 元のリガンドを読み込み...")
original_ligand = next(Chem.SDMolSupplier("data/atom_matching/4hw3_A_lig.sdf"))
original_ligand_no_h = Chem.RemoveHs(original_ligand)
original_coords = original_ligand_no_h.GetConformer().GetPositions()
print(f"   原子数（水素なし）: {original_ligand_no_h.GetNumAtoms()}")

# 置換後のリガンドを読み込み
print("\n2. 置換後のリガンドを読み込み...")
replaced_ligand = next(Chem.SDMolSupplier("test_output/integrated_simple/pattern_0_ligand_replaced.sdf"))
replaced_coords = replaced_ligand.GetConformer().GetPositions()
print(f"   原子数: {replaced_ligand.GetNumAtoms()}")

# E23のマッチ部分を特定
print("\n3. E23のマッチ部分を特定...")
from inverse_msmd.utils.mol_utils import read_mol_from_pdb_smi
e23_mol = read_mol_from_pdb_smi("data/sample_probes/E23.pdb", "data/sample_probes/E23.smi")
e23_no_h = Chem.RemoveHs(e23_mol)

matches = original_ligand_no_h.GetSubstructMatches(e23_no_h)
if len(matches) == 0:
    print("   エラー: マッチが見つかりません")
    exit(1)

match = matches[0]
print(f"   マッチインデックス: {match}")
print(f"   マッチ原子数: {len(match)}")

# 置換されていない原子を特定
match_set = set(match)
non_replaced_indices = [i for i in range(original_ligand_no_h.GetNumAtoms()) if i not in match_set]
print(f"\n4. 置換されていない原子: {len(non_replaced_indices)}個")
print(f"   インデックス: {non_replaced_indices}")

# 座標の比較
print("\n5. 座標の比較...")
print("   注意: 元のリガンドと置換後のリガンドで原子数が異なるため、")
print("   インデックスの対応を確認する必要があります。")

# 元のリガンドの各原子について、置換後のリガンドでの対応を見つける
# 簡易的な方法: 最初のN個の非置換原子について座標を比較
print("\n   非置換原子の座標比較（最初の5個）:")
for i, idx in enumerate(non_replaced_indices[:5]):
    orig_coord = original_coords[idx]
    
    # 置換後のリガンドでは、マッチ部分が削除され、
    # 新しい部分が追加されているため、インデックスが変わっている
    # 単純な仮定: 非置換原子は順序が保持されている場合
    # しかし、実際には削除・追加の影響でインデックスがずれる
    
    atom = original_ligand_no_h.GetAtomWithIdx(idx)
    print(f"\n   元のリガンド 原子{idx} ({atom.GetSymbol()}):")
    print(f"     座標: [{orig_coord[0]:.4f}, {orig_coord[1]:.4f}, {orig_coord[2]:.4f}]")

# より正確な検証: SMARTSパターンで対応を見つける
print("\n6. より正確な検証...")
print("   元のリガンドと置換後のリガンドを重ね合わせて、")
print("   非置換部分の座標が保持されているかを確認します。")

# E24の原子数を確認
e24_mol = read_mol_from_pdb_smi("data/sample_probes/E24.pdb", "data/sample_probes/E24.smi")
e24_no_h = Chem.RemoveHs(e24_mol)
atom_count_diff = e24_no_h.GetNumAtoms() - e23_no_h.GetNumAtoms()

print(f"\n   元のリガンド原子数: {original_ligand_no_h.GetNumAtoms()}")
print(f"   E23原子数: {e23_no_h.GetNumAtoms()}")
print(f"   E24原子数: {e24_no_h.GetNumAtoms()}")
print(f"   原子数の差: {atom_count_diff}")
print(f"   期待される置換後の原子数: {original_ligand_no_h.GetNumAtoms() - e23_no_h.GetNumAtoms() + e24_no_h.GetNumAtoms()}")
print(f"   実際の置換後の原子数: {replaced_ligand.GetNumAtoms()}")

# 非置換部分の座標が保持されているかを確認
# インデックス変換マップを作成
# match内の原子は削除される
# マッチより前の原子はインデックスが保持される
# マッチより後の原子は削除された分だけ前にずれる

print("\n7. インデックス変換マップの作成...")
# ソートされたマッチインデックス
sorted_match = sorted(match)

# 元のインデックス -> 置換後のインデックス
index_map = {}
deleted_count = 0
for i in range(original_ligand_no_h.GetNumAtoms()):
    if i in match_set:
        deleted_count += 1
        # この原子は削除される
        continue
    else:
        # この原子は保持される
        # 削除された原子の分だけ前にずれる
        index_map[i] = i - deleted_count

print(f"   インデックスマップ（最初の10個）:")
for orig_idx in sorted(index_map.keys())[:10]:
    new_idx = index_map[orig_idx]
    atom = original_ligand_no_h.GetAtomWithIdx(orig_idx)
    print(f"     {orig_idx} -> {new_idx} ({atom.GetSymbol()})")

# 座標の一致度を確認
print("\n8. 座標の一致度を確認...")
max_diff = 0.0
total_diff = 0.0
count = 0

for orig_idx, new_idx in index_map.items():
    orig_coord = original_coords[orig_idx]
    new_coord = replaced_coords[new_idx]
    
    diff = np.linalg.norm(orig_coord - new_coord)
    total_diff += diff
    max_diff = max(max_diff, diff)
    count += 1

avg_diff = total_diff / count if count > 0 else 0.0

print(f"\n   非置換原子数: {count}")
print(f"   平均座標差: {avg_diff:.6f} Å")
print(f"   最大座標差: {max_diff:.6f} Å")

if max_diff < 1e-4:
    print("\n   ✓ 座標が完全に保持されています！")
elif max_diff < 1e-2:
    print("\n   ⚠ 座標がほぼ保持されています（許容範囲内）")
else:
    print("\n   ✗ 座標が保持されていません！")
    print("\n   詳細確認（差が大きい原子）:")
    for orig_idx, new_idx in sorted(index_map.items(), key=lambda x: x[0])[:5]:
        orig_coord = original_coords[orig_idx]
        new_coord = replaced_coords[new_idx]
        diff = np.linalg.norm(orig_coord - new_coord)
        atom = original_ligand_no_h.GetAtomWithIdx(orig_idx)
        print(f"     原子{orig_idx} ({atom.GetSymbol()}): 差={diff:.6f} Å")
        print(f"       元: [{orig_coord[0]:.4f}, {orig_coord[1]:.4f}, {orig_coord[2]:.4f}]")
        print(f"       新: [{new_coord[0]:.4f}, {new_coord[1]:.4f}, {new_coord[2]:.4f}]")

print("\n" + "=" * 70)
print("検証完了")
print("=" * 70)