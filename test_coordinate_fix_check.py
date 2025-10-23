#!/usr/bin/env python
"""
座標修正の検証スクリプト
E24の座標が保持され、リガンドとタンパク質が変換されていることを確認
"""

import numpy as np
from rdkit import Chem
from inverse_msmd.utils.bio_utils import PDB
from inverse_msmd.utils.mol_utils import read_mol_from_pdb_smi

print("=" * 70)
print("座標修正の検証")
print("=" * 70)

# 1. 元のE24の座標を取得
print("\n1. 元のE24の座標を読み込み...")
e24_original = read_mol_from_pdb_smi(
    "data/sample_probes/E24.pdb",
    "data/sample_probes/E24.smi"
)
e24_original_no_h = Chem.RemoveHs(e24_original)
e24_original_coords = e24_original_no_h.GetConformer().GetPositions()
print(f"  E24原子数: {len(e24_original_coords)}")
print(f"  E24の最初の原子座標: {e24_original_coords[0]}")

# 2. 出力されたリガンドのE24部分の座標を確認
print("\n2. 出力されたリガンドを読み込み...")
ligand_replaced = next(Chem.SDMolSupplier("test_output/fix_validation/pattern_0_ligand_replaced.sdf"))
ligand_replaced_coords = ligand_replaced.GetConformer().GetPositions()
print(f"  置換後リガンド原子数: {len(ligand_replaced_coords)}")

# リガンド中のE24部分（後半12原子がE24と仮定）
# 元のリガンドは25原子、E23は9原子、E24は12原子
# 置換後は25-9+12=28原子になるはず
expected_atoms = 25 - 9 + 12
print(f"  期待される原子数: {expected_atoms}")
print(f"  実際の原子数: {len(ligand_replaced_coords)}")

# E24部分の座標を抽出（最後の12原子）
e24_in_ligand_coords = ligand_replaced_coords[-12:]
print(f"  リガンド中のE24部分の最初の原子座標: {e24_in_ligand_coords[0]}")

# 3. E24の座標が保持されているか確認
print("\n3. E24の座標が保持されているか確認...")
# 座標の差を計算
coords_diff = np.abs(e24_original_coords - e24_in_ligand_coords)
max_diff = np.max(coords_diff)
mean_diff = np.mean(coords_diff)

print(f"  座標の最大差: {max_diff:.6f} Å")
print(f"  座標の平均差: {mean_diff:.6f} Å")

if max_diff < 0.01:  # 0.01 Å以下の差を許容
    print("  ✓ E24の座標が正しく保持されています")
else:
    print("  ✗ E24の座標に大きな差があります")
    print("\n  元のE24座標（最初の3原子）:")
    for i in range(min(3, len(e24_original_coords))):
        print(f"    原子{i}: {e24_original_coords[i]}")
    print("\n  リガンド中のE24座標（最初の3原子）:")
    for i in range(min(3, len(e24_in_ligand_coords))):
        print(f"    原子{i}: {e24_in_ligand_coords[i]}")

# 4. 元のタンパク質との位置関係を確認
print("\n4. タンパク質の変換を確認...")
protein_original = PDB.get_structure("data/sample_proteins/4hw3_A.pdb")
protein_original_coords = PDB.get_attr(protein_original, "coord")
print(f"  元のタンパク質の最初の原子座標: {protein_original_coords[0]}")

protein_transformed = PDB.get_structure("test_output/fix_validation/pattern_0_protein_aligned.pdb")
protein_transformed_coords = PDB.get_attr(protein_transformed, "coord")
print(f"  変換後タンパク質の最初の原子座標: {protein_transformed_coords[0]}")

# タンパク質が変換されているか確認
protein_diff = np.abs(protein_original_coords[0] - protein_transformed_coords[0])
protein_moved = np.linalg.norm(protein_diff)
print(f"  タンパク質の移動距離: {protein_moved:.6f} Å")

if protein_moved > 0.01:
    print("  ✓ タンパク質が変換されています")
else:
    print("  ✗ タンパク質が変換されていません")

# 5. 元のリガンドとの位置関係を確認
print("\n5. 元のリガンドとの位置関係を確認...")
ligand_original = next(Chem.SDMolSupplier("data/atom_matching/4hw3_A_lig.sdf"))
ligand_original_coords = ligand_original.GetConformer().GetPositions()
print(f"  元のリガンドの最初の原子座標: {ligand_original_coords[0]}")
print(f"  置換後リガンドの最初の原子座標: {ligand_replaced_coords[0]}")

ligand_diff = np.abs(ligand_original_coords[0] - ligand_replaced_coords[0])
ligand_moved = np.linalg.norm(ligand_diff)
print(f"  リガンドの移動距離: {ligand_moved:.6f} Å")

if ligand_moved > 0.01:
    print("  ✓ リガンドが変換されています")
else:
    print("  ✗ リガンドが変換されていません")

print("\n" + "=" * 70)
print("検証完了")
print("=" * 70)

# 最終判定
if max_diff < 0.01 and protein_moved > 0.01 and ligand_moved > 0.01:
    print("\n✓ 全ての検証に合格:")
    print("  - E24の座標が保持されている")
    print("  - タンパク質が変換されている")
    print("  - リガンドが変換されている")
    print("\n修正が正しく適用されています！")
else:
    print("\n✗ 検証に問題があります")
    if max_diff >= 0.01:
        print("  - E24の座標が保持されていない")
    if protein_moved <= 0.01:
        print("  - タンパク質が変換されていない")
    if ligand_moved <= 0.01:
        print("  - リガンドが変換されていない")