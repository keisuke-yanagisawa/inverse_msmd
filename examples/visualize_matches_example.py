#!/usr/bin/env python3
"""A01を使った複数マッチ可視化のテスト"""

from rdkit import Chem
from inverse_msmd.substructure_replacement import (
    find_substructure_in_ligand,
    visualize_multiple_matches
)
from inverse_msmd.utils.mol_utils import read_mol_from_pdb_smi
from pathlib import Path

print("=" * 60)
print("A01を使った複数マッチ可視化テスト")
print("=" * 60)

# テストデータの読み込み
print("\n1. テストデータを読み込み中...")
ligand_mol = next(Chem.SDMolSupplier("data/atom_matching/4hw3_A_lig.sdf"))
a01_mol = read_mol_from_pdb_smi(
    "data/sample_probes/A01.pdb",
    "data/sample_probes/A01.smi"
)

print(f"   リガンド原子数: {ligand_mol.GetNumAtoms()}")
print(f"   A01原子数: {a01_mol.GetNumAtoms()}")

# 水素を除去した後の原子数も確認
ligand_no_h = Chem.RemoveHs(ligand_mol)
a01_no_h = Chem.RemoveHs(a01_mol)
print(f"   リガンド原子数（水素なし）: {ligand_no_h.GetNumAtoms()}")
print(f"   A01原子数（水素なし）: {a01_no_h.GetNumAtoms()}")

# 部分構造探索
print("\n2. 部分構造探索を実行中...")
matches = find_substructure_in_ligand(ligand_mol, a01_mol)
print(f"   ✓ マッチ数: {len(matches)}")

if len(matches) == 0:
    print("   ✗ マッチが見つかりませんでした")
    exit(1)

for i, match in enumerate(matches):
    print(f"   マッチ {i}: 原子インデックス = {match}")
    print(f"     原子数 = {len(match)}")

# PNG画像を生成
output_path = "test_output/substructure_matches_A01_test.png"
Path("test_output").mkdir(exist_ok=True)

print(f"\n3. PNG画像を生成中: {output_path}")
visualize_multiple_matches(
    ligand_mol, 
    a01_mol, 
    matches, 
    output_path
)

# 結果の確認
assert Path(output_path).exists(), "PNG画像が生成されていません"
assert Path(output_path).stat().st_size > 0, "PNG画像のサイズが0です"

print(f"   ✓ PNG画像を生成: {output_path}")
print(f"   ✓ ファイルサイズ: {Path(output_path).stat().st_size} bytes")

print("\n" + "=" * 60)
print("テスト完了！")
print("=" * 60)
print(f"\n📊 結果サマリー:")
print(f"   - マッチ数: {len(matches)}")
print(f"   - 出力画像: {output_path}")
print(f"\n👀 次のステップ:")
print(f"   画像を開いて以下を確認してください：")
print(f"   1. {len(matches)}個のマッチがそれぞれ表示されている")
print(f"   2. 各マッチでハイライト部分が異なる")
print(f"   3. 2D構造が見やすく表示されている")
print(f"\n   コマンド: code {output_path}")
print()