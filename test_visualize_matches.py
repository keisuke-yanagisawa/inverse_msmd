#!/usr/bin/env python3
"""複数マッチ可視化のテスト"""

from rdkit import Chem
from inverse_msmd.substructure_replacement import (
    find_substructure_in_ligand,
    visualize_multiple_matches
)
from inverse_msmd.utils.mol_utils import read_mol_from_pdb_smi
from pathlib import Path

# テストデータの読み込み
ligand_mol = next(Chem.SDMolSupplier("data/atom_matching/4hw3_A_lig.sdf"))
substructure_mol = read_mol_from_pdb_smi(
    "data/sample_probes/E23.pdb",
    "data/sample_probes/E23.smi"
)

# 部分構造探索
matches = find_substructure_in_ligand(ligand_mol, substructure_mol)
print(f"マッチ数: {len(matches)}")

# PNG画像を生成
output_path = "test_output/substructure_matches_test.png"
Path("test_output").mkdir(exist_ok=True)

print(f"\nPNG画像を生成中: {output_path}")
visualize_multiple_matches(
    ligand_mol, 
    substructure_mol, 
    matches, 
    output_path
)

# 結果の確認
assert Path(output_path).exists(), "PNG画像が生成されていません"
assert Path(output_path).stat().st_size > 0, "PNG画像のサイズが0です"

print(f"✓ PNG画像を生成: {output_path}")
print(f"  ファイルサイズ: {Path(output_path).stat().st_size} bytes")
print("\n手動確認: 画像を開いてマッチ部分がハイライトされているか確認してください")
print(f"  コマンド: code {output_path}")