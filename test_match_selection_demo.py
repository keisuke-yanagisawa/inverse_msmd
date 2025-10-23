#!/usr/bin/env python3
"""マッチ選択のデモンストレーション"""

from rdkit import Chem
from inverse_msmd.substructure_replacement import (
    find_substructure_in_ligand,
    visualize_multiple_matches
)
from inverse_msmd.utils.mol_utils import read_mol_from_pdb_smi
from pathlib import Path

print("=" * 70)
print("マッチ選択デモンストレーション")
print("=" * 70)

# テストデータの読み込み
print("\n📂 ステップ1: データ読み込み")
ligand_mol = next(Chem.SDMolSupplier("data/atom_matching/4hw3_A_lig.sdf"))
a01_mol = read_mol_from_pdb_smi(
    "data/sample_probes/A01.pdb",
    "data/sample_probes/A01.smi"
)

# 部分構造探索
print("\n🔍 ステップ2: 部分構造探索")
matches = find_substructure_in_ligand(ligand_mol, a01_mol)
print(f"   検出されたマッチ数: {len(matches)}")

# 可視化画像を生成
output_path = "test_output/match_selection_guide.png"
Path("test_output").mkdir(exist_ok=True)
visualize_multiple_matches(ligand_mol, a01_mol, matches, output_path)
print(f"   可視化画像: {output_path}")

print("\n" + "=" * 70)
print("✅ 複数マッチが見つかりました！")
print("=" * 70)

print(f"""
🎯 マッチを選択する方法

【方法1】Pythonコードから選択
────────────────────────────────────────
以下のようにmatch_indexパラメータを指定します：

```python
from inverse_msmd.substructure_replacement import integrated_substructure_replacement

# Match 0 を選択する場合
results = integrated_substructure_replacement(
    ligand_file="data/atom_matching/4hw3_A_lig.sdf",
    protein_file="data/sample_proteins/4hw3_A.pdb",
    from_file="data/sample_probes/A01",
    to_file="data/sample_probes/E24",
    output_dir="output/integrated/",
    match_index=0  # ← ここで選択（0始まり）
)

# Match 1 を選択する場合
results = integrated_substructure_replacement(
    ligand_file="data/atom_matching/4hw3_A_lig.sdf",
    protein_file="data/sample_proteins/4hw3_A.pdb",
    from_file="data/sample_probes/A01",
    to_file="data/sample_probes/E24",
    output_dir="output/integrated/",
    match_index=1  # ← Match 1を選択
)
```

【方法2】CLIから選択（将来実装予定）
────────────────────────────────────────
コマンドラインから --match-index オプションで指定：

# Match 0 を選択
$ python scripts/integrated_replacement.py \\
    --ligand data/atom_matching/4hw3_A_lig.sdf \\
    --protein data/sample_proteins/4hw3_A.pdb \\
    --from-file data/sample_probes/A01 \\
    --to-file data/sample_probes/E24 \\
    --output output/integrated/ \\
    --match-index 0

# Match 1 を選択
$ python scripts/integrated_replacement.py \\
    --ligand data/atom_matching/4hw3_A_lig.sdf \\
    --protein data/sample_proteins/4hw3_A.pdb \\
    --from-file data/sample_probes/A01 \\
    --to-file data/sample_probes/E24 \\
    --output output/integrated/ \\
    --match-index 1

【方法3】match_indexを指定しない場合（デフォルト）
────────────────────────────────────────
match_index=None の場合：
- 1つのマッチのみ → 自動的にそれを使用
- 複数マッチ → PNG画像を出力し、ユーザーに選択を促す

""")

print("=" * 70)
print("📊 現在検出されているマッチ:")
print("=" * 70)

for i, match in enumerate(matches):
    print(f"\n[Match {i}]")
    print(f"  原子インデックス: {match}")
    print(f"  原子数: {len(match)}")
    print(f"  選択方法: match_index={i} または --match-index {i}")

print("\n" + "=" * 70)
print("💡 ヒント:")
print("=" * 70)
print(f"""
1. まず可視化画像を確認: {output_path}
2. 使いたいマッチを決定（Match 0, Match 1, ...）
3. match_indexパラメータでそのインデックスを指定

例：Match 1を使いたい場合
  → match_index=1 と指定
""")