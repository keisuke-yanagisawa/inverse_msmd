# 統合部分構造置換機能 テストチェックリスト

本ドキュメントは、実装の各段階で行うべき動作チェック方法を詳細に記述します。
各機能を実装したら、必ずこのチェックリストに従って動作確認を行ってください。

## テスト環境の準備

### 必要なテストデータ
```bash
# 既存のサンプルデータを確認
ls -la data/atom_matching/4hw3_A_lig.sdf
ls -la data/sample_proteins/4hw3_A.pdb
ls -la data/sample_probes/E23.{pdb,smi}
ls -la data/sample_probes/E24.{pdb,smi}
```

---

## タスク1: モジュール基本構造の作成

### 実装内容
`inverse_msmd/substructure_replacement.py` の作成と基本インポート

### テストコード
```python
# test_imports.py
"""基本インポートのテスト"""

# モジュールがインポートできることを確認
try:
    from inverse_msmd import substructure_replacement
    print("✓ モジュールのインポート成功")
except ImportError as e:
    print(f"✗ インポートエラー: {e}")
    exit(1)

# 必要な依存関係の確認
try:
    from rdkit import Chem
    from Bio.PDB import Structure
    import numpy as np
    print("✓ 全ての依存関係がインポート可能")
except ImportError as e:
    print(f"✗ 依存関係エラー: {e}")
    exit(1)

print("\n全てのインポートチェック完了")
```

### 実行方法
```bash
cd /workspaces/inverse_msmd
python test_imports.py
```

### 期待される結果
```
✓ モジュールのインポート成功
✓ 全ての依存関係がインポート可能

全てのインポートチェック完了
```

### 確認ポイント
- [ ] モジュールがエラーなくインポートできる
- [ ] 必要な依存パッケージが全て利用可能

---

## タスク2: 部分構造探索関数の実装

### 実装内容
`find_substructure_in_ligand()` - リガンド中の部分構造を探索し、複数マッチを返す

### テストコード
```python
# test_find_substructure.py
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
print(f"部分構造原子数: {substructure_mol.GetNumAtoms()}")

# 部分構造探索の実行
print("\n部分構造探索を実行中...")
matches = find_substructure_in_ligand(ligand_mol, substructure_mol)

# 結果の確認
print(f"\n✓ マッチ数: {len(matches)}")
for i, match in enumerate(matches):
    print(f"  マッチ {i}: 原子インデックス = {match}")
    print(f"    原子数 = {len(match)}")

# 検証
assert len(matches) > 0, "少なくとも1つのマッチが見つかる必要があります"
for match in matches:
    assert len(match) == substructure_mol.GetNumAtoms(), "マッチの原子数が一致しません"
    
print("\n✓ 全ての検証に合格")
```

### 実行方法
```bash
python test_find_substructure.py
```

### 期待される結果
```
テストデータを読み込み中...
リガンド原子数: XX
部分構造原子数: YY

部分構造探索を実行中...

✓ マッチ数: N (1以上の整数)
  マッチ 0: 原子インデックス = (1, 3, 5, ...)
    原子数 = YY
  ...

✓ 全ての検証に合格
```

### 確認ポイント
- [ ] 少なくとも1つのマッチが見つかる
- [ ] マッチした原子数が部分構造の原子数と一致
- [ ] 複数マッチがある場合、全て正しく返される
- [ ] 各マッチの原子インデックスが有効範囲内

---

## タスク3: 複数マッチ可視化関数の実装

### 実装内容
`visualize_multiple_matches()` - 複数マッチをPNG画像として可視化

### テストコード
```python
# test_visualize_matches.py
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
```

### 実行方法
```bash
python test_visualize_matches.py
# 生成された画像を確認
# Linux: xdg-open test_output/substructure_matches_test.png
# または VSCode で画像を開く
```

### 期待される結果
```
マッチ数: N

PNG画像を生成中: test_output/substructure_matches_test.png
✓ PNG画像を生成: test_output/substructure_matches_test.png
  ファイルサイズ: XXXXX bytes

手動確認: 画像を開いてマッチ部分がハイライトされているか確認してください
```

### 確認ポイント
- [ ] PNG画像ファイルが正常に生成される
- [ ] ファイルサイズが0より大きい
- [ ] **手動確認**: 画像を開いて以下を確認
  - [ ] リガンド構造が表示されている
  - [ ] マッチした部分構造がハイライト（色付け）されている
  - [ ] 複数マッチがある場合、それぞれが識別可能
  - [ ] 画像が読みやすい解像度・サイズ

---

## タスク4: Atom Matching関数の実装

### 実装内容
`match_substructures()` - E23部分構造とE24のatom matchingを実行

### テストコード
```python
# test_atom_matching.py
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

# リガンド中のE23を探索
ligand_e23_matches = find_substructure_in_ligand(ligand_mol, e23_mol)
print(f"\nリガンド中のE23マッチ数: {len(ligand_e23_matches)}")

# 最初のマッチを使用してテスト
first_match = ligand_e23_matches[0]
print(f"使用するマッチ: {first_match}")

# リガンドのE23部分とE24のatom matchingを実行
# （リガンド部分構造の座標とE24分子を使用）
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
    assert np.all(pairs[0] < e23_mol.GetNumAtoms()), "E23のインデックスが範囲外"
    assert np.all(pairs[1] < e24_mol.GetNumAtoms()), "E24のインデックスが範囲外"

print("\n✓ 全ての検証に合格")
```

### 実行方法
```bash
python test_atom_matching.py
```

### 期待される結果
```
E23原子数: XX
E24原子数: YY

リガンド中のE23マッチ数: N
使用するマッチ: (1, 3, 5, ...)

Atom Matchingを実行中...

✓ Atom Matchingパターン数: M
パターン 0:
  形状: (2, Z)
  E23側インデックス: [0 1 2 ...]
  E24側インデックス: [0 1 2 ...]
  
✓ 全ての検証に合格
```

### 確認ポイント
- [ ] 少なくとも1つのマッチングパターンが見つかる
- [ ] atom pairsの形状が正しい（2行 x N列）
- [ ] インデックスが有効範囲内
- [ ] 対応する原子の元素が一致している（可能であれば）

---

## タスク5: Superimpose計算関数の実装

### 実装内容
`calculate_transformation()` - 変換行列（回転+並進）を計算

### テストコード
```python
# test_superimpose.py
"""Superimpose計算機能のテスト"""

from rdkit import Chem
import numpy as np
from inverse_msmd.substructure_replacement import (
    find_substructure_in_ligand,
    match_substructures,
    calculate_transformation
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

# リガンドとE23/E24の座標を取得
ligand_coords = ligand_mol.GetConformer().GetPositions()
e23_coords = e23_mol.GetConformer().GetPositions()
e24_coords = e24_mol.GetConformer().GetPositions()

print(f"リガンド座標形状: {ligand_coords.shape}")
print(f"E23座標形状: {e23_coords.shape}")
print(f"E24座標形状: {e24_coords.shape}")

# マッチング情報を取得
ligand_e23_matches = find_substructure_in_ligand(ligand_mol, e23_mol)
atom_pair_patterns = match_substructures(e23_mol, e24_mol)

# 最初のパターンでテスト
match = ligand_e23_matches[0]
atom_pairs = atom_pair_patterns[0]

print(f"\n使用するマッチ: {match}")
print(f"Atom pairs形状: {atom_pairs.shape}")

# リガンドのE23部分の座標を抽出
ligand_e23_coords = ligand_coords[list(match)]

# 変換行列を計算
print("\nSuperimpose計算を実行中...")
rot, tran = calculate_transformation(
    ligand_e23_coords,
    e24_coords,
    atom_pairs
)

print(f"\n✓ 回転行列 (rot):")
print(rot)
print(f"  形状: {rot.shape}")

print(f"\n✓ 並進ベクトル (tran):")
print(tran)
print(f"  形状: {tran.shape}")

# 検証
assert rot.shape == (3, 3), "回転行列は3x3である必要があります"
assert tran.shape == (3,), "並進ベクトルは長さ3である必要があります"

# 回転行列の直交性をチェック
rot_test = np.dot(rot, rot.T)
identity = np.eye(3)
assert np.allclose(rot_test, identity, atol=1e-6), "回転行列が直交行列ではありません"

# 行列式が1に近いことをチェック（回転行列の性質）
det = np.linalg.det(rot)
assert np.abs(det - 1.0) < 1e-6, f"回転行列の行列式が1ではありません: {det}"

print("\n✓ 全ての検証に合格")
print("  - 回転行列は正しい形状")
print("  - 回転行列は直交行列")
print("  - 回転行列の行列式は1")
```

### 実行方法
```bash
python test_superimpose.py
```

### 期待される結果
```
リガンド座標形状: (XX, 3)
E23座標形状: (YY, 3)
E24座標形状: (ZZ, 3)

使用するマッチ: (1, 3, 5, ...)
Atom pairs形状: (2, N)

Superimpose計算を実行中...

✓ 回転行列 (rot):
[[ 0.xxx  0.xxx  0.xxx]
 [ 0.xxx  0.xxx  0.xxx]
 [ 0.xxx  0.xxx  0.xxx]]
  形状: (3, 3)

✓ 並進ベクトル (tran):
[ X.xxx  Y.xxx  Z.xxx]
  形状: (3,)

✓ 全ての検証に合格
  - 回転行列は正しい形状
  - 回転行列は直交行列
  - 回転行列の行列式は1
```

### 確認ポイント
- [ ] 回転行列の形状が3x3
- [ ] 並進ベクトルの形状が(3,)
- [ ] 回転行列が直交行列（R^T * R = I）
- [ ] 回転行列の行列式が1に近い

---

## タスク6: タンパク質変換関数の実装

### 実装内容
`apply_transformation_to_protein()` - 変換行列をタンパク質に適用

### テストコード
```python
# test_protein_transformation.py
"""タンパク質変換機能のテスト"""

from rdkit import Chem
import numpy as np
from inverse_msmd.substructure_replacement import (
    find_substructure_in_ligand,
    match_substructures,
    calculate_transformation,
    apply_transformation_to_protein
)
from inverse_msmd.utils.mol_utils import read_mol_from_pdb_smi
from inverse_msmd.utils.bio_utils import PDB

# テストデータの読み込み
print("テストデータを読み込み中...")
ligand_mol = next(Chem.SDMolSupplier("data/atom_matching/4hw3_A_lig.sdf"))
protein = PDB.get_structure("data/sample_proteins/4hw3_A.pdb")
e23_mol = read_mol_from_pdb_smi(
    "data/sample_probes/E23.pdb",
    "data/sample_probes/E23.smi"
)
e24_mol = read_mol_from_pdb_smi(
    "data/sample_probes/E24.pdb",
    "data/sample_probes/E24.smi"
)

# 変換行列を計算
ligand_coords = ligand_mol.GetConformer().GetPositions()
e23_coords = e23_mol.GetConformer().GetPositions()
e24_coords = e24_mol.GetConformer().GetPositions()

ligand_e23_matches = find_substructure_in_ligand(ligand_mol, e23_mol)
atom_pair_patterns = match_substructures(e23_mol, e24_mol)

match = ligand_e23_matches[0]
atom_pairs = atom_pair_patterns[0]
ligand_e23_coords = ligand_coords[list(match)]

rot, tran = calculate_transformation(ligand_e23_coords, e24_coords, atom_pairs)

# 元のタンパク質座標を取得
original_coords = PDB.get_attr(protein, "coord")
print(f"タンパク質原子数: {len(original_coords)}")
print(f"元の座標の最初の3原子:\n{original_coords[:3]}")

# 変換を適用
print("\nタンパク質に変換を適用中...")
transformed_protein = apply_transformation_to_protein(protein, rot, tran)

# 変換後の座標を取得
new_coords = PDB.get_attr(transformed_protein, "coord")
print(f"\n変換後の座標の最初の3原子:\n{new_coords[:3]}")

# 検証
assert new_coords.shape == original_coords.shape, "座標の形状が変わっています"
assert not np.allclose(new_coords, original_coords), "座標が変換されていません"

# 変換が正しく適用されているか確認（手動計算と比較）
expected_coords = np.dot(original_coords, rot) + tran
assert np.allclose(new_coords, expected_coords, atol=1e-6), "変換が正しく適用されていません"

print("\n✓ 全ての検証に合格")
print("  - 座標の形状は保持されている")
print("  - 座標が実際に変換されている")
print("  - 変換式（rot @ coords + tran）が正しく適用されている")

# 変換後のタンパク質を保存してテスト
test_output = "test_output/protein_transformed_test.pdb"
PDB.save(transformed_protein, test_output)
print(f"\n✓ 変換後のタンパク質を保存: {test_output}")
print("  手動確認: PyMOLやChimeraで元のタンパク質と比較してください")
```

### 実行方法
```bash
python test_protein_transformation.py
```

### 期待される結果
```
テストデータを読み込み中...
タンパク質原子数: XXXX
元の座標の最初の3原子:
[[x1 y1 z1]
 [x2 y2 z2]
 [x3 y3 z3]]

タンパク質に変換を適用中...

変換後の座標の最初の3原子:
[[x1' y1' z1']
 [x2' y2' z2']
 [x3' y3' z3']]

✓ 全ての検証に合格
  - 座標の形状は保持されている
  - 座標が実際に変換されている
  - 変換式（rot @ coords + tran）が正しく適用されている

✓ 変換後のタンパク質を保存: test_output/protein_transformed_test.pdb
  手動確認: PyMOLやChimeraで元のタンパク質と比較してください
```

### 確認ポイント
- [ ] 座標の形状が保持されている
- [ ] 座標が実際に変換されている
- [ ] 変換式が正しく適用されている
- [ ] **手動確認**: PyMOLで元のPDBと変換後のPDBを重ねて表示
  - [ ] タンパク質が回転・並進している
  - [ ] 構造が破壊されていない

---

## タスク7: リガンド置換関数の実装

### 実装内容
`replace_ligand_substructure()` - リガンドの部分構造をE24で置換

### テストコード
```python
# test_ligand_replacement.py
"""リガンド部分構造置換のテスト"""

from rdkit import Chem
from inverse_msmd.substructure_replacement import (
    find_substructure_in_ligand,
    match_substructures,
    replace_ligand_substructure
)
from inverse_msmd.utils.mol_utils import read_mol_from_pdb_smi

# テストデータの読み込み
print("テストデータを読み込み中...")
ligand_mol = next(Chem.SDMolSupplier("data/atom_matching/4hw3_A_lig.sdf"))
e23_mol = read_mol_from_pdb_smi(
    "data/sample_probes/E23.pdb",
    "data/sample_probes/E23.smi"
)
e24_mol = read_mol_from_pdb_smi(
    "data/sample_probes/E24.pdb",
    "data/sample_probes/E24.smi"
)

print(f"元のリガンド原子数: {ligand_mol.GetNumAtoms()}")
print(f"E23原子数: {e23_mol.GetNumAtoms()}")
print(f"E24原子数: {e24_mol.GetNumAtoms()}")

# マッチング情報を取得
ligand_e23_matches = find_substructure_in_ligand(ligand_mol, e23_mol)
atom_pair_patterns = match_substructures(e23_mol, e24_mol)

match = ligand_e23_matches[0]
atom_pairs = atom_pair_patterns[0]

print(f"\n使用するマッチ: {match}")
print(f"Atom pairs: {atom_pairs}")

# 部分構造を置換
print("\nリガンドの部分構造を置換中...")
replaced_ligand = replace_ligand_substructure(
    ligand_mol,
    match,
    e24_mol,
    atom_pairs
)

print(f"\n✓ 置換後のリガンド原子数: {replaced_ligand.GetNumAtoms()}")

# 期待される原子数を計算
expected_atoms = (ligand_mol.GetNumAtoms() - 
                  e23_mol.GetNumAtoms() + 
                  e24_mol.GetNumAtoms())
print(f"  期待される原子数: {expected_atoms}")

# 検証
assert replaced_ligand is not None, "置換に失敗しました"

# Sanitizeチェック
try:
    Chem.SanitizeMol(replaced_ligand)
    print("✓ 置換後の分子は化学的に妥当")
except:
    print("✗ 警告: 置換後の分子のSanitizeに失敗")

# SDFファイルとして保存
test_output = "test_output/ligand_replaced_test.sdf"
writer = Chem.SDWriter(test_output)
writer.SetKekulize(False)
writer.write(replaced_ligand)
writer.close()

print(f"\n✓ 置換後のリガンドを保存: {test_output}")
print("  手動確認: ChemDrawやRDKitで元のリガンドと比較してください")

# SMILES比較（参考）
original_smiles = Chem.MolToSmiles(Chem.RemoveHs(ligand_mol))
replaced_smiles = Chem.MolToSmiles(Chem.RemoveHs(replaced_ligand))
print(f"\n元のSMILES: {original_smiles}")
print(f"置換後のSMILES: {replaced_smiles}")
```

### 実行方法
```bash
python test_ligand_replacement.py
```

### 期待される結果
```
テストデータを読み込み中...
元のリガンド原子数: XX
E23原子数: YY
E24原子数: ZZ

使用するマッチ: (1, 3, 5, ...)
Atom pairs: [[...] [...]]

リガンドの部分構造を置換中...

✓ 置換後のリガンド原子数: WW
  期待される原子数: WW
✓ 置換後の分子は化学的に妥当

✓ 置換後のリガンドを保存: test_output/ligand_replaced_test.sdf
  手動確認: ChemDrawやRDKitで元のリガンドと比較してください

元のSMILES: ...
置換後のSMILES: ...
```

### 確認ポイント
- [ ] 置換後の原子数が期待値と一致
- [ ] 分子がSanitizeを通過（化学的に妥当）
- [ ] SDFファイルが正常に保存される
- [ ] **手動確認**: 分子構造ビューアで確認
  - [ ] E23部分がE24に置換されている
  - [ ] 結合が適切に保持されている
  - [ ] 不自然な結合切断がない

---

## タスク8-11: 統合ワークフロー・CLI・出力機能

### 統合テスト
全ての機能を組み合わせた統合テストを実行

### テストコード
```python
# test_integrated_workflow.py
"""統合ワークフロー全体

## タスク8-11: 統合ワークフロー・CLI・出力機能

### 統合テスト
全ての機能を組み合わせた統合テストを実行

### テストコード
```python
# test_integrated_workflow.py
"""統合ワークフロー全体のテスト"""

from pathlib import Path
from inverse_msmd.substructure_replacement import integrated_substructure_replacement

# 出力ディレクトリを作成
output_dir = "test_output/integrated"
Path(output_dir).mkdir(parents=True, exist_ok=True)

print("統合ワークフローを実行中...")
print("=" * 60)

# 統合関数を実行
results = integrated_substructure_replacement(
    ligand_file="data/atom_matching/4hw3_A_lig.sdf",
    protein_file="data/sample_proteins/4hw3_A.pdb",
    from_file="data/sample_probes/E23",
    to_file="data/sample_probes/E24",
    output_dir=output_dir,
    match_index=None  # 自動選択またはユーザー選択
)

print("\n" + "=" * 60)
print(f"✓ 処理完了: {len(results)} パターンの結果を生成")

# 結果の検証
for i, result in enumerate(results):
    print(f"\nパターン {i}:")
    print(f"  リガンドファイル: {result['ligand_file']}")
    print(f"  タンパク質ファイル: {result['protein_file']}")
    
    # ファイルの存在確認
    assert Path(result['ligand_file']).exists(), f"リガンドファイルが存在しません"
    assert Path(result['protein_file']).exists(), f"タンパク質ファイルが存在しません"
    
    # ファイルサイズ確認
    ligand_size = Path(result['ligand_file']).stat().st_size
    protein_size = Path(result['protein_file']).stat().st_size
    
    assert ligand_size > 0, "リガンドファイルが空です"
    assert protein_size > 0, "タンパク質ファイルが空です"
    
    print(f"  ✓ ファイル存在確認OK")
    print(f"  ✓ リガンドサイズ: {ligand_size} bytes")
    print(f"  ✓ タンパク質サイズ: {protein_size} bytes")

print("\n✓ 全ての検証に合格")
```

### 実行方法
```bash
python test_integrated_workflow.py
```

### 確認ポイント
- [ ] 全てのパターンで出力ファイルが生成される
- [ ] ファイルサイズが0より大きい
- [ ] **手動確認**: 各パターンの結果を可視化して正しく変換されているか確認

---

## タスク10: CLIスクリプトの実装

### CLIテストコマンド

#### ケース1: 基本的な使用
```bash
python scripts/integrated_replacement.py \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --from-file data/sample_probes/E23 \
    --to-file data/sample_probes/E24 \
    --output test_output/cli_test1/
```

#### ケース2: ヘルプ表示
```bash
python scripts/integrated_replacement.py --help
```

### 確認ポイント
- [ ] 全てのコマンドラインオプションが機能する
- [ ] ヘルプメッセージが適切に表示される
- [ ] エラー時に分かりやすいメッセージが出る

---

## 最終統合テスト

### 実際のデータでのエンドツーエンドテスト

```bash
python scripts/integrated_replacement.py \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --from-file data/sample_probes/E23 \
    --to-file data/sample_probes/E24 \
    --output output/final_test/ \
    --verbose
```

### 最終確認チェックリスト

#### 1. 機能確認
- [ ] 部分構造探索が正しく動作する
- [ ] Atom matchingが複数パターンを正しく見つける
- [ ] Superimposeの計算が正確
- [ ] タンパク質の座標変換が正しい
- [ ] リガンドの部分構造置換が正しい
- [ ] 全パターンで出力ファイルが生成される

#### 2. 出力ファイル確認
- [ ] `pattern_N_ligand_replaced.sdf` が存在する
- [ ] `pattern_N_protein_aligned.pdb` が存在する
- [ ] 全てのファイルサイズが0より大きい

#### 3. 構造の視覚的確認
- [ ] PyMOLで元のPDBと変換後のPDBを比較
- [ ] 分子ビューアで元のリガンドと置換後のリガンドを比較
- [ ] E23部分がE24に置き換わっている
- [ ] 結合が適切に保持されている

---

## テスト実行の推奨順序

1. **タスク1**: モジュール基本構造
2. **タスク2**: 部分構造探索
3. **タスク3**: 可視化
4. **タスク4**: Atom Matching
5. **タスク5**: Superimpose
6. **タスク6**: タンパク質変換
7. **タスク7**: リガンド置換
8. **タスク8-11**: 統合ワークフロー
9. **タスク10**: CLI
10. **タスク13**: エラーハンドリング
11. **最終**: 統合テスト

各テストをパスしてから次のステップに進むことで、問題の早期発見と修正が可能になります。