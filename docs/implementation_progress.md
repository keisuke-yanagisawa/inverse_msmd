# 統合部分構造置換機能 実装進捗記録

このドキュメントは、実装の進捗状況と各タスクの詳細な記録を提供します。

**最終更新**: 2025-10-23 (MCS重ね合わせ問題修正完了)

---

## 📊 進捗サマリー

| フェーズ | 完了 | 進行中 | 未着手 | 進捗率 |
|---------|------|--------|--------|--------|
| Phase 1: 基本機能実装 | 4 | 0 | 0 | 100% ✅ |
| Phase 2: 座標変換機能 | 3 | 0 | 0 | 100% ✅ |
| Phase 3: 統合とインターフェース | 3 | 0 | 0 | 100% ✅ |
| Phase 4: 品質保証 | 3 | 0 | 0 | 100% ✅ |
| **合計** | **13** | **0** | **0** | **100%** ✅ |

---

## ✅ Phase 1: 基本機能実装 (完了)

### T1: モジュール基本構造作成

**実装日**: 2025-10-23
**ファイル**: [`inverse_msmd/substructure_replacement.py`](../inverse_msmd/substructure_replacement.py)
**テストファイル**: [`tests/unit/test_imports.py`](../tests/unit/test_imports.py)

#### 実装内容
- 全7関数のスタブ定義完了
- 完全な型ヒント（`typing.List`, `numpy.typing.npt`）
- 詳細なdocstring（Google形式）
- 必要なインポート文の整理

#### 実装された関数スタブ
1. `find_substructure_in_ligand()`
2. `visualize_multiple_matches()`
3. `match_substructures()`
4. `calculate_transformation()`
5. `apply_transformation_to_protein()`
6. `replace_ligand_substructure()`
7. `integrated_substructure_replacement()`

#### テスト結果
```
✓ モジュールのインポート成功
✓ 全ての依存関係がインポート可能
✓ 全ての関数スタブが定義されている
```

---

### T2: 部分構造探索関数実装

**実装日**: 2025-10-23
**関数**: [`find_substructure_in_ligand()`](../inverse_msmd/substructure_replacement.py:46)
**テストファイル**: [`tests/unit/test_substructure_search.py`](../tests/unit/test_substructure_search.py)

#### 実装内容
```python
def find_substructure_in_ligand(
    ligand_mol: Chem.Mol,
    substructure_mol: Chem.Mol
) -> List[Tuple[int, ...]]:
    # 水素を除いた分子で処理
    ligand_no_h = Chem.RemoveHs(ligand_mol)
    substructure_no_h = Chem.RemoveHs(substructure_mol)
    
    # 部分構造を検索（全てのマッチを返す）
    matches = ligand_no_h.GetSubstructMatches(substructure_no_h)
    
    return list(matches)
```

#### 技術的詳細
- RDKitの`GetSubstructMatches()`を使用
- 水素原子を除去してマッチング精度を向上
- 複数マッチを全てリストとして返却

#### テスト結果
```
リガンド原子数: 25
部分構造原子数（水素あり）: 18
リガンド原子数（水素なし）: 25
部分構造原子数（水素なし）: 9

✓ マッチ数: 1
  マッチ 0: 原子インデックス = (0, 12, 6, 14, 7, 13, 1, 15, 23)
    原子数 = 9

✓ 全ての検証に合格
```

#### 参考コード
- [`scripts/replace_substructure.py:186`](../scripts/replace_substructure.py) - 部分構造探索の実装例

---

### T3: 複数マッチ可視化関数実装

**実装日**: 2025-10-23
**関数**: [`visualize_multiple_matches()`](../inverse_msmd/substructure_replacement.py:81)
**テストファイル**: [`tests/unit/test_visualization.py`](../tests/unit/test_visualization.py)

#### 実装内容
- matplotlibとRDKit Drawを使用したPNG画像生成
- マッチ部分の自動ハイライト
- グリッドレイアウト（最大4列）
- 出力ディレクトリの自動作成

#### 技術的詳細
```python
# マッチした原子をハイライト
img = Draw.MolToImage(
    ligand_mol,
    size=(400, 400),
    highlightAtoms=list(match)
)
```

#### テスト結果
```
マッチ数: 1

PNG画像を生成中: test_output/substructure_matches_test.png
✓ PNG画像を生成: test_output/substructure_matches_test.png
  ファイルサイズ: 47974 bytes
```

#### 出力例
- ファイル: [`test_output/substructure_matches_test.png`](../test_output/substructure_matches_test.png)
- サイズ: 47,974 bytes
- 内容: リガンド構造とハイライトされたマッチ部分

#### 注意事項
⚠️ **視覚的確認が必要**  
このタスクは自動テストのみでは不十分です。生成された画像を開いて以下を確認してください：
- リガンド構造が明瞭に表示されているか
- マッチ部分が適切にハイライトされているか
- 画像の解像度は十分か

#### 参考コード
- [`scripts/replace_substructure.py:404`](../scripts/replace_substructure.py) - 分子構造の描画例

---

### T4: Atom Matching関数実装

**実装日**: 2025-10-23
**関数**: [`match_substructures()`](../inverse_msmd/substructure_replacement.py:159)
**テストファイル**: [`tests/unit/test_atom_matching.py`](../tests/unit/test_atom_matching.py)

#### 実装内容
```python
def match_substructures(
    mol1: Chem.Mol,
    mol2: Chem.Mol
) -> List[npt.NDArray[np.int_]]:
    from rdkit.Chem import rdFMCS
    
    # 水素を除去
    mol1_no_h = Chem.RemoveHs(mol1)
    mol2_no_h = Chem.RemoveHs(mol2)
    
    # MCS検索
    mcs_result = rdFMCS.FindMCS([mol1_no_h, mol2_no_h])
    if mcs_result.numAtoms == 0:
        return []
    
    mcs = Chem.MolFromSmarts(mcs_result.smartsString)
    
    # すべてのマッチングを取得
    mol1_matches = mol1_no_h.GetSubstructMatches(mcs, uniquify=False)
    mol2_matches = mol2_no_h.GetSubstructMatches(mcs, uniquify=False)
    
    # 重複除外
    seen_pairings = set()
    matches = []
    
    for mol1_match in mol1_matches:
        for mol2_match in mol2_matches:
            pairing = frozenset(zip(mol1_match, mol2_match))
            
            if pairing in seen_pairings:
                continue
            seen_pairings.add(pairing)
            
            atom_pairs = np.array([mol1_match, mol2_match], dtype=np.int_)
            matches.append(atom_pairs)
    
    return matches
```

#### 技術的詳細
- MCS（Maximum Common Substructure）検索を使用
- 重複パターンの自動除外
- 複数のマッチングパターンを全て返却
- NumPy配列形式（shape: 2 x n_atoms）

#### テスト結果
```
E23原子数（水素なし）: 9
E24原子数（水素なし）: 12

✓ Atom Matchingパターン数: 16

パターン 0:
  形状: (2, 8)
  E23側インデックス: [1 0 2 3 4 7 5 6]
  E24側インデックス: [ 3  2  4  5  6 10 11  0]

... (16パターン)

✓ 全ての検証に合格
```

#### データ構造
各マッチングパターンは以下の形式：
```
atom_pairs[0, :] = mol1の原子インデックス
atom_pairs[1, :] = mol2の原子インデックス
```

#### 参考コード
- [`inverse_msmd/alignment.py:102`](../inverse_msmd/alignment.py) - MCS検索とatom matchingの実装例

---

## ✅ Phase 2: 座標変換機能 (完了)

### T5: Superimpose計算関数実装

**実装日**: 2025-10-23
**関数**: [`calculate_transformation()`](../inverse_msmd/substructure_replacement.py:232)
**テストファイル**: [`tests/unit/test_transformation.py`](../tests/unit/test_transformation.py)

#### 実装内容
```python
def calculate_transformation(
    source_coords: npt.NDArray[np.float64],
    target_coords: npt.NDArray[np.float64],
    atom_pairs: npt.NDArray[np.int_]
) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    from inverse_msmd.utils.bio_utils import SuperImposer
    
    # atom_pairsに基づいて対応する座標を抽出
    source_matched = source_coords[atom_pairs[0]]
    target_matched = target_coords[atom_pairs[1]]
    
    # SuperImposerで変換行列を計算
    si = SuperImposer()
    si.fit(source_matched, target_matched)
    
    return si.rot_, si.tran_
```

#### 技術的詳細
- `SuperImposer`クラスを使用してSVD法による最適化
- atom_pairsに基づいて対応する座標を自動抽出
- 回転行列（3x3）と並進ベクトル（3,）を返却
- 変換式: `new_coords = rot @ coords + tran`

#### テスト結果
```
✓ 回転行列の形状が3x3
✓ 並進ベクトルの形状が(3,)
✓ 回転行列が直交行列（R^T * R = I）
✓ 回転行列の行列式が1
✓ 複数のマッチングパターンで動作確認
```

**テスト統計**: 5個のテストケース全てパス

---

### T6: タンパク質変換関数実装

**実装日**: 2025-10-23
**関数**: [`apply_transformation_to_protein()`](../inverse_msmd/substructure_replacement.py:272)
**テストファイル**: [`tests/unit/test_protein_transformation.py`](../tests/unit/test_protein_transformation.py)

#### 実装内容
```python
def apply_transformation_to_protein(
    protein: Structure,
    rot: npt.NDArray[np.float64],
    tran: npt.NDArray[np.float64]
) -> Structure:
    from inverse_msmd.utils.bio_utils import PDB
    
    # タンパク質の座標を取得
    protein_coords = PDB.get_attr(protein, "coord")
    
    # 変換を適用: new_coords = rot @ coords + tran
    transformed_coords = np.dot(protein_coords, rot) + tran
    
    # 変換後の座標を設定
    PDB.set_attr(protein, "coord", transformed_coords)
    
    return protein
```

#### 技術的詳細
- BioPythonのStructureオブジェクトを直接操作
- NumPy配列で効率的な座標変換
- 変換式を正確に適用
- 元のタンパク質構造を変更（in-place）

#### テスト結果
```
✓ 座標の形状が保持される
✓ 座標が実際に変換される
✓ 変換式（rot @ coords + tran）が正確に適用される
✓ PDBファイルとして保存・読み込み可能
```

**テスト統計**: 5個のテストケース全てパス

---

### T7: リガンド置換関数実装

**実装日**: 2025-10-23
**関数**: [`replace_ligand_substructure()`](../inverse_msmd/substructure_replacement.py:305)
**テストファイル**: [`tests/unit/test_ligand_replacement.py`](../tests/unit/test_ligand_replacement.py)

#### 実装内容
- RDKitの`RWMol`を使用した分子編集
- 接続点の自動検出
- 水素を除去した分子での処理
- 置換部分構造の原子と結合を追加
- 古い部分構造の原子を削除

#### 技術的詳細
```python
def replace_ligand_substructure(...) -> Chem.Mol:
    # 水素を除去
    ligand_no_h = Chem.RemoveHs(ligand_mol)
    replacement_no_h = Chem.RemoveHs(replacement_mol)
    
    # 接続点を検出
    connections = []  # [(replacement_atom_idx, ligand_neighbor_idx, bond_type), ...]
    
    # RWMolで分子を編集
    rwmol = Chem.RWMol(copy.deepcopy(ligand_no_h))
    
    # 新しい原子と結合を追加
    # 接続を追加
    # 古い原子を削除（逆順）
    
    return rwmol.GetMol()
```

#### テスト結果
```
✓ 基本的な置換が機能
✓ 原子数が期待値に近い
✓ Sanitizeチェックをパス（化学的に妥当）
✓ SDFファイルとして保存可能
✓ 複数のマッチングパターンで動作
✓ SMILESが変化（置換が実行されている）
✓ コンフォーマー情報が保持される
```

**テスト統計**: 7個のテストケース全てパス

#### 🐛 座標保持バグの修正 (2025-10-23)

**問題**: 置換後のリガンドSDFファイルで、置換部分（E24など）の座標が全て原点(0, 0, 0)になっていた。

**原因**:
- `RWMol.GetMol()`で変換した際、元のリガンドのコンフォーマーが残存
- 新しいコンフォーマーを`AddConformer()`で追加すると、2つのコンフォーマーが存在
- デフォルトでは最初のコンフォーマー（座標が初期化されていない）が使用される

**解決策**:
[`inverse_msmd/substructure_replacement.py:431`](../inverse_msmd/substructure_replacement.py:431)に以下の修正を適用：

```python
# RWMolをMolに変換
new_mol = rwmol.GetMol()

# 既存のコンフォーマーを全て削除  ← 追加した重要な行
new_mol.RemoveAllConformers()

# 座標をコンフォーマーとして設定
from rdkit.Geometry import Point3D
conf = Chem.Conformer(new_mol.GetNumAtoms())
for i, coords in enumerate(new_coords_list):
    conf.SetAtomPosition(i, Point3D(float(coords[0]), float(coords[1]), float(coords[2])))
new_mol.AddConformer(conf)
```

**検証結果**:
```
# 修正前（全て原点）
   0.0000    0.0000    0.0000 C   # E24部分
   0.0000    0.0000    0.0000 C
   ...

# 修正後（正しい座標）
   3.5640   -0.0000    0.0000 C   # E24部分の正しい座標
   2.8600    1.1390    0.3940 C
   1.4660    1.1390    0.3940 C
   ...
```


---

## 🐛 構造重ね合わせの座標系バグ修正 (2025-10-23)

### 問題の発見

統合ワークフロー関数[`integrated_substructure_replacement()`](../inverse_msmd/substructure_replacement.py:583-622)において、構造重ね合わせの座標系の基準が設計仕様と逆になっていることが判明しました。

**問題点**:
- 当初の設計: E24の座標系を基準とし、リガンド+タンパク質を変換
- 実装されていた内容: リガンドの座標系を基準とし、E24のみを変換
- 結果: E24の本来の向き・位置情報が失われ、MSMDプロファイルとの対応が崩れる

### 原因分析

[`inverse_msmd/substructure_replacement.py:583-618`](../inverse_msmd/substructure_replacement.py:583-618)での処理:

```python
# 修正前（誤った実装）
# 1. E24をリガンドに合わせる変換を計算
rot, tran = calculate_transformation(
    to_coords,              # E24（変換元）
    ligand_match_coords,    # リガンド（変換先）
    atom_pairs[::-1]
)

# 2. E24のみを変換
to_mol_transformed_coords = np.dot(to_coords, rot) + tran

# 3. タンパク質は変換なし
protein_copy = copy.deepcopy(protein)
transformed_protein = protein_copy  # 変換なし
```

この実装では:
- ✗ E24がリガンドに合わせて変形される
- ✗ E24の本来の向き・位置が失われる
- ✗ タンパク質が変換されない
- ✗ MSMDプロファイルとの対応関係が崩れる

### 修正内容

[`inverse_msmd/substructure_replacement.py:583-622`](../inverse_msmd/substructure_replacement.py:583-622)を以下のように修正:

```python
# 修正後（正しい実装）
# 1. リガンドをE24に合わせる変換を計算
rot, tran = calculate_transformation(
    ligand_match_coords,    # リガンド部分（変換元）
    to_coords,              # E24（変換先・基準）
    atom_pairs              # 順序も修正
)

# 2. リガンド全体を変換
ligand_transformed_coords = np.dot(ligand_coords, rot) + tran
# コンフォーマーに設定...

# 3. E24は元の座標で使用
replaced_ligand = replace_ligand_substructure(
    ligand_transformed,     # 変換後のリガンド
    selected_match,
    to_no_h,               # 元のE24（変換なし）
    atom_pairs
)

# 4. タンパク質にも同じ変換を適用
protein_copy = copy.deepcopy(protein)
protein_coords = PDB.get_attr(protein_copy, "coord")
protein_transformed_coords = np.dot(protein_coords, rot) + tran
PDB.set_attr(protein_copy, "coord", protein_transformed_coords)
transformed_protein = protein_copy
```

### 修正箇所の詳細

**変更1**: 変換の方向を逆転（行585-589）
- 変換元: E24 → リガンド部分
- 変換先: リガンド → E24
- atom_pairsの順序: [::-1]を削除

**変更2**: 変換対象をE24からリガンドへ変更（行591-604）
- E24の座標変換を削除
- リガンド全体の座標を変換

**変更3**: タンパク質にも変換を適用（行616-622）
- タンパク質座標を変換
- PDB.set_attr()で座標を更新

**変更4**: 置換に使うE24を元の座標に変更（行609-614）
- 変換後のE24ではなく元のE24（to_no_h）を使用
- 変換後のリガンド（ligand_transformed）を使用

### 検証結果

修正後の動作確認:

```bash
python test_coordinate_fix_check.py
```

**検証項目**:
- ✅ E24の座標が保持されている（差: 0.000000 Å）
- ✅ タンパク質が変換されている（移動: 63.7 Å）
- ✅ リガンドが変換されている（移動: 52.7 Å）
- ✅ 16パターン全て正常に生成

### 影響

**修正前の動作**:
- リガンドとタンパク質の相対位置関係は保持

---

## 🐛 MCS重ね合わせ問題の修正 (2025-10-23)

### 問題の発見

統合ワークフロー関数において、E23（置換前）とE24（置換後）の芳香環が正しく重ならない問題が発覚しました。

**症状**:
- E24の2つの芳香環の中央部分にリガンドのE23部分が配置される
- MCS（最大共通部分構造）ベースの重ね合わせがうまく機能していない
- 重ね合わせ後のRMSDが1.96Åと非常に大きい

### 原因分析

**初期のMCS設定（問題あり）**:
```python
# 環内と環外の区別なし
mcs_result = rdFMCS.FindMCS([mol1_no_h, mol2_no_h])
```

この設定では：
- MCS原子数: 8個
- しかし、環内原子と環外原子が混在してマッチング
- E23とE24の3D座標が大きく異なるため、RMSD 1.96Å

### 解決策

[`inverse_msmd/substructure_replacement.py:200`](../inverse_msmd/substructure_replacement.py:200)に修正を適用：

```python
# MCS検索（環内と環外のマッチングを制限）
mcs_result = rdFMCS.FindMCS([mol1_no_h, mol2_no_h], ringMatchesRingOnly=True)
```

**重要なポイント**:
- `ringMatchesRingOnly=True`をキーワード引数として直接渡す
- MCSParametersオブジェクトではなく、キーワード引数を使用

### 修正結果

**修正前**:
- MCS原子数: 8個
- Atom matchingパターン数: 16個
- RMSD: **1.958979 Å**（✗ 問題あり）
- E23とE24の芳香環が重ならない

**修正後**:
- MCS原子数: 6個（より正確な芳香環マッチング）
- Atom matchingパターン数: 24個
- RMSD: **0.021014 Å**（✅ ほぼ完璧）
- E23とE24の芳香環が正確に重なる

### 影響

この修正により：
- ✅ 芳香環の正確な重ね合わせが実現
- ✅ 化学的に意味のある構造置換が可能に
- ✅ RMSDが劇的に改善（1.96Å → 0.021Å）
- ✅ より多くのマッチングパターン（16 → 24）で柔軟性が向上

### 検証

修正後の動作確認：
```bash
python scripts/integrated_replacement.py \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --from-file data/sample_probes/E23 \
    --to-file data/sample_probes/E24 \
    --output output/final/ \
    --match-index 0
```

出力構造をPyMOLで確認すると、E23とE24の芳香環が正確に重なっていることが確認できます。

- しかしE24の向き・位置情報が失われる
- MSMDワークフローでの後続処理に影響

**修正後の動作**:
- E24の座標系を基準とする（設計通り）
- リガンドとタンパク質がE24の位置に移動
- MSMDプロファイルとの対応が維持される
- 科学的に意味のある出力構造

### テストファイル

- [`test_fix_validation.py`](../test_fix_validation.py) - 基本動作確認
- [`test_coordinate_fix_check.py`](../test_coordinate_fix_check.py) - 座標検証

**結論**: 修正により、設計仕様通りの動作となり、E24の座標系を基準とした構造重ね合わせが正しく実行されるようになりました。

**影響**: 全40個のテストが引き続き成功し、座標が正しく保持されることを確認。


---

## 🔧 立体障害チェック機能の追加 (2025-10-23)

### 問題の背景

統合ワークフロー関数において、置換後のリガンドが立体化学的に不適切な構造になる可能性がありました。

**具体的な問題**:
- 置換基の向きによっては、分子内で原子間距離が異常に近くなる
- 立体障害が発生しているパターンが出力に含まれる
- 化学的に不自然な構造が生成される

### 実装内容

#### 1. 立体障害チェック関数の追加

[`inverse_msmd/substructure_replacement.py:447`](../inverse_msmd/substructure_replacement.py:447)に新関数を追加：

```python
def check_steric_clash(
    mol: Chem.Mol,
    min_distance: float = 2.0,
    exclude_bonded: bool = True
) -> Tuple[bool, List[Tuple[int, int, float]]]:
    """
    分子内の立体障害（原子間距離が異常に近い）をチェックします。
    
    Parameters
    ----------
    mol : Chem.Mol
        チェックする分子
    min_distance : float, default=2.0
        許容最小原子間距離（Å）
    exclude_bonded : bool, default=True
        結合している原子対を除外するかどうか
    
    Returns
    -------
    is_valid : bool
        立体障害がない場合True
    clashes : List[Tuple[int, int, float]]
        立体障害がある原子対のリスト
    """
```

**機能**:
- 全原子間距離を計算（scipy.spatial.distance.pdist使用）
- 結合している原子対を除外
- 指定した最小距離（デフォルト2.0Å）より近い原子対を検出
- 詳細な立体障害情報を返す

#### 2. 統合ワークフローへの組み込み

[`inverse_msmd/substructure_replacement.py:682`](../inverse_msmd/substructure_replacement.py:682)で、リガンド置換後に立体障害チェックを実行：

```python
# 5. 立体障害チェック
print(f"  立体障害チェック中...")
is_valid, clashes = check_steric_clash(replaced_ligand, min_distance=2.0)

if not is_valid:
    print(f"  ⚠ 警告: 立体障害を検出しました（{len(clashes)}箇所）")
    for atom_i, atom_j, dist in clashes:
        atom_i_symbol = replaced_ligand.GetAtomWithIdx(atom_i).GetSymbol()
        atom_j_symbol = replaced_ligand.GetAtomWithIdx(atom_j).GetSymbol()
        print(f"    原子{atom_i}({atom_i_symbol}) - 原子{atom_j}({atom_j_symbol}): {dist:.3f}Å")
    print(f"  → このパターンをスキップします")
    continue

print(f"  ✓ 立体障害なし")
```

**動作**:
- 立体障害が検出された場合、そのパターンをスキップ
- 詳細な警告メッセージを表示（原子番号、元素記号、距離）
- 有効なパターンのみを出力

### 検証結果

**テスト実行**:
```bash
python test_integrated_simple.py
```

**結果サマリー**:
- 入力パターン数: 24パターン
- 立体障害検出: 8パターン
- 有効な出力: 16パターン（33%削減）

**立体障害の検出例**:

**パターン3** (13箇所の立体障害):
```
⚠ 警告: 立体障害を検出しました（13箇所）
  原子4(C) - 原子21(C): 1.676Å
  原子4(C) - 原子22(C): 0.981Å
  原子4(C) - 原子23(C): 1.068Å
  原子14(O) - 原子20(C): 0.152Å  ← 深刻な衝突
```

**パターン5** (軽微な立体障害):
```
⚠ 警告: 立体障害を検出しました（1箇所）
  原子5(C) - 原子25(C): 1.876Å
```

### 効果

この機能により：
- ✅ 化学的に不適切な構造が自動的に除外される
- ✅ 出力される構造の品質が向上
- ✅ 後続の計算やシミュレーションでのエラーを防止
- ✅ ユーザーが手動で構造を検証する手間が削減

### 技術詳細

**距離計算方法**:
- `scipy.spatial.distance.pdist`: 全原子対の距離を効率的に計算
- `squareform`: 距離行列に変換
- 結合情報を考慮して、結合原子対を除外

**閾値設定**:
- デフォルト: 2.0Å
- 一般的な原子間の最小距離（van der Waals半径の和）より保守的
- 化学的に妥当な構造のみを許容

**依存関係**:
- `scipy.spatial.distance`: 距離計算ライブラリ（既存の依存関係）

---

## ✅ Phase 3: 統合とインターフェース (完了)

### T8-9: 統合ワークフロー関数実装

**実装日**: 2025-10-23
**関数**: [`integrated_substructure_replacement()`](../inverse_msmd/substructure_replacement.py:423)
**テストファイル**: [`tests/integration/test_workflow.py`](../tests/integration/test_workflow.py)

#### 実装内容
```python
def integrated_substructure_replacement(
    ligand_file: str,
    protein_file: str,
    from_file: str,
    to_file: str,
    output_dir: str,
    match_index: Optional[int] = None
) -> List[Dict[str, str]]:
    # 全ての既存機能を統合
    # 1. ファイル読み込み
    # 2. 部分構造探索
    # 3. マッチ選択（自動/手動/可視化）
    # 4. Atom matching
    # 5. 各パターンについて:
    #    - Superimpose計算
    #    - タンパク質変換
    #    - リガンド置換
    #    - ファイル出力
```

#### 主要機能
- ファイルの自動読み込み（SDF, PDB, PDB+SMI）
- 部分構造の自動探索とマッチ選択
- 複数マッチ時の可視化画像出力
- 全atom matchingパターンの処理
- 結果の自動保存（命名規則に従う）
- 詳細な進捗表示

#### テスト結果
```
✓ 16パターンの結果を生成
✓ 全てのファイルが正常に作成
✓ リガンド: pattern_N_ligand_replaced.sdf (2.4KB)
✓ タンパク質: pattern_N_protein_aligned.pdb (97KB)
✓ 統合テストをパス
```

**テスト統計**: 5個のテストケース全てパス

---

### T10: CLIスクリプト作成

**実装日**: 2025-10-23
**ファイル**: [`scripts/integrated_replacement.py`](../scripts/integrated_replacement.py)

#### 実装内容
完全なコマンドラインインターフェースを実装：

```bash
python scripts/integrated_replacement.py \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --from-file data/sample_probes/E23 \
    --to-file data/sample_probes/E24 \
    --output test_output/cli_test/ \
    --match-index 0 \
    --verbose
```

#### 機能
- **必須オプション**: ligand, protein, from-file, to-file, output
- **オプション**: match-index, verbose, version
- **入力検証**: 全ての必須ファイルの存在確認
- **エラーハンドリング**: 適切なエラーメッセージ
- **ヘルプ**: 詳細な使用例とドキュメント

#### テスト結果
```
✓ ヘルプメッセージが適切に表示
✓ 全てのオプションが機能
✓ ファイル検証が正常動作
✓ 詳細出力（--verbose）が正常動作
✓ 16パターンの結果を正常生成
```

---

### T11: 出力機能とテスト

**実装日**: 2025-10-23
**成果物**: 完全な出力機能とテストスイート

#### 実装内容

**ファイル命名規則**:
- リガンド: `pattern_N_ligand_replaced.sdf`
- タンパク質: `pattern_N_protein_aligned.pdb`
- 可視化: `substructure_matches.png` (複数マッチ時)

**出力ディレクトリ管理**:
- 自動作成（存在しない場合）
- 相対/絶対パス両対応

#### テスト結果
```
✓ 出力ディレクトリが自動作成される
✓ ファイル名が規則に従う
✓ 全パターンで出力ファイル生成
✓ ファイルサイズが適切
✓ ファイル形式が正しい（SDF, PDB）
✓ 内容が読み込み可能
```

**生成ファイル例**:
```
test_output/cli_test/
├── pattern_0_ligand_replaced.sdf    (2.4KB)
├── pattern_0_protein_aligned.pdb    (97KB)
├── pattern_1_ligand_replaced.sdf
├── pattern_1_protein_aligned.pdb
...
└── pattern_15_protein_aligned.pdb
```

---

## ✅ Phase 4: 品質保証 (完了)

### T12: 総合テスト実行

**実施日**: 2025-10-23
**テストスイート**: pytest (40テストケース)

#### テスト結果

**全テスト成功**: 40/40 (100%)

```bash
./run_tests.sh
========================================
  inverse_msmd テストスイート
========================================

全てのテストを実行中...
===================================================== 40 passed in 2.28s ======================================================
```

#### テストカバレッジ

**単体テスト (30テスト)**:
- ✅ インポートテスト (3テスト)
- ✅ 部分構造探索テスト (4テスト)
- ✅ 可視化テスト (3テスト)
- ✅ Atom Matchingテスト (5テスト)
- ✅ 座標変換テスト (5テスト)
- ✅ タンパク質変換テスト (5テスト)
- ✅ リガンド置換テスト (7テスト)

**統合テスト (10テスト)**:
- ✅ 部分構造探索とマッチングワークフロー
- ✅ 変換行列計算ワークフロー
- ✅ 可視化統合テスト
- ✅ 複数プローブワークフロー
- ✅ 統合部分構造置換（基本）
- ✅ match_index指定テスト
- ✅ match_index未指定テスト
- ✅ 無効なmatch_indexエラーハンドリング
- ✅ その他統合テスト

#### 動作確認済み機能

**コア機能**:
- [x] 部分構造の自動探索
- [x] 複数マッチの可視化
- [x] Atom matchingパターン生成（16パターン）
- [x] Superimpose計算
- [x] タンパク質座標変換
- [x] リガンド部分構造置換
- [x] ファイル自動出力

**CLIツール**:
- [x] 全オプション動作確認
- [x] 入力ファイル検証
- [x] エラーメッセージ表示
- [x] ヘルプ表示
- [x] 詳細出力モード

---

### T13: エラーハンドリングとバリデーション

**実施日**: 2025-10-23

#### 実装済みエラーハンドリング

**入力検証**:
- ✅ ファイル存在確認（リガンド、タンパク質、部分構造）
- ✅ match_indexの範囲チェック
- ✅ 部分構造マッチの存在確認
- ✅ Atom matchingパターンの存在確認

**エラーメッセージ**:
```python
# 無効なmatch_index
ValueError: 無効なmatch_index: 999。有効範囲は 0 から 0 です

# 部分構造が見つからない
ValueError: リガンド中に部分構造が見つかりませんでした

# Atom matching失敗
ValueError: Atom matchingに失敗しました
```

**CLIでの入力検証**:
```python
# 必須ファイルの存在確認
errors = []
if not ligand_path.exists():
    errors.append(f"リガンドファイルが見つかりません: {ligand_path}")
# ... 全ファイルチェック

if errors:
    print("エラー: 必要なファイルが見つかりません\n", file=sys.stderr)
    for error in errors:
        print(f"  - {error}", file=sys.stderr)
    sys.exit(1)
```

#### テスト確認済み

- ✅ 無効なmatch_indexでValueErrorが発生
- ✅ 存在しないファイルで適切なエラーメッセージ
- ✅ 部分構造未検出時のエラーハンドリング
- ✅ Atom matching失敗時のエラーハンドリング

---

### T14: ドキュメント整備

**実施日**: 2025-10-23

#### 更新されたドキュメント

**メインドキュメント**:
- ✅ [`README.md`](../README.md) - 統合機能の使用例とAPI説明
- ✅ [`scripts/README.md`](../scripts/README.md) - CLIツールの詳細ガイド
- ✅ [`docs/implementation_progress.md`](implementation_progress.md) - 本ドキュメント（完全な実装記録）

**既存ドキュメント**:
- ✅ [`docs/integrated_replacement_plan.md`](integrated_replacement_plan.md) - 設計仕様
- ✅ [`docs/task_handoff_guide.md`](task_handoff_guide.md) - タスク引き継ぎ
- ✅ [`docs/testing_checklist.md`](testing_checklist.md) - テスト手順
- ✅ [`docs/testing_responsibility.md`](testing_responsibility.md) - テスト責任分担

#### ドキュメントの内容

**使用例**:
- Python APIからの使用方法
- CLIツールの使用方法
- 出力ファイルの説明
- トラブルシューティング

**技術詳細**:
- 各関数の実装内容
- データフロー
- テスト結果
- 参考コード

---

## 🎉 実装完了サマリー

### 達成された目標

**全14タスクを完了** (100%)

1. ✅ **Phase 1** (4タスク): 基本機能実装
2. ✅ **Phase 2** (3タスク): 座標変換機能
3. ✅ **Phase 3** (3タスク): 統合とインターフェース
4. ✅ **Phase 4** (3タスク): 品質保証
5. ✅ **T14**: ドキュメント整備

### 主要成果物

**実装ファイル**:
- [`inverse_msmd/substructure_replacement.py`](../inverse_msmd/substructure_replacement.py) - 7関数、約600行
- [`scripts/integrated_replacement.py`](../scripts/integrated_replacement.py) - CLIツール、約180行

**テストスイート**:
- 単体テスト: 30テストケース
- 統合テスト: 10テストケース
- **合計**: 40テストケース（全て成功）

**ドキュメント**:
- README更新
- scripts/README作成
- 実装進捗記録（本ドキュメント）

### 動作確認

**実行例**:
```bash
python scripts/integrated_replacement.py \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --from-file data/sample_probes/E23 \
    --to-file data/sample_probes/E24 \
    --output output/integrated/ \
    --verbose
```

**出力**:
- 16パターンのリガンド（各2.4KB）
- 16パターンのタンパク質（各97KB）
- 全ファイルが正常に生成

### 次のステップ

統合部分構造置換機能の実装は完了しました。今後の改善案：

1. **パフォーマンス最適化**: 大規模分子での処理速度向上
2. **追加機能**:
   - 複数の部分構造を一度に置換
   - カスタムマッチングルール
3. **ユーザーインターフェース**:
   - GUIツールの開発
   - 可視化機能の強化

---

## 📚 参考リンク

### プロジェクト内部
- [全体設計](integrated_replacement_plan.md)
- [タスク引き継ぎ](task_handoff_guide.md)
- [テストチェックリスト](testing_checklist.md)
- [メインREADME](../README.md)
- [CLIツールガイド](../scripts/README.md)

### 実装ファイル
- [`inverse_msmd/substructure_replacement.py`](../inverse_msmd/substructure_replacement.py)
- [`scripts/integrated_replacement.py`](../scripts/integrated_replacement.py)
- [`inverse_msmd/alignment.py`](../inverse_msmd/alignment.py)
- [`inverse_msmd/utils/bio_utils.py`](../inverse_msmd/utils/bio_utils.py)
- [`inverse_msmd/utils/mol_utils.py`](../inverse_msmd/utils/mol_utils.py)

### テストスイート
- [`tests/`](../tests/) - pytest形式の単体テスト・統合テスト
- [`run_tests.sh`](../run_tests.sh) - テスト実行スクリプト

---

**最終更新**: 2025-10-23
**実装者**: AI Assistant (Roo)
**ステータス**: ✅ 完了 (100%)

---

## 🔔 重要な未実装機能テスト

### match_index機能の動作チェック（Phase 3完了後）

**テスト時期**: T8-9（統合ワークフロー）完成後

**目的**: 複数マッチが見つかった場合に、ユーザーが特定のマッチを選択できることを確認

#### テストケース準備完了
- **テストデータ**: A01を使用（リガンド中に2箇所マッチ）
- **可視化画像**: [`test_output/substructure_matches_A01_test.png`](../test_output/substructure_matches_A01_test.png)
- **マッチ情報**:
  - Match 0: 原子インデックス (2, 3, 5, 19, 18, 4)
  - Match 1: 原子インデックス (6, 12, 15, 13, 7, 14)

#### 実施予定のテスト

**テスト1: Match 0を選択**
```python
results_match0 = integrated_substructure_replacement(
    ligand_file="data/atom_matching/4hw3_A_lig.sdf",
    protein_file="data/sample_proteins/4hw3_A.pdb",
    from_file="data/sample_probes/A01",
    to_file="data/sample_probes/E24",
    output_dir="test_output/match0/",
    match_index=0  # Match 0を指定
)
```

**テスト2: Match 1を選択**
```python
results_match1 = integrated_substructure_replacement(
    ligand_file="data/atom_matching/4hw3_A_lig.sdf",
    protein_file="data/sample_proteins/4hw3_A.pdb",
    from_file="data/sample_probes/A01",
    to_file="data/sample_probes/E24",
    output_dir="test_output/match1/",
    match_index=1  # Match 1を指定
)
```

#### 確認すべき違い

1. **リガンド構造**
   - Match 0とMatch 1で置換される部分構造の位置が異なる
   - 出力SDFファイルの原子座標が異なる

2. **タンパク質構造**
   - 座標変換の基準点が異なるため、タンパク質の配置が異なる
   - 出力PDBファイルの座標が異なる

3. **視覚的確認（PyMOL等）**
   - 2つの出力を重ねて表示
   - リガンドの置換位置が異なることを確認
   - タンパク質の配置が異なることを確認

#### テスト実行スクリプト（作成予定）
- `test_match_selection_functionality.py` - match_index機能の包括的テスト
- 両方のmatch_indexで実行し、出力を比較
- 差分が適切であることを検証

#### 参考資料
- [`test_match_selection_demo.py`](../test_match_selection_demo.py) - マッチ選択方法のデモ
- [`docs/visual_check_guide.md`](visual_check_guide.md) - 視覚的確認ガイド

---


### T5: Superimpose計算関数

**状態**: 未実装  
**予定関数**: `calculate_transformation()`

#### 実装予定
- `SuperImposer`クラスを使用
- 回転行列（3x3）と並進ベクトル（3,）を計算
- 変換式: `new_coords = rot @ coords + tran`

#### 参考コード
- [`inverse_msmd/utils/bio_utils.py:39`](../inverse_msmd/utils/bio_utils.py) - SuperImposerクラス

---

### T6: タンパク質変換関数

**状態**: 未実装  
**予定関数**: `apply_transformation_to_protein()`

#### 実装予定
- 変換行列をタンパク質全体に適用
- BioPythonのPDB操作を使用
- 座標の保持と更新

#### 参考コード
- [`inverse_msmd/alignment.py:191`](../inverse_msmd/alignment.py) - タンパク質座標変換の例

---

### T7: リガンド置換関数

**状態**: 未実装  
**予定関数**: `replace_ligand_substructure()`

#### 実装予定
- 既存の`create_replacement()`を参考
- 結合情報の保持
- 分子のSanitize

#### 参考コード
- [`scripts/replace_substructure.py:106`](../scripts/replace_substructure.py) - 部分構造置換の実装例

---

## ⏳ Phase 3: 統合とインターフェース (未着手)

### T8-9: 統合ワークフロー関数

**状態**: 未実装  
**予定関数**: `integrated_substructure_replacement()`

#### 実装予定
- T1-T7の全機能を統合
- ファイル入出力の管理
- 複数パターンの処理

---

### T10: CLIスクリプト

**状態**: 未実装  
**予定ファイル**: `scripts/integrated_replacement.py`

#### 実装予定
- argparseによるコマンドライン引数処理
- ヘルプメッセージ
- エラーハンドリング

---

### T11: 出力機能

**状態**: 未実装

#### 実装予定
- ファイル命名規則の実装
- SDF/PDB出力
- ログ出力

---

## ⏳ Phase 4: 品質保証 (未着手)

### T12: 総合テスト

**状態**: 未実装

#### 予定
- エンドツーエンドテスト
- 実データでの検証
- PyMOLでの可視化確認

---

### T13: エラーハンドリング

**状態**: 未実装

#### 予定
- 入力検証
- 例外処理
- エラーメッセージ

---

### T14: ドキュメント

**状態**: 一部完了

#### 完了
- [`docs/integrated_replacement_plan.md`](integrated_replacement_plan.md) - 全体設計
- [`docs/task_handoff_guide.md`](task_handoff_guide.md) - タスク引き継ぎ
- [`docs/testing_checklist.md`](testing_checklist.md) - テスト方法
- [`docs/testing_responsibility.md`](testing_responsibility.md) - テスト責任分担
- 本ドキュメント - 進捗記録

#### 未完了
- README更新
- 使用例の追加
- FAQ作成

---

## 📝 技術メモ

### 使用したライブラリとバージョン
- RDKit
- BioPython
- NumPy
- Matplotlib

### コーディング規約
- 型ヒントを必ず付与
- Google形式のdocstring
- PEP 8準拠

### ファイル命名規則
- テストファイル: `test_*.py`
- 出力ファイル: `test_output/`

---

## 🔜 次回作業の推奨事項

### 優先順位1: Phase 2の完了
1. **T5: calculate_transformation()** - Superimpose計算
2. **T6: apply_transformation_to_protein()** - タンパク質変換
3. **T7: replace_ligand_substructure()** - リガンド置換

これら3つは依存関係がないため、並行実装可能です。

### 優先順位2: Phase 3の開始
4. **T8-9: integrated_substructure_replacement()** - 統合ワークフロー
   - Phase 2の完了が必要

### 準備事項
- Phase 2の各関数用テストデータの確認
- PyMOL等の可視化ツールの準備（手動確認用）

---

## 📚 参考リンク

### プロジェクト内部
- [全体設計](integrated_replacement_plan.md)
- [タスク引き継ぎ](task_handoff_guide.md)
- [テストチェックリスト](testing_checklist.md)

### 実装ファイル
- [`inverse_msmd/substructure_replacement.py`](../inverse_msmd/substructure_replacement.py)
- [`inverse_msmd/alignment.py`](../inverse_msmd/alignment.py)
- [`scripts/replace_substructure.py`](../scripts/replace_substructure.py)

### テストスイート
- [`tests/`](../tests/) - pytest形式の単体テスト・統合テスト
- [`tests/unit/test_imports.py`](../tests/unit/test_imports.py)
- [`tests/unit/test_substructure_search.py`](../tests/unit/test_substructure_search.py)
- [`tests/unit/test_visualization.py`](../tests/unit/test_visualization.py)
- [`tests/unit/test_atom_matching.py`](../tests/unit/test_atom_matching.py)
- [`tests/integration/test_workflow.py`](../tests/integration/test_workflow.py)
- [`run_tests.sh`](../run_tests.sh) - テスト実行スクリプト

### テスト実行方法
```bash
# 全テスト実行
./run_tests.sh all

# 単体テストのみ
./run_tests.sh unit

# 詳細はtests/README.mdを参照
```