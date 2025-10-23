# 統合部分構造置換機能 実装進捗記録

このドキュメントは、実装の進捗状況と各タスクの詳細な記録を提供します。

**最終更新**: 2025-10-23

---

## 📊 進捗サマリー

| フェーズ | 完了 | 進行中 | 未着手 | 進捗率 |
|---------|------|--------|--------|--------|
| Phase 1: 基本機能実装 | 4 | 0 | 0 | 100% ✅ |
| Phase 2: 座標変換機能 | 0 | 0 | 3 | 0% ⏳ |
| Phase 3: 統合とインターフェース | 0 | 0 | 3 | 0% ⏳ |
| Phase 4: 品質保証 | 0 | 0 | 3 | 0% ⏳ |
| **合計** | **4** | **0** | **10** | **29%** |

---

## ✅ Phase 1: 基本機能実装 (完了)

### T1: モジュール基本構造作成

**実装日**: 2025-10-23  
**ファイル**: [`inverse_msmd/substructure_replacement.py`](../inverse_msmd/substructure_replacement.py)  
**テストファイル**: [`test_imports.py`](../test_imports.py)

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
**テストファイル**: [`test_find_substructure.py`](../test_find_substructure.py)

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
**テストファイル**: [`test_visualize_matches.py`](../test_visualize_matches.py)

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
**テストファイル**: [`test_atom_matching.py`](../test_atom_matching.py)

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

## ⏳ Phase 2: 座標変換機能 (未着手)

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

### テストファイル
- [`test_imports.py`](../test_imports.py)
- [`test_find_substructure.py`](../test_find_substructure.py)
- [`test_visualize_matches.py`](../test_visualize_matches.py)
- [`test_atom_matching.py`](../test_atom_matching.py)