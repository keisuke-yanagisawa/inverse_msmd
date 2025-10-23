# プロファイルスコア計算統合機能 実装計画

> **注意:** このドキュメントは統合計画の概要と要求仕様を提供します。
> 
> - 詳細なアーキテクチャ → [`architecture.md`](architecture.md)
> - 実装進捗とタスクリスト → [`implementation_progress.md`](implementation_progress.md)
> - テスト手順 → [`testing_guide.md`](testing_guide.md)
> - タスク引き継ぎ → [`task_handoff.md`](task_handoff.md)

## 概要

`integrated_substructure_replacement`機能を拡張し、プロファイルマッチングスコア計算を統合して一括処理を可能にします。

**主な目標:**
- 部分構造置換とスコア計算を一つのワークフローで実行
- 複数のatom matchingパターンをスコアでランク付け
- 後方互換性を維持（スコア計算はオプション）
- 個別機能のエクスポートでカスタムワークフローに対応

## 要求仕様

### 入力

既存の`integrated_substructure_replacement`の入力に加えて：

| パラメータ | 型 | デフォルト | 説明 |
|-----------|-----|----------|------|
| `profile_dir` | Optional[str] | None | プロファイルディレクトリパス |
| `probe_id` | Optional[str] | None | プローブID（例: "E24"） |
| `gamma` | float | 0.0 | 距離重み付けパラメータ |

**プロファイルファイル形式:**
- ファイル名: `{probe_id}_{残基名}_profile.dx.gz`
- 例: `E24_ALA_profile.dx.gz`, `E24_ARG_profile.dx.gz`, ...
- 19種類（GLYを除く標準アミノ酸）

### 処理フロー

```
1. 部分構造置換処理（既存機能）
   ├─ リガンド中の部分構造探索
   ├─ Atom matching
   ├─ 各パターンの座標変換と置換
   └─ 立体障害チェック

2. プロファイルスコア計算（新規機能、オプション）
   ├─ 各パターンの置換後タンパク質構造を評価
   ├─ Cβ原子位置でプロファイル値を3D補間
   ├─ 距離重み付け（gamma > 0の場合）
   └─ 対数スケールでスコア統合

3. 結果の整理と出力
   ├─ 全パターンにスコアを付与（計算時のみ）
   ├─ スコアで降順ソート（計算時のみ）
   └─ ファイル出力とスコア情報の返却
```

### 出力

各atom matchingパターンごとに：

**ファイル出力:**
- `pattern_N_ligand_replaced.sdf`: 置換後のリガンド構造
- `pattern_N_protein_aligned.pdb`: 座標変換後のタンパク質構造

**戻り値:**
```python
[
    {
        'ligand_file': 'path/to/pattern_0_ligand_replaced.sdf',
        'protein_file': 'path/to/pattern_0_protein_aligned.pdb',
        'pattern_index': 0,
        'score': -125.43  # スコア計算時のみ
    },
    # ... more patterns (スコア降順ソート済み)
]
```

## 実装フェーズ

プロジェクトは5つのフェーズに分かれています：

1. **Phase 1**: プロファイルスコア計算モジュール作成
   - `inverse_msmd/profile_scoring.py`の実装

2. **Phase 2**: 統合ワークフローの拡張
   - `integrated_substructure_replacement()`にスコア計算を統合

3. **Phase 3**: パッケージエクスポートの更新
   - `inverse_msmd/__init__.py`で個別機能をエクスポート

4. **Phase 4**: CLIスクリプトの拡張
   - `scripts/integrated_replacement.py`に新規オプション追加

5. **Phase 5**: テストとドキュメント
   - テストスイートの作成
   - ドキュメントの更新

詳細は[`implementation_progress.md`](implementation_progress.md)を参照してください。

## 使用例

### 1. 統合ワークフロー（プロファイルスコアあり）

```python
from inverse_msmd.substructure_replacement import integrated_substructure_replacement

results = integrated_substructure_replacement(
    ligand_file="data/atom_matching/4hw3_A_lig.sdf",
    protein_file="data/sample_proteins/4hw3_A.pdb",
    from_file="data/sample_probes/E23",
    to_file="data/sample_probes/E24",
    output_dir="output/integrated/",
    profile_dir="data/profiles/",  # プロファイル計算を有効化
    probe_id="E24",                # 必須（profile_dir指定時）
    gamma=0.0                      # 距離重み付けなし
)

# 結果はスコアで降順ソート済み
for i, result in enumerate(results):
    print(f"パターン {i}: スコア={result['score']:.2f}")
    print(f"  リガンド: {result['ligand_file']}")
    print(f"  タンパク質: {result['protein_file']}")
```

### 2. CLIから実行

```bash
# スコア計算なし（既存の動作）
python scripts/integrated_replacement.py \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --from-file data/sample_probes/E23 \
    --to-file data/sample_probes/E24 \
    --output output/integrated/

# スコア計算あり
python scripts/integrated_replacement.py \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --from-file data/sample_probes/E23 \
    --to-file data/sample_probes/E24 \
    --output output/integrated/ \
    --profile-dir data/profiles/ \
    --probe-id E24 \
    --gamma 0.0 \
    --verbose
```

### 3. 個別機能の使用（カスタムワークフロー）

```python
from inverse_msmd import (
    find_substructure_in_ligand,
    match_substructures,
    calculate_transformation,
    calculate_profile_score
)

# カスタムワークフロー例
matches = find_substructure_in_ligand(ligand_mol, from_mol)
atom_pairs = match_substructures(from_mol, to_mol)
rot, tran = calculate_transformation(from_mol, to_mol, atom_pairs[0])

# プロファイルスコアの個別計算
score = calculate_profile_score(
    protein, probe_center,
    "data/profiles/", "E24", gamma=0.0
)
```

## 設計上の重要な決定

### 1. 後方互換性の維持

```python
# profile_dir=Noneで既存の動作を維持
results = integrated_substructure_replacement(
    ...,
    profile_dir=None  # スコア計算をスキップ
)
# -> 'score'キーは含まれない
```

### 2. オプションパラメータ

- `profile_dir`が`None`の場合、スコア計算をスキップ
- `profile_dir`指定時は`probe_id`が必須
- `gamma`はデフォルト0.0（重み付けなし）

### 3. エラーハンドリング

```python
# パラメータバリデーション
if profile_dir is not None and probe_id is None:
    raise ValueError("profile_dirが指定されている場合、probe_idも必須です")

# プロファイルファイルの存在確認
if not profiles:
    raise ValueError(f"プロファイルファイルが見つかりません: {profile_dir}")

# Cβ原子の存在確認
if not cb_atoms:
    raise ValueError("Cβ原子が見つかりません")
```

## プロファイルスコア計算のアルゴリズム

**基本的な流れ:**

1. 19種類のアミノ酸プロファイルを読み込み（GLY除く）
2. タンパク質からCβ原子を抽出
3. 各Cβ原子について：
   - 残基タイプに対応するプロファイルを選択
   - Cβ座標でプロファイル値を3D線形補間
   - プローブ中心からの距離を計算
   - 距離重み付け: `w = exp(-gamma * d^2)`
   - 重み付き対数スコアを累積: `score += log(value) * w`

**数式:**

```
score = Σ log(profile_value_i) * exp(-gamma * d_i^2)
```

ここで：
- `profile_value_i`: i番目のCβ原子位置でのプロファイル値（3D補間）
- `d_i`: i番目のCβ原子とプローブ中心の距離
- `gamma`: 距離重み付けパラメータ（0.0で重み付けなし）

詳細は[`architecture.md`](architecture.md)を参照してください。

## 技術スタック

| コンポーネント | ライブラリ | 用途 |
|--------------|-----------|------|
| プロファイル読み込み | gridData | .dx.gz形式の読み込みと3D補間 |
| タンパク質処理 | Bio.PDB | PDB構造の読み書き |
| 分子構造処理 | RDKit | 部分構造探索、分子操作 |
| 数値計算 | NumPy | 行列演算、座標変換 |

## 参考リンク

### ドキュメント
- [`README.md`](README.md) - ドキュメント体系の概要
- [`architecture.md`](architecture.md) - 詳細なアーキテクチャ設計
- [`implementation_progress.md`](implementation_progress.md) - タスクリストと進捗
- [`testing_guide.md`](testing_guide.md) - テスト手順とチェックリスト
- [`task_handoff.md`](task_handoff.md) - AI間の引き継ぎガイド

### コード
- [`inverse_msmd/substructure_replacement.py`](../../inverse_msmd/substructure_replacement.py) - 既存の統合ワークフロー
- [`examples/calculate_matching.py`](../../examples/calculate_matching.py) - プロファイルスコア計算の元実装
- [`scripts/integrated_replacement.py`](../../scripts/integrated_replacement.py) - CLIスクリプト

### プロジェクト全体
- [`../../README.md`](../../README.md) - プロジェクトメインドキュメント
- [`../documentation_best_practices.md`](../documentation_best_practices.md) - ドキュメントベストプラクティス

---

**最終更新:** 2025-10-23  
**ステータス:** Phase 0完了（ドキュメント整備）、Phase 1未着手