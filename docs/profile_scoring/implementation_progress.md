# プロファイルスコア計算機能 実装進捗記録

最終更新: 2025-10-23

## 📊 進捗サマリー

| フェーズ | ステータス | 進捗率 | 備考 |
|---------|----------|-------|------|
| Phase 1: プロファイルスコア計算モジュール | ⏸️ 未着手 | 0% | - |
| Phase 2: 統合ワークフローの拡張 | ⏸️ 未着手 | 0% | - |
| Phase 3: パッケージエクスポートの更新 | ⏸️ 未着手 | 0% | - |
| Phase 4: CLIスクリプトの拡張 | ⏸️ 未着手 | 0% | - |
| Phase 5: テストとドキュメント | ⏸️ 未着手 | 0% | - |
| **全体進捗** | **⏸️ 未着手** | **0%** | - |

## Phase 1: プロファイルスコア計算モジュール作成

### T1: `inverse_msmd/profile_scoring.py`の基本構造作成

**ステータス:** ⏸️ 未着手

**タスク内容:**
- [ ] モジュールdocstringの作成
- [ ] 必要なインポート（gridData, numpy, Bio.PDB等）の追加
- [ ] ファイル作成と基本構造の準備

**実装詳細:**
```python
"""
プロファイルマッチングスコア計算モジュール

このモジュールは、タンパク質構造とプロファイルデータから
マッチングスコアを計算する機能を提供します。
"""

from gridData import Grid
import numpy as np
from pathlib import Path
from typing import Dict
from Bio.PDB.Structure import Structure
```

**完了条件:**
- [ ] ファイルが正しく作成されている
- [ ] 必要なライブラリがインポートされている
- [ ] docstringが適切に記載されている

**注意事項:**
- gridDataライブラリが依存関係に含まれているか確認
- Bio.PDBのインポートパスを確認

---

### T2: `calculate_profile_score()`関数の実装

**ステータス:** ⏸️ 未着手

**タスク内容:**
- [ ] 関数シグネチャとdocstringの作成
- [ ] プロファイル読み込みロジックの実装
- [ ] Cβ原子抽出ロジックの実装
- [ ] スコア計算ロジック（3D補間、重み付け）の実装
- [ ] エラーハンドリングの追加

**関数シグネチャ:**
```python
def calculate_profile_score(
    protein: Structure,
    probe_center: np.ndarray,
    profile_dir: str,
    probe_id: str,
    gamma: float = 0.0
) -> float:
    """
    タンパク質構造とプローブ中心からマッチングスコアを計算
    
    Parameters
    ----------
    protein : Bio.PDB.Structure.Structure
        評価対象のタンパク質構造
    probe_center : np.ndarray, shape (3,)
        プローブ分子の中心座標
    profile_dir : str
        プロファイルファイルのディレクトリパス
    probe_id : str
        プローブID（ファイル名プレフィックス、例: "E24"）
    gamma : float, default=0.0
        距離重み付けパラメータ（0.0=重み付けなし）
    
    Returns
    -------
    float
        対数マッチングスコア
    """
```

**実装のポイント:**
1. プロファイル読み込み
   - 19種類のアミノ酸（GLY除く）のプロファイルを読み込み
   - ファイル名: `{probe_id}_{resname}_profile.dx.gz`
   
2. Cβ原子抽出
   - `atom.get_name() == "CB"`でフィルタ
   - GLY残基は自動的にスキップ

3. スコア計算
   - 3D線形補間: `profiles[resname].interpolated([x], [y], [z])[0]`
   - 距離重み: `w = exp(-gamma * d^2)`
   - スコア累積: `score += log(value) * weight`

**完了条件:**
- [ ] 関数が正しく動作する
- [ ] プロファイルファイルが正しく読み込まれる
- [ ] スコア計算が正確に行われる
- [ ] エラーハンドリングが適切に機能する

**テストケース:**
- [ ] E24プロファイルで4hw3_Aタンパク質のスコアを計算
- [ ] gamma=0.0とgamma=0.003で異なる結果が得られる
- [ ] GLY残基のみのタンパク質でエラーが発生しない

**技術的課題:**
- プロファイル値が負になる場合の処理
- グリッド外の座標の処理

**解決策:**
- 負の値は`profiles[resname].grid.min()`で置換
- グリッド外は補間関数が自動的に処理

---

## Phase 2: 統合ワークフローの拡張

### T3: `integrated_substructure_replacement()`の拡張

**ステータス:** ⏸️ 未着手

**タスク内容:**
- [ ] 新規パラメータ追加（profile_dir, probe_id, gamma）
- [ ] パラメータバリデーションの実装
- [ ] スコア計算の統合（オプション）
- [ ] 結果にスコア情報を追加
- [ ] スコアでソート機能の実装

**新規パラメータ:**
```python
profile_dir: Optional[str] = None,  # プロファイルディレクトリパス
probe_id: Optional[str] = None,     # プローブID
gamma: float = 0.0                  # 距離重み付けパラメータ
```

**実装のポイント:**
1. パラメータバリデーション
   ```python
   calculate_scores = profile_dir is not None
   if calculate_scores and probe_id is None:
       raise ValueError("profile_dirが指定されている場合、probe_idも必須です")
   ```

2. プローブ中心座標の計算
   ```python
   if calculate_scores:
       to_no_h = Chem.RemoveHs(to_mol)
       to_center = to_no_h.GetConformer().GetPositions().mean(axis=0)
   ```

3. 各パターンでのスコア計算
   ```python
   if calculate_scores:
       transformed_center = np.dot(to_center, rot) + tran
       score = calculate_profile_score(
           protein_copy, transformed_center,
           profile_dir, probe_id, gamma
       )
       result['score'] = score
   ```

4. スコアソート
   ```python
   if calculate_scores and results:
       results.sort(key=lambda x: x['score'], reverse=True)
   ```

**完了条件:**
- [ ] 後方互換性が維持されている（profile_dir=Noneで従来通り動作）
- [ ] スコア計算が正しく統合されている
- [ ] 結果がスコアで降順ソートされる
- [ ] 適切なエラーメッセージが表示される

**テストケース:**
- [ ] スコア計算なし（profile_dir=None）で従来通り動作
- [ ] スコア計算あり（profile_dir指定）で全パターンにスコアが付与
- [ ] 結果が正しくソートされる
- [ ] probe_idなしでエラーが発生

---

## Phase 3: パッケージエクスポートの更新

### T4: `inverse_msmd/__init__.py`の更新

**ステータス:** ⏸️ 未着手

**タスク内容:**
- [ ] 個別機能のインポート追加
- [ ] `calculate_profile_score`のインポート
- [ ] `__all__`リストの更新
- [ ] docstringの更新
- [ ] バージョン番号の更新（0.1.0 → 0.2.0）

**実装内容:**
```python
from .substructure_replacement import (
    find_substructure_in_ligand,
    visualize_multiple_matches,
    match_substructures,
    calculate_transformation,
    apply_transformation_to_protein,
    replace_ligand_substructure,
    check_steric_clash,
    integrated_substructure_replacement
)
from .profile_scoring import calculate_profile_score

__version__ = "0.2.0"
```

**完了条件:**
- [ ] 全ての個別機能がエクスポートされている
- [ ] `from inverse_msmd import calculate_profile_score`が動作する
- [ ] バージョン番号が更新されている
- [ ] docstringが最新の内容を反映している

---

## Phase 4: CLIスクリプトの拡張

### T5: `scripts/integrated_replacement.py`の更新

**ステータス:** ⏸️ 未着手

**タスク内容:**
- [ ] 新規コマンドライン引数の追加
- [ ] パラメータバリデーションの実装
- [ ] 統合ワークフロー呼び出しの更新
- [ ] 結果表示の改善（スコア情報含む）
- [ ] ヘルプメッセージの更新

**新規引数:**
```python
parser.add_argument("--profile-dir", default=None, 
    help="プロファイルディレクトリのパス")
parser.add_argument("--probe-id", default=None,
    help="プローブID（例: E24）")
parser.add_argument("--gamma", type=float, default=0.0,
    help="距離重み付けパラメータ")
```

**完了条件:**
- [ ] 新規引数が正しく動作する
- [ ] パラメータバリデーションが機能する
- [ ] スコア情報が適切に表示される
- [ ] ヘルプメッセージが分かりやすい

---

## Phase 5: テストとドキュメント

### T6: テストコードの作成

**ステータス:** ⏸️ 未着手

**タスク内容:**
- [ ] `tests/unit/test_profile_scoring.py`: プロファイルスコア計算のテスト
- [ ] `tests/integration/test_profile_integration.py`: 統合ワークフローのテスト
- [ ] エッジケースのテスト

**テストケース:**
1. **ユニットテスト（test_profile_scoring.py）**
   - [ ] 正常なスコア計算
   - [ ] GLY残基のスキップ
   - [ ] プロファイルファイルが見つからない場合のエラー
   - [ ] Cβ原子がない場合のエラー
   - [ ] 負のプロファイル値の処理

2. **統合テスト（test_profile_integration.py）**
   - [ ] スコア計算なしの統合ワークフロー
   - [ ] スコア計算ありの統合ワークフロー
   - [ ] 結果のスコアソート
   - [ ] 複数パターンでのスコア計算

**完了条件:**
- [ ] 全てのテストがパスする
- [ ] カバレッジが80%以上
- [ ] エッジケースが適切にテストされている

---

### T7: README.mdの更新

**ステータス:** ⏸️ 未着手

**タスク内容:**
- [ ] 新機能の使用例を追加
- [ ] プロファイルスコア計算の説明
- [ ] APIリファレンスの更新

**追加内容:**
- プロファイルスコア計算機能の説明
- 統合ワークフローでの使用例
- 個別機能の使用例

**完了条件:**
- [ ] 使用例が実際に動作する
- [ ] 説明が分かりやすい
- [ ] リンクが正しく設定されている

---

### T8: サンプルスクリプトの作成

**ステータス:** ⏸️ 未着手

**タスク内容:**
- [ ] `examples/profile_scoring_workflow.py`: 完全な使用例

**サンプル内容:**
- 統合ワークフローの実行例
- 個別機能の組み合わせ例
- スコア計算のカスタマイズ例

**完了条件:**
- [ ] サンプルが実際に動作する
- [ ] コメントが充実している
- [ ] 初心者でも理解できる

---

## 🐛 バグと技術的課題

### 発見されたバグ

現在、バグは報告されていません。

### 技術的課題

| 課題 | 影響 | 優先度 | ステータス | 解決策 |
|------|------|--------|-----------|--------|
| プロファイル値が負になる場合の処理 | スコア計算の精度 | 高 | ⏸️ 検討中 | 最小値で置換 |
| 大量のパターンでのパフォーマンス | 処理速度 | 中 | ⏸️ 検討中 | 並列処理の検討 |

---

## 📝 技術的決定の記録

### 決定事項

**決定日:** 2025-10-23

**タイトル:** プロファイルスコア計算の統合方法

**決定内容:**
`integrated_substructure_replacement()`にオプションパラメータとして追加し、後方互換性を維持する。

**理由:**
- 既存のワークフローを壊さない
- 柔軟性が高い（スコア計算の有無を選択可能）
- 段階的な移行が可能

**代替案:**
1. 別関数として実装 → 統合ワークフローの利点が失われる
2. 常にスコア計算を実行 → 不要な場合もコストがかかる

---

## 📋 次のステップ

### 直近のタスク（優先度順）

1. **Phase 1を完了させる**
   - [ ] T1: profile_scoring.pyの基本構造作成
   - [ ] T2: calculate_profile_score()の実装

2. **Phase 2を開始する**
   - [ ] T3: integrated_substructure_replacement()の拡張

3. **簡易テストで動作確認**
   - [ ] 手動テストでスコア計算を確認

### 今後の拡張案

- 並列処理対応（マルチプロセッシング）
- キャッシング機能
- 他のスコアリング関数の統合
- バッチ処理機能

---

## 🔄 変更履歴

### 2025-10-23
- **ドキュメント作成:** 実装進捗記録を初期化
- **Phase 1-5の詳細タスク定義**
- **技術的課題の洗い出し**

---

## 📚 参考リンク

- [`architecture.md`](architecture.md) - システムアーキテクチャ
- [`integration_plan.md`](integration_plan.md) - 統合計画の詳細
- [`testing_guide.md`](testing_guide.md) - テスト手順
- [`task_handoff.md`](task_handoff.md) - タスク引き継ぎ

---

## 💡 開発のヒント

### コーディング規約
- docstringはNumPy形式を使用
- 型ヒントを必ず追加
- エラーメッセージは日本語で分かりやすく

### テストの書き方
- 各関数に対して最低1つのテストケース
- エッジケースを必ずテスト
- テストデータは`tests/data/`に配置

### コミットメッセージ
- 形式: `[Phase X] タスクの説明`
- 例: `[Phase 1] profile_scoring.pyの基本構造作成`