# プロファイルスコア計算機能 タスク引き継ぎガイド

最終更新: 2025-10-23

## 🎯 このドキュメントの目的

このドキュメントは、AI間でのタスク引き継ぎを円滑にするために作成されています。前回のセッションから継続して作業を行う際に、このドキュメントを読むことで迅速に状況を把握できます。

---

## 📊 現在の状況

### プロジェクトステータス

**全体進捗:** 0% (未着手)

**現在のフェーズ:** Phase 0 - ドキュメント整備完了

**最終更新日:** 2025-10-23

### 完了した作業

✅ **ドキュメント体系の構築**
- [`README.md`](README.md) - ドキュメントインデックス
- [`architecture.md`](architecture.md) - システムアーキテクチャ設計
- [`implementation_progress.md`](implementation_progress.md) - 実装進捗記録
- [`testing_guide.md`](testing_guide.md) - テスト手順書
- [`task_handoff.md`](task_handoff.md) - このドキュメント

---

## 🚀 次に取り組むべきタスク

### 優先度1: Phase 1の開始（プロファイルスコア計算モジュール）

**タスク番号:** T1

**タスク名:** `inverse_msmd/profile_scoring.py`の基本構造作成

**所要時間（推定）:** 30分

**詳細:** [`implementation_progress.md`](implementation_progress.md#t1-inverse_msmdprofile_scoringpyの基本構造作成) を参照

**実装手順:**
1. `inverse_msmd/profile_scoring.py`ファイルを作成
2. モジュールdocstringを追加
3. 必要なインポートを追加
4. 基本的なコメント構造を作成

**成功基準:**
- [ ] ファイルが作成されている
- [ ] インポートエラーがない
- [ ] `from inverse_msmd.profile_scoring import calculate_profile_score`でエラーが出ない（関数はまだ未実装でOK）

**次のステップ:** T2に進む

---

### 優先度2: T2の実装

**タスク番号:** T2

**タスク名:** `calculate_profile_score()`関数の実装

**所要時間（推定）:** 2-3時間

**詳細:** [`implementation_progress.md`](implementation_progress.md#t2-calculate_profile_score関数の実装) を参照

**実装のポイント:**
- プロファイル読み込みロジック
- Cβ原子抽出ロジック
- 3D補間とスコア計算
- エラーハンドリング

**成功基準:**
- [ ] 関数が正しく動作する
- [ ] 簡易テストがパスする
- [ ] エラーケースが適切に処理される

---

## 📚 重要なコンテキスト

### プロジェクトの背景

このプロジェクトは、`integrated_substructure_replacement`機能を拡張し、プロファイルマッチングスコア計算を統合するものです。

**目的:**
- リガンド部分構造の置換と同時にスコア計算を実行
- 複数パターンをスコアでランク付け
- 後方互換性を維持

### 設計上の重要な決定事項

1. **後方互換性の維持**
   - `profile_dir`パラメータをオプションにする
   - `None`の場合はスコア計算をスキップ

2. **モジュラー設計**
   - スコア計算機能を独立モジュールとして実装
   - 個別機能をパッケージレベルでエクスポート

3. **エラーハンドリング**
   - プロファイルファイルが見つからない場合の適切なエラー
   - GLY残基（Cβなし）の自動スキップ

詳細は [`architecture.md`](architecture.md) を参照してください。

---

## 🔍 知っておくべき技術的詳細

### プロファイルスコア計算のアルゴリズム

```python
# 基本的な流れ
1. プロファイルファイルを読み込み（19種類の残基タイプ）
2. Cβ原子を抽出（GLYを除く）
3. 各Cβ位置でプロファイル値を3D補間
4. 距離重み付けを適用（オプション）
5. 対数スケールで統合
```

### 使用するライブラリ

- **gridData**: プロファイル読み込みと3D補間
- **Bio.PDB**: タンパク質構造の処理
- **NumPy**: 数値計算
- **RDKit**: 分子構造の処理

### ファイル構成

```
inverse_msmd/
├── profile_scoring.py         # 新規作成予定
├── substructure_replacement.py # 更新予定
└── __init__.py                # 更新予定
```

---

## ⚠️ 注意事項

### 実装時の注意点

1. **プロファイル値の処理**
   - 負の値が出た場合は`profiles[resname].grid.min()`で置換
   - NaNやInfが発生しないよう注意

2. **座標変換**
   - プローブ中心座標は変換後の座標を使用
   - `transformed_center = np.dot(to_center, rot) + tran`

3. **パフォーマンス**
   - プロファイルは各残基タイプにつき1回のみ読み込み
   - 不要な計算を避ける

### よくある落とし穴

❌ **間違い:** 全ての原子でスコア計算
✅ **正解:** Cβ原子のみでスコア計算

❌ **間違い:** プロファイルを毎回読み込み
✅ **正解:** 最初に一度だけ読み込んで辞書に保存

❌ **間違い:** 変換前の座標でスコア計算
✅ **正解:** 変換後の座標でスコア計算

---

## 🧪 テスト戦略

### 実装時のテスト手順

1. **最小限の動作確認**
   ```python
   # T1完了後
   from inverse_msmd.profile_scoring import calculate_profile_score
   print("インポート成功")
   ```

2. **簡易機能テスト**
   ```python
   # T2完了後
   from inverse_msmd.profile_scoring import calculate_profile_score
   from inverse_msmd.utils.bio_utils import PDB
   import numpy as np
   
   protein = PDB.get_structure("data/sample_proteins/4hw3_A.pdb")
   center = np.array([10.0, 15.0, 20.0])
   score = calculate_profile_score(protein, center, "data/profiles/", "E24")
   print(f"Score: {score}")
   ```

3. **統合テスト**
   - Phase 2完了後に実施
   - 詳細は [`testing_guide.md`](testing_guide.md) を参照

---

## 📋 チェックリスト：作業開始前

新しいAIセッションで作業を開始する前に、以下を確認してください：

- [ ] [`README.md`](README.md) を読んでプロジェクト概要を把握
- [ ] [`architecture.md`](architecture.md) でシステム設計を理解
- [ ] [`implementation_progress.md`](implementation_progress.md) で現在の進捗を確認
- [ ] このドキュメントで次のタスクを確認
- [ ] 必要なデータファイルが存在するか確認
  - [ ] `data/profiles/E24_*_profile.dx.gz`
  - [ ] `data/sample_proteins/4hw3_A.pdb`
  - [ ] `data/sample_probes/E24.*`

---

## 🔄 作業完了後の手順

タスクを完了したら、以下を実行してください：

1. **進捗記録の更新**
   - [`implementation_progress.md`](implementation_progress.md) のステータスを更新
   - 完了したタスクにチェックマークを付ける
   - 技術的な決定事項や課題があれば記録

2. **このドキュメントの更新**
   - 「現在の状況」セクションを更新
   - 「次に取り組むべきタスク」を更新
   - 新しく発見した注意事項があれば追加

3. **テストの実行**
   - 実装した機能の簡易テストを実行
   - 結果を記録

4. **コミット**
   - 変更をgitにコミット
   - コミットメッセージの形式: `[Phase X] タスクの説明`

---

## 💬 引き継ぎメッセージ（テンプレート）

```markdown
## 作業完了報告

**実施日:** YYYY-MM-DD
**完了したタスク:** TX - タスク名
**所要時間:** X時間
**ステータス:** ✅ 完了 / ⚠️ 部分完了 / ❌ 問題あり

### 実施内容
- 箇条書きで実施した内容を記載

### 技術的な決定事項
- 重要な設計判断や実装方法の選択理由

### 発見した問題
- バグや課題があれば記載

### 次のセッションへの引き継ぎ事項
- 次に取り組むべきタスク
- 注意すべき点
- 未解決の課題

### テスト結果
- [ ] インポートテスト: 成功/失敗
- [ ] 簡易機能テスト: 成功/失敗
- [ ] 特記事項: あれば記載
```

---

## 🔗 クイックリンク

### ドキュメント
- [README](README.md) - ドキュメントインデックス
- [Architecture](architecture.md) - システム設計
- [Progress](implementation_progress.md) - 進捗記録
- [Testing](testing_guide.md) - テスト手順
- [Integration Plan](integration_plan.md) - 統合計画の詳細

### コードファイル
- [`inverse_msmd/__init__.py`](../../inverse_msmd/__init__.py) - パッケージエントリポイント
- [`inverse_msmd/substructure_replacement.py`](../../inverse_msmd/substructure_replacement.py) - 統合ワークフロー
- [`scripts/integrated_replacement.py`](../../scripts/integrated_replacement.py) - CLIスクリプト

### データ
- [`data/profiles/`](../../data/profiles/) - プロファイルファイル
- [`data/sample_proteins/`](../../data/sample_proteins/) - サンプルタンパク質
- [`data/sample_probes/`](../../data/sample_probes/) - サンプルプローブ

### テスト
- [`tests/unit/`](../../tests/unit/) - ユニットテスト
- [`tests/integration/`](../../tests/integration/) - 統合テスト

---

## 📞 サポート情報

### 困ったときは

1. **エラーが発生した場合**
   - エラーメッセージを確認
   - [`testing_guide.md`](testing_guide.md) のトラブルシューティングを参照
   - [`implementation_progress.md`](implementation_progress.md) の技術的課題を確認

2. **設計判断に迷った場合**
   - [`architecture.md`](architecture.md) の設計原則を確認
   - [`integration_plan.md`](integration_plan.md) の実装詳細を参照
   - 不明な点は質問として記録

3. **テストに失敗した場合**
   - [`testing_guide.md`](testing_guide.md) で期待される動作を確認
   - 簡易テストから段階的に確認

---

## 📊 プロジェクト全体のマイルストーン

```
Phase 0: ドキュメント整備 ✅ 完了
  └─ ドキュメント体系の構築

Phase 1: プロファイルスコア計算モジュール ⏸️ 未着手
  ├─ T1: 基本構造作成
  └─ T2: スコア計算関数の実装

Phase 2: 統合ワークフローの拡張 ⏸️ 未着手
  └─ T3: integrated_substructure_replacement()の拡張

Phase 3: パッケージエクスポート ⏸️ 未着手
  └─ T4: __init__.pyの更新

Phase 4: CLIスクリプトの拡張 ⏸️ 未着手
  └─ T5: integrated_replacement.pyの更新

Phase 5: テストとドキュメント ⏸️ 未着手
  ├─ T6: テストコードの作成
  ├─ T7: READMEの更新
  └─ T8: サンプルスクリプトの作成
```

---

## 🎓 学習リソース

### 関連技術の理解

- **gridData**: [MDAnalysis GridDataFormats](https://www.mdanalysis.org/GridDataFormats/)
- **Bio.PDB**: [Biopython Tutorial](https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ)
- **RDKit**: [RDKit Documentation](https://www.rdkit.org/docs/)

### プロジェクト固有の知識

- [`examples/calculate_matching.py`](../../examples/calculate_matching.py) - 元のスコア計算実装
- [`../../README.md`](../../README.md) - プロジェクト全体の説明

---

**最終更新:** 2025-10-23
**次回レビュー予定:** Phase 1完了時