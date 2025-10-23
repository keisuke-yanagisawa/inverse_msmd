# プロファイルスコア計算機能 ドキュメント

このディレクトリは、プロファイルマッチングスコア計算機能の統合に関するドキュメントを含んでいます。

## 📚 ドキュメント一覧

### 必須ドキュメント

- **[`integration_plan.md`](integration_plan.md)** ⭐⭐⭐
  - プロファイルスコア計算機能の統合計画
  - アーキテクチャ設計と実装仕様
  - タスクリストと処理フロー

### 開発中のドキュメント

- **[`architecture.md`](architecture.md)** ⭐⭐
  - システムアーキテクチャの詳細
  - モジュール間の依存関係
  - データフロー図

- **[`implementation_progress.md`](implementation_progress.md)** ⭐⭐
  - 実装進捗の詳細記録
  - 完了したタスクと残タスク
  - バグ修正と技術的決定の記録

- **[`testing_guide.md`](testing_guide.md)** ⭐⭐
  - テスト手順とチェックリスト
  - 手動確認が必要な項目
  - テストデータの説明

- **[`task_handoff.md`](task_handoff.md)** ⭐
  - AI間のタスク引き継ぎガイド
  - 現在の作業状況
  - 次のステップと注意事項

## 🎯 クイックスタート

### このドキュメント群の目的

プロファイルスコア計算機能の統合は、複数のAIセッションにわたって開発される可能性があります。これらのドキュメントは以下を目的としています：

1. **一貫性の維持**: 複数のAIセッション間で設計思想を共有
2. **効率的な引き継ぎ**: 前回のセッションからスムーズに再開
3. **進捗の可視化**: 何が完了し、何が残っているかを明確化
4. **技術的決定の記録**: なぜその実装方法を選んだかの記録

### ドキュメントの読み方

**初めて参加する場合:**
1. [`integration_plan.md`](integration_plan.md) で全体像を把握
2. [`architecture.md`](architecture.md) でアーキテクチャを理解
3. [`implementation_progress.md`](implementation_progress.md) で現在の進捗を確認
4. [`task_handoff.md`](task_handoff.md) で次のアクションを確認

**作業を再開する場合:**
1. [`task_handoff.md`](task_handoff.md) で前回の作業を確認
2. [`implementation_progress.md`](implementation_progress.md) で進捗を確認
3. 必要に応じて [`testing_guide.md`](testing_guide.md) でテスト方法を確認

## 📋 プロジェクト概要

### 機能概要

`integrated_substructure_replacement`機能を拡張し、プロファイルマッチングスコア計算を統合します。

**主な機能:**
- リガンド部分構造の置換
- タンパク質構造の座標変換
- **プロファイルスコア計算**（新規）
- 結果のスコアソート（新規）

### 開発フェーズ

```
Phase 1: プロファイルスコア計算モジュール作成
Phase 2: 統合ワークフローの拡張
Phase 3: パッケージエクスポートの更新
Phase 4: CLIスクリプトの拡張
Phase 5: テストとドキュメント
```

## 🔄 ドキュメントのライフサイクル

### 開発完了後の整理計画

このプロジェクトが完了したら、以下のドキュメント整理を実施します：

**統合先:**
- 重要な設計思想 → [`../../README.md`](../../README.md)に統合
- API使用例 → [`../../README.md`](../../README.md)と[`../../examples/README.md`](../../examples/README.md)
- テスト方法 → [`../../tests/README.md`](../../tests/README.md)
- 実装詳細 → コードのdocstring

**削除対象:**
- [`implementation_progress.md`](implementation_progress.md) - 進捗記録は開発完了後は不要
- [`task_handoff.md`](task_handoff.md) - タスク管理は開発完了後は不要
- [`testing_guide.md`](testing_guide.md) - テストスイートに移行後は不要

**保持対象:**
- [`architecture.md`](architecture.md) - 設計思想は長期的に有用（簡潔化）
- [`integration_plan.md`](integration_plan.md) - 技術的な詳細リファレンスとして保持（簡潔化）

詳細は [`../documentation_best_practices.md`](../documentation_best_practices.md) を参照してください。

## 🔗 関連ドキュメント

- [`../../README.md`](../../README.md) - プロジェクトメインドキュメント
- [`../../examples/README.md`](../../examples/README.md) - 使用例
- [`../../tests/README.md`](../../tests/README.md) - テスト方法
- [`../documentation_best_practices.md`](../documentation_best_practices.md) - ドキュメントベストプラクティス

## 📝 更新履歴

| 日付 | 更新内容 |
|------|---------|
| 2025-10-23 | 初版作成 - ドキュメント体系の整備 |