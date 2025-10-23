# inverse_msmd ドキュメント

このディレクトリには、inverse_msmdパッケージの各機能に関する詳細なドキュメントが含まれています。

## 📚 ドキュメント一覧

### プロファイルスコア統合機能
- **[profile_scoring_integration_plan.md](profile_scoring_integration_plan.md)** - プロファイルマッチングスコア計算を統合する機能の詳細な実装計画

## 🎯 機能概要

### プロファイルスコア統合
`integrated_substructure_replacement`機能を拡張し、部分構造置換の各パターンをプロファイルマッチングスコアで評価できるようにします。

**主な特徴:**
- プロファイルベースのスコア計算
- 個別機能のパッケージエクスポート
- 統合ワークフローへのスコア計算組み込み
- スコアによる結果の自動ソート
- 完全な後方互換性

**使用例:**
```python
from inverse_msmd.substructure_replacement import integrated_substructure_replacement

results = integrated_substructure_replacement(
    ligand_file="ligand.sdf",
    protein_file="protein.pdb",
    from_file="data/sample_probes/E23",
    to_file="data/sample_probes/E24",
    output_dir="output/",
    profile_dir="data/profiles/",  # プロファイルスコア計算
    probe_id="E24"
)

# 結果はスコアで降順ソート済み
for i, result in enumerate(results):
    print(f"Pattern {i}: score={result['score']:.2f}")
```

## 📖 ドキュメント構成

各ドキュメントには以下の情報が含まれています：

1. **概要と要求仕様** - 機能の目的と入出力
2. **アーキテクチャ設計** - モジュール構成と主要関数
3. **処理フロー** - 詳細なワークフロー図
4. **実装タスクリスト** - 段階的な実装計画
5. **使用例** - 実践的なコード例
6. **テスト計画** - 品質保証のための計画
7. **技術的考慮事項** - パフォーマンスやエラーハンドリング

## 🚀 実装の進め方

各機能の実装は、ドキュメント内の実装タスクリストに従って進めます：

1. Phase 1: 基本モジュールの作成
2. Phase 2: コア機能の実装
3. Phase 3: パッケージ統合
4. Phase 4: CLI拡張
5. Phase 5: テストとドキュメント

## 🔗 関連リンク

- [README.md](../README.md) - パッケージの概要
- [examples/](../examples/) - 使用例
- [tests/](../tests/) - テストコード