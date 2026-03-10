# テスト

## テストの実行

### run_tests.sh（推奨）

```bash
./run_tests.sh all         # 全テスト
./run_tests.sh unit        # 単体テストのみ
./run_tests.sh integration # 統合テストのみ
./run_tests.sh fast        # 高速テスト（visual/slowを除く）
./run_tests.sh coverage    # カバレッジレポート付き
```

### pytest を直接使用

```bash
pytest tests/                                     # 全テスト
pytest tests/ -v                                  # 詳細出力
pytest tests/ -m "unit"                           # 単体テストのみ
pytest tests/ -m "integration"                    # 統合テストのみ
pytest tests/unit/test_imports.py                 # 特定ファイル
pytest tests/ --cov=inverse_msmd --cov-report=html  # カバレッジ
```

## テスト構造

```
tests/
├── conftest.py              # 共通フィクスチャ
├── unit/                    # 単体テスト
│   ├── test_imports.py
│   ├── test_substructure_search.py
│   ├── test_visualization.py
│   └── test_atom_matching.py
└── integration/             # 統合テスト
    └── test_workflow.py
```

## テストマーカー

| マーカー | 説明 |
|---------|------|
| `@pytest.mark.unit` | 単体テスト |
| `@pytest.mark.integration` | 統合テスト |
| `@pytest.mark.slow` | 実行時間が長いテスト |
| `@pytest.mark.visual` | 視覚的な確認が必要なテスト |

## 開発時のワークフロー

```bash
# 1. 高速テストで基本確認
./run_tests.sh fast

# 2. 変更した機能のテスト
pytest tests/unit/test_substructure_search.py -v

# 3. 全テストで最終確認
./run_tests.sh all
```
