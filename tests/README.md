# inverse_msmd テストスイート

このディレクトリには、`inverse_msmd`パッケージの包括的なテストスイートが含まれています。

## 概要

テストは以下のように構成されています：

- **単体テスト** (`unit/`): 個々の機能やモジュールの動作を検証
- **統合テスト** (`integration/`): 複数のコンポーネントを組み合わせた動作を検証
- **共通フィクスチャ** (`conftest.py`): 全テストで共有されるテストデータとユーティリティ

## テストの実行

### 基本的な使い方

```bash
# プロジェクトルートから全テストを実行
pytest tests/

# 詳細な出力で実行
pytest tests/ -v

# 特定のディレクトリのテストのみ実行
pytest tests/unit/
pytest tests/integration/
```

### テストスクリプトの使用（推奨）

プロジェクトルートにある[`run_tests.sh`](../run_tests.sh)スクリプトを使用すると便利です：

```bash
# 全テスト実行
./run_tests.sh all

# 単体テストのみ
./run_tests.sh unit

# 統合テストのみ
./run_tests.sh integration

# 高速テスト（slowとvisualマーカーを除く）
./run_tests.sh fast

# カバレッジレポート生成
./run_tests.sh coverage
```

## テストマーカー

テストには以下のマーカーが使用されています：

- `@pytest.mark.unit` - 単体テスト（個別の機能をテスト）
- `@pytest.mark.integration` - 統合テスト（複数機能の連携をテスト）
- `@pytest.mark.slow` - 実行時間が長いテスト
- `@pytest.mark.visual` - 視覚的な確認が必要なテスト（画像生成など）

### マーカーを使ったテスト実行

```bash
# 単体テストのみ実行
pytest tests/ -m "unit"

# 統合テストのみ実行
pytest tests/ -m "integration"

# slowとvisualを除いた高速テスト
pytest tests/ -m "not slow and not visual"

# visualテストのみ実行（画像生成を確認したい場合）
pytest tests/ -m "visual"
```

## テスト構造

### 単体テスト (`unit/`)

個別の機能をテストします：

- [`test_imports.py`](unit/test_imports.py) - モジュールインポートのテスト
  - モジュールが正しくインポートできることを確認
  - 必要な依存関係が利用可能か確認
  - 関数スタブの存在確認

- [`test_substructure_search.py`](unit/test_substructure_search.py) - 部分構造探索のテスト
  - 基本的な部分構造探索機能
  - マッチした原子数の検証
  - 原子インデックスの範囲確認
  - 複数プローブでの動作確認

- [`test_visualization.py`](unit/test_visualization.py) - 可視化機能のテスト
  - PNG画像生成の確認
  - 複数マッチの可視化
  - 単一マッチの可視化

- [`test_atom_matching.py`](unit/test_atom_matching.py) - Atom Matching機能のテスト
  - 基本的なマッチング機能
  - atom pairsの形状検証
  - インデックスの範囲確認
  - 元素の一致性確認

### 統合テスト (`integration/`)

複数の機能を組み合わせたワークフローをテストします：

- [`test_workflow.py`](integration/test_workflow.py) - 統合ワークフローのテスト
  - 部分構造探索→atom matching の流れ
  - 変換行列計算までの統合処理
  - 可視化との統合
  - 複数プローブでのワークフロー

## フィクスチャ

[`conftest.py`](conftest.py)には、全テストで使用可能な共通フィクスチャが定義されています：

### パス関連
- `project_root` - プロジェクトルートディレクトリ
- `data_dir` - dataディレクトリ
- `output_dir` - 各テスト用の一時出力ディレクトリ

### テストデータ（分子）
- `ligand_mol` - 4hw3_A_lig.sdfのリガンド分子
- `e23_mol` - E23プローブ分子
- `e24_mol` - E24プローブ分子
- `a01_mol` - A01プローブ分子
- `a08_mol` - A08プローブ分子

### テストデータ（タンパク質）
- `protein_pdb_path` - タンパク質PDBファイルのパス

### 使用例

```python
def test_example(ligand_mol, e23_mol, output_dir):
    """フィクスチャを使用したテスト例"""
    # ligand_mol, e23_mol は自動的に読み込まれる
    # output_dir は各テスト用の一時ディレクトリ
    matches = find_substructure_in_ligand(ligand_mol, e23_mol)
    assert len(matches) > 0
```

## カバレッジ

テストカバレッジを確認するには：

```bash
# カバレッジレポート生成
pytest tests/ --cov=inverse_msmd --cov-report=html --cov-report=term

# HTMLレポートを開く
# htmlcov/index.html をブラウザで開く
```

カバレッジレポートは以下の情報を提供します：
- 各モジュールのカバレッジ率
- テストされていないコード行
- 分岐カバレッジ

## テスト駆動開発（TDD）

新機能を追加する際の推奨ワークフロー：

1. **テストを先に書く**
   ```python
   def test_new_feature(ligand_mol, e23_mol):
       """新機能のテスト"""
       result = new_feature(ligand_mol, e23_mol)
       assert result is not None
   ```

2. **テストを実行して失敗を確認**
   ```bash
   pytest tests/unit/test_new_feature.py -v
   ```

3. **機能を実装**
   ```python
   def new_feature(ligand, substructure):
       # 実装
       return result
   ```

4. **テストを再実行して成功を確認**
   ```bash
   pytest tests/unit/test_new_feature.py -v
   ```

5. **全テストを実行して回帰がないことを確認**
   ```bash
   ./run_tests.sh all
   ```

## CI/CD

継続的インテグレーションでテストを実行する場合：

```bash
# GitHubActionsなどでの使用例
pip install -e ".[dev]"
pytest tests/ --cov=inverse_msmd --cov-report=xml
```

## トラブルシューティング

### テストデータが見つからない

```
pytest.skip: Test data not found: data/...
```

→ `data/`ディレクトリが正しく存在することを確認してください

### ImportError

```
ImportError: No module named 'inverse_msmd'
```

→ パッケージを開発モードでインストールしてください：
```bash
pip install -e .
```

### RDKitエラー

```
ImportError: No module named 'rdkit'
```

→ RDKitをインストールしてください：
```bash
conda install -c conda-forge rdkit
```

## 既存のテストスクリプトとの関係

プロジェクトルートにある既存のテストスクリプト（`test_*.py`）は、開発中の手動テストとして残されています：

- `test_imports.py`
- `test_find_substructure.py`
- `test_visualize_matches.py`
- `test_atom_matching.py`
- `test_match_selection_demo.py`

これらは以下のように使い分けます：

- **`tests/`ディレクトリのテスト**: 自動テスト、CI/CD、品質管理用
- **ルートの`test_*.py`**: 手動実行、開発中の動作確認用

将来的には、ルートのテストスクリプトは削除または`examples/`に移動することを推奨します。

## 参考資料

- [pytest公式ドキュメント](https://docs.pytest.org/)
- [テストチェックリスト](../docs/testing_checklist.md) - 詳細なテスト手順
- [テスト責任分担](../docs/testing_responsibility.md) - AIと人間の役割分担