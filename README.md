# inverse_msmd

混合溶媒分子動力学法（MSMD）の逆解析ツール

## TL;DR

**リガンド中のプローブ部分構造を入れ替えて、その置換がどれだけ良いかスコアリングするツール**です。

### やること

タンパク質-リガンド複合体に対して、リガンド中のプローブ部分構造を入れ替え、MSMDプロファイルマップと重ね合わせてスコアリングします。

| (A) 複合体構造 | (B) プローブ＋プロファイルマップ | (C) 統合結果 |
|:---:|:---:|:---:|
| <img src="docs/figures/panel_a_complex.png" width="260" alt="タンパク質-リガンド複合体"> | <img src="docs/figures/panel_b_probe_map.png" width="260" alt="プローブ+プロファイルマップ"> | <img src="docs/figures/panel_c_combined.png" width="260" alt="統合結果"> |
| タンパク質＋リガンド | プローブ（緑）＋相互作用マップ | (A)を(B)の座標系に変換して統合 |

### 最小実行例

```bash
# 1件だけ試す（部分構造置換 + スコア計算）
python scripts/integrated_replacement.py \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --from-file data/sample_probes/E23 \
    --to-file data/sample_probes/E24 \
    --output output/integrated/ \
    --profile-dir data/profiles --probe-id E24

# たくさん一括で回す（バッチ処理）
python scripts/run_batch.py \
    --batch-csv examples/batch_config_sample.csv \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --probe-dir data/sample_probes \
    --profile-dir data/profiles \
    --output output/batch_results
```

**入力**: タンパク質PDB + リガンドSDF + プローブ分子ペア（置換前/後）
**出力**: 置換後リガンドSDF + 変換済みタンパク質PDB + スコア

### 同梱プローブ分子

<img src="docs/figures/probe_molecules.png" width="600" alt="プローブ分子一覧: A01, A08, E23, E24">

---

## 概要

`inverse_msmd`は、混合溶媒分子動力学（Mixed-Solvent Molecular Dynamics, MSMD）シミュレーションの逆解析を行うためのPythonパッケージです。タンパク質構造の重ね合わせと、相互作用プロファイルに基づくマッチングスコア計算の機能を提供します。

## 特徴

- **統合部分構造置換**: リガンドの部分構造置換とタンパク質座標変換を統合
- **構造重ね合わせ**: BioPythonベースの構造重ね合わせ
- **原子マッチング**: MCS（最大共通部分構造）ベースの原子ペア検出（アイソトープラベル対応）
- **バッチ処理**: 複数の置換パターンを一括処理し、スコアリング・3D描画まで自動実行
- **3D構造描画**: PyMOLによるタンパク質-リガンド複合体、プローブマップの描画

## インストール

### 必要な環境

- Python >= 3.8
- numpy >= 1.20.0
- biopython >= 1.79
- scikit-learn >= 0.24.0
- scipy >= 1.7.0
- GridDataFormats >= 0.6.0
- libcoffee >= 0.4.3
- matplotlib >= 3.10.1
- rdkit（原子マッチング機能を使用する場合）

### インストール方法

#### 開発モード（推奨）

```bash
# リポジトリをクローン
git clone https://github.com/akiyamalab/inverse_msmd.git
cd inverse_msmd

# 開発モードでインストール（開発用依存関係も含む）
pip install -e .[dev]
```

#### 通常インストール

```bash
# 基本的な依存関係のみ
pip install .

# または、開発用依存関係も含めてインストール
pip install .[dev]
```

## 使い方

### 統合部分構造置換

リガンドの部分構造を別の部分構造で置換し、タンパク質構造も適切に座標変換します。

**重要な仕様**:
- MCS（最大共通部分構造）ベースの正確な重ね合わせ（RMSD < 0.1 Å）
- `ringMatchesRingOnly=True`により、芳香環が正確に重なります
- リガンドとタンパク質の元の座標系を保持します
- **立体障害の自動検出**: 置換後に原子間距離が2.0Å未満の構造を自動的に除外

```bash
# 基本的な使用方法
python scripts/integrated_replacement.py \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --from-file data/sample_probes/E23 \
    --to-file data/sample_probes/E24 \
    --output output/integrated/

# 詳細な進捗表示と特定のマッチを指定
python scripts/integrated_replacement.py \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --from-file data/sample_probes/E23 \
    --to-file data/sample_probes/E24 \
    --output output/integrated/ \
    --match-index 0 \
    --verbose
```

**コマンドラインオプション:**

| オプション | 必須 | 説明 |
|-----------|------|------|
| `--ligand` | ✓ | リガンドSDFファイルのパス |
| `--protein` | ✓ | タンパク質PDBファイルのパス |
| `--from-file` | ✓ | 置換前の部分構造（拡張子なし、.pdbと.smiを自動読込） |
| `--to-file` | ✓ | 置換後の部分構造（拡張子なし、.pdbと.smiを自動読込） |
| `--output` | ✓ | 出力ディレクトリのパス |
| `--match-index` | | 特定のマッチパターンを指定（0始まり） |
| `--verbose` | | 詳細な進捗情報を表示 |
| `--version` | | バージョン情報を表示 |

**出力ファイル:**
- `pattern_N_ligand_replaced.sdf`: 部分構造が置換されたリガンド
- `pattern_N_protein_aligned.pdb`: 座標変換されたタンパク質
- `substructure_matches.png`: 複数マッチ時の可視化画像（match_index未指定時）

**注意**: 立体障害チェックにより、化学的に不適切な構造は自動的に除外されます。24パターン生成される場合、通常16パターン程度の有効な構造が出力されます。

**実行結果の確認:**
```bash
# 生成されたファイル数を確認
ls output/integrated/*.sdf | wc -l

# ファイル一覧を表示
ls -lh output/integrated/
```

### サンプルスクリプト

パッケージには、サンプルデータ（[`data/`](data/)ディレクトリ）付きのサンプルスクリプト（[`examples/`](examples/)ディレクトリ）とツール（[`scripts/`](scripts/)ディレクトリ）が含まれています：

**サンプルスクリプト:**
- [`examples/integrated_alignment.py`](examples/integrated_alignment.py): ⭐ 統合アライメントAPIの使用例（推奨）
- [`examples/calculate_matching.py`](examples/calculate_matching.py): マッチングスコア計算の例
- [`examples/add_isotope_labels.py`](examples/add_isotope_labels.py): アイソトープラベル付与ツール

**コマンドラインツール:**
- [`scripts/integrated_replacement.py`](scripts/integrated_replacement.py): ⭐ 統合部分構造置換ツール（新機能）
- [`scripts/replace_substructure.py`](scripts/replace_substructure.py): 部分構造置換ツール

詳細な使用方法については、各スクリプト内のコメントと[`examples/README.md`](examples/README.md)、[`scripts/README.md`](scripts/README.md)を参照してください。

## ディレクトリ構造

```
inverse_msmd/
├── pyproject.toml          # パッケージ設定（依存関係を含む）
├── README.md               # このファイル
├── LICENSE                 # ライセンスファイル
├── inverse_msmd/           # メインパッケージ
│   ├── __init__.py        # パブリックAPI
│   ├── alignment.py       # アライメント機能
│   ├── substructure_replacement.py  # 統合部分構造置換
│   └── utils/             # ユーティリティモジュール
│       ├── __init__.py
│       ├── bio_utils.py   # BioPythonユーティリティ
│       ├── mol_utils.py   # RDKitユーティリティ
│       ├── spatial_utils.py  # 空間計算ユーティリティ
│       └── path_utils.py  # パス処理ユーティリティ
├── data/                   # サンプルデータ
│   ├── sample_proteins/   # タンパク質構造ファイル（PDB）
│   ├── sample_probes/     # プローブ分子ファイル（PDB）
│   ├── atom_matching/     # 原子マッチングデータ
│   └── profiles/          # 相互作用プロファイル（.dx.gz）
├── examples/               # サンプルスクリプト
│   ├── README.md          # サンプルの説明
│   ├── integrated_alignment.py  # 統合API使用例
│   ├── calculate_matching.py  # スコア計算例
│   └── add_isotope_labels.py  # アイソトープラベル付与
├── scripts/                # コマンドラインツール
│   ├── README.md          # ツールの説明
│   ├── integrated_replacement.py  # 統合部分構造置換CLI
│   ├── replace_substructure.py  # 部分構造置換ツール
│   └── add_isotope_labels.py  # アイソトープラベル付与
├── docs/                   # ドキュメント
│   └── figures/            # README用の図
└── tests/                  # テストスイート
    ├── conftest.py        # 共通フィクスチャ
    ├── unit/              # 単体テスト
    ├── integration/       # 統合テスト
    └── data/              # テストデータ（最小限）
```

## ライセンス

MIT License

## 著者

Keisuke Yanagisawa (yanagisawa@comp.isct.ac.jp)

## テスト

このプロジェクトには包括的なテストスイートが含まれており、品質管理を確実に行うことができます。

### テストの実行

#### 簡単な方法（推奨）

提供されているテスト実行スクリプトを使用します：

```bash
# 全てのテストを実行
./run_tests.sh all

# 単体テストのみ実行
./run_tests.sh unit

# 統合テストのみ実行
./run_tests.sh integration

# 高速テスト（visualとslowマーカーを除く）
./run_tests.sh fast

# カバレッジレポート付きで実行
./run_tests.sh coverage
```

#### pytestを直接使用

```bash
# 全てのテストを実行
pytest tests/

# 詳細な出力で実行
pytest tests/ -v

# 特定のマーカーのみ実行
pytest tests/ -m "unit"
pytest tests/ -m "integration"

# 特定のテストファイルを実行
pytest tests/unit/test_imports.py

# カバレッジレポートを生成
pytest tests/ --cov=inverse_msmd --cov-report=html
```

### テストの構造

```
tests/
├── conftest.py              # 共通フィクスチャ定義
├── unit/                    # 単体テスト
│   ├── test_imports.py     # インポートテスト
│   ├── test_substructure_search.py  # 部分構造探索
│   ├── test_visualization.py        # 可視化機能
│   └── test_atom_matching.py        # Atom Matching
└── integration/             # 統合テスト
    └── test_workflow.py    # ワークフロー全体のテスト
```

### テストマーカー

- `@pytest.mark.unit` - 単体テスト
- `@pytest.mark.integration` - 統合テスト
- `@pytest.mark.slow` - 実行時間が長いテスト
- `@pytest.mark.visual` - 視覚的な確認が必要なテスト

### 開発時のテスト

開発時は以下のワークフローを推奨します：

```bash
# 1. 高速テストを実行して基本的な動作を確認
./run_tests.sh fast

# 2. 変更した機能に関連するテストを実行
pytest tests/unit/test_substructure_search.py -v

# 3. 全テストを実行して問題がないことを確認
./run_tests.sh all

# 4. カバレッジを確認
./run_tests.sh coverage
# htmlcov/index.html をブラウザで開く
```

## 引用

このパッケージを使用する場合は、適切な引用を行ってください。

## 貢献

バグ報告や機能要望は、GitHub Issuesでお願いします。

## トラブルシューティング

### RDKitのインストールエラー

RDKitは`pip`ではなく`conda`を使用することを推奨します：

```bash
conda install -c conda-forge rdkit
```

### ImportError: No module named 'inverse_msmd'

パッケージが正しくインストールされているか確認してください：

```bash
pip install -e .
# または開発用依存関係も含めて
pip install -e .[dev]
```

### ファイルが見つからない

- サンプルスクリプトは`examples/`ディレクトリから実行してください
- サンプルデータが`data/`ディレクトリに存在することを確認してください
