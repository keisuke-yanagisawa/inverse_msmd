# inverse_msmd

混合溶媒分子動力学法（MSMD）の逆解析ツール

## 概要

`inverse_msmd`は、混合溶媒分子動力学（Mixed-Solvent Molecular Dynamics, MSMD）シミュレーションの逆解析を行うためのPythonパッケージです。タンパク質構造の重ね合わせと、相互作用プロファイルに基づくマッチングスコア計算の機能を提供します。

## 特徴

- **統合アライメントAPI**: 原子マッチングと構造重ね合わせを一括で実行する統合API
- **構造重ね合わせ**: BioPythonベースの構造重ね合わせ（scikit-learn互換インターフェース）
- **原子マッチング**: MCS（最大共通部分構造）ベースの原子ペア検出（アイソトープラベル対応）
- **統合部分構造置換** ✅: リガンドの部分構造置換とタンパク質座標変換を統合
- **PDB操作**: PDBファイルの読み込み、属性の取得・設定、保存機能
- **空間計算**: 球体の体積推定などの空間計算ユーティリティ
- **パス処理**: 環境変数とチルダ展開に対応したパス処理

## インストール

### 必要な環境

- Python >= 3.8
- numpy >= 1.20.0
- biopython >= 1.79
- scikit-learn >= 0.24.0
- scipy >= 1.7.0
- gridData >= 0.6.0
- rdkit（原子マッチング機能を使用する場合）

### インストール方法

#### 開発モード（推奨）

```bash
# リポジトリをクローン
git clone https://github.com/akiyamalab/inverse_msmd.git
cd inverse_msmd

# 開発モードでインストール
pip install -e .
```

#### 通常インストール

```bash
pip install -r requirements.txt
pip install .
```

## 使い方

### 統合アライメントAPI（推奨）

統合APIを使用すると、原子マッチングと構造重ね合わせを一括で実行できます：

```python
from inverse_msmd import align_structures

# プローブ分子とリガンドのマッチングに基づいてタンパク質を重ね合わせ
results = align_structures(
    protein_file="data/sample_proteins/4hw3_A.pdb",
    ligand_file="data/atom_matching/4hw3_A_lig.sdf",
    probe_files={
        "A08": "data/sample_probes/A08.pdb",
        "E24": "data/sample_probes/E24.pdb"
    },
    output_dir="./aligned_structures"
)

# 結果の確認
for result in results:
    print(f"{result.probe_id}_{result.match_id}: "
          f"{result.atom_pairs.shape[1]} 個の原子がマッチ")
```

### 統合部分構造置換API

リガンドの部分構造を別の部分構造で置換し、タンパク質構造も適切に座標変換します：

```python
from inverse_msmd.substructure_replacement import integrated_substructure_replacement

# E23をE24に置換
results = integrated_substructure_replacement(
    ligand_file="data/atom_matching/4hw3_A_lig.sdf",
    protein_file="data/sample_proteins/4hw3_A.pdb",
    from_file="data/sample_probes/E23",  # 拡張子なし
    to_file="data/sample_probes/E24",    # 拡張子なし
    output_dir="output/integrated/",
    match_index=0  # オプション: 特定のマッチを指定
)

# 結果の確認
for i, result in enumerate(results):
    print(f"パターン {i}:")
    print(f"  リガンド: {result['ligand_file']}")

**重要な仕様**:
- 統合部分構造置換では、**置換後の部分構造（to_file）の座標系を基準**とします
- リガンドとタンパク質が置換後の部分構造に合わせて変換されます
- これにより、MSMDプロファイルとの対応関係が維持されます
    print(f"  タンパク質: {result['protein_file']}")
```

**CLIからの使用:**

```bash
python scripts/integrated_replacement.py \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --from-file data/sample_probes/E23 \
    --to-file data/sample_probes/E24 \
    --output output/integrated/ \
    --match-index 0 \
    --verbose
```

**出力ファイル:**
- `pattern_N_ligand_replaced.sdf`: 部分構造が置換されたリガンド
- `pattern_N_protein_aligned.pdb`: 座標変換されたタンパク質
- `substructure_matches.png`: 複数マッチ時の可視化画像（match_index未指定時）

### 低レベルAPI

より細かい制御が必要な場合は、低レベルAPIを使用できます：

```python
from inverse_msmd import SuperImposer, PDB
import numpy as np

# PDBファイルを読み込み
protein = PDB.get_structure("protein.pdb")
probe = PDB.get_structure("probe.pdb")

# 座標を取得
protein_coords = PDB.get_attr(protein, "coord")
probe_coords = PDB.get_attr(probe, "coord")

# 原子ペアに基づいて座標を選択
atom_pairs = np.loadtxt("atom_matching.txt", int)
probe_coords_target = probe_coords[atom_pairs[0]]
protein_coords_target = protein_coords[atom_pairs[1]]

# 構造を重ね合わせ
si = SuperImposer()
si.fit(protein_coords_target, probe_coords_target)

# 変換後の座標を設定
transformed_coords = si.transform(protein_coords)
PDB.set_attr(protein, "coord", transformed_coords)

# 結果を保存
PDB.save(protein, "aligned_protein.pdb")
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
├── pyproject.toml          # パッケージ設定
├── README.md               # このファイル
├── LICENSE                 # ライセンスファイル
├── requirements.txt        # 依存関係リスト
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
│   ├── implementation_progress.md  # 実装進捗記録
│   ├── integrated_replacement_plan.md  # 統合機能設計
│   └── testing_checklist.md  # テスト手順
└── tests/                  # テストスイート
    ├── conftest.py        # 共通フィクスチャ
    ├── unit/              # 単体テスト
    ├── integration/       # 統合テスト
    └── data/              # テストデータ（最小限）
```

## API リファレンス

### アライメントモジュール

#### align_structures()

プローブ分子とリガンドのマッチングに基づいてタンパク質を重ね合わせる統合関数（推奨）

```python
align_structures(protein_file, ligand_file, probe_files, output_dir, iso_value=1)
```

**パラメータ:**
- `protein_file` (str): タンパク質PDBファイルのパス
- `ligand_file` (str): 共結晶化リガンドSDFファイルのパス  
- `probe_files` (str | list | dict): プローブPDBファイル
  - str: 単一プローブのパス
  - list: プローブパスのリスト
  - dict: `{probe_id: filepath}` の辞書
- `output_dir` (str): 出力ディレクトリのパス
- `iso_value` (int): アイソトープラベル値（デフォルト: 1）

**戻り値:** `List[AlignmentResult]` - 重ね合わせ結果のリスト

#### AlignmentResult

重ね合わせ結果を保持するデータクラス

**属性:**
- `probe_id` (str): プローブ分子のID
- `match_id` (int): マッチングID
- `aligned_protein` (Structure): 重ね合わせ後のタンパク質構造
- `atom_pairs` (np.ndarray): 原子ペアのインデックス配列

#### find_atom_matches()

プローブとリファレンス分子間のMCS検索と原子ペアマッチングを実行

```python
find_atom_matches(probe_mol, ref_mol, iso_value=1)
```

**パラメータ:**
- `probe_mol` (Chem.Mol): プローブ分子（RDKit Mol）
- `ref_mol` (Chem.Mol): リファレンス分子
- `iso_value` (int): アイソトープラベル値（デフォルト: 1）

**戻り値:** `List[np.ndarray]` - 原子ペアマッチングのリスト

#### align_structure()

原子ペアに基づいてタンパク質構造を重ね合わせ

```python
align_structure(protein, ligand_coords, probe_coords, atom_pairs)
```

**パラメータ:**
- `protein` (Structure): タンパク質構造（BioPython Structure）
- `ligand_coords` (np.ndarray): リガンド座標
- `probe_coords` (np.ndarray): プローブ座標
- `atom_pairs` (np.ndarray): 原子ペア配列

**戻り値:** `Structure` - 重ね合わせ後の構造

### SuperImposer

構造重ね合わせクラス（scikit-learn互換インターフェース）

**メソッド:**
- `fit(coords, reference_coords)`: 重ね合わせパラメータを計算
- `transform(coords)`: 座標を変換
- `inverse_transform(coords)`: 逆変換を実行

**使用例:**
```python
si = SuperImposer()
si.fit(moving_coords, target_coords)
transformed = si.transform(coords_to_move)
```

### PDB

PDB構造操作のためのユーティリティ関数

**主要関数:**
- `get_structure(filepath)`: PDBファイルを読み込み
- `get_attr(model, attr, sele=None)`: 属性を取得
- `set_attr(model, attr, lst, sele=None)`: 属性を設定
- `save(structs, path)`: 構造をPDBファイルに保存

**使用例:**
```python
from inverse_msmd import PDB

# 読み込み
protein = PDB.get_structure("protein.pdb")

# 座標取得
coords = PDB.get_attr(protein, "coord")

# 座標設定
PDB.set_attr(protein, "coord", new_coords)

# 保存
PDB.save(protein, "output.pdb")
```

### その他のユーティリティ

#### estimate_volume()

球体の集合の体積を推定

```python
estimate_volume(points, radii, granularity=10)
```

**パラメータ:**
- `points` (array-like): 点の座標リスト
- `radii` (array-like): 各球の半径
- `granularity` (int): グリッドの粒度（デフォルト: 10）

**戻り値:** `float` - 推定体積

#### expandpath()

環境変数とチルダを展開したパスを返す

```python
expandpath(path)
```

**パラメータ:**
- `path` (str): 展開するパス文字列

**戻り値:** `str` - 展開されたパス

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

### gridDataが見つからない

```bash
pip install gridData
```

### ImportError: No module named 'inverse_msmd'

パッケージが正しくインストールされているか確認してください：

```bash
pip install -e .
```

### ファイルが見つからない

- サンプルスクリプトは`examples/`ディレクトリから実行してください
- サンプルデータが`data/`ディレクトリに存在することを確認してください
