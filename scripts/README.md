# PDB + SMI → SDF 変換スクリプト

このスクリプトは、PDBファイル（座標情報）とSMIファイル（結合情報）から、両方の情報を含むSDFファイルを作成します。

## 必要な環境

- Python 3.7以上
- RDKit
- inverse_msmdパッケージ

### インストール

```bash
# inverse_msmdパッケージをインストール
pip install -e .

# RDKitをインストール
conda install -c conda-forge rdkit
```

## 使い方

### コマンドラインから使用

```bash
python scripts/create_sdf_from_pdb_smi.py <PDBファイル> <SMIファイル> <出力SDFファイル>
```

#### 例

```bash
# A08プローブのSDFファイルを作成
python scripts/create_sdf_from_pdb_smi.py \
    data/sample_probes/A08.pdb \
    data/sample_probes/A08.smi \
    output/A08.sdf

# E24プローブのSDFファイルを作成
python scripts/create_sdf_from_pdb_smi.py \
    data/sample_probes/E24.pdb \
    data/sample_probes/E24.smi \
    output/E24.sdf

# 分子名を指定する場合
python scripts/create_sdf_from_pdb_smi.py \
    data/sample_probes/A08.pdb \
    data/sample_probes/A08.smi \
    output/A08.sdf \
    --name "Toluene"
```

#### オプション

- `--name`, `-n`: 分子名を指定（省略時は出力ファイル名から自動設定）

#### ヘルプの表示

```bash
python scripts/create_sdf_from_pdb_smi.py --help
```

### Pythonコードから使用

```python
from inverse_msmd import create_sdf_from_pdb_smi

# 基本的な使用法
create_sdf_from_pdb_smi(
    pdb_file="data/sample_probes/A08.pdb",
    smi_file="data/sample_probes/A08.smi",
    output_file="output/A08.sdf"
)

# 分子名を指定
create_sdf_from_pdb_smi(
    pdb_file="molecule.pdb",
    smi_file="molecule.smi",
    output_file="output.sdf",
    molecule_name="MyMolecule"
)

# 進捗メッセージを非表示
create_sdf_from_pdb_smi(
    pdb_file="molecule.pdb",
    smi_file="molecule.smi",
    output_file="output.sdf",
    verbose=False
)
```

## 入力ファイルの要件

### PDBファイル
- 標準的なPDB形式であること
- ATOM または HETATM レコードを含むこと
- 元素記号が76-78列目に記載されているか、atom nameから推定可能であること

### SMIファイル
- 1行目にSMILES文字列が記載されていること
- 形式: `SMILES` または `SMILES NAME`

### 重要な注意点

1. **原子数の一致**: PDBファイルとSMILESから生成される分子の原子数が一致している必要があります
2. **原子の順序**: PDBファイルの原子順序とSMILESから生成される原子順序が一致している必要があります
3. **水素原子**: SMILESに水素が明示的に含まれていない場合、スクリプトが自動的に追加します

## 出力ファイル

SDFファイル（Structure Data File）形式で以下の情報を含みます：

- **分子名**: 指定された名前または出力ファイル名
- **3D座標**: PDBファイルから取得
- **元素記号**: PDBファイルから取得
- **結合情報**: SMILESから取得
- **結合次数**: 単結合、二重結合、芳香族結合など

## トラブルシューティング

### エラー: 原子数が一致しません

```
警告: 原子数が一致しません (SMILES: 15, PDB: 14)
```

**原因**: PDBファイルの原子数とSMILESから生成される分子の原子数が異なります。

**解決方法**:
- PDBファイルに水素原子が含まれているか確認してください
- SMILESが正しいか確認してください
- 必要に応じてPDBファイルまたはSMILESを修正してください

### エラー: 元素の種類または数が一致しません

```
エラー: 元素の種類または数が一致しません
  SMILES由来: ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'H', ...]
  PDB由来:    ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', ...]
```

**原因**: 元素の種類または数が一致していません。

**解決方法**:
- PDBファイルの元素記号が正しいか確認してください
- SMILESが正しい分子を表しているか確認してください

### 警告: 元素の順序が一致しない可能性があります

この警告は、元素の種類と数は一致しているが順序が異なる場合に表示されます。多くの場合、この警告が出ても処理は続行されますが、結果を確認することをお勧めします。

## 動作確認済みの環境

- Python 3.8+
- RDKit 2023.03+
- Linux / macOS / Windows

## サンプルデータ

このリポジトリには以下のサンプルデータが含まれています：

```
data/sample_probes/
├── A08.pdb  # トルエン（メチルベンゼン）のPDBファイル
├── A08.smi  # トルエンのSMILES: Cc1ccccc1
├── E24.pdb  # ビフェニルのPDBファイル
└── E24.smi  # ビフェニルのSMILES: c1ccc(-c2ccccc2)cc1
```

これらのファイルを使用して動作確認ができます。

## ライセンス


---

# アイソトープラベル付与スクリプト

このスクリプトは、SDFファイルの分子に対して、SMARTS記法で指定した部分構造にアイソトープ番号を付与します。

## 使い方

### コマンドラインから使用

```bash
python scripts/add_isotope_labels.py <入力SDF> <出力SDF> <SMARTSパターン> <同位体番号> [オプション]
```

#### 例

```bash
# ベンゼン環の炭素に同位体番号13を付与（全マッチ）
python scripts/add_isotope_labels.py input.sdf output.sdf "c1ccccc1" 13

# 対話的モードで複数マッチから選択
python scripts/add_isotope_labels.py input.sdf output.sdf "c1ccccc1" 13 --interactive

# 2番目のマッチのみに適用
python scripts/add_isotope_labels.py input.sdf output.sdf "c1ccccc1" 13 --match-index 1

# 可視化オプション付き
python scripts/add_isotope_labels.py input.sdf output.sdf "c1ccccc1" 13 --visualize

# カルボニル基の炭素と酸素に異なる同位体番号を付与（2段階）
python scripts/add_isotope_labels.py input.sdf temp.sdf "[C]=[O]" 13
python scripts/add_isotope_labels.py temp.sdf output.sdf "[O]=[C]" 18

# アミン基の窒素に同位体番号15を付与
python scripts/add_isotope_labels.py input.sdf output.sdf "[NH2]" 15
```

#### オプション

- `--interactive`, `-i`: 複数マッチ時に対話的にマッチを選択
- `--match-index N`, `-m N`: N番目のマッチのみに適用（0始まり）
- `--visualize`, `-v`: 同位体ラベル付与前後の分子構造を画像で出力
- `--image-dir DIR`, `-d DIR`: 画像ファイルの出力ディレクトリ

#### ヘルプの表示

```bash
python scripts/add_isotope_labels.py --help
```

## 主な機能

- SMARTS記法による部分構造指定
- 複数マッチの自動検出
- 対話的マッチ選択モード
- ラベル付与前後の分子構造可視化
- 複数分子の一括処理

## 使用する場面

- 特定の官能基のみを原子マッチングの対象にしたい場合
- 複数の部分構造候補から特定のものを選択したい場合
- カスタムマッチング戦略を実装する場合

## トラブルシューティング

### 警告: SMARTSパターンにマッチする部分構造が見つかりませんでした

**原因**: 指定したSMARTSパターンが分子内に存在しません。

**解決方法**:
- SMARTSパターンが正しいか確認してください
- 分子構造を確認してください

### 警告: match_index が範囲外です

**原因**: 指定したマッチインデックスが実際のマッチ数を超えています。

**解決方法**:
- `--interactive`モードで実際のマッチ数を確認してください
- または、`--match-index`を指定せずに実行してマッチ数を確認してください
このスクリプトは、inverse_msmdプロジェクトのライセンスに従います。