# サンプルスクリプト

このディレクトリには、`inverse_msmd`パッケージの使用方法を示すサンプルスクリプトが含まれています。

## 前提条件

サンプルスクリプトを実行する前に、以下を確認してください：

1. パッケージがインストールされていること：
```bash
pip install -e .
```

2. サンプルデータが[`../data/`](../data/)ディレクトリに配置されていること

3. RDKitがインストールされていること（原子マッチング機能を使用する場合）：
```bash
conda install -c conda-forge rdkit
```

## サンプルデータ構造

```
data/
├── sample_proteins/     # タンパク質構造ファイル（PDB形式）
│   ├── 4hw3_A.pdb      # サンプルタンパク質（PDB ID: 4hw3, Chain A）
│   └── 4hw3.pdb        # フルタンパク質
├── sample_probes/       # プローブ分子ファイル（PDB形式）
│   ├── A08.pdb         # アセトン（プローブA08）
│   └── E24.pdb         # エタノール（プローブE24）
├── atom_matching/       # 原子マッチングファイル
│   ├── 4hw3_A_lig.sdf  # 共結晶化リガンド（SDF形式）
│   ├── 4hw3_A_lig_with_subst_label.sdf  # アイソトープラベル付きリガンド
│   └── atom_matching_* # 原子マッチング結果ファイル
└── profiles/            # 相互作用プロファイル（gzip圧縮済み）
    ├── A08_ALA_profile.dx.gz  # プローブA08のアラニンプロファイル
    ├── E24_ALA_profile.dx.gz  # プローブE24のアラニンプロファイル
    └── ...                     # 他のアミノ酸プロファイル
```

**注意**: プロファイルファイルは効率的なストレージのためgzip圧縮（.dx.gz）されています。gridDataライブラリは圧縮ファイルを直接読み込めるため、解凍する必要はありません。

## サンプルスクリプト一覧

### 1. 統合アライメント（integrated_alignment.py）⭐ 推奨

**説明:**  
統合APIを使用して、原子マッチングと構造重ね合わせを一括で実行します。最もシンプルで使いやすいAPIです。

**実行方法:**
```bash
cd examples
python integrated_alignment.py
```

**処理内容:**
1. プローブ分子（A08, E24）と共結晶化リガンドをロード
2. MCS（最大共通部分構造）を自動検索
3. 原子ペアを自動特定
4. タンパク質構造を各プローブに重ね合わせ
5. 結果を`./aligned_structures/`に自動保存

**主な機能:**
- ワンステップでアライメント完了
- 複数プローブの一括処理
- アイソトープラベル対応
- 結果の自動保存

**出力:**
```
aligned_structures/
├── aligned_to_A08_0.pdb
├── aligned_to_A08_1.pdb
├── aligned_to_E24_0.pdb
└── ...
```

**コード例:**
```python
from inverse_msmd import align_structures

results = align_structures(
    protein_file="../data/sample_proteins/4hw3_A.pdb",
    ligand_file="../data/atom_matching/4hw3_A_lig_with_subst_label.sdf",
    probe_files={
        "A08": "../data/sample_probes/A08.pdb",
        "E24": "../data/sample_probes/E24.pdb",
    },
    output_dir="./aligned_structures",
    iso_value=1  # ISO=1ラベルのある原子のみを使用
)
```

---

### 2. マッチングスコア計算（calculate_matching.py）

**説明:**  
相互作用プロファイルに基づいてマッチングスコアを計算するサンプルです。重ね合わせ済みのタンパク質構造とプローブのプロファイルを使用します。

**実行方法:**
```bash
cd examples
python calculate_matching.py
```

**処理内容:**
1. 重ね合わせ済みのタンパク質構造を読み込み
2. プローブ分子の相互作用プロファイルを読み込み
3. 各残基のCβ原子位置でプロファイル値を補間
4. プローブ中心からの距離で重み付け
5. 対数スケールでマッチングスコアを計算

**主な機能:**
- gridDataライブラリを使用したプロファイルの読み込み
- 3D空間での線形補間を使用したスコア計算
- 距離ベースの重み付け（ガウシアン）
- 複数残基タイプの統合評価

**パラメータ:**
- `GAMMA`: 距離重み付けパラメータ（デフォルト: 0.00）
  - 0.00: 重み付けなし（全残基均等）
  - 0.003: 距離に基づく重み付け

**出力例:**
```
A08 -245.67
E24 -189.23
E24_1 -201.45
E24_2 -195.89
```

**使用する場面:**
- MSMD解析結果の評価
- タンパク質-リガンド相互作用の定量化
- 複数プローブ配置の比較

---

### 3. アイソトープラベル付与ツール

**注意:** このツールは[`../scripts/add_isotope_labels.py`](../scripts/add_isotope_labels.py)に移動しました。

**説明:**
SDFファイルの分子に対して、SMARTS記法で指定した部分構造にアイソトープ番号を付与する実用的なCLIツールです。原子マッチングで特定の部分構造のみを考慮したい場合に使用します。

**詳細なドキュメント:**
[`../scripts/README.md`](../scripts/README.md)の「アイソトープラベル付与スクリプト」セクションを参照してください。

**クイックスタート:**
```bash
# 基本的な使用法
python scripts/add_isotope_labels.py input.sdf output.sdf "c1ccccc1" 13

# 対話的モードで選択
python scripts/add_isotope_labels.py input.sdf output.sdf "c1ccccc1" 13 --interactive

# ヘルプを表示
python scripts/add_isotope_labels.py --help
```

**使用する場面:**
- 特定の官能基のみを原子マッチングの対象にしたい場合
- 複数の部分構造候補から特定のものを選択したい場合
- カスタムマッチング戦略を実装する場合

---

## 実行順序

### 基本的なワークフロー（統合API使用）

```bash
# 1. 統合アライメントを実行（推奨）
python integrated_alignment.py

# 2. マッチングスコアを計算
python calculate_matching.py
```

## カスタマイズ

### 使用するマッチングIDの変更

各スクリプト内のマッチングIDリストを編集：

```python
# superimposition.py
matching_ids = [f"A08_{i}" for i in range(12)] + [f"E24_{i}" for i in range(24)]

# calculate_matching.py
matching_ids = "A08 E24 E24_1 E24_2".split(" ")
```

### パラメータの調整

距離重み付けパラメータを調整（calculate_matching.py）：

```python
GAMMA = 0.003  # 距離重み付けパラメータ
# 0.00: 重み付けなし
# 0.003: 距離に基づく重み付け
```

### プローブの追加

新しいプローブを追加する場合：

1. プローブPDBファイルを`../data/sample_probes/`に配置
2. プローブのプロファイルを`../data/profiles/`に配置
3. スクリプトのプローブリストに追加

## 注意事項

### 実行ディレクトリ

- **必ず`examples/`ディレクトリから実行してください**
- 相対パスが`examples/`ディレクトリを基準に設定されています

### データファイル

- サンプルデータがない場合、スクリプトの実行に失敗します
- プロファイルファイル（.dx形式）は、すべての残基タイプに対して必要です
- アイソトープラベル付きリガンドがない場合は、先に作成してください

### プロファイルファイル

- プロファイルファイル（.dx.gz）はgzip圧縮されています
- 約83%のサイズ削減（非圧縮時: 約339MB → 圧縮後: 約57MB）
- gridDataライブラリが.dx.gzファイルを直接読み込めます
- 独自のプロファイルを追加する場合：
  ```bash
  gzip -9 your_profile.dx
  ```

## トラブルシューティング

### ImportError が発生する場合

パッケージが正しくインストールされているか確認：
```bash
pip install -e .
```

### FileNotFoundError が発生する場合

1. 現在のディレクトリが`examples/`であることを確認
   ```bash
   pwd  # /path/to/inverse_msmd/examples であるべき
   ```

2. サンプルデータが`../data/`に配置されていることを確認
   ```bash
   ls ../data/sample_proteins/
   ls ../data/sample_probes/
   ```

### RDKit が見つからない場合

RDKitをインストール：
```bash
conda install -c conda-forge rdkit
```

### gridData が見つからない場合

gridDataをインストール：
```bash
pip install gridData
```

### メモリ不足エラー

大きなプロファイルファイルを扱う場合、メモリ不足になることがあります：
- `granularity`パラメータを小さくする
- 処理するプローブ数を減らす
- より多くのメモリを持つマシンで実行

### 重ね合わせ結果が期待と異なる

1. 原子マッチングデータを確認
2. アイソトープラベルが正しく設定されているか確認
3. MCS検出結果を確認（atom_matching.pyの出力）
4. 複数のマッチングがある場合、適切なものを選択

## さらなる情報

- パッケージAPI: [`../inverse_msmd/`](../inverse_msmd/)
- メインREADME: [`../README.md`](../README.md)
- テストコード: [`../tests/`](../tests/)