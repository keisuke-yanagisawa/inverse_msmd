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

### 2. 原子マッチング（atom_matching.py）

**説明:**  
プローブ分子と参照リガンドの間で最大共通部分構造（MCS）を検出するリファレンス実装です。統合APIの内部動作を理解したい場合に参照してください。

**実行方法:**
```bash
cd examples
python atom_matching.py
```

**処理内容:**
1. プローブ分子（PDB形式）と参照リガンド（SDF形式）を読み込み
2. RDKitを使用してMCS（最大共通部分構造）を検出
3. アイソトープラベル（ISO=1）を持つ原子のみを考慮
4. すべての可能なマッチングを列挙（重複除外）
5. 一致する原子のインデックスペアを保存

**主な機能:**
- RDKitライブラリを使用した分子構造の読み込み
- アイソトープラベルに基づく部分構造抽出
- MCS検出とマッチング列挙
- 原子ペアデータの生成と保存

**出力:**
```
../data/atom_matching/
├── atom_matching_A08_0    # マッチング0
├── atom_matching_A08_1    # マッチング1
├── atom_matching_E24_0
└── ...
```

各ファイルの形式：
```
# 1行目: プローブ側の原子インデックス
# 2行目: リガンド側の原子インデックス
0 1 2 3 4 5
10 11 12 13 14 15
```

**使用する場面:**
- カスタムマッチングアルゴリズムを開発する場合
- MCS検索の詳細を確認したい場合
- 統合APIとは異なるマッチング戦略を試したい場合

---

### 3. 構造重ね合わせ（superimposition.py）

**説明:**  
事前に計算された原子マッチングデータを使用して、タンパク質構造とプローブ分子を重ね合わせるリファレンス実装です。

**実行方法:**
```bash
cd examples
python superimposition.py
```

**処理内容:**
1. タンパク質構造（4hw3_A.pdb）を読み込み
2. 共結晶化リガンドとプローブ分子を読み込み
3. 各プローブとのマッチングデータ（原子ペア）を読み込み
4. プローブとリガンドの対応原子に基づいて重ね合わせを計算
5. 重ね合わせ変換をタンパク質全体に適用
6. 結果を保存（`4hw3_aligned_to_*.pdb`）

**主な機能:**
- [`SuperImposer`](../inverse_msmd/utils/bio_utils.py)クラスを使用した構造重ね合わせ
- [`PDB`](../inverse_msmd/utils/bio_utils.py)ユーティリティを使用したPDBファイルの読み書き
- 原子ペアに基づく座標マッチング
- RDKitとBioPythonの組み合わせ

**設計原則:**
- タンパク質、共結晶化リガンド、プローブは別ファイルからロード
- 各ファイルの原子インデックスは独立（0始まり）
- オフセット調整は不要
- リガンドとプローブはRDKit、タンパク質はBioPythonで処理

**使用する場面:**
- 既存の原子マッチングデータを使用する場合
- 重ね合わせの詳細をカスタマイズしたい場合
- 統合APIとは異なる処理フローを構築する場合

---

### 4. マッチングスコア計算（calculate_matching.py）

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

### 5. アイソトープラベル付与（add_isotope_labels.py）

**説明:**  
SDFファイルの分子に対して、SMARTS記法で指定した部分構造にアイソトープ番号を付与するツールです。原子マッチングで特定の部分構造のみを考慮したい場合に使用します。

**実行方法:**
```bash
cd examples

# 基本的な使用法（全マッチに適用）
python add_isotope_labels.py input.sdf output.sdf "c1ccccc1" 13

# 対話的モードで選択
python add_isotope_labels.py input.sdf output.sdf "c1ccccc1" 13 --interactive

# 特定のマッチのみに適用
python add_isotope_labels.py input.sdf output.sdf "c1ccccc1" 13 --match-index 1

# 可視化オプション付き
python add_isotope_labels.py input.sdf output.sdf "c1ccccc1" 13 --visualize
```

**主な機能:**
- SMARTS記法による部分構造指定
- 複数マッチの自動検出
- 対話的マッチ選択モード
- ラベル付与前後の分子構造可視化
- 複数分子の一括処理

**引数:**
- `input`: 入力SDFファイル
- `output`: 出力SDFファイル
- `smarts`: SMARTS記法による部分構造パターン
- `isotope`: 付与する同位体番号
- `--interactive, -i`: 対話的マッチ選択モード
- `--match-index N, -m N`: N番目のマッチのみに適用（0始まり）
- `--visualize, -v`: 分子構造を画像で出力
- `--image-dir DIR, -d DIR`: 画像出力ディレクトリ

**使用例:**

```bash
# ベンゼン環の炭素に同位体番号13を付与
python add_isotope_labels.py ligand.sdf labeled.sdf "c1ccccc1" 13

# カルボニル基の炭素に13、酸素に18を付与（2段階）
python add_isotope_labels.py input.sdf temp.sdf "[C]=[O]" 13
python add_isotope_labels.py temp.sdf output.sdf "[O]=[C]" 18

# アミン基の窒素に同位体番号15を付与
python add_isotope_labels.py input.sdf output.sdf "[NH2]" 15
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

### 詳細な制御が必要な場合

```bash
# 1. アイソトープラベルを付与（オプション）
python add_isotope_labels.py \
    ../data/atom_matching/4hw3_A_lig.sdf \
    ../data/atom_matching/4hw3_A_lig_with_subst_label.sdf \
    "your_smarts_pattern" 1

# 2. 原子マッチングを実行
python atom_matching.py

# 3. 構造重ね合わせを実行
python superimposition.py

# 4. マッチングスコアを計算
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