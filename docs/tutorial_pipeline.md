# チュートリアル: MSMDプロファイル生成から逆解析まで

本チュートリアルでは、MSMDシミュレーションの結果から残基相互作用プロファイルを生成し、
inverse_msmdで逆解析を行うまでの全パイプラインを解説します。

## 全体像

```
[Step 0] MSMDシミュレーション (exprorer_msmd / GROMACS)
    入力: タンパク質PDB + プローブ分子
    出力: MDトラジェクトリ (.xtc + .pdb)
              ↓
[Step 1] プロファイル生成 (inverse_msmd.profile_generator)  ← 本チュートリアルで解説
    入力: MDトラジェクトリ群 + 参照プローブPDB
    出力: 残基相互作用プロファイル (.dx.gz)
              ↓
[Step 2] 逆解析 (inverse_msmd)
    入力: プロファイル + リガンドSDF + タンパク質PDB + プローブ分子ペア
    出力: 置換後リガンド + スコア + 3D描画
```

---

## 前提条件

```bash
# 基本インストール（逆解析のみ）
pip install -e .

# プロファイル生成も行う場合（pytraj必須）
pip install -e ".[profile]"
# または
conda install -c conda-forge pytraj
pip install joblib
```

---

## Step 0: MSMDシミュレーション（外部ツール）

MSMDシミュレーションは [exprorer_msmd](https://github.com/keisuke-yanagisawa/exprorer_msmd) で実行します。
本パッケージのスコープ外ですが、出力構造を理解するために概要を示します。

### exprorer_msmdの出力構造

```
working_dir/
├── system0/
│   ├── simulation/
│   │   ├── protein_probe.xtc    # MDトラジェクトリ
│   │   └── ...
│   └── prep/
│       └── protein_probe.pdb    # トポロジー（初期構造）
├── system1/
│   └── ...
├── ...
└── system19/
    └── ...
```

各 `systemN/` は異なる初期配置でのMDシミュレーション結果です。
統計的精度のため、通常10〜20システムを実行します。

---

## Step 1: プロファイル生成

### 概要

MDトラジェクトリから、プローブ周辺の各アミノ酸残基Cβ原子の3D密度分布を計算します。

処理の流れ:
1. **個別プロファイル生成**: 各トラジェクトリについて、プローブを原点に固定し、
   周辺残基のCβ密度を100×100×100のグリッド（1Å間隔）上に計算
2. **統合・正規化**: 全トラジェクトリの密度を合算し、バルク溶媒での期待値で割ることで
   「バルク比」に変換（>1: バルクより高頻度、<1: バルクより低頻度）
3. **圧縮**: gzip圧縮して `.dx.gz` 形式で保存

### CLIでの実行

```bash
# 基本的な使用法
python scripts/generate_profiles.py \
    --trajectories system0/simulation/protein_probe.xtc \
                   system1/simulation/protein_probe.xtc \
                   system2/simulation/protein_probe.xtc \
    --topologies   system0/prep/protein_probe.pdb \
                   system1/prep/protein_probe.pdb \
                   system2/prep/protein_probe.pdb \
    --ref-probe    probe/A17/A17.pdb \
    --probe-resname A17 \
    --probe-id     A17 \
    --output       profiles/A17

# 並列実行（4ジョブ）
python scripts/generate_profiles.py \
    --trajectories system*/simulation/protein_probe.xtc \
    --topologies   system*/prep/protein_probe.pdb \
    --ref-probe    probe/A17/A17.pdb \
    --probe-resname A17 \
    --probe-id     A17 \
    --output       profiles/A17 \
    --n-jobs 4

# 特定のアミノ酸のみ生成
python scripts/generate_profiles.py \
    --trajectories system*/simulation/protein_probe.xtc \
    --topologies   system*/prep/protein_probe.pdb \
    --ref-probe    probe/A17/A17.pdb \
    --probe-resname A17 \
    --probe-id     A17 \
    --output       profiles/A17 \
    --amino-acids ALA LEU VAL ILE PHE
```

### CLIオプション一覧

```
python scripts/generate_profiles.py --help
```

| オプション | 必須 | 説明 |
|---|---|---|
| `--trajectories` | ○ | トラジェクトリファイル (.xtc)。複数指定可 |
| `--topologies` | ○ | トポロジーファイル (.pdb)。`--trajectories` と同数 |
| `--ref-probe` | ○ | 参照プローブPDB |
| `--probe-resname` | ○ | トラジェクトリ中のプローブ残基名 |
| `--probe-id` | ○ | 出力ファイル名用のID |
| `--output` | ○ | 出力ディレクトリ |
| `--amino-acids` | | 対象アミノ酸（省略時: 全20種） |
| `--n-jobs` | | 並列ジョブ数（デフォルト: 1） |
| `--grid-size` | | グリッド座標数（デフォルト: 100） |
| `--grid-pitch` | | グリッド間隔 Å（デフォルト: 1.0） |
| `--no-latter-half` | | トラジェクトリ全体を使用 |
| `--no-compress` | | gzip圧縮しない |
| `--stop-on-error` | | エラー時に処理を中断 |
| `-v` | | 詳細ログ |

### 出力ファイル

```
profiles/A17/
├── A17_ALA_profile.dx.gz    # アラニン相互作用プロファイル
├── A17_ARG_profile.dx.gz    # アルギニン
├── A17_ASN_profile.dx.gz    # アスパラギン
├── ...                       # 全20アミノ酸
├── A17_TYR_profile.dx.gz
└── single_profiles/          # 中間ファイル（個別プロファイル）
    ├── sys0_A17_ALA_environment.dx
    ├── sys0_A17_ALA_bulk_cnt.txt
    └── ...
```

### パラメータの説明

### 注意事項

- **pytraj が必要**: `create_single_profile()` は pytraj（CPPTRAJ Python wrapper）に依存します。
  `conda install -c conda-forge pytraj` でインストールしてください。
- **メモリ**: 1トラジェクトリあたり数GB必要。`n_jobs` を上げすぎるとOOMになります。
- **レジューム**: 既に存在するファイルは自動的にスキップされます。
  中断後に再実行すると未完了分から処理を再開します。
- **HPC環境**: 大規模データではStep 1aをジョブスケジューラ（qsub等）で並列投入し、
  完了後にStep 1b, 1cを実行するのが効率的です。

---

## Step 2: 逆解析（部分構造置換 + スコアリング）

Step 1 で生成したプロファイルを使って、リガンドの部分構造置換とスコアリングを行います。

### 単件実行

```bash
python scripts/integrated_replacement.py \
    --ligand 4hw3_A_lig.sdf \
    --protein 4hw3_A.pdb \
    --from-file probe/A38 \
    --to-file probe/A17 \
    --output output/A38_to_A17/ \
    --profile-dir profiles/A17 \
    --probe-id A17
```

### バッチ実行

```bash
python scripts/run_batch.py \
    --batch-csv batch_config.csv \
    --ligand 4hw3_A_lig.sdf \
    --protein 4hw3_A.pdb \
    --probe-dir probe/ \
    --profile-dir profiles/ \
    --output output/batch_results \
    --render-figures
```

### バッチCSV形式

```csv
job_id,from_probe,to_probe,match_index,comment,enabled
1-01,A38,A38,0,自己置換（ベースライン）,yes
1-02,A38,A01,0,A38→A01,yes
1-03,A38,A17,0,A38→A17,yes
```

---

## 実践例: TYK2阻害剤のプローブ置換解析

`20251024_tyk2liang/` ディレクトリに実際の解析例があります。

### データ構造

```
20251024_tyk2liang/
├── target_protein/
│   ├── 4GIH.pdb              # TYK2タンパク質
│   └── 4gih_B_0X5.sdf        # 共結晶化リガンド
├── probe/
│   ├── A01.pdb, A01.smi      # プローブ分子群（66種）
│   ├── A08.pdb, A08.smi
│   └── ...
├── profiles/
│   ├── A01_ALA_profile.dx.gz  # 事前生成済みプロファイル
│   ├── A01_ARG_profile.dx.gz
│   └── ...
├── batch_config.csv           # バッチ設定
└── output/                    # 解析結果
```

### 実行

```bash
cd 20251024_tyk2liang

# バッチ実行
python ../scripts/run_batch.py \
    --batch-csv batch_config.csv \
    --ligand target_protein/4gih_B_0X5.sdf \
    --protein target_protein/4GIH.pdb \
    --probe-dir probe \
    --profile-dir profiles \
    --output output \
    --render-figures
```

---

## スコアの解釈

マッチングスコアは、置換後のリガンド周辺のアミノ酸配置が
MSMDプロファイル（プローブ周辺の期待される残基分布）と
どの程度一致するかを定量化します。

```
score = Σ log(profile_value_at_Cβ)
```

- **スコアが高い（0に近い）**: プロファイルとの一致が良い → 有利な置換
- **スコアが低い（大きな負値）**: プロファイルとの一致が悪い → 不利な置換
- **同じプローブへの自己置換**: ベースラインスコアとして使用

### 注意

- スコアの絶対値はプローブごとに異なるため、同一プローブ内での相対比較に使用
- 実験値（IC50等）との相関は系に依存。参考指標として使用

---

## トラブルシューティング

### pytraj が見つからない

```bash
# conda 環境で実行（pip では入らない）
conda install -c conda-forge pytraj
```

### メモリ不足 (OOM)

```bash
# 並列数を下げる
python scripts/generate_profiles.py ... --n-jobs 1

# 特定のアミノ酸だけ先に処理
python scripts/generate_profiles.py ... --amino-acids ALA LEU
```

### 中断からの再開

`generate_profiles` はレジューム機能を持っています。
既に生成済みのファイル（サイズ > 0）は自動スキップされるため、
同じコマンドを再実行するだけで未完了分から処理を再開します。

### プロファイルの確認

生成されたプロファイルは PyMOL で可視化できます:

```python
# PyMOL
load profiles/A17/A17_ALA_profile.dx, ALA_profile
isomesh ALA_mesh, ALA_profile, level=1.5
```

または inverse_msmd の描画機能:

```bash
python scripts/render_figures.py \
    --protein protein.pdb \
    --ligand ligand.sdf \
    --probe probe/A17/A17.pdb \
    --profile-dir profiles/A17 \
    --probe-id A17 \
    --output figure.png \
    --save-pse  # PyMOLセッションファイルも保存
```
