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
# exprorer_msmd出力ディレクトリを直接指定（推奨）
python scripts/generate_profiles.py \
    --msmd-dir /path/to/working_dir \
    --ref-probe probe/A17/A17.pdb \
    --probe-id A17 \
    --output profiles/A17

# exprorer_msmd出力ディレクトリ + 並列実行
python scripts/generate_profiles.py \
    --msmd-dir /path/to/working_dir \
    --ref-probe probe/A17/A17.pdb \
    --probe-id A17 \
    --output profiles/A17 \
    --n-jobs 4

# トラジェクトリ/トポロジーを個別指定（非標準のディレクトリ構造の場合）
python scripts/generate_profiles.py \
    --trajectories system*/simulation/protein_probe.xtc \
    --topologies   system*/prep/protein_probe.pdb \
    --ref-probe    probe/A17/A17.pdb \
    --probe-id     A17 \
    --output       profiles/A17 \
    --n-jobs 4

# 特定のアミノ酸のみ生成
python scripts/generate_profiles.py \
    --msmd-dir /path/to/working_dir \
    --ref-probe probe/A17/A17.pdb \
    --probe-id A17 \
    --output profiles/A17 \
    --amino-acids ALA LEU VAL ILE PHE

# プローブ残基名がprobe-idと異なる場合のみ--probe-resnameを指定
python scripts/generate_profiles.py \
    --msmd-dir /path/to/working_dir \
    --ref-probe probe/E14/E14.pdb \
    --probe-id E14 \
    --probe-resname E14X \
    --output profiles/E14
```

### CLIオプション一覧

```
python scripts/generate_profiles.py --help
```

| オプション | 必須 | 説明 |
|---|---|---|
| `--msmd-dir` | ※1 | exprorer_msmd出力ディレクトリ。`--trajectories/--topologies` の代替 |
| `--traj-pattern` | | `--msmd-dir` 使用時のトラジェクトリglob（デフォルト: `system*/simulation/*.xtc`） |
| `--topo-pattern` | | `--msmd-dir` 使用時のトポロジーglob（デフォルト: `system*/prep/*.pdb`） |
| `--trajectories` | ※1 | トラジェクトリファイル (.xtc)。複数指定可 |
| `--topologies` | ※1 | トポロジーファイル (.pdb)。`--trajectories` と同数 |
| `--ref-probe` | ○ | 参照プローブPDB |
| `--probe-id` | ○ | 出力ファイル名用のID |
| `--probe-resname` | | トラジェクトリ中のプローブ残基名（省略時: `--probe-id` と同じ） |
| `--output` | ○ | 出力ディレクトリ |
| `--amino-acids` | | 対象アミノ酸（省略時: 全20種） |
| `--n-jobs` | | 並列ジョブ数（デフォルト: 1） |
| `--grid-size` | | グリッド座標数（デフォルト: 100） |
| `--grid-pitch` | | グリッド間隔 Å（デフォルト: 1.0） |
| `--eps` | | バルク正規化時のε（デフォルト: 0.1） |
| `--no-latter-half` | | トラジェクトリ全体を使用 |
| `--no-compress` | | gzip圧縮しない |
| `--stop-on-error` | | エラー時に処理を中断 |
| `-v` | | 詳細ログ |

※1: `--msmd-dir` または `--trajectories`+`--topologies` のいずれか一方が必須（同時指定不可）

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

### HPC環境での個別ステップ実行

大規模データでは、パイプラインを個別ステップに分割してジョブスケジューラ（qsub等）で
並列投入するのが効率的です。

#### Step 1a: 個別プロファイル生成

各トラジェクトリ × 各アミノ酸の組み合わせを個別に実行:

```bash
python scripts/create_single_profile.py \
    --trajectory system0/simulation/protein_probe.xtc \
    --topology system0/prep/protein_probe.pdb \
    --ref-probe probe/A17/A17.pdb \
    --probe-resname A17 \
    --amino-acid ALA \
    --output-dx single_profiles/sys0_A17_ALA_environment.dx \
    --output-bulk-cnt single_profiles/sys0_A17_ALA_bulk_cnt.txt
```

#### Step 1b+1c: 統合・正規化・圧縮

全システムの個別プロファイルを統合し、バルク正規化 + gzip圧縮:

```bash
python scripts/combine_profiles.py \
    --profile-dxs single_profiles/sys*_A17_ALA_environment.dx \
    --bulk-cnts single_profiles/sys*_A17_ALA_bulk_cnt.txt \
    --output profiles/A17_ALA_profile.dx \
    --eps 0.1 \
    --compress
```

### 注意事項

- **pytraj が必要**: `create_single_profile.py` は pytraj（CPPTRAJ Python wrapper）に依存します。
  `conda install -c conda-forge pytraj` でインストールしてください。
- **メモリ**: 1トラジェクトリあたり数GB必要。`n_jobs` を上げすぎるとOOMになります。
- **レジューム**: 既に存在するファイルは自動的にスキップされます。
  中断後に再実行すると未完了分から処理を再開します。

---

## Step 2: 逆解析（部分構造置換 + スコアリング）

Step 1 で生成したプロファイルを使って、リガンドの部分構造置換とスコアリングを行います。

### 単件実行

```bash
inverse-msmd-run \
    --ligand target_protein/4gih_B_0X5.sdf \
    --protein target_protein/4GIH.pdb \
    --from-probe probe/A38 \
    --to-probe probe/A17 \
    --output output/A38_to_A17/ \
    --profile-dir profiles \
    --probe-id A17
```

### バッチ実行

```bash
inverse-msmd-batch \
    --batch-csv batch_config.csv \
    --ligand target_protein/4gih_B_0X5.sdf \
    --protein target_protein/4GIH.pdb \
    --probe-dir probe \
    --profile-dir profiles \
    --output output/batch_results
```

### profile-dir のディレクトリ構造

`--profile-dir` にはプロファイルファイルをフラットに配置します:

```
profiles/
├── A01_ALA_profile.dx.gz
├── A01_ARG_profile.dx.gz
├── ...
├── A17_ALA_profile.dx.gz
└── ...
```

ファイル名の規則: `{probe_id}_{アミノ酸3文字}_profile.dx.gz`

`generate_profiles.py` の出力をそのままバッチ処理に使用する場合:
```bash
# 各プローブのプロファイルを共通ディレクトリに生成
python scripts/generate_profiles.py ... --probe-id A01 --output profiles/
python scripts/generate_profiles.py ... --probe-id A17 --output profiles/
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

`20251024_tyk2liang/` ディレクトリにTYK2阻害剤を対象とした解析例があります。

### データ構造

```
20251024_tyk2liang/
├── target_protein/
│   ├── 4GIH.pdb              # TYK2タンパク質
│   └── 4gih_B_0X5.sdf        # 共結晶化リガンド（0X5）
├── probe/
│   ├── A01.pdb, A01.smi      # プローブ分子群
│   ├── A17.pdb, A17.smi
│   ├── A38.pdb, A38.smi
│   └── ...
├── profiles/
│   ├── A01_ALA_profile.dx.gz  # 事前生成済みプロファイル（フラット配置）
│   ├── A01_ARG_profile.dx.gz
│   ├── A17_ALA_profile.dx.gz
│   └── ...
├── batch_config.csv           # バッチ設定
└── output/                    # 解析結果
```

### 単件実行例

```bash
cd 20251024_tyk2liang

# A38→A17 プローブ置換
inverse-msmd-run \
    --ligand target_protein/4gih_B_0X5.sdf \
    --protein target_protein/4GIH.pdb \
    --from-probe probe/A38 \
    --to-probe probe/A17 \
    --output output/A38_to_A17/ \
    --profile-dir profiles \
    --probe-id A17
```

出力例:
```
パターン 0: スコア=28.35, SMILES=O=C(Nc1ccnc(NC(=O)C2CC2)c1)c1ccccc1Cl
```

### バッチ実行例

```bash
cd 20251024_tyk2liang

inverse-msmd-batch \
    --batch-csv batch_config.csv \
    --ligand target_protein/4gih_B_0X5.sdf \
    --protein target_protein/4GIH.pdb \
    --probe-dir probe \
    --profile-dir profiles \
    --output output
```

バッチ実行後、`output/batch_summary.csv` にスコア一覧が出力されます:

```csv
job_id,from_probe,to_probe,match_index,status,num_patterns,best_score,...
1-01,A38,A38,0,success,1,64.70,...    # 自己置換（ベースライン）
1-02,A38,A01,0,success,1,23.82,...
1-03,A38,A17,0,success,1,28.35,...
```

自己置換（A38→A38）のスコア 64.70 がベースラインとなり、
各プローブへの置換スコアとの比較で置換の有利・不利を判断します。

---

## スコアの解釈

マッチングスコアは、置換後のリガンド周辺のアミノ酸配置が
MSMDプロファイル（プローブ周辺の期待される残基分布）と
どの程度一致するかを定量化します。

```
score = Σ log(profile_value_at_Cβ)
```

- **スコアが高い**: プロファイルとの一致が良い → 有利な置換
- **スコアが低い（小さい値）**: プロファイルとの一致が悪い → 不利な置換
- **同じプローブへの自己置換**: ベースラインスコアとして使用
