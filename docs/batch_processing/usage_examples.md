# バッチ処理機能 使用例

このドキュメントでは、バッチ処理機能の具体的な使用例を示します。

## サンプルCSVファイル

### 基本的なバッチ設定

**ファイル名:** `batch_config_basic.csv`

```csv
job_id,from_probe,to_probe,match_index,comment,enabled
exp_001,E23,E24,0,基本パターン - E23からE24への置換,yes
exp_002,E23,E24,1,基本パターン - 別の置換位置,yes
exp_003,E23,A08,0,プローブ変更 - E23からA08へ,yes
exp_004,A08,E24,0,逆方向パターン - A08からE24へ,yes
exp_005,A08,A08,0,同一プローブ置換 - match位置0,yes
exp_006,A08,A08,1,同一プローブ置換 - match位置1,yes
```

### 複数プローブの組み合わせ

**ファイル名:** `batch_config_comprehensive.csv`

```csv
job_id,from_probe,to_probe,match_index,comment,enabled
job_E23_E24_0,E23,E24,0,E23→E24 位置0,yes
job_E23_E24_1,E23,E24,1,E23→E24 位置1,yes
job_E23_E24_2,E23,E24,2,E23→E24 位置2,yes
job_E23_A08_0,E23,A08,0,E23→A08 位置0,yes
job_E23_A08_1,E23,A08,1,E23→A08 位置1,yes
job_A08_E24_0,A08,E24,0,A08→E24 位置0,yes
job_A08_E24_1,A08,E24,1,A08→E24 位置1,yes
job_A08_A08_0,A08,A08,0,A08→A08 位置0,yes
job_E24_E24_0,E24,E24,0,E24→E24 位置0 - 同一プローブ,yes
job_test_disabled,E23,E24,0,このジョブは無効化されています,no
```

### 系統的なスクリーニング

**ファイル名:** `batch_config_screening.csv`

```csv
job_id,from_probe,to_probe,match_index,comment,enabled
screen_001,E23,E24,0,スクリーニング実験1,yes
screen_002,E23,E24,0,スクリーニング実験2 - 重複テスト,yes
screen_003,A01,E24,0,A01プローブテスト,yes
screen_004,A08,E24,0,A08プローブテスト,yes
screen_005,E23,A01,0,E24以外のターゲット,yes
screen_006,E23,A08,0,A08ターゲット,yes
```

## Python API使用例

### 例1: 基本的な使用

```python
from inverse_msmd.batch_processing import run_batch_processing

# バッチ処理を実行
result = run_batch_processing(
    batch_csv="experiments/batch_config_basic.csv",
    ligand_file="data/atom_matching/4hw3_A_lig.sdf",
    protein_file="data/sample_proteins/4hw3_A.pdb",
    probe_base_dir="data/sample_probes",
    profile_base_dir="data/profiles",
    output_base_dir="output/batch_results_basic"
)

# 結果のサマリーを表示
print(f"\n{'='*60}")
print(f"バッチ処理完了")
print(f"{'='*60}")
print(f"総ジョブ数: {result.total_jobs}")
print(f"成功: {result.num_success}")
print(f"失敗: {result.num_failed}")
print(f"スキップ: {result.num_skipped}")
print(f"実行時間: {result.total_execution_time:.2f}秒")

# 各ジョブの詳細を表示
print(f"\n{'='*60}")
print(f"ジョブ詳細")
print(f"{'='*60}")
for job_result in result.job_results:
    if job_result.status == "success":
        print(f"\n[成功] {job_result.job_id}")
        print(f"  {job_result.from_probe} → {job_result.to_probe} (index={job_result.match_index})")
        print(f"  パターン数: {job_result.num_patterns}")
        if job_result.best_score is not None:
            print(f"  最高スコア: {job_result.best_score:.2f} (パターン{job_result.best_pattern_index})")
        print(f"  実行時間: {job_result.execution_time:.2f}秒")
    elif job_result.status == "failed":
        print(f"\n[失敗] {job_result.job_id}")
        print(f"  エラー: {job_result.error_message}")
```

### 例2: 並列処理を使用

```python
from inverse_msmd.batch_processing import run_batch_processing

# 4つの並列ワーカーで実行
result = run_batch_processing(
    batch_csv="experiments/batch_config_comprehensive.csv",
    ligand_file="data/atom_matching/4hw3_A_lig.sdf",
    protein_file="data/sample_proteins/4hw3_A.pdb",
    probe_base_dir="data/sample_probes",
    profile_base_dir="data/profiles",
    output_base_dir="output/batch_results_parallel",
    parallel=True,
    max_workers=4
)

print(f"並列処理完了: {result.num_success}/{result.total_jobs} 成功")
```

### 例3: エラー時に停止

```python
from inverse_msmd.batch_processing import run_batch_processing

# エラーが発生したら即座に停止
try:
    result = run_batch_processing(
        batch_csv="experiments/batch_config_basic.csv",
        ligand_file="data/atom_matching/4hw3_A_lig.sdf",
        protein_file="data/sample_proteins/4hw3_A.pdb",
        probe_base_dir="data/sample_probes",
        profile_base_dir="data/profiles",
        output_base_dir="output/batch_results_strict",
        continue_on_error=False  # エラー時に停止
    )
except Exception as e:
    print(f"バッチ処理がエラーで停止しました: {e}")
```

### 例4: ログファイルを指定

```python
from inverse_msmd.batch_processing import run_batch_processing

# ログファイルを指定
result = run_batch_processing(
    batch_csv="experiments/batch_config_screening.csv",
    ligand_file="data/atom_matching/4hw3_A_lig.sdf",
    protein_file="data/sample_proteins/4hw3_A.pdb",
    probe_base_dir="data/sample_probes",
    profile_base_dir="data/profiles",
    output_base_dir="output/batch_results_screening",
    log_file="output/batch_results_screening/execution.log"
)

print(f"ログファイル: output/batch_results_screening/execution.log")
```

### 例5: 結果の分析

```python
from inverse_msmd.batch_processing import run_batch_processing
import pandas as pd
import matplotlib.pyplot as plt

# バッチ処理を実行
result = run_batch_processing(
    batch_csv="experiments/batch_config_comprehensive.csv",
    ligand_file="data/atom_matching/4hw3_A_lig.sdf",
    protein_file="data/sample_proteins/4hw3_A.pdb",
    probe_base_dir="data/sample_probes",
    profile_base_dir="data/profiles",
    output_base_dir="output/batch_analysis"
)

# 結果サマリーCSVを読み込み
summary_df = pd.read_csv("output/batch_analysis/batch_summary.csv")

# 成功したジョブのみフィルタ
success_df = summary_df[summary_df['status'] == 'success']

# スコアが最も高いジョブを特定
best_job = success_df.loc[success_df['best_score'].idxmax()]
print(f"\n最高スコアのジョブ:")
print(f"  Job ID: {best_job['job_id']}")
print(f"  {best_job['from_probe']} → {best_job['to_probe']}")
print(f"  スコア: {best_job['best_score']:.2f}")
print(f"  パターン数: {best_job['num_patterns']}")

# スコア分布をプロット
plt.figure(figsize=(10, 6))
plt.bar(success_df['job_id'], success_df['best_score'])
plt.xticks(rotation=45, ha='right')
plt.xlabel('Job ID')
plt.ylabel('Best Score')
plt.title('Batch Processing Results - Score Distribution')
plt.tight_layout()
plt.savefig('output/batch_analysis/score_distribution.png')
print(f"\nスコア分布グラフ: output/batch_analysis/score_distribution.png")

# プローブの組み合わせごとにグループ化
grouped = success_df.groupby(['from_probe', 'to_probe']).agg({
    'best_score': 'mean',
    'num_patterns': 'sum'
}).reset_index()

print(f"\nプローブ組み合わせごとの平均スコア:")
print(grouped.to_string(index=False))
```

## CLIスクリプト使用例

### 基本的な使用

```bash
# scripts/run_batch.py として実装予定
python scripts/run_batch.py \
    --batch-csv experiments/batch_config_basic.csv \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --probe-dir data/sample_probes \
    --profile-dir data/profiles \
    --output output/batch_cli
```

### 並列処理を有効化

```bash
python scripts/run_batch.py \
    --batch-csv experiments/batch_config_comprehensive.csv \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --probe-dir data/sample_probes \
    --profile-dir data/profiles \
    --output output/batch_parallel \
    --parallel \
    --max-workers 4
```

### ログファイルを指定

```bash
python scripts/run_batch.py \
    --batch-csv experiments/batch_config_screening.csv \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --probe-dir data/sample_probes \
    --profile-dir data/profiles \
    --output output/batch_screening \
    --log-file output/batch_screening/execution.log
```

### ヘルプを表示

```bash
python scripts/run_batch.py --help
```

## 出力ファイルの構造

バッチ処理を実行すると、以下のディレクトリ構造が生成されます：

```
output/batch_results/
├── batch_summary.csv          # 全ジョブのサマリー
├── batch_summary.json         # JSON形式のサマリー（オプション）
├── batch_execution.log        # 実行ログ
│
├── exp_001/                   # ジョブ exp_001 の結果
│   ├── pattern_0_ligand_replaced.sdf
│   ├── pattern_0_protein_aligned.pdb
│   ├── pattern_1_ligand_replaced.sdf
│   ├── pattern_1_protein_aligned.pdb
│   └── results.csv           # このジョブの詳細結果
│
├── exp_002/                   # ジョブ exp_002 の結果
│   ├── pattern_0_ligand_replaced.sdf
│   ├── pattern_0_protein_aligned.pdb
│   └── results.csv
│
└── exp_003/                   # ジョブ exp_003 の結果
    ├── pattern_0_ligand_replaced.sdf
    ├── pattern_0_protein_aligned.pdb
    └── results.csv
```

### batch_summary.csv の例

```csv
job_id,from_probe,to_probe,match_index,status,num_patterns,best_score,best_pattern_index,execution_time,error_message
exp_001,E23,E24,0,success,3,-125.45,0,12.34,
exp_002,E23,E24,1,success,2,-138.92,1,10.21,
exp_003,E23,A08,0,success,4,-142.67,2,15.89,
exp_004,A08,E24,0,failed,0,,,8.56,部分構造が見つかりませんでした
exp_005,A08,A08,0,success,5,-118.34,1,14.23,
exp_006,A08,A08,1,skipped,0,,,0.00,enabled=False
```

### 各ジョブの results.csv の例

```csv
pattern_index,score,ligand_smiles,ligand_file,protein_file
0,-125.45,CC(=O)Nc1ccc(O)cc1,exp_001/pattern_0_ligand_replaced.sdf,exp_001/pattern_0_protein_aligned.pdb
1,-128.23,CC(=O)Nc1ccc(O)cc1,exp_001/pattern_1_ligand_replaced.sdf,exp_001/pattern_1_protein_aligned.pdb
2,-131.89,CC(=O)Nc1ccc(O)cc1,exp_001/pattern_2_ligand_replaced.sdf,exp_001/pattern_2_protein_aligned.pdb
```

## よくある使用パターン

### パターン1: 複数の置換位置を試す

同じプローブ組み合わせで、すべての可能な置換位置を試す：

```csv
job_id,from_probe,to_probe,match_index,comment,enabled
pos_0,E23,E24,0,置換位置0,yes
pos_1,E23,E24,1,置換位置1,yes
pos_2,E23,E24,2,置換位置2,yes
pos_3,E23,E24,3,置換位置3,yes
```

### パターン2: 複数のプローブ組み合わせ

すべてのプローブの組み合わせを試す（同一位置）：

```csv
job_id,from_probe,to_probe,match_index,comment,enabled
E23_to_E24,E23,E24,0,E23→E24,yes
E23_to_A08,E23,A08,0,E23→A08,yes
A08_to_E24,A08,E24,0,A08→E24,yes
A08_to_A01,A08,A01,0,A08→A01,yes
E24_to_A01,E24,A01,0,E24→A01,yes
```

### パターン3: グリッドサーチ

プローブと位置の全組み合わせ：

```csv
job_id,from_probe,to_probe,match_index,comment,enabled
E23_E24_0,E23,E24,0,E23→E24 位置0,yes
E23_E24_1,E23,E24,1,E23→E24 位置1,yes
E23_A08_0,E23,A08,0,E23→A08 位置0,yes
E23_A08_1,E23,A08,1,E23→A08 位置1,yes
A08_E24_0,A08,E24,0,A08→E24 位置0,yes
A08_E24_1,A08,E24,1,A08→E24 位置1,yes
```

## トラブルシューティング

### エラー: "部分構造が見つかりませんでした"

**原因:** 指定した`from_probe`がリガンド中に存在しない

**解決策:**
1. `from_probe`が正しいか確認
2. `match_index`が有効な範囲内か確認
3. 手動で部分構造探索を実行して、マッチ数を確認

```python
from inverse_msmd.substructure_replacement import find_substructure_in_ligand
from rdkit import Chem
from inverse_msmd.utils.mol_utils import read_mol_from_pdb_smi

ligand = next(Chem.SDMolSupplier("data/atom_matching/4hw3_A_lig.sdf"))
probe = read_mol_from_pdb_smi("data/sample_probes/E23.pdb", "data/sample_probes/E23.smi")

matches = find_substructure_in_ligand(ligand, probe)
print(f"マッチ数: {len(matches)}")
print(f"有効なmatch_index範囲: 0-{len(matches)-1}")
```

### エラー: "プロファイルファイルが見つかりません"

**原因:** 指定した`to_probe`のプロファイルファイルが欠落

**解決策:**
1. `profile_base_dir`が正しいか確認
2. 必要なプロファイルファイルが存在するか確認

```bash
ls data/profiles/E24_*.dx.gz
```

### 並列処理でメモリ不足

**原因:** 大量のジョブを並列実行するとメモリ不足になる

**解決策:**
1. `max_workers`を減らす（例: 2 または 1）
2. バッチを小さなチャンクに分割して実行

### ジョブの一部が失敗する

**原因:** 特定のパラメータ組み合わせで問題が発生

**解決策:**
1. `batch_summary.csv`でエラーメッセージを確認
2. 失敗したジョブを個別に実行してデバッグ
3. 失敗したジョブの`enabled`を`no`に設定して再実行

## パフォーマンス最適化のヒント

### 1. 並列処理の適切な設定

```python
# CPUコア数に応じて設定
import multiprocessing
max_workers = max(1, multiprocessing.cpu_count() - 1)

result = run_batch_processing(
    ...,
    parallel=True,
    max_workers=max_workers
)
```

### 2. バッチサイズの調整

大量のジョブは小さなバッチに分割：

```bash
# バッチ1: ジョブ1-50
python scripts/run_batch.py --batch-csv batch_part1.csv ...

# バッチ2: ジョブ51-100
python scripts/run_batch.py --batch-csv batch_part2.csv ...
```

### 3. プログレス監視

`tqdm`によるプログレスバーで進捗を確認できます。

## 関連ドキュメント

- [設計ドキュメント](README.md)
- [実装仕様](implementation_spec.md)
- [既存ワークフロー](../../inverse_msmd/substructure_replacement.py)