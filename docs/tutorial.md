# チュートリアル

このチュートリアルでは、リポジトリ同梱のサンプルデータを使って一通りの処理を実行します。

### サンプルデータの配置

以下のコマンド例で使用するファイルは全てリポジトリに含まれています。

| 用途 | パス |
|------|------|
| リガンド（SDF） | `data/atom_matching/4hw3_A_lig.sdf` |
| タンパク質（PDB） | `data/sample_proteins/4hw3_A.pdb` |
| プローブ分子 | `data/sample_probes/` （E23, E24, A08 等。各プローブに `.pdb` と `.smi`） |
| プロファイルマップ | `data/profiles/` （`{probe}_{AA}_profile.dx.gz` 形式） |
| バッチ設定CSV例 | `examples/batch_config_sample.csv` |

## バッチ処理

CSVファイルで複数の置換パターンを定義し、一括処理します。

```bash
python scripts/run_batch.py \
    --batch-csv examples/batch_config_sample.csv \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --probe-dir data/sample_probes \
    --profile-dir data/profiles \
    --output output/batch_results \
    --parallel --render-figures
```

### バッチ設定CSV

```csv
job_id,from_probe,to_probe,comment,enabled
exp_001,E23,E24,E23→E24,yes
exp_002,E23,A08,E23→A08,yes
```

| 列名 | 必須 | 説明 |
|------|------|------|
| `job_id` | ✓ | ジョブの一意識別子。出力ディレクトリ名に使用される |
| `from_probe` | ✓ | 置換前のプローブID。`--probe-dir` 内の `{ID}.pdb` / `{ID}.smi` を参照 |
| `to_probe` | ✓ | 置換後のプローブID。同上 |
| `match_index` | | MCSマッチの位置インデックス（0始まり）。省略時は全マッチを試行し最良スコアを自動選択 |
| `comment` | | メモ（処理には影響しない） |
| `enabled` | | `yes`（デフォルト）/ `no`。`no` にするとそのジョブをスキップ |

`from_probe` と `to_probe` には、`--probe-dir` に配置されたプローブ分子のID（拡張子なし）を指定します。各プローブには `.pdb`（3D座標）と `.smi`（SMILES）の2ファイルが必要です。

### 出力構造

```
output/batch_results/
├── batch_summary.csv          # 全ジョブのサマリー
├── batch_summary.json
├── exp_001/
│   ├── results.csv
│   ├── pattern_N_ligand_replaced.sdf
│   ├── pattern_N_protein_aligned.pdb
│   ├── pattern_N_complex.png      # Panel A（複合体）
│   ├── pattern_N_combined.png     # Panel C（統合図）
│   └── probe_map.png             # Panel B（プローブ+マップ）
├── exp_002/
│   └── ...
```

## 単一ジョブ（CSVなしで1件だけ試す）

```bash
python scripts/integrated_replacement.py \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --from-file data/sample_probes/E23 \
    --to-file data/sample_probes/E24 \
    --output output/integrated/ \
    --profile-dir data/profiles --probe-id E24
```

バッチ処理と同じオプション（`--render-figures`, `--skip-steric-clash-check` 等）が使えます。
