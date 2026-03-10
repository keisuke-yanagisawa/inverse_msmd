# チュートリアル

## 単一ジョブ: 部分構造置換 + スコア計算

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

# スコア計算付き
python scripts/integrated_replacement.py \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --from-file data/sample_probes/E23 \
    --to-file data/sample_probes/E24 \
    --output output/integrated/ \
    --profile-dir data/profiles --probe-id E24

# 特定のマッチパターンを指定
python scripts/integrated_replacement.py \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --from-file data/sample_probes/E23 \
    --to-file data/sample_probes/E24 \
    --output output/integrated/ \
    --match-index 0 \
    --verbose
```

### コマンドラインオプション

| オプション | 必須 | 説明 |
|-----------|------|------|
| `--ligand` | ✓ | リガンドSDFファイルのパス |
| `--protein` | ✓ | タンパク質PDBファイルのパス |
| `--from-file` | ✓ | 置換前の部分構造（拡張子なし、.pdbと.smiを自動読込） |
| `--to-file` | ✓ | 置換後の部分構造（拡張子なし、.pdbと.smiを自動読込） |
| `--output` | ✓ | 出力ディレクトリのパス |
| `--profile-dir` | | プロファイルディレクトリ（スコア計算時） |
| `--probe-id` | | プローブID（スコア計算時） |
| `--match-index` | | 特定のマッチパターンを指定（0始まり） |
| `--verbose` | | 詳細な進捗情報を表示 |

### 出力ファイル

- `pattern_N_ligand_replaced.sdf`: 部分構造が置換されたリガンド
- `pattern_N_protein_aligned.pdb`: 座標変換されたタンパク質
- `results.csv`: スコア計算結果（`--profile-dir` 指定時）

立体障害チェックにより、化学的に不適切な構造は自動的に除外されます。

### 実行結果の確認

```bash
ls -lh output/integrated/
```

## バッチ処理: 複数パターンの一括実行

CSVファイルで複数の置換パターンを定義し、一括処理できます。

```bash
python scripts/run_batch.py \
    --batch-csv examples/batch_config_sample.csv \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --probe-dir data/sample_probes \
    --profile-dir data/profiles \
    --output output/batch_results \
    --parallel --max-workers 12 \
    --render-figures
```

### バッチ設定CSV

```csv
job_id,from_probe,to_probe,match_index,comment,enabled
1-01,A38,A38,0,自己置換（ベースライン）,yes
1-02,A38,A01,0,A38→A01,yes
1-03,A38,A17,0,A38→A17,yes
```

### バッチ出力構造

```
output/batch_results/
├── batch_summary.csv          # 全ジョブのサマリー
├── batch_summary.json
├── panel_b_probe_map.png      # Panel B（プローブ共通）
├── 1-01/
│   ├── results.csv
│   ├── pattern_0_ligand_replaced.sdf
│   ├── pattern_0_protein_aligned.pdb
│   ├── panel_a_complex.png    # Panel A（複合体）
│   └── panel_c_combined.png   # Panel C（統合図）
├── 1-02/
│   └── ...
```
