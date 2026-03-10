# inverse_msmd

Inverse MSMD（共溶媒分子動力学, mixed-solvent molecular dynamics）

## 概要

**リガンド中のプローブ部分構造を入れ替えて、その置換がどれだけ良いかスコアリングする手法**です。

タンパク質-リガンド複合体に対して、リガンド中のプローブ部分構造を入れ替え、MSMDプロファイルマップと重ね合わせてスコアリングします。

| (A) 複合体構造 | (B) プローブ＋プロファイルマップ | (C) 統合結果 |
|:---:|:---:|:---:|
| <img src="docs/figures/panel_a_complex.png" width="260" alt="タンパク質-リガンド複合体"> | <img src="docs/figures/panel_b_probe_map.png" width="260" alt="プローブ+プロファイルマップ"> | <img src="docs/figures/panel_c_combined.png" width="260" alt="統合結果"> |
| タンパク質＋リガンド | プローブ（緑）＋相互作用マップ | (A)を(B)の座標系に変換して統合 |

## 特徴

- **統合部分構造置換**: MCSベースの原子マッチング + 構造重ね合わせ + 立体障害チェック
- **バッチ処理**: 複数の置換パターンを並列処理し、スコアリング・3D描画まで自動実行
- **3D構造描画**: PyMOLによるタンパク質-リガンド複合体、プローブマップの自動描画

## 実行環境

**VS Code + Dev Containers** を想定しています。リポジトリを開くだけで依存関係が自動インストールされます。

```bash
# 1. リポジトリをクローン
git clone https://github.com/keisuke-yanagisawa/inverse_msmd.git

# 2. VS Code で開き、「Reopen in Container」を実行
#    → Dockerfile に基づくコンテナが起動し、pip install -e .[dev] が自動実行されます
```

Dev Containersを使わない場合は、手動でセットアップしてください：

```bash
cd inverse_msmd
pip install -e .[dev]
```

## クイックスタート

```bash
# 1件だけ試す（部分構造置換 + スコア計算）
python scripts/integrated_replacement.py \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --from-file data/sample_probes/E23 \
    --to-file data/sample_probes/E24 \
    --output output/integrated/ \
    --profile-dir data/profiles --probe-id E24

# バッチ処理（並列 + 3D描画）
python scripts/run_batch.py \
    --batch-csv examples/batch_config_sample.csv \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --probe-dir data/sample_probes \
    --profile-dir data/profiles \
    --output output/batch_results \
    --parallel --render-figures
```

**入力**: タンパク質PDB + リガンドSDF + プローブ分子ペア（置換前/後）

**出力**: 置換後リガンドSDF + 変換済みタンパク質PDB + スコア + 3Dパネル画像

## ドキュメント

- [チュートリアル](docs/tutorial.md) — 単一ジョブ・バッチ処理の詳細な使い方
- [テスト](docs/testing.md) — テストの実行方法

## 引用

本手法を使用する場合は、以下の論文を引用してください。

柳澤 渓甫, 吉野 龍ノ介, 工藤 玄己, 広川 貴次. Inverse MSMDによる化合物部分構造プロファイリングと結合親和性予測への応用. 情報処理学会研究報告, 2026-BIO-84, pp.1-8, 2026.

```bibtex
@inproceedings{yanagisawa2026a,
  author       = {{柳澤 渓甫，吉野 龍ノ介，工藤 玄己，広川 貴次}},
  title        = {Inverse MSMDによる化合物部分構造プロファイリングと結合親和性予測への応用},
  booktitle    = {情報処理学会研究報告},
  volume       = {2026-BIO-84},
  pages        = {1--8},
  year         = {2026},
  month        = {2},
  day          = {18}
}
```

## ライセンス

MIT License — [LICENSE](LICENSE)

## 著者

Keisuke Yanagisawa (yanagisawa@comp.isct.ac.jp)
