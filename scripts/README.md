# コマンドラインツール

このディレクトリには、inverse_msmdパッケージで使用できるコマンドラインツールが含まれています。

## 統合部分構造置換ツール ⭐ NEW

### integrated_replacement.py

リガンドの部分構造を別の部分構造で置換し、タンパク質構造も適切に座標変換する統合ツールです。

#### 基本的な使用方法

```bash
python scripts/integrated_replacement.py \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --from-file data/sample_probes/E23 \
    --to-file data/sample_probes/E24 \
    --output output/integrated/
```

#### オプション

**必須オプション:**
- `--ligand LIGAND` - リガンドSDFファイルのパス
- `--protein PROTEIN` - タンパク質PDBファイルのパス
- `--from-file FROM_FILE` - 置換前の部分構造のベースパス（拡張子なし）
- `--to-file TO_FILE` - 置換後の部分構造のベースパス（拡張子なし）
- `--output OUTPUT` - 出力ディレクトリのパス

**オプション:**
- `--match-index N` - 部分構造マッチのインデックス指定（0始まり）
- `--verbose` - 詳細な出力を表示
- `--version` - バージョン情報を表示
- `--help` - ヘルプメッセージを表示

#### 使用例

**例1: 基本的な使用（詳細出力付き）**

```bash
python scripts/integrated_replacement.py \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --from-file data/sample_probes/E23 \
    --to-file data/sample_probes/E24 \
    --output output/e23_to_e24/ \
    --verbose
```

**例2: 特定のマッチを指定**

```bash
python scripts/integrated_replacement.py \
    --ligand data/atom_matching/4hw3_A_lig.sdf \
    --protein data/sample_proteins/4hw3_A.pdb \
    --from-file data/sample_probes/A01 \
    --to-file data/sample_probes/E24 \
    --output output/a01_to_e24/ \
    --match-index 1
```

#### 出力ファイル

以下のファイルが出力ディレクトリに生成されます：

- `pattern_N_ligand_replaced.sdf` - 部分構造が置換されたリガンド（最良パターン）
- `pattern_N_protein_aligned.pdb` - 座標変換されたタンパク質（最良パターン）

#### 処理フロー

1. **ファイル読み込み**: リガンド、タンパク質、部分構造ファイルを読み込み
2. **部分構造探索**: リガンド中の置換前部分構造を検索
3. **Atom Matching**: 置換前後の部分構造間で原子対応付け
4. **全パターンについて処理**: Superimpose計算、座標変換、部分構造置換、スコア計算
5. **最良パターンを出力**: スコアが最も高いパターンのみファイル出力

#### 注意事項

- `--from-file`と`--to-file`には拡張子なしのベースパスを指定してください
- `.pdb`と`.smi`ファイルが自動的に読み込まれます
- 複数のatom matchingパターンが見つかった場合、全てのパターンで結果が生成されます

---

## 部分構造置換ツール

### replace_substructure.py

リガンドの部分構造を置換するツールです（従来版）。

#### 基本的な使用方法

詳細は[`replace_substructure.py`](replace_substructure.py)のヘルプを参照してください：

```bash
python scripts/replace_substructure.py --help
```

---

## アイソトープラベル付与ツール

### add_isotope_labels.py

分子にアイソトープラベルを付与するツールです。

#### 基本的な使用方法

詳細は[`add_isotope_labels.py`](add_isotope_labels.py)のヘルプを参照してください：

```bash
python scripts/add_isotope_labels.py --help
```

---

## トラブルシューティング

### ファイルが見つからないエラー

```
エラー: 必要なファイルが見つかりません
  - リガンドファイルが見つかりません: ...
```

**解決方法:**
- ファイルパスが正しいか確認してください
- カレントディレクトリを確認してください（プロジェクトルートから実行することを推奨）
- 相対パスを使用している場合は、絶対パスを試してください

### RDKitのインポートエラー

```
エラー: モジュールのインポートに失敗しました: No module named 'rdkit'
```

**解決方法:**
```bash
conda install -c conda-forge rdkit
```

### 部分構造が見つからないエラー

```
ValueError: リガンド中に部分構造が見つかりませんでした
```

**解決方法:**
- 置換前の部分構造（--from-file）がリガンド中に実際に存在するか確認してください
- 可視化ツールで確認することをお勧めします

---

## 関連ドキュメント

- [統合部分構造置換機能の設計](../docs/integrated_replacement_plan.md)
- [実装進捗記録](../docs/implementation_progress.md)
- [テストチェックリスト](../docs/testing_checklist.md)
- [メインREADME](../README.md)