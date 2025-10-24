# プロファイルスコア計算機能 テスト手順書

最終更新: 2025-10-23

## 📋 テスト概要

このドキュメントは、プロファイルスコア計算機能の統合に関する全テスト項目と手順を記載しています。

### テスト戦略

- **ユニットテスト**: 個別関数の動作確認
- **統合テスト**: 複数モジュールの連携確認
- **回帰テスト**: 既存機能への影響確認
- **手動テスト**: 視覚的確認が必要な項目

## 🧪 テスト環境

### 必要なデータ

| データ | パス | 用途 |
|--------|------|------|
| サンプルリガンド | `data/atom_matching/4hw3_A_lig.sdf` | 部分構造置換のテスト |
| サンプルタンパク質 | `data/sample_proteins/4hw3_A.pdb` | スコア計算のテスト |
| E23プローブ | `data/sample_probes/E23.*` | 置換元構造 |
| E24プローブ | `data/sample_probes/E24.*` | 置換先構造 |
| E24プロファイル | `data/profiles/E24_*_profile.dx.gz` | スコア計算 |

### テスト実行環境

```bash
# 仮想環境の有効化
source venv/bin/activate

# 依存関係のインストール
pip install -r requirements.txt
pip install pytest pytest-cov

# テストディレクトリの準備
mkdir -p test_output
```

## 📝 テストチェックリスト

### Phase 1: プロファイルスコア計算モジュール

#### T1-1: モジュール基本構造のテスト

**ファイル:** `tests/unit/test_profile_scoring.py`

**テスト項目:**

- [ ] **インポートテスト**
  ```python
  def test_import_profile_scoring():
      """profile_scoringモジュールがインポートできる"""
      from inverse_msmd.profile_scoring import calculate_profile_score
      assert callable(calculate_profile_score)
  ```
  
  **期待結果:** モジュールが正常にインポートされる

---

#### T1-2: プロファイル読み込みテスト

- [ ] **正常系: プロファイルファイルの読み込み**
  ```python
  def test_load_profiles():
      """プロファイルファイルが正しく読み込まれる"""
      # E24プロファイルの読み込みテスト
      from gridData import Grid
      profile = Grid("data/profiles/E24_ALA_profile.dx.gz")
      assert profile.grid.shape[0] > 0
  ```
  
  **期待結果:** 全19種類のプロファイルが読み込まれる

- [ ] **異常系: プロファイルディレクトリが存在しない**
  ```python
  def test_profile_dir_not_found():
      """存在しないディレクトリでエラーが発生する"""
      with pytest.raises(ValueError):
          calculate_profile_score(
              protein, probe_center,
              "nonexistent_dir", "E24"
          )
  ```
  
  **期待結果:** `ValueError`が発生し、適切なエラーメッセージが表示される

---

#### T1-3: Cβ原子抽出テスト

- [ ] **正常系: Cβ原子の抽出**
  ```python
  def test_extract_cb_atoms():
      """Cβ原子が正しく抽出される"""
      from inverse_msmd.utils.bio_utils import PDB
      protein = PDB.get_structure("data/sample_proteins/4hw3_A.pdb")
      cb_atoms = [a for a in protein.get_atoms() if a.get_name() == "CB"]
      assert len(cb_atoms) > 0
  ```
  
  **期待結果:** Cβ原子が抽出される（GLYを除く残基数と一致）

- [ ] **異常系: Cβ原子が存在しない**
  ```python
  def test_no_cb_atoms():
      """Cβ原子がない場合のエラー処理"""
      # GLYのみのタンパク質でテスト
      with pytest.raises(ValueError):
          calculate_profile_score(gly_only_protein, ...)
  ```
  
  **期待結果:** `ValueError`が発生

---

#### T1-4: スコア計算テスト

- [ ] **正常系: スコア計算**
  ```python
  def test_calculate_score():
      """スコアが正しく計算される"""
      score = calculate_profile_score(
          protein, probe_center,
          "data/profiles/", "E24"
      )
      assert isinstance(score, float)
      assert score < 0  # 対数スコアは通常負
  ```
  
  **期待結果:** floatの負の値が返される

- [ ] **エッジケース: 負のプロファイル値の処理**
  ```python
  def test_negative_profile_values():
      """負のプロファイル値が適切に処理される"""
      # 最小値で置換されることを確認
      score = calculate_profile_score(protein, probe_center, ...)
      assert not np.isnan(score)
      assert not np.isinf(score)
  ```
  
  **期待結果:** NaNやInfが発生しない

---

### Phase 2: 統合ワークフローの拡張

#### T2-1: 後方互換性テスト

**ファイル:** `tests/integration/test_profile_integration.py`

- [ ] **スコア計算なし（既存機能）**
  ```python
  def test_backward_compatibility():
      """profile_dir=Noneで既存の動作が維持される"""
      results = integrated_substructure_replacement(
          ligand_file="data/atom_matching/4hw3_A_lig.sdf",
          protein_file="data/sample_proteins/4hw3_A.pdb",
          from_file="data/sample_probes/E23",
          to_file="data/sample_probes/E24",
          output_dir="test_output/backward/",
          profile_dir=None  # スコア計算なし
      )
      assert len(results) > 0
      assert 'score' not in results[0]  # スコアが含まれていない
  ```
  
  **期待結果:** 
  - 既存の動作と同じ結果が得られる
  - スコア情報が含まれていない

---

#### T2-2: スコア計算統合テスト

- [ ] **正常系: スコア計算あり**
  ```python
  def test_with_profile_scoring():
      """スコア計算が正しく統合される"""
      results = integrated_substructure_replacement(
          ligand_file="data/atom_matching/4hw3_A_lig.sdf",
          protein_file="data/sample_proteins/4hw3_A.pdb",
          from_file="data/sample_probes/E23",
          to_file="data/sample_probes/E24",
          output_dir="test_output/with_score/",
          profile_dir="data/profiles/",
          probe_id="E24"
      )
      assert len(results) > 0
      assert 'score' in results[0]
      assert isinstance(results[0]['score'], float)
  ```
  
  **期待結果:**
  - 全パターンにスコアが付与される
  - スコアはfloat型

- [ ] **スコアソートテスト**
  ```python
  def test_score_sorting():
      """結果がスコアで降順ソートされる"""
      results = integrated_substructure_replacement(
          ...,
          profile_dir="data/profiles/",
          probe_id="E24"
      )
      if len(results) > 1:
          for i in range(len(results) - 1):
              assert results[i]['score'] >= results[i+1]['score']
  ```
  
  **期待結果:** スコアが高い順に並んでいる

---

#### T2-3: パラメータバリデーションテスト

- [ ] **異常系: probe_idが指定されていない**
  ```python
  def test_missing_probe_id():
      """profile_dir指定時にprobe_idがないとエラー"""
      with pytest.raises(ValueError) as exc_info:
          integrated_substructure_replacement(
              ...,
              profile_dir="data/profiles/",
              probe_id=None  # エラー
          )
      assert "probe_id" in str(exc_info.value).lower()
  ```
  
  **期待結果:** 適切なエラーメッセージが表示される

---

### Phase 3: パッケージエクスポート

#### T3-1: インポートテスト

**ファイル:** `tests/unit/test_imports.py`

- [ ] **個別機能のインポート**
  ```python
  def test_import_individual_functions():
      """個別機能がインポートできる"""
      from inverse_msmd import (
          find_substructure_in_ligand,
          match_substructures,
          calculate_profile_score,
          integrated_substructure_replacement
      )
      assert all(callable(f) for f in [
          find_substructure_in_ligand,
          match_substructures,
          calculate_profile_score,
          integrated_substructure_replacement
      ])
  ```
  
  **期待結果:** 全ての関数がインポートできる

- [ ] **バージョン確認**
  ```python
  def test_version():
      """バージョンが更新されている"""
      import inverse_msmd
      assert inverse_msmd.__version__ == "0.2.0"
  ```
  
  **期待結果:** バージョンが0.2.0

---

### Phase 4: CLIスクリプト

#### T4-1: コマンドライン引数テスト

- [ ] **ヘルプメッセージ**
  ```bash
  python scripts/integrated_replacement.py --help
  ```
  
  **期待結果:**
  - `--profile-dir`オプションが表示される
  - `--probe-id`オプションが表示される

- [ ] **スコア計算なしの実行**
  ```bash
  python scripts/integrated_replacement.py \
      --ligand data/atom_matching/4hw3_A_lig.sdf \
      --protein data/sample_proteins/4hw3_A.pdb \
      --from-file data/sample_probes/E23 \
      --to-file data/sample_probes/E24 \
      --output test_output/cli_no_score/
  ```
  
  **期待結果:** 正常に実行され、ファイルが出力される

- [ ] **スコア計算ありの実行**
  ```bash
  python scripts/integrated_replacement.py \
      --ligand data/atom_matching/4hw3_A_lig.sdf \
      --protein data/sample_proteins/4hw3_A.pdb \
      --from-file data/sample_probes/E23 \
      --to-file data/sample_probes/E24 \
      --output test_output/cli_with_score/ \
      --profile-dir data/profiles/ \
      --probe-id E24 \
      --verbose
  ```
  
  **期待結果:**
  - スコアが表示される
  - スコアで降順ソート済みと表示される

---

### Phase 5: 統合テスト

#### T5-1: エンドツーエンドテスト

- [ ] **完全なワークフロー**
  ```python
  def test_complete_workflow():
      """完全なワークフローが正常に動作する"""
      # 1. 統合ワークフローの実行
      results = integrated_substructure_replacement(
          ligand_file="data/atom_matching/4hw3_A_lig.sdf",
          protein_file="data/sample_proteins/4hw3_A.pdb",
          from_file="data/sample_probes/E23",
          to_file="data/sample_probes/E24",
          output_dir="test_output/e2e/",
          profile_dir="data/profiles/",
          probe_id="E24"
      )
      
      # 2. 結果の検証
      assert len(results) > 0
      assert all('score' in r for r in results)
      
      # 3. ファイルの存在確認
      for result in results:
          assert Path(result['ligand_file']).exists()
          assert Path(result['protein_file']).exists()
      
      # 4. スコアの範囲確認
      scores = [r['score'] for r in results]
      assert all(isinstance(s, float) for s in scores)
      assert all(s < 0 for s in scores)  # 対数スコア
  ```
  
  **期待結果:** 全ての処理が正常に完了する

---

## 🔍 手動確認項目

以下の項目は自動テストでカバーできないため、手動で確認が必要です。

### 視覚的確認

- [ ] **出力ファイルの確認**
  - 生成されたSDFファイルをChimeraなどで開く
  - リガンドの部分構造が正しく置換されているか確認
  - タンパク質の座標が正しく変換されているか確認

- [ ] **スコアの妥当性確認**
  - スコアが高いパターンほど、視覚的に良好な構造か確認
  - 立体障害がないか確認

### パフォーマンス確認

- [ ] **処理時間の測定**
  ```python
  import time
  start = time.time()
  results = integrated_substructure_replacement(...)
  elapsed = time.time() - start
  print(f"処理時間: {elapsed:.2f}秒")
  ```
  
  **期待結果:** 合理的な時間で完了する（目安: 1パターンあたり1秒以内）

---

## 🧹 テスト実行手順

### 1. 全テストの実行

```bash
# 全テストを実行
pytest tests/ -v

# カバレッジ付きで実行
pytest tests/ --cov=inverse_msmd --cov-report=html

# 特定のテストのみ実行
pytest tests/unit/test_profile_scoring.py -v
```

### 2. 統合テストのみ実行

```bash
pytest tests/integration/ -v
```

### 3. CLIテスト

```bash
# ヘルプの確認
python scripts/integrated_replacement.py --help

# サンプルデータで実行
bash run_tests.sh  # テストスクリプトがある場合
```

---

## ✅ 完了条件

### テスト成功基準

- [ ] 全ユニットテストがパスする
- [ ] 全統合テストがパスする
- [ ] コードカバレッジが80%以上
- [ ] 手動確認項目が全て完了
- [ ] パフォーマンスが許容範囲内
- [ ] ドキュメントが最新

### テスト完了時のチェックリスト

- [ ] テストコードがリポジトリにコミットされている
- [ ] テスト結果がドキュメント化されている
- [ ] 発見されたバグが全て修正されている
- [ ] CI/CDパイプラインが正常に動作する（設定されている場合）

---

## 📊 テスト結果の記録

### テスト実行ログ

```
実行日: YYYY-MM-DD
実行者: [AI/人間の名前]
結果: [成功/失敗]
カバレッジ: XX%
備考: [特記事項]
```

### 発見されたバグ

現在、バグは報告されていません。

---

## 🔗 関連ドキュメント

- [`implementation_progress.md`](implementation_progress.md) - 実装進捗
- [`architecture.md`](architecture.md) - アーキテクチャ
- [`../../tests/README.md`](../../tests/README.md) - テスト全般の説明

---

## 💡 トラブルシューティング

### よくある問題

**問題1:** `gridData`のインポートエラー
```
解決策: pip install GridDataFormats
```

**問題2:** テストデータが見つからない
```
解決策: プロジェクトルートから実行しているか確認
cd /workspaces/inverse_msmd
pytest tests/
```

**問題3:** スコアがNaNになる
```
原因: プロファイル値が負またはゼロ
解決策: 最小値での置換ロジックを確認