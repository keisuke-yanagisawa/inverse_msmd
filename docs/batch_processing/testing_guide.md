# バッチ処理機能 テスト設計

このドキュメントは、バッチ処理機能の包括的なテスト設計を定義します。

## テスト戦略

### テストレベル

1. **ユニットテスト**: 個別関数・クラスの動作確認
2. **統合テスト**: モジュール間の連携確認
3. **End-to-Endテスト**: 実際のユースケースの動作確認

### テストカバレッジ目標

- ユニットテスト: 80%以上
- 統合テスト: 主要ワークフロー100%
- エッジケースとエラーハンドリング: 100%

## テストデータ準備

### 必要なテストデータ

```
tests/data/batch_processing/
├── valid_batch.csv              # 正常なバッチ設定
├── invalid_batch_*.csv          # 異常系のバッチ設定
├── test_ligand.sdf             # テスト用リガンド
├── test_protein.pdb            # テスト用タンパク質
├── test_probes/                # テスト用プローブ
│   ├── E23.pdb
│   ├── E23.smi
│   ├── E24.pdb
│   ├── E24.smi
│   ├── A08.pdb
│   └── A08.smi
└── test_profiles/              # テスト用プロファイル
    ├── E24_ALA_profile.dx.gz
    ├── E24_ARG_profile.dx.gz
    └── ...
```

### テスト用CSVファイル

#### 1. valid_batch.csv（正常系）

```csv
job_id,from_probe,to_probe,match_index,comment,enabled
test_001,E23,E24,0,Basic test case,yes
test_002,E23,E24,1,Different match index,yes
test_003,E23,A08,0,Different probe,yes
```

#### 2. invalid_batch_missing_column.csv（異常系：必須列欠落）

```csv
job_id,from_probe,to_probe,comment
test_001,E23,E24,Missing match_index column
```

#### 3. invalid_batch_duplicate_id.csv（異常系：重複ID）

```csv
job_id,from_probe,to_probe,match_index,comment,enabled
test_001,E23,E24,0,First entry,yes
test_001,E23,A08,0,Duplicate ID,yes
```

#### 4. batch_with_disabled.csv（無効化ジョブ含む）

```csv
job_id,from_probe,to_probe,match_index,comment,enabled
test_001,E23,E24,0,Enabled job,yes
test_002,E23,A08,0,Disabled job,no
test_003,A08,E24,0,Enabled job,yes
```

## ユニットテスト設計

### tests/unit/test_batch_processing.py

#### テストクラス1: TestBatchJob

```python
import unittest
from inverse_msmd.batch_processing import BatchJob

class TestBatchJob(unittest.TestCase):
    """BatchJobデータクラスのテスト"""
    
    def test_create_valid_batch_job(self):
        """正常なBatchJobの作成"""
        job = BatchJob(
            job_id="test_001",
            from_probe="E23",
            to_probe="E24",
            match_index=0,
            comment="Test job"
        )
        self.assertEqual(job.job_id, "test_001")
        self.assertEqual(job.from_probe, "E23")
        self.assertEqual(job.to_probe, "E24")
        self.assertEqual(job.match_index, 0)
        self.assertEqual(job.comment, "Test job")
        self.assertTrue(job.enabled)
    
    def test_empty_job_id_raises_error(self):
        """空のjob_idでエラー"""
        with self.assertRaises(ValueError) as context:
            BatchJob(
                job_id="",
                from_probe="E23",
                to_probe="E24",
                match_index=0
            )
        self.assertIn("job_id", str(context.exception))
    
    def test_negative_match_index_raises_error(self):
        """負のmatch_indexでエラー"""
        with self.assertRaises(ValueError) as context:
            BatchJob(
                job_id="test_001",
                from_probe="E23",
                to_probe="E24",
                match_index=-1
            )
        self.assertIn("match_index", str(context.exception))
    
    def test_default_enabled_is_true(self):
        """デフォルトでenabledがTrue"""
        job = BatchJob(
            job_id="test_001",
            from_probe="E23",
            to_probe="E24",
            match_index=0
        )
        self.assertTrue(job.enabled)
    
    def test_optional_comment_is_none(self):
        """commentを省略した場合None"""
        job = BatchJob(
            job_id="test_001",
            from_probe="E23",
            to_probe="E24",
            match_index=0
        )
        self.assertIsNone(job.comment)
```

#### テストクラス2: TestJobResult

```python
from inverse_msmd.batch_processing import JobResult
from datetime import datetime

class TestJobResult(unittest.TestCase):
    """JobResultデータクラスのテスト"""
    
    def test_create_success_result(self):
        """成功結果の作成"""
        result = JobResult(
            job_id="test_001",
            from_probe="E23",
            to_probe="E24",
            match_index=0,
            status="success",
            num_patterns=3,
            best_score=-125.45,
            best_pattern_index=0,
            execution_time=12.34
        )
        self.assertEqual(result.status, "success")
        self.assertEqual(result.num_patterns, 3)
        self.assertAlmostEqual(result.best_score, -125.45)
    
    def test_create_failed_result(self):
        """失敗結果の作成"""
        result = JobResult(
            job_id="test_001",
            from_probe="E23",
            to_probe="E24",
            match_index=0,
            status="failed",
            error_message="部分構造が見つかりませんでした"
        )
        self.assertEqual(result.status, "failed")
        self.assertEqual(result.error_message, "部分構造が見つかりませんでした")
        self.assertEqual(result.num_patterns, 0)
    
    def test_to_dict_conversion(self):
        """辞書変換のテスト"""
        result = JobResult(
            job_id="test_001",
            from_probe="E23",
            to_probe="E24",
            match_index=0,
            status="success",
            num_patterns=2
        )
        result_dict = result.to_dict()
        
        self.assertIsInstance(result_dict, dict)
        self.assertEqual(result_dict['job_id'], "test_001")
        self.assertEqual(result_dict['status'], "success")
        self.assertEqual(result_dict['num_patterns'], 2)
```

#### テストクラス3: TestLoadBatchConfig

```python
from inverse_msmd.batch_processing import load_batch_config
import tempfile
import os

class TestLoadBatchConfig(unittest.TestCase):
    """load_batch_config関数のテスト"""
    
    def setUp(self):
        """テスト前の準備"""
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """テスト後のクリーンアップ"""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_load_valid_csv(self):
        """正常なCSVの読み込み"""
        csv_content = """job_id,from_probe,to_probe,match_index,comment,enabled
test_001,E23,E24,0,Test comment,yes
test_002,E23,A08,1,Another test,yes
"""
        csv_path = os.path.join(self.temp_dir, "valid.csv")
        with open(csv_path, 'w') as f:
            f.write(csv_content)
        
        jobs = load_batch_config(csv_path)
        
        self.assertEqual(len(jobs), 2)
        self.assertEqual(jobs[0].job_id, "test_001")
        self.assertEqual(jobs[0].from_probe, "E23")
        self.assertEqual(jobs[0].to_probe, "E24")
        self.assertEqual(jobs[0].match_index, 0)
        self.assertTrue(jobs[0].enabled)
    
    def test_missing_required_column(self):
        """必須列欠落でエラー"""
        csv_content = """job_id,from_probe,to_probe,comment
test_001,E23,E24,Missing match_index
"""
        csv_path = os.path.join(self.temp_dir, "invalid.csv")
        with open(csv_path, 'w') as f:
            f.write(csv_content)
        
        with self.assertRaises(ValueError) as context:
            load_batch_config(csv_path)
        self.assertIn("必須列", str(context.exception))
    
    def test_duplicate_job_id(self):
        """重複job_idでエラー"""
        csv_content = """job_id,from_probe,to_probe,match_index
test_001,E23,E24,0
test_001,E23,A08,0
"""
        csv_path = os.path.join(self.temp_dir, "duplicate.csv")
        with open(csv_path, 'w') as f:
            f.write(csv_content)
        
        with self.assertRaises(ValueError) as context:
            load_batch_config(csv_path)
        self.assertIn("重複", str(context.exception))
    
    def test_enabled_no_is_disabled(self):
        """enabled=noでジョブが無効化される"""
        csv_content = """job_id,from_probe,to_probe,match_index,enabled
test_001,E23,E24,0,no
"""
        csv_path = os.path.join(self.temp_dir, "disabled.csv")
        with open(csv_path, 'w') as f:
            f.write(csv_content)
        
        jobs = load_batch_config(csv_path)
        self.assertFalse(jobs[0].enabled)
    
    def test_file_not_found(self):
        """存在しないファイルでエラー"""
        with self.assertRaises(FileNotFoundError):
            load_batch_config("nonexistent.csv")
```

#### テストクラス4: TestSaveBatchSummary

```python
from inverse_msmd.batch_processing import save_batch_summary, BatchResult, JobResult
import csv
import json

class TestSaveBatchSummary(unittest.TestCase):
    """save_batch_summary関数のテスト"""
    
    def setUp(self):
        """テスト前の準備"""
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """テスト後のクリーンアップ"""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_save_csv_summary(self):
        """CSVサマリーの保存"""
        # テストデータ作成
        batch_result = BatchResult(
            batch_csv="test.csv",
            total_jobs=2,
            num_success=1,
            num_failed=1
        )
        batch_result.job_results = [
            JobResult(
                job_id="test_001",
                from_probe="E23",
                to_probe="E24",
                match_index=0,
                status="success",
                num_patterns=3,
                best_score=-125.45,
                best_pattern_index=0,
                execution_time=12.34
            ),
            JobResult(
                job_id="test_002",
                from_probe="E23",
                to_probe="A08",
                match_index=0,
                status="failed",
                error_message="Error occurred"
            )
        ]
        
        # CSV保存
        csv_path = os.path.join(self.temp_dir, "summary.csv")
        save_batch_summary(batch_result, csv_path)
        
        # CSV読み込みと検証
        with open(csv_path, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]['job_id'], "test_001")
        self.assertEqual(rows[0]['status'], "success")
        self.assertEqual(rows[1]['status'], "failed")
    
    def test_save_json_summary(self):
        """JSONサマリーの保存"""
        batch_result = BatchResult(
            batch_csv="test.csv",
            total_jobs=1,
            num_success=1
        )
        batch_result.job_results = [
            JobResult(
                job_id="test_001",
                from_probe="E23",
                to_probe="E24",
                match_index=0,
                status="success"
            )
        ]
        
        # JSON保存
        csv_path = os.path.join(self.temp_dir, "summary.csv")
        json_path = os.path.join(self.temp_dir, "summary.json")
        save_batch_summary(batch_result, csv_path, json_path)
        
        # JSON読み込みと検証
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        self.assertIn('summary', data)
        self.assertIn('jobs', data)
        self.assertEqual(data['summary']['total_jobs'], 1)
```

## 統合テスト設計

### tests/integration/test_batch_workflow.py

```python
import unittest
import tempfile
import shutil
from pathlib import Path
from inverse_msmd.batch_processing import run_batch_processing

class TestBatchProcessingWorkflow(unittest.TestCase):
    """バッチ処理の統合テスト"""
    
    @classmethod
    def setUpClass(cls):
        """テストクラス全体の準備"""
        # テストデータディレクトリのパス
        cls.test_data_dir = Path("tests/data/batch_processing")
        
        # テストデータの存在確認
        required_files = [
            cls.test_data_dir / "test_ligand.sdf",
            cls.test_data_dir / "test_protein.pdb",
            cls.test_data_dir / "test_probes" / "E23.pdb",
            cls.test_data_dir / "test_probes" / "E24.pdb",
        ]
        
        for file_path in required_files:
            if not file_path.exists():
                raise FileNotFoundError(
                    f"テストデータが見つかりません: {file_path}\n"
                    f"tests/data/batch_processing/配下にテストデータを配置してください"
                )
    
    def setUp(self):
        """各テストの前準備"""
        self.temp_output_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """各テスト後のクリーンアップ"""
        shutil.rmtree(self.temp_output_dir)
    
    def test_small_batch_sequential(self):
        """小規模バッチの順次処理"""
        # テスト用CSV作成
        csv_content = """job_id,from_probe,to_probe,match_index,comment,enabled
test_001,E23,E24,0,Test job 1,yes
test_002,E23,E24,1,Test job 2,yes
"""
        csv_path = Path(self.temp_output_dir) / "batch_config.csv"
        csv_path.write_text(csv_content)
        
        # バッチ処理実行
        result = run_batch_processing(
            batch_csv=str(csv_path),
            ligand_file=str(self.test_data_dir / "test_ligand.sdf"),
            protein_file=str(self.test_data_dir / "test_protein.pdb"),
            probe_base_dir=str(self.test_data_dir / "test_probes"),
            profile_base_dir=str(self.test_data_dir / "test_profiles"),
            output_base_dir=self.temp_output_dir,
            parallel=False
        )
        
        # 結果検証
        self.assertEqual(result.total_jobs, 2)
        self.assertGreaterEqual(result.num_success, 0)
        
        # サマリーCSVの存在確認
        summary_csv = Path(self.temp_output_dir) / "batch_summary.csv"
        self.assertTrue(summary_csv.exists())
    
    def test_batch_with_errors(self):
        """エラーを含むバッチ処理"""
        csv_content = """job_id,from_probe,to_probe,match_index,comment,enabled
test_001,E23,E24,0,Valid job,yes
test_002,E23,E24,99,Invalid match_index,yes
"""
        csv_path = Path(self.temp_output_dir) / "batch_config.csv"
        csv_path.write_text(csv_content)
        
        # continue_on_error=Trueで実行
        result = run_batch_processing(
            batch_csv=str(csv_path),
            ligand_file=str(self.test_data_dir / "test_ligand.sdf"),
            protein_file=str(self.test_data_dir / "test_protein.pdb"),
            probe_base_dir=str(self.test_data_dir / "test_probes"),
            profile_base_dir=str(self.test_data_dir / "test_profiles"),
            output_base_dir=self.temp_output_dir,
            continue_on_error=True
        )
        
        # 一部成功、一部失敗することを確認
        self.assertEqual(result.total_jobs, 2)
        self.assertGreater(result.num_failed, 0)
    
    def test_batch_with_disabled_jobs(self):
        """無効化ジョブを含むバッチ処理"""
        csv_content = """job_id,from_probe,to_probe,match_index,enabled
test_001,E23,E24,0,yes
test_002,E23,E24,1,no
test_003,E23,A08,0,yes
"""
        csv_path = Path(self.temp_output_dir) / "batch_config.csv"
        csv_path.write_text(csv_content)
        
        result = run_batch_processing(
            batch_csv=str(csv_path),
            ligand_file=str(self.test_data_dir / "test_ligand.sdf"),
            protein_file=str(self.test_data_dir / "test_protein.pdb"),
            probe_base_dir=str(self.test_data_dir / "test_probes"),
            profile_base_dir=str(self.test_data_dir / "test_profiles"),
            output_base_dir=self.temp_output_dir
        )
        
        # スキップされたジョブがあることを確認
        self.assertEqual(result.num_skipped, 1)
    
    def test_output_directory_structure(self):
        """出力ディレクトリ構造の確認"""
        csv_content = """job_id,from_probe,to_probe,match_index,enabled
test_001,E23,E24,0,yes
"""
        csv_path = Path(self.temp_output_dir) / "batch_config.csv"
        csv_path.write_text(csv_content)
        
        result = run_batch_processing(
            batch_csv=str(csv_path),
            ligand_file=str(self.test_data_dir / "test_ligand.sdf"),
            protein_file=str(self.test_data_dir / "test_protein.pdb"),
            probe_base_dir=str(self.test_data_dir / "test_probes"),
            profile_base_dir=str(self.test_data_dir / "test_profiles"),
            output_base_dir=self.temp_output_dir
        )
        
        # 必須ファイルの存在確認
        output_path = Path(self.temp_output_dir)
        self.assertTrue((output_path / "batch_summary.csv").exists())
        
        # 成功したジョブのディレクトリ確認
        if result.num_success > 0:
            job_dir = output_path / "test_001"
            self.assertTrue(job_dir.exists())
```

## パフォーマンステスト

### tests/performance/test_batch_performance.py

```python
import unittest
import time
from inverse_msmd.batch_processing import run_batch_processing

class TestBatchPerformance(unittest.TestCase):
    """パフォーマンステスト"""
    
    @unittest.skipUnless(os.getenv('RUN_PERFORMANCE_TESTS'), "Performance tests disabled")
    def test_parallel_vs_sequential(self):
        """並列処理と順次処理の速度比較"""
        # 同じバッチを順次と並列で実行して時間を比較
        
        # 順次処理
        start_time = time.time()
        result_seq = run_batch_processing(
            batch_csv="tests/data/performance_batch.csv",
            ligand_file="tests/data/test_ligand.sdf",
            protein_file="tests/data/test_protein.pdb",
            probe_base_dir="tests/data/test_probes",
            profile_base_dir="tests/data/test_profiles",
            output_base_dir="output/perf_sequential",
            parallel=False
        )
        seq_time = time.time() - start_time
        
        # 並列処理
        start_time = time.time()
        result_par = run_batch_processing(
            batch_csv="tests/data/performance_batch.csv",
            ligand_file="tests/data/test_ligand.sdf",
            protein_file="tests/data/test_protein.pdb",
            probe_base_dir="tests/data/test_probes",
            profile_base_dir="tests/data/test_profiles",
            output_base_dir="output/perf_parallel",
            parallel=True,
            max_workers=4
        )
        par_time = time.time() - start_time
        
        # 並列処理の方が速いことを確認（少なくとも1.5倍）
        print(f"\n順次処理: {seq_time:.2f}秒")
        print(f"並列処理: {par_time:.2f}秒")
        print(f"スピードアップ: {seq_time/par_time:.2f}x")
        
        self.assertLess(par_time, seq_time / 1.5)
```

## テスト実行方法

### 全テスト実行

```bash
# すべてのテストを実行
pytest tests/

# カバレッジレポート付き
pytest --cov=inverse_msmd.batch_processing --cov-report=html tests/
```

### ユニットテストのみ

```bash
pytest tests/unit/test_batch_processing.py
```

### 統合テストのみ

```bash
pytest tests/integration/test_batch_workflow.py
```

### 特定のテストクラス

```bash
pytest tests/unit/test_batch_processing.py::TestBatchJob
```

### 特定のテストメソッド

```bash
pytest tests/unit/test_batch_processing.py::TestBatchJob::test_create_valid_batch_job
```

## CI/CD統合

### GitHub Actions設定例

```.github/workflows/test_batch_processing.yml
name: Batch Processing Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    
    steps:
    - uses: actions/checkout@v2
    
    - name: Set up Python
      uses: actions/setup-python@v2
      with:
        python-version: '3.9'
    
    - name: Install dependencies
      run: |
        pip install -e .
        pip install pytest pytest-cov
    
    - name: Run tests
      run: |
        pytest tests/unit/test_batch_processing.py --cov=inverse_msmd.batch_processing
    
    - name: Upload coverage
      uses: codecov/codecov-action@v2
```

## テストチェックリスト

実装完了前に以下を確認：

### ユニットテスト
- [ ] BatchJobの作成と検証
- [ ] JobResultの作成と辞書変換
- [ ] BatchResultのサマリー生成
- [ ] load_batch_config()の正常系・異常系
- [ ] save_batch_summary()のCSV/JSON出力
- [ ] process_single_job()の正常系・異常系
- [ ] エラーハンドリングの各ケース

### 統合テスト
- [ ] 小規模バッチ（3-5ジョブ）の実行
- [ ] エラーを含むバッチの処理
- [ ] 無効化ジョブの処理
- [ ] 出力ディレクトリ構造の確認
- [ ] サマリーCSVの内容確認

### End-to-Endテスト
- [ ] 実際のデータでのバッチ処理
- [ ] 並列処理の動作確認
- [ ] CLIスクリプトの動作確認

### パフォーマンステスト（オプション）
- [ ] 並列処理の速度向上確認
- [ ] メモリ使用量の監視

## トラブルシューティング

### テストが失敗する場合

1. **テストデータの確認**
   ```bash
   ls -la tests/data/batch_processing/
   ```

2. **依存パッケージの確認**
   ```bash
   pip list | grep -E "pytest|rdkit|biopython"
   ```

3. **詳細なログ出力**
   ```bash
   pytest -v -s tests/unit/test_batch_processing.py
   ```

## 関連ドキュメント

- [設計ドキュメント](README.md)
- [実装仕様](implementation_spec.md)
- [実装計画](implementation_plan.md)