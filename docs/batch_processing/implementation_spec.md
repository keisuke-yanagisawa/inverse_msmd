# バッチ処理機能 実装仕様書

## モジュール設計

### ファイル構成

```
inverse_msmd/
├── batch_processing.py          # バッチ処理モジュール（新規作成）
├── substructure_replacement.py  # 既存（変更なし）
└── profile_scoring.py           # 既存（変更なし）

scripts/
├── run_batch.py                 # CLIスクリプト（新規作成）
└── README.md                    # 更新

examples/
├── batch_processing_example.py  # サンプルスクリプト（新規作成）
└── batch_config_sample.csv      # サンプル設定ファイル（新規作成）
```

## データ構造

### 1. BatchJob（データクラス）

```python
from dataclasses import dataclass
from typing import Optional

@dataclass
class BatchJob:
    """バッチ処理の1つのジョブを表すデータクラス"""
    
    job_id: str
    from_probe: str
    to_probe: str
    match_index: int
    comment: Optional[str] = None
    enabled: bool = True
    
    def __post_init__(self):
        """バリデーション"""
        if not self.job_id:
            raise ValueError("job_idは空にできません")
        if self.match_index < 0:
            raise ValueError("match_indexは0以上である必要があります")
```

### 2. JobResult（データクラス）

```python
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
from datetime import datetime

@dataclass
class JobResult:
    """ジョブ実行結果を表すデータクラス"""
    
    job_id: str
    from_probe: str
    to_probe: str
    match_index: int
    status: str  # "success", "failed", "skipped"
    
    # 成功時の情報
    num_patterns: int = 0
    best_score: Optional[float] = None
    best_pattern_index: Optional[int] = None
    patterns: List[Dict[str, Any]] = field(default_factory=list)
    
    # 実行情報
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    execution_time: Optional[float] = None
    
    # エラー情報
    error_message: Optional[str] = None
    error_type: Optional[str] = None
    
    def to_dict(self) -> Dict[str, Any]:
        """辞書形式に変換"""
        return {
            'job_id': self.job_id,
            'from_probe': self.from_probe,
            'to_probe': self.to_probe,
            'match_index': self.match_index,
            'status': self.status,
            'num_patterns': self.num_patterns,
            'best_score': self.best_score,
            'best_pattern_index': self.best_pattern_index,
            'execution_time': self.execution_time,
            'error_message': self.error_message
        }
```

### 3. BatchResult（データクラス）

```python
@dataclass
class BatchResult:
    """バッチ処理全体の結果を表すデータクラス"""
    
    batch_csv: str
    total_jobs: int
    num_success: int = 0
    num_failed: int = 0
    num_skipped: int = 0
    
    job_results: List[JobResult] = field(default_factory=list)
    
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    total_execution_time: Optional[float] = None
    
    def summary(self) -> Dict[str, Any]:
        """サマリー辞書を返す"""
        return {
            'batch_csv': self.batch_csv,
            'total_jobs': self.total_jobs,
            'num_success': self.num_success,
            'num_failed': self.num_failed,
            'num_skipped': self.num_skipped,
            'total_execution_time': self.total_execution_time
        }
```

## 主要関数の詳細仕様

### 1. load_batch_config()

```python
def load_batch_config(
    batch_csv: str,
    validate: bool = True
) -> List[BatchJob]:
    """
    バッチ設定CSVを読み込みます。
    
    Parameters
    ----------
    batch_csv : str
        バッチ設定CSVファイルのパス
    validate : bool, default=True
        読み込み時にバリデーションを実行するか
    
    Returns
    -------
    List[BatchJob]
        バッチジョブのリスト
    
    Raises
    ------
    FileNotFoundError
        CSVファイルが見つからない場合
    ValueError
        CSVの形式が不正な場合
    
    Notes
    -----
    必須列: job_id, from_probe, to_probe, match_index
    オプション列: comment, enabled
    
    enabled列が'no'または'false'の場合、そのジョブはスキップされます
    
    Examples
    --------
    >>> jobs = load_batch_config("batch_config.csv")
    >>> print(f"読み込んだジョブ数: {len(jobs)}")
    """
```

**実装例:**

```python
import csv
from pathlib import Path

def load_batch_config(batch_csv: str, validate: bool = True) -> List[BatchJob]:
    batch_path = Path(batch_csv)
    
    if not batch_path.exists():
        raise FileNotFoundError(f"バッチCSVが見つかりません: {batch_csv}")
    
    jobs = []
    seen_job_ids = set()
    
    with open(batch_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        
        # 必須列チェック
        required_columns = {'job_id', 'from_probe', 'to_probe', 'match_index'}
        if not required_columns.issubset(reader.fieldnames):
            missing = required_columns - set(reader.fieldnames)
            raise ValueError(f"必須列が不足: {missing}")
        
        for row_num, row in enumerate(reader, start=2):  # ヘッダー行の次から
            # enabled列のチェック
            enabled_str = row.get('enabled', 'yes').lower()
            enabled = enabled_str not in ['no', 'false', '0']
            
            # BatchJobオブジェクトを作成
            try:
                job = BatchJob(
                    job_id=row['job_id'].strip(),
                    from_probe=row['from_probe'].strip(),
                    to_probe=row['to_probe'].strip(),
                    match_index=int(row['match_index']),
                    comment=row.get('comment', '').strip() or None,
                    enabled=enabled
                )
            except (ValueError, KeyError) as e:
                raise ValueError(f"行{row_num}のパースエラー: {e}")
            
            # job_idの重複チェック
            if validate:
                if job.job_id in seen_job_ids:
                    raise ValueError(f"job_idが重複: {job.job_id}")
                seen_job_ids.add(job.job_id)
            
            jobs.append(job)
    
    return jobs
```

### 2. process_single_job()

```python
def process_single_job(
    job: BatchJob,
    ligand_file: str,
    protein_file: str,
    probe_base_dir: str,
    profile_base_dir: str,
    output_base_dir: str,
    logger: Optional[logging.Logger] = None
) -> JobResult:
    """
    単一ジョブを処理します。
    
    Parameters
    ----------
    job : BatchJob
        処理するジョブ
    ligand_file : str
        リガンドファイルパス
    protein_file : str
        タンパク質ファイルパス
    probe_base_dir : str
        プローブベースディレクトリ
    profile_base_dir : str
        プロファイルベースディレクトリ
    output_base_dir : str
        出力ベースディレクトリ
    logger : Optional[logging.Logger]
        ロガーインスタンス
    
    Returns
    -------
    JobResult
        ジョブ実行結果
    
    Examples
    --------
    >>> job = BatchJob("job_001", "E23", "E24", 0)
    >>> result = process_single_job(
    ...     job, "ligand.sdf", "protein.pdb",
    ...     "data/probes", "data/profiles", "output"
    ... )
    >>> print(f"ステータス: {result.status}")
    """
```

**実装例:**

```python
import time
from datetime import datetime
from pathlib import Path
from inverse_msmd.substructure_replacement import integrated_substructure_replacement

def process_single_job(
    job: BatchJob,
    ligand_file: str,
    protein_file: str,
    probe_base_dir: str,
    profile_base_dir: str,
    output_base_dir: str,
    logger: Optional[logging.Logger] = None
) -> JobResult:
    
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # ジョブが無効化されている場合
    if not job.enabled:
        logger.info(f"ジョブ {job.job_id} はスキップされました（enabled=False）")
        return JobResult(
            job_id=job.job_id,
            from_probe=job.from_probe,
            to_probe=job.to_probe,
            match_index=job.match_index,
            status="skipped"
        )
    
    logger.info(f"ジョブ {job.job_id} を開始: {job.from_probe} → {job.to_probe} (index={job.match_index})")
    
    result = JobResult(
        job_id=job.job_id,
        from_probe=job.from_probe,
        to_probe=job.to_probe,
        match_index=job.match_index,
        status="failed",  # デフォルトは失敗、成功したら更新
        start_time=datetime.now()
    )
    
    start_time = time.time()
    
    try:
        # プローブファイルパスを構築
        from_file = str(Path(probe_base_dir) / job.from_probe)
        to_file = str(Path(probe_base_dir) / job.to_probe)
        
        # 出力ディレクトリ
        job_output_dir = str(Path(output_base_dir) / job.job_id)
        
        # CSV出力パス
        csv_output = str(Path(job_output_dir) / "results.csv")
        
        # 統合ワークフローを実行
        patterns = integrated_substructure_replacement(
            ligand_file=ligand_file,
            protein_file=protein_file,
            from_file=from_file,
            to_file=to_file,
            output_dir=job_output_dir,
            match_index=job.match_index,
            profile_dir=profile_base_dir,
            probe_id=job.to_probe,
            csv_output=csv_output
        )
        
        # 結果を更新
        result.status = "success"
        result.num_patterns = len(patterns)
        result.patterns = patterns
        
        # スコアがある場合、最高スコアを記録
        if patterns and 'score' in patterns[0]:
            best_pattern = max(patterns, key=lambda x: x.get('score', float('-inf')))
            result.best_score = best_pattern['score']
            result.best_pattern_index = best_pattern['pattern_index']
        
        logger.info(f"ジョブ {job.job_id} 完了: {len(patterns)} パターン生成")
        
    except Exception as e:
        result.error_message = str(e)
        result.error_type = type(e).__name__
        logger.error(f"ジョブ {job.job_id} 失敗: {e}")
    
    finally:
        end_time = time.time()
        result.end_time = datetime.now()
        result.execution_time = end_time - start_time
    
    return result
```

### 3. run_batch_processing()

```python
def run_batch_processing(
    batch_csv: str,
    ligand_file: str,
    protein_file: str,
    probe_base_dir: str,
    profile_base_dir: str,
    output_base_dir: str,
    parallel: bool = False,
    max_workers: int = 4,
    continue_on_error: bool = True,
    log_file: Optional[str] = None
) -> BatchResult:
    """
    バッチ処理を実行します。
    
    Parameters
    ----------
    batch_csv : str
        バッチ設定CSVファイルのパス
    ligand_file : str
        リガンドSDFファイルのパス
    protein_file : str
        タンパク質PDBファイルのパス
    probe_base_dir : str
        プローブファイルのベースディレクトリ
    profile_base_dir : str
        プロファイルファイルのベースディレクトリ
    output_base_dir : str
        出力ベースディレクトリ
    parallel : bool, default=False
        並列処理を有効にするか
    max_workers : int, default=4
        並列処理時の最大ワーカー数
    continue_on_error : bool, default=True
        エラー発生時も処理を継続するか
    log_file : Optional[str], default=None
        ログファイルのパス（Noneの場合は標準出力のみ）
    
    Returns
    -------
    BatchResult
        バッチ処理結果
    
    Examples
    --------
    >>> result = run_batch_processing(
    ...     batch_csv="batch_config.csv",
    ...     ligand_file="ligand.sdf",
    ...     protein_file="protein.pdb",
    ...     probe_base_dir="data/probes",
    ...     profile_base_dir="data/profiles",
    ...     output_base_dir="output"
    ... )
    >>> print(f"成功: {result.num_success}, 失敗: {result.num_failed}")
    """
```

## ユーティリティ関数

### 1. save_batch_summary()

```python
def save_batch_summary(
    batch_result: BatchResult,
    output_csv: str,
    output_json: Optional[str] = None
) -> None:
    """
    バッチ処理結果をCSVおよびJSON形式で保存します。
    
    Parameters
    ----------
    batch_result : BatchResult
        バッチ処理結果
    output_csv : str
        出力CSVファイルパス
    output_json : Optional[str]
        出力JSONファイルパス（Noneの場合はJSON出力なし）
    """
    import csv
    import json
    
    # CSV出力
    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        fieldnames = [
            'job_id', 'from_probe', 'to_probe', 'match_index',
            'status', 'num_patterns', 'best_score', 'best_pattern_index',
            'execution_time', 'error_message'
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        
        for job_result in batch_result.job_results:
            writer.writerow(job_result.to_dict())
    
    # JSON出力（オプション）
    if output_json:
        data = {
            'summary': batch_result.summary(),
            'jobs': [r.to_dict() for r in batch_result.job_results]
        }
        with open(output_json, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
```

### 2. setup_logger()

```python
def setup_logger(
    name: str,
    log_file: Optional[str] = None,
    level: int = logging.INFO
) -> logging.Logger:
    """
    ロガーを設定します。
    
    Parameters
    ----------
    name : str
        ロガー名
    log_file : Optional[str]
        ログファイルパス（Noneの場合は標準出力のみ）
    level : int
        ログレベル
    
    Returns
    -------
    logging.Logger
        設定されたロガー
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    # フォーマット
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # 標準出力ハンドラ
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # ファイルハンドラ（オプション）
    if log_file:
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger
```

## エラーハンドリング

### カスタム例外クラス

```python
class BatchProcessingError(Exception):
    """バッチ処理エラーの基底クラス"""
    pass

class BatchConfigError(BatchProcessingError):
    """バッチ設定エラー"""
    pass

class JobExecutionError(BatchProcessingError):
    """ジョブ実行エラー"""
    pass
```

### エラー処理パターン

```python
# continue_on_error=True の場合
for job in jobs:
    try:
        result = process_single_job(job, ...)
    except Exception as e:
        if not continue_on_error:
            raise JobExecutionError(f"ジョブ {job.job_id} が失敗しました") from e
        # エラーを記録して継続
        logger.error(f"ジョブ {job.job_id} でエラー: {e}")
        result = JobResult(
            job_id=job.job_id,
            status="failed",
            error_message=str(e)
        )
    results.append(result)
```

## パフォーマンス最適化

### 1. プローブファイルのキャッシング

同じプローブを複数回使用する場合、キャッシュで高速化：

```python
from functools import lru_cache

@lru_cache(maxsize=32)
def load_probe_cached(probe_file: str):
    """プローブ分子をキャッシュ付きで読み込む"""
    return read_mol_from_pdb_smi(f"{probe_file}.pdb", f"{probe_file}.smi")
```

### 2. 並列処理のメモリ管理

```python
# 大量のジョブを処理する場合、チャンクに分割
def process_in_chunks(jobs, chunk_size=10):
    for i in range(0, len(jobs), chunk_size):
        chunk = jobs[i:i+chunk_size]
        yield process_chunk_parallel(chunk)
```

## テストケース

### 1. ユニットテスト例

```python
import unittest
from inverse_msmd.batch_processing import load_batch_config, BatchJob

class TestBatchProcessing(unittest.TestCase):
    
    def test_load_batch_config_valid(self):
        """正常なCSVファイルの読み込み"""
        jobs = load_batch_config("test_data/valid_batch.csv")
        self.assertEqual(len(jobs), 3)
        self.assertEqual(jobs[0].job_id, "job_001")
    
    def test_load_batch_config_duplicate_job_id(self):
        """重複job_idのエラー検出"""
        with self.assertRaises(ValueError):
            load_batch_config("test_data/duplicate_job_id.csv")
    
    def test_batch_job_validation(self):
        """BatchJobのバリデーション"""
        with self.assertRaises(ValueError):
            BatchJob("", "E23", "E24", 0)  # 空のjob_id
        
        with self.assertRaises(ValueError):
            BatchJob("job_001", "E23", "E24", -1)  # 負のmatch_index
```

## デプロイメント

### 依存関係

`pyproject.toml`に追加する依存パッケージ：

```toml
dependencies = [
    # 既存の依存関係...
    "tqdm>=4.65.0",  # プログレスバー
]

[project.optional-dependencies]
dev = [
    # 既存の開発依存関係...
]
```

### パッケージ構造

```python
# inverse_msmd/__init__.py に追加
from .batch_processing import (
    run_batch_processing,
    load_batch_config,
    BatchJob,
    JobResult,
    BatchResult
)
```

## まとめ

この実装仕様に従って、以下の順序で実装を進めます：

1. データクラスの定義（`BatchJob`, `JobResult`, `BatchResult`）
2. ユーティリティ関数（`load_batch_config`, `save_batch_summary`, `setup_logger`）
3. コア処理関数（`process_single_job`, `run_batch_processing`）
4. CLIスクリプト（`scripts/run_batch.py`）
5. サンプルとドキュメント
6. テストコード