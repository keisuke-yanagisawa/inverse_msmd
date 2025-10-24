"""
バッチ処理モジュール

複数の置換パターン（from_probe, to_probe, match_index）に対して、
matching scoreを一括計算するバッチ処理機能を提供します。

主要機能:
- バッチ設定CSVの読み込み
- 複数ジョブの順次または並列処理
- 結果のCSV出力
- エラーハンドリング

使用例:
    >>> from inverse_msmd.batch_processing import run_batch_processing
    >>> 
    >>> result = run_batch_processing(
    ...     batch_csv="batch_config.csv",
    ...     ligand_file="data/ligand.sdf",
    ...     protein_file="data/protein.pdb",
    ...     probe_base_dir="data/probes",
    ...     profile_base_dir="data/profiles",
    ...     output_base_dir="output/batch"
    ... )
    >>> print(f"成功: {result.num_success}, 失敗: {result.num_failed}")
"""

from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
from datetime import datetime
import logging
import csv
import json
import time
from pathlib import Path
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

from inverse_msmd.substructure_replacement import integrated_substructure_replacement


# データクラス定義

@dataclass
class BatchJob:
    """
    バッチ処理の1つのジョブを表すデータクラス
    
    Attributes
    ----------
    job_id : str
        ジョブの一意識別子（出力ディレクトリ名に使用）
    from_probe : str
        置換前プローブID（例: "E23"）
    to_probe : str
        置換後プローブID（例: "E24"）
    match_index : int
        置換位置インデックス（0始まり）
    comment : Optional[str]
        説明（オプション、ログ出力に使用）
    enabled : bool
        ジョブの有効/無効フラグ（デフォルト: True）
    
    Raises
    ------
    ValueError
        job_idが空文字列の場合
        match_indexが負の値の場合
    
    Examples
    --------
    >>> job = BatchJob(
    ...     job_id="exp_001",
    ...     from_probe="E23",
    ...     to_probe="E24",
    ...     match_index=0,
    ...     comment="基本パターン"
    ... )
    >>> print(job.job_id)
    'exp_001'
    """
    job_id: str
    from_probe: str
    to_probe: str
    match_index: int
    comment: Optional[str] = None
    enabled: bool = True
    
    def __post_init__(self):
        """バリデーション"""
        if not self.job_id or not self.job_id.strip():
            raise ValueError("job_idは空にできません")
        if self.match_index < 0:
            raise ValueError("match_indexは0以上である必要があります")


@dataclass
class JobResult:
    """
    ジョブ実行結果を表すデータクラス
    
    Attributes
    ----------
    job_id : str
        ジョブID
    from_probe : str
        置換前プローブID
    to_probe : str
        置換後プローブID
    match_index : int
        置換位置インデックス
    status : str
        実行ステータス（"success", "failed", "skipped"）
    num_patterns : int
        生成されたパターン数（デフォルト: 0）
    best_score : Optional[float]
        最高スコア（スコア計算時のみ）
    best_pattern_index : Optional[int]
        最高スコアのパターン番号
    patterns : List[Dict[str, Any]]
        各パターンの詳細情報
    start_time : Optional[datetime]
        開始時刻
    end_time : Optional[datetime]
        終了時刻
    execution_time : Optional[float]
        実行時間（秒）
    error_message : Optional[str]
        エラーメッセージ（失敗時のみ）
    error_type : Optional[str]
        エラータイプ（例: "ValueError"）
    """
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
        """
        辞書形式に変換
        
        Returns
        -------
        Dict[str, Any]
            結果の辞書表現
        """
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


@dataclass
class BatchResult:
    """
    バッチ処理全体の結果を表すデータクラス
    
    Attributes
    ----------
    batch_csv : str
        バッチ設定CSVファイルパス
    total_jobs : int
        総ジョブ数
    num_success : int
        成功したジョブ数
    num_failed : int
        失敗したジョブ数
    num_skipped : int
        スキップされたジョブ数
    job_results : List[JobResult]
        各ジョブの実行結果
    start_time : Optional[datetime]
        バッチ処理開始時刻
    end_time : Optional[datetime]
        バッチ処理終了時刻
    total_execution_time : Optional[float]
        総実行時間（秒）
    """
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
        """
        サマリー辞書を返す
        
        Returns
        -------
        Dict[str, Any]
            バッチ処理のサマリー情報
        """
        return {
            'batch_csv': self.batch_csv,
            'total_jobs': self.total_jobs,
            'num_success': self.num_success,
            'num_failed': self.num_failed,
            'num_skipped': self.num_skipped,
            'total_execution_time': self.total_execution_time
        }


# ユーティリティ関数

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
        CSVの形式が不正な場合（必須列の欠落、job_idの重複など）
    
    Notes
    -----
    必須列: job_id, from_probe, to_probe, match_index
    オプション列: comment, enabled
    
    enabled列が'no'、'false'、'0'の場合、そのジョブのenabledはFalseになります
    
    Examples
    --------
    >>> jobs = load_batch_config("batch_config.csv")
    >>> print(f"読み込んだジョブ数: {len(jobs)}")
    読み込んだジョブ数: 5
    """
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
            raise ValueError(
                f"必須列が不足しています: {missing}\n"
                f"必須列: {required_columns}"
            )
        
        for row_num, row in enumerate(reader, start=2):  # ヘッダー行の次から
            # enabled列のチェック
            enabled_str = row.get('enabled', 'yes').lower().strip()
            enabled = enabled_str not in ['no', 'false', '0', '']
            
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
                raise ValueError(
                    f"CSVファイルの行{row_num}でエラーが発生しました: {e}\n"
                    f"行の内容: {row}"
                )
            
            # job_idの重複チェック
            if validate:
                if job.job_id in seen_job_ids:
                    raise ValueError(
                        f"job_idが重複しています: '{job.job_id}'\n"
                        f"各job_idはユニークである必要があります"
                    )
                seen_job_ids.add(job.job_id)
            
            jobs.append(job)
    
    if len(jobs) == 0:
        raise ValueError(f"CSVファイルにジョブが定義されていません: {batch_csv}")
    
    return jobs


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
    
    Examples
    --------
    >>> save_batch_summary(
    ...     batch_result,
    ...     "output/batch_summary.csv",
    ...     "output/batch_summary.json"
    ... )
    """
    # 出力ディレクトリを作成
    csv_path = Path(output_csv)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    
    # CSV出力
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
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
        json_path = Path(output_json)
        json_path.parent.mkdir(parents=True, exist_ok=True)
        
        data = {
            'summary': batch_result.summary(),
            'jobs': [r.to_dict() for r in batch_result.job_results]
        }
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False, default=str)


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
        ログレベル（デフォルト: logging.INFO）
    
    Returns
    -------
    logging.Logger
        設定されたロガー
    
    Examples
    --------
    >>> logger = setup_logger("batch_processing", "output/batch.log")
    >>> logger.info("Processing started")
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    # 既存のハンドラをクリア（重複を防ぐ）
    logger.handlers.clear()
    
    # フォーマット
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # 標準出力ハンドラ
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # ファイルハンドラ（オプション）
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger


# コア処理関数

def process_single_job(
    job: BatchJob,
    ligand_file: str,
    protein_file: str,
    probe_base_dir: str,
    profile_base_dir: str,
    output_base_dir: str,
    logger: Optional[logging.Logger] = None,
    skip_steric_clash_check: bool = False
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
    ステータス: success
    """
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
    
    # ログ出力
    comment_str = f" ({job.comment})" if job.comment else ""
    logger.info(
        f"ジョブ {job.job_id} を開始: "
        f"{job.from_probe} → {job.to_probe} "
        f"(index={job.match_index}){comment_str}"
    )
    
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
        # deduplicate_by_smiles=Trueで、異なる構造式のパターンのみを保持
        patterns = integrated_substructure_replacement(
            ligand_file=ligand_file,
            protein_file=protein_file,
            from_file=from_file,
            to_file=to_file,
            output_dir=job_output_dir,
            match_index=job.match_index,
            profile_dir=profile_base_dir,
            probe_id=job.to_probe,
            csv_output=csv_output,
            deduplicate_by_smiles=True,  # 異なるSMILESのパターンをすべて保持
            skip_steric_clash_check=skip_steric_clash_check
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
            logger.info(
                f"ジョブ {job.job_id} 完了: {len(patterns)} パターン生成、"
                f"最高スコア={result.best_score:.2f}"
            )
        else:
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


def _process_job_wrapper(args):
    """
    並列処理用のジョブ処理ラッパー関数
    
    ProcessPoolExecutorでは関数の引数をpickle化する必要があるため、
    ロガーは各プロセスで再作成します。
    
    Parameters
    ----------
    args : tuple
        (job, ligand_file, protein_file, probe_base_dir,
         profile_base_dir, output_base_dir, skip_steric_clash_check) のタプル
    
    Returns
    -------
    JobResult
        ジョブ実行結果
    """
    job, ligand_file, protein_file, probe_base_dir, profile_base_dir, output_base_dir, skip_steric_clash_check = args
    
    # 各プロセスで独自のロガーを作成（pickle化の問題を回避）
    logger = logging.getLogger(f"batch_processing.worker.{job.job_id}")
    logger.setLevel(logging.INFO)
    
    # ハンドラが既に存在する場合はクリア
    logger.handlers.clear()
    
    return process_single_job(
        job,
        ligand_file,
        protein_file,
        probe_base_dir,
        profile_base_dir,
        output_base_dir,
        logger,
        skip_steric_clash_check
    )


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
    log_file: Optional[str] = None,
    skip_steric_clash_check: bool = False
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
        並列処理時の最大ワーカー数（parallelがTrueの場合に使用）
    continue_on_error : bool, default=True
        エラー発生時も処理を継続するか
    log_file : Optional[str], default=None
        ログファイルのパス（Noneの場合は標準出力のみ）
    
    Returns
    -------
    BatchResult
        バッチ処理結果
    
    Raises
    ------
    FileNotFoundError
        バッチCSVファイルが見つからない場合
    ValueError
        バッチCSVの形式が不正な場合
    
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
    成功: 3, 失敗: 0
    """
    # ロガーのセットアップ
    if log_file is None:
        log_file = str(Path(output_base_dir) / "batch_execution.log")
    
    logger = setup_logger("batch_processing", log_file)
    
    logger.info("=" * 70)
    logger.info("バッチ処理を開始")
    logger.info("=" * 70)
    logger.info(f"バッチCSV: {batch_csv}")
    logger.info(f"出力先: {output_base_dir}")
    
    # 出力ディレクトリを作成
    output_path = Path(output_base_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # バッチ設定を読み込み
    logger.info(f"\nバッチ設定を読み込み中...")
    try:
        jobs = load_batch_config(batch_csv)
    except Exception as e:
        logger.error(f"バッチCSVの読み込みに失敗しました: {e}")
        raise
    
    logger.info(f"読み込み完了: {len(jobs)} ジョブ")
    
    # BatchResultオブジェクトを作成
    batch_result = BatchResult(
        batch_csv=batch_csv,
        total_jobs=len(jobs),
        start_time=datetime.now()
    )
    
    batch_start_time = time.time()
    
    # 各ジョブを処理
    logger.info(f"\n{'='*70}")
    logger.info(f"ジョブ処理を開始")
    if parallel:
        logger.info(f"並列処理モード: 最大{max_workers}ワーカー")
    else:
        logger.info(f"順次処理モード")
    logger.info(f"{'='*70}\n")
    
    if parallel:
        # 並列処理
        logger.info(f"{max_workers}個のワーカーで並列処理を実行します")
        
        # 各ジョブの引数をタプルとして準備
        job_args = [
            (job, ligand_file, protein_file, probe_base_dir, profile_base_dir, output_base_dir, skip_steric_clash_check)
            for job in jobs
        ]
        
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # ジョブを投入
            future_to_job = {
                executor.submit(_process_job_wrapper, args): args[0]
                for args in job_args
            }
            
            # 完了したジョブから結果を収集
            with tqdm(total=len(jobs), desc="Processing jobs") as pbar:
                for future in as_completed(future_to_job):
                    job = future_to_job[future]
                    try:
                        result = future.result()
                        batch_result.job_results.append(result)
                        
                        # ステータスカウント
                        if result.status == "success":
                            batch_result.num_success += 1
                        elif result.status == "failed":
                            batch_result.num_failed += 1
                        elif result.status == "skipped":
                            batch_result.num_skipped += 1
                        
                        pbar.update(1)
                        
                    except Exception as e:
                        logger.error(f"ジョブ {job.job_id} の処理中に予期しないエラー: {e}")
                        
                        if not continue_on_error:
                            logger.error("continue_on_error=Falseのため処理を中断します")
                            raise
                        
                        # エラー結果を記録
                        error_result = JobResult(
                            job_id=job.job_id,
                            from_probe=job.from_probe,
                            to_probe=job.to_probe,
                            match_index=job.match_index,
                            status="failed",
                            error_message=str(e),
                            error_type=type(e).__name__
                        )
                        batch_result.job_results.append(error_result)
                        batch_result.num_failed += 1
                        pbar.update(1)
    else:
        # 順次処理
        for job in tqdm(jobs, desc="Processing jobs"):
            try:
                result = process_single_job(
                    job,
                    ligand_file,
                    protein_file,
                    probe_base_dir,
                    profile_base_dir,
                    output_base_dir,
                    logger,
                    skip_steric_clash_check
                )
                
                batch_result.job_results.append(result)
                
                # ステータスカウント
                if result.status == "success":
                    batch_result.num_success += 1
                elif result.status == "failed":
                    batch_result.num_failed += 1
                elif result.status == "skipped":
                    batch_result.num_skipped += 1
                    
            except Exception as e:
                logger.error(f"ジョブ {job.job_id} の処理中に予期しないエラー: {e}")
                
                if not continue_on_error:
                    logger.error("continue_on_error=Falseのため処理を中断します")
                    raise
                
                # エラー結果を記録
                error_result = JobResult(
                    job_id=job.job_id,
                    from_probe=job.from_probe,
                    to_probe=job.to_probe,
                    match_index=job.match_index,
                    status="failed",
                    error_message=str(e),
                    error_type=type(e).__name__
                )
                batch_result.job_results.append(error_result)
                batch_result.num_failed += 1
    
    # バッチ処理終了
    batch_end_time = time.time()
    batch_result.end_time = datetime.now()
    batch_result.total_execution_time = batch_end_time - batch_start_time
    
    # 結果サマリーを保存
    logger.info(f"\n{'='*70}")
    logger.info(f"結果サマリーを保存中...")
    logger.info(f"{'='*70}")
    
    summary_csv = str(output_path / "batch_summary.csv")
    summary_json = str(output_path / "batch_summary.json")
    save_batch_summary(batch_result, summary_csv, summary_json)
    
    logger.info(f"CSV: {summary_csv}")
    logger.info(f"JSON: {summary_json}")
    
    # サマリー表示
    logger.info(f"\n{'='*70}")
    logger.info(f"バッチ処理完了")
    logger.info(f"{'='*70}")
    logger.info(f"総ジョブ数: {batch_result.total_jobs}")
    logger.info(f"成功: {batch_result.num_success}")
    logger.info(f"失敗: {batch_result.num_failed}")
    logger.info(f"スキップ: {batch_result.num_skipped}")
    logger.info(f"実行時間: {batch_result.total_execution_time:.2f}秒")
    
    return batch_result


# モジュールのエクスポート
__all__ = [
    'BatchJob',
    'JobResult',
    'BatchResult',
    'load_batch_config',
    'save_batch_summary',
    'setup_logger',
    'process_single_job',
    'run_batch_processing',
]