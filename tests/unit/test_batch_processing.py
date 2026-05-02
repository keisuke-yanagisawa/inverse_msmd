"""
バッチ処理モジュールのユニットテスト
"""

import pytest
import tempfile
import csv
import json
from pathlib import Path
from datetime import datetime

from inverse_msmd.batch_processing import (
    BatchJob,
    JobResult,
    BatchResult,
    load_batch_config,
    save_batch_summary,
    setup_logger,
)


class TestBatchJob:
    """BatchJobデータクラスのテスト"""
    
    def test_valid_job(self):
        """正常なジョブの作成"""
        job = BatchJob(
            job_id="test_001",
            from_probe="E23",
            to_probe="E24",
            match_index=0,
            comment="テストジョブ",
            enabled=True
        )
        assert job.job_id == "test_001"
        assert job.from_probe == "E23"
        assert job.to_probe == "E24"
        assert job.match_index == 0
        assert job.comment == "テストジョブ"
        assert job.enabled is True
    
    def test_empty_job_id(self):
        """空のjob_idでエラー"""
        with pytest.raises(ValueError, match="job_idは空にできません"):
            BatchJob(
                job_id="",
                from_probe="E23",
                to_probe="E24",
                match_index=0
            )
    
    def test_whitespace_job_id(self):
        """空白のみのjob_idでエラー"""
        with pytest.raises(ValueError, match="job_idは空にできません"):
            BatchJob(
                job_id="   ",
                from_probe="E23",
                to_probe="E24",
                match_index=0
            )
    
    def test_negative_match_index(self):
        """負のmatch_indexでエラー"""
        with pytest.raises(ValueError, match="match_indexは0以上である必要があります"):
            BatchJob(
                job_id="test_001",
                from_probe="E23",
                to_probe="E24",
                match_index=-1
            )
    
    def test_default_values(self):
        """デフォルト値の確認"""
        job = BatchJob(
            job_id="test_001",
            from_probe="E23",
            to_probe="E24",
            match_index=0
        )
        assert job.comment is None
        assert job.enabled is True


class TestJobResult:
    """JobResultデータクラスのテスト"""
    
    def test_to_dict(self):
        """to_dict()メソッドのテスト"""
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
        
        result_dict = result.to_dict()
        
        assert result_dict['job_id'] == "test_001"
        assert result_dict['from_probe'] == "E23"
        assert result_dict['to_probe'] == "E24"
        assert result_dict['match_index'] == 0
        assert result_dict['status'] == "success"
        assert result_dict['num_patterns'] == 3
        assert result_dict['best_score'] == -125.45
        assert result_dict['best_pattern_index'] == 0
        assert result_dict['execution_time'] == 12.34
        assert result_dict['error_message'] is None


class TestBatchResult:
    """BatchResultデータクラスのテスト"""
    
    def test_summary(self):
        """summary()メソッドのテスト"""
        result = BatchResult(
            batch_csv="test.csv",
            total_jobs=5,
            num_success=3,
            num_failed=1,
            num_skipped=1,
            total_execution_time=45.67
        )
        
        summary = result.summary()
        
        assert summary['batch_csv'] == "test.csv"
        assert summary['total_jobs'] == 5
        assert summary['num_success'] == 3
        assert summary['num_failed'] == 1
        assert summary['num_skipped'] == 1
        assert summary['total_execution_time'] == 45.67


class TestLoadBatchConfig:
    """load_batch_config()関数のテスト"""
    
    def test_valid_csv(self, tmp_path):
        """正常なCSVファイルの読み込み"""
        csv_file = tmp_path / "test_batch.csv"
        csv_file.write_text(
            "job_id,from_probe,to_probe,match_index,comment,enabled\n"
            "job_001,E23,E24,0,基本パターン,yes\n"
            "job_002,E23,E24,1,別の位置,yes\n"
            "job_003,E23,A08,0,プローブ変更,no\n",
            encoding='utf-8'
        )
        
        jobs = load_batch_config(str(csv_file))
        
        assert len(jobs) == 3
        assert jobs[0].job_id == "job_001"
        assert jobs[0].from_probe == "E23"
        assert jobs[0].to_probe == "E24"
        assert jobs[0].match_index == 0
        assert jobs[0].comment == "基本パターン"
        assert jobs[0].enabled is True
        
        assert jobs[2].job_id == "job_003"
        assert jobs[2].enabled is False
    
    def test_missing_file(self):
        """存在しないファイルでエラー"""
        with pytest.raises(FileNotFoundError):
            load_batch_config("nonexistent.csv")
    
    def test_missing_required_column(self, tmp_path):
        """必須列の欠落でエラー"""
        # to_probe 列を欠落させた CSV
        csv_file = tmp_path / "test_batch.csv"
        csv_file.write_text(
            "job_id,from_probe\n"
            "job_001,E23\n",
            encoding='utf-8'
        )

        with pytest.raises(ValueError, match="必須列が不足しています"):
            load_batch_config(str(csv_file))
    
    def test_duplicate_job_id(self, tmp_path):
        """job_idの重複でエラー"""
        csv_file = tmp_path / "test_batch.csv"
        csv_file.write_text(
            "job_id,from_probe,to_probe,match_index\n"
            "job_001,E23,E24,0\n"
            "job_001,E23,A08,0\n",
            encoding='utf-8'
        )
        
        with pytest.raises(ValueError, match="job_idが重複しています"):
            load_batch_config(str(csv_file))
    
    def test_empty_csv(self, tmp_path):
        """空のCSVファイルでエラー"""
        csv_file = tmp_path / "test_batch.csv"
        csv_file.write_text(
            "job_id,from_probe,to_probe,match_index\n",
            encoding='utf-8'
        )
        
        with pytest.raises(ValueError, match="CSVファイルにジョブが定義されていません"):
            load_batch_config(str(csv_file))
    
    def test_invalid_match_index(self, tmp_path):
        """不正なmatch_indexでエラー"""
        csv_file = tmp_path / "test_batch.csv"
        csv_file.write_text(
            "job_id,from_probe,to_probe,match_index\n"
            "job_001,E23,E24,abc\n",
            encoding='utf-8'
        )
        
        with pytest.raises(ValueError, match="行2でエラーが発生しました"):
            load_batch_config(str(csv_file))
    
    def test_enabled_variations(self, tmp_path):
        """enabled列の様々な値のテスト"""
        csv_file = tmp_path / "test_batch.csv"
        csv_file.write_text(
            "job_id,from_probe,to_probe,match_index,enabled\n"
            "job_001,E23,E24,0,yes\n"
            "job_002,E23,E24,0,no\n"
            "job_003,E23,E24,0,true\n"
            "job_004,E23,E24,0,false\n"
            "job_005,E23,E24,0,1\n"
            "job_006,E23,E24,0,0\n",
            encoding='utf-8'
        )
        
        jobs = load_batch_config(str(csv_file))
        
        assert jobs[0].enabled is True   # yes
        assert jobs[1].enabled is False  # no
        assert jobs[2].enabled is True   # true
        assert jobs[3].enabled is False  # false
        assert jobs[4].enabled is True   # 1
        assert jobs[5].enabled is False  # 0
    
    def test_missing_optional_columns(self, tmp_path):
        """オプション列がない場合"""
        csv_file = tmp_path / "test_batch.csv"
        csv_file.write_text(
            "job_id,from_probe,to_probe,match_index\n"
            "job_001,E23,E24,0\n",
            encoding='utf-8'
        )
        
        jobs = load_batch_config(str(csv_file))
        
        assert len(jobs) == 1
        assert jobs[0].comment is None
        assert jobs[0].enabled is True


class TestSaveBatchSummary:
    """save_batch_summary()関数のテスト"""
    
    def test_save_csv_only(self, tmp_path):
        """CSVのみ保存"""
        batch_result = BatchResult(
            batch_csv="test.csv",
            total_jobs=2,
            num_success=1,
            num_failed=1
        )
        
        batch_result.job_results = [
            JobResult(
                job_id="job_001",
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
                job_id="job_002",
                from_probe="E23",
                to_probe="A08",
                match_index=0,
                status="failed",
                error_message="部分構造が見つかりませんでした"
            )
        ]
        
        output_csv = tmp_path / "summary.csv"
        save_batch_summary(batch_result, str(output_csv))
        
        assert output_csv.exists()
        
        # CSVの内容を確認
        with open(output_csv, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        
        assert len(rows) == 2
        assert rows[0]['job_id'] == "job_001"
        assert rows[0]['status'] == "success"
        assert rows[0]['best_score'] == "-125.45"
        assert rows[1]['job_id'] == "job_002"
        assert rows[1]['status'] == "failed"
        assert rows[1]['error_message'] == "部分構造が見つかりませんでした"
    
    def test_save_csv_and_json(self, tmp_path):
        """CSVとJSON両方保存"""
        batch_result = BatchResult(
            batch_csv="test.csv",
            total_jobs=1,
            num_success=1,
            total_execution_time=12.34
        )
        
        batch_result.job_results = [
            JobResult(
                job_id="job_001",
                from_probe="E23",
                to_probe="E24",
                match_index=0,
                status="success",
                num_patterns=3
            )
        ]
        
        output_csv = tmp_path / "summary.csv"
        output_json = tmp_path / "summary.json"
        save_batch_summary(batch_result, str(output_csv), str(output_json))
        
        assert output_csv.exists()
        assert output_json.exists()
        
        # JSONの内容を確認
        with open(output_json, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        assert 'summary' in data
        assert 'jobs' in data
        assert data['summary']['total_jobs'] == 1
        assert data['summary']['num_success'] == 1
        assert len(data['jobs']) == 1
    
    def test_create_parent_directory(self, tmp_path):
        """親ディレクトリの自動作成"""
        output_csv = tmp_path / "subdir" / "summary.csv"
        
        batch_result = BatchResult(
            batch_csv="test.csv",
            total_jobs=0
        )
        
        save_batch_summary(batch_result, str(output_csv))
        
        assert output_csv.exists()
        assert output_csv.parent.exists()


class TestSetupLogger:
    """setup_logger()関数のテスト"""
    
    def test_logger_creation(self):
        """ロガーの作成"""
        logger = setup_logger("test_logger")
        
        assert logger.name == "test_logger"
        assert len(logger.handlers) > 0
    
    def test_logger_with_file(self, tmp_path):
        """ファイル出力付きロガー"""
        log_file = tmp_path / "test.log"
        logger = setup_logger("test_logger_file", str(log_file))
        
        logger.info("テストメッセージ")
        
        assert log_file.exists()
        content = log_file.read_text(encoding='utf-8')
        assert "テストメッセージ" in content
    
    def test_logger_creates_directory(self, tmp_path):
        """ログディレクトリの自動作成"""
        log_file = tmp_path / "logs" / "test.log"
        logger = setup_logger("test_logger_dir", str(log_file))
        
        logger.info("テスト")
        
        assert log_file.exists()
        assert log_file.parent.exists()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])