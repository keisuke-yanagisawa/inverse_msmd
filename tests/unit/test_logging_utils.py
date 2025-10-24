"""
logging_utils モジュールのテスト
"""

import pytest
from io import StringIO
import sys

from inverse_msmd.utils.logging_utils import (
    verbose_print,
    format_progress_message,
    format_section_header,
    format_key_value
)


class TestVerbosePrint:
    """verbose_print関数のテスト"""
    
    def test_verbose_true_with_prefix(self, capsys):
        """verboseがTrueでプレフィックス付きの出力テスト"""
        verbose_print("テストメッセージ", verbose=True, level="INFO", prefix=True)
        captured = capsys.readouterr()
        assert "INFO: テストメッセージ" in captured.out
    
    def test_verbose_true_without_prefix(self, capsys):
        """verboseがTrueでプレフィックスなしの出力テスト"""
        verbose_print("テストメッセージ", verbose=True, prefix=False)
        captured = capsys.readouterr()
        assert "テストメッセージ" in captured.out
        assert "INFO:" not in captured.out
    
    def test_verbose_false_no_output(self, capsys):
        """verboseがFalseの場合に何も出力されないことをテスト"""
        verbose_print("このメッセージは表示されない", verbose=False)
        captured = capsys.readouterr()
        assert captured.out == ""
    
    def test_different_log_levels(self, capsys):
        """異なるログレベルのテスト"""
        levels = ["INFO", "WARNING", "ERROR", "DEBUG"]
        
        for level in levels:
            verbose_print(f"{level}メッセージ", verbose=True, level=level)
            captured = capsys.readouterr()
            assert f"{level}:" in captured.out
    
    def test_warning_level(self, capsys):
        """WARNINGレベルの出力テスト"""
        verbose_print("警告メッセージ", verbose=True, level="WARNING")
        captured = capsys.readouterr()
        assert "WARNING: 警告メッセージ" in captured.out
    
    def test_error_level(self, capsys):
        """ERRORレベルの出力テスト"""
        verbose_print("エラーメッセージ", verbose=True, level="ERROR")
        captured = capsys.readouterr()
        assert "ERROR: エラーメッセージ" in captured.out


class TestFormatProgressMessage:
    """format_progress_message関数のテスト"""
    
    def test_half_progress(self):
        """50%進捗のフォーマットテスト"""
        msg = format_progress_message(5, 10, "処理中")
        assert "50.0%" in msg
        assert "(5/10)" in msg
        assert "処理中" in msg
    
    def test_complete_progress(self):
        """100%完了のフォーマットテスト"""
        msg = format_progress_message(10, 10, "完了")
        assert "100.0%" in msg
        assert "(10/10)" in msg
        assert "完了" in msg
    
    def test_progress_bar_visual(self):
        """プログレスバーの視覚的表示テスト"""
        msg = format_progress_message(5, 10)
        # プログレスバーに'='と'.'が含まれることを確認
        assert "=" in msg
        assert "." in msg
        assert "[" in msg
        assert "]" in msg
    
    def test_zero_progress(self):
        """0%進捗のテスト"""
        msg = format_progress_message(0, 10)
        assert "0.0%" in msg
        assert "(0/10)" in msg
    
    def test_custom_width(self):
        """カスタム幅のテスト"""
        msg = format_progress_message(5, 10, width=20)
        # 幅20の場合、バーの長さが変わる
        assert "[" in msg
        assert "]" in msg
    
    def test_no_message(self):
        """追加メッセージなしのテスト"""
        msg = format_progress_message(3, 10)
        assert "30.0%" in msg
        assert "(3/10)" in msg


class TestFormatSectionHeader:
    """format_section_header関数のテスト"""
    
    def test_default_header(self):
        """デフォルトのヘッダーフォーマットテスト"""
        header = format_section_header("テストセクション")
        lines = header.split("\n")
        
        # 3行であることを確認
        assert len(lines) == 3
        # 中央行がタイトル
        assert lines[1] == "テストセクション"
        # 上下が区切り線
        assert "=" in lines[0]
        assert "=" in lines[2]
    
    def test_custom_char(self):
        """カスタム区切り文字のテスト"""
        header = format_section_header("サブセクション", char="-")
        lines = header.split("\n")
        
        assert "-" in lines[0]
        assert "-" in lines[2]
        assert lines[1] == "サブセクション"
    
    def test_custom_width(self):
        """カスタム幅のテスト"""
        header = format_section_header("タイトル", width=50)
        lines = header.split("\n")
        
        # 区切り線の長さが50であることを確認
        assert len(lines[0]) == 50
        assert len(lines[2]) == 50
    
    def test_long_title(self):
        """長いタイトルのテスト"""
        long_title = "これは非常に長いタイトルでデフォルト幅を超えます" * 2
        header = format_section_header(long_title, width=50)
        lines = header.split("\n")
        
        # 区切り線がタイトルの長さに合わせられていることを確認
        assert len(lines[0]) >= len(long_title)


class TestFormatKeyValue:
    """format_key_value関数のテスト"""
    
    def test_basic_format(self):
        """基本的なキーと値のフォーマットテスト"""
        formatted = format_key_value("キー", "値")
        assert "キー" in formatted
        assert "値" in formatted
        assert ": " in formatted
    
    def test_custom_key_width(self):
        """カスタムキー幅のテスト"""
        formatted = format_key_value("短いキー", "値", key_width=30)
        # キーが30文字幅で左詰めされていることを確認
        parts = formatted.split(": ")
        assert len(parts[0]) == 30
    
    def test_custom_separator(self):
        """カスタム区切り文字のテスト"""
        formatted = format_key_value("キー", "値", separator=" = ")
        assert " = " in formatted
        assert ": " not in formatted
    
    def test_long_key(self):
        """長いキーのテスト"""
        long_key = "これは非常に長いキー名です"
        formatted = format_key_value(long_key, "値", key_width=10)
        # キーが幅より長い場合、そのまま表示される
        assert long_key in formatted
    
    def test_numeric_values(self):
        """数値の値のテスト"""
        formatted = format_key_value("原子数", "42")
        assert "原子数" in formatted
        assert "42" in formatted
    
    def test_empty_value(self):
        """空の値のテスト"""
        formatted = format_key_value("キー", "")
        assert "キー" in formatted
        assert ": " in formatted


class TestIntegration:
    """統合テスト"""
    
    def test_combined_output(self, capsys):
        """複数の関数を組み合わせた出力テスト"""
        # ヘッダーを出力
        header = format_section_header("処理開始")
        print(header)
        
        # 詳細情報を出力
        print(format_key_value("ファイル名", "test.sdf"))
        print(format_key_value("原子数", "100"))
        
        # 進捗を出力
        for i in range(1, 11):
            verbose_print(
                format_progress_message(i, 10, "処理中"),
                verbose=True,
                prefix=False
            )
        
        captured = capsys.readouterr()
        
        # 全ての出力が含まれていることを確認
        assert "処理開始" in captured.out
        assert "ファイル名" in captured.out
        assert "test.sdf" in captured.out
        assert "100.0%" in captured.out
        assert "(10/10)" in captured.out
    
    def test_verbose_control(self, capsys):
        """verbose制御のテスト"""
        # verbose=Trueの場合
        verbose_print("表示される", verbose=True)
        captured1 = capsys.readouterr()
        assert "表示される" in captured1.out
        
        # verbose=Falseの場合
        verbose_print("表示されない", verbose=False)
        captured2 = capsys.readouterr()
        assert captured2.out == ""