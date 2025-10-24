"""
file_utils モジュールのテスト
"""

import pytest
from pathlib import Path

from inverse_msmd.utils.file_utils import (
    ensure_output_dir,
    validate_file_exists,
    find_pdb_smi_pair,
    validate_output_path,
    get_file_extension,
    list_files_with_extension
)


class TestEnsureOutputDir:
    """ensure_output_dir関数のテスト"""
    
    def test_creates_directory(self, tmp_path):
        """ディレクトリが作成されることをテスト"""
        output_file = tmp_path / "output" / "results" / "data.csv"
        
        parent_dir = ensure_output_dir(output_file)
        
        assert parent_dir.exists()
        assert parent_dir.is_dir()
        assert parent_dir == tmp_path / "output" / "results"
    
    def test_nested_directories(self, tmp_path):
        """ネストされたディレクトリが作成されることをテスト"""
        output_file = tmp_path / "a" / "b" / "c" / "d" / "file.txt"
        
        parent_dir = ensure_output_dir(output_file)
        
        assert parent_dir.exists()
        assert (tmp_path / "a" / "b" / "c" / "d").exists()
    
    def test_existing_directory(self, tmp_path):
        """既存のディレクトリでエラーにならないことをテスト"""
        existing_dir = tmp_path / "existing"
        existing_dir.mkdir()
        
        output_file = existing_dir / "file.txt"
        
        parent_dir = ensure_output_dir(output_file)
        
        assert parent_dir.exists()
        assert parent_dir == existing_dir
    
    def test_current_directory(self, tmp_path):
        """カレントディレクトリの場合は何も作成しないことをテスト"""
        output_file = Path("file.txt")
        
        parent_dir = ensure_output_dir(output_file)
        
        assert parent_dir == Path('.')


class TestValidateFileExists:
    """validate_file_exists関数のテスト"""
    
    def test_existing_file(self, tmp_path):
        """存在するファイルが検証されることをテスト"""
        test_file = tmp_path / "test.txt"
        test_file.write_text("test content")
        
        validated_path = validate_file_exists(test_file, "Test file")
        
        assert validated_path.exists()
        assert validated_path == test_file
    
    def test_missing_file_raises_error(self, tmp_path):
        """存在しないファイルでエラーが発生することをテスト"""
        missing_file = tmp_path / "missing.txt"
        
        with pytest.raises(FileNotFoundError) as exc_info:
            validate_file_exists(missing_file, "Missing file")
        
        assert "Missing file not found" in str(exc_info.value)
        assert str(missing_file) in str(exc_info.value)
    
    def test_default_description(self, tmp_path):
        """デフォルトの説明が使用されることをテスト"""
        missing_file = tmp_path / "missing.txt"
        
        with pytest.raises(FileNotFoundError) as exc_info:
            validate_file_exists(missing_file)
        
        assert "File not found" in str(exc_info.value)


class TestFindPdbSmiPair:
    """find_pdb_smi_pair関数のテスト"""
    
    def test_both_files_exist(self, tmp_path):
        """両方のファイルが存在する場合のテスト"""
        base_path = tmp_path / "probe"
        pdb_file = base_path.with_suffix('.pdb')
        smi_file = base_path.with_suffix('.smi')
        
        pdb_file.write_text("PDB content")
        smi_file.write_text("SMI content")
        
        pdb, smi = find_pdb_smi_pair(base_path)
        
        assert pdb.exists()
        assert smi.exists()
        assert pdb == pdb_file
        assert smi == smi_file
    
    def test_missing_pdb_raises_error(self, tmp_path):
        """PDBファイルが存在しない場合のテスト"""
        base_path = tmp_path / "probe"
        smi_file = base_path.with_suffix('.smi')
        
        smi_file.write_text("SMI content")
        
        with pytest.raises(FileNotFoundError) as exc_info:
            find_pdb_smi_pair(base_path)
        
        assert "PDB file not found" in str(exc_info.value)
    
    def test_missing_smi_raises_error(self, tmp_path):
        """SMIファイルが存在しない場合のテスト"""
        base_path = tmp_path / "probe"
        pdb_file = base_path.with_suffix('.pdb')
        
        pdb_file.write_text("PDB content")
        
        with pytest.raises(FileNotFoundError) as exc_info:
            find_pdb_smi_pair(base_path)
        
        assert "SMI file not found" in str(exc_info.value)
    
    def test_both_missing_lists_all_errors(self, tmp_path):
        """両方のファイルが存在しない場合、全エラーが列挙されることをテスト"""
        base_path = tmp_path / "probe"
        
        with pytest.raises(FileNotFoundError) as exc_info:
            find_pdb_smi_pair(base_path)
        
        error_message = str(exc_info.value)
        assert "PDB file not found" in error_message
        assert "SMI file not found" in error_message


class TestValidateOutputPath:
    """validate_output_path関数のテスト"""
    
    def test_creates_parent_directory(self, tmp_path):
        """親ディレクトリが作成されることをテスト"""
        output_file = tmp_path / "output" / "result.txt"
        
        validated_path = validate_output_path(output_file)
        
        assert validated_path.parent.exists()
        assert validated_path == output_file
    
    def test_allows_overwrite_by_default(self, tmp_path):
        """デフォルトで上書きが許可されることをテスト"""
        existing_file = tmp_path / "existing.txt"
        existing_file.write_text("existing content")
        
        validated_path = validate_output_path(existing_file)
        
        assert validated_path.exists()
        assert validated_path == existing_file
    
    def test_prevents_overwrite_when_requested(self, tmp_path):
        """must_not_exist=Trueで上書きが防止されることをテスト"""
        existing_file = tmp_path / "existing.txt"
        existing_file.write_text("existing content")
        
        with pytest.raises(FileExistsError) as exc_info:
            validate_output_path(existing_file, must_not_exist=True)
        
        assert "already exists" in str(exc_info.value)
    
    def test_allows_new_file_with_must_not_exist(self, tmp_path):
        """must_not_exist=Trueでも新規ファイルは許可されることをテスト"""
        new_file = tmp_path / "new.txt"
        
        validated_path = validate_output_path(new_file, must_not_exist=True)
        
        assert validated_path == new_file
        assert validated_path.parent.exists()


class TestGetFileExtension:
    """get_file_extension関数のテスト"""
    
    def test_simple_extension(self):
        """シンプルな拡張子のテスト"""
        ext = get_file_extension("file.txt")
        assert ext == ".txt"
    
    def test_no_extension(self):
        """拡張子なしのテスト"""
        ext = get_file_extension("file_without_extension")
        assert ext == ""
    
    def test_multiple_dots(self):
        """複数のドットを含むファイル名のテスト"""
        ext = get_file_extension("archive.tar.gz")
        assert ext == ".gz"
    
    def test_hidden_file(self):
        """隠しファイルのテスト"""
        ext = get_file_extension(".bashrc")
        # 隠しファイルは拡張子なしとして扱われる
        assert ext == ""
    
    def test_path_object(self):
        """Pathオブジェクトのテスト"""
        ext = get_file_extension(Path("data/molecule.sdf"))
        assert ext == ".sdf"


class TestListFilesWithExtension:
    """list_files_with_extension関数のテスト"""
    
    def test_finds_files_in_directory(self, tmp_path):
        """ディレクトリ内のファイルを見つけることをテスト"""
        # テストファイルを作成
        (tmp_path / "file1.sdf").write_text("content1")
        (tmp_path / "file2.sdf").write_text("content2")
        (tmp_path / "file3.txt").write_text("content3")
        
        files = list_files_with_extension(tmp_path, ".sdf")
        
        assert len(files) == 2
        assert all(f.suffix == ".sdf" for f in files)
    
    def test_extension_without_dot(self, tmp_path):
        """ドットなしの拡張子でも動作することをテスト"""
        (tmp_path / "test.pdb").write_text("content")
        
        files = list_files_with_extension(tmp_path, "pdb")
        
        assert len(files) == 1
        assert files[0].suffix == ".pdb"
    
    def test_case_insensitive(self, tmp_path):
        """大文字小文字を区別しないことをテスト"""
        (tmp_path / "file1.PDB").write_text("content1")
        (tmp_path / "file2.pdb").write_text("content2")
        
        # .pdbで検索
        files_lower = list_files_with_extension(tmp_path, ".pdb")
        # .PDBで検索
        files_upper = list_files_with_extension(tmp_path, ".PDB")
        
        # どちらの検索でも対応するファイルが見つかる
        assert len(files_lower) >= 1  # file2.pdb
        assert len(files_upper) >= 1  # file1.PDB
        
        # 両方を合わせると2つ
        all_pdb_files = list_files_with_extension(tmp_path, ".pdb") + \
                        list_files_with_extension(tmp_path, ".PDB")
        # 重複を除去して確認
        unique_files = set(all_pdb_files)
        assert len(unique_files) == 2
    
    def test_recursive_search(self, tmp_path):
        """再帰的検索のテスト"""
        # ネストされた構造を作成
        (tmp_path / "file1.sdf").write_text("content1")
        subdir = tmp_path / "subdir"
        subdir.mkdir()
        (subdir / "file2.sdf").write_text("content2")
        
        # 非再帰的検索
        files_non_recursive = list_files_with_extension(tmp_path, ".sdf", recursive=False)
        assert len(files_non_recursive) == 1
        
        # 再帰的検索
        files_recursive = list_files_with_extension(tmp_path, ".sdf", recursive=True)
        assert len(files_recursive) == 2
    
    def test_nonexistent_directory(self, tmp_path):
        """存在しないディレクトリで空リストを返すことをテスト"""
        files = list_files_with_extension(tmp_path / "nonexistent", ".sdf")
        
        assert files == []
    
    def test_sorted_output(self, tmp_path):
        """ソートされた出力のテスト"""
        (tmp_path / "c.txt").write_text("content")
        (tmp_path / "a.txt").write_text("content")
        (tmp_path / "b.txt").write_text("content")
        
        files = list_files_with_extension(tmp_path, ".txt")
        
        filenames = [f.name for f in files]
        assert filenames == ["a.txt", "b.txt", "c.txt"]


class TestIntegration:
    """統合テスト"""
    
    def test_workflow_with_real_structure(self, tmp_path):
        """実際のファイル構造を使ったワークフローテスト"""
        # ディレクトリ構造を作成
        data_dir = tmp_path / "data" / "probes"
        data_dir.mkdir(parents=True)
        
        # プローブファイルを作成
        (data_dir / "E24.pdb").write_text("PDB content")
        (data_dir / "E24.smi").write_text("SMI content")
        
        # 出力ディレクトリを確保
        output_dir = tmp_path / "output" / "results"
        output_file = output_dir / "result.txt"
        ensure_output_dir(output_file)
        
        # ファイルペアを検索
        pdb, smi = find_pdb_smi_pair(data_dir / "E24")
        
        # 検証
        assert pdb.exists()
        assert smi.exists()
        assert output_dir.exists()
        
        # ファイルリストを取得
        pdb_files = list_files_with_extension(data_dir, ".pdb")
        assert len(pdb_files) == 1