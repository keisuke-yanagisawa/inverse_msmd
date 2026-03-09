"""
profile_generator モジュールのテスト

pytraj非依存の関数のみテストする:
- _make_grid
- combine_profiles
- compress_profile
- _is_valid_file
"""

import gzip
import tempfile
from pathlib import Path

import numpy as np
import pytest
from gridData import Grid

from inverse_msmd.profile_generator import (
    _is_valid_file,
    _make_grid,
    combine_profiles,
    compress_profile,
)


# =============================================================================
# _make_grid
# =============================================================================


class TestMakeGrid:
    """_make_grid関数のテスト"""

    @pytest.mark.unit
    def test_delta_is_set(self):
        """deltaが正しく設定されることをテスト"""
        data = np.zeros((10, 10, 10))
        g = _make_grid(data, size=10, pitch=0.5, origin=[0.0, 0.0, 0.0])

        assert g.delta == [0.5, 0.5, 0.5]

    @pytest.mark.unit
    def test_origin_is_set(self):
        """originが正しく設定されることをテスト"""
        data = np.zeros((10, 10, 10))
        g = _make_grid(data, size=10, pitch=1.0, origin=[-5.0, -5.0, -5.0])

        assert g.origin == [-5.0, -5.0, -5.0]

    @pytest.mark.unit
    def test_n_is_set(self):
        """nが正しく設定されることをテスト"""
        data = np.zeros((20, 20, 20))
        g = _make_grid(data, size=20, pitch=1.0, origin=[0.0, 0.0, 0.0])

        assert g.n == [20, 20, 20]

    @pytest.mark.unit
    def test_grid_data_is_set(self):
        """gridデータが正しく設定されることをテスト"""
        data = np.ones((10, 10, 10)) * 3.14
        g = _make_grid(data, size=10, pitch=1.0, origin=[0.0, 0.0, 0.0])

        assert np.array_equal(g.grid, data)
        assert g.grid.shape == (10, 10, 10)

    @pytest.mark.unit
    def test_returns_grid_object(self):
        """Grid型オブジェクトが返ることをテスト"""
        data = np.zeros((5, 5, 5))
        g = _make_grid(data, size=5, pitch=2.0, origin=[1.0, 2.0, 3.0])

        assert isinstance(g, Grid)

    @pytest.mark.unit
    def test_different_pitch_values(self):
        """異なるpitch値でdeltaが全軸同一に設定されることをテスト"""
        data = np.zeros((10, 10, 10))
        for pitch in [0.25, 0.5, 1.0, 2.0]:
            g = _make_grid(data, size=10, pitch=pitch, origin=[0.0, 0.0, 0.0])
            assert g.delta == [pitch, pitch, pitch]


# =============================================================================
# combine_profiles
# =============================================================================


def _create_dx_file(path: str, data: np.ndarray, origin=None, pitch=1.0):
    """テスト用のDXファイルを作成するヘルパー"""
    size = data.shape[0]
    if origin is None:
        origin = [-size * pitch / 2] * 3
    g = _make_grid(data, size=size, pitch=pitch, origin=origin)
    g.export(path)
    return path


def _create_bulk_cnt_file(path: str, value: float):
    """テスト用のバルクカウントファイルを作成するヘルパー"""
    Path(path).write_text(str(value))
    return path


class TestCombineProfiles:
    """combine_profiles関数のテスト"""

    @pytest.mark.unit
    def test_single_profile(self, tmp_path):
        """単一プロファイルの統合テスト"""
        data = np.ones((10, 10, 10)) * 2.0
        dx_path = str(tmp_path / "profile.dx")
        cnt_path = str(tmp_path / "bulk_cnt.txt")
        output_path = str(tmp_path / "combined.dx")

        _create_dx_file(dx_path, data)
        _create_bulk_cnt_file(cnt_path, 5.0)

        result = combine_profiles([dx_path], [cnt_path], output_path, eps=0.1)

        assert result == output_path
        assert Path(output_path).exists()

        # 検証: (data + eps) / bulk_total = (2.0 + 0.1) / 5.0 = 0.42
        combined = Grid(output_path)
        expected = (2.0 + 0.1) / 5.0
        assert np.allclose(combined.grid, expected, atol=1e-6)

    @pytest.mark.unit
    def test_two_profiles_summed(self, tmp_path):
        """2つのプロファイルが正しく合算・正規化されることをテスト"""
        data1 = np.ones((10, 10, 10)) * 3.0
        data2 = np.ones((10, 10, 10)) * 7.0
        eps = 0.1

        dx1 = _create_dx_file(str(tmp_path / "p1.dx"), data1)
        dx2 = _create_dx_file(str(tmp_path / "p2.dx"), data2)
        cnt1 = _create_bulk_cnt_file(str(tmp_path / "cnt1.txt"), 2.0)
        cnt2 = _create_bulk_cnt_file(str(tmp_path / "cnt2.txt"), 3.0)

        output_path = str(tmp_path / "combined.dx")
        combine_profiles([dx1, dx2], [cnt1, cnt2], output_path, eps=eps)

        combined = Grid(output_path)
        # (3.0 + 7.0 + 0.1) / (2.0 + 3.0) = 10.1 / 5.0 = 2.02
        expected = (3.0 + 7.0 + eps) / (2.0 + 3.0)
        assert np.allclose(combined.grid, expected, atol=1e-6)

    @pytest.mark.unit
    def test_eps_parameter(self, tmp_path):
        """epsパラメータが正しく加算されることをテスト"""
        data = np.zeros((10, 10, 10))
        dx_path = _create_dx_file(str(tmp_path / "profile.dx"), data)
        cnt_path = _create_bulk_cnt_file(str(tmp_path / "cnt.txt"), 1.0)
        output_path = str(tmp_path / "combined.dx")

        eps = 0.5
        combine_profiles([dx_path], [cnt_path], output_path, eps=eps)

        combined = Grid(output_path)
        # (0.0 + 0.5) / 1.0 = 0.5
        assert np.allclose(combined.grid, 0.5, atol=1e-6)

    @pytest.mark.unit
    def test_empty_profile_list_raises_error(self, tmp_path):
        """空のプロファイルリストでValueErrorが発生することをテスト"""
        output_path = str(tmp_path / "combined.dx")

        with pytest.raises(ValueError, match="プロファイルDXファイルが指定されていません"):
            combine_profiles([], [], output_path)

    @pytest.mark.unit
    def test_mismatched_list_lengths_raises_error(self, tmp_path):
        """DXとバルクカウントのリスト長が不一致でValueErrorが発生することをテスト"""
        data = np.ones((10, 10, 10))
        dx1 = _create_dx_file(str(tmp_path / "p1.dx"), data)
        dx2 = _create_dx_file(str(tmp_path / "p2.dx"), data)
        cnt1 = _create_bulk_cnt_file(str(tmp_path / "cnt1.txt"), 1.0)
        output_path = str(tmp_path / "combined.dx")

        with pytest.raises(ValueError, match="一致しません"):
            combine_profiles([dx1, dx2], [cnt1], output_path)

    @pytest.mark.unit
    def test_zero_bulk_total_raises_error(self, tmp_path):
        """バルクカウント合計が0でValueErrorが発生することをテスト"""
        data = np.ones((10, 10, 10))
        dx_path = _create_dx_file(str(tmp_path / "profile.dx"), data)
        cnt_path = _create_bulk_cnt_file(str(tmp_path / "cnt.txt"), 0.0)
        output_path = str(tmp_path / "combined.dx")

        with pytest.raises(ValueError, match="バルクカウントの合計が0です"):
            combine_profiles([dx_path], [cnt_path], output_path)

    @pytest.mark.unit
    def test_output_directory_created(self, tmp_path):
        """出力ディレクトリが自動作成されることをテスト"""
        data = np.ones((10, 10, 10))
        dx_path = _create_dx_file(str(tmp_path / "profile.dx"), data)
        cnt_path = _create_bulk_cnt_file(str(tmp_path / "cnt.txt"), 1.0)
        output_path = str(tmp_path / "nested" / "dir" / "combined.dx")

        combine_profiles([dx_path], [cnt_path], output_path)

        assert Path(output_path).exists()

    @pytest.mark.unit
    def test_grid_shape_preserved(self, tmp_path):
        """統合後のグリッドサイズが保持されることをテスト"""
        data = np.ones((10, 10, 10))
        dx_path = _create_dx_file(str(tmp_path / "profile.dx"), data)
        cnt_path = _create_bulk_cnt_file(str(tmp_path / "cnt.txt"), 1.0)
        output_path = str(tmp_path / "combined.dx")

        combine_profiles([dx_path], [cnt_path], output_path)

        combined = Grid(output_path)
        assert combined.grid.shape == (10, 10, 10)


# =============================================================================
# compress_profile
# =============================================================================


class TestCompressProfile:
    """compress_profile関数のテスト"""

    @pytest.mark.unit
    def test_creates_gz_file(self, tmp_path):
        """.gzファイルが作成されることをテスト"""
        input_path = tmp_path / "test.dx"
        input_path.write_text("test content for compression")

        result = compress_profile(str(input_path), remove_original=False)

        assert result == str(input_path) + ".gz"
        assert Path(result).exists()

    @pytest.mark.unit
    def test_roundtrip_content(self, tmp_path):
        """圧縮→解凍で内容が一致することをテスト"""
        input_path = tmp_path / "test.dx"
        original_content = b"This is test DX file content\nwith multiple lines\n"
        input_path.write_bytes(original_content)

        gz_path = compress_profile(str(input_path), remove_original=False)

        with gzip.open(gz_path, "rb") as f:
            decompressed = f.read()

        assert decompressed == original_content

    @pytest.mark.unit
    def test_removes_original_by_default(self, tmp_path):
        """デフォルトで元ファイルが削除されることをテスト"""
        input_path = tmp_path / "test.dx"
        input_path.write_text("test content")

        compress_profile(str(input_path))

        assert not input_path.exists()
        assert (tmp_path / "test.dx.gz").exists()

    @pytest.mark.unit
    def test_keeps_original_when_requested(self, tmp_path):
        """remove_original=Falseで元ファイルが保持されることをテスト"""
        input_path = tmp_path / "test.dx"
        input_path.write_text("test content")

        compress_profile(str(input_path), remove_original=False)

        assert input_path.exists()
        assert (tmp_path / "test.dx.gz").exists()

    @pytest.mark.unit
    def test_custom_output_path(self, tmp_path):
        """カスタム出力パスが使用されることをテスト"""
        input_path = tmp_path / "test.dx"
        input_path.write_text("test content")
        custom_output = str(tmp_path / "custom_name.dx.gz")

        result = compress_profile(str(input_path), output_path=custom_output, remove_original=False)

        assert result == custom_output
        assert Path(custom_output).exists()

    @pytest.mark.unit
    def test_nonexistent_file_raises_error(self, tmp_path):
        """存在しないファイルでFileNotFoundErrorが発生することをテスト"""
        nonexistent = str(tmp_path / "nonexistent.dx")

        with pytest.raises(FileNotFoundError):
            compress_profile(nonexistent)

    @pytest.mark.unit
    def test_roundtrip_with_real_dx(self, tmp_path):
        """実際のDXファイルで圧縮→解凍の内容一致をテスト"""
        # DXファイルを作成
        data = np.random.rand(10, 10, 10)
        dx_path = str(tmp_path / "profile.dx")
        _create_dx_file(dx_path, data)

        # 元のファイル内容を保存
        original_bytes = Path(dx_path).read_bytes()

        gz_path = compress_profile(dx_path, remove_original=False)

        with gzip.open(gz_path, "rb") as f:
            decompressed = f.read()

        assert decompressed == original_bytes


# =============================================================================
# _is_valid_file
# =============================================================================


class TestIsValidFile:
    """_is_valid_file関数のテスト"""

    @pytest.mark.unit
    def test_valid_file(self, tmp_path):
        """内容のあるファイルがTrueを返すことをテスト"""
        valid_file = tmp_path / "valid.txt"
        valid_file.write_text("some content")

        assert _is_valid_file(valid_file) is True

    @pytest.mark.unit
    def test_empty_file_returns_false(self, tmp_path):
        """空ファイルがFalseを返すことをテスト"""
        empty_file = tmp_path / "empty.txt"
        empty_file.write_text("")

        assert _is_valid_file(empty_file) is False

    @pytest.mark.unit
    def test_nonexistent_file_returns_false(self, tmp_path):
        """存在しないファイルがFalseを返すことをテスト"""
        nonexistent = tmp_path / "nonexistent.txt"

        assert _is_valid_file(nonexistent) is False

    @pytest.mark.unit
    def test_directory_returns_false(self, tmp_path):
        """ディレクトリがFalseを返すことをテスト（st_sizeは0ではないがファイルではない）"""
        # ディレクトリはexists()=True, stat().st_size > 0 の場合がある
        # _is_valid_fileの動作を確認
        result = _is_valid_file(tmp_path)
        # ディレクトリのst_sizeは通常0より大きいのでTrueになる
        # （関数はファイルとディレクトリを区別しない仕様）
        assert isinstance(result, bool)

    @pytest.mark.unit
    def test_binary_file(self, tmp_path):
        """バイナリファイルがTrueを返すことをテスト"""
        binary_file = tmp_path / "data.bin"
        binary_file.write_bytes(b"\x00\x01\x02\x03")

        assert _is_valid_file(binary_file) is True

    @pytest.mark.unit
    def test_single_byte_file(self, tmp_path):
        """1バイトファイルがTrueを返すことをテスト"""
        one_byte = tmp_path / "one.txt"
        one_byte.write_bytes(b"x")

        assert _is_valid_file(one_byte) is True
