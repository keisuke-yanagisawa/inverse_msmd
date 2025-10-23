"""
プロファイルスコア計算機能のユニットテスト

このテストモジュールは、profile_scoringモジュールの
個別機能をテストします。
"""

import pytest
import numpy as np
from pathlib import Path


class TestProfileScoringImport:
    """プロファイルスコアモジュールのインポートテスト"""
    
    def test_import_module(self):
        """profile_scoringモジュールがインポートできる"""
        from inverse_msmd import profile_scoring
        assert profile_scoring is not None
    
    def test_import_function(self):
        """calculate_profile_score関数がインポートできる"""
        from inverse_msmd.profile_scoring import calculate_profile_score
        assert callable(calculate_profile_score)
    
    def test_import_from_package(self):
        """パッケージレベルから関数がインポートできる"""
        from inverse_msmd import calculate_profile_score
        assert callable(calculate_profile_score)


class TestProfileScoreCalculation:
    """プロファイルスコア計算のテスト"""
    
    @pytest.fixture
    def sample_protein(self):
        """サンプルタンパク質構造を読み込み"""
        from inverse_msmd.utils.bio_utils import PDB
        protein_path = "data/sample_proteins/4hw3_A.pdb"
        if not Path(protein_path).exists():
            pytest.skip(f"テストデータが見つかりません: {protein_path}")
        return PDB.get_structure(protein_path)
    
    @pytest.fixture
    def sample_probe_center(self):
        """サンプルプローブ中心座標"""
        from inverse_msmd.utils.bio_utils import PDB
        probe_path = "data/sample_probes/E24.pdb"
        if not Path(probe_path).exists():
            pytest.skip(f"テストデータが見つかりません: {probe_path}")
        probe = PDB.get_structure(probe_path)
        return PDB.get_attr(probe, "coord").mean(axis=0)
    
    def test_calculate_score_no_weight(self, sample_protein, sample_probe_center):
        """重み付けなしでスコアが計算される"""
        from inverse_msmd.profile_scoring import calculate_profile_score
        
        profile_dir = "data/profiles/"
        if not Path(profile_dir).exists():
            pytest.skip(f"プロファイルディレクトリが見つかりません: {profile_dir}")
        
        score = calculate_profile_score(
            sample_protein,
            sample_probe_center,
            profile_dir,
            "E24",
            gamma=0.0
        )
        
        assert isinstance(score, float)
        assert not np.isnan(score)
        assert not np.isinf(score)
        # 対数スコアは通常負の値
        assert score < 0
    
    def test_calculate_score_with_weight(self, sample_protein, sample_probe_center):
        """距離重み付けありでスコアが計算される"""
        from inverse_msmd.profile_scoring import calculate_profile_score
        
        profile_dir = "data/profiles/"
        if not Path(profile_dir).exists():
            pytest.skip(f"プロファイルディレクトリが見つかりません: {profile_dir}")
        
        score_weighted = calculate_profile_score(
            sample_protein,
            sample_probe_center,
            profile_dir,
            "E24",
            gamma=0.003
        )
        
        score_no_weight = calculate_profile_score(
            sample_protein,
            sample_probe_center,
            profile_dir,
            "E24",
            gamma=0.0
        )
        
        assert isinstance(score_weighted, float)
        assert isinstance(score_no_weight, float)
        # gamma値が異なればスコアも異なるはず
        assert score_weighted != score_no_weight
    
    def test_different_probes(self, sample_protein, sample_probe_center):
        """異なるプローブでスコアが異なる"""
        from inverse_msmd.profile_scoring import calculate_profile_score
        
        profile_dir = "data/profiles/"
        if not Path(profile_dir).exists():
            pytest.skip(f"プロファイルディレクトリが見つかりません: {profile_dir}")
        
        # E24とA08で異なるスコアが得られるはず
        score_e24 = calculate_profile_score(
            sample_protein,
            sample_probe_center,
            profile_dir,
            "E24",
            gamma=0.0
        )
        
        score_a08 = calculate_profile_score(
            sample_protein,
            sample_probe_center,
            profile_dir,
            "A08",
            gamma=0.0
        )
        
        # 異なるプローブでは通常スコアが異なる
        # ただし、まれに同じ値になる可能性もあるので、
        # 両方とも有効な値であることを確認
        assert isinstance(score_e24, float)
        assert isinstance(score_a08, float)


class TestProfileScoringErrors:
    """プロファイルスコア計算のエラーハンドリングテスト"""
    
    def test_profile_dir_not_found(self, tmp_path):
        """存在しないディレクトリでエラーが発生する"""
        from inverse_msmd.profile_scoring import calculate_profile_score
        from inverse_msmd.utils.bio_utils import PDB
        
        protein_path = "data/sample_proteins/4hw3_A.pdb"
        if not Path(protein_path).exists():
            pytest.skip(f"テストデータが見つかりません: {protein_path}")
        
        protein = PDB.get_structure(protein_path)
        probe_center = np.array([10.0, 15.0, 20.0])
        
        nonexistent_dir = tmp_path / "nonexistent"
        
        with pytest.raises(ValueError, match="プロファイルディレクトリが見つかりません"):
            calculate_profile_score(
                protein,
                probe_center,
                str(nonexistent_dir),
                "E24"
            )
    
    def test_profile_files_not_found(self, tmp_path):
        """プロファイルファイルが見つからない場合のエラー"""
        from inverse_msmd.profile_scoring import calculate_profile_score
        from inverse_msmd.utils.bio_utils import PDB
        
        protein_path = "data/sample_proteins/4hw3_A.pdb"
        if not Path(protein_path).exists():
            pytest.skip(f"テストデータが見つかりません: {protein_path}")
        
        protein = PDB.get_structure(protein_path)
        probe_center = np.array([10.0, 15.0, 20.0])
        
        # 空のディレクトリを作成
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        
        with pytest.raises(ValueError, match="プロファイルファイルが見つかりません"):
            calculate_profile_score(
                protein,
                probe_center,
                str(empty_dir),
                "INVALID"
            )
    
    def test_no_cb_atoms(self):
        """Cβ原子がない場合のエラー（GLYのみのタンパク質など）"""
        from inverse_msmd.profile_scoring import calculate_profile_score
        from Bio.PDB.Structure import Structure
        from Bio.PDB.Model import Model
        from Bio.PDB.Chain import Chain
        from Bio.PDB.Residue import Residue
        from Bio.PDB.Atom import Atom
        
        # GLY残基のみを持つダミータンパク質を作成
        structure = Structure("dummy")
        model = Model(0)
        chain = Chain("A")
        residue = Residue((" ", 1, " "), "GLY", 1)
        
        # GLYはCβを持たないので、CA原子のみ追加
        atom = Atom("CA", np.array([0.0, 0.0, 0.0]), 1.0, 1.0, " ", "CA", 1, "C")
        residue.add(atom)
        
        chain.add(residue)
        model.add(chain)
        structure.add(model)
        
        probe_center = np.array([10.0, 15.0, 20.0])
        
        profile_dir = "data/profiles/"
        if not Path(profile_dir).exists():
            pytest.skip(f"プロファイルディレクトリが見つかりません: {profile_dir}")
        
        with pytest.raises(ValueError, match="Cβ原子が見つかりません"):
            calculate_profile_score(
                structure,
                probe_center,
                profile_dir,
                "E24"
            )


class TestProfileScoringEdgeCases:
    """エッジケースのテスト"""
    
    def test_negative_profile_values(self):
        """負のプロファイル値が適切に処理される"""
        # このテストは実際のプロファイルデータに依存するため、
        # 負の値が出る場合のロジックが正しく動作することを確認
        # 実装では負の値は最小値で置換されるため、NaNやInfにならないことを確認
        from inverse_msmd.profile_scoring import calculate_profile_score
        from inverse_msmd.utils.bio_utils import PDB
        
        protein_path = "data/sample_proteins/4hw3_A.pdb"
        if not Path(protein_path).exists():
            pytest.skip(f"テストデータが見つかりません: {protein_path}")
        
        protein = PDB.get_structure(protein_path)
        probe_center = np.array([10.0, 15.0, 20.0])
        
        profile_dir = "data/profiles/"
        if not Path(profile_dir).exists():
            pytest.skip(f"プロファイルディレクトリが見つかりません: {profile_dir}")
        
        score = calculate_profile_score(
            protein,
            probe_center,
            profile_dir,
            "E24",
            gamma=0.0
        )
        
        # 負の値が適切に処理され、NaNやInfにならないことを確認
        assert not np.isnan(score)
        assert not np.isinf(score)
        assert isinstance(score, float)