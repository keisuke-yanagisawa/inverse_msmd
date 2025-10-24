"""
プロファイルスコア計算統合テスト

このテストモジュールは、プロファイルスコア計算機能と
統合ワークフローの連携をテストします。
"""

import pytest
import numpy as np
from pathlib import Path


class TestIntegratedWorkflowWithoutScoring:
    """スコア計算なしの統合ワークフロー（後方互換性）テスト"""
    
    def test_backward_compatibility(self, tmp_path):
        """profile_dir=Noneで既存の動作が維持される"""
        from inverse_msmd.substructure_replacement import integrated_substructure_replacement
        
        # テストデータの存在確認
        ligand_file = "data/atom_matching/4hw3_A_lig.sdf"
        protein_file = "data/sample_proteins/4hw3_A.pdb"
        from_file = "data/sample_probes/E23"
        to_file = "data/sample_probes/E24"
        
        if not all(Path(f).exists() for f in [ligand_file, protein_file, 
                                                f"{from_file}.pdb", f"{from_file}.smi",
                                                f"{to_file}.pdb", f"{to_file}.smi"]):
            pytest.skip("テストデータが見つかりません")
        
        output_dir = tmp_path / "backward_compat"
        
        results = integrated_substructure_replacement(
            ligand_file=ligand_file,
            protein_file=protein_file,
            from_file=from_file,
            to_file=to_file,
            output_dir=str(output_dir),
            match_index=0,
            profile_dir=None  # スコア計算なし
        )
        
        assert len(results) > 0
        # スコアが含まれていないことを確認
        assert 'score' not in results[0]
        # 必須キーが含まれていることを確認
        assert 'ligand_file' in results[0]
        assert 'protein_file' in results[0]
        assert 'pattern_index' in results[0]
        
        # ファイルが実際に生成されていることを確認
        for result in results:
            assert Path(result['ligand_file']).exists()
            assert Path(result['protein_file']).exists()


class TestIntegratedWorkflowWithScoring:
    """スコア計算ありの統合ワークフローテスト"""
    
    def test_with_profile_scoring(self, tmp_path):
        """スコア計算が正しく統合される"""
        from inverse_msmd.substructure_replacement import integrated_substructure_replacement
        
        # テストデータの存在確認
        ligand_file = "data/atom_matching/4hw3_A_lig.sdf"
        protein_file = "data/sample_proteins/4hw3_A.pdb"
        from_file = "data/sample_probes/E23"
        to_file = "data/sample_probes/E24"
        profile_dir = "data/profiles/"
        
        if not all(Path(f).exists() for f in [ligand_file, protein_file,
                                                f"{from_file}.pdb", f"{from_file}.smi",
                                                f"{to_file}.pdb", f"{to_file}.smi"]):
            pytest.skip("テストデータが見つかりません")
        
        if not Path(profile_dir).exists():
            pytest.skip(f"プロファイルディレクトリが見つかりません: {profile_dir}")
        
        output_dir = tmp_path / "with_score"
        
        results = integrated_substructure_replacement(
            ligand_file=ligand_file,
            protein_file=protein_file,
            from_file=from_file,
            to_file=to_file,
            output_dir=str(output_dir),
            match_index=0,
            profile_dir=profile_dir,
            probe_id="E24"
        )
        
        assert len(results) > 0
        
        # 全パターンにスコアが付与されていることを確認
        for result in results:
            assert 'score' in result
            assert isinstance(result['score'], float)
            assert not np.isnan(result['score'])
            assert not np.isinf(result['score'])
            
            # その他のキーも確認
            assert 'ligand_file' in result
            assert 'protein_file' in result
            assert 'pattern_index' in result
            
            # ファイルが生成されていることを確認
            assert Path(result['ligand_file']).exists()
            assert Path(result['protein_file']).exists()
    
    def test_score_sorting(self, tmp_path):
        """結果がスコアで降順ソートされる"""
        from inverse_msmd.substructure_replacement import integrated_substructure_replacement
        
        # テストデータの存在確認
        ligand_file = "data/atom_matching/4hw3_A_lig.sdf"
        protein_file = "data/sample_proteins/4hw3_A.pdb"
        from_file = "data/sample_probes/E23"
        to_file = "data/sample_probes/E24"
        profile_dir = "data/profiles/"
        
        if not all(Path(f).exists() for f in [ligand_file, protein_file,
                                                f"{from_file}.pdb", f"{from_file}.smi",
                                                f"{to_file}.pdb", f"{to_file}.smi"]):
            pytest.skip("テストデータが見つかりません")
        
        if not Path(profile_dir).exists():
            pytest.skip(f"プロファイルディレクトリが見つかりません: {profile_dir}")
        
        output_dir = tmp_path / "sorted"
        
        results = integrated_substructure_replacement(
            ligand_file=ligand_file,
            protein_file=protein_file,
            from_file=from_file,
            to_file=to_file,
            output_dir=str(output_dir),
            match_index=0,
            profile_dir=profile_dir,
            probe_id="E24"
        )
        
        # 複数パターンがある場合、スコアで降順ソートされているか確認
        if len(results) > 1:
            scores = [r['score'] for r in results]
            # 降順ソートされていることを確認
            assert scores == sorted(scores, reverse=True)
    

class TestParameterValidation:
    """パラメータバリデーションのテスト"""
    
    def test_missing_probe_id(self, tmp_path):
        """profile_dir指定時にprobe_idがないとエラー"""
        from inverse_msmd.substructure_replacement import integrated_substructure_replacement
        
        # テストデータの存在確認
        ligand_file = "data/atom_matching/4hw3_A_lig.sdf"
        protein_file = "data/sample_proteins/4hw3_A.pdb"
        from_file = "data/sample_probes/E23"
        to_file = "data/sample_probes/E24"
        profile_dir = "data/profiles/"
        
        if not all(Path(f).exists() for f in [ligand_file, protein_file,
                                                f"{from_file}.pdb", f"{from_file}.smi",
                                                f"{to_file}.pdb", f"{to_file}.smi"]):
            pytest.skip("テストデータが見つかりません")
        
        if not Path(profile_dir).exists():
            pytest.skip(f"プロファイルディレクトリが見つかりません: {profile_dir}")
        
        output_dir = tmp_path / "error_test"
        
        with pytest.raises(ValueError, match="probe_idも必須です"):
            integrated_substructure_replacement(
                ligand_file=ligand_file,
                protein_file=protein_file,
                from_file=from_file,
                to_file=to_file,
                output_dir=str(output_dir),
                match_index=0,
                profile_dir=profile_dir,
                probe_id=None  # エラーを引き起こす
            )


class TestEndToEndWorkflow:
    """エンドツーエンドの完全なワークフローテスト"""
    
    def test_complete_workflow(self, tmp_path):
        """完全なワークフローが正常に動作する"""
        from inverse_msmd.substructure_replacement import integrated_substructure_replacement
        
        # テストデータの存在確認
        ligand_file = "data/atom_matching/4hw3_A_lig.sdf"
        protein_file = "data/sample_proteins/4hw3_A.pdb"
        from_file = "data/sample_probes/E23"
        to_file = "data/sample_probes/E24"
        profile_dir = "data/profiles/"
        
        if not all(Path(f).exists() for f in [ligand_file, protein_file,
                                                f"{from_file}.pdb", f"{from_file}.smi",
                                                f"{to_file}.pdb", f"{to_file}.smi"]):
            pytest.skip("テストデータが見つかりません")
        
        if not Path(profile_dir).exists():
            pytest.skip(f"プロファイルディレクトリが見つかりません: {profile_dir}")
        
        output_dir = tmp_path / "e2e"
        
        # 1. 統合ワークフローの実行
        results = integrated_substructure_replacement(
            ligand_file=ligand_file,
            protein_file=protein_file,
            from_file=from_file,
            to_file=to_file,
            output_dir=str(output_dir),
            match_index=0,
            profile_dir=profile_dir,
            probe_id="E24"
        )
        
        # 2. 結果の検証
        assert len(results) > 0
        assert all('score' in r for r in results)
        
        # 3. ファイルの存在確認
        for result in results:
            ligand_path = Path(result['ligand_file'])
            protein_path = Path(result['protein_file'])
            
            assert ligand_path.exists()
            assert protein_path.exists()
            
            # ファイルサイズが0より大きいことを確認
            assert ligand_path.stat().st_size > 0
            assert protein_path.stat().st_size > 0
        
        # 4. スコアの範囲確認
        scores = [r['score'] for r in results]
        assert all(isinstance(s, float) for s in scores)
        # スコアが有効な値であることを確認（NaNやInfでない）
        assert all(not np.isnan(s) and not np.isinf(s) for s in scores)
        
        # 5. スコアで降順ソートされていることを確認
        if len(scores) > 1:
            assert scores == sorted(scores, reverse=True)