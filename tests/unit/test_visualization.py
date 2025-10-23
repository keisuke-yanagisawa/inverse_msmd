"""可視化機能の単体テスト"""

import pytest
from pathlib import Path
from inverse_msmd.substructure_replacement import (
    find_substructure_in_ligand,
    visualize_multiple_matches
)


class TestVisualization:
    """可視化機能のテスト"""

    @pytest.mark.unit
    @pytest.mark.visual
    def test_png_generation(self, ligand_mol, e23_mol, output_dir):
        """PNG画像が正常に生成されることを確認"""
        # 部分構造探索
        matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        
        # PNG画像を生成
        output_path = output_dir / "substructure_matches_test.png"
        visualize_multiple_matches(
            ligand_mol,
            e23_mol,
            matches,
            str(output_path)
        )
        
        # ファイルの存在確認
        assert output_path.exists(), "PNG画像が生成されていません"
        assert output_path.stat().st_size > 0, "PNG画像のサイズが0です"

    @pytest.mark.unit
    @pytest.mark.visual
    def test_multiple_matches_visualization(self, ligand_mol, a01_mol, output_dir):
        """複数マッチの可視化が正常に動作することを確認"""
        matches = find_substructure_in_ligand(ligand_mol, a01_mol)
        
        # 複数マッチがある場合のみテスト
        if len(matches) > 1:
            output_path = output_dir / "multiple_matches_test.png"
            visualize_multiple_matches(
                ligand_mol,
                a01_mol,
                matches,
                str(output_path)
            )
            
            assert output_path.exists(), "複数マッチのPNG画像が生成されていません"
            assert output_path.stat().st_size > 0, "PNG画像のサイズが0です"
        else:
            pytest.skip("このテストには複数マッチが必要です")

    @pytest.mark.unit
    @pytest.mark.visual
    def test_single_match_visualization(self, ligand_mol, e23_mol, output_dir):
        """単一マッチの可視化が正常に動作することを確認"""
        matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        
        if len(matches) > 0:
            # 最初のマッチのみを使用
            single_match = [matches[0]]
            output_path = output_dir / "single_match_test.png"
            
            visualize_multiple_matches(
                ligand_mol,
                e23_mol,
                single_match,
                str(output_path)
            )
            
            assert output_path.exists(), "単一マッチのPNG画像が生成されていません"
            assert output_path.stat().st_size > 0, "PNG画像のサイズが0です"