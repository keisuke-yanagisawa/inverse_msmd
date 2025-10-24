#!/usr/bin/env python3
"""
スコア付き可視化機能のユニットテスト
"""

import pytest
import numpy as np
from pathlib import Path
from rdkit import Chem
from inverse_msmd.utils.visualization_utils import draw_scored_molecules_grid


class TestDrawScoredMoleculesGrid:
    """draw_scored_molecules_grid関数のテスト"""
    
    def test_basic_visualization(self, tmp_path):
        """基本的な可視化のテスト"""
        # テスト用の分子を作成
        molecules = [
            Chem.MolFromSmiles("CCO"),
            Chem.MolFromSmiles("CCN"),
            Chem.MolFromSmiles("CCC")
        ]
        scores = [-123.45, -145.67, -156.89]
        output_path = tmp_path / "test_scored_grid.png"
        
        # 関数を実行
        draw_scored_molecules_grid(
            molecules,
            scores,
            str(output_path),
            align_molecules=False
        )
        
        # 出力ファイルが生成されたことを確認
        assert output_path.exists()
        assert output_path.stat().st_size > 0
    
    def test_with_alignment(self, tmp_path):
        """アライメント付き可視化のテスト"""
        # ベンゼン環を含む分子を作成（アライメントが効果的）
        molecules = [
            Chem.MolFromSmiles("c1ccccc1CCO"),
            Chem.MolFromSmiles("c1ccccc1CCN"),
            Chem.MolFromSmiles("c1ccccc1CCC")
        ]
        scores = [-100.0, -120.0, -140.0]
        output_path = tmp_path / "test_aligned_grid.png"
        
        # 関数を実行（アライメント有効）
        draw_scored_molecules_grid(
            molecules,
            scores,
            str(output_path),
            align_molecules=True
        )
        
        # 出力ファイルが生成されたことを確認
        assert output_path.exists()
        assert output_path.stat().st_size > 0
    
    def test_with_titles(self, tmp_path):
        """タイトル付き可視化のテスト"""
        molecules = [
            Chem.MolFromSmiles("CCO"),
            Chem.MolFromSmiles("CCN")
        ]
        scores = [-123.45, -145.67]
        titles = ["Pattern 0", "Pattern 1"]
        output_path = tmp_path / "test_titled_grid.png"
        
        # 関数を実行
        draw_scored_molecules_grid(
            molecules,
            scores,
            str(output_path),
            titles=titles
        )
        
        # 出力ファイルが生成されたことを確認
        assert output_path.exists()
    
    def test_single_molecule(self, tmp_path):
        """単一分子の可視化のテスト"""
        molecules = [Chem.MolFromSmiles("CCO")]
        scores = [-123.45]
        output_path = tmp_path / "test_single.png"
        
        # 関数を実行
        draw_scored_molecules_grid(
            molecules,
            scores,
            str(output_path)
        )
        
        # 出力ファイルが生成されたことを確認
        assert output_path.exists()
    
    def test_many_molecules(self, tmp_path):
        """多数の分子の可視化のテスト（グリッドレイアウトの確認）"""
        # 8個の分子を作成
        smiles_list = ["CCO", "CCN", "CCC", "CCCO", "CCCN", "CCCC", "CC(C)O", "CC(C)N"]
        molecules = [Chem.MolFromSmiles(smi) for smi in smiles_list]
        scores = [-100.0 - i * 10 for i in range(len(molecules))]
        output_path = tmp_path / "test_many_molecules.png"
        
        # 関数を実行（最大4列）
        draw_scored_molecules_grid(
            molecules,
            scores,
            str(output_path),
            max_cols=4
        )
        
        # 出力ファイルが生成されたことを確認
        assert output_path.exists()
    
    def test_empty_molecules_raises_error(self):
        """空の分子リストでエラーが発生することを確認"""
        with pytest.raises(ValueError, match="molecules リストは空にできません"):
            draw_scored_molecules_grid([], [], "output.png")
    
    def test_mismatched_lengths_raises_error(self):
        """分子とスコアの数が一致しない場合にエラーが発生することを確認"""
        molecules = [Chem.MolFromSmiles("CCO"), Chem.MolFromSmiles("CCN")]
        scores = [-123.45]  # 1つだけ
        
        with pytest.raises(ValueError, match="moleculesの数.*scoresの数.*一致しません"):
            draw_scored_molecules_grid(molecules, scores, "output.png")
    
    def test_custom_image_size(self, tmp_path):
        """カスタム画像サイズのテスト"""
        molecules = [Chem.MolFromSmiles("CCO")]
        scores = [-123.45]
        output_path = tmp_path / "test_custom_size.png"
        
        # カスタムサイズで実行
        draw_scored_molecules_grid(
            molecules,
            scores,
            str(output_path),
            image_size=(600, 600),
            dpi=200
        )
        
        # 出力ファイルが生成されたことを確認
        assert output_path.exists()
    
    def test_output_directory_creation(self, tmp_path):
        """出力ディレクトリが自動作成されることを確認"""
        molecules = [Chem.MolFromSmiles("CCO")]
        scores = [-123.45]
        output_path = tmp_path / "nested" / "dir" / "test.png"
        
        # ディレクトリが存在しないことを確認
        assert not output_path.parent.exists()
        
        # 関数を実行
        draw_scored_molecules_grid(molecules, scores, str(output_path))
        
        # ディレクトリと出力ファイルが生成されたことを確認
        assert output_path.parent.exists()
        assert output_path.exists()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])