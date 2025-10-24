"""
visualization_utils モジュールのテスト
"""

import pytest
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from PIL import Image

from inverse_msmd.utils.visualization_utils import (
    draw_molecule_comparison,
    draw_molecule_grid,
    draw_molecule_with_highlights
)


@pytest.fixture
def simple_molecules():
    """テスト用の簡単な分子を作成"""
    mol1 = Chem.MolFromSmiles("CCO")  # エタノール
    mol2 = Chem.MolFromSmiles("CCN")  # エチルアミン
    mol3 = Chem.MolFromSmiles("CCC")  # プロパン
    return [mol1, mol2, mol3]


@pytest.fixture
def benzene_molecule():
    """ベンゼン分子を作成"""
    return Chem.MolFromSmiles("c1ccccc1")


@pytest.fixture
def temp_output_dir(tmp_path):
    """一時的な出力ディレクトリを作成"""
    output_dir = tmp_path / "test_output"
    output_dir.mkdir()
    return output_dir


class TestDrawMoleculeComparison:
    """draw_molecule_comparison関数のテスト"""
    
    def test_basic_comparison(self, simple_molecules, temp_output_dir):
        """基本的な分子比較描画のテスト"""
        output_path = temp_output_dir / "comparison.png"
        
        draw_molecule_comparison(
            simple_molecules[:2],
            output_path
        )
        
        # ファイルが作成されたことを確認
        assert output_path.exists()
        
        # 画像として読み込めることを確認
        img = Image.open(output_path)
        assert img.size[0] > 0
        assert img.size[1] > 0
    
    def test_comparison_with_titles(self, simple_molecules, temp_output_dir):
        """タイトル付き分子比較描画のテスト"""
        output_path = temp_output_dir / "comparison_with_titles.png"
        titles = ["エタノール", "エチルアミン"]
        
        draw_molecule_comparison(
            simple_molecules[:2],
            output_path,
            titles=titles
        )
        
        assert output_path.exists()
    
    def test_comparison_with_highlights(self, simple_molecules, temp_output_dir):
        """ハイライト付き分子比較描画のテスト"""
        output_path = temp_output_dir / "comparison_with_highlights.png"
        highlight_atoms_list = [[0, 1], [0, 1, 2]]
        
        draw_molecule_comparison(
            simple_molecules[:2],
            output_path,
            highlight_atoms_list=highlight_atoms_list
        )
        
        assert output_path.exists()
    
    def test_empty_molecules_raises_error(self, temp_output_dir):
        """空の分子リストでエラーが発生することを確認"""
        output_path = temp_output_dir / "empty.png"
        
        with pytest.raises(ValueError, match="molecules リストは空にできません"):
            draw_molecule_comparison([], output_path)
    
    def test_mismatched_titles_raises_error(self, simple_molecules, temp_output_dir):
        """タイトル数が一致しない場合のエラーを確認"""
        output_path = temp_output_dir / "mismatched.png"
        
        with pytest.raises(ValueError, match="titlesの数"):
            draw_molecule_comparison(
                simple_molecules[:2],
                output_path,
                titles=["タイトル1"]  # 2つの分子に対して1つのタイトル
            )
    
    def test_mismatched_highlights_raises_error(self, simple_molecules, temp_output_dir):
        """ハイライトリスト数が一致しない場合のエラーを確認"""
        output_path = temp_output_dir / "mismatched_highlights.png"
        
        with pytest.raises(ValueError, match="highlight_atoms_listの数"):
            draw_molecule_comparison(
                simple_molecules[:2],
                output_path,
                highlight_atoms_list=[[0, 1]]  # 2つの分子に対して1つのリスト
            )
    
    def test_single_molecule(self, simple_molecules, temp_output_dir):
        """単一分子の描画テスト"""
        output_path = temp_output_dir / "single.png"
        
        draw_molecule_comparison(
            simple_molecules[:1],
            output_path
        )
        
        assert output_path.exists()
    
    def test_auto_create_directory(self, simple_molecules, tmp_path):
        """出力ディレクトリの自動作成をテスト"""
        output_path = tmp_path / "new_dir" / "subdir" / "comparison.png"
        
        draw_molecule_comparison(
            simple_molecules[:2],
            output_path
        )
        
        assert output_path.exists()
        assert output_path.parent.exists()


class TestDrawMoleculeGrid:
    """draw_molecule_grid関数のテスト"""
    
    def test_basic_grid(self, benzene_molecule, temp_output_dir):
        """基本的なグリッド描画のテスト"""
        output_path = temp_output_dir / "grid.png"
        
        # ベンゼン環全体をマッチとして使用
        matches = [(0, 1, 2, 3, 4, 5)]
        
        draw_molecule_grid(
            benzene_molecule,
            matches,
            output_path
        )
        
        assert output_path.exists()
    
    def test_multiple_matches(self, temp_output_dir):
        """複数マッチのグリッド描画テスト"""
        # トルエン（メチルベンゼン）
        mol = Chem.MolFromSmiles("Cc1ccccc1")
        # ベンゼン環をサブストラクチャとして検索
        pattern = Chem.MolFromSmiles("c1ccccc1")
        matches = mol.GetSubstructMatches(pattern)
        
        output_path = temp_output_dir / "multiple_matches.png"
        
        draw_molecule_grid(
            mol,
            matches,
            output_path,
            title="トルエンのベンゼン環"
        )
        
        assert output_path.exists()
    
    def test_many_matches_grid_layout(self, temp_output_dir):
        """多数のマッチでグリッドレイアウトをテスト"""
        mol = Chem.MolFromSmiles("c1ccccc1")  # ベンゼン
        # ベンゼンの有効な原子インデックスでマッチを作成（0-5の範囲内）
        matches = [(i, (i+1)%6, (i+2)%6) for i in range(6)]
        
        output_path = temp_output_dir / "many_matches.png"
        
        draw_molecule_grid(
            mol,
            matches,
            output_path,
            max_cols=3
        )
        
        assert output_path.exists()
    
    def test_empty_matches_raises_error(self, benzene_molecule, temp_output_dir):
        """空のマッチリストでエラーが発生することを確認"""
        output_path = temp_output_dir / "empty_matches.png"
        
        with pytest.raises(ValueError, match="matches リストは空にできません"):
            draw_molecule_grid(benzene_molecule, [], output_path)
    
    def test_custom_image_size(self, benzene_molecule, temp_output_dir):
        """カスタム画像サイズのテスト"""
        output_path = temp_output_dir / "custom_size.png"
        matches = [(0, 1, 2)]
        
        draw_molecule_grid(
            benzene_molecule,
            matches,
            output_path,
            image_size=(600, 600),
            dpi=200
        )
        
        assert output_path.exists()


class TestDrawMoleculeWithHighlights:
    """draw_molecule_with_highlights関数のテスト"""
    
    def test_basic_highlight(self, benzene_molecule):
        """基本的なハイライト描画のテスト"""
        img = draw_molecule_with_highlights(benzene_molecule)
        
        # PIL Imageオブジェクトが返されることを確認
        assert isinstance(img, Image.Image)
        assert img.size == (400, 400)
    
    def test_with_highlight_atoms(self, benzene_molecule):
        """原子ハイライト付き描画のテスト"""
        img = draw_molecule_with_highlights(
            benzene_molecule,
            highlight_atoms=[0, 1, 2]
        )
        
        assert isinstance(img, Image.Image)
    
    def test_custom_size(self, benzene_molecule):
        """カスタムサイズのテスト"""
        img = draw_molecule_with_highlights(
            benzene_molecule,
            size=(600, 600)
        )
        
        assert isinstance(img, Image.Image)
        assert img.size == (600, 600)
    
    def test_save_image(self, benzene_molecule, temp_output_dir):
        """画像保存のテスト"""
        img = draw_molecule_with_highlights(benzene_molecule)
        output_path = temp_output_dir / "highlighted.png"
        
        img.save(output_path)
        
        assert output_path.exists()


class TestIntegration:
    """統合テスト"""
    
    def test_workflow_comparison_and_grid(self, simple_molecules, temp_output_dir):
        """比較描画とグリッド描画を組み合わせたワークフローテスト"""
        # 比較描画
        comparison_path = temp_output_dir / "workflow_comparison.png"
        draw_molecule_comparison(
            simple_molecules[:2],
            comparison_path,
            titles=["分子1", "分子2"]
        )
        
        # グリッド描画
        grid_path = temp_output_dir / "workflow_grid.png"
        mol = simple_molecules[0]
        # 単純なマッチを作成
        matches = [(0, 1), (1, 2)]
        draw_molecule_grid(
            mol,
            matches,
            grid_path
        )
        
        assert comparison_path.exists()
        assert grid_path.exists()