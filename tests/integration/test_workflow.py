"""統合ワークフローのテスト"""

import pytest
from pathlib import Path
from rdkit import Chem
import numpy as np
from inverse_msmd.substructure_replacement import (
    find_substructure_in_ligand,
    match_substructures,
    calculate_transformation,
    visualize_multiple_matches,
    integrated_substructure_replacement
)
from inverse_msmd.utils.bio_utils import PDB


class TestIntegratedWorkflow:
    """統合ワークフローのテスト"""

    @pytest.mark.integration
    def test_search_and_match_workflow(self, ligand_mol, e23_mol, e24_mol):
        """部分構造探索からatom matchingまでの一連の流れをテスト"""
        # ステップ1: 部分構造探索
        matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        assert len(matches) > 0, "部分構造が見つかりません"
        
        # ステップ2: Atom matching
        atom_pair_patterns = match_substructures(e23_mol, e24_mol)
        assert len(atom_pair_patterns) > 0, "Atom matchingパターンが見つかりません"
        
        # ステップ3: 座標の取得と検証
        ligand_coords = ligand_mol.GetConformer().GetPositions()
        e23_coords = e23_mol.GetConformer().GetPositions()
        e24_coords = e24_mol.GetConformer().GetPositions()
        
        assert ligand_coords.shape[0] == ligand_mol.GetNumAtoms()
        assert e23_coords.shape[0] == e23_mol.GetNumAtoms()
        assert e24_coords.shape[0] == e24_mol.GetNumAtoms()

    @pytest.mark.integration
    def test_transformation_calculation(self, ligand_mol, e23_mol, e24_mol):
        """変換行列計算までの一連の流れをテスト"""
        # 部分構造探索
        matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        match = matches[0]
        
        # Atom matching
        atom_pair_patterns = match_substructures(e23_mol, e24_mol)
        atom_pairs = atom_pair_patterns[0]
        
        # 座標取得
        ligand_coords = ligand_mol.GetConformer().GetPositions()
        e24_coords = e24_mol.GetConformer().GetPositions()
        ligand_e23_coords = ligand_coords[list(match)]
        
        # 変換行列計算
        rot, tran = calculate_transformation(
            ligand_e23_coords,
            e24_coords,
            atom_pairs
        )
        
        # 回転行列の検証
        assert rot.shape == (3, 3), "回転行列の形状が不正です"
        assert tran.shape == (3,), "並進ベクトルの形状が不正です"
        
        # 回転行列が直交行列であることを確認
        rot_test = np.dot(rot, rot.T)
        identity = np.eye(3)
        assert np.allclose(rot_test, identity, atol=1e-6), \
            "回転行列が直交行列ではありません"
        
        # 行列式が1に近いことを確認
        det = np.linalg.det(rot)
        assert np.abs(det - 1.0) < 1e-6, \
            f"回転行列の行列式が1ではありません: {det}"

    @pytest.mark.integration
    @pytest.mark.slow
    def test_visualization_with_search(self, ligand_mol, e23_mol, output_dir):
        """部分構造探索と可視化の統合テスト"""
        # 部分構造探索
        matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        
        # 可視化
        output_path = output_dir / "integrated_visualization.png"
        visualize_multiple_matches(
            ligand_mol,
            e23_mol,
            matches,
            str(output_path)
        )
        
        # 検証
        assert output_path.exists(), "可視化画像が生成されていません"
        assert output_path.stat().st_size > 0, "画像ファイルが空です"

    @pytest.mark.integration
    def test_multiple_probe_workflow(self, ligand_mol, e23_mol, e24_mol, a01_mol):
        """複数のプローブを使用したワークフローをテスト"""
        # E23での処理
        e23_matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        assert len(e23_matches) > 0, "E23のマッチが見つかりません"
        
        # A01での処理
        a01_matches = find_substructure_in_ligand(ligand_mol, a01_mol)
        assert len(a01_matches) > 0, "A01のマッチが見つかりません"
        
        # E23-E24のatom matching
        e23_e24_pairs = match_substructures(e23_mol, e24_mol)
        assert len(e23_e24_pairs) > 0, "E23-E24のマッチングが見つかりません"

    @pytest.mark.integration
    @pytest.mark.slow
    def test_integrated_substructure_replacement(self, output_dir):
        """integrated_substructure_replacement()の完全なワークフローをテスト"""
        # テスト用の出力ディレクトリ
        test_output = output_dir / "integrated_test"
        
        # 統合ワークフローを実行
        results = integrated_substructure_replacement(
            ligand_file="data/atom_matching/4hw3_A_lig.sdf",
            protein_file="data/sample_proteins/4hw3_A.pdb",
            from_file="data/sample_probes/E23",
            to_file="data/sample_probes/E24",
            output_dir=str(test_output),
            match_index=0  # 最初のマッチを使用
        )
        
        # 結果の検証
        assert len(results) > 0, "結果が返されていません"
        
        for i, result in enumerate(results):
            # 必須キーの存在確認
            assert 'ligand_file' in result, f"パターン{i}: ligand_fileキーが存在しません"
            assert 'protein_file' in result, f"パターン{i}: protein_fileキーが存在しません"
            
            # ファイルの存在確認
            ligand_path = Path(result['ligand_file'])
            protein_path = Path(result['protein_file'])
            
            assert ligand_path.exists(), f"パターン{i}: リガンドファイルが存在しません: {ligand_path}"
            assert protein_path.exists(), f"パターン{i}: タンパク質ファイルが存在しません: {protein_path}"
            
            # ファイルサイズ確認
            assert ligand_path.stat().st_size > 0, f"パターン{i}: リガンドファイルが空です"
            assert protein_path.stat().st_size > 0, f"パターン{i}: タンパク質ファイルが空です"
            
            # リガンドファイルの読み込み確認
            ligand_mol = next(Chem.SDMolSupplier(str(ligand_path)))
            assert ligand_mol is not None, f"パターン{i}: リガンドファイルの読み込みに失敗しました"
            assert ligand_mol.GetNumAtoms() > 0, f"パターン{i}: リガンドに原子がありません"
            
            # タンパク質ファイルの読み込み確認
            protein = PDB.get_structure(str(protein_path))
            coords = PDB.get_attr(protein, "coord")
            assert len(coords) > 0, f"パターン{i}: タンパク質に原子がありません"
    
    @pytest.mark.integration
    def test_integrated_with_match_index(self, output_dir):
        """match_indexパラメータの動作をテスト"""
        test_output = output_dir / "integrated_match_index_test"
        
        # match_index=0で実行
        results = integrated_substructure_replacement(
            ligand_file="data/atom_matching/4hw3_A_lig.sdf",
            protein_file="data/sample_proteins/4hw3_A.pdb",
            from_file="data/sample_probes/E23",
            to_file="data/sample_probes/E24",
            output_dir=str(test_output),
            match_index=0
        )
        
        assert len(results) > 0, "結果が生成されていません"
        
        # ファイルが生成されていることを確認
        for result in results:
            assert Path(result['ligand_file']).exists()
            assert Path(result['protein_file']).exists()
    
    @pytest.mark.integration
    def test_integrated_without_match_index(self, output_dir):
        """match_indexを指定しない場合の動作をテスト（複数マッチ時は可視化）"""
        test_output = output_dir / "integrated_no_index_test"
        
        # match_indexなしで実行（複数マッチがある場合は可視化画像も生成される）
        results = integrated_substructure_replacement(
            ligand_file="data/atom_matching/4hw3_A_lig.sdf",
            protein_file="data/sample_proteins/4hw3_A.pdb",
            from_file="data/sample_probes/E23",
            to_file="data/sample_probes/E24",
            output_dir=str(test_output),
            match_index=None
        )
        
        assert len(results) > 0, "結果が生成されていません"
        
        # 結果ファイルの確認
        for result in results:
            assert Path(result['ligand_file']).exists()
            assert Path(result['protein_file']).exists()
    
    @pytest.mark.integration
    def test_integrated_invalid_match_index(self, output_dir):
        """無効なmatch_indexでエラーが発生することをテスト"""
        test_output = output_dir / "integrated_invalid_index_test"
        
        # 無効なインデックスでエラーが発生することを確認
        with pytest.raises(ValueError, match="無効なmatch_index"):
            integrated_substructure_replacement(
                ligand_file="data/atom_matching/4hw3_A_lig.sdf",
                protein_file="data/sample_proteins/4hw3_A.pdb",
                from_file="data/sample_probes/E23",
                to_file="data/sample_probes/E24",
                output_dir=str(test_output),
                match_index=999  # 明らかに無効なインデックス
            )