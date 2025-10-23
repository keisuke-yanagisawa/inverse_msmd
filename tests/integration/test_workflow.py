"""統合ワークフローのテスト"""

import pytest
from pathlib import Path
from rdkit import Chem
import numpy as np
from inverse_msmd.substructure_replacement import (
    find_substructure_in_ligand,
    match_substructures,
    calculate_transformation,
    visualize_multiple_matches
)


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