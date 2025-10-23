"""Superimpose計算機能の単体テスト"""

import pytest
import numpy as np
from rdkit import Chem
from inverse_msmd.substructure_replacement import (
    find_substructure_in_ligand,
    match_substructures,
    calculate_transformation
)


class TestTransformation:
    """変換行列計算機能のテスト"""

    @pytest.mark.unit
    def test_basic_transformation(self, ligand_mol, e23_mol, e24_mol):
        """基本的な変換行列計算が機能することを確認"""
        # リガンド座標とプローブ座標を取得
        ligand_coords = ligand_mol.GetConformer().GetPositions()
        e23_coords = e23_mol.GetConformer().GetPositions()
        e24_coords = e24_mol.GetConformer().GetPositions()
        
        # マッチング情報を取得
        ligand_e23_matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        assert len(ligand_e23_matches) > 0, "E23のマッチが見つかりません"
        
        atom_pair_patterns = match_substructures(e23_mol, e24_mol)
        assert len(atom_pair_patterns) > 0, "Atom matchingパターンが見つかりません"
        
        # 最初のパターンでテスト
        match = ligand_e23_matches[0]
        atom_pairs = atom_pair_patterns[0]
        
        # リガンドのE23部分の座標を抽出
        ligand_e23_coords = ligand_coords[list(match)]
        
        # 変換行列を計算
        rot, tran = calculate_transformation(
            ligand_e23_coords,
            e24_coords,
            atom_pairs
        )
        
        assert rot is not None, "回転行列がNoneです"
        assert tran is not None, "並進ベクトルがNoneです"

    @pytest.mark.unit
    def test_rotation_matrix_shape(self, ligand_mol, e23_mol, e24_mol):
        """回転行列の形状が正しいことを確認"""
        ligand_coords = ligand_mol.GetConformer().GetPositions()
        e23_coords = e23_mol.GetConformer().GetPositions()
        e24_coords = e24_mol.GetConformer().GetPositions()
        
        ligand_e23_matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        atom_pair_patterns = match_substructures(e23_mol, e24_mol)
        
        match = ligand_e23_matches[0]
        atom_pairs = atom_pair_patterns[0]
        ligand_e23_coords = ligand_coords[list(match)]
        
        rot, tran = calculate_transformation(
            ligand_e23_coords,
            e24_coords,
            atom_pairs
        )
        
        assert rot.shape == (3, 3), \
            f"回転行列は3x3である必要があります（実際: {rot.shape}）"
        assert tran.shape == (3,), \
            f"並進ベクトルは長さ3である必要があります（実際: {tran.shape}）"

    @pytest.mark.unit
    def test_rotation_matrix_orthogonal(self, ligand_mol, e23_mol, e24_mol):
        """回転行列が直交行列であることを確認"""
        ligand_coords = ligand_mol.GetConformer().GetPositions()
        e23_coords = e23_mol.GetConformer().GetPositions()
        e24_coords = e24_mol.GetConformer().GetPositions()
        
        ligand_e23_matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        atom_pair_patterns = match_substructures(e23_mol, e24_mol)
        
        match = ligand_e23_matches[0]
        atom_pairs = atom_pair_patterns[0]
        ligand_e23_coords = ligand_coords[list(match)]
        
        rot, tran = calculate_transformation(
            ligand_e23_coords,
            e24_coords,
            atom_pairs
        )
        
        # 回転行列の直交性をチェック: R^T * R = I
        rot_test = np.dot(rot, rot.T)
        identity = np.eye(3)
        
        assert np.allclose(rot_test, identity, atol=1e-6), \
            "回転行列が直交行列ではありません"

    @pytest.mark.unit
    def test_rotation_matrix_determinant(self, ligand_mol, e23_mol, e24_mol):
        """回転行列の行列式が1に近いことを確認"""
        ligand_coords = ligand_mol.GetConformer().GetPositions()
        e23_coords = e23_mol.GetConformer().GetPositions()
        e24_coords = e24_mol.GetConformer().GetPositions()
        
        ligand_e23_matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        atom_pair_patterns = match_substructures(e23_mol, e24_mol)
        
        match = ligand_e23_matches[0]
        atom_pairs = atom_pair_patterns[0]
        ligand_e23_coords = ligand_coords[list(match)]
        
        rot, tran = calculate_transformation(
            ligand_e23_coords,
            e24_coords,
            atom_pairs
        )
        
        # 行列式が1に近いことをチェック（回転行列の性質）
        det = np.linalg.det(rot)
        
        assert np.abs(det - 1.0) < 1e-6, \
            f"回転行列の行列式が1ではありません: {det}"

    @pytest.mark.unit
    def test_multiple_patterns(self, ligand_mol, e23_mol, e24_mol):
        """複数のマッチングパターンで変換行列が計算できることを確認"""
        ligand_coords = ligand_mol.GetConformer().GetPositions()
        e23_coords = e23_mol.GetConformer().GetPositions()
        e24_coords = e24_mol.GetConformer().GetPositions()
        
        ligand_e23_matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        atom_pair_patterns = match_substructures(e23_mol, e24_mol)
        
        match = ligand_e23_matches[0]
        ligand_e23_coords = ligand_coords[list(match)]
        
        # 全てのパターンでテスト
        for i, atom_pairs in enumerate(atom_pair_patterns[:3]):  # 最初の3パターン
            rot, tran = calculate_transformation(
                ligand_e23_coords,
                e24_coords,
                atom_pairs
            )
            
            assert rot.shape == (3, 3), \
                f"パターン{i}: 回転行列の形状が不正"
            assert tran.shape == (3,), \
                f"パターン{i}: 並進ベクトルの形状が不正"
            
            # 直交性チェック
            rot_test = np.dot(rot, rot.T)
            identity = np.eye(3)
            assert np.allclose(rot_test, identity, atol=1e-6), \
                f"パターン{i}: 回転行列が直交行列ではありません"