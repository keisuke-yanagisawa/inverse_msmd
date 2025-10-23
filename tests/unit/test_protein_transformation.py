"""タンパク質変換機能の単体テスト"""

import pytest
import numpy as np
from rdkit import Chem
from inverse_msmd.substructure_replacement import (
    find_substructure_in_ligand,
    match_substructures,
    calculate_transformation,
    apply_transformation_to_protein
)
from inverse_msmd.utils.bio_utils import PDB


class TestProteinTransformation:
    """タンパク質変換機能のテスト"""

    @pytest.mark.unit
    def test_basic_protein_transformation(self, protein_pdb_path, ligand_mol, e23_mol, e24_mol):
        """基本的なタンパク質変換が機能することを確認"""
        # タンパク質を読み込み
        protein = PDB.get_structure(protein_pdb_path)
        
        # 変換行列を計算
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
        
        # 元の座標を保存
        original_coords = PDB.get_attr(protein, "coord").copy()
        
        # 変換を適用
        transformed_protein = apply_transformation_to_protein(protein, rot, tran)
        
        assert transformed_protein is not None, "変換後のタンパク質がNoneです"

    @pytest.mark.unit
    def test_coordinates_shape_preserved(self, protein_pdb_path, ligand_mol, e23_mol, e24_mol):
        """座標の形状が保持されることを確認"""
        protein = PDB.get_structure(protein_pdb_path)
        original_coords = PDB.get_attr(protein, "coord").copy()
        
        # 変換行列を計算
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
        
        # 変換を適用
        transformed_protein = apply_transformation_to_protein(protein, rot, tran)
        new_coords = PDB.get_attr(transformed_protein, "coord")
        
        assert new_coords.shape == original_coords.shape, \
            f"座標の形状が変わっています（元: {original_coords.shape}, 変換後: {new_coords.shape}）"

    @pytest.mark.unit
    def test_coordinates_actually_transformed(self, protein_pdb_path, ligand_mol, e23_mol, e24_mol):
        """座標が実際に変換されていることを確認"""
        protein = PDB.get_structure(protein_pdb_path)
        original_coords = PDB.get_attr(protein, "coord").copy()
        
        # 変換行列を計算
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
        
        # 変換を適用
        transformed_protein = apply_transformation_to_protein(protein, rot, tran)
        new_coords = PDB.get_attr(transformed_protein, "coord")
        
        # 座標が変化していることを確認
        assert not np.allclose(new_coords, original_coords), \
            "座標が変換されていません"

    @pytest.mark.unit
    def test_transformation_formula(self, protein_pdb_path, ligand_mol, e23_mol, e24_mol):
        """変換式（rot @ coords + tran）が正しく適用されていることを確認"""
        protein = PDB.get_structure(protein_pdb_path)
        original_coords = PDB.get_attr(protein, "coord").copy()
        
        # 変換行列を計算
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
        
        # 変換を適用
        transformed_protein = apply_transformation_to_protein(protein, rot, tran)
        new_coords = PDB.get_attr(transformed_protein, "coord")
        
        # 期待される座標を手動計算
        expected_coords = np.dot(original_coords, rot) + tran
        
        assert np.allclose(new_coords, expected_coords, atol=1e-6), \
            "変換式（rot @ coords + tran）が正しく適用されていません"

    @pytest.mark.unit
    def test_pdb_save_and_load(self, protein_pdb_path, ligand_mol, e23_mol, e24_mol, output_dir):
        """変換後のタンパク質をPDBファイルとして保存・読み込みできることを確認"""
        protein = PDB.get_structure(protein_pdb_path)
        
        # 変換行列を計算
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
        
        # 変換を適用
        transformed_protein = apply_transformation_to_protein(protein, rot, tran)
        
        # PDBファイルとして保存
        output_path = output_dir / "protein_transformed_test.pdb"
        PDB.save(transformed_protein, str(output_path))
        
        # ファイルが生成されたことを確認
        assert output_path.exists(), "PDBファイルが生成されていません"
        assert output_path.stat().st_size > 0, "PDBファイルのサイズが0です"
        
        # 保存したファイルを読み込めることを確認
        loaded_protein = PDB.get_structure(str(output_path))
        loaded_coords = PDB.get_attr(loaded_protein, "coord")
        transformed_coords = PDB.get_attr(transformed_protein, "coord")
        
        assert np.allclose(loaded_coords, transformed_coords, atol=1e-3), \
            "保存・読み込み後の座標が一致しません"