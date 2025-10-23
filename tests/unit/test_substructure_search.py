"""部分構造探索機能の単体テスト"""

import pytest
from rdkit import Chem
from inverse_msmd.substructure_replacement import find_substructure_in_ligand


class TestSubstructureSearch:
    """部分構造探索のテスト"""

    @pytest.mark.unit
    def test_find_substructure_basic(self, ligand_mol, e23_mol):
        """基本的な部分構造探索が機能することを確認"""
        matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        
        assert len(matches) > 0, "少なくとも1つのマッチが見つかる必要があります"

    @pytest.mark.unit
    def test_match_atom_count(self, ligand_mol, e23_mol):
        """マッチした原子数が部分構造の原子数と一致することを確認"""
        matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        
        # 水素を除去した部分構造の原子数
        substructure_no_h = Chem.RemoveHs(e23_mol)
        expected_atom_count = substructure_no_h.GetNumAtoms()
        
        for match in matches:
            assert len(match) == expected_atom_count, \
                f"マッチの原子数（{len(match)}）が期待値（{expected_atom_count}）と一致しません"

    @pytest.mark.unit
    def test_match_indices_valid(self, ligand_mol, e23_mol):
        """マッチした原子インデックスが有効範囲内であることを確認"""
        matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        ligand_atom_count = ligand_mol.GetNumAtoms()
        
        for match in matches:
            for atom_idx in match:
                assert 0 <= atom_idx < ligand_atom_count, \
                    f"原子インデックス {atom_idx} が範囲外です（0-{ligand_atom_count-1}）"

    @pytest.mark.unit
    def test_multiple_probes(self, ligand_mol, e23_mol, a01_mol):
        """異なるプローブでの部分構造探索が機能することを確認"""
        # E23プローブ
        e23_matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        assert len(e23_matches) > 0, "E23のマッチが見つかりません"
        
        # A01プローブ
        a01_matches = find_substructure_in_ligand(ligand_mol, a01_mol)
        assert len(a01_matches) > 0, "A01のマッチが見つかりません"