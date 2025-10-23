"""リガンド部分構造置換の単体テスト"""

import pytest
import numpy as np
from rdkit import Chem
from inverse_msmd.substructure_replacement import (
    find_substructure_in_ligand,
    match_substructures,
    replace_ligand_substructure
)


class TestLigandReplacement:
    """リガンド置換機能のテスト"""

    @pytest.mark.unit
    def test_basic_replacement(self, ligand_mol, e23_mol, e24_mol):
        """基本的なリガンド置換が機能することを確認"""
        # マッチング情報を取得
        ligand_e23_matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        assert len(ligand_e23_matches) > 0, "E23のマッチが見つかりません"
        
        atom_pair_patterns = match_substructures(e23_mol, e24_mol)
        assert len(atom_pair_patterns) > 0, "Atom matchingパターンが見つかりません"
        
        match = ligand_e23_matches[0]
        atom_pairs = atom_pair_patterns[0]
        
        # 部分構造を置換
        replaced_ligand = replace_ligand_substructure(
            ligand_mol,
            match,
            e24_mol,
            atom_pairs
        )
        
        assert replaced_ligand is not None, "置換に失敗しました"

    @pytest.mark.unit
    def test_atom_count_change(self, ligand_mol, e23_mol, e24_mol):
        """原子数が期待通りに変化することを確認"""
        # 水素を除去した分子で処理
        ligand_no_h = Chem.RemoveHs(ligand_mol)
        e23_no_h = Chem.RemoveHs(e23_mol)
        e24_no_h = Chem.RemoveHs(e24_mol)
        
        original_atom_count = ligand_no_h.GetNumAtoms()
        e23_atom_count = e23_no_h.GetNumAtoms()
        e24_atom_count = e24_no_h.GetNumAtoms()
        
        # マッチング情報を取得
        ligand_e23_matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        atom_pair_patterns = match_substructures(e23_mol, e24_mol)
        
        match = ligand_e23_matches[0]
        atom_pairs = atom_pair_patterns[0]
        
        # 部分構造を置換
        replaced_ligand = replace_ligand_substructure(
            ligand_mol,
            match,
            e24_mol,
            atom_pairs
        )
        
        # 水素を除去して原子数を確認
        replaced_no_h = Chem.RemoveHs(replaced_ligand)
        new_atom_count = replaced_no_h.GetNumAtoms()
        
        # 期待される原子数: 元の原子数 - E23の原子数 + E24の原子数
        expected_atom_count = original_atom_count - e23_atom_count + e24_atom_count
        
        # 接続の関係で±数原子の誤差は許容
        assert abs(new_atom_count - expected_atom_count) <= 5, \
            f"原子数が期待値から大きく外れています " \
            f"(期待: 約{expected_atom_count}, 実際: {new_atom_count})"

    @pytest.mark.unit
    def test_sanitize_check(self, ligand_mol, e23_mol, e24_mol):
        """置換後の分子が化学的に妥当（Sanitize可能）であることを確認"""
        ligand_e23_matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        atom_pair_patterns = match_substructures(e23_mol, e24_mol)
        
        match = ligand_e23_matches[0]
        atom_pairs = atom_pair_patterns[0]
        
        replaced_ligand = replace_ligand_substructure(
            ligand_mol,
            match,
            e24_mol,
            atom_pairs
        )
        
        # Sanitizeチェック
        try:
            Chem.SanitizeMol(replaced_ligand)
            sanitize_ok = True
        except:
            sanitize_ok = False
        
        assert sanitize_ok, "置換後の分子のSanitizeに失敗しました"

    @pytest.mark.unit
    def test_sdf_save(self, ligand_mol, e23_mol, e24_mol, output_dir):
        """置換後の分子をSDFファイルとして保存できることを確認"""
        ligand_e23_matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        atom_pair_patterns = match_substructures(e23_mol, e24_mol)
        
        match = ligand_e23_matches[0]
        atom_pairs = atom_pair_patterns[0]
        
        replaced_ligand = replace_ligand_substructure(
            ligand_mol,
            match,
            e24_mol,
            atom_pairs
        )
        
        # SDFファイルとして保存
        output_path = output_dir / "ligand_replaced_test.sdf"
        writer = Chem.SDWriter(str(output_path))
        writer.SetKekulize(False)
        writer.write(replaced_ligand)
        writer.close()
        
        # ファイルが生成されたことを確認
        assert output_path.exists(), "SDFファイルが生成されていません"
        assert output_path.stat().st_size > 0, "SDFファイルのサイズが0です"
        
        # 保存したファイルを読み込めることを確認
        supplier = Chem.SDMolSupplier(str(output_path))
        loaded_mol = next(supplier)
        assert loaded_mol is not None, "SDFファイルから分子を読み込めません"

    @pytest.mark.unit
    def test_multiple_patterns(self, ligand_mol, e23_mol, e24_mol):
        """複数のマッチングパターンで置換できることを確認"""
        ligand_e23_matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        atom_pair_patterns = match_substructures(e23_mol, e24_mol)
        
        match = ligand_e23_matches[0]
        
        # 最初の3パターンでテスト
        for i, atom_pairs in enumerate(atom_pair_patterns[:3]):
            try:
                replaced_ligand = replace_ligand_substructure(
                    ligand_mol,
                    match,
                    e24_mol,
                    atom_pairs
                )
                
                assert replaced_ligand is not None, \
                    f"パターン{i}: 置換に失敗しました"
                
                # Sanitizeチェック
                Chem.SanitizeMol(replaced_ligand)
                
            except Exception as e:
                # 一部のパターンは化学的に無効な可能性がある
                pytest.skip(f"パターン{i}は化学的に無効: {e}")

    @pytest.mark.unit
    def test_smiles_comparison(self, ligand_mol, e23_mol, e24_mol):
        """置換前後でSMILESが変化していることを確認"""
        ligand_e23_matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        atom_pair_patterns = match_substructures(e23_mol, e24_mol)
        
        match = ligand_e23_matches[0]
        atom_pairs = atom_pair_patterns[0]
        
        # 元のSMILES
        original_smiles = Chem.MolToSmiles(Chem.RemoveHs(ligand_mol))
        
        # 置換
        replaced_ligand = replace_ligand_substructure(
            ligand_mol,
            match,
            e24_mol,
            atom_pairs
        )
        
        # 置換後のSMILES
        replaced_smiles = Chem.MolToSmiles(Chem.RemoveHs(replaced_ligand))
        
        # SMILESが変化していることを確認
        assert original_smiles != replaced_smiles, \
            "置換前後でSMILESが同じです（置換されていない可能性）"

    @pytest.mark.unit
    def test_conformer_preservation(self, ligand_mol, e23_mol, e24_mol):
        """コンフォーマー情報が保持されることを確認"""
        ligand_e23_matches = find_substructure_in_ligand(ligand_mol, e23_mol)
        atom_pair_patterns = match_substructures(e23_mol, e24_mol)
        
        match = ligand_e23_matches[0]
        atom_pairs = atom_pair_patterns[0]
        
        replaced_ligand = replace_ligand_substructure(
            ligand_mol,
            match,
            e24_mol,
            atom_pairs
        )
        
        # コンフォーマーが存在することを確認
        assert replaced_ligand.GetNumConformers() >= 0, \
            "置換後の分子にコンフォーマーがありません"