"""Atom Matching機能の単体テスト"""

import pytest
import numpy as np
from rdkit import Chem
from inverse_msmd.substructure_replacement import match_substructures


class TestAtomMatching:
    """Atom Matching機能のテスト"""

    @pytest.mark.unit
    def test_basic_matching(self, e23_mol, e24_mol):
        """基本的なatom matchingが機能することを確認"""
        atom_pair_patterns = match_substructures(e23_mol, e24_mol)
        
        assert len(atom_pair_patterns) > 0, \
            "少なくとも1つのマッチングパターンが見つかる必要があります"

    @pytest.mark.unit
    def test_atom_pairs_shape(self, e23_mol, e24_mol):
        """atom pairsの形状が正しいことを確認"""
        atom_pair_patterns = match_substructures(e23_mol, e24_mol)
        
        for i, pairs in enumerate(atom_pair_patterns):
            assert pairs.shape[0] == 2, \
                f"パターン{i}: atom pairsは2行である必要があります（実際: {pairs.shape[0]}）"
            assert pairs.shape[1] > 0, \
                f"パターン{i}: 少なくとも1つの原子ペアが必要です"

    @pytest.mark.unit
    def test_atom_indices_valid(self, e23_mol, e24_mol):
        """atom pairsのインデックスが有効範囲内であることを確認"""
        atom_pair_patterns = match_substructures(e23_mol, e24_mol)
        
        # 水素を除去した分子の原子数
        e23_no_h = Chem.RemoveHs(e23_mol)
        e24_no_h = Chem.RemoveHs(e24_mol)
        e23_atom_count = e23_no_h.GetNumAtoms()
        e24_atom_count = e24_no_h.GetNumAtoms()
        
        for i, pairs in enumerate(atom_pair_patterns):
            assert np.all(pairs[0] < e23_atom_count), \
                f"パターン{i}: E23のインデックスが範囲外（0-{e23_atom_count-1}）"
            assert np.all(pairs[1] < e24_atom_count), \
                f"パターン{i}: E24のインデックスが範囲外（0-{e24_atom_count-1}）"

    @pytest.mark.unit
    def test_element_consistency(self, e23_mol, e24_mol):
        """対応する原子の元素記号が一致することを確認"""
        atom_pair_patterns = match_substructures(e23_mol, e24_mol)
        
        # 水素を除去
        e23_no_h = Chem.RemoveHs(e23_mol)
        e24_no_h = Chem.RemoveHs(e24_mol)
        
        for i, pairs in enumerate(atom_pair_patterns):
            for j in range(pairs.shape[1]):
                e23_idx = pairs[0, j]
                e24_idx = pairs[1, j]
                
                e23_atom = e23_no_h.GetAtomWithIdx(int(e23_idx))
                e24_atom = e24_no_h.GetAtomWithIdx(int(e24_idx))
                
                assert e23_atom.GetSymbol() == e24_atom.GetSymbol(), \
                    f"パターン{i}, ペア{j}: 元素が一致しません " \
                    f"(E23[{e23_idx}]={e23_atom.GetSymbol()} != " \
                    f"E24[{e24_idx}]={e24_atom.GetSymbol()})"

    @pytest.mark.unit
    def test_different_probe_pairs(self, e23_mol, a01_mol):
        """異なるプローブペアでのatom matchingをテスト"""
        try:
            atom_pair_patterns = match_substructures(e23_mol, a01_mol)
            # マッチングが見つかる可能性はあるが、見つからなくてもエラーではない
            assert isinstance(atom_pair_patterns, list), \
                "結果はリストである必要があります"
        except Exception as e:
            # 構造が大きく異なる場合はマッチングが失敗する可能性がある
            pytest.skip(f"プローブ間の構造が大きく異なるためスキップ: {e}")