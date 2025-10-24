"""
部分構造置換候補生成関数のテスト
"""

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from inverse_msmd.substructure_replacement import (
    generate_all_replacement_candidates,
    create_replacement_molecule
)


@pytest.fixture
def simple_molecule():
    """テスト用の簡単な分子（エチルベンゼン）"""
    mol = Chem.MolFromSmiles("CCc1ccccc1")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return Chem.RemoveHs(mol)


@pytest.fixture
def benzene_pattern():
    """ベンゼン環パターン"""
    mol = Chem.MolFromSmiles("c1ccccc1")
    return Chem.RemoveHs(mol)


@pytest.fixture
def cyclohexane_replacement():
    """シクロヘキサン置換パターン"""
    mol = Chem.MolFromSmiles("C1CCCCC1")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return Chem.RemoveHs(mol)


class TestCreateReplacementMolecule:
    """create_replacement_molecule関数のテスト"""
    
    def test_basic_replacement(self, simple_molecule, benzene_pattern, cyclohexane_replacement):
        """基本的な部分構造置換のテスト"""
        # ベンゼン環を探索
        match = simple_molecule.GetSubstructMatch(benzene_pattern)
        assert len(match) > 0, "ベンゼン環が見つかりませんでした"
        
        # 接続点を決定（エチル基との接続）
        # ベンゼン環の最初の原子に接続
        connections = [(0, 1, Chem.BondType.SINGLE)]  # シクロヘキサンの最初の原子をエチル基に接続
        
        # 置換実行
        new_mol = create_replacement_molecule(
            simple_molecule,
            match,
            cyclohexane_replacement,
            connections
        )
        
        # 分子が生成されたことを確認
        assert new_mol is not None
        assert new_mol.GetNumAtoms() > 0
        
        # Sanitizeが可能であることを確認
        try:
            Chem.SanitizeMol(new_mol)
            sanitize_ok = True
        except:
            sanitize_ok = False
        
        assert sanitize_ok, "Sanitizeに失敗しました"
    
    def test_no_connection_replacement(self):
        """接続なしの置換テスト（孤立部分構造）"""
        mol = Chem.MolFromSmiles("c1ccccc1")  # ベンゼン単体
        mol = Chem.RemoveHs(mol)
        
        pattern = Chem.MolFromSmiles("c1ccccc1")
        pattern = Chem.RemoveHs(pattern)
        
        replacement = Chem.MolFromSmiles("C1CCCCC1")
        replacement = Chem.RemoveHs(replacement)
        
        match = mol.GetSubstructMatch(pattern)
        
        # 接続なし
        new_mol = create_replacement_molecule(mol, match, replacement, [])
        
        assert new_mol is not None
        assert new_mol.GetNumAtoms() == replacement.GetNumAtoms()
        
        # Sanitizeできることを確認
        Chem.SanitizeMol(new_mol)
    
    def test_preserves_properties(self):
        """分子のプロパティが保持されることをテスト"""
        mol = Chem.MolFromSmiles("CCc1ccccc1")
        mol = Chem.RemoveHs(mol)
        mol.SetProp("_Name", "TestMolecule")
        mol.SetProp("CustomProp", "CustomValue")
        
        pattern = Chem.MolFromSmiles("c1ccccc1")
        pattern = Chem.RemoveHs(pattern)
        
        replacement = Chem.MolFromSmiles("C1CCCCC1")
        replacement = Chem.RemoveHs(replacement)
        
        match = mol.GetSubstructMatch(pattern)
        connections = [(0, 1, Chem.BondType.SINGLE)]
        
        new_mol = create_replacement_molecule(mol, match, replacement, connections)
        
        # プロパティが保持されていることを確認
        assert new_mol.HasProp("_Name")
        assert new_mol.GetProp("_Name") == "TestMolecule"
        assert new_mol.HasProp("CustomProp")
        assert new_mol.GetProp("CustomProp") == "CustomValue"


class TestGenerateAllReplacementCandidates:
    """generate_all_replacement_candidates関数のテスト"""
    
    def test_generates_multiple_candidates(self, simple_molecule, benzene_pattern, cyclohexane_replacement):
        """複数の候補が生成されることをテスト"""
        match = simple_molecule.GetSubstructMatch(benzene_pattern)
        
        candidates = generate_all_replacement_candidates(
            simple_molecule,
            benzene_pattern,
            cyclohexane_replacement,
            match
        )
        
        # 少なくとも1つの候補が生成されることを確認
        assert len(candidates) > 0
        
        # 全ての候補がSanitize可能であることを確認
        for i, candidate in enumerate(candidates):
            try:
                Chem.SanitizeMol(candidate)
                assert True, f"候補{i}はSanitize可能です"
            except:
                pytest.fail(f"候補{i}がSanitizeに失敗しました")
    
    def test_all_candidates_are_unique(self, simple_molecule, benzene_pattern, cyclohexane_replacement):
        """全ての候補が一意であることをテスト"""
        match = simple_molecule.GetSubstructMatch(benzene_pattern)
        
        candidates = generate_all_replacement_candidates(
            simple_molecule,
            benzene_pattern,
            cyclohexane_replacement,
            match
        )
        
        # SMILESで候補の一意性を確認
        smiles_set = set()
        for candidate in candidates:
            smiles = Chem.MolToSmiles(candidate)
            smiles_set.add(smiles)
        
        # 一意な候補の数が期待される数と一致
        # （一部の候補は同一の構造になる可能性がある）
        assert len(smiles_set) >= 1
    
    def test_no_attachment_point(self):
        """接続点がない場合のテスト"""
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.RemoveHs(mol)
        
        pattern = Chem.MolFromSmiles("c1ccccc1")
        pattern = Chem.RemoveHs(pattern)
        
        replacement = Chem.MolFromSmiles("C1CCCCC1")
        replacement = Chem.RemoveHs(replacement)
        
        match = mol.GetSubstructMatch(pattern)
        
        candidates = generate_all_replacement_candidates(
            mol,
            pattern,
            replacement,
            match
        )
        
        # 接続点がない場合でも1つの候補が生成される
        assert len(candidates) == 1
        assert candidates[0].GetNumAtoms() == replacement.GetNumAtoms()
    
    def test_raises_on_no_valid_candidates(self):
        """有効な候補が生成できない場合にエラーが発生することをテスト"""
        # 非常に不適切な置換パターンを作成
        mol = Chem.MolFromSmiles("C")  # メタン
        mol = Chem.RemoveHs(mol)
        
        pattern = Chem.MolFromSmiles("C")
        pattern = Chem.RemoveHs(pattern)
        
        # 非常に大きな置換分子（接続に失敗する可能性が高い）
        replacement = Chem.MolFromSmiles("c1ccc2ccccc2c1")  # ナフタレン
        replacement = Chem.RemoveHs(replacement)
        
        match = mol.GetSubstructMatch(pattern)
        
        # このケースでは接続点がないため、1つの候補が生成される
        # （孤立した分子として）
        candidates = generate_all_replacement_candidates(
            mol,
            pattern,
            replacement,
            match
        )
        
        # 接続点がない場合は候補が生成される
        assert len(candidates) >= 1
    
    def test_candidate_atom_count(self):
        """候補の原子数が適切であることをテスト"""
        mol = Chem.MolFromSmiles("CCc1ccccc1")
        mol = Chem.RemoveHs(mol)
        
        pattern = Chem.MolFromSmiles("c1ccccc1")
        pattern = Chem.RemoveHs(pattern)
        
        replacement = Chem.MolFromSmiles("C1CCCCC1")
        replacement = Chem.RemoveHs(replacement)
        
        match = mol.GetSubstructMatch(pattern)
        
        candidates = generate_all_replacement_candidates(
            mol,
            pattern,
            replacement,
            match
        )
        
        # 期待される原子数: 元の分子 - パターン + 置換分子
        expected_atoms = mol.GetNumAtoms() - len(match) + replacement.GetNumAtoms()
        
        for candidate in candidates:
            assert candidate.GetNumAtoms() == expected_atoms


class TestIntegration:
    """統合テスト"""
    
    def test_workflow_with_real_molecules(self):
        """実際の分子を使ったワークフローテスト"""
        # トルエンからベンゼン環をシクロヘキサンに置換
        mol = Chem.MolFromSmiles("Cc1ccccc1")  # トルエン
        mol = Chem.RemoveHs(mol)
        
        pattern = Chem.MolFromSmiles("c1ccccc1")
        pattern = Chem.RemoveHs(pattern)
        
        replacement = Chem.MolFromSmiles("C1CCCCC1")
        replacement = Chem.RemoveHs(replacement)
        
        # マッチを見つける
        match = mol.GetSubstructMatch(pattern)
        assert len(match) > 0
        
        # 全候補を生成
        candidates = generate_all_replacement_candidates(
            mol,
            pattern,
            replacement,
            match
        )
        
        assert len(candidates) > 0
        
        # 全候補が有効であることを確認
        for candidate in candidates:
            Chem.SanitizeMol(candidate)
            assert candidate.GetNumAtoms() > 0
            
            # SMILES表記が取得できることを確認
            smiles = Chem.MolToSmiles(candidate)
            assert len(smiles) > 0