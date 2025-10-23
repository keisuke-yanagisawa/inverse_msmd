"""pytest共通フィクスチャ定義

このファイルはpytestが自動的に読み込み、全てのテストで使用可能なフィクスチャを提供します。
"""

import pytest
from pathlib import Path
from rdkit import Chem
from inverse_msmd.utils.mol_utils import read_mol_from_pdb_smi


# プロジェクトルートパス
@pytest.fixture(scope="session")
def project_root():
    """プロジェクトのルートディレクトリパスを返す"""
    return Path(__file__).parent.parent


# データディレクトリパス
@pytest.fixture(scope="session")
def data_dir(project_root):
    """dataディレクトリのパスを返す"""
    return project_root / "data"


# 出力ディレクトリ
@pytest.fixture(scope="function")
def output_dir(tmp_path):
    """各テスト用の一時出力ディレクトリを返す"""
    output = tmp_path / "test_output"
    output.mkdir(exist_ok=True)
    return output


# テストデータ: リガンド
@pytest.fixture(scope="session")
def ligand_mol(data_dir):
    """4hw3_A_lig.sdfからリガンド分子を読み込む"""
    sdf_path = data_dir / "atom_matching" / "4hw3_A_lig.sdf"
    if not sdf_path.exists():
        pytest.skip(f"Test data not found: {sdf_path}")
    return next(Chem.SDMolSupplier(str(sdf_path)))


# テストデータ: E23プローブ
@pytest.fixture(scope="session")
def e23_mol(data_dir):
    """E23プローブ分子を読み込む"""
    pdb_path = data_dir / "sample_probes" / "E23.pdb"
    smi_path = data_dir / "sample_probes" / "E23.smi"
    if not pdb_path.exists() or not smi_path.exists():
        pytest.skip(f"Test data not found: {pdb_path} or {smi_path}")
    return read_mol_from_pdb_smi(str(pdb_path), str(smi_path))


# テストデータ: E24プローブ
@pytest.fixture(scope="session")
def e24_mol(data_dir):
    """E24プローブ分子を読み込む"""
    pdb_path = data_dir / "sample_probes" / "E24.pdb"
    smi_path = data_dir / "sample_probes" / "E24.smi"
    if not pdb_path.exists() or not smi_path.exists():
        pytest.skip(f"Test data not found: {pdb_path} or {smi_path}")
    return read_mol_from_pdb_smi(str(pdb_path), str(smi_path))


# テストデータ: A01プローブ
@pytest.fixture(scope="session")
def a01_mol(data_dir):
    """A01プローブ分子を読み込む"""
    pdb_path = data_dir / "sample_probes" / "A01.pdb"
    smi_path = data_dir / "sample_probes" / "A01.smi"
    if not pdb_path.exists() or not smi_path.exists():
        pytest.skip(f"Test data not found: {pdb_path} or {smi_path}")
    return read_mol_from_pdb_smi(str(pdb_path), str(smi_path))


# テストデータ: A08プローブ
@pytest.fixture(scope="session")
def a08_mol(data_dir):
    """A08プローブ分子を読み込む"""
    pdb_path = data_dir / "sample_probes" / "A08.pdb"
    smi_path = data_dir / "sample_probes" / "A08.smi"
    if not pdb_path.exists() or not smi_path.exists():
        pytest.skip(f"Test data not found: {pdb_path} or {smi_path}")
    return read_mol_from_pdb_smi(str(pdb_path), str(smi_path))


# テストデータ: タンパク質構造
@pytest.fixture(scope="session")
def protein_pdb_path(data_dir):
    """タンパク質PDBファイルのパスを返す"""
    pdb_path = data_dir / "sample_proteins" / "4hw3_A.pdb"
    if not pdb_path.exists():
        pytest.skip(f"Test data not found: {pdb_path}")
    return str(pdb_path)