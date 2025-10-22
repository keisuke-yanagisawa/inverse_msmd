# coding: utf-8

"""
分子データ操作ユーティリティモジュール

分子構造データの読み込み、変換、出力機能を提供します。
PDBファイルとSMIファイルから、座標と結合情報の両方を含む
分子構造データを作成できます。

このモジュールはRDKitを使用します。
version: 1.0.0
Authors: Keisuke Yanagisawa
"""

from typing import Optional, Tuple
import numpy as np
import numpy.typing as npt
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem


def _create_mol_from_pdb_with_smiles_template(pdb_file: str, smiles: str) -> 'Chem.Mol':
    """
    PDBファイルとSMILESテンプレートから分子を作成（内部関数）
    
    この関数は、PDBファイルから座標を読み込み、SMILESから結合情報を取得して、
    両方の情報を統合した分子オブジェクトを作成します。原子の順序が異なる場合でも
    自動的にマッチングを行います。
    
    注意: これは内部実装の詳細です。通常はcreate_sdf_from_pdb_smi()を使用してください。
    
    Parameters
    ----------
    pdb_file : str
        PDBファイルのパス
    smiles : str
        SMILES文字列（結合情報のテンプレート）
        
    Returns
    -------
    Chem.Mol
        座標と結合情報を持つ分子オブジェクト
        
    Raises
    ------
    ValueError
        分子の作成に失敗した場合
    """
    
    # PDBファイルから分子を読み込む（結合情報なし）
    pdb_mol = Chem.MolFromPDBFile(pdb_file, removeHs=False, sanitize=False)
    if pdb_mol is None:
        raise ValueError(f"PDBファイルから分子を作成できませんでした: {pdb_file}")
    
    # SMILESから分子を作成（結合情報あり、座標なし）
    template_mol = Chem.MolFromSmiles(smiles)
    if template_mol is None:
        raise ValueError(f"SMILESから分子を作成できませんでした: {smiles}")
    
    # 水素を追加
    template_mol = Chem.AddHs(template_mol)
    
    # 原子数のチェック
    if pdb_mol.GetNumAtoms() != template_mol.GetNumAtoms():
        raise ValueError(
            f"原子数が一致しません (PDB: {pdb_mol.GetNumAtoms()}, SMILES: {template_mol.GetNumAtoms()})"
        )
    
    try:
        # テンプレートから結合情報を割り当て
        # この関数は自動的に原子のマッチングを行います
        mol_with_bonds = AllChem.AssignBondOrdersFromTemplate(template_mol, pdb_mol)
        return mol_with_bonds
    except Exception as e:
        raise ValueError(f"テンプレートマッチングに失敗しました: {e}")


def read_mol_from_pdb_smi(pdb_file: str, smi_file: str, verbose: bool = False) -> 'Chem.Mol':
    """
    PDBファイルとSMIファイルから分子オブジェクトを読み込む
    
    PDBファイルから3D座標情報を、SMIファイルから結合情報を取得し、
    両方の情報を含むRDKit分子オブジェクトを返します。
    テンプレートマッチングにより、原子の順序が異なる場合でも自動的にマッチングします。
    
    Parameters
    ----------
    pdb_file : str
        PDBファイルのパス（座標情報）
    smi_file : str
        SMIファイルのパス（結合情報、SMILES文字列）
    verbose : bool, optional
        進捗メッセージを表示するか（デフォルト: False）
        
    Returns
    -------
    Chem.Mol
        座標と結合情報を持つRDKit分子オブジェクト
        
    Raises
    ------
    FileNotFoundError
        入力ファイルが存在しない場合
    ValueError
        原子数または元素が一致しない場合、分子の作成に失敗した場合
        
    Notes
    -----
    - PDBファイルとSMILESから生成される分子の原子数が一致している必要があります
    - テンプレートマッチングにより、原子の順序が異なっていても自動的にマッチングします
    - SMILESに水素が明示的に含まれていない場合、自動的に追加されます
    - 返された分子オブジェクトは、ユーザーが自由に加工・保存できます
    
    Examples
    --------
    >>> # 分子を読み込む
    >>> mol = read_mol_from_pdb_smi(
    ...     "data/sample_probes/A08.pdb",
    ...     "data/sample_probes/A08.smi"
    ... )
    
    >>> # 分子名を設定
    >>> mol.SetProp("_Name", "Toluene")
    
    >>> # SDFファイルとして保存
    >>> from rdkit import Chem
    >>> writer = Chem.SDWriter("output/A08.sdf")
    >>> writer.write(mol)
    >>> writer.close()
    
    >>> # または、さらに加工してから保存
    >>> mol.SetProp("MyProperty", "MyValue")
    >>> # 他の処理...
    """
    # ファイルの存在確認
    if not Path(pdb_file).exists():
        raise FileNotFoundError(f"PDBファイルが見つかりません: {pdb_file}")
    if not Path(smi_file).exists():
        raise FileNotFoundError(f"SMIファイルが見つかりません: {smi_file}")
    
    # SMIファイルからSMILESを読み込み
    if verbose:
        print(f"SMIファイルを読み込み中: {smi_file}")
    
    with open(smi_file, 'r') as f:
        # 最初の行を読み取る（SMILES文字列のみ、または SMILES NAME 形式）
        line = f.readline().strip()
        # スペースで分割して最初の要素を取得
        smiles = line.split()[0]
    
    if verbose:
        print(f"  SMILES: {smiles}")
    
    # テンプレートマッチングを使用して分子を作成
    if verbose:
        print(f"PDBファイルを読み込み中: {pdb_file}")
        print("  テンプレートマッチングを使用（原子順序が異なってもOK）")
    
    mol = _create_mol_from_pdb_with_smiles_template(pdb_file, smiles)
        
    if verbose:
        print(f"  原子数: {mol.GetNumAtoms()}")
        print("✓ 分子の読み込み成功")
    
    return mol


__all__ = [
    "read_mol_from_pdb_smi",
]