"""
Structure alignment module.

このモジュールは、プローブ分子と共結晶化リガンドのatom matchingに基づいて
タンパク質構造の重ね合わせを行う機能を提供します。

主要な機能:
- MCS（最大共通部分構造）検索
- 原子ペアの特定
- 構造の重ね合わせ
- バッチ処理のサポート
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Union
import numpy as np
import numpy.typing as npt
from rdkit import Chem
from rdkit.Chem import rdFMCS
from Bio.PDB.Structure import Structure

from .utils.bio_utils import SuperImposer, PDB
from .utils.path_utils import expandpath


@dataclass
class AlignmentResult:
    """
    重ね合わせ結果を保持するデータクラス。
    
    Attributes
    ----------
    probe_id : str
        プローブ分子のID
    match_id : int
        マッチングID（同一プローブで複数マッチがある場合の識別子）
    aligned_protein : Structure
        重ね合わせ後のタンパク質構造
    atom_pairs : np.ndarray
        原子ペアのインデックス配列 (shape: (2, n_atoms))
        atom_pairs[0] がプローブ側、atom_pairs[1] がリガンド側のインデックス
    """
    probe_id: str
    match_id: int
    aligned_protein: Structure
    atom_pairs: npt.NDArray[np.int_]


def _get_iso_labeled_atoms(mol: Chem.Mol, iso_value: int = 1) -> List[int]:
    """
    指定されたアイソトープラベル値を持つ原子のインデックスを取得します。
    
    Parameters
    ----------
    mol : Chem.Mol
        RDKit分子オブジェクト
    iso_value : int, optional
        抽出するアイソトープラベル値（デフォルト: 1）
    
    Returns
    -------
    List[int]
        アイソトープラベルを持つ原子のインデックスリスト
    """
    iso_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetIsotope() == iso_value:
            iso_atoms.append(atom.GetIdx())
    return iso_atoms


def _create_submol_from_atoms(mol: Chem.Mol, atom_indices: List[int]) -> Chem.Mol:
    """
    指定された原子インデックスから部分構造を作成します。
    
    Parameters
    ----------
    mol : Chem.Mol
        元のRDKit分子オブジェクト
    atom_indices : List[int]
        抽出する原子のインデックスリスト
    
    Returns
    -------
    Chem.Mol
        部分構造分子
    """
    atom_set = set(atom_indices)
    
    # これらの原子間の結合を抽出
    bonds_to_include = []
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in atom_set and bond.GetEndAtomIdx() in atom_set:
            bonds_to_include.append(bond.GetIdx())
    
    # 部分構造を作成
    submol = Chem.PathToSubmol(mol, bonds_to_include)
    return submol


def find_atom_matches(
    probe_mol: Chem.Mol,
    ref_mol: Chem.Mol,
    iso_value: int = 1
) -> List[npt.NDArray[np.int_]]:
    """
    プローブ分子とリファレンス分子間のMCS（最大共通部分構造）を検索し、
    すべての原子ペアマッチングを返します。
    
    リファレンス分子にアイソトープラベルがある場合、
    指定されたアイソトープ値を持つ原子のみを考慮します。
    
    Parameters
    ----------
    probe_mol : Chem.Mol
        プローブ分子（RDKit Mol オブジェクト）
    ref_mol : Chem.Mol
        リファレンス分子（共結晶化リガンド等）
    iso_value : int, optional
        リファレンス分子のアイソトープラベル値（デフォルト: 1）
        0の場合、すべての原子を考慮します
    
    Returns
    -------
    List[np.ndarray]
        原子ペアマッチングのリスト。各要素は shape (2, n_atoms) の配列で、
        [0, :] がプローブ側、[1, :] がリファレンス側の原子インデックス
    """
    # アイソトープラベルがある場合、部分構造を作成
    if iso_value > 0:
        iso_atoms = _get_iso_labeled_atoms(ref_mol, iso_value)
        if not iso_atoms:
            # アイソトープラベルがない場合は全体を使用
            ref_submol = ref_mol
            submol_to_mol_map = {i: i for i in range(ref_mol.GetNumAtoms())}
        else:
            ref_submol = _create_submol_from_atoms(ref_mol, iso_atoms)
            # 部分構造のインデックスを元の分子のインデックスにマッピング
            submol_to_mol_map = {}
            for submol_idx, atom in enumerate(ref_submol.GetAtoms()):
                original_idx = atom.GetAtomMapNum()
                if original_idx == 0:
                    if submol_idx < len(iso_atoms):
                        original_idx = iso_atoms[submol_idx]
                else:
                    original_idx -= 1  # 1-indexed to 0-indexed
                submol_to_mol_map[submol_idx] = original_idx
    else:
        ref_submol = ref_mol
        submol_to_mol_map = {i: i for i in range(ref_mol.GetNumAtoms())}
    
    # MCS検索
    mcs_result = rdFMCS.FindMCS([probe_mol, ref_submol])
    if mcs_result.numAtoms == 0:
        return []
    
    mcs = Chem.MolFromSmarts(mcs_result.smartsString)
    
    # すべてのマッチングを取得
    probe_matches = probe_mol.GetSubstructMatches(mcs, uniquify=False)
    ref_submol_matches = ref_submol.GetSubstructMatches(mcs, uniquify=False)
    
    # 重複除外のため、原子ペアのセットを記録
    seen_pairings = set()
    matches = []
    
    for probe_match in probe_matches:
        for ref_submol_match in ref_submol_matches:
            # 部分構造のマッチを元の分子のインデックスに変換
            ref_match = tuple(
                submol_to_mol_map.get(idx, idx)
                for idx in ref_submol_match
            )
            
            # 原子ペアのセットを作成（順序に依存しない）
            pairing = frozenset(zip(probe_match, ref_match))
            
            # このペアリングが既に見たことがあるかチェック
            if pairing in seen_pairings:
                continue
            seen_pairings.add(pairing)
            
            # 配列として保存
            atom_pairs = np.array([probe_match, ref_match], dtype=np.int_)
            matches.append(atom_pairs)
    
    return matches


def align_structure(
    protein: Structure,
    ligand_coords: npt.NDArray[np.float_],
    probe_coords: npt.NDArray[np.float_],
    atom_pairs: npt.NDArray[np.int_]
) -> Structure:
    """
    原子ペアに基づいてタンパク質構造を重ね合わせます。
    
    プローブ分子と共結晶化リガンドの対応する原子座標を使用して、
    タンパク質全体を重ね合わせた新しい構造を返します。
    
    Parameters
    ----------
    protein : Structure
        重ね合わせるタンパク質構造（BioPython Structure オブジェクト）
    ligand_coords : np.ndarray
        共結晶化リガンドの座標配列 (shape: (n_atoms, 3))
    probe_coords : np.ndarray
        プローブ分子の座標配列 (shape: (n_atoms, 3))
    atom_pairs : np.ndarray
        原子ペアのインデックス配列 (shape: (2, n_atoms))
        atom_pairs[0] がプローブ側、atom_pairs[1] がリガンド側
    
    Returns
    -------
    Structure
        重ね合わせ後のタンパク質構造（新しいオブジェクト）
    """
    # タンパク質の座標を取得
    protein_coords = PDB.get_attr(protein, "coord")
    
    # マッチした原子の座標を取得
    probe_coords_target = probe_coords[atom_pairs[0]]
    ligand_coords_target = ligand_coords[atom_pairs[1]]
    
    # 重ね合わせを実行
    si = SuperImposer()
    si.fit(ligand_coords_target, probe_coords_target)
    
    # タンパク質に変換を適用
    PDB.set_attr(protein, "coord", si.transform(protein_coords))
    
    return protein


def align_structures(
    protein_file: str,
    ligand_file: str,
    probe_files: Union[str, List[str], Dict[str, str]],
    output_dir: str,
    iso_value: int = 1
) -> List[AlignmentResult]:
    """
    プローブ分子とリガンドのマッチングに基づいてタンパク質を重ね合わせます。
    
    この関数は、atom matchingとstructure superimpositionを統合した
    高レベルAPIです。複数のプローブに対してバッチ処理を行います。
    
    Parameters
    ----------
    protein_file : str
        タンパク質PDBファイルのパス
    ligand_file : str
        共結晶化リガンドSDFファイルのパス
    probe_files : str, List[str], or Dict[str, str]
        プローブPDBファイルのパス。以下の形式を受け付けます：
        - str: 単一のプローブファイルパス
        - List[str]: プローブファイルパスのリスト
        - Dict[str, str]: {probe_id: filepath} の辞書
    output_dir : str
        出力ディレクトリのパス
    iso_value : int, optional
        リファレンス分子のアイソトープラベル値（デフォルト: 1）
        0の場合、すべての原子を考慮します
    
    Returns
    -------
    List[AlignmentResult]
        重ね合わせ結果のリスト。各結果には以下が含まれます：
        - probe_id: プローブ分子のID
        - match_id: マッチングID
        - aligned_protein: 重ね合わせ後のタンパク質構造
        - atom_pairs: 原子ペア情報
    
    Examples
    --------
    >>> from inverse_msmd import align_structures
    >>> 
    >>> # 単一プローブの場合
    >>> results = align_structures(
    ...     protein_file="protein.pdb",
    ...     ligand_file="ligand.sdf",
    ...     probe_files="probe.pdb",
    ...     output_dir="./output"
    ... )
    >>> 
    >>> # 複数プローブの場合（辞書形式）
    >>> results = align_structures(
    ...     protein_file="protein.pdb",
    ...     ligand_file="ligand.sdf",
    ...     probe_files={
    ...         "A08": "data/probes/A08.pdb",
    ...         "E24": "data/probes/E24.pdb"
    ...     },
    ...     output_dir="./output"
    ... )
    """
    # 出力ディレクトリの作成
    output_path = Path(expandpath(output_dir))
    output_path.mkdir(parents=True, exist_ok=True)
    
    # プローブファイルを辞書形式に正規化
    if isinstance(probe_files, str):
        # 単一ファイルの場合
        probe_id = Path(probe_files).stem
        probe_dict = {probe_id: probe_files}
    elif isinstance(probe_files, list):
        # リストの場合
        probe_dict = {Path(f).stem: f for f in probe_files}
    else:
        # 辞書の場合
        probe_dict = probe_files
    
    # 共結晶化リガンドを読み込み
    ligand_file = expandpath(ligand_file)
    ref_lig_mol = next(Chem.SDMolSupplier(ligand_file))
    ligand_coords = ref_lig_mol.GetConformer().GetPositions()
    
    results = []
    
    # 各プローブについて処理
    for probe_id, probe_file in probe_dict.items():
        probe_file = expandpath(probe_file)
        
        # プローブ分子を読み込み
        probe_mol = Chem.MolFromPDBFile(probe_file)
        if probe_mol is None:
            print(f"警告: プローブファイル {probe_file} を読み込めませんでした")
            continue
        
        probe_coords = probe_mol.GetConformer().GetPositions()
        
        # 原子マッチングを実行
        matches = find_atom_matches(probe_mol, ref_lig_mol, iso_value)
        
        if not matches:
            print(f"警告: {probe_id} のマッチングが見つかりませんでした")
            continue
        
        print(f"{probe_id}: {len(matches)} 個のマッチングを発見")
        
        # 各マッチングについて重ね合わせを実行
        for match_id, atom_pairs in enumerate(matches):
            # タンパク質を読み込み（各マッチごとに新しいインスタンス）
            protein = PDB.get_structure(expandpath(protein_file))
            
            # 重ね合わせを実行
            aligned_protein = align_structure(
                protein,
                ligand_coords,
                probe_coords,
                atom_pairs
            )
            
            # 結果を保存
            output_file = output_path / f"aligned_to_{probe_id}_{match_id}.pdb"
            PDB.save(aligned_protein, str(output_file))
            
            # 結果をリストに追加
            result = AlignmentResult(
                probe_id=probe_id,
                match_id=match_id,
                aligned_protein=aligned_protein,
                atom_pairs=atom_pairs
            )
            results.append(result)
            
            print(f"  マッチ {match_id}: {output_file.name} に保存")
    
    return results


__all__ = [
    "AlignmentResult",
    "find_atom_matches",
    "align_structure",
    "align_structures",
]