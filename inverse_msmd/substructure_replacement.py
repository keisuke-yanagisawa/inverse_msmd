"""
統合部分構造置換モジュール

このモジュールは、リガンド中の部分構造を別の部分構造で置換しつつ、
タンパク質構造も適切に座標変換する統合機能を提供します。

主要な機能:
- リガンド中の部分構造探索と可視化
- 部分構造間のatom matching
- Superimposeによる変換行列計算
- タンパク質への変換適用
- リガンドの部分構造置換
- 統合ワークフロー

使用例:
    >>> from inverse_msmd.substructure_replacement import integrated_substructure_replacement
    >>> 
    >>> results = integrated_substructure_replacement(
    ...     ligand_file="data/atom_matching/4hw3_A_lig.sdf",
    ...     protein_file="data/sample_proteins/4hw3_A.pdb",
    ...     from_file="data/sample_probes/E23",
    ...     to_file="data/sample_probes/E24",
    ...     output_dir="output/integrated/",
    ...     match_index=None
    ... )
    >>> 
    >>> for i, result in enumerate(results):
    ...     print(f"Pattern {i}: {result['ligand_file']}, {result['protein_file']}")
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np
import numpy.typing as npt
from rdkit import Chem
from rdkit.Chem import Draw
from Bio.PDB.Structure import Structure

from .utils.bio_utils import SuperImposer, PDB
from .utils.mol_utils import read_mol_from_pdb_smi


def find_substructure_in_ligand(
    ligand_mol: Chem.Mol,
    substructure_mol: Chem.Mol
) -> List[Tuple[int, ...]]:
    """
    リガンド中の部分構造を探索し、複数のマッチを全て返します。
    
    Parameters
    ----------
    ligand_mol : Chem.Mol
        対象のリガンド分子
    substructure_mol : Chem.Mol
        検索する部分構造
    
    Returns
    -------
    List[Tuple[int, ...]]
        マッチした原子インデックスのタプルのリスト
        各タプルは部分構造の原子に対応するリガンドの原子インデックスを含む
    
    Examples
    --------
    >>> ligand = Chem.SDMolSupplier("ligand.sdf")[0]
    >>> substructure = read_mol_from_pdb_smi("E23.pdb", "E23.smi")
    >>> matches = find_substructure_in_ligand(ligand, substructure)
    >>> print(f"Found {len(matches)} matches")
    """
    # 水素を除いた分子で処理（SDFファイルには通常水素がないため）
    ligand_no_h = Chem.RemoveHs(ligand_mol)
    substructure_no_h = Chem.RemoveHs(substructure_mol)
    
    # 部分構造を検索（全てのマッチを返す）
    matches = ligand_no_h.GetSubstructMatches(substructure_no_h)
    
    # タプルのリストとして返す
    return list(matches)


def visualize_multiple_matches(
    ligand_mol: Chem.Mol,
    substructure_mol: Chem.Mol,
    matches: List[Tuple[int, ...]],
    output_path: str
) -> None:
    """
    複数の部分構造マッチをPNG画像として可視化します。
    
    マッチした部分をハイライト表示し、ユーザーが選択しやすいようにします。
    
    Parameters
    ----------
    ligand_mol : Chem.Mol
        対象のリガンド分子
    substructure_mol : Chem.Mol
        検索された部分構造
    matches : List[Tuple[int, ...]]
        マッチした原子インデックスのリスト
    output_path : str
        出力PNG画像のパス
    
    Examples
    --------
    >>> visualize_multiple_matches(
    ...     ligand_mol,
    ...     substructure_mol,
    ...     matches,
    ...     "output/substructure_matches.png"
    ... )
    """
    from rdkit.Chem import AllChem
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    # 出力ディレクトリを作成
    output_path_obj = Path(output_path)
    output_path_obj.parent.mkdir(parents=True, exist_ok=True)
    
    # マッチ数に応じてレイアウトを決定
    n_matches = len(matches)
    if n_matches == 0:
        return
    
    # グリッドレイアウト（最大4列）
    n_cols = min(4, n_matches)
    n_rows = (n_matches + n_cols - 1) // n_cols
    
    # 図のサイズを設定
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
    if n_matches == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    # リガンドのコピーを作成して2D座標を生成
    ligand_2d = Chem.Mol(ligand_mol)
    AllChem.Compute2DCoords(ligand_2d)
    
    # 各マッチを描画
    for i, match in enumerate(matches):
        # マッチした原子をハイライト
        img = Draw.MolToImage(
            ligand_2d,
            size=(400, 400),
            highlightAtoms=list(match)
        )
        
        axes[i].imshow(img)
        axes[i].axis('off')
        axes[i].set_title(f'Match {i}', fontsize=12, pad=10)
    
    # 使用しない軸を非表示
    for i in range(n_matches, len(axes)):
        axes[i].axis('off')
    
    plt.suptitle(f'Found {n_matches} substructure match(es)', fontsize=14, y=0.98)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def match_substructures(
    mol1: Chem.Mol,
    mol2: Chem.Mol
) -> List[npt.NDArray[np.int_]]:
    """
    2つの部分構造間のatom matchingを実行します。
    
    MCS（最大共通部分構造）検索に基づいて、複数のマッチングパターンを見つけます。
    
    Parameters
    ----------
    mol1 : Chem.Mol
        第一の部分構造
    mol2 : Chem.Mol
        第二の部分構造
    
    Returns
    -------
    List[np.ndarray]
        原子ペアのリスト
        各要素は shape (2, n_atoms) の配列で、
        [0, :] が mol1 のインデックス、[1, :] が mol2 のインデックス
    
    Examples
    --------
    >>> e23_mol = read_mol_from_pdb_smi("E23.pdb", "E23.smi")
    >>> e24_mol = read_mol_from_pdb_smi("E24.pdb", "E24.smi")
    >>> atom_pair_patterns = match_substructures(e23_mol, e24_mol)
    >>> for i, pairs in enumerate(atom_pair_patterns):
    ...     print(f"Pattern {i}: {pairs.shape[1]} atom pairs")
    """
    from rdkit.Chem import rdFMCS
    
    # 水素を除去
    mol1_no_h = Chem.RemoveHs(mol1)
    mol2_no_h = Chem.RemoveHs(mol2)
    
    # MCS検索
    mcs_result = rdFMCS.FindMCS([mol1_no_h, mol2_no_h])
    if mcs_result.numAtoms == 0:
        return []
    
    mcs = Chem.MolFromSmarts(mcs_result.smartsString)
    
    # すべてのマッチングを取得
    mol1_matches = mol1_no_h.GetSubstructMatches(mcs, uniquify=False)
    mol2_matches = mol2_no_h.GetSubstructMatches(mcs, uniquify=False)
    
    # 重複除外のため、原子ペアのセットを記録
    seen_pairings = set()
    matches = []
    
    for mol1_match in mol1_matches:
        for mol2_match in mol2_matches:
            # 原子ペアのセットを作成（順序に依存しない）
            pairing = frozenset(zip(mol1_match, mol2_match))
            
            # このペアリングが既に見たことがあるかチェック
            if pairing in seen_pairings:
                continue
            seen_pairings.add(pairing)
            
            # 配列として保存
            atom_pairs = np.array([mol1_match, mol2_match], dtype=np.int_)
            matches.append(atom_pairs)
    
    return matches


def calculate_transformation(
    source_coords: npt.NDArray[np.float64],
    target_coords: npt.NDArray[np.float64],
    atom_pairs: npt.NDArray[np.int_]
) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """
    Superimposeによる変換行列（回転+並進）を計算します。
    
    Parameters
    ----------
    source_coords : np.ndarray
        変換元の座標配列 (shape: n_atoms x 3)
    target_coords : np.ndarray
        変換先の座標配列 (shape: m_atoms x 3)
    atom_pairs : np.ndarray
        原子ペアのインデックス配列 (shape: 2 x k_pairs)
        atom_pairs[0] が source 側、atom_pairs[1] が target 側のインデックス
    
    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (rot, tran) のタプル
        - rot: 回転行列 (3 x 3)
        - tran: 並進ベクトル (3,)
        
        変換式: new_coords = rot @ coords + tran
    
    Examples
    --------
    >>> rot, tran = calculate_transformation(
    ...     ligand_coords[match_indices],
    ...     e24_coords,
    ...     atom_pairs
    ... )
    >>> transformed = rot @ protein_coords.T + tran[:, np.newaxis]
    """
    # TODO: 実装予定
    pass


def apply_transformation_to_protein(
    protein: Structure,
    rot: npt.NDArray[np.float64],
    tran: npt.NDArray[np.float64]
) -> Structure:
    """
    変換行列をタンパク質構造に適用します。
    
    Parameters
    ----------
    protein : Structure
        変換するタンパク質構造
    rot : np.ndarray
        回転行列 (3 x 3)
    tran : np.ndarray
        並進ベクトル (3,)
    
    Returns
    -------
    Structure
        変換後のタンパク質構造
    
    Examples
    --------
    >>> protein = PDB.get_structure("protein.pdb")
    >>> rot, tran = calculate_transformation(...)
    >>> transformed_protein = apply_transformation_to_protein(protein, rot, tran)
    >>> PDB.save(transformed_protein, "transformed_protein.pdb")
    """
    # TODO: 実装予定
    pass


def replace_ligand_substructure(
    ligand_mol: Chem.Mol,
    match: Tuple[int, ...],
    replacement_mol: Chem.Mol,
    atom_pairs: npt.NDArray[np.int_]
) -> Chem.Mol:
    """
    リガンドの部分構造を新しい部分構造で置換します。
    
    Parameters
    ----------
    ligand_mol : Chem.Mol
        対象のリガンド分子
    match : Tuple[int, ...]
        置換する部分のリガンド原子インデックス
    replacement_mol : Chem.Mol
        置換後の部分構造
    atom_pairs : np.ndarray
        原子対応のインデックス配列 (shape: 2 x n_atoms)
        atom_pairs[0] が元の部分構造、atom_pairs[1] が置換後の部分構造
    
    Returns
    -------
    Chem.Mol
        部分構造が置換された新しいリガンド分子
    
    Examples
    --------
    >>> replaced_ligand = replace_ligand_substructure(
    ...     ligand_mol,
    ...     match,
    ...     e24_mol,
    ...     atom_pairs
    ... )
    >>> Chem.SanitizeMol(replaced_ligand)
    """
    # TODO: 実装予定
    pass


def integrated_substructure_replacement(
    ligand_file: str,
    protein_file: str,
    from_file: str,
    to_file: str,
    output_dir: str,
    match_index: Optional[int] = None
) -> List[Dict[str, str]]:
    """
    統合部分構造置換ワークフローを実行します。
    
    リガンド中の部分構造を別の部分構造で置換し、
    タンパク質構造も適切に座標変換します。
    
    Parameters
    ----------
    ligand_file : str
        リガンドSDFファイルのパス
    protein_file : str
        タンパク質PDBファイルのパス
    from_file : str
        置換前の部分構造のベースパス（拡張子なし）
        .pdbと.smiファイルを自動的に読み込みます
    to_file : str
        置換後の部分構造のベースパス（拡張子なし）
        .pdbと.smiファイルを自動的に読み込みます
    output_dir : str
        出力ディレクトリのパス
    match_index : Optional[int], default=None
        部分構造マッチのインデックス指定（0始まり）
        Noneの場合、複数マッチ時は画像出力してユーザーに選択を促す
    
    Returns
    -------
    List[Dict[str, str]]
        各atom matchingパターンの結果リスト
        各辞書は以下のキーを含む:
        - 'ligand_file': 置換後のリガンドSDFファイルパス
        - 'protein_file': 座標変換後のタンパク質PDBファイルパス
    
    Examples
    --------
    >>> results = integrated_substructure_replacement(
    ...     ligand_file="data/atom_matching/4hw3_A_lig.sdf",
    ...     protein_file="data/sample_proteins/4hw3_A.pdb",
    ...     from_file="data/sample_probes/E23",
    ...     to_file="data/sample_probes/E24",
    ...     output_dir="output/integrated/",
    ...     match_index=0
    ... )
    >>> for i, result in enumerate(results):
    ...     print(f"Pattern {i}:")
    ...     print(f"  Ligand: {result['ligand_file']}")
    ...     print(f"  Protein: {result['protein_file']}")
    """
    # TODO: 実装予定
    pass