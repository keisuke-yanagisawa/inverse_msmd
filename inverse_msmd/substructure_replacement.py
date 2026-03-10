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
from typing import Dict, List, Optional, Tuple, Union
import copy
import numpy as np
import numpy.typing as npt
from rdkit import Chem
from rdkit.Chem import Draw
from Bio.PDB.Structure import Structure
from scipy.spatial.distance import pdist, squareform
import csv

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
    
    # MCS検索（環内と環外のマッチングを制限）
    mcs_result = rdFMCS.FindMCS([mol1_no_h, mol2_no_h], ringMatchesRingOnly=True)
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
    from inverse_msmd.utils.bio_utils import SuperImposer
    
    # atom_pairsに基づいて対応する座標を抽出
    source_matched = source_coords[atom_pairs[0]]
    target_matched = target_coords[atom_pairs[1]]
    
    # SuperImposerで変換行列を計算
    si = SuperImposer()
    si.fit(source_matched, target_matched)
    
    return si.rot_, si.tran_


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
    from inverse_msmd.utils.bio_utils import PDB
    
    # タンパク質の座標を取得
    protein_coords = PDB.get_attr(protein, "coord")
    
    # 変換を適用: new_coords = rot @ coords + tran
    transformed_coords = np.dot(protein_coords, rot) + tran
    
    # 変換後の座標を設定
    PDB.set_attr(protein, "coord", transformed_coords)
    
    return protein


def replace_ligand_substructure(
    ligand_mol: Chem.Mol,
    match: Tuple[int, ...],
    replacement_mol: Chem.Mol,
    atom_pairs: npt.NDArray[np.int_]
) -> Tuple[Chem.Mol, Dict[int, int]]:
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
    Tuple[Chem.Mol, Dict[int, int]]
        - 部分構造が置換された新しいリガンド分子
        - 置換部分のマッピング辞書 {置換後リガンド内インデックス: 元のreplacement_mol内インデックス}
    
    Examples
    --------
    >>> replaced_ligand, replacement_map = replace_ligand_substructure(
    ...     ligand_mol,
    ...     match,
    ...     e24_mol,
    ...     atom_pairs
    ... )
    >>> Chem.SanitizeMol(replaced_ligand)
    """
    import copy
    
    # 水素を除去（RemoveHsは座標を保持する）
    ligand_no_h = Chem.RemoveHs(ligand_mol)
    replacement_no_h = Chem.RemoveHs(replacement_mol)
    
    # 座標を取得
    ligand_coords = ligand_no_h.GetConformer().GetPositions()
    replacement_coords = replacement_no_h.GetConformer().GetPositions()
    
    # マッチした部分構造の原子インデックスをセットに変換
    match_set = set(match)
    
    # 接続点を見つける（MCSマッピング経由で対応付け）
    # atom_pairs[0]はE23のMCS原子、atom_pairs[1]はE24のMCS原子
    connections = []
    
    # E23->E24のMCSマッピングを作成
    e23_to_e24_map = {}
    for e23_idx, e24_idx in zip(atom_pairs[0], atom_pairs[1]):
        e23_to_e24_map[e23_idx] = e24_idx
    
    for i, ligand_atom_idx in enumerate(match):
        atom = ligand_no_h.GetAtomWithIdx(ligand_atom_idx)
        for bond in atom.GetBonds():
            neighbor_idx = bond.GetOtherAtomIdx(ligand_atom_idx)
            if neighbor_idx not in match_set:
                # この原子は置換部分の外側にある
                # match内でのligand_atom_idxの位置を見つける
                match_list = list(match)
                pos_in_match = match_list.index(ligand_atom_idx)
                
                # MCSマッピングを使用して対応するE24原子を見つける
                if pos_in_match in e23_to_e24_map:
                    replacement_atom_idx = e23_to_e24_map[pos_in_match]
                    connections.append((replacement_atom_idx, neighbor_idx, bond.GetBondType()))
    
    # RWMolオブジェクトを作成（編集可能な分子）
    rwmol = Chem.RWMol(copy.deepcopy(ligand_no_h))
    
    # 新しい座標リストを作成（削除前の全原子分）
    new_coords_list = []
    for i in range(ligand_no_h.GetNumAtoms()):
        new_coords_list.append(ligand_coords[i])
    
    # 新しい部分構造の原子を追加
    new_atom_map = {}  # replacement_molの原子インデックス -> rwmolの原子インデックス
    for atom in replacement_no_h.GetAtoms():
        atom_idx = atom.GetIdx()
        new_idx = rwmol.AddAtom(atom)
        new_atom_map[atom_idx] = new_idx
        # 置換後の部分構造の元の座標を追加
        new_coords_list.append(replacement_coords[atom_idx])
    
    # 新しい部分構造内の結合を追加
    for bond in replacement_no_h.GetBonds():
        begin_idx = new_atom_map[bond.GetBeginAtomIdx()]
        end_idx = new_atom_map[bond.GetEndAtomIdx()]
        rwmol.AddBond(begin_idx, end_idx, bond.GetBondType())
    
    # 接続を追加
    for replacement_atom_idx, ligand_neighbor_idx, bond_type in connections:
        new_atom_idx = new_atom_map[replacement_atom_idx]
        rwmol.AddBond(new_atom_idx, ligand_neighbor_idx, bond_type)
    
    # 古い部分構造の原子を削除（逆順で削除してインデックスのずれを防ぐ）
    # 座標も同時に削除
    sorted_match = sorted(match)
    for atom_idx in sorted(match, reverse=True):
        rwmol.RemoveAtom(atom_idx)
        del new_coords_list[atom_idx]
    
    # RWMolをMolに変換
    new_mol = rwmol.GetMol()
    
    # 既存のコンフォーマーを全て削除
    new_mol.RemoveAllConformers()
    
    # 座標をコンフォーマーとして設定
    from rdkit.Geometry import Point3D
    conf = Chem.Conformer(new_mol.GetNumAtoms())
    for i, coords in enumerate(new_coords_list):
        conf.SetAtomPosition(i, Point3D(float(coords[0]), float(coords[1]), float(coords[2])))
    new_mol.AddConformer(conf)
    
    # 置換部分のマッピング辞書を作成
    # new_atom_mapは {replacement_mol内インデックス: 追加時のrwmol内インデックス}
    # 削除後の最終的なインデックスに変換する必要がある
    replacement_map = {}
    for replacement_idx, added_idx in new_atom_map.items():
        # added_idxは削除前のrwmol内インデックス
        # 削除前に、matchに含まれる原子でadded_idxより小さいものの数を数える
        num_deleted_before = sum(1 for m in sorted_match if m < added_idx)
        final_idx = added_idx - num_deleted_before
        replacement_map[final_idx] = replacement_idx
    
    return new_mol, replacement_map
def generate_all_replacement_candidates(
    mol: Chem.Mol,
    from_mol: Chem.Mol,
    to_mol: Chem.Mol,
    match: Tuple[int, ...]
) -> List[Chem.Mol]:
    """
    すべての可能な置換候補を生成します。
    
    1箇所の接続点のみを前提とし、置換先の各原子について
    置換候補を生成します。
    
    Parameters
    ----------
    mol : Chem.Mol
        対象の分子
    from_mol : Chem.Mol
        検索する部分構造（水素除去済み）
    to_mol : Chem.Mol
        置き換える部分構造（水素除去済み）
    match : Tuple[int, ...]
        マッチした原子インデックスのタプル
    
    Returns
    -------
    List[Chem.Mol]
        すべての有効な置換候補のリスト
        無効な候補（Sanitizeに失敗）は含まれない
    
    Raises
    ------
    ValueError
        有効な置換候補が1つも生成できない場合
    
    Examples
    --------
    >>> ligand = next(Chem.SDMolSupplier("ligand.sdf"))
    >>> from_mol = read_mol_from_pdb_smi("E23.pdb", "E23.smi")
    >>> to_mol = read_mol_from_pdb_smi("E24.pdb", "E24.smi")
    >>> from_no_h = Chem.RemoveHs(from_mol)
    >>> to_no_h = Chem.RemoveHs(to_mol)
    >>> ligand_no_h = Chem.RemoveHs(ligand)
    >>> matches = ligand_no_h.GetSubstructMatches(from_no_h)
    >>> candidates = generate_all_replacement_candidates(
    ...     ligand_no_h, from_no_h, to_no_h, matches[0]
    ... )
    >>> print(f"生成された候補数: {len(candidates)}")
    
    Notes
    -----
    - 接続点は自動検出されます
    - 化学的に妥当な候補のみが返されます
    - 各候補はSanitizeを通過しています
    """
    # マッチした部分構造の原子インデックスをセットに変換
    match_set = set(match)
    
    # 接続点を見つける（1箇所のみを想定）
    attachment_info = None  # (neighbor_idx, bond_type)
    
    for atom_idx in match:
        atom = mol.GetAtomWithIdx(atom_idx)
        for bond in atom.GetBonds():
            neighbor_idx = bond.GetOtherAtomIdx(atom_idx)
            if neighbor_idx not in match_set:
                attachment_info = (neighbor_idx, bond.GetBondType())
                break
        if attachment_info:
            break
    
    if attachment_info is None:
        # 接続点がない場合（孤立した部分構造）
        return [create_replacement_molecule(mol, match, to_mol, [])]
    
    neighbor_idx, bond_type = attachment_info
    
    # 新しい部分構造の各原子について置換候補を生成
    n_to_atoms = to_mol.GetNumAtoms()
    
    valid_mols = []
    
    for to_atom_idx in range(n_to_atoms):
        connections = [(to_atom_idx, neighbor_idx, bond_type)]
        
        try:
            new_mol = create_replacement_molecule(mol, match, to_mol, connections)
            Chem.SanitizeMol(new_mol)
            valid_mols.append(new_mol)
        except:
            # 無効な候補はスキップ
            pass
    
    if not valid_mols:
        raise ValueError("有効な置換候補を生成できませんでした")
    
    return valid_mols


def create_replacement_molecule(
    mol: Chem.Mol,
    match: Tuple[int, ...],
    to_mol: Chem.Mol,
    connections: List[Tuple[int, int, Chem.BondType]]
) -> Chem.Mol:
    """
    部分構造を置き換えた新しい分子を作成します。
    
    Parameters
    ----------
    mol : Chem.Mol
        元の分子
    match : Tuple[int, ...]
        マッチした原子インデックス
    to_mol : Chem.Mol
        置き換える部分構造
    connections : List[Tuple[int, int, Chem.BondType]]
        接続情報のリスト
        [(to_atom_idx, mol_neighbor_idx, bond_type), ...]
    
    Returns
    -------
    Chem.Mol
        置換後の分子（元のプロパティを保持）
    
    Examples
    --------
    >>> from rdkit import Chem
    >>> mol = Chem.MolFromSmiles("CCc1ccccc1")
    >>> from_mol = Chem.MolFromSmiles("c1ccccc1")
    >>> to_mol = Chem.MolFromSmiles("C1CCCCC1")
    >>> mol_no_h = Chem.RemoveHs(mol)
    >>> from_no_h = Chem.RemoveHs(from_mol)
    >>> to_no_h = Chem.RemoveHs(to_mol)
    >>> match = mol_no_h.GetSubstructMatch(from_no_h)
    >>> # 接続点を見つける
    >>> connections = [(0, 1, Chem.BondType.SINGLE)]
    >>> new_mol = create_replacement_molecule(mol_no_h, match, to_no_h, connections)
    >>> Chem.SanitizeMol(new_mol)
    >>> print(Chem.MolToSmiles(new_mol))
    
    Notes
    -----
    - 元の分子のプロパティはdeep copyで保持されます
    - 座標情報も保持されます（コンフォーマーがある場合）
    """
    # RWMolオブジェクトを作成（編集可能な分子）
    # deep copyで元のプロパティを保持
    rwmol = Chem.RWMol(copy.deepcopy(mol))
    
    # 新しい部分構造の原子を追加
    new_atom_map = {}  # to_molの原子インデックス -> rwmolの原子インデックス
    for atom in to_mol.GetAtoms():
        new_idx = rwmol.AddAtom(atom)
        new_atom_map[atom.GetIdx()] = new_idx
    
    # 新しい部分構造内の結合を追加
    for bond in to_mol.GetBonds():
        begin_idx = new_atom_map[bond.GetBeginAtomIdx()]
        end_idx = new_atom_map[bond.GetEndAtomIdx()]
        rwmol.AddBond(begin_idx, end_idx, bond.GetBondType())
    
    # 接続を追加
    for to_atom_idx, mol_neighbor_idx, bond_type in connections:
        new_atom_idx = new_atom_map[to_atom_idx]
        rwmol.AddBond(new_atom_idx, mol_neighbor_idx, bond_type)
    
    # 古い部分構造の原子を削除（逆順で削除してインデックスのずれを防ぐ）
    for atom_idx in sorted(match, reverse=True):
        rwmol.RemoveAtom(atom_idx)
    
    # RWMolをMolに変換
    new_mol = rwmol.GetMol()
    
    return new_mol



def check_steric_clash(
    mol: Chem.Mol,
    min_distance: float = 2.0,
    exclude_bonded: bool = True
) -> Tuple[bool, List[Tuple[int, int, float]]]:
    """
    分子内の立体障害（原子間距離が異常に近い）をチェックします。
    
    Parameters
    ----------
    mol : Chem.Mol
        チェックする分子
    min_distance : float, default=2.0
        許容最小原子間距離（Å）。これより近い非結合原子対は立体障害とみなされます
    exclude_bonded : bool, default=True
        結合している原子対を除外するかどうか
    
    Returns
    -------
    is_valid : bool
        立体障害がない場合True、ある場合False
    clashes : List[Tuple[int, int, float]]
        立体障害がある原子対のリスト [(atom_idx1, atom_idx2, distance), ...]
    
    Examples
    --------
    >>> is_valid, clashes = check_steric_clash(mol, min_distance=2.0)
    >>> if not is_valid:
    ...     print(f"立体障害検出: {len(clashes)}箇所")
    ...     for i, j, dist in clashes:
    ...         print(f"  原子{i}-原子{j}: {dist:.2f}Å")
    """
    if mol.GetNumConformers() == 0:
        raise ValueError("分子にコンフォーマー（3D座標）がありません")
    
    # 座標を取得
    coords = mol.GetConformer().GetPositions()
    n_atoms = mol.GetNumAtoms()
    
    # 結合している原子対のセットを作成
    bonded_pairs = set()
    if exclude_bonded:
        for bond in mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            bonded_pairs.add((min(i, j), max(i, j)))
    
    # 全原子間距離を計算
    distances = squareform(pdist(coords))
    
    # 立体障害を検出
    clashes = []
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            # 結合している原子対はスキップ
            if exclude_bonded and (i, j) in bonded_pairs:
                continue
            
            dist = distances[i, j]
            if dist < min_distance:
                clashes.append((i, j, dist))
    
    is_valid = len(clashes) == 0
    return is_valid, clashes



def integrated_substructure_replacement(
    ligand_file: str,
    protein_file: str,
    from_file: str,
    to_file: str,
    output_dir: str,
    match_index: Optional[int] = None,
    profile_dir: Optional[str] = None,
    probe_id: Optional[str] = None,
    csv_output: Optional[str] = None,
    image_output: Optional[str] = None,
    deduplicate_by_smiles: bool = False,
    skip_steric_clash_check: bool = False,
    render_figures: bool = False
) -> List[Dict[str, Union[str, float, int]]]:
    """
    統合部分構造置換ワークフローを実行します。
    
    リガンド中の部分構造を別の部分構造で置換し、
    タンパク質構造も適切に座標変換します。
    オプションでプロファイルマッチングスコアも計算できます。
    
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
    profile_dir : Optional[str], default=None
        プロファイルファイルのディレクトリパス
        指定した場合、各パターンのプロファイルスコアを計算します
        Noneの場合、スコア計算をスキップ（後方互換性維持）
    probe_id : Optional[str], default=None
        プローブID（例: "E24"）
        profile_dirが指定されている場合は必須
        プロファイルファイル名: {probe_id}_{残基名}_profile.dx.gz
    csv_output : Optional[str], default=None
        CSV出力ファイルのパス
        指定した場合、結果をCSVファイルとして保存します
        CSVには以下の列が含まれます：
        - pattern_index: パターン番号
        - score: マッチングスコア（profile_dir指定時のみ）
        - ligand_smiles: 置換後のリガンドのSMILES表記
        - ligand_file: 出力されたリガンドSDFファイルパス
        - protein_file: 出力されたタンパク質PDBファイルパス
    image_output : Optional[str], default=None
        画像出力ファイルのパス
        指定した場合、全パターンの置換後リガンドとスコアを
        グリッド形式で可視化した画像（PNG）を保存します
        スコアが計算されている場合のみ有効
        2D構造は自動的にアライメントされ、向きが統一されます
    deduplicate_by_smiles : bool, default=False
        SMILES表記が同じリガンドの重複を除去するかどうか
        Trueの場合、同じSMILES構造を持つパターンの中で
        最高スコア（最も正の値が大きい）のものだけを保持します
        スコア計算が有効な場合のみ機能します
    skip_steric_clash_check : bool, default=False
        立体障害チェックをスキップするかどうか
        Trueの場合、原子間距離が近すぎるパターンもそのまま出力します
        Falseの場合（デフォルト）、立体障害があるパターンは除外されます
    render_figures : bool, default=False
        3D構造図（PyMOL）を自動生成するかどうか。
        Trueの場合、各パターンの複合体図・統合図を出力ディレクトリに保存します。
        視点はプローブ分子のPCA（compute_probe_view）で自動計算されます。
        PyMOLがインストールされていない場合は警告を表示してスキップします。
    
    Returns
    -------
    List[Dict[str, Union[str, float, int]]]
        各atom matchingパターンの結果リスト（スコア計算時は降順ソート済み）
        各辞書は以下のキーを含む:
        - 'ligand_file': 置換後のリガンドSDFファイルパス
        - 'protein_file': 座標変換後のタンパク質PDBファイルパス
        - 'pattern_index': パターンインデックス
        - 'score': プロファイルスコア（profile_dir指定時のみ）
    
    Raises
    ------
    ValueError
        profile_dirが指定されているがprobe_idがNoneの場合
        その他の入力パラメータが不正な場合
    
    Examples
    --------
    基本的な使用（スコア計算なし）:
    
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
    
    プロファイルスコア計算を含む使用:
    
    >>> results = integrated_substructure_replacement(
    ...     ligand_file="data/atom_matching/4hw3_A_lig.sdf",
    ...     protein_file="data/sample_proteins/4hw3_A.pdb",
    ...     from_file="data/sample_probes/E23",
    ...     to_file="data/sample_probes/E24",
    ...     output_dir="output/integrated/",
    ...     profile_dir="data/profiles/",
    ...     probe_id="E24",
    ...     gamma=0.0
    ... )
    >>> # 結果はスコアで降順ソート済み
    >>> for i, result in enumerate(results):
    ...     print(f"Pattern {i}: スコア={result['score']:.2f}")
    ...     print(f"  Ligand: {result['ligand_file']}")
    """
    from pathlib import Path
    from inverse_msmd.utils.mol_utils import read_mol_from_pdb_smi
    from inverse_msmd.utils.bio_utils import PDB
    
    # パラメータバリデーション
    calculate_scores = profile_dir is not None
    if calculate_scores and probe_id is None:
        raise ValueError(
            "profile_dirが指定されている場合、probe_idも必須です。\n"
            "例: probe_id='E24'"
        )
    
    # 出力ディレクトリを作成
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # ファイルの読み込み
    print(f"ファイルを読み込み中...")
    ligand_mol = next(Chem.SDMolSupplier(ligand_file))
    protein = PDB.get_structure(protein_file)
    from_mol = read_mol_from_pdb_smi(f"{from_file}.pdb", f"{from_file}.smi")
    to_mol = read_mol_from_pdb_smi(f"{to_file}.pdb", f"{to_file}.smi")
    
    print(f"  リガンド原子数: {ligand_mol.GetNumAtoms()}")
    print(f"  タンパク質原子数: {len(PDB.get_attr(protein, 'coord'))}")
    print(f"  置換前部分構造原子数: {from_mol.GetNumAtoms()}")
    print(f"  置換後部分構造原子数: {to_mol.GetNumAtoms()}")
    
    # リガンド中の部分構造を探索
    print(f"\nリガンド中の部分構造を探索中...")
    matches = find_substructure_in_ligand(ligand_mol, from_mol)
    print(f"  マッチ数: {len(matches)}")
    
    if len(matches) == 0:
        raise ValueError("リガンド中に部分構造が見つかりませんでした")
    
    # マッチの選択
    selected_match = None
    if match_index is not None:
        # インデックスが指定されている場合（マッチ数に関わらず検証）
        if match_index < 0 or match_index >= len(matches):
            raise ValueError(
                f"無効なmatch_index: {match_index}。"
                f"有効範囲は 0 から {len(matches)-1} です"
            )
        selected_match = matches[match_index]
        print(f"  指定されたインデックス {match_index} のマッチを使用: {selected_match}")
    elif len(matches) == 1:
        # 1つしかない場合は自動選択
        selected_match = matches[0]
        print(f"  マッチが1つのみのため自動選択: {selected_match}")
    else:
        # 複数マッチがあり、インデックスが指定されていない場合
        print(f"  複数のマッチが見つかりました。可視化画像を生成します...")
        vis_path = output_path / "substructure_matches.png"
        visualize_multiple_matches(ligand_mol, from_mol, matches, str(vis_path))
        print(f"  可視化画像: {vis_path}")
        
        # デフォルトで最初のマッチを使用（ユーザーは画像を見て再実行できる）
        selected_match = matches[0]
        print(f"  警告: 最初のマッチ（インデックス0）を使用します")
        print(f"  他のマッチを使用する場合は、match_indexパラメータを指定して再実行してください")
    
    # Atom matchingを実行
    print(f"\nAtom matchingを実行中...")
    atom_pair_patterns = match_substructures(from_mol, to_mol)
    print(f"  Atom matchingパターン数: {len(atom_pair_patterns)}")
    
    if len(atom_pair_patterns) == 0:
        raise ValueError("Atom matchingに失敗しました")
    
    # 水素を除去してから座標を取得
    ligand_no_h = Chem.RemoveHs(ligand_mol)
    from_no_h = Chem.RemoveHs(from_mol)
    to_no_h = Chem.RemoveHs(to_mol)
    
    # 座標を取得
    ligand_coords = ligand_no_h.GetConformer().GetPositions()
    from_coords = from_no_h.GetConformer().GetPositions()
    to_coords = to_no_h.GetConformer().GetPositions()
    
    # プローブ中心座標を事前計算（スコア計算が必要な場合）
    to_center_e24 = None
    if calculate_scores:
        to_center_e24 = to_no_h.GetConformer().GetPositions().mean(axis=0)
    
    # 各atom matchingパターンについて処理
    results = []
    clashed_patterns = []  # 立体障害でスキップされたパターンを記録
    for pattern_idx, atom_pairs in enumerate(atom_pair_patterns):
        print(f"\nパターン {pattern_idx} を処理中...")
        
        # 1. MCS対応原子の座標を抽出
        # atom_pairs[0]はE23内のインデックス、atom_pairs[1]はE24内のインデックス
        # リガンド中のE23部分のMCS原子の座標を取得
        ligand_mcs_indices = [selected_match[i] for i in atom_pairs[0]]
        ligand_mcs_coords = ligand_coords[ligand_mcs_indices]
        
        # TODO: このコードが、E23, 24専用になっていないか確認
        # E24のMCS原子の座標を取得
        e24_mcs_coords = to_coords[atom_pairs[1]]
        
        # 2. E24をリガンドに合わせる変換を計算（MCS原子のみ使用）
        print(f"  Superimpose計算中（MCS原子: {len(atom_pairs[0])}個）...")
        from inverse_msmd.utils.bio_utils import SuperImposer
        si = SuperImposer()
        si.fit(e24_mcs_coords, ligand_mcs_coords)  # E24をリガンドに合わせる
        rot, tran = si.rot_, si.tran_

        # プローブ中心をE24空間→タンパク質空間に順変換
        if calculate_scores:
            to_center_e24_coord = to_no_h.GetConformer().GetPositions().mean(axis=0)
            transformed_center = np.dot(to_center_e24_coord, rot) + tran

        # 3. E24は変換しない（元の座標のまま使用）
        print(f"  E24を元の座標で使用...")
        import copy
        to_mol_copy = copy.deepcopy(to_no_h)
        
        # 4. リガンドの部分構造をE24で置換
        print(f"  リガンド部分構造を置換中...")
        replaced_ligand, replacement_map = replace_ligand_substructure(
            ligand_no_h,
            selected_match,
            to_mol_copy,
            atom_pairs
        )
        
        # ステップ5-NEW: 置換部分のみ順変換（E24空間→タンパク質空間）
        print(f"  置換部分の座標を修正中...")
        final_coords = replaced_ligand.GetConformer().GetPositions().copy()
        # 非置換部分は既にタンパク質空間にある（変換不要）
        # 置換部分のみE24空間→タンパク質空間に順変換
        for new_ligand_idx, replacement_idx in replacement_map.items():
            e24_coord = to_coords[replacement_idx]
            final_coords[new_ligand_idx] = np.dot(e24_coord, rot) + tran
        
        # 最終座標をリガンドに設定
        replaced_ligand.RemoveAllConformers()
        from rdkit.Geometry import Point3D
        conf = Chem.Conformer(replaced_ligand.GetNumAtoms())
        for i in range(replaced_ligand.GetNumAtoms()):
            conf.SetAtomPosition(i, Point3D(float(final_coords[i, 0]),
                                             float(final_coords[i, 1]),
                                             float(final_coords[i, 2])))
        replaced_ligand.AddConformer(conf)
        print(f"  ✓ {len(replacement_map)}個の原子座標を修正しました")
        
        # 7. 立体障害チェック
        if not skip_steric_clash_check:
            print(f"  立体障害チェック中...")
            is_valid, clashes = check_steric_clash(replaced_ligand, min_distance=2.0)
            
            if not is_valid:
                min_clash_dist = min(d for _, _, d in clashes)
                print(f"  ⚠ 警告: 立体障害を検出しました（{len(clashes)}箇所）")
                for atom_i, atom_j, dist in clashes:
                    atom_i_symbol = replaced_ligand.GetAtomWithIdx(atom_i).GetSymbol()
                    atom_j_symbol = replaced_ligand.GetAtomWithIdx(atom_j).GetSymbol()
                    print(f"    原子{atom_i}({atom_i_symbol}) - 原子{atom_j}({atom_j_symbol}): {dist:.3f}Å")
                clashed_patterns.append({
                    'pattern_idx': pattern_idx,
                    'min_clash_dist': min_clash_dist,
                    'clashes': clashes,
                    'replaced_ligand': replaced_ligand,
                    'rot': rot,
                    'tran': tran,
                    'transformed_protein': copy.deepcopy(protein),
                    'transformed_center': transformed_center if calculate_scores else None,
                })
                print(f"  → このパターンをスキップします")
                continue
            
            print(f"  ✓ 立体障害なし")
        else:
            print(f"  立体障害チェックをスキップ")
        
        # ステップ8-NEW: タンパク質は元座標のまま（変換不要）
        protein_copy = copy.deepcopy(protein)
        transformed_protein = protein_copy
        
        # 9. ファイルを出力
        ligand_output = output_path / f"pattern_{pattern_idx}_ligand_replaced.sdf"
        protein_output = output_path / f"pattern_{pattern_idx}_protein_aligned.pdb"
        
        print(f"  出力中...")
        # リガンドを保存
        writer = Chem.SDWriter(str(ligand_output))
        writer.SetKekulize(False)
        writer.write(replaced_ligand)
        writer.close()
        
        # タンパク質を保存
        PDB.save(transformed_protein, str(protein_output))
        
        print(f"  ✓ リガンド: {ligand_output}")
        print(f"  ✓ タンパク質: {protein_output}")
        
        # リガンドのSMILES表記を取得
        ligand_smiles = Chem.MolToSmiles(replaced_ligand)
        
        # 結果辞書を作成
        result = {
            'ligand_file': str(ligand_output),
            'protein_file': str(protein_output),
            'pattern_index': pattern_idx,
            'ligand_smiles': ligand_smiles,
            'transform_rot': rot,    # 追加: E24→タンパク質空間の回転行列
            'transform_tran': tran,  # 追加: E24→タンパク質空間の並進ベクトル
        }
        
        # プロファイルスコア計算（オプション）
        if calculate_scores:
            print(f"  プロファイルスコア計算中...")
            from .profile_scoring import calculate_profile_score
            
            # スコアを計算
            try:
                score = calculate_profile_score(
                    transformed_protein,
                    transformed_center,
                    profile_dir,
                    probe_id
                )
                result['score'] = score
                print(f"  ✓ スコア: {score:.2f}")
            except Exception as e:
                print(f"  ⚠ 警告: スコア計算に失敗しました: {e}")
                # スコア計算失敗時もパターンは保存
        
        results.append(result)

    # 全パターンが立体障害でスキップされた場合、最も衝突距離が遠いものを採用
    if not results and clashed_patterns:
        best_clash = max(clashed_patterns, key=lambda x: x['min_clash_dist'])
        pat_idx = best_clash['pattern_idx']
        print(f"\n⚠ 警告: 全パターンが立体障害でスキップされました。")
        print(f"  最も衝突距離が遠いパターン {pat_idx}（最小衝突距離: {best_clash['min_clash_dist']:.3f}Å）を採用します。")
        print(f"  ⚠ このパターンのスコアは信頼性が低い可能性があります。")

        replaced_ligand = best_clash['replaced_ligand']
        rot = best_clash['rot']
        tran = best_clash['tran']
        transformed_protein = best_clash['transformed_protein']
        transformed_center = best_clash['transformed_center']

        ligand_output = output_path / f"pattern_{pat_idx}_ligand_replaced.sdf"
        protein_output = output_path / f"pattern_{pat_idx}_protein_aligned.pdb"
        writer = Chem.SDWriter(str(ligand_output))
        writer.SetKekulize(False)
        writer.write(replaced_ligand)
        writer.close()
        PDB.save(transformed_protein, str(protein_output))

        ligand_smiles = Chem.MolToSmiles(replaced_ligand)
        result = {
            'ligand_file': str(ligand_output),
            'protein_file': str(protein_output),
            'pattern_index': pat_idx,
            'ligand_smiles': ligand_smiles,
            'transform_rot': rot,
            'transform_tran': tran,
            'steric_clash_warning': True,
            'min_clash_distance': best_clash['min_clash_dist'],
        }

        if calculate_scores:
            from .profile_scoring import calculate_profile_score
            try:
                score = calculate_profile_score(
                    transformed_protein, transformed_center,
                    profile_dir, probe_id
                )
                result['score'] = score
                print(f"  スコア: {score:.2f} ⚠（立体障害あり、信頼性低）")
            except Exception as e:
                print(f"  ⚠ 警告: スコア計算に失敗しました: {e}")

        results.append(result)

    # スコア計算時は降順ソート
    if calculate_scores and results:
        # スコアがあるもののみソート
        results_with_score = [r for r in results if 'score' in r]
        if results_with_score:
            results_with_score.sort(key=lambda x: x['score'], reverse=True)
            print(f"\n✓ 結果をスコアで降順ソートしました")
            results = results_with_score
            
            # SMILES重複除去（オプション）
            if deduplicate_by_smiles:
                print(f"\nSMILES重複除去を実行中...")
                seen_smiles = {}
                unique_results = []
                
                for result in results:
                    smiles = result.get('ligand_smiles', '')
                    if smiles not in seen_smiles:
                        # 初めて見るSMILES - 保持
                        seen_smiles[smiles] = result
                        unique_results.append(result)
                    # else: すでに見たSMILES - より高スコアのものが既に追加されているのでスキップ
                
                removed_count = len(results) - len(unique_results)
                results = unique_results
                print(f"  ✓ {removed_count} パターンを重複として除去しました")
                print(f"  ✓ 残り {len(results)} パターン（ユニークなSMILES）")
    
    # CSV出力（オプション）
    if csv_output and results:
        print(f"\nCSVファイルを出力中...")
        csv_path = Path(csv_output)
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(csv_path, 'w', newline='', encoding='utf-8') as f:
            # CSVの列を決定（スコアの有無で変わる）
            fieldnames = ['pattern_index']
            if calculate_scores and 'score' in results[0]:
                fieldnames.append('score')
            fieldnames.extend(['ligand_smiles', 'ligand_file', 'protein_file'])
            
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            
            for result in results:
                # 出力する辞書を作成（fielднamesに含まれるキーのみ）
                row = {key: result.get(key, '') for key in fieldnames}
                writer.writerow(row)
        
        print(f"  ✓ CSV出力: {csv_path}")

    # 画像出力（オプション）
    if image_output and results and calculate_scores:
        # スコアがあるパターンのみを収集
        results_with_score = [r for r in results if 'score' in r]
        
        if results_with_score:
            print(f"\n画像ファイルを生成中...")
            from .utils.visualization_utils import draw_scored_molecules_grid
            
            # 各パターンの置換後リガンド分子を読み込み
            molecules = []
            scores = []
            titles = []
            
            for result in results_with_score:
                # SDFファイルから分子を読み込み
                ligand_path = result['ligand_file']
                mol = next(Chem.SDMolSupplier(ligand_path))
                if mol is not None:
                    molecules.append(mol)
                    scores.append(result['score'])
                    titles.append(f"Pattern {result['pattern_index']}")
            
            if molecules:
                image_path = Path(image_output)
                image_path.parent.mkdir(parents=True, exist_ok=True)
                
                draw_scored_molecules_grid(
                    molecules,
                    scores,
                    str(image_path),
                    titles=titles,
                    max_cols=4,
                    align_molecules=True
                )
                print(f"  ✓ 画像出力: {image_path}")
                print(f"  ✓ {len(molecules)} パターンの構造式を可視化しました")
        else:
            print(f"  ⚠ 警告: スコア情報がないため画像を生成できませんでした")
    
    # 3D構造図の生成（オプション）
    if render_figures and results:
        print(f"\n3D構造図を生成中...")
        try:
            from .pymol_visualization import (
                compute_ligand_pca_view,
                render_complex, render_combined,
                render_probe_with_maps, _find_profile_files
            )

            probe_pdb = f"{to_file}.pdb"

            # プロファイルファイルの検索
            profile_files = {}
            if calculate_scores and profile_dir and probe_id:
                profile_files = _find_profile_files(profile_dir, probe_id)

            # プローブ+マップ図（1回のみ、プロファイルがある場合）
            if profile_files:
                panel_b = str(output_path / "probe_map.png")
                render_probe_with_maps(
                    probe_pdb=probe_pdb,
                    profile_files=profile_files,
                    output_png=panel_b,
                )
                print(f"  ✓ プローブ+マップ: {panel_b}")

            # ベストパターンのみ描画（スコア降順ソート済みの先頭）
            best = results[0]
            pat_idx = best['pattern_index']

            # リガンドPCAで視点を計算（原子分散最大化）
            protein_view = compute_ligand_pca_view(best['ligand_file'])

            # 複合体図（タンパク質+リガンド）
            panel_a = str(output_path / f"pattern_{pat_idx}_complex.png")
            render_complex(
                protein_pdb=best['protein_file'],
                ligand_sdf=best['ligand_file'],
                output_png=panel_a,
                view=protein_view,
                original_ligand_sdf=ligand_file,
            )
            print(f"  ✓ パターン {pat_idx} 複合体: {panel_a}")

            # 統合図（タンパク質+リガンド+プローブ+マップ）
            if profile_files:
                panel_c = str(output_path / f"pattern_{pat_idx}_combined.png")
                render_combined(
                    protein_pdb=best['protein_file'],
                    ligand_sdf=best['ligand_file'],
                    probe_pdb=probe_pdb,
                    profile_files=profile_files,
                    output_png=panel_c,
                    view=protein_view,
                    transform_rot=best['transform_rot'],
                    transform_tran=best['transform_tran'],
                    original_ligand_sdf=ligand_file,
                )
                print(f"  ✓ パターン {pat_idx} 統合図: {panel_c}")

            print(f"  ✓ 3D描画完了")
        except ImportError as e:
            print(f"  ⚠ 警告: PyMOLが利用できないため3D描画をスキップします: {e}")
        except Exception as e:
            print(f"  ⚠ 警告: 3D描画中にエラーが発生しました: {e}")

    print(f"\n完了: {len(results)} パターンの結果を生成しました")
    return results