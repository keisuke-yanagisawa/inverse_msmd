"""
分子構造の可視化ユーティリティ

このモジュールは、RDKitとmatplotlibを使用した分子構造の可視化機能を提供します。
複数の分子の比較表示、部分構造マッチのグリッド表示などの機能があります。
"""

from typing import List, Optional, Tuple, Union
from pathlib import Path
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def draw_molecule_comparison(
    molecules: List[Chem.Mol],
    output_path: Union[str, Path],
    titles: Optional[List[str]] = None,
    highlight_atoms_list: Optional[List[List[int]]] = None,
    image_size: Tuple[int, int] = (400, 400),
    dpi: int = 150
) -> None:
    """
    複数の分子を横並びで比較描画します。
    
    Parameters
    ----------
    molecules : List[Chem.Mol]
        描画する分子のリスト（2個以上推奨）
    output_path : str or Path
        出力PNG画像のパス
    titles : List[str], optional
        各分子のタイトル。Noneの場合は分子名を使用
    highlight_atoms_list : List[List[int]], optional
        各分子でハイライトする原子インデックスのリスト
        Noneの場合はハイライトなし
    image_size : Tuple[int, int], default=(400, 400)
        各分子画像のサイズ (width, height)
    dpi : int, default=150
        出力画像の解像度
    
    Raises
    ------
    ValueError
        molecules が空リストの場合
    FileNotFoundError
        出力ディレクトリが作成できない場合
    
    Examples
    --------
    >>> from rdkit import Chem
    >>> mol1 = Chem.MolFromSmiles("CCO")
    >>> mol2 = Chem.MolFromSmiles("CCN")
    >>> draw_molecule_comparison(
    ...     [mol1, mol2],
    ...     "comparison.png",
    ...     titles=["エタノール", "エチルアミン"]
    ... )
    
    Notes
    -----
    - 自動的に2D座標を生成します
    - matplotlibのバックエンドはAggに設定されます
    - 出力ディレクトリが存在しない場合は自動作成されます
    """
    # 入力検証
    if not molecules:
        raise ValueError("molecules リストは空にできません")
    if titles is not None and len(titles) != len(molecules):
        raise ValueError(f"titlesの数({len(titles)})が分子数({len(molecules)})と一致しません")
    if highlight_atoms_list is not None and len(highlight_atoms_list) != len(molecules):
        raise ValueError(
            f"highlight_atoms_listの数({len(highlight_atoms_list)})が分子数({len(molecules)})と一致しません"
        )
    
    # 2D座標生成
    mol_2d_list = []
    for mol in molecules:
        mol_2d = Chem.Mol(mol)  # コピーを作成
        AllChem.Compute2DCoords(mol_2d)
        mol_2d_list.append(mol_2d)
    
    # 描画レイアウト
    n_mols = len(molecules)
    fig, axes = plt.subplots(
        1, n_mols,
        figsize=(image_size[0]/100 * n_mols, image_size[1]/100)
    )
    if n_mols == 1:
        axes = [axes]
    
    # 各分子の描画
    for i, (mol_2d, ax) in enumerate(zip(mol_2d_list, axes)):
        highlight_atoms = highlight_atoms_list[i] if highlight_atoms_list else None
        img = Draw.MolToImage(
            mol_2d,
            size=image_size,
            highlightAtoms=highlight_atoms
        )
        ax.imshow(img)
        ax.axis('off')
        
        # タイトル設定
        if titles:
            title = titles[i]
        elif mol_2d.HasProp("_Name"):
            title = mol_2d.GetProp("_Name")
        else:
            title = f"Molecule {i+1}"
        ax.set_title(title, fontsize=12, pad=10)
    
    # 保存処理
    output_path_obj = Path(output_path)
    output_path_obj.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def draw_molecule_grid(
    molecule: Chem.Mol,
    matches: List[Tuple[int, ...]],
    output_path: Union[str, Path],
    title: str = "Substructure Matches",
    max_cols: int = 4,
    image_size: Tuple[int, int] = (400, 400),
    dpi: int = 150
) -> None:
    """
    同一分子の複数マッチをグリッド形式で描画します。
    
    Parameters
    ----------
    molecule : Chem.Mol
        対象の分子
    matches : List[Tuple[int, ...]]
        マッチした原子インデックスのリスト
        各要素はタプル（原子インデックスの組）
    output_path : str or Path
        出力PNG画像のパス
    title : str, default="Substructure Matches"
        図全体のタイトル
    max_cols : int, default=4
        グリッドの最大列数
    image_size : Tuple[int, int], default=(400, 400)
        各グリッドセルのサイズ (width, height)
    dpi : int, default=150
        出力画像の解像度
    
    Raises
    ------
    ValueError
        matches が空リストの場合
    
    Examples
    --------
    >>> ligand = Chem.SDMolSupplier("ligand.sdf")[0]
    >>> substructure = Chem.MolFromSmiles("c1ccccc1")
    >>> matches = ligand.GetSubstructMatches(substructure)
    >>> draw_molecule_grid(
    ...     ligand,
    ...     matches,
    ...     "matches.png",
    ...     title="ベンゼン環マッチ"
    ... )
    
    Notes
    -----
    - 各マッチは別のグリッドセルに描画されます
    - マッチした原子はハイライト表示されます
    - グリッドレイアウトは自動的に計算されます
    """
    # 入力検証とレイアウト計算
    if not matches:
        raise ValueError("matches リストは空にできません")
    
    n_matches = len(matches)
    n_cols = min(max_cols, n_matches)
    n_rows = (n_matches + n_cols - 1) // n_cols
    
    # 2D座標生成
    mol_2d = Chem.Mol(molecule)
    AllChem.Compute2DCoords(mol_2d)
    
    # グリッド描画
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
    if n_matches == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    for i, match in enumerate(matches):
        img = Draw.MolToImage(
            mol_2d,
            size=image_size,
            highlightAtoms=list(match)
        )
        axes[i].imshow(img)
        axes[i].axis('off')
        axes[i].set_title(f'Match {i}', fontsize=12, pad=10)
    
    # 未使用の軸を非表示
    for i in range(n_matches, len(axes)):
        axes[i].axis('off')
    
    # タイトルと保存
    plt.suptitle(f'{title} ({n_matches} match(es))', fontsize=14, y=0.98)
    output_path_obj = Path(output_path)
    output_path_obj.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def draw_molecule_with_highlights(
    molecule: Chem.Mol,
    highlight_atoms: Optional[List[int]] = None,
    size: Tuple[int, int] = (400, 400)
) -> Image.Image:
    """
    ハイライト付きで分子を描画します（PIL Imageを返す）。
    
    この関数は主に他の描画関数から呼び出される内部用関数です。
    
    Parameters
    ----------
    molecule : Chem.Mol
        描画する分子
    highlight_atoms : List[int], optional
        ハイライトする原子インデックス
    size : Tuple[int, int], default=(400, 400)
        画像サイズ (width, height)
    
    Returns
    -------
    PIL.Image.Image
        描画された分子画像
    
    Examples
    --------
    >>> mol = Chem.MolFromSmiles("c1ccccc1")
    >>> img = draw_molecule_with_highlights(mol, highlight_atoms=[0, 1, 2])
    >>> img.save("benzene_highlighted.png")
    """
    # 2D座標を持つコピーを作成
    mol_2d = Chem.Mol(molecule)
    AllChem.Compute2DCoords(mol_2d)
    
    # 描画
    img = Draw.MolToImage(
        mol_2d,
        size=size,
        highlightAtoms=highlight_atoms
    )
    


def draw_scored_molecules_grid(
    molecules: List[Chem.Mol],
    scores: List[float],
    output_path: Union[str, Path],
    titles: Optional[List[str]] = None,
    max_cols: int = 4,
    image_size: Tuple[int, int] = (400, 400),
    dpi: int = 150,
    align_molecules: bool = True
) -> None:
    """
    複数の分子をマッチングスコア付きでグリッド形式で描画します。
    
    各分子の下にスコアをキャプションとして表示し、2D構造の向きを
    統一することができます。
    
    Parameters
    ----------
    molecules : List[Chem.Mol]
        描画する分子のリスト
    scores : List[float]
        各分子のマッチングスコア
    output_path : str or Path
        出力PNG画像のパス
    titles : List[str], optional
        各分子の追加タイトル（パターン番号など）
        Noneの場合は "Pattern N" を使用
    max_cols : int, default=4
        グリッドの最大列数
    image_size : Tuple[int, int], default=(400, 400)
        各グリッドセルのサイズ (width, height)
    dpi : int, default=150
        出力画像の解像度
    align_molecules : bool, default=True
        2D構造の向きを統一するかどうか
        Trueの場合、最初の分子を基準に他の分子をアライメント
    
    Raises
    ------
    ValueError
        molecules または scores が空リストの場合
        molecules と scores の長さが異なる場合
    
    Examples
    --------
    >>> from rdkit import Chem
    >>> molecules = [Chem.MolFromSmiles("CCO"), Chem.MolFromSmiles("CCN")]
    >>> scores = [-123.45, -145.67]
    >>> draw_scored_molecules_grid(
    ...     molecules,
    ...     scores,
    ...     "output/scored_molecules.png",
    ...     titles=["Pattern 0", "Pattern 1"]
    ... )
    
    Notes
    -----
    - 自動的に2D座標を生成します
    - align_molecules=Trueの場合、RDKitのAlignMolを使用して
      最初の分子を基準に他の分子をアライメントします
    - スコアは小数点以下2桁で表示されます
    - matplotlibのバックエンドはAggに設定されます
    """
    # 入力検証
    if not molecules:
        raise ValueError("molecules リストは空にできません")
    if not scores:
        raise ValueError("scores リストは空にできません")
    if len(molecules) != len(scores):
        raise ValueError(
            f"moleculesの数({len(molecules)})とscoresの数({len(scores)})が一致しません"
        )
    
    # グリッドレイアウトを計算
    n_mols = len(molecules)
    n_cols = min(max_cols, n_mols)
    n_rows = (n_mols + n_cols - 1) // n_cols
    
    # 2D座標を生成
    mol_2d_list = []
    for mol in molecules:
        mol_2d = Chem.Mol(mol)  # コピーを作成
        AllChem.Compute2DCoords(mol_2d)
        mol_2d_list.append(mol_2d)
    
    # 分子のアライメント（オプション）
    if align_molecules and len(mol_2d_list) > 1:
        # 最初の分子を基準として使用
        reference_mol = mol_2d_list[0]
        
        for i in range(1, len(mol_2d_list)):
            try:
                # MCS（最大共通部分構造）を使用してアライメント
                from rdkit.Chem import rdFMCS
                
                mcs_result = rdFMCS.FindMCS(
                    [reference_mol, mol_2d_list[i]],
                    ringMatchesRingOnly=True,
                    completeRingsOnly=True
                )
                
                if mcs_result.numAtoms > 0:
                    # MCSパターンを取得
                    mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
                    
                    # 各分子でMCSにマッチする原子を取得
                    ref_match = reference_mol.GetSubstructMatch(mcs_mol)
                    mol_match = mol_2d_list[i].GetSubstructMatch(mcs_mol)
                    
                    if ref_match and mol_match:
                        # アライメントを実行
                        from rdkit.Chem import rdMolAlign
                        rdMolAlign.AlignMol(
                            mol_2d_list[i],
                            reference_mol,
                            atomMap=list(zip(mol_match, ref_match))
                        )
            except Exception:
                # アライメントに失敗した場合は元の座標を使用
                pass
    
    # グリッド描画
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4.5 * n_rows))
    if n_mols == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    for i, (mol_2d, score) in enumerate(zip(mol_2d_list, scores)):
        # 分子画像を生成
        img = Draw.MolToImage(
            mol_2d,
            size=image_size
        )
        
        axes[i].imshow(img)
        axes[i].axis('off')
        
        # タイトルとスコアを設定
        if titles and i < len(titles):
            title = titles[i]
        else:
            title = f'Pattern {i}'
        
        # スコアを含むキャプションを作成
        caption = f'{title}\nScore: {score:.2f}'
        axes[i].set_title(caption, fontsize=12, pad=10)
    
    # 未使用の軸を非表示
    for i in range(n_mols, len(axes)):
        axes[i].axis('off')
    
    # 全体タイトル
    plt.suptitle(
        f'Replacement Results ({n_mols} pattern(s))',
        fontsize=14,
        y=0.98
    )
    
    # 保存処理
    output_path_obj = Path(output_path)
    output_path_obj.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    return img