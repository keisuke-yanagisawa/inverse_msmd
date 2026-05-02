"""
プロファイルマッチングスコア計算モジュール

このモジュールは、タンパク質構造とプロファイルデータから
マッチングスコアを計算する機能を提供します。

入力 dx は log 形式の RIprofile (= log(occupancy / bulk_occupancy)) を想定し、
スコアは各 Cβ 原子位置の RIprofile 値の単純和として計算されます
(論文 Eq. (2) と一致)。

主要な機能:
- 相互作用プロファイルの読み込み
- Cβ原子位置でのプロファイル値の3D補間
- RIprofile 値の単純和によるスコア統合

使用例:
    >>> from inverse_msmd.profile_scoring import calculate_profile_score
    >>> from inverse_msmd.utils.bio_utils import PDB
    >>> import numpy as np
    >>> 
    >>> protein = PDB.get_structure("data/sample_proteins/4hw3_A.pdb")
    >>> probe_center = np.array([10.0, 15.0, 20.0])
    >>> score = calculate_profile_score(
    ...     protein, probe_center,
    ...     "data/profiles/", "E24"
    ... )
    >>> print(f"スコア: {score:.2f}")
"""

from gridData import Grid
import numpy as np
import numpy.typing as npt
from pathlib import Path
from typing import Dict
from Bio.PDB.Structure import Structure


def calculate_profile_score(
    protein: Structure,
    probe_center: npt.NDArray[np.float64],
    profile_dir: str,
    probe_id: str,
    inverse_rot: npt.NDArray[np.float64] = None,
    inverse_tran: npt.NDArray[np.float64] = None,
) -> float:
    """
    タンパク質構造とプローブ中心からマッチングスコアを計算します。
    
    各残基のCβ原子位置での相互作用プロファイル値を3D補間により取得し、
    統合スコアを算出します。
    
    Parameters
    ----------
    protein : Bio.PDB.Structure.Structure
        評価対象のタンパク質構造
    probe_center : np.ndarray, shape (3,)
        プローブ分子の中心座標（Å単位）
    profile_dir : str
        プロファイルファイル（.dx.gz形式）のディレクトリパス
    probe_id : str
        プローブID（ファイル名プレフィックス、例: "E24"）
        プロファイルファイル名: {probe_id}_{残基名}_profile.dx.gz
    
    Returns
    -------
    float
        マッチングスコア (RIprofile 値の単純和)。
        より大きい値ほど良いマッチングを示します。

    Raises
    ------
    ValueError
        プロファイルファイルが見つからない場合
        Cβ原子が見つからない場合

    Notes
    -----
    - GLY残基はCβ原子を持たないため、自動的にスキップされます
    - 各原子位置でのプロファイル値は3D線形補間により取得されます
    - グリッド外 (補間結果が NaN) のCβ原子はバルク参照 (寄与 0) として扱われます
    - 入力 dx は log 形式の RIprofile を想定しています
    
    Examples
    --------
    >>> from inverse_msmd.profile_scoring import calculate_profile_score
    >>> from inverse_msmd.utils.bio_utils import PDB
    >>> import numpy as np
    >>> 
    >>> # タンパク質構造を読み込み
    >>> protein = PDB.get_structure("data/sample_proteins/4hw3_A.pdb")
    >>> 
    >>> # プローブ中心座標（例）
    >>> probe_center = np.array([10.0, 15.0, 20.0])
    >>>
    >>> # スコア計算
    >>> score = calculate_profile_score(
    ...     protein, probe_center,
    ...     "data/profiles/", "E24"
    ... )
    >>> print(f"スコア: {score:.2f}")
    """
    # プロファイルディレクトリのパスオブジェクトを作成
    profile_path = Path(profile_dir)
    
    if not profile_path.exists():
        raise ValueError(f"プロファイルディレクトリが見つかりません: {profile_dir}")
    
    # 利用可能な残基タイプのリスト（GLYを除く標準アミノ酸）
    profile_residues = [
        "ALA", "CYS", "GLU", "ASP", "PHE", "HIS", "ILE", "LYS", "LEU",
        "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"
    ]
    
    # 各残基タイプのプロファイルを読み込み
    profiles: Dict[str, Grid] = {}
    for res in profile_residues:
        profile_file = profile_path / f"{probe_id}_{res}_profile.dx.gz"
        if not profile_file.exists():
            raise ValueError(
                f"プロファイルファイルが見つかりません: {profile_file}\n"
                f"プローブID '{probe_id}' のプロファイルファイルが全て揃っているか確認してください"
            )
        profiles[res] = Grid(str(profile_file))
    
    # Cβ原子のみを選択（GLYにはCβがないので自動的に除外される）
    atoms_of_interest = [atom for atom in protein.get_atoms() if atom.get_name() == "CB"]
    
    if len(atoms_of_interest) == 0:
        raise ValueError(
            "Cβ原子が見つかりません。タンパク質構造にCβ原子を持つ残基が含まれているか確認してください"
        )
    
    # マッチングスコアを計算
    score = 0.0
    for atom in atoms_of_interest:
        resname = atom.get_parent().get_resname()

        # プロファイルがない残基はスキップ（例: GLY）
        if resname not in profile_residues:
            continue

        # 原子座標を取得
        coord = atom.get_coord()
        # 案B: タンパク質空間 → プローブ空間に逆変換してからサンプリング
        if inverse_rot is not None:
            coord = np.dot(coord - inverse_tran, inverse_rot.T)

        # 3D補間でRIprofile値 (log値) を取得
        profile_value = profiles[resname].interpolated(
            [coord[0]], [coord[1]], [coord[2]]
        )[0]

        # グリッド外 → 補間結果が NaN になる場合はバルク参照 (log(1)=0) として扱う
        if np.isnan(profile_value):
            continue

        # RIprofile 値はすでに log 変換済みのため単純和で累積
        score += float(profile_value)

    return score