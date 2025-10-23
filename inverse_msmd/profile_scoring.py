"""
プロファイルマッチングスコア計算モジュール

このモジュールは、タンパク質構造とプロファイルデータから
マッチングスコアを計算する機能を提供します。

主要な機能:
- 相互作用プロファイルの読み込み
- Cβ原子位置でのプロファイル値の3D補間
- 距離に基づく重み付け（オプション）
- 対数スケールでのスコア統合

使用例:
    >>> from inverse_msmd.profile_scoring import calculate_profile_score
    >>> from inverse_msmd.utils.bio_utils import PDB
    >>> import numpy as np
    >>> 
    >>> protein = PDB.get_structure("data/sample_proteins/4hw3_A.pdb")
    >>> probe_center = np.array([10.0, 15.0, 20.0])
    >>> score = calculate_profile_score(
    ...     protein, probe_center,
    ...     "data/profiles/", "E24", gamma=0.0
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
    gamma: float = 0.0
) -> float:
    """
    タンパク質構造とプローブ中心からマッチングスコアを計算します。
    
    各残基のCβ原子位置での相互作用プロファイル値を3D補間により取得し、
    プローブ中心からの距離で重み付けして統合スコアを算出します。
    
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
    gamma : float, default=0.0
        距離重み付けパラメータ
        - 0.0: 重み付けなし（全残基を均等に扱う）
        - 0.003: 距離に基づく重み付け（ガウシアン減衰）
        値が大きいほど、プローブ中心に近い残基の寄与が大きくなります
    
    Returns
    -------
    float
        対数マッチングスコア
        より大きい（絶対値が小さい）値ほど良いマッチングを示します
    
    Raises
    ------
    ValueError
        プロファイルファイルが見つからない場合
        Cβ原子が見つからない場合
    
    Notes
    -----
    - GLY残基はCβ原子を持たないため、自動的にスキップされます
    - 各原子位置でのプロファイル値は3D線形補間により取得されます
    - 補間値が負の場合、プロファイルの最小値で置き換えられます
    - 距離重み付けはガウシアン関数で計算されます: exp(-gamma * distance^2)
    
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
    >>> # スコア計算（重み付けなし）
    >>> score = calculate_profile_score(
    ...     protein, probe_center,
    ...     "data/profiles/", "E24", gamma=0.0
    ... )
    >>> print(f"スコア: {score:.2f}")
    
    >>> # スコア計算（距離重み付けあり）
    >>> score_weighted = calculate_profile_score(
    ...     protein, probe_center,
    ...     "data/profiles/", "E24", gamma=0.003
    ... )
    >>> print(f"重み付きスコア: {score_weighted:.2f}")
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
    log_score = 0.0
    for atom in atoms_of_interest:
        resname = atom.get_parent().get_resname()
        
        # プロファイルがない残基はスキップ（例: GLY）
        if resname not in profile_residues:
            continue
        
        # 原子座標とプローブ中心からの距離を計算
        coord = atom.get_coord()
        distance = np.linalg.norm(probe_center - coord)
        
        # 距離重み付けを計算（ガウシアン関数）
        weight = np.exp(-gamma * (distance ** 2))
        
        # 3D補間でプロファイル値を取得
        # interpolatedメソッドは各座標をリストで渡す必要がある
        profile_value = profiles[resname].interpolated(
            [coord[0]], [coord[1]], [coord[2]]
        )[0]
        
        # 補間値が負の場合はプロファイルの最小値を使用
        # これにより log(負の値) でエラーが発生するのを防ぐ
        if profile_value < 0:
            profile_value = profiles[resname].grid.min()
        
        # さらに安全のため、非常に小さい値の場合は最小値を使用
        profile_value = max(profile_value, profiles[resname].grid.min())
        
        # 対数スケールで重み付きスコアを累積
        log_score += np.log(profile_value) * weight
    
    return log_score