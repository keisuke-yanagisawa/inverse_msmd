"""
マッチングスコア計算スクリプト

このスクリプトは、inverse_msmdパッケージを使用して、
相互作用プロファイルに基づいてマッチングスコアを計算する方法を示します。

概要
----
重ね合わせ済みのタンパク質構造とプローブ分子の相互作用プロファイルを使用して、
タンパク質-プローブ間のマッチングスコアを計算します。各残基のCβ原子位置での
プロファイル値を補間し、プローブ中心からの距離で重み付けして統合スコアを算出します。

処理の流れ
----------
1. 重ね合わせ済みのタンパク質構造を読み込み
2. プローブ分子を読み込み、中心座標を計算
3. タンパク質のCβ原子を選択
4. 各残基タイプの相互作用プロファイルを読み込み
5. 各Cβ原子の位置でプロファイル値を3D補間
6. 距離に基づく重み付けを適用（オプション）
7. 対数スケールでスコアを統合

パラメータ
----------
GAMMA : float
    距離重み付けパラメータ
    - 0.00: 重み付けなし（すべての残基を均等に扱う）
    - 0.003: 距離に基づく重み付け（ガウシアン減衰）
    値が大きいほど、プローブ中心に近い残基の寄与が大きくなります

使用方法
--------
このスクリプトを実行する前に、integrated_alignment.pyを実行して
重ね合わせ済みのタンパク質構造を生成してください。

    $ cd examples
    $ python integrated_alignment.py
    $ python calculate_matching.py

出力
----
各マッチングについて、以下の形式でスコアが出力されます：
    matching_id score

例：
    A08 -245.67
    E24 -189.23

スコアの解釈
------------
- スコアは対数スケールです
- より大きい（絶対値が小さい）スコアほど良いマッチングを示します
- 負の値が大きいほどマッチングが悪いことを意味します

注意事項
--------
- プロファイルファイル（.dx.gz）がすべての残基タイプについて必要です
- GLY残基はCβ原子を持たないため、自動的にスキップされます
- 補間時に負の値が出る場合、プロファイルの最小値で置き換えられます
"""
from inverse_msmd.utils.bio_utils import SuperImposer, PDB
from gridData import Grid
import numpy as np

# パラメータ設定
GAMMA = 0.00  # 距離重み付けパラメータ（0.00=重み付けなし、0.003=距離重み付けあり）

def calculate_matching_score(atoms_of_interest, profile_residues, probe_center, profiles):
    """
    相互作用プロファイルに基づいてマッチングスコアを計算します。
    
    Parameters
    ----------
    atoms_of_interest : list of Bio.PDB.Atom.Atom
        評価対象の原子リスト（通常はCβ原子）
    profile_residues : list of str
        プロファイルが利用可能な残基名のリスト
    probe_center : np.ndarray, shape (3,)
        プローブ分子の中心座標
    profiles : dict
        各残基タイプのGridオブジェクトの辞書
        キー: 残基名（3文字コード）、値: gridData.Grid
        
    Returns
    -------
    float
        対数マッチングスコア
        
    Notes
    -----
    - GLY残基など、profile_residuesに含まれない残基は自動的にスキップされます
    - 各原子位置でのプロファイル値は3D線形補間により取得されます
    - 距離重み付けはガウシアン関数で計算されます: exp(-GAMMA * distance^2)
    - 補間値が負の場合、プロファイルの最小値で置き換えられます
    """
    log_score = 0
    for atom in atoms_of_interest:
        resname = atom.get_parent().get_resname()
        if resname not in profile_residues:
            continue  # GLYなど、プロファイルがない残基はスキップ
        coord = atom.get_coord()
        distance = np.linalg.norm(probe_center - coord)
        weight = np.exp(-GAMMA * (distance**2))
        # 補間値が負の場合はプロファイルの最小値を使用
        profile_value = max(
            profiles[resname].interpolated([coord[0]], [coord[1]], [coord[2]])[0],
            profiles[resname].grid.min()
        )
        log_score += np.log(profile_value) * weight
    return log_score


# メイン処理
if __name__ == "__main__":
    # 評価するマッチングIDのリスト
    matching_ids = "A08 E24 E24_1 E24_2".split(" ")
    
    for matching in matching_ids:
        probe_id = matching[:3]
        
        # 重ね合わせ済みのタンパク質構造を読み込み
        # integrated_alignment.pyで生成されたファイルを使用
        protein = PDB.get_structure(f"./aligned_structures/aligned_to_{matching}.pdb")
        
        # プローブ分子を読み込み、中心座標を計算
        probe = PDB.get_structure(f"../data/sample_probes/{probe_id}.pdb")
        probe_center = PDB.get_attr(probe, "coord").mean(axis=0)
        
        # Cβ原子のみを選択（GLYにはCβがないので自動的に除外される）
        atoms_of_interest = [a for a in protein.get_atoms() if a.get_name() == "CB"]

        # 利用可能な残基タイプのリスト（GLYを除く標準アミノ酸）
        profile_residues = "ALA CYS GLU ASP PHE HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split(" ")
        
        # 各残基タイプのプロファイルを読み込み
        profiles = {res: Grid(f"../data/profiles/{probe_id}_{res}_profile.dx.gz") for res in profile_residues}
        
        # マッチングスコアを計算
        score = calculate_matching_score(atoms_of_interest, profile_residues, probe_center, profiles)
        
        # 結果を出力
        print(matching, f"{score:.2f}")