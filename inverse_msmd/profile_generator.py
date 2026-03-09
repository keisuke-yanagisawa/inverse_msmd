"""
プロファイル生成モジュール

MSMDトラジェクトリから残基相互作用プロファイル（DXファイル）を生成する。

パイプライン:
    1. create_single_profile: 個別トラジェクトリからプローブ周辺の残基密度を計算
    2. combine_profiles: 複数プロファイルを統合し、バルク比で正規化
    3. compress_profile: DXファイルをgzip圧縮

依存ライブラリ:
    - pytraj (Step 1のみ): CPPTRAJ Python wrapper。conda install -c conda-forge pytraj
    - gridData: DXファイルの読み書き
"""

import gzip
import logging
import shutil
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
from gridData import Grid
from tqdm import tqdm

logger = logging.getLogger(__name__)


# =============================================================================
# Step 1: 個別プロファイル生成 (pytraj依存)
# =============================================================================

def _make_grid(grid_data: np.ndarray, size: int, pitch: float, origin: List[float]) -> Grid:
    """numpy配列からgridData.Gridオブジェクトを作成する。

    Args:
        grid_data: 3次元のnumpy配列
        size: グリッドの座標数（各軸同一の立方体を仮定）
        pitch: グリッドの間隔（Å）
        origin: グリッドの開始点座標 [x, y, z]
    """
    g = Grid()
    g.delta = [pitch, pitch, pitch]
    g.origin = origin
    g.n = [size, size, size]
    g.grid = grid_data
    return g


def _preprocess_trajectory(pt, traj, ref_probe, probe_resi: int):
    """トラジェクトリをプローブ中心に変換する。

    プローブの向きを参照構造に合わせ、プローブを原点に移動し、
    周期境界条件を処理してタンパク質原子をプローブ近傍に配置する。

    Args:
        pt: pytraj module
        traj: pytraj trajectory
        ref_probe: 参照プローブ構造
        probe_resi: プローブ残基のインデックス
    """
    tmp = pt.strip(traj, ":WAT|:HOH|@/H")
    pt.fiximagedbonds(tmp)
    pt.rmsd(tmp, mask=f":{probe_resi}&!@VIS", ref=ref_probe, ref_mask=":*&!@/H")
    tmp = pt.center(tmp, mask=f":{probe_resi}", mass=True, center="origin")
    tmp = pt.image(tmp, "origin center byatom")
    return tmp


def create_single_profile(
    trajectory_path: str,
    topology_path: str,
    ref_probe_path: str,
    probe_resname: str,
    amino_acid: str,
    output_dx: str,
    output_bulk_cnt: str,
    grid_size: int = 100,
    grid_pitch: float = 1.0,
    use_latter_half: bool = True,
) -> Tuple[Optional[str], Optional[str]]:
    """単一トラジェクトリから残基相互作用プロファイルを生成する。

    MDトラジェクトリを読み込み、各プローブ分子について:
    1. プローブを参照構造に重ね合わせ・原点に移動
    2. 指定アミノ酸のCβ原子の3D密度をグリッド上に計算
    3. バルクカウント（ボックス体積あたりの存在数）を計算

    Args:
        trajectory_path: トラジェクトリファイルのパス (.xtc)
        topology_path: トポロジーファイルのパス (.pdb)
        ref_probe_path: 参照プローブ構造のパス (.pdb)
        probe_resname: プローブの残基名 (例: "A17", "E14")
        amino_acid: 対象アミノ酸の3文字コード (例: "ALA", "LEU")
        output_dx: 出力DXファイルのパス
        output_bulk_cnt: バルクカウント出力ファイルのパス
        grid_size: グリッドの座標数 (デフォルト: 100)
        grid_pitch: グリッドの間隔 Å (デフォルト: 1.0)
        use_latter_half: トラジェクトリの後半のみ使用 (デフォルト: True)

    Returns:
        (output_dx, output_bulk_cnt) のタプル。エラー時は (None, None)。

    Raises:
        ImportError: pytrajがインストールされていない場合
    """
    try:
        import pytraj as pt
    except ImportError:
        raise ImportError(
            "pytrajが必要です。conda install -c conda-forge pytraj でインストールしてください。"
        )

    # トラジェクトリの読み込み
    traj = pt.iterload(trajectory_path, topology_path)
    if use_latter_half:
        traj = traj[traj.n_frames // 2:]

    # ヒスチジンの変異体を考慮
    aa_sele = "HIE,HID,HIP" if amino_acid == "HIS" else amino_acid

    # 対象アミノ酸がトラジェクトリに存在するか確認
    try:
        traj[f"((@/C&@CB)&(!:{probe_resname})&(:{aa_sele}))"]
    except IndexError:
        logger.warning(f"{amino_acid} はこのトラジェクトリに存在しません。スキップします。")
        return None, None

    # 参照プローブの読み込み
    ref_probe = pt.load(ref_probe_path)

    box_volume = np.prod(traj.unitcells[0][:3])

    # プローブ残基のインデックスを取得
    probe_resis = set()
    for residue in traj[f":{probe_resname}"].topology.residues:
        probe_resis.add(residue.index)

    grid_data = None
    bulk_cnt = 0.0
    for probe_resi in tqdm(probe_resis, desc=f"{amino_acid}", leave=False):
        ret = _preprocess_trajectory(pt, traj, ref_probe, probe_resi)
        data = pt.grid(
            ret,
            command=f"{grid_size} {grid_pitch} {grid_size} {grid_pitch} {grid_size} {grid_pitch} "
                    f"((@/C&@CB)&(!:{probe_resname})&(:{aa_sele}))",
            dtype="dataset",
        ).to_ndarray()
        bulk_cnt += data.sum() / box_volume
        if grid_data is None:
            grid_data = data
        else:
            grid_data += data

    # 出力
    origin = [-grid_size * grid_pitch / 2] * 3
    grid = _make_grid(grid_data, grid_size, grid_pitch, origin)

    Path(output_dx).parent.mkdir(parents=True, exist_ok=True)
    grid.export(output_dx)

    Path(output_bulk_cnt).parent.mkdir(parents=True, exist_ok=True)
    with open(output_bulk_cnt, "w") as f:
        f.write(f"{bulk_cnt}")

    logger.info(f"生成完了: {output_dx} (bulk_cnt={bulk_cnt:.4f})")
    return output_dx, output_bulk_cnt


# =============================================================================
# Step 2: プロファイル統合・正規化
# =============================================================================

def combine_profiles(
    profile_dx_paths: List[str],
    bulk_cnt_paths: List[str],
    output_path: str,
    eps: float = 0.1,
) -> str:
    """複数のプロファイルを統合し、バルク比で正規化する。

    各DXファイルの密度を合計し、バルクカウントの合計で割ることで、
    バルク溶媒中での期待値に対する比率に変換する。
    log(0) を回避するため、微小値εを加算する。

    Args:
        profile_dx_paths: 個別プロファイルDXファイルのパスリスト
        bulk_cnt_paths: バルクカウントファイルのパスリスト
        output_path: 出力DXファイルのパス
        eps: log(0) 回避用の微小値 (デフォルト: 0.1)

    Returns:
        出力ファイルのパス
    """
    if len(profile_dx_paths) == 0:
        raise ValueError("プロファイルDXファイルが指定されていません")
    if len(profile_dx_paths) != len(bulk_cnt_paths):
        raise ValueError(
            f"DXファイル数({len(profile_dx_paths)})とバルクカウントファイル数"
            f"({len(bulk_cnt_paths)})が一致しません"
        )

    # グリッドの合計
    logger.info(f"{len(profile_dx_paths)} 個のプロファイルを統合中...")
    grid = Grid(profile_dx_paths[0])
    for path in tqdm(profile_dx_paths[1:], desc="統合", leave=False):
        try:
            tmp = Grid(path)
        except (AttributeError, IOError) as e:
            logger.error(f"DXファイルの読み込みに失敗: {path} - {e}")
            raise
        grid.grid += tmp.grid
    grid.grid += eps

    # バルクカウントの合計
    bulk_total = 0.0
    for path in bulk_cnt_paths:
        val = float(Path(path).read_text().strip())
        bulk_total += val
    if bulk_total == 0:
        raise ValueError("バルクカウントの合計が0です。データに問題があります。")

    # 正規化
    grid.grid /= bulk_total

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    grid.export(output_path)
    logger.info(f"統合完了: {output_path} (bulk_total={bulk_total:.4f})")
    return output_path


# =============================================================================
# Step 3: gzip圧縮
# =============================================================================

def compress_profile(input_path: str, output_path: Optional[str] = None, remove_original: bool = True) -> str:
    """DXファイルをgzip圧縮する。

    Args:
        input_path: 入力DXファイルのパス
        output_path: 出力パス（デフォルト: input_path + ".gz"）
        remove_original: 圧縮後に元ファイルを削除するか

    Returns:
        圧縮ファイルのパス
    """
    input_p = Path(input_path)
    if not input_p.exists():
        raise FileNotFoundError(f"ファイルが見つかりません: {input_path}")

    if output_path is None:
        output_path = str(input_p) + ".gz" if not str(input_p).endswith(".gz") else str(input_p)

    with open(input_p, "rb") as f_in:
        with gzip.open(output_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    if remove_original and str(input_p) != output_path:
        input_p.unlink()

    logger.info(f"圧縮完了: {output_path}")
    return output_path


# =============================================================================
# 統合パイプライン
# =============================================================================

AMINO_ACIDS = [
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE",
    "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER",
    "THR", "VAL", "TRP", "TYR",
]


def _is_valid_file(path: Path) -> bool:
    """ファイルが存在し、サイズが0より大きいか確認する。"""
    return path.exists() and path.stat().st_size > 0


def generate_profiles(
    trajectory_topology_pairs: List[Tuple[str, str]],
    ref_probe_path: str,
    probe_resname: str,
    probe_id: str,
    output_dir: str,
    amino_acids: Optional[List[str]] = None,
    grid_size: int = 100,
    grid_pitch: float = 1.0,
    use_latter_half: bool = True,
    compress: bool = True,
    n_jobs: int = 1,
    eps: float = 0.1,
    continue_on_error: bool = True,
) -> List[str]:
    """トラジェクトリ群からプロファイルを一括生成する統合パイプライン。

    Args:
        trajectory_topology_pairs: (trajectory_path, topology_path) のリスト
        ref_probe_path: 参照プローブPDBのパス
        probe_resname: プローブの残基名 (例: "A17")
        probe_id: 出力ファイル名用のプローブID (例: "A17")
        output_dir: 出力ディレクトリ
        amino_acids: 対象アミノ酸リスト (デフォルト: 全20種)
        grid_size: グリッドサイズ (デフォルト: 100)
        grid_pitch: グリッド間隔 Å (デフォルト: 1.0)
        use_latter_half: トラジェクトリの後半のみ使用
        compress: gzip圧縮するか
        n_jobs: 並列ジョブ数 (Step 1のみ)
        eps: バルク正規化時のεパラメータ
        continue_on_error: エラー時に継続するか (デフォルト: True)

    Returns:
        生成されたプロファイルファイルのパスリスト
    """
    if amino_acids is None:
        amino_acids = AMINO_ACIDS

    out_dir = Path(output_dir)
    single_dir = out_dir / "single_profiles"
    single_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: 個別プロファイル生成
    logger.info(f"Step 1: 個別プロファイル生成 ({len(trajectory_topology_pairs)} トラジェクトリ × {len(amino_acids)} アミノ酸)")

    if n_jobs > 1:
        from joblib import Parallel, delayed

        def _run_single(traj_path, topo_path, aa, sys_idx):
            dx_path = str(single_dir / f"sys{sys_idx}_{probe_id}_{aa}_environment.dx")
            cnt_path = str(single_dir / f"sys{sys_idx}_{probe_id}_{aa}_bulk_cnt.txt")
            if _is_valid_file(Path(dx_path)) and _is_valid_file(Path(cnt_path)):
                logger.info(f"スキップ（既存）: {dx_path}")
                return dx_path, cnt_path
            return create_single_profile(
                traj_path, topo_path, ref_probe_path, probe_resname, aa,
                dx_path, cnt_path, grid_size, grid_pitch, use_latter_half,
            )

        results = Parallel(n_jobs=n_jobs, backend="multiprocessing")(
            delayed(_run_single)(traj, topo, aa, i)
            for i, (traj, topo) in enumerate(trajectory_topology_pairs)
            for aa in amino_acids
        )
    else:
        results = []
        for i, (traj_path, topo_path) in enumerate(trajectory_topology_pairs):
            for aa in amino_acids:
                dx_path = str(single_dir / f"sys{i}_{probe_id}_{aa}_environment.dx")
                cnt_path = str(single_dir / f"sys{i}_{probe_id}_{aa}_bulk_cnt.txt")
                if _is_valid_file(Path(dx_path)) and _is_valid_file(Path(cnt_path)):
                    logger.info(f"スキップ（既存）: {dx_path}")
                    results.append((dx_path, cnt_path))
                    continue
                result = create_single_profile(
                    traj_path, topo_path, ref_probe_path, probe_resname, aa,
                    dx_path, cnt_path, grid_size, grid_pitch, use_latter_half,
                )
                results.append(result)

    # Step 2 & 3: アミノ酸ごとに統合・圧縮
    output_files = []
    for aa in tqdm(amino_acids, desc="Step 2-3: 統合・圧縮"):
        # この AA に関する全ての single profile を収集
        dx_paths = []
        cnt_paths = []
        for i in range(len(trajectory_topology_pairs)):
            dx_path = single_dir / f"sys{i}_{probe_id}_{aa}_environment.dx"
            cnt_path = single_dir / f"sys{i}_{probe_id}_{aa}_bulk_cnt.txt"
            if _is_valid_file(dx_path) and _is_valid_file(cnt_path):
                dx_paths.append(str(dx_path))
                cnt_paths.append(str(cnt_path))

        if len(dx_paths) == 0:
            logger.warning(f"{aa}: プロファイルが見つかりません。スキップします。")
            continue

        try:
            # 統合
            combined_path = str(out_dir / f"{probe_id}_{aa}_profile.dx")
            combine_profiles(dx_paths, cnt_paths, combined_path, eps=eps)

            # 圧縮
            if compress:
                final_path = compress_profile(combined_path)
            else:
                final_path = combined_path

            output_files.append(final_path)
        except Exception as e:
            logger.error(f"{aa}: 統合失敗 - {e}")
            if not continue_on_error:
                raise

    logger.info(f"完了: {len(output_files)} プロファイルを生成しました → {output_dir}")
    return output_files
