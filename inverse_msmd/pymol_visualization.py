"""
PyMOLによる3D構造可視化モジュール

バッチ処理結果に対して、タンパク質-リガンド複合体、プローブ+プロファイルマップ、
および統合図を自動生成する機能を提供します。

依存パッケージ:
    - pymol-open-source (PyMOL headless)
    - numpy

主要関数:
    - render_complex: タンパク質-リガンド複合体の描画
    - render_probe_with_maps: プローブ+プロファイルマップの描画
    - render_combined: 統合図（複合体+プローブ+マップ）の描画
    - compute_probe_view: プローブ分子のPCA視点を自動計算
    - render_batch_results: バッチ処理結果に対する一括描画

使用例:
    >>> from inverse_msmd.pymol_visualization import render_combined, compute_probe_view
    >>>
    >>> view = compute_probe_view("probe/A38.pdb")
    >>> render_combined(
    ...     protein_pdb="output/pattern_0_protein_aligned.pdb",
    ...     ligand_sdf="output/pattern_0_ligand_replaced.sdf",
    ...     probe_pdb="probe/A38.pdb",
    ...     profile_files={"LEU": "profiles/A38_LEU_profile.dx.gz"},
    ...     output_png="output/combined.png",
    ...     view=view,
    ... )
"""

import gzip
import shutil
import tempfile
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


def _ensure_pymol():
    """PyMOLを初期化して cmd を返す（既に起動済みなら何もしない）"""
    import pymol
    from pymol import cmd
    if not hasattr(pymol, '_inverse_msmd_launched'):
        pymol.finish_launching(["pymol", "-cq"])
        pymol._inverse_msmd_launched = True
    return cmd


def compute_probe_view(
    probe_pdb: str,
    distance: float = 50.0,
    slab_near: float = 20.0,
    slab_far: float = 80.0,
    tilt_deg: float = 45.0,
) -> Tuple[float, ...]:
    """
    プローブ分子の原子座標からPCA視点を自動計算します。

    プローブの平面法線方向にカメラを置き、tilt_deg分だけX軸回転を加えます。

    Parameters
    ----------
    probe_pdb : str
        プローブのPDBファイルパス
    distance : float
        カメラ距離
    slab_near : float
        スラブ前方距離
    slab_far : float
        スラブ後方距離
    tilt_deg : float
        X軸方向の傾斜角度（度）

    Returns
    -------
    Tuple[float, ...]
        PyMOLの set_view に渡す18要素のタプル
    """
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("probe", probe_pdb)
    coords = np.array([atom.get_vector().get_array() for atom in structure.get_atoms()])
    centroid = coords.mean(axis=0)

    # PCA
    centered = coords - centroid
    cov = np.cov(centered.T)
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    # 固有値の大きい順にソート
    idx = np.argsort(eigenvalues)[::-1]
    eigenvectors = eigenvectors[:, idx]

    # PC1, PC2が面内、PC3が法線
    pc1 = eigenvectors[:, 0]
    pc2 = eigenvectors[:, 1]
    normal = eigenvectors[:, 2]

    # 右手系を保証
    if np.dot(np.cross(pc1, pc2), normal) < 0:
        normal = -normal

    # X軸方向に tilt_deg 回転
    tilt_rad = np.radians(tilt_deg)
    cos_t, sin_t = np.cos(tilt_rad), np.sin(tilt_rad)
    rot_x = np.array([
        [1, 0, 0],
        [0, cos_t, -sin_t],
        [0, sin_t, cos_t],
    ])

    # 視線行列: 行 = [right, up, -look]
    # PC1(長軸)を右、PC2を上、法線方向から見下ろす（分子平面を正面視）
    view_matrix = np.array([pc1, pc2, normal])
    view_matrix = rot_x @ view_matrix

    # プローブ座標系の重心（通常は原点近傍）
    center = centroid

    return (
        view_matrix[0, 0], view_matrix[0, 1], view_matrix[0, 2],
        view_matrix[1, 0], view_matrix[1, 1], view_matrix[1, 2],
        view_matrix[2, 0], view_matrix[2, 1], view_matrix[2, 2],
        0.0, 0.0, -distance,
        center[0], center[1], center[2],
        slab_near, slab_far, -20.0,
    )


def compute_protein_view(
    protein_pdb: str,
    ligand_sdf: str = None,
    distance: float = 60.0,
    slab_near: float = 20.0,
    slab_far: float = 80.0,
    tilt_deg: float = 30.0,
) -> Tuple[float, ...]:
    """
    タンパク質構造の視点を計算する。

    ligand_sdfが指定された場合、タンパク質中心→リガンド中心方向を視線とし、
    ポケットを外側から覗く視点を返す。リガンド中心にカメラを向ける。
    ligand_sdfが未指定の場合、CA原子PCAに基づく固定視点を返す。
    """
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", protein_pdb)

    # CA原子座標
    ca_coords = []
    for atom in structure.get_atoms():
        if atom.get_name() == "CA":
            ca_coords.append(atom.get_vector().get_array())

    if len(ca_coords) < 3:
        ca_coords = [atom.get_vector().get_array() for atom in structure.get_atoms()]

    coords = np.array(ca_coords)
    protein_centroid = coords.mean(axis=0)

    if ligand_sdf is not None:
        # ポケット視点: タンパク質中心→リガンド中心を視線方向にする
        from rdkit import Chem
        supplier = Chem.SDMolSupplier(ligand_sdf, removeHs=True)
        lig_mol = next(supplier)
        lig_coords = lig_mol.GetConformer().GetPositions()
        lig_centroid = lig_coords.mean(axis=0)

        # 視線方向: タンパク質中心→リガンド中心（カメラはこの逆方向に置かれる）
        look = lig_centroid - protein_centroid
        look = look / np.linalg.norm(look)

        # up方向: PCA第2主軸をベースにGram-Schmidt直交化
        centered = coords - protein_centroid
        cov = np.cov(centered.T)
        eigenvalues, eigenvectors = np.linalg.eigh(cov)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvectors = eigenvectors[:, idx]

        # PCA符号安定化
        for i in range(3):
            max_abs_idx = np.argmax(np.abs(eigenvectors[:, i]))
            if eigenvectors[max_abs_idx, i] < 0:
                eigenvectors[:, i] = -eigenvectors[:, i]

        up_candidate = eigenvectors[:, 1]  # PC2をup候補に
        # lookと直交化
        up = up_candidate - np.dot(up_candidate, look) * look
        up = up / np.linalg.norm(up)

        right = np.cross(look, up)
        right = right / np.linalg.norm(right)

        # PyMOL view matrix: 行 = [right, up, -look]
        view_matrix = np.array([right, up, -look])

        center = lig_centroid
    else:
        # PCAベースの視点（従来動作）
        centered = coords - protein_centroid
        cov = np.cov(centered.T)
        eigenvalues, eigenvectors = np.linalg.eigh(cov)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvectors = eigenvectors[:, idx]

        for i in range(3):
            max_abs_idx = np.argmax(np.abs(eigenvectors[:, i]))
            if eigenvectors[max_abs_idx, i] < 0:
                eigenvectors[:, i] = -eigenvectors[:, i]

        pc1 = eigenvectors[:, 0]
        pc2 = eigenvectors[:, 1]
        normal = eigenvectors[:, 2]

        if np.dot(np.cross(pc1, pc2), normal) < 0:
            normal = -normal

        tilt_rad = np.radians(tilt_deg)
        cos_t, sin_t = np.cos(tilt_rad), np.sin(tilt_rad)
        rot_x = np.array([
            [1, 0, 0],
            [0, cos_t, -sin_t],
            [0, sin_t, cos_t],
        ])

        view_matrix = np.array([pc1, pc2, normal])
        view_matrix = rot_x @ view_matrix

        center = protein_centroid

    return (
        view_matrix[0, 0], view_matrix[0, 1], view_matrix[0, 2],
        view_matrix[1, 0], view_matrix[1, 1], view_matrix[1, 2],
        view_matrix[2, 0], view_matrix[2, 1], view_matrix[2, 2],
        0.0, 0.0, -distance,
        center[0], center[1], center[2],
        slab_near, slab_far, -20.0,
    )


def _build_pymol_transform(rot, tran):
    """
    SuperImposerのrot, tranからPyMOL用4x4変換行列を生成。

    SuperImposerの順変換: y = x @ rot + tran （行ベクトル規約）
    PyMOLのtransform_object: p' = M @ p （列ベクトル規約）

    変換: M = [[rot^T, tran^T], [0, 0, 0, 1]]
    """
    R = rot.T  # 行ベクトル → 列ベクトル規約変換
    return [
        R[0, 0], R[0, 1], R[0, 2], float(tran[0]),
        R[1, 0], R[1, 1], R[1, 2], float(tran[1]),
        R[2, 0], R[2, 1], R[2, 2], float(tran[2]),
        0.0,     0.0,     0.0,     1.0,
    ]


def _resample_dx_file(dx_path_in, dx_path_out, rot, tran):
    """
    DXファイルのボリュームデータをタンパク質空間に再サンプリングする。

    PyMOLのDXリーダーは軸整列直交グリッドのみ対応するため、
    delta vectorsの回転では正しく読み込めない。
    代わりに、新しい軸整列グリッドを作成し、各グリッド点について
    逆変換でプローブ空間の座標を求め、元データから補間する。

    Parameters
    ----------
    dx_path_in : str
        入力DXファイルパス（プローブ空間）
    dx_path_out : str
        出力DXファイルパス（タンパク質空間、軸整列）
    rot : np.ndarray (3,3)
        行ベクトル規約の回転行列: y = x @ rot + tran
    tran : np.ndarray (3,)
        並進ベクトル
    """
    from scipy.ndimage import map_coordinates
    from gridData import Grid

    g = Grid(dx_path_in)
    data = g.grid  # (nx, ny, nz)
    origin = np.array(g.origin, dtype=np.float64)
    delta = float(g.delta[0])  # isotropic spacing assumed
    shape = np.array(data.shape)

    # プローブ空間の8頂点をタンパク質空間に変換し、バウンディングボックスを求める
    corners_probe = np.array([
        origin + np.array([i * (shape[0] - 1), j * (shape[1] - 1), k * (shape[2] - 1)]) * delta
        for i in (0, 1) for j in (0, 1) for k in (0, 1)
    ])
    corners_protein = corners_probe @ rot + tran

    new_origin = corners_protein.min(axis=0) - delta  # 1グリッド分のマージン
    new_end = corners_protein.max(axis=0) + delta
    new_shape = np.ceil((new_end - new_origin) / delta).astype(int) + 1

    # 新しい軸整列グリッドの各点座標（タンパク質空間）
    xi = np.arange(new_shape[0]) * delta + new_origin[0]
    yi = np.arange(new_shape[1]) * delta + new_origin[1]
    zi = np.arange(new_shape[2]) * delta + new_origin[2]
    gx, gy, gz = np.meshgrid(xi, yi, zi, indexing='ij')
    protein_coords = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)

    # タンパク質空間 → プローブ空間への逆変換: x = (y - tran) @ rot^T
    probe_coords = (protein_coords - tran) @ rot.T

    # プローブ空間でのグリッドインデックス（連続値）
    grid_indices = (probe_coords - origin) / delta

    # scipy map_coordinates で補間（範囲外は0）
    resampled = map_coordinates(
        data,
        grid_indices.T,  # (3, N)
        order=0,  # nearest neighbor interpolation
        mode='constant',
        cval=0.0,
    )
    new_data = resampled.reshape(new_shape)

    # 新しいDXファイルを書き出し
    new_grid = Grid(new_data, origin=new_origin,
                    delta=np.array([delta, delta, delta]))
    new_grid.export(dx_path_out, file_format="DX")


def _common_settings(cmd):
    """共通のレンダリング設定を適用"""
    cmd.bg_color("white")
    cmd.set("depth_cue", 1)
    cmd.set("fog_start", 0.3)
    cmd.set("fog", 0.8)
    cmd.set("ray_opaque_background", 1)
    cmd.set("ray_shadow", 0)
    cmd.set("antialias", 2)


def _show_protein(cmd, obj_name="prot", ligand_obj="lig"):
    """タンパク質の表示設定"""
    cmd.hide("everything", obj_name)
    cmd.show("cartoon", obj_name)
    cmd.color("palecyan", obj_name)
    cmd.set("cartoon_transparency", 0.4, obj_name)

    # リガンド周辺残基
    cmd.select("near_lig", f"{obj_name} and byres ({ligand_obj} around 5.0)")
    cmd.show("sticks", "near_lig")
    cmd.color("gray70", "near_lig and elem C")
    cmd.color("red", "near_lig and elem O")
    cmd.color("blue", "near_lig and elem N")
    cmd.set("stick_radius", 0.1, "near_lig")


def _show_ligand(cmd, obj_name="lig", stick_radius=0.2, sphere_scale=0.25):
    """リガンドの表示設定"""
    cmd.show("sticks", obj_name)
    cmd.color("tv_orange", f"{obj_name} and elem C")
    cmd.color("red", f"{obj_name} and elem O")
    cmd.color("blue", f"{obj_name} and elem N")
    cmd.color("yellow", f"{obj_name} and elem S")
    cmd.color("splitpea", f"{obj_name} and elem Cl")
    cmd.set("stick_radius", stick_radius, obj_name)
    cmd.show("spheres", obj_name)
    cmd.set("sphere_scale", sphere_scale, obj_name)


def _show_probe(cmd, obj_name="probe", stick_radius=0.2, sphere_scale=0.3):
    """プローブの表示設定"""
    cmd.show("sticks", obj_name)
    cmd.color("tv_green", f"{obj_name} and elem C")
    cmd.color("splitpea", f"{obj_name} and elem Cl")
    cmd.color("red", f"{obj_name} and elem O")
    cmd.color("blue", f"{obj_name} and elem N")
    cmd.color("yellow", f"{obj_name} and elem S")
    cmd.set("stick_radius", stick_radius, obj_name)
    cmd.show("spheres", obj_name)
    cmd.set("sphere_scale", sphere_scale, obj_name)


def _compute_adaptive_isomesh_level(
    probe_pdb: str,
    profile_files: Dict[str, str],
    radius: float = 15.0,
    percentile: float = 99.95,
    min_level: float = 2.0,
) -> float:
    """
    プローブ近傍のプロファイル値からisomesh閾値を自動計算する。

    Parameters
    ----------
    probe_pdb : str
        プローブPDBファイルパス
    profile_files : Dict[str, str]
        {アミノ酸名: dx.gzファイルパス}
    radius : float
        プローブ原子からの近傍半径（Å）
    percentile : float
        近傍ボクセル値のパーセンタイル（閾値として使用）
    min_level : float
        閾値の最小値

    Returns
    -------
    float
        計算されたisomesh_level
    """
    from Bio.PDB import PDBParser

    # プローブ原子座標
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("probe", probe_pdb)
    probe_coords = np.array([a.get_vector().get_array() for a in structure.get_atoms()])

    nearby_values = []

    for aa, gz_path in profile_files.items():
        # DXファイルを読み込み
        with gzip.open(str(gz_path), "rt") as f:
            origin = None
            delta = np.zeros(3)
            grid_shape = None
            vals = []
            reading = False
            delta_idx = 0
            for line in f:
                if line.startswith("origin"):
                    origin = np.array([float(x) for x in line.split()[1:4]])
                elif line.startswith("delta"):
                    parts = line.split()
                    delta[delta_idx] = float(parts[1 + delta_idx])
                    delta_idx += 1
                elif line.startswith("object 1"):
                    parts = line.split()
                    grid_shape = (int(parts[5]), int(parts[6]), int(parts[7]))
                elif "follows" in line:
                    reading = True
                    continue
                elif reading:
                    if line.startswith("attribute") or line.startswith("object"):
                        break
                    vals.extend(float(x) for x in line.split())

        if origin is None or grid_shape is None or len(vals) == 0:
            continue

        grid = np.array(vals).reshape(grid_shape)

        # 各プローブ原子の近傍ボクセルを収集
        for atom_coord in probe_coords:
            # 原子座標 → グリッドインデックス
            idx_float = (atom_coord - origin) / delta
            # 半径内のインデックス範囲
            r_voxels = int(np.ceil(radius / delta[0]))
            center_idx = np.round(idx_float).astype(int)

            i_min = max(0, center_idx[0] - r_voxels)
            i_max = min(grid_shape[0], center_idx[0] + r_voxels + 1)
            j_min = max(0, center_idx[1] - r_voxels)
            j_max = min(grid_shape[1], center_idx[1] + r_voxels + 1)
            k_min = max(0, center_idx[2] - r_voxels)
            k_max = min(grid_shape[2], center_idx[2] + r_voxels + 1)

            subgrid = grid[i_min:i_max, j_min:j_max, k_min:k_max]

            # 実際の距離フィルタ
            ii, jj, kk = np.mgrid[i_min:i_max, j_min:j_max, k_min:k_max]
            coords = origin + np.stack([ii, jj, kk], axis=-1) * delta
            dist = np.linalg.norm(coords - atom_coord, axis=-1)
            mask = dist <= radius

            values = subgrid[mask]
            finite_values = values[np.isfinite(values) & (values > 1.0)]
            nearby_values.extend(finite_values)

    if len(nearby_values) == 0:
        logger.warning("プローブ近傍に有意なプロファイル値がありません。min_levelを使用します。")
        return min_level

    nearby_values = np.array(nearby_values)
    level = float(np.percentile(nearby_values, percentile))
    level = max(level, min_level)
    logger.info(f"適応的isomesh閾値: {level:.1f} (p{percentile}, {len(nearby_values)} voxels)")
    return level


def _load_profiles(cmd, profile_files, tmpdir, isomesh_level=9.0,
                   transform_rot=None, transform_tran=None):
    """
    プロファイルマップを読み込んでisomesh表示

    Parameters
    ----------
    profile_files : Dict[str, str]
        {アミノ酸名: dx.gzファイルパス} の辞書
    tmpdir : str
        一時ファイルの展開先
    isomesh_level : float
        isomeshの等値面レベル
    transform_rot : np.ndarray (3,3), optional
        行ベクトル規約の回転行列。指定時にDX座標系を変換。
    transform_tran : np.ndarray (3,), optional
        並進ベクトル
    """
    # アミノ酸 → 性質に基づく5色マッピング
    # 疎水性=tv_green, 親水性=cyan, 酸性=tv_red, 塩基性=tv_blue, 芳香性=gray60
    default_colors = {
        # 疎水性（緑）
        "ALA": "tv_green", "VAL": "tv_green", "LEU": "tv_green",
        "ILE": "tv_green", "MET": "tv_green", "PRO": "tv_green",
        # 親水性（シアン）
        "SER": "cyan", "THR": "cyan", "ASN": "cyan",
        "GLN": "cyan", "CYS": "cyan", "GLY": "cyan",
        # 酸性（赤）
        "ASP": "tv_red", "GLU": "tv_red",
        # 塩基性（青）
        "LYS": "tv_blue", "ARG": "tv_blue", "HIS": "tv_blue",
        # 芳香性（濃い灰色）
        "PHE": "gray60", "TRP": "gray60", "TYR": "gray60",
    }

    for aa, gz_path in profile_files.items():
        gz_path = str(gz_path)
        dx_path = Path(tmpdir) / f"{aa}_profile.dx"

        with gzip.open(gz_path, "rb") as f_in:
            with open(str(dx_path), "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        # DXデータをタンパク質空間に再サンプリング（指定時）
        # プローブPDBとDXグリッドは同じ座標空間で生成されているため、
        # プローブと同じ変換行列をそのまま適用する
        if transform_rot is not None:
            dx_transformed = Path(tmpdir) / f"{aa}_profile_transformed.dx"
            _resample_dx_file(str(dx_path), str(dx_transformed), transform_rot, transform_tran)
            dx_path = dx_transformed

        obj_name = f"profile_{aa.lower()}"
        cmd.load(str(dx_path), obj_name)

        # PyMOL fails to create map objects from DX files containing inf values
        if obj_name not in cmd.get_names("objects"):
            logger.warning(f"プロファイル {aa} の読み込みに失敗（データにinfが含まれている可能性）、スキップ")
            continue

        mesh_name = f"map_{aa.lower()}"
        cmd.isomesh(mesh_name, obj_name, level=isomesh_level)
        color = default_colors.get(aa.upper(), "gray50")
        cmd.color(color, mesh_name)
        cmd.set("mesh_color", color, mesh_name)
        cmd.set("transparency", 0.4, mesh_name)


def render_complex(
    protein_pdb: str,
    ligand_sdf: str,
    output_png: str,
    view: Optional[Tuple[float, ...]] = None,
    ray_size: Tuple[int, int] = (800, 800),
    dpi: int = 150,
    save_pse: bool = False,
    zoom_target: str = "lig",
    zoom_buffer: float = 8.0,
) -> str:
    """
    タンパク質-リガンド複合体を描画します。

    Parameters
    ----------
    protein_pdb : str
        タンパク質PDBファイルパス
    ligand_sdf : str
        リガンドSDFファイルパス
    output_png : str
        出力PNG画像パス
    view : Tuple[float, ...], optional
        PyMOL視点（18要素タプル）。Noneの場合はリガンドにzoom
    ray_size : Tuple[int, int]
        レイトレーシング画像サイズ
    dpi : int
        出力画像のDPI
    zoom_target : str
        ズーム対象のPyMOLオブジェクト名
    zoom_buffer : float
        ズーム時の余白（Å）

    Returns
    -------
    str
        出力画像のパス
    """
    cmd = _ensure_pymol()
    cmd.reinitialize()
    _common_settings(cmd)

    cmd.load(protein_pdb, "prot")
    cmd.load(ligand_sdf, "lig")
    cmd.remove("hydrogens")

    _show_protein(cmd)
    _show_ligand(cmd)

    if view is not None:
        cmd.set_view(view)
    cmd.zoom(zoom_target, zoom_buffer)

    Path(output_png).parent.mkdir(parents=True, exist_ok=True)
    cmd.ray(*ray_size)
    cmd.png(output_png, dpi=dpi)
    logger.info(f"Saved: {output_png}")

    if save_pse:
        pse_path = str(Path(output_png).with_suffix(".pse"))
        cmd.save(pse_path)
        logger.info(f"Saved: {pse_path}")

    return output_png


def render_probe_with_maps(
    probe_pdb: str,
    profile_files: Dict[str, str],
    output_png: str,
    view: Optional[Tuple[float, ...]] = None,
    isomesh_level: Optional[float] = None,
    ray_size: Tuple[int, int] = (800, 800),
    dpi: int = 150,
    transform_rot=None,
    transform_tran=None,
    save_pse: bool = False,
) -> str:
    """
    プローブ+プロファイルマップを描画します。

    Parameters
    ----------
    probe_pdb : str
        プローブPDBファイルパス
    profile_files : Dict[str, str]
        {アミノ酸名: dx.gzファイルパス} の辞書
    output_png : str
        出力PNG画像パス
    view : Tuple[float, ...], optional
        PyMOL視点。Noneの場合はプローブにzoom
    isomesh_level : float, optional
        isomeshの等値面レベル。Noneの場合はプローブ近傍値から自動計算
    ray_size : Tuple[int, int]
        レイトレーシング画像サイズ
    dpi : int
        出力画像のDPI

    Returns
    -------
    str
        出力画像のパス
    """
    if isomesh_level is None:
        isomesh_level = _compute_adaptive_isomesh_level(probe_pdb, profile_files)

    cmd = _ensure_pymol()
    cmd.reinitialize()
    _common_settings(cmd)

    cmd.load(probe_pdb, "probe")
    cmd.remove("hydrogens")
    _show_probe(cmd)

    tmpdir = tempfile.mkdtemp()
    try:
        _load_profiles(cmd, profile_files, tmpdir, isomesh_level,
                       transform_rot=transform_rot, transform_tran=transform_tran)

        # プローブをタンパク質空間に変換（DXは_load_profiles内で変換済み）
        if transform_rot is not None:
            transform_matrix = _build_pymol_transform(transform_rot, transform_tran)
            cmd.transform_object("probe", transform_matrix)

        if view is not None:
            cmd.set_view(view)
        else:
            # プローブ+isomesh全体が見える視点を自動計算し、斜め視点に回転
            cmd.orient("probe")
            cmd.zoom("probe", 4)
            cmd.turn("x", -30)
            cmd.turn("y", 30)

        Path(output_png).parent.mkdir(parents=True, exist_ok=True)
        cmd.ray(*ray_size)
        cmd.png(output_png, dpi=dpi)
        logger.info(f"Saved: {output_png}")

        if save_pse:
            pse_path = str(Path(output_png).with_suffix(".pse"))
            cmd.save(pse_path)
            logger.info(f"Saved: {pse_path}")
    finally:
        shutil.rmtree(tmpdir)

    return output_png


def render_combined(
    protein_pdb: str,
    ligand_sdf: str,
    probe_pdb: str,
    profile_files: Dict[str, str],
    output_png: str,
    view: Optional[Tuple[float, ...]] = None,
    isomesh_level: Optional[float] = None,
    ray_size: Tuple[int, int] = (800, 800),
    dpi: int = 150,
    transform_rot=None,
    transform_tran=None,
    save_pse: bool = False,
    zoom_target: str = "lig",
    zoom_buffer: float = 8.0,
) -> str:
    """
    統合図（タンパク質+リガンド+プローブ+プロファイルマップ）を描画します。

    Parameters
    ----------
    protein_pdb : str
        変換後タンパク質PDBファイルパス
    ligand_sdf : str
        置換後リガンドSDFファイルパス
    probe_pdb : str
        プローブPDBファイルパス
    profile_files : Dict[str, str]
        {アミノ酸名: dx.gzファイルパス} の辞書
    output_png : str
        出力PNG画像パス
    view : Tuple[float, ...], optional
        PyMOL視点。Noneの場合はリガンドにzoom
    isomesh_level : float, optional
        isomeshの等値面レベル。Noneの場合はプローブ近傍値から自動計算
    ray_size : Tuple[int, int]
        レイトレーシング画像サイズ
    dpi : int
        出力画像のDPI
    zoom_target : str
        ズーム対象のPyMOLオブジェクト名
    zoom_buffer : float
        ズーム時の余白（Å）

    Returns
    -------
    str
        出力画像のパス
    """
    if isomesh_level is None:
        isomesh_level = _compute_adaptive_isomesh_level(probe_pdb, profile_files)

    cmd = _ensure_pymol()
    cmd.reinitialize()
    _common_settings(cmd)

    cmd.load(protein_pdb, "prot")
    cmd.load(ligand_sdf, "lig")
    cmd.load(probe_pdb, "probe")
    cmd.remove("hydrogens")

    _show_protein(cmd)
    _show_ligand(cmd, stick_radius=0.15, sphere_scale=0.2)
    _show_probe(cmd)

    tmpdir = tempfile.mkdtemp()
    try:
        _load_profiles(cmd, profile_files, tmpdir, isomesh_level,
                       transform_rot=transform_rot, transform_tran=transform_tran)

        # プローブをタンパク質空間に変換（DXは_load_profiles内で変換済み）
        if transform_rot is not None:
            transform_matrix = _build_pymol_transform(transform_rot, transform_tran)
            cmd.transform_object("probe", transform_matrix)

        if view is not None:
            cmd.set_view(view)
        cmd.zoom(zoom_target, zoom_buffer)

        Path(output_png).parent.mkdir(parents=True, exist_ok=True)
        cmd.ray(*ray_size)
        cmd.png(output_png, dpi=dpi)
        logger.info(f"Saved: {output_png}")

        if save_pse:
            pse_path = str(Path(output_png).with_suffix(".pse"))
            cmd.save(pse_path)
            logger.info(f"Saved: {pse_path}")
    finally:
        shutil.rmtree(tmpdir)

    return output_png


def _find_profile_files(
    profile_dir: str,
    probe_id: str,
    amino_acids: Optional[List[str]] = None,
) -> Dict[str, str]:
    """
    プロファイルディレクトリからdx.gzファイルを検索します。

    Parameters
    ----------
    profile_dir : str
        プロファイルディレクトリ
    probe_id : str
        プローブID（例: "A38"）
    amino_acids : List[str], optional
        描画するアミノ酸のリスト。Noneの場合は全て

    Returns
    -------
    Dict[str, str]
        {アミノ酸名: dx.gzファイルパス} の辞書
    """
    profile_path = Path(profile_dir)
    result = {}

    for gz_file in sorted(profile_path.glob(f"{probe_id}_*_profile.dx.gz")):
        # ファイル名からアミノ酸名を抽出: A38_LEU_profile.dx.gz → LEU
        parts = gz_file.stem.replace("_profile.dx", "").split("_")
        if len(parts) >= 2:
            aa = parts[1]
            if amino_acids is None or aa in amino_acids:
                result[aa] = str(gz_file)

    return result


def render_batch_results(
    output_base_dir: str,
    probe_pdb: str,
    profile_dir: str,
    probe_id: str,
    profile_amino_acids: Optional[List[str]] = None,
    view: Optional[Tuple[float, ...]] = None,
    protein_view: Optional[Tuple[float, ...]] = None,
    isomesh_level: float = 9.0,
    ray_size: Tuple[int, int] = (800, 800),
    dpi: int = 150,
    job_ids: Optional[List[str]] = None,
    transform_rot=None,
    transform_tran=None,
) -> List[Dict]:
    """
    バッチ処理結果に対して3パネル図を一括生成します。

    各ジョブの出力ディレクトリ内に以下のファイルを生成:
      - panel_a_complex.png (タンパク質+リガンド)
      - panel_b_probe_map.png (プローブ+マップ) ※全ジョブ共通
      - panel_c_combined.png (統合図)

    Parameters
    ----------
    output_base_dir : str
        バッチ処理の出力ベースディレクトリ
    probe_pdb : str
        プローブPDBファイルパス
    profile_dir : str
        プロファイルディレクトリパス
    probe_id : str
        プローブID（例: "A38"）
    profile_amino_acids : List[str], optional
        描画するアミノ酸のリスト（例: ["LEU", "PHE"]）。
        Noneの場合はプロファイルディレクトリ内の全アミノ酸
    view : Tuple[float, ...], optional
        PyMOL視点。Noneの場合はprobe_pdbからPCAで自動計算
    isomesh_level : float
        isomeshの等値面レベル
    ray_size : Tuple[int, int]
        レイトレーシング画像サイズ
    dpi : int
        出力画像のDPI
    job_ids : List[str], optional
        描画対象のジョブIDリスト。Noneの場合は全ジョブ

    Returns
    -------
    List[Dict]
        各ジョブの描画結果情報のリスト
        [{"job_id": ..., "panel_a": ..., "panel_b": ..., "panel_c": ...}, ...]
    """
    output_path = Path(output_base_dir)

    # 視点の自動計算
    if view is None:
        logger.info(f"プローブ {probe_pdb} からPCA視点を計算中...")
        view = compute_probe_view(probe_pdb)

    # ジョブディレクトリの検索（protein_view計算のために先に確定させる）
    if job_ids is not None:
        job_dirs = [Path(output_base_dir) / jid for jid in job_ids]
    else:
        job_dirs = sorted(
            [d for d in Path(output_base_dir).iterdir() if d.is_dir()],
            key=lambda d: d.name,
        )

    if protein_view is None:
        # 最初のジョブのタンパク質PDBを探す（job_dirsの最初の要素を参照）
        import glob
        first_job_dir = job_dirs[0] if job_dirs else None
        if first_job_dir:
            protein_candidates = (
                list(Path(first_job_dir).glob("*protein_aligned.pdb")) +
                list(Path(first_job_dir).glob("*protein.pdb"))
            )
            if protein_candidates:
                # リガンドSDFも検索してポケット視点を計算
                ligand_candidates = list(Path(first_job_dir).glob("*ligand_replaced.sdf"))
                ligand_sdf = str(ligand_candidates[0]) if ligand_candidates else None
                protein_view = compute_protein_view(
                    str(protein_candidates[0]),
                    ligand_sdf=ligand_sdf,
                )
        if protein_view is None:
            protein_view = view  # フォールバック

    # プロファイルファイルの検索
    profile_files = _find_profile_files(profile_dir, probe_id, profile_amino_acids)
    if not profile_files:
        logger.warning(f"プロファイルファイルが見つかりません: {profile_dir}/{probe_id}_*")
    else:
        logger.info(f"プロファイル: {list(profile_files.keys())}")

    # Panel B（プローブ+マップ）は全ジョブ共通なので1回だけ描画
    panel_b_path = str(output_path / "panel_b_probe_map.png")
    if profile_files:
        render_probe_with_maps(
            probe_pdb=probe_pdb,
            profile_files=profile_files,
            output_png=panel_b_path,
            view=view,
            isomesh_level=isomesh_level,
            ray_size=ray_size,
            dpi=dpi,
            transform_rot=transform_rot,
            transform_tran=transform_tran,
        )
    else:
        panel_b_path = None

    results = []

    for job_dir in job_dirs:
        if not job_dir.is_dir():
            logger.warning(f"ジョブディレクトリが見つかりません: {job_dir}")
            continue

        job_id = job_dir.name

        # パターンファイルを検索（pattern_0を優先、results.csvがあればbest_scoreのパターンを使用）
        best_pattern = _find_best_pattern(job_dir)
        if best_pattern is None:
            logger.warning(f"ジョブ {job_id}: パターンファイルが見つかりません")
            continue

        protein_pdb_path = str(job_dir / f"pattern_{best_pattern}_protein_aligned.pdb")
        ligand_sdf_path = str(job_dir / f"pattern_{best_pattern}_ligand_replaced.sdf")

        if not Path(protein_pdb_path).exists() or not Path(ligand_sdf_path).exists():
            logger.warning(f"ジョブ {job_id}: pattern_{best_pattern} のファイルが不完全です")
            continue

        logger.info(f"ジョブ {job_id} (pattern {best_pattern}) を描画中...")

        result_info = {"job_id": job_id, "pattern_index": best_pattern}

        # Panel A: 複合体
        panel_a = str(job_dir / "panel_a_complex.png")
        render_complex(
            protein_pdb=protein_pdb_path,
            ligand_sdf=ligand_sdf_path,
            output_png=panel_a,
            view=protein_view,
            ray_size=ray_size,
            dpi=dpi,
        )
        result_info["panel_a"] = panel_a

        # Panel B: 共通
        result_info["panel_b"] = panel_b_path

        # Panel C: 統合（各ジョブ固有の変換行列を使用）
        if profile_files:
            # ジョブ固有の変換行列を読み込み
            job_rot, job_tran = transform_rot, transform_tran
            transform_json = job_dir / "transform.json"
            if transform_json.exists():
                import json as _json
                import numpy as np
                with open(transform_json) as tf:
                    tdata = _json.load(tf)
                job_rot = np.array(tdata['transform_rot'])
                job_tran = np.array(tdata['transform_tran'])

            panel_c = str(job_dir / "panel_c_combined.png")
            render_combined(
                protein_pdb=protein_pdb_path,
                ligand_sdf=ligand_sdf_path,
                probe_pdb=probe_pdb,
                profile_files=profile_files,
                output_png=panel_c,
                view=protein_view,
                isomesh_level=isomesh_level,
                ray_size=ray_size,
                dpi=dpi,
                transform_rot=job_rot,
                transform_tran=job_tran,
            )
            result_info["panel_c"] = panel_c

        results.append(result_info)

    logger.info(f"描画完了: {len(results)} ジョブ")
    return results


def _find_best_pattern(job_dir: Path) -> Optional[int]:
    """
    ジョブディレクトリからベストパターンのインデックスを取得します。

    results.csvがある場合は最高スコアのパターンを、
    ない場合はpattern_0を返します。
    """
    import csv

    results_csv = job_dir / "results.csv"
    if results_csv.exists():
        best_score = float("-inf")
        best_idx = None
        with open(results_csv, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    score = float(row["score"])
                    idx = int(row["pattern_index"])
                    if score > best_score:
                        best_score = score
                        best_idx = idx
                except (ValueError, KeyError):
                    continue
        if best_idx is not None:
            return best_idx

    # results.csvがない場合、pattern_0を探す
    if (job_dir / "pattern_0_protein_aligned.pdb").exists():
        return 0

    # 任意のパターンを探す
    for f in sorted(job_dir.glob("pattern_*_protein_aligned.pdb")):
        try:
            idx = int(f.stem.split("_")[1])
            return idx
        except (ValueError, IndexError):
            continue

    return None


__all__ = [
    "compute_probe_view",
    "render_complex",
    "render_probe_with_maps",
    "render_combined",
    "render_batch_results",
]
