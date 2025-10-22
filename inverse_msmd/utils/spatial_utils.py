"""
空間計算ユーティリティモジュール

3D空間における幾何学的計算のためのユーティリティ関数を提供します。
"""
from scipy.spatial import KDTree
import numpy as np


def estimate_volume(points, radii, granularity=10):
    """
    球体の集合の体積を推定します。
    
    複数の球体が占める体積を、ボクセルベースの近似により推定します。
    3D空間をグリッドに分割し、各球体の内部にあるボクセルをカウントして
    体積を計算します。
    
    推定誤差は通常約10%以内です（警告が出ない場合）。

    Parameters
    ----------
    points : array-like, shape (n_spheres, 3)
        各球体の中心座標のリスト
    radii : array-like, shape (n_spheres,)
        各球体の半径
    granularity : int, optional
        グリッドの粒度（デフォルト: 10）
        大きいほど精度が上がりますが、計算時間が増加します。
        
    Returns
    -------
    float
        球体の集合が占める体積（Å³）
        
    Warnings
    --------
    - 球体の半径がグリッド間隔の2倍より小さい場合、小さな球体の体積が
      過小評価される可能性があります。この場合、granularityを大きくすることを
      推奨する警告が出力されます。
    - ボクセル数が10^7を超える場合、計算時間が長くなる可能性があります。
    
    Notes
    -----
    この関数は、球体が重なっている場合でも正しく体積を計算します。
    重複部分は一度だけカウントされます。
    
    アルゴリズム:
    1. すべての球体を含む境界ボックスを計算
    2. 境界ボックスを granularity に基づいてグリッド分割
    3. KDTreeを使用して、各球体内にあるボクセルを高速に検索
    4. 占有されているボクセル数をカウント
    5. ボクセルサイズを乗じて体積を計算
    
    Examples
    --------
    >>> points = np.array([[0, 0, 0], [2, 0, 0], [0, 2, 0]])
    >>> radii = np.array([1.0, 1.0, 1.0])
    >>> volume = estimate_volume(points, radii, granularity=20)
    >>> print(f"体積: {volume:.2f} Å³")
    
    >>> # 高精度計算
    >>> volume_high = estimate_volume(points, radii, granularity=50)
    """

    # グリッドを生成
    x_max, y_max, z_max = np.max(points, axis=0)
    x_min, y_min, z_min = np.min(points, axis=0)
    x_pitch = np.min([(x_max - x_min) / granularity, 1])
    y_pitch = np.min([(y_max - y_min) / granularity, 1])
    z_pitch = np.min([(z_max - z_min) / granularity, 1])

    if np.min(radii) < np.max([x_pitch, y_pitch, z_pitch]) * 2:  # 2 はマジックナンバー
        print("警告: 小さな球体の体積が過小評価される可能性があります。")
        print("      granularityを大きくすることを検討してください。")

    # 占有される可能性のある点を列挙
    xs = np.arange(x_min - np.max(radii), x_max + np.max(radii) + 1, x_pitch)
    ys = np.arange(y_min - np.max(radii), y_max + np.max(radii) + 1, y_pitch)
    zs = np.arange(z_min - np.max(radii), z_max + np.max(radii) + 1, z_pitch)
    if len(xs) * len(ys) * len(zs) > 1e7:
        print(f"警告: ボクセル数が非常に多くなっています: {len(xs)*len(ys)*len(zs):.0f}")
        print("      計算に時間がかかる可能性があります。")
    xx, yy, zz = np.meshgrid(xs, ys, zs)
    grid_points = np.vstack((xx.ravel(), yy.ravel(), zz.ravel())).T
    # 高速な最近傍探索のためのKDTreeを作成
    tree = KDTree(grid_points)

    # 占有されているボクセルをカウント
    occupied = np.zeros(len(grid_points), dtype=bool)
    for point, radius in zip(points, radii):
        occupied[tree.query_ball_point(point, radius)] = True
    return occupied.sum() * x_pitch * y_pitch * z_pitch