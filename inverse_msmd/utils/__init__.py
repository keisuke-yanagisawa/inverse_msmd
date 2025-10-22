"""
inverse_msmdパッケージのユーティリティモジュール

このモジュールは、inverse_msmdパッケージで使用される
各種ユーティリティ関数とクラスを提供します。

モジュール構成
------------
bio_utils : BioPythonユーティリティ
    PDB構造の読み込み、操作、保存、および構造重ね合わせ機能
    - SuperImposer: 構造重ね合わせクラス
    - PDB: PDBファイル操作の名前空間クラス

spatial_utils : 空間計算ユーティリティ
    3D空間における幾何学的計算
    - estimate_volume: 球体集合の体積推定

path_utils : パス処理ユーティリティ
    ファイルパスの展開と処理
    - expandpath: 環境変数とチルダを展開

mol_utils : 分子データ操作ユーティリティ
    分子構造データの読み込み、変換、出力機能
    - create_sdf_from_pdb_smi: 異なる形式の構造データを統合
    - create_mol_from_pdb_with_smiles_template: テンプレートマッチング
    - read_smiles: SMILES文字列の読み込み
    - read_pdb_coordinates: 座標データの読み込み

使用例
------
>>> from inverse_msmd.utils import SuperImposer, PDB, estimate_volume, expandpath
>>> 
>>> # PDBファイルの読み込み
>>> protein = PDB.get_structure("protein.pdb")
>>> 
>>> # 構造の重ね合わせ
>>> si = SuperImposer()
>>> si.fit(moving_coords, target_coords)
>>> transformed = si.transform(coords)
>>> 
>>> # 体積推定
>>> volume = estimate_volume(points, radii)
>>> 
>>> # パス展開
>>> full_path = expandpath("~/data/protein.pdb")
"""

from .bio_utils import SuperImposer, PDB
from .spatial_utils import estimate_volume
from .path_utils import expandpath
from .mol_utils import read_mol_from_pdb_smi

__all__ = [
    "SuperImposer",
    "PDB",
    "estimate_volume",
    "expandpath",
    "read_mol_from_pdb_smi",
]