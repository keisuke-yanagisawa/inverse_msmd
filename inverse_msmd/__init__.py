"""
inverse_msmd - 混合溶媒分子動力学法（MSMD）の逆解析パッケージ

このパッケージは、混合溶媒分子動力学（Mixed-Solvent Molecular Dynamics, MSMD）
シミュレーションの逆解析を行うためのユーティリティを提供します。
タンパク質構造の重ね合わせと、相互作用プロファイルに基づくマッチングスコア計算が
主な機能です。

主要モジュール
------------
alignment : アライメント機能
    統合的な原子マッチングと構造重ね合わせ機能
    - align_structures: 統合アライメント関数（推奨）
    - find_atom_matches: 原子マッチング検索
    - align_structure: 構造重ね合わせ
    - AlignmentResult: 結果を保持するデータクラス

utils.bio_utils : BioPythonユーティリティ
    PDB構造の操作と重ね合わせ
    - SuperImposer: 構造重ね合わせクラス
    - PDB: PDBファイル操作の名前空間

utils.spatial_utils : 空間計算ユーティリティ
    3D幾何学的計算
    - estimate_volume: 球体集合の体積推定

utils.path_utils : パス処理ユーティリティ
    ファイルパスの展開
    - expandpath: 環境変数とチルダを展開

utils.mol_utils : 分子データ操作ユーティリティ
    分子構造データの読み込み、変換、出力機能
    - create_sdf_from_pdb_smi: 異なる形式の構造データを統合

使用例
------
統合APIを使用した基本的なワークフロー：

>>> from inverse_msmd import align_structures
>>> 
>>> # プローブとリガンドのマッチングに基づいてタンパク質を重ね合わせ
>>> results = align_structures(
...     protein_file="data/sample_proteins/4hw3_A.pdb",
...     ligand_file="data/atom_matching/4hw3_A_lig.sdf",
...     probe_files={
...         "A08": "data/sample_probes/A08.pdb",
...         "E24": "data/sample_probes/E24.pdb"
...     },
...     output_dir="./aligned_structures"
... )
>>> 
>>> # 結果の確認
>>> for result in results:
...     print(f"{result.probe_id}_{result.match_id}: "
...           f"{result.atom_pairs.shape[1]} 個の原子がマッチ")

低レベルAPIを使用した詳細な制御：

>>> from inverse_msmd import SuperImposer, PDB
>>> import numpy as np
>>> 
>>> # PDBファイルの読み込み
>>> protein = PDB.get_structure("protein.pdb")
>>> probe = PDB.get_structure("probe.pdb")
>>> 
>>> # 座標の取得
>>> protein_coords = PDB.get_attr(protein, "coord")
>>> probe_coords = PDB.get_attr(probe, "coord")
>>> 
>>> # 原子ペアに基づく座標選択
>>> atom_pairs = np.loadtxt("atom_matching.txt", int)
>>> probe_coords_target = probe_coords[atom_pairs[0]]
>>> protein_coords_target = protein_coords[atom_pairs[1]]
>>> 
>>> # 構造の重ね合わせ
>>> si = SuperImposer()
>>> si.fit(protein_coords_target, probe_coords_target)
>>> transformed_coords = si.transform(protein_coords)
>>> 
>>> # 結果の保存
>>> PDB.set_attr(protein, "coord", transformed_coords)
>>> PDB.save(protein, "aligned_protein.pdb")

参照
----
詳細な使用方法については、以下を参照してください：
- README.md: パッケージの概要とインストール方法
- examples/: サンプルスクリプトと詳細な使用例
- examples/README.md: サンプルの詳細な説明
"""

from .utils.bio_utils import SuperImposer, PDB
from .utils.spatial_utils import estimate_volume
from .utils.path_utils import expandpath
from .utils.mol_utils import read_mol_from_pdb_smi
from .alignment import AlignmentResult, align_structures, find_atom_matches, align_structure

__version__ = "0.1.0"
__author__ = "Keisuke Yanagisawa"
__email__ = "yanagisawa@comp.isct.ac.jp"

__all__ = [
    "SuperImposer",
    "PDB",
    "estimate_volume",
    "expandpath",
    "read_mol_from_pdb_smi",
    "AlignmentResult",
    "align_structures",
    "find_atom_matches",
    "align_structure",
]