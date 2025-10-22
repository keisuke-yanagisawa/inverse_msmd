"""
構造重ね合わせスクリプト

このスクリプトは、inverse_msmdパッケージを使用して、
原子マッチングデータに基づいてタンパク質構造を重ね合わせる方法を示します。

設計原則
--------
- タンパク質、共結晶化リガンド、プローブは別々のファイルから読み込まれます
- 各ファイルの原子インデックスは独立しています（0始まり）
- 異なる分子実体間でのオフセット調整は不要です
- リガンドとプローブはRDKit（SDF/PDB形式）で読み込まれます
- タンパク質はBioPython（PDB形式）で読み込まれます
- 重ね合わせは共結晶化リガンドとプローブの間で行われます
- その変換がタンパク質構造全体に適用されます

処理の流れ
----------
1. 共結晶化リガンドをRDKitで読み込み（SDF形式）
2. 各マッチングについて：
   a. タンパク質構造をBioPythonで読み込み
   b. プローブ分子をRDKitで読み込み（PDB形式）
   c. 原子ペアのインデックスを読み込み
   d. 対応する原子の座標を取得
   e. SuperImposerで重ね合わせパラメータを計算
   f. タンパク質全体に変換を適用
   g. 結果を保存

使用方法
--------
このスクリプトを実行する前に、atom_matching.pyを実行して
原子マッチングデータを生成してください。

    $ cd examples
    $ python atom_matching.py
    $ python superimposition.py

出力
----
各マッチングについて、重ね合わせ後のタンパク質構造が
4hw3_aligned_to_{matching_id}.pdb として保存されます。
"""
from inverse_msmd.utils.bio_utils import SuperImposer, PDB
from rdkit import Chem
import numpy as np

# ファイルパス
protein_file = "../data/sample_proteins/4hw3_A.pdb"
ligand_file = "../data/atom_matching/4hw3_A_lig.sdf"
matching_ids = [f"A08_{i}" for i in range(12)] + [f"E24_{i}" for i in range(24)]
probe_pdb = "../data/sample_probes/{probe}.pdb"

# RDKitを使用して共結晶化リガンドを読み込み
ref_lig_mol = [mol for mol in Chem.SDMolSupplier(ligand_file)][0]
ligand_coords = ref_lig_mol.GetConformer().GetPositions()

for matching in matching_ids:
    cid = matching[:3]
    
    # BioPythonを使用してタンパク質構造を読み込み（各マッチング毎に新規）
    protein = PDB.get_structure(protein_file)
    protein_coords = PDB.get_attr(protein, "coord")
    
    # RDKitを使用してプローブを読み込み
    probe_mol = Chem.rdmolfiles.MolFromPDBFile(probe_pdb.format(probe=cid))
    probe_coords = probe_mol.GetConformer().GetPositions()

    # 原子ペアのインデックスを読み込み
    # 各行: [プローブの原子インデックス, リガンドの原子インデックス]
    # インデックスは0始まりで、それぞれの構造に対して相対的
    atom_pairs = np.loadtxt(f"../data/atom_matching/atom_matching_{matching}", int)

    # マッチした原子の座標を取得
    probe_coords_target = probe_coords[atom_pairs[0]]
    ligand_coords_target = ligand_coords[atom_pairs[1]]

    # プローブと共結晶化リガンドの間で重ね合わせを実行
    si = SuperImposer()
    si.fit(ligand_coords_target, probe_coords_target)

    # タンパク質構造全体に変換を適用
    PDB.set_attr(protein, "coord", si.transform(protein_coords))
    PDB.save(protein, f"4hw3_aligned_to_{matching}.pdb")