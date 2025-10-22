"""
原子マッチングスクリプト

このスクリプトは、RDKitを使用してプローブ分子と参照リガンドの間で
最大共通部分構造（MCS）を検出する方法を示します。
すべての可能な部分構造マッチを列挙し、保存します。

アイソトープラベルについて
------------------------
参照リガンドにアイソトープラベル（ISO=1）がある場合、
ラベル付けされた原子のみが考慮されます。これにより、
リガンドの特定の部分のみをマッチング対象にできます。

処理の流れ
----------
1. プローブ分子（PDB形式）を読み込み
2. 参照リガンド（SDF形式）を読み込み
3. ISO=1ラベルを持つ原子を抽出
4. ラベル付き原子から部分構造を作成
5. プローブと部分構造の間でMCS検索
6. すべての可能なマッチングを列挙（重複除外）
7. 各マッチングの原子ペアを保存

出力形式
--------
各マッチングは以下の形式で保存されます：
    ../data/atom_matching/atom_matching_{probe_id}_{match_id}

ファイルは2行からなり：
    1行目: プローブ側の原子インデックス（0始まり、スペース区切り）
    2行目: リファレンス側の原子インデックス（0始まり、スペース区切り）

例：
    0 1 2 3 4 5
    10 11 12 13 14 15

使用方法
--------
    $ cd examples
    $ python atom_matching.py

注意事項
--------
- サンプルデータ（A08とE24プローブ）のみを処理します
- RDKitがインストールされている必要があります
- アイソトープラベル付きリガンドファイルが必要です
"""
from rdkit import Chem
from rdkit.Chem import rdFMCS
import numpy as np

def get_iso_labeled_atoms(mol, iso_value=1):
    """
    指定されたアイソトープラベル値を持つ原子のインデックスを取得します。
    
    Parameters
    ----------
    mol : rdkit.Chem.Mol
        分子オブジェクト
    iso_value : int, optional
        抽出するアイソトープラベル値（デフォルト: 1）
    
    Returns
    -------
    list of int
        アイソトープラベルを持つ原子のインデックスリスト
    """
    iso_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetIsotope() == iso_value:
            iso_atoms.append(atom.GetIdx())
    return iso_atoms

def create_submol_from_atoms(mol, atom_indices):
    """
    指定された原子インデックスから部分構造を作成します。
    
    Parameters
    ----------
    mol : rdkit.Chem.Mol
        元の分子
    atom_indices : list of int
        抽出する原子のインデックスリスト
    
    Returns
    -------
    rdkit.Chem.Mol
        部分構造分子
        
    Notes
    -----
    指定された原子間の結合のみが部分構造に含まれます。
    """
    # 原子のセットを作成
    atom_set = set(atom_indices)
    
    # これらの原子間の結合を抽出
    bonds_to_include = []
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in atom_set and bond.GetEndAtomIdx() in atom_set:
            bonds_to_include.append(bond.GetIdx())
    
    # 部分構造を作成
    submol = Chem.PathToSubmol(mol, bonds_to_include)
    return submol

# ファイルパス
probe_pdb = "../data/sample_probes/{probe}.pdb"
ref_lig_sdf_labeled = "../data/atom_matching/4hw3_A_lig_with_subst_label.sdf"

# サンプルデータがあるプローブのみを処理
for probe_id in "A08 E24".split(" "):
    # プローブ分子を読み込み
    probe_mol = Chem.rdmolfiles.MolFromPDBFile(probe_pdb.format(probe=probe_id))
    # 参照リガンドを読み込み
    ref_lig_mol = [mol for mol in Chem.SDMolSupplier(ref_lig_sdf_labeled)][0]
    
    # ISO=1でラベル付けされた原子を取得
    iso1_atoms = get_iso_labeled_atoms(ref_lig_mol, iso_value=1)
    print(f"\n{probe_id}: 参照リガンドのISO 1ラベル原子: {sorted(iso1_atoms)}")
    
    # ISO=1の原子のみから部分構造を作成
    ref_submol = create_submol_from_atoms(ref_lig_mol, iso1_atoms)
    
    # プローブ分子全体とリファレンスのISO 1部分のMCSを計算
    mcs_result = rdFMCS.FindMCS([probe_mol, ref_submol])
    print(f"  MCS SMARTS: {mcs_result.smartsString}")
    print(f"  MCS サイズ: {mcs_result.numAtoms} 個の原子")
    
    mcs = Chem.MolFromSmarts(mcs_result.smartsString)
    
    # GetSubstructMatches（複数形）を使用してすべてのマッチングを取得
    # uniquify=Falseで重複を含むすべてのマッチングを取得
    probe_matches = probe_mol.GetSubstructMatches(mcs, uniquify=False)
    ref_submol_matches = ref_submol.GetSubstructMatches(mcs, uniquify=False)
    
    print(f"  発見: プローブマッチ {len(probe_matches)} 個、参照部分構造マッチ {len(ref_submol_matches)} 個")
    
    # 部分構造のインデックスを元の分子のインデックスにマッピング
    # ref_submolの原子インデックスをref_lig_molの原子インデックスに変換
    submol_to_mol_map = {}
    for submol_idx, atom in enumerate(ref_submol.GetAtoms()):
        # 原子マッピング番号を使用して元のインデックスを取得
        original_idx = atom.GetAtomMapNum()
        if original_idx == 0:
            # マッピング番号が設定されていない場合、ISO 1原子リストから取得
            if submol_idx < len(iso1_atoms):
                original_idx = iso1_atoms[submol_idx]
        else:
            original_idx -= 1  # 1-indexed から 0-indexed へ変換
        submol_to_mol_map[submol_idx] = original_idx
    
    # 重複除外のため、原子ペアのセットを記録
    seen_pairings = set()
    
    # すべての可能な組み合わせを保存（重複を除外）
    match_count = 0
    for i, probe_match in enumerate(probe_matches):
        for j, ref_submol_match in enumerate(ref_submol_matches):
            # 部分構造のマッチを元の分子のインデックスに変換
            ref_match = tuple(submol_to_mol_map.get(idx, iso1_atoms[idx] if idx < len(iso1_atoms) else idx)
                            for idx in ref_submol_match)
            
            # 原子ペアのセットを作成（順序に依存しない）
            pairing = frozenset(zip(probe_match, ref_match))
            
            # このペアリングが既に見たことがあるかチェック
            if pairing in seen_pairings:
                continue
            seen_pairings.add(pairing)
            
            arr = np.array([probe_match, ref_match])
            
            # 複数のマッチングがある場合は番号を付けて保存
            if len(probe_matches) > 1 or len(ref_submol_matches) > 1:
                output_file = f"../data/atom_matching/atom_matching_{probe_id}_{match_count}"
            else:
                output_file = f"../data/atom_matching/atom_matching_{probe_id}"
            
            np.savetxt(output_file, arr, fmt='%d')
            print(f"  マッチ {match_count}: {output_file} に保存")
            print(f"    プローブ原子:  {probe_match}")
            print(f"    参照原子:      {ref_match}")
            match_count += 1
    
    print(f"  保存した一意のマッチ総数: {match_count}")