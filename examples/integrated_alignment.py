"""
統合されたアライメントAPIの使用例。

このスクリプトは、inverse_msmdパッケージの統合アライメント機能を使用して、
プローブ分子と共結晶化リガンドのマッチングに基づいて
タンパク質構造を重ね合わせます。

従来のatom_matching.pyとsuperimposition.pyを統合した
シンプルなAPIを使用します。
"""

from inverse_msmd import align_structures

# ファイルパス
protein_file = "../data/sample_proteins/4hw3_A.pdb"
ligand_file = "../data/atom_matching/4hw3_A_lig_with_subst_label.sdf"

# プローブファイル（辞書形式で指定）
probe_files = {
    "A08": "../data/sample_probes/A08.pdb",
    "E24": "../data/sample_probes/E24.pdb",
}

# 統合APIを使用してアライメントを実行
# - atom matching（MCS検索、原子ペア特定）
# - structure superimposition（重ね合わせ）
# を一括で実行
results = align_structures(
    protein_file=protein_file,
    ligand_file=ligand_file,
    probe_files=probe_files,
    output_dir="./aligned_structures",
    iso_value=1  # ISO=1ラベルのある原子のみを使用
)

# 結果の表示
print(f"\n総計: {len(results)} 個のアライメント結果を生成")
print("\n結果サマリー:")
for result in results:
    n_atoms = result.atom_pairs.shape[1]
    print(f"  {result.probe_id}_{result.match_id}: "
          f"マッチした原子数 = {n_atoms}")

print("\n出力ファイル: ./aligned_structures/aligned_to_*.pdb")