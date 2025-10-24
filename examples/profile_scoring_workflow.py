"""
プロファイルスコア計算統合ワークフロー - サンプルスクリプト

このスクリプトは、プロファイルスコア計算機能を統合した
部分構造置換ワークフローの使用例を示します。

使用例:
    python examples/profile_scoring_workflow.py
"""

from pathlib import Path
from inverse_msmd.substructure_replacement import integrated_substructure_replacement


def main():
    """メイン処理"""
    
    print("=" * 70)
    print("プロファイルスコア計算統合ワークフロー - サンプル")
    print("=" * 70)
    print()
    
    # 入力ファイルのパス
    ligand_file = "data/atom_matching/4hw3_A_lig.sdf"
    protein_file = "data/sample_proteins/4hw3_A.pdb"
    from_file = "data/sample_probes/E23"
    to_file = "data/sample_probes/E24"
    profile_dir = "data/profiles/"
    output_dir = "output/profile_scoring_example/"
    
    # ファイルの存在確認
    required_files = [
        ligand_file,
        protein_file,
        f"{from_file}.pdb",
        f"{from_file}.smi",
        f"{to_file}.pdb",
        f"{to_file}.smi",
    ]
    
    missing_files = [f for f in required_files if not Path(f).exists()]
    if missing_files:
        print("エラー: 以下のファイルが見つかりません:")
        for f in missing_files:
            print(f"  - {f}")
        return
    
    if not Path(profile_dir).exists():
        print(f"エラー: プロファイルディレクトリが見つかりません: {profile_dir}")
        return
    
    print("入力ファイル:")
    print(f"  リガンド      : {ligand_file}")
    print(f"  タンパク質    : {protein_file}")
    print(f"  置換前部分構造: {from_file}")
    print(f"  置換後部分構造: {to_file}")
    print(f"  プロファイル  : {profile_dir}")
    print(f"  出力先        : {output_dir}")
    print()
    
    # ===================================================================
    # 例1: スコア計算なしの基本的な使用（後方互換性）
    # ===================================================================
    print("-" * 70)
    print("例1: スコア計算なし（既存の動作）")
    print("-" * 70)
    
    output_dir_no_score = f"{output_dir}/no_score/"
    csv_output_no_score = f"{output_dir}/no_score/results.csv"
    
    results_no_score = integrated_substructure_replacement(
        ligand_file=ligand_file,
        protein_file=protein_file,
        from_file=from_file,
        to_file=to_file,
        output_dir=output_dir_no_score,
        match_index=0,
        profile_dir=None,  # スコア計算をスキップ
        csv_output=csv_output_no_score
    )
    
    print(f"\n結果: {len(results_no_score)} パターン生成")
    for i, result in enumerate(results_no_score):
        print(f"  パターン {i}:")
        print(f"    SMILES    : {result['ligand_smiles']}")
        print(f"    リガンド  : {result['ligand_file']}")
        print(f"    タンパク質: {result['protein_file']}")
    print(f"\nCSV出力: {csv_output_no_score}")
    
    print()
    
    # ===================================================================
    # 例2: プロファイルスコア計算あり
    # ===================================================================
    print("-" * 70)
    print("例2: プロファイルスコア計算")
    print("-" * 70)
    
    output_dir_score = f"{output_dir}/with_score/"
    csv_output_score = f"{output_dir}/with_score/results.csv"
    
    results_score = integrated_substructure_replacement(
        ligand_file=ligand_file,
        protein_file=protein_file,
        from_file=from_file,
        to_file=to_file,
        output_dir=output_dir_score,
        match_index=0,
        profile_dir=profile_dir,
        probe_id="E24",
        csv_output=csv_output_score
    )
    
    print(f"\n結果: {len(results_score)} パターン生成（スコア降順ソート済み）")
    for i, result in enumerate(results_score):
        print(f"  パターン {i}:")
        print(f"    スコア    : {result['score']:.2f}")
        print(f"    SMILES    : {result['ligand_smiles']}")
        print(f"    リガンド  : {result['ligand_file']}")
        print(f"    タンパク質: {result['protein_file']}")
    print(f"\nCSV出力: {csv_output_score}")
    
    print()
    print("=" * 70)
    print("完了")
    print("=" * 70)
    print()
    print(f"出力ディレクトリ: {output_dir}")
    print()
    print("次のステップ:")
    print("  1. CSVファイルで結果を確認")
    print("  2. 生成されたSDFファイルを分子可視化ソフト（PyMOL、Chimera等）で確認")
    print("  3. スコアが高いパターンを優先的に検討")
    print("  4. 必要に応じてgamma値を調整して再実行")


if __name__ == "__main__":
    main()