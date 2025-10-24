#!/usr/bin/env python3
"""
スコア付きリガンド可視化のサンプルスクリプト

このスクリプトは、integrated_substructure_replacement関数を使用して
複数パターンの置換後リガンドとマッチングスコアを可視化する方法を示します。

使用方法
--------
まず、プロファイルファイルが必要です。以下のコマンドを実行してください：

    $ cd examples
    $ python visualize_scored_ligands.py

出力
----
- output/scored_visualization/pattern_N_ligand_replaced.sdf: 各パターンの置換後リガンド
- output/scored_visualization/pattern_N_protein_aligned.pdb: 各パターンの変換後タンパク質
- output/scored_visualization/results.csv: スコア付き結果CSV
- output/scored_visualization/scored_ligands.png: スコア付き可視化画像
"""

from inverse_msmd.substructure_replacement import integrated_substructure_replacement
from pathlib import Path

def main():
    """メイン関数"""
    print("=" * 70)
    print("スコア付きリガンド可視化のサンプル")
    print("=" * 70)
    print()
    
    # 出力ディレクトリ
    output_dir = "output/scored_visualization"
    
    # データファイルのパス
    ligand_file = "../data/atom_matching/4hw3_A_lig.sdf"
    protein_file = "../data/sample_proteins/4hw3_A.pdb"
    from_file = "../data/sample_probes/E23"  # 拡張子なし
    to_file = "../data/sample_probes/E24"    # 拡張子なし
    
    # プロファイルディレクトリ（スコア計算用）
    profile_dir = "../data/profiles"
    probe_id = "E24"
    
    # CSV出力パス
    csv_output = f"{output_dir}/results.csv"
    
    # 画像出力パス
    image_output = f"{output_dir}/scored_ligands.png"
    
    print(f"入力ファイル:")
    print(f"  リガンド      : {ligand_file}")
    print(f"  タンパク質    : {protein_file}")
    print(f"  置換前プローブ: {from_file}")
    print(f"  置換後プローブ: {to_file}")
    print(f"\nプロファイル:")
    print(f"  ディレクトリ  : {profile_dir}")
    print(f"  プローブID    : {probe_id}")
    print(f"\n出力:")
    print(f"  ディレクトリ  : {output_dir}")
    print(f"  CSV           : {csv_output}")
    print(f"  画像          : {image_output}")
    print()
    
    # 統合置換処理を実行
    print("処理を開始します...")
    print("-" * 70)
    
    try:
        results = integrated_substructure_replacement(
            ligand_file=ligand_file,
            protein_file=protein_file,
            from_file=from_file,
            to_file=to_file,
            output_dir=output_dir,
            match_index=0,  # 最初のマッチを使用
            profile_dir=profile_dir,
            probe_id=probe_id,
            csv_output=csv_output,
            image_output=image_output,
            deduplicate_by_smiles=True  # SMILES重複除去を有効化
        )
        
        print("-" * 70)
        print()
        print("=" * 70)
        print("✓ 処理完了")
        print("=" * 70)
        print()
        print(f"生成されたパターン数: {len(results)}")
        print()
        
        if results:
            print("スコア順の結果:")
            for i, result in enumerate(results):
                score = result.get('score', 'N/A')
                print(f"  {i+1}. Pattern {result['pattern_index']}: スコア = {score:.2f}" if score != 'N/A' else f"  {i+1}. Pattern {result['pattern_index']}")
        
        print()
        print("出力ファイル:")
        print(f"  📊 CSV : {csv_output}")
        print(f"  🖼️  画像: {image_output}")
        print()
        print("次のステップ:")
        print(f"  画像を開いて結果を確認してください：")
        print(f"    $ code {image_output}")
        print()
        
    except Exception as e:
        print()
        print("=" * 70)
        print("✗ エラーが発生しました")
        print("=" * 70)
        print(f"エラー: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())