#!/usr/bin/env python
"""
統合部分構造置換 - コマンドラインインターフェース

リガンド中の部分構造を別の部分構造で置換し、
タンパク質構造も適切に座標変換します。

使用例:
    python scripts/integrated_replacement.py \
        --ligand data/atom_matching/4hw3_A_lig.sdf \
        --protein data/sample_proteins/4hw3_A.pdb \
        --from-file data/sample_probes/E23 \
        --to-file data/sample_probes/E24 \
        --output output/integrated/
"""

import argparse
import sys
from pathlib import Path


def main():
    """メイン関数"""
    parser = argparse.ArgumentParser(
        description="統合部分構造置換: リガンドの部分構造を置換し、タンパク質も座標変換します",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用例:
  基本的な使用:
    %(prog)s --ligand ligand.sdf --protein protein.pdb \\
             --from-file probes/E23 --to-file probes/E24 \\
             --output output/

  マッチインデックスを指定:
    %(prog)s --ligand ligand.sdf --protein protein.pdb \\
             --from-file probes/E23 --to-file probes/E24 \\
             --output output/ --match-index 1

注意:
  --from-fileと--to-fileには拡張子なしのベースパスを指定してください。
  .pdbと.smiファイルが自動的に読み込まれます。
        """
    )
    
    # 必須引数
    parser.add_argument(
        "--ligand",
        required=True,
        help="リガンドSDFファイルのパス"
    )
    parser.add_argument(
        "--protein",
        required=True,
        help="タンパク質PDBファイルのパス"
    )
    parser.add_argument(
        "--from-file",
        required=True,
        help="置換前の部分構造のベースパス（拡張子なし）"
    )
    parser.add_argument(
        "--to-file",
        required=True,
        help="置換後の部分構造のベースパス（拡張子なし）"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="出力ディレクトリのパス"
    )
    
    # オプション引数
    parser.add_argument(
        "--match-index",
        type=int,
        default=None,
        help="部分構造マッチのインデックス指定（0始まり）。"
             "指定しない場合、複数マッチ時は可視化画像を出力します"
    )
    parser.add_argument(
        "--profile-dir",
        default=None,
        help="プロファイルディレクトリのパス。"
             "指定した場合、各パターンのプロファイルスコアを計算します"
    )
    parser.add_argument(
        "--probe-id",
        default=None,
        help="プローブID（例: E24）。"
             "profile-dirが指定されている場合は必須です"
    )
    parser.add_argument(
        "--csv-output",
        default=None,
        help="CSV出力ファイルのパス。"
             "指定した場合、結果をCSVファイルとして保存します"
    )
    parser.add_argument(
        "--image-output",
        default=None,
        help="画像出力ファイルのパス（PNG形式）。"
             "指定した場合、全パターンの置換後リガンドとスコアを可視化した画像を生成します。"
             "profile-dirが指定されている場合のみ有効です"
    )
    parser.add_argument(
        "--deduplicate-by-smiles",
        action="store_true",
        help="SMILES表記が同じリガンドの重複を除去します。"
             "同じSMILES構造を持つパターンの中で最高スコアのものだけを保持します。"
             "profile-dirが指定されている場合のみ有効です"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="詳細な出力を表示"
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s 2.0.0"
    )
    
    args = parser.parse_args()
    
    # ファイルの存在確認
    ligand_path = Path(args.ligand)
    protein_path = Path(args.protein)
    from_pdb = Path(f"{args.from_file}.pdb")
    from_smi = Path(f"{args.from_file}.smi")
    to_pdb = Path(f"{args.to_file}.pdb")
    to_smi = Path(f"{args.to_file}.smi")
    
    errors = []
    if not ligand_path.exists():
        errors.append(f"リガンドファイルが見つかりません: {ligand_path}")
    if not protein_path.exists():
        errors.append(f"タンパク質ファイルが見つかりません: {protein_path}")
    if not from_pdb.exists():
        errors.append(f"置換前部分構造PDBが見つかりません: {from_pdb}")
    if not from_smi.exists():
        errors.append(f"置換前部分構造SMIが見つかりません: {from_smi}")
    if not to_pdb.exists():
        errors.append(f"置換後部分構造PDBが見つかりません: {to_pdb}")
    if not to_smi.exists():
        errors.append(f"置換後部分構造SMIが見つかりません: {to_smi}")
    
    if errors:
        print("エラー: 必要なファイルが見つかりません\n", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        sys.exit(1)
    
    # 統合ワークフローをインポートして実行
    try:
        from inverse_msmd.substructure_replacement import integrated_substructure_replacement
    except ImportError as e:
        print(f"エラー: モジュールのインポートに失敗しました: {e}", file=sys.stderr)
        print("inverse_msmdパッケージがインストールされているか確認してください", file=sys.stderr)
        sys.exit(1)
    
    # プロファイル計算のパラメータ検証
    if args.profile_dir is not None and args.probe_id is None:
        print("エラー: --profile-dirを指定する場合、--probe-idも必須です", file=sys.stderr)
        sys.exit(1)
    
    if args.verbose:
        print("=" * 70)
        print("統合部分構造置換を開始します")
        print("=" * 70)
        print(f"リガンド      : {args.ligand}")
        print(f"タンパク質    : {args.protein}")
        print(f"置換前部分構造: {args.from_file}")
        print(f"置換後部分構造: {args.to_file}")
        print(f"出力先        : {args.output}")
        print(f"マッチ指定    : {args.match_index if args.match_index is not None else '自動'}")
        if args.profile_dir:
            print(f"プロファイル  : {args.profile_dir}")
            print(f"プローブID    : {args.probe_id}")
        if args.csv_output:
            print(f"CSV出力       : {args.csv_output}")
        if args.image_output:
            print(f"画像出力      : {args.image_output}")
        if args.deduplicate_by_smiles:
            print(f"SMILES重複除去: 有効")
        print("=" * 70)
        print()
    
    try:
        # 統合ワークフローを実行
        results = integrated_substructure_replacement(
            ligand_file=str(ligand_path),
            protein_file=str(protein_path),
            from_file=args.from_file,
            to_file=args.to_file,
            output_dir=args.output,
            match_index=args.match_index,
            profile_dir=args.profile_dir,
            probe_id=args.probe_id,
            csv_output=args.csv_output,
            image_output=args.image_output,
            deduplicate_by_smiles=args.deduplicate_by_smiles
        )
        
        # 結果のサマリーを表示
        print()
        print("=" * 70)
        print(f"✓ 処理完了: {len(results)} パターンの結果を生成しました")
        print("=" * 70)
        
        if args.verbose:
            for i, result in enumerate(results):
                print(f"\nパターン {i}:")
                if 'score' in result:
                    print(f"  スコア    : {result['score']:.2f}")
                print(f"  リガンド  : {result['ligand_file']}")
                print(f"  タンパク質: {result['protein_file']}")
        else:
            print(f"\n出力ディレクトリ: {args.output}")
            print(f"  pattern_N_ligand_replaced.sdf")
            print(f"  pattern_N_protein_aligned.pdb")
            if args.profile_dir and results:
                print(f"\n✓ 結果はプロファイルスコアで降順ソート済みです")
                if results:
                    best_score = results[0].get('score')
                    if best_score is not None:
                        print(f"  最高スコア: {best_score:.2f} (パターン {results[0]['pattern_index']})")
        
        print()
        
    except ValueError as e:
        print(f"\nエラー: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"\n予期しないエラーが発生しました: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()