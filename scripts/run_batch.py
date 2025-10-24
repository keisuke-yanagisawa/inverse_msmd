#!/usr/bin/env python
"""
バッチ処理CLIスクリプト

複数の置換パターン（from_probe, to_probe, match_index）に対して、
matching scoreを一括計算するコマンドラインツールです。

使用例:
    # 基本的な使用方法
    python scripts/run_batch.py \\
        --batch-csv experiments/batch_config.csv \\
        --ligand data/atom_matching/4hw3_A_lig.sdf \\
        --protein data/sample_proteins/4hw3_A.pdb \\
        --probe-dir data/sample_probes \\
        --profile-dir data/profiles \\
        --output output/batch_results
    
    # 並列処理を使用
    python scripts/run_batch.py \\
        --batch-csv experiments/batch_config.csv \\
        --ligand data/atom_matching/4hw3_A_lig.sdf \\
        --protein data/sample_proteins/4hw3_A.pdb \\
        --probe-dir data/sample_probes \\
        --profile-dir data/profiles \\
        --output output/batch_results \\
        --parallel \\
        --max-workers 4
    
    # エラー時に処理を中断
    python scripts/run_batch.py \\
        --batch-csv experiments/batch_config.csv \\
        ... \\
        --no-continue-on-error
"""

import argparse
import sys
from pathlib import Path

from inverse_msmd.batch_processing import run_batch_processing


def parse_args():
    """コマンドライン引数をパースします"""
    parser = argparse.ArgumentParser(
        description="複数の置換パターンに対してmatching scoreを一括計算します",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用例:
  %(prog)s --batch-csv batch_config.csv \\
      --ligand ligand.sdf --protein protein.pdb \\
      --probe-dir data/probes --profile-dir data/profiles \\
      --output output/results

  %(prog)s --batch-csv batch_config.csv \\
      --ligand ligand.sdf --protein protein.pdb \\
      --probe-dir data/probes --profile-dir data/profiles \\
      --output output/results --parallel --max-workers 8

詳細は docs/batch_processing/README.md を参照してください。
        """
    )
    
    # 必須引数
    required = parser.add_argument_group('必須引数')
    required.add_argument(
        '--batch-csv',
        required=True,
        help='バッチ設定CSVファイルのパス'
    )
    required.add_argument(
        '--ligand',
        required=True,
        help='リガンドSDFファイルのパス'
    )
    required.add_argument(
        '--protein',
        required=True,
        help='タンパク質PDBファイルのパス'
    )
    required.add_argument(
        '--probe-dir',
        required=True,
        help='プローブファイルのベースディレクトリ'
    )
    required.add_argument(
        '--profile-dir',
        required=True,
        help='プロファイルファイルのベースディレクトリ'
    )
    required.add_argument(
        '--output',
        required=True,
        help='出力ベースディレクトリ'
    )
    
    # オプション引数
    optional = parser.add_argument_group('オプション引数')
    optional.add_argument(
        '--parallel',
        action='store_true',
        help='並列処理を有効にする'
    )
    optional.add_argument(
        '--max-workers',
        type=int,
        default=4,
        help='並列処理時の最大ワーカー数（デフォルト: 4）'
    )
    optional.add_argument(
        '--no-continue-on-error',
        action='store_true',
        help='エラー発生時に処理を中断する（デフォルトは継続）'
    )
    optional.add_argument(
        '--log-file',
        help='ログファイルのパス（指定しない場合は<output>/batch_execution.log）'
    )
    optional.add_argument(
        '--visualize',
        action='store_true',
        help='バッチ処理完了後に結果を可視化する（成功したジョブのリガンドとスコアをグリッド表示）'
    )
    optional.add_argument(
        '--visualization-output',
        help='可視化画像の出力パス（デフォルト: <output>/batch_visualization.png）'
    )
    
    return parser.parse_args()


def validate_paths(args):
    """入力ファイルとディレクトリの存在を検証します"""
    errors = []
    
    # 必須ファイルの存在確認
    required_files = {
        'バッチCSV': args.batch_csv,
        'リガンドファイル': args.ligand,
        'タンパク質ファイル': args.protein,
    }
    
    for name, filepath in required_files.items():
        if not Path(filepath).exists():
            errors.append(f"{name}が見つかりません: {filepath}")
    
    # ディレクトリの存在確認
    required_dirs = {
        'プローブディレクトリ': args.probe_dir,
        'プロファイルディレクトリ': args.profile_dir,
    }
    
    for name, dirpath in required_dirs.items():
        if not Path(dirpath).exists():
            errors.append(f"{name}が見つかりません: {dirpath}")
    
    return errors


def main():
    """メイン処理"""
    args = parse_args()
    
    # パスの検証
    print("入力ファイルを検証中...")
    errors = validate_paths(args)
    if errors:
        print("\nエラー: 以下のファイル/ディレクトリが見つかりません:")
        for error in errors:
            print(f"  - {error}")
        return 1
    
    print("検証完了\n")
    
    # パラメータの表示
    print("=" * 70)
    print("バッチ処理の設定")
    print("=" * 70)
    print(f"バッチCSV: {args.batch_csv}")
    print(f"リガンドファイル: {args.ligand}")
    print(f"タンパク質ファイル: {args.protein}")
    print(f"プローブディレクトリ: {args.probe_dir}")
    print(f"プロファイルディレクトリ: {args.profile_dir}")
    print(f"出力ディレクトリ: {args.output}")
    print(f"並列処理: {'有効' if args.parallel else '無効'}")
    if args.parallel:
        print(f"最大ワーカー数: {args.max_workers}")
    print(f"エラー時の継続: {'無効' if args.no_continue_on_error else '有効'}")
    print()
    
    # バッチ処理を実行
    try:
        result = run_batch_processing(
            batch_csv=args.batch_csv,
            ligand_file=args.ligand,
            protein_file=args.protein,
            probe_base_dir=args.probe_dir,
            profile_base_dir=args.profile_dir,
            output_base_dir=args.output,
            parallel=args.parallel,
            max_workers=args.max_workers,
            continue_on_error=not args.no_continue_on_error,
            log_file=args.log_file
        )
        
        # 結果のサマリーを表示
        print("\n" + "=" * 70)
        print("バッチ処理完了")
        print("=" * 70)
        print(f"総ジョブ数: {result.total_jobs}")
        print(f"成功: {result.num_success}")
        print(f"失敗: {result.num_failed}")
        print(f"スキップ: {result.num_skipped}")
        print(f"実行時間: {result.total_execution_time:.2f}秒")
        
        # 成功率
        if result.total_jobs > 0:
            success_rate = (result.num_success / result.total_jobs) * 100
            print(f"成功率: {success_rate:.1f}%")
        
        # 出力ファイルの案内
        print("\n" + "=" * 70)
        print("出力ファイル")
        print("=" * 70)
        output_path = Path(args.output)
        print(f"サマリーCSV: {output_path / 'batch_summary.csv'}")
        print(f"サマリーJSON: {output_path / 'batch_summary.json'}")
        log_file = args.log_file or str(output_path / 'batch_execution.log')
        print(f"実行ログ: {log_file}")
        print()
        
        # 可視化処理（オプション）
        if args.visualize and result.num_success > 0:
            print("=" * 70)
            print("結果を可視化中...")
            print("=" * 70)
            try:
                from inverse_msmd.utils.visualization_utils import draw_scored_molecules_grid
                from rdkit import Chem
                import pandas as pd
                
                # サマリーCSVから成功したジョブの情報を読み込み
                summary_df = pd.read_csv(output_path / 'batch_summary.csv')
                success_df = summary_df[summary_df['status'] == 'success'].copy()
                
                if len(success_df) > 0:
                    # スコアでソート
                    success_df = success_df.sort_values('best_score', ascending=False)
                    
                    molecules = []
                    scores = []
                    titles = []
                    
                    for _, row in success_df.iterrows():
                        job_id = row['job_id']
                        job_dir = output_path / job_id
                        
                        # best_pattern_indexのSDFファイルを読み込み
                        pattern_index = int(row['best_pattern_index']) if pd.notna(row['best_pattern_index']) else 0
                        sdf_file = job_dir / f"pattern_{pattern_index}_ligand_replaced.sdf"
                        
                        if sdf_file.exists():
                            suppl = Chem.SDMolSupplier(str(sdf_file))
                            mol = next(suppl)
                            if mol is not None:
                                molecules.append(mol)
                                scores.append(float(row['best_score']))
                                titles.append(f"{job_id}")
                    
                    if molecules:
                        # 可視化出力パス
                        viz_output = args.visualization_output or str(output_path / 'batch_visualization.png')
                        
                        # グリッド描画
                        draw_scored_molecules_grid(
                            molecules,
                            scores,
                            viz_output,
                            titles=titles,
                            max_cols=4,
                            align_molecules=True
                        )
                        
                        print(f"可視化画像を保存しました: {viz_output}")
                        print()
                
            except Exception as e:
                print(f"警告: 可視化中にエラーが発生しました: {e}")
                print("バッチ処理は正常に完了しています。")
                import traceback
                traceback.print_exc()
        
        # 失敗したジョブがある場合は警告
        if result.num_failed > 0:
            print("警告: 一部のジョブが失敗しました。")
            print("詳細はサマリーCSVまたはログファイルを確認してください。")
            return 1
        
        return 0
        
    except KeyboardInterrupt:
        print("\n\n処理が中断されました")
        return 130
    except Exception as e:
        print(f"\nエラー: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())