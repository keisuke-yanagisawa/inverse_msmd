"""
inverse_msmd CLIエントリーポイント

コマンドラインからinverse_msmdの主要機能を実行するためのエントリーポイントを提供します。

コマンド:
    inverse-msmd-run   : 単体の部分構造置換+スコア計算
    inverse-msmd-batch : バッチ処理（複数パターン一括計算）
"""

import argparse
import sys
from pathlib import Path


def run_single_main():
    """単体部分構造置換+スコア計算のCLIエントリーポイント"""
    from inverse_msmd.substructure_replacement import integrated_substructure_replacement

    parser = argparse.ArgumentParser(
        description="リガンドの部分構造置換とプロファイルマッチングスコア計算を行います",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用例:
  # スコア計算なし
  %(prog)s --ligand ligand.sdf --protein protein.pdb \\
      --from-probe data/probes/E23 --to-probe data/probes/E24 \\
      --output output/results

  # スコア計算あり + CSV出力
  %(prog)s --ligand ligand.sdf --protein protein.pdb \\
      --from-probe data/probes/E23 --to-probe data/probes/E24 \\
      --output output/results \\
      --profile-dir data/profiles --probe-id E24 \\
      --csv-output output/results/results.csv

  # 全オプション指定
  %(prog)s --ligand ligand.sdf --protein protein.pdb \\
      --from-probe data/probes/E23 --to-probe data/probes/E24 \\
      --output output/results \\
      --profile-dir data/profiles --probe-id E24 \\
      --match-index 0 \\
      --csv-output output/results/results.csv \\
      --image-output output/results/results.png \\
      --deduplicate
        """
    )

    # 必須引数
    required = parser.add_argument_group('必須引数')
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
        '--from-probe',
        required=True,
        help='置換前の部分構造のベースパス（拡張子なし、.pdbと.smiを自動読み込み）'
    )
    required.add_argument(
        '--to-probe',
        required=True,
        help='置換後の部分構造のベースパス（拡張子なし、.pdbと.smiを自動読み込み）'
    )
    required.add_argument(
        '--output',
        required=True,
        help='出力ディレクトリのパス'
    )

    # スコア計算オプション
    scoring = parser.add_argument_group('スコア計算オプション')
    scoring.add_argument(
        '--profile-dir',
        help='プロファイルファイル（.dx.gz）のディレクトリパス（指定するとスコア計算を実行）'
    )
    scoring.add_argument(
        '--probe-id',
        help='プローブID（例: "E24"）。--profile-dir指定時は必須'
    )

    # その他のオプション
    optional = parser.add_argument_group('その他のオプション')
    optional.add_argument(
        '--match-index',
        type=int,
        default=None,
        help='部分構造マッチのインデックス（0始まり）。'
             '未指定の場合、複数マッチ時は可視化画像を生成して最初のマッチを使用'
    )
    optional.add_argument(
        '--csv-output',
        help='結果CSVファイルの出力パス'
    )
    optional.add_argument(
        '--image-output',
        help='結果可視化画像（PNG）の出力パス（スコア計算時のみ有効）'
    )
    optional.add_argument(
        '--deduplicate',
        action='store_true',
        help='SMILES表記が同じリガンドの重複を除去する（スコア計算時のみ有効）'
    )
    optional.add_argument(
        '--skip-steric-clash-check',
        action='store_true',
        help='立体障害チェックをスキップする'
    )

    args = parser.parse_args()

    # パラメータバリデーション
    if args.profile_dir and not args.probe_id:
        print("エラー: --profile-dir を指定した場合、--probe-id も必須です")
        return 1

    # パスの検証
    print("入力ファイルを検証中...")
    errors = []
    required_files = {
        'リガンドファイル': args.ligand,
        'タンパク質ファイル': args.protein,
        '置換前部分構造PDB': f"{args.from_probe}.pdb",
        '置換前部分構造SMI': f"{args.from_probe}.smi",
        '置換後部分構造PDB': f"{args.to_probe}.pdb",
        '置換後部分構造SMI': f"{args.to_probe}.smi",
    }
    for name, filepath in required_files.items():
        if not Path(filepath).exists():
            errors.append(f"{name}が見つかりません: {filepath}")
    if args.profile_dir and not Path(args.profile_dir).exists():
        errors.append(f"プロファイルディレクトリが見つかりません: {args.profile_dir}")

    if errors:
        print("\nエラー: 以下のファイル/ディレクトリが見つかりません:")
        for error in errors:
            print(f"  - {error}")
        return 1

    print("検証完了\n")

    # パラメータの表示
    print("=" * 70)
    print("実行設定")
    print("=" * 70)
    print(f"リガンドファイル    : {args.ligand}")
    print(f"タンパク質ファイル  : {args.protein}")
    print(f"置換前部分構造      : {args.from_probe}")
    print(f"置換後部分構造      : {args.to_probe}")
    print(f"出力ディレクトリ    : {args.output}")
    if args.match_index is not None:
        print(f"マッチインデックス  : {args.match_index}")
    if args.profile_dir:
        print(f"プロファイルディレクトリ: {args.profile_dir}")
        print(f"プローブID          : {args.probe_id}")
    print(f"立体障害チェック    : {'スキップ' if args.skip_steric_clash_check else '有効'}")
    if args.deduplicate:
        print(f"SMILES重複除去      : 有効")
    print()

    # 統合ワークフローを実行
    try:
        results = integrated_substructure_replacement(
            ligand_file=args.ligand,
            protein_file=args.protein,
            from_file=args.from_probe,
            to_file=args.to_probe,
            output_dir=args.output,
            match_index=args.match_index,
            profile_dir=args.profile_dir,
            probe_id=args.probe_id,
            csv_output=args.csv_output,
            image_output=args.image_output,
            deduplicate_by_smiles=args.deduplicate,
            skip_steric_clash_check=args.skip_steric_clash_check
        )

        # 結果サマリー
        print("\n" + "=" * 70)
        print("結果サマリー")
        print("=" * 70)
        print(f"生成パターン数: {len(results)}")

        if results:
            has_score = 'score' in results[0]
            for result in results:
                pattern_idx = result['pattern_index']
                smiles = result.get('ligand_smiles', 'N/A')
                if has_score:
                    score = result.get('score', float('nan'))
                    print(f"  パターン {pattern_idx}: スコア={score:.2f}, SMILES={smiles}")
                else:
                    print(f"  パターン {pattern_idx}: SMILES={smiles}")

        # 出力ファイルの案内
        print("\n" + "=" * 70)
        print("出力ファイル")
        print("=" * 70)
        for result in results:
            pattern_idx = result['pattern_index']
            print(f"  パターン {pattern_idx}:")
            print(f"    リガンド  : {result['ligand_file']}")
            print(f"    タンパク質: {result['protein_file']}")
        if args.csv_output:
            print(f"  CSV: {args.csv_output}")
        if args.image_output:
            print(f"  画像: {args.image_output}")
        print()

        return 0

    except KeyboardInterrupt:
        print("\n\n処理が中断されました")
        return 130
    except Exception as e:
        print(f"\nエラー: {e}")
        import traceback
        traceback.print_exc()
        return 1


def run_batch_main():
    """バッチ処理のCLIエントリーポイント（scripts/run_batch.pyからの移行）"""
    # 既存のrun_batch.pyのmain関数を再利用
    from inverse_msmd.batch_processing import run_batch_processing

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
        help='バッチ処理完了後に結果を可視化する'
    )
    optional.add_argument(
        '--visualization-output',
        help='可視化画像の出力パス（デフォルト: <output>/batch_visualization.png）'
    )
    optional.add_argument(
        '--skip-steric-clash-check',
        action='store_true',
        help='立体障害チェックをスキップする'
    )

    args = parser.parse_args()

    # パスの検証
    print("入力ファイルを検証中...")
    errors = []
    required_files = {
        'バッチCSV': args.batch_csv,
        'リガンドファイル': args.ligand,
        'タンパク質ファイル': args.protein,
    }
    for name, filepath in required_files.items():
        if not Path(filepath).exists():
            errors.append(f"{name}が見つかりません: {filepath}")
    required_dirs = {
        'プローブディレクトリ': args.probe_dir,
        'プロファイルディレクトリ': args.profile_dir,
    }
    for name, dirpath in required_dirs.items():
        if not Path(dirpath).exists():
            errors.append(f"{name}が見つかりません: {dirpath}")

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
    print(f"立体障害チェック: {'スキップ' if args.skip_steric_clash_check else '有効'}")
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
            log_file=args.log_file,
            skip_steric_clash_check=args.skip_steric_clash_check
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

                summary_df = pd.read_csv(output_path / 'batch_summary.csv')
                success_df = summary_df[summary_df['status'] == 'success'].copy()

                if len(success_df) > 0:
                    molecules = []
                    scores = []
                    titles = []

                    for _, row in success_df.iterrows():
                        job_id = row['job_id']
                        job_dir = output_path / job_id
                        results_csv = job_dir / 'results.csv'

                        if results_csv.exists():
                            job_results_df = pd.read_csv(results_csv)
                            for _, pattern_row in job_results_df.iterrows():
                                pattern_index = int(pattern_row['pattern_index'])
                                pattern_score = float(pattern_row['score'])
                                sdf_file = job_dir / f"pattern_{pattern_index}_ligand_replaced.sdf"

                                if sdf_file.exists():
                                    suppl = Chem.SDMolSupplier(str(sdf_file))
                                    mol = next(suppl)
                                    if mol is not None:
                                        molecules.append(mol)
                                        scores.append(pattern_score)
                                        titles.append(f"{job_id}_p{pattern_index}")

                    if molecules:
                        sorted_indices = sorted(
                            range(len(scores)), key=lambda i: scores[i], reverse=True
                        )
                        molecules = [molecules[i] for i in sorted_indices]
                        scores = [scores[i] for i in sorted_indices]
                        titles = [titles[i] for i in sorted_indices]

                        viz_output = args.visualization_output or str(
                            output_path / 'batch_visualization.png'
                        )
                        draw_scored_molecules_grid(
                            molecules, scores, viz_output,
                            titles=titles, max_cols=4, align_molecules=True
                        )
                        print(f"可視化画像を保存しました: {viz_output}")

            except Exception as e:
                print(f"警告: 可視化中にエラーが発生しました: {e}")
                print("バッチ処理は正常に完了しています。")

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
