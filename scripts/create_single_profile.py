#!/usr/bin/env python3
"""
単一プロファイル生成CLIスクリプト

単一トラジェクトリから1アミノ酸の残基相互作用プロファイルを生成する。
HPC環境でのジョブスケジューラ連携を想定した個別ステップCLI。

使用例:
    python scripts/create_single_profile.py \
        --trajectory system0/simulation/protein_probe.xtc \
        --topology system0/prep/protein_probe.pdb \
        --ref-probe probe/A17/A17.pdb \
        --probe-resname A17 \
        --amino-acid ALA \
        --output-dx single_profiles/sys0_A17_ALA_environment.dx \
        --output-bulk-cnt single_profiles/sys0_A17_ALA_bulk_cnt.txt
"""

import argparse
import logging
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="単一トラジェクトリから残基相互作用プロファイルを生成",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--trajectory", required=True,
        help="トラジェクトリファイル (.xtc) のパス",
    )
    parser.add_argument(
        "--topology", required=True,
        help="トポロジーファイル (.pdb) のパス",
    )
    parser.add_argument(
        "--ref-probe", required=True,
        help="参照プローブPDBファイルのパス",
    )
    parser.add_argument(
        "--probe-resname", required=True,
        help="トラジェクトリ中のプローブ残基名 (例: A17, E14)",
    )
    parser.add_argument(
        "--amino-acid", required=True,
        help="対象アミノ酸の3文字コード (例: ALA, LEU)",
    )
    parser.add_argument(
        "--output-dx", required=True,
        help="出力DXファイルのパス",
    )
    parser.add_argument(
        "--output-bulk-cnt", required=True,
        help="バルクカウント出力ファイルのパス",
    )
    parser.add_argument(
        "--grid-size", type=int, default=100,
        help="グリッドの座標数 (デフォルト: 100)",
    )
    parser.add_argument(
        "--grid-pitch", type=float, default=1.0,
        help="グリッドの間隔 Å (デフォルト: 1.0)",
    )
    parser.add_argument(
        "--no-latter-half", action="store_true",
        help="トラジェクトリ全体を使用（デフォルトは後半のみ）",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="詳細ログを表示",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    # 入力ファイルの存在確認
    missing = []
    if not Path(args.trajectory).exists():
        missing.append(f"トラジェクトリ: {args.trajectory}")
    if not Path(args.topology).exists():
        missing.append(f"トポロジー: {args.topology}")
    if not Path(args.ref_probe).exists():
        missing.append(f"参照プローブ: {args.ref_probe}")
    if missing:
        print("エラー: 以下のファイルが見つかりません:")
        for m in missing:
            print(f"  - {m}")
        return 1

    print("=" * 70)
    print("単一プロファイル生成")
    print("=" * 70)
    print(f"トラジェクトリ   : {args.trajectory}")
    print(f"トポロジー       : {args.topology}")
    print(f"参照プローブ     : {args.ref_probe}")
    print(f"プローブ残基名   : {args.probe_resname}")
    print(f"対象アミノ酸     : {args.amino_acid}")
    print(f"出力DX           : {args.output_dx}")
    print(f"出力バルクカウント: {args.output_bulk_cnt}")
    print(f"グリッドサイズ   : {args.grid_size}")
    print(f"グリッド間隔     : {args.grid_pitch} Å")
    print()

    from inverse_msmd.profile_generator import create_single_profile

    try:
        output_dx, output_bulk_cnt = create_single_profile(
            trajectory_path=args.trajectory,
            topology_path=args.topology,
            ref_probe_path=args.ref_probe,
            probe_resname=args.probe_resname,
            amino_acid=args.amino_acid,
            output_dx=args.output_dx,
            output_bulk_cnt=args.output_bulk_cnt,
            grid_size=args.grid_size,
            grid_pitch=args.grid_pitch,
            use_latter_half=not args.no_latter_half,
        )

        if output_dx is None:
            print(f"\n警告: {args.amino_acid} はこのトラジェクトリに存在しません。スキップしました。")
            return 0

        print("\n" + "=" * 70)
        print("完了")
        print("=" * 70)
        print(f"  DX           : {output_dx}")
        print(f"  バルクカウント: {output_bulk_cnt}")
        return 0

    except ImportError as e:
        print(f"\nエラー: {e}")
        return 1
    except Exception as e:
        print(f"\nエラー: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
