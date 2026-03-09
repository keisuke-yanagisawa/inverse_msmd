#!/usr/bin/env python3
"""
プロファイル生成CLIスクリプト

MSMDトラジェクトリ群から残基相互作用プロファイルを一括生成する。

使用例:
    # 基本的な使用法（単一トラジェクトリ）
    python scripts/generate_profiles.py \
        --trajectories system0/traj.xtc \
        --topologies system0/topology.pdb \
        --ref-probe probe/A17/A17.pdb \
        --probe-resname A17 \
        --probe-id A17 \
        --output profiles/A17

    # 複数トラジェクトリ（system0〜system19）
    python scripts/generate_profiles.py \
        --trajectories system*/traj.xtc \
        --topologies system*/topology.pdb \
        --ref-probe probe/A17/A17.pdb \
        --probe-resname A17 \
        --probe-id A17 \
        --output profiles/A17 \
        --n-jobs 4

    # 特定のアミノ酸のみ
    python scripts/generate_profiles.py \
        --trajectories system*/traj.xtc \
        --topologies system*/topology.pdb \
        --ref-probe probe/A17/A17.pdb \
        --probe-resname A17 \
        --probe-id A17 \
        --output profiles/A17 \
        --amino-acids ALA LEU VAL ILE

    # Step別の実行も可能（後述の個別スクリプト参照）
"""

import argparse
import logging
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="MSMDトラジェクトリから残基相互作用プロファイルを生成",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--trajectories", nargs="+", required=True,
        help="トラジェクトリファイル (.xtc) のパス（複数指定可）",
    )
    parser.add_argument(
        "--topologies", nargs="+", required=True,
        help="トポロジーファイル (.pdb) のパス（--trajectoriesと同数必要）",
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
        "--probe-id", required=True,
        help="出力ファイル名用のプローブID (例: A17)",
    )
    parser.add_argument(
        "--output", required=True,
        help="出力ディレクトリ",
    )
    parser.add_argument(
        "--amino-acids", nargs="+", default=None,
        help="対象アミノ酸（デフォルト: 全20種）",
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
        "--no-compress", action="store_true",
        help="gzip圧縮を行わない",
    )
    parser.add_argument(
        "--n-jobs", type=int, default=1,
        help="並列ジョブ数 (デフォルト: 1)",
    )
    parser.add_argument(
        "--eps", type=float, default=0.1,
        help="バルク正規化時のε (デフォルト: 0.1)",
    )
    parser.add_argument(
        "--stop-on-error", action="store_true",
        help="エラー発生時に処理を中断（デフォルトは継続）",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="詳細ログを表示",
    )

    args = parser.parse_args()

    # ログ設定
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    # バリデーション
    if len(args.trajectories) != len(args.topologies):
        print(
            f"エラー: トラジェクトリ数({len(args.trajectories)})と"
            f"トポロジー数({len(args.topologies)})が一致しません"
        )
        return 1

    pairs = list(zip(args.trajectories, args.topologies))

    # 入力ファイルの存在確認
    missing = []
    for traj, topo in pairs:
        if not Path(traj).exists():
            missing.append(f"トラジェクトリ: {traj}")
        if not Path(topo).exists():
            missing.append(f"トポロジー: {topo}")
    if not Path(args.ref_probe).exists():
        missing.append(f"参照プローブ: {args.ref_probe}")
    if missing:
        print("エラー: 以下のファイルが見つかりません:")
        for m in missing:
            print(f"  - {m}")
        return 1

    print("=" * 70)
    print("プロファイル生成")
    print("=" * 70)
    print(f"トラジェクトリ数 : {len(pairs)}")
    print(f"参照プローブ     : {args.ref_probe}")
    print(f"プローブ残基名   : {args.probe_resname}")
    print(f"プローブID       : {args.probe_id}")
    print(f"出力ディレクトリ : {args.output}")
    print(f"グリッドサイズ   : {args.grid_size}")
    print(f"グリッド間隔     : {args.grid_pitch} Å")
    print(f"並列ジョブ数     : {args.n_jobs}")
    if args.amino_acids:
        print(f"対象アミノ酸     : {', '.join(args.amino_acids)}")
    print()

    from inverse_msmd.profile_generator import generate_profiles

    try:
        output_files = generate_profiles(
            trajectory_topology_pairs=pairs,
            ref_probe_path=args.ref_probe,
            probe_resname=args.probe_resname,
            probe_id=args.probe_id,
            output_dir=args.output,
            amino_acids=args.amino_acids,
            grid_size=args.grid_size,
            grid_pitch=args.grid_pitch,
            use_latter_half=not args.no_latter_half,
            compress=not args.no_compress,
            n_jobs=args.n_jobs,
            eps=args.eps,
            continue_on_error=not args.stop_on_error,
        )

        print("\n" + "=" * 70)
        print(f"完了: {len(output_files)} プロファイルを生成しました")
        print("=" * 70)
        for f in output_files:
            print(f"  {f}")
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
