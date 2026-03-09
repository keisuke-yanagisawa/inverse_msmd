#!/usr/bin/env python3
"""
プロファイル統合CLIスクリプト

複数の単一プロファイルDXファイルを統合・正規化する。
HPC環境でのジョブスケジューラ連携を想定した個別ステップCLI。

使用例:
    python scripts/combine_profiles.py \
        --profile-dxs single_profiles/sys*_A17_ALA_environment.dx \
        --bulk-cnts single_profiles/sys*_A17_ALA_bulk_cnt.txt \
        --output profiles/A17_ALA_profile.dx \
        --eps 0.1 \
        --compress
"""

import argparse
import logging
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="複数の単一プロファイルを統合・正規化する",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--profile-dxs", nargs="+", required=True,
        help="個別プロファイルDXファイルのパス（複数指定可）",
    )
    parser.add_argument(
        "--bulk-cnts", nargs="+", required=True,
        help="バルクカウントファイルのパス（--profile-dxsと同数必要）",
    )
    parser.add_argument(
        "--output", required=True,
        help="出力DXファイルのパス",
    )
    parser.add_argument(
        "--eps", type=float, default=0.1,
        help="バルク正規化時のε (デフォルト: 0.1)",
    )
    parser.add_argument(
        "--compress", action="store_true",
        help="出力ファイルをgzip圧縮する",
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

    # バリデーション
    if len(args.profile_dxs) != len(args.bulk_cnts):
        print(
            f"エラー: DXファイル数({len(args.profile_dxs)})と"
            f"バルクカウントファイル数({len(args.bulk_cnts)})が一致しません"
        )
        return 1

    # 入力ファイルの存在確認
    missing = []
    for path in args.profile_dxs:
        if not Path(path).exists():
            missing.append(f"DX: {path}")
    for path in args.bulk_cnts:
        if not Path(path).exists():
            missing.append(f"バルクカウント: {path}")
    if missing:
        print("エラー: 以下のファイルが見つかりません:")
        for m in missing:
            print(f"  - {m}")
        return 1

    print("=" * 70)
    print("プロファイル統合")
    print("=" * 70)
    print(f"DXファイル数     : {len(args.profile_dxs)}")
    print(f"出力ファイル     : {args.output}")
    print(f"ε                : {args.eps}")
    print(f"gzip圧縮         : {'あり' if args.compress else 'なし'}")
    print()

    from inverse_msmd.profile_generator import combine_profiles, compress_profile

    try:
        output_path = combine_profiles(
            profile_dx_paths=args.profile_dxs,
            bulk_cnt_paths=args.bulk_cnts,
            output_path=args.output,
            eps=args.eps,
        )

        if args.compress:
            output_path = compress_profile(output_path)

        print("\n" + "=" * 70)
        print("完了")
        print("=" * 70)
        print(f"  出力: {output_path}")
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
