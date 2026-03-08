#!/usr/bin/env python
"""
3D構造図の生成スクリプト

部分構造置換の結果ディレクトリに対して、PyMOLによる3D構造図を生成します。
既に置換処理が完了した出力ディレクトリを指定して、後から図だけを生成できます。

使用例:
    # 基本的な使用方法（出力ディレクトリ内の全パターンを描画）
    python scripts/render_figures.py \
        --output-dir output/test_render \
        --probe-pdb data/sample_probes/E24.pdb \
        --profile-dir data/profiles \
        --probe-id E24

    # 特定のパターンのみ描画
    python scripts/render_figures.py \
        --output-dir output/test_render \
        --probe-pdb data/sample_probes/E24.pdb \
        --profile-dir data/profiles \
        --probe-id E24 \
        --patterns 0 2 16

    # isomeshレベルとtilt角度を調整
    python scripts/render_figures.py \
        --output-dir output/test_render \
        --probe-pdb data/sample_probes/E24.pdb \
        --profile-dir data/profiles \
        --probe-id E24 \
        --isomesh-level 3.0 \
        --tilt-deg 30
"""

import argparse
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="部分構造置換結果に対してPyMOLによる3D構造図を生成",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    required = parser.add_argument_group("必須引数")
    required.add_argument(
        "--output-dir", required=True,
        help="置換結果の出力ディレクトリ（pattern_N_*.pdb/sdf があるディレクトリ）"
    )
    required.add_argument(
        "--probe-pdb", required=True,
        help="プローブPDBファイルのパス"
    )

    optional = parser.add_argument_group("オプション引数")
    optional.add_argument(
        "--profile-dir",
        help="プロファイルディレクトリのパス（指定時にプローブ+マップ図と統合図を生成）"
    )
    optional.add_argument(
        "--probe-id",
        help="プローブID（例: E24）。--profile-dir と併用"
    )
    optional.add_argument(
        "--patterns", nargs="+", type=int,
        help="描画するパターン番号（例: 0 2 16）。未指定時は全パターン"
    )
    optional.add_argument(
        "--isomesh-level", type=float, default=5.0,
        help="isomeshの等値面レベル（デフォルト: 5.0）"
    )
    optional.add_argument(
        "--tilt-deg", type=float, default=45.0,
        help="視点のX軸方向傾斜角度（デフォルト: 45.0）"
    )
    optional.add_argument(
        "--distance", type=float, default=50.0,
        help="カメラ距離（デフォルト: 50.0）"
    )
    optional.add_argument(
        "--ray-size", type=int, nargs=2, default=[800, 800],
        metavar=("W", "H"),
        help="レイトレーシング画像サイズ（デフォルト: 800 800）"
    )
    optional.add_argument(
        "--dpi", type=int, default=150,
        help="出力画像のDPI（デフォルト: 150）"
    )

    args = parser.parse_args()

    # 入力検証
    output_dir = Path(args.output_dir)
    if not output_dir.is_dir():
        print(f"エラー: 出力ディレクトリが存在しません: {output_dir}", file=sys.stderr)
        return 1

    probe_pdb = Path(args.probe_pdb)
    if not probe_pdb.exists():
        print(f"エラー: プローブPDBが存在しません: {probe_pdb}", file=sys.stderr)
        return 1

    if args.profile_dir and not args.probe_id:
        print("エラー: --profile-dir を指定する場合は --probe-id も必要です", file=sys.stderr)
        return 1

    # PyMOLのインポート確認
    try:
        from inverse_msmd.pymol_visualization import (
            compute_probe_view, render_complex, render_combined,
            render_probe_with_maps, _find_profile_files,
        )
    except ImportError as e:
        print(f"エラー: 必要なモジュールが見つかりません: {e}", file=sys.stderr)
        print("PyMOLがインストールされた環境で実行してください", file=sys.stderr)
        return 1

    # 視点の計算
    print(f"プローブ {probe_pdb} から視点を計算中...")
    view = compute_probe_view(
        str(probe_pdb),
        distance=args.distance,
        tilt_deg=args.tilt_deg,
    )

    ray_size = tuple(args.ray_size)

    # プロファイルファイルの検索
    profile_files = {}
    if args.profile_dir and args.probe_id:
        profile_files = _find_profile_files(args.profile_dir, args.probe_id)
        if profile_files:
            print(f"プロファイル: {list(profile_files.keys())}")
        else:
            print(f"警告: プロファイルファイルが見つかりません: {args.profile_dir}/{args.probe_id}_*")

    # Panel B: プローブ+マップ（プロファイルがある場合のみ、1回だけ）
    if profile_files:
        panel_b = str(output_dir / "probe_map.png")
        render_probe_with_maps(
            probe_pdb=str(probe_pdb),
            profile_files=profile_files,
            output_png=panel_b,
            view=view,
            isomesh_level=args.isomesh_level,
            ray_size=ray_size,
            dpi=args.dpi,
        )
        print(f"  probe_map.png")

    # パターンファイルを検索
    if args.patterns is not None:
        pdb_files = [
            output_dir / f"pattern_{p}_protein_aligned.pdb"
            for p in args.patterns
        ]
    else:
        pdb_files = sorted(output_dir.glob("pattern_*_protein_aligned.pdb"))

    rendered = 0
    for pdb_path in pdb_files:
        if not pdb_path.exists():
            print(f"  警告: {pdb_path.name} が見つかりません、スキップ")
            continue

        pat = pdb_path.stem.replace("_protein_aligned", "")
        sdf_path = output_dir / f"{pat}_ligand_replaced.sdf"
        if not sdf_path.exists():
            print(f"  警告: {sdf_path.name} が見つかりません、スキップ")
            continue

        # Panel A: 複合体
        render_complex(
            protein_pdb=str(pdb_path),
            ligand_sdf=str(sdf_path),
            output_png=str(output_dir / f"{pat}_complex.png"),
            view=view,
            ray_size=ray_size,
            dpi=args.dpi,
        )

        # Panel C: 統合図（プロファイルがある場合のみ）
        if profile_files:
            render_combined(
                protein_pdb=str(pdb_path),
                ligand_sdf=str(sdf_path),
                probe_pdb=str(probe_pdb),
                profile_files=profile_files,
                output_png=str(output_dir / f"{pat}_combined.png"),
                view=view,
                isomesh_level=args.isomesh_level,
                ray_size=ray_size,
                dpi=args.dpi,
            )

        rendered += 1
        print(f"  {pat} done")

    print(f"\n完了: {rendered} パターンを描画しました")
    return 0


if __name__ == "__main__":
    sys.exit(main())
