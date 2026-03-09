#!/usr/bin/env python
"""
3D構造図の生成スクリプト

部分構造置換の結果ディレクトリに対して、PyMOLによる3D構造図を生成します。
既に置換処理が完了した出力ディレクトリを指定して、後から図だけを生成できます。

使用例:
    # 基本的な使用方法（ベストパターンのみ描画、変換行列を自動計算）
    python scripts/render_figures.py \
        --output-dir output/batch_results/1-05 \
        --probe-pdb probe/E21.pdb \
        --from-probe probe/A38 \
        --ligand-sdf target_protein/4gih_B_0X5.sdf \
        --profile-dir profiles \
        --probe-id E21

    # 特定のパターンのみ描画
    python scripts/render_figures.py \
        --output-dir output/batch_results/1-05 \
        --probe-pdb probe/E21.pdb \
        --from-probe probe/A38 \
        --ligand-sdf target_protein/4gih_B_0X5.sdf \
        --profile-dir profiles \
        --probe-id E21 \
        --patterns 0 2 16

    # isomeshレベルを調整
    python scripts/render_figures.py \
        --output-dir output/batch_results/1-05 \
        --probe-pdb probe/E21.pdb \
        --from-probe probe/A38 \
        --ligand-sdf target_protein/4gih_B_0X5.sdf \
        --profile-dir profiles \
        --probe-id E21 \
        --isomesh-level 3.0
"""

import argparse
import sys
from pathlib import Path


def _compute_transform(ligand_sdf, from_probe_base, probe_pdb,
                       match_index=0, pattern_index=0):
    """
    リガンド・置換前プローブ・置換後プローブから変換行列を計算する。

    Parameters
    ----------
    ligand_sdf : str
        元のリガンドSDFファイルパス
    from_probe_base : str
        置換前プローブのベースパス（拡張子なし）
    probe_pdb : str
        置換後プローブのPDBファイルパス
    match_index : int
        リガンド中の部分構造マッチのインデックス
    pattern_index : int
        atom matchingパターンのインデックス

    Returns
    -------
    rot : np.ndarray (3,3)
    tran : np.ndarray (3,)
    """
    from inverse_msmd.substructure_replacement import (
        find_substructure_in_ligand,
        match_substructures,
    )
    from inverse_msmd.utils.mol_utils import read_mol_from_pdb_smi
    from inverse_msmd.utils.bio_utils import SuperImposer
    from rdkit import Chem

    ligand = next(Chem.SDMolSupplier(ligand_sdf))
    from_mol = read_mol_from_pdb_smi(
        f"{from_probe_base}.pdb", f"{from_probe_base}.smi"
    )
    to_mol = read_mol_from_pdb_smi(
        probe_pdb,
        str(Path(probe_pdb).with_suffix(".smi")),
    )

    # 水素除去（integrated_substructure_replacementと同じ前処理）
    ligand_no_h = Chem.RemoveHs(ligand)
    from_no_h = Chem.RemoveHs(from_mol)
    to_no_h = Chem.RemoveHs(to_mol)

    matches = find_substructure_in_ligand(ligand_no_h, from_no_h)
    atom_pairs_list = match_substructures(from_no_h, to_no_h)

    if match_index >= len(matches):
        raise ValueError(f"match_index={match_index} が範囲外（マッチ数: {len(matches)}）")
    if pattern_index >= len(atom_pairs_list):
        raise ValueError(f"pattern_index={pattern_index} が範囲外（パターン数: {len(atom_pairs_list)}）")

    selected_match = matches[match_index]
    atom_pairs = atom_pairs_list[pattern_index]

    # integrated_substructure_replacementと同じSuperImpose計算
    ligand_coords = ligand_no_h.GetConformer().GetPositions()
    to_coords = to_no_h.GetConformer().GetPositions()

    ligand_mcs_indices = [selected_match[i] for i in atom_pairs[0]]
    ligand_mcs_coords = ligand_coords[ligand_mcs_indices]
    e24_mcs_coords = to_coords[atom_pairs[1]]

    si = SuperImposer()
    si.fit(e24_mcs_coords, ligand_mcs_coords)
    return si.rot_, si.tran_


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
        help="置換後プローブのPDBファイルパス（例: probe/E21.pdb）"
    )

    transform = parser.add_argument_group("座標変換引数（プローブ・マップをタンパク質空間に変換）")
    transform.add_argument(
        "--from-probe",
        help="置換前プローブのベースパス（拡張子なし、例: probe/A38）"
    )
    transform.add_argument(
        "--ligand-sdf",
        help="元のリガンドSDFファイルパス"
    )
    transform.add_argument(
        "--match-index", type=int, default=0,
        help="リガンド中の部分構造マッチのインデックス（デフォルト: 0）"
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
        help="描画するパターン番号（例: 0 2 16）。未指定時はベストスコアのパターン"
    )
    optional.add_argument(
        "--isomesh-level", type=float, default=9.0,
        help="isomeshの等値面レベル（デフォルト: 5.0）"
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
    optional.add_argument(
        "--save-pse", action="store_true",
        help="PyMOLセッションファイル(.pse)も保存する"
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
            compute_protein_view,
            render_complex, render_combined,
            render_probe_with_maps, _find_profile_files,
        )
    except ImportError as e:
        print(f"エラー: 必要なモジュールが見つかりません: {e}", file=sys.stderr)
        print("PyMOLがインストールされた環境で実行してください", file=sys.stderr)
        return 1

    ray_size = tuple(args.ray_size)

    # プロファイルファイルの検索
    profile_files = {}
    if args.profile_dir and args.probe_id:
        profile_files = _find_profile_files(args.profile_dir, args.probe_id)
        if profile_files:
            print(f"プロファイル: {list(profile_files.keys())}")
        else:
            print(f"警告: プロファイルファイルが見つかりません: {args.profile_dir}/{args.probe_id}_*")

    # パターンファイルを検索（ベストパターンのインデックスも取得）
    best_pattern_idx = 0
    if args.patterns is not None:
        pdb_files = [
            output_dir / f"pattern_{p}_protein_aligned.pdb"
            for p in args.patterns
        ]
        best_pattern_idx = args.patterns[0]
    else:
        # results.csvがあればベストスコアのパターンのみ描画
        results_csv = output_dir / "results.csv"
        if results_csv.exists():
            import csv
            best_score = float("-inf")
            best_idx = None
            with open(results_csv, "r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    try:
                        score = float(row["score"])
                        idx = int(row["pattern_index"])
                        if score > best_score:
                            best_score = score
                            best_idx = idx
                    except (ValueError, KeyError):
                        continue
            if best_idx is not None:
                pdb_files = [output_dir / f"pattern_{best_idx}_protein_aligned.pdb"]
                best_pattern_idx = best_idx
                print(f"results.csvからベストパターン {best_idx} (score={best_score:.3f}) を選択")
            else:
                pdb_files = sorted(output_dir.glob("pattern_*_protein_aligned.pdb"))
        else:
            pdb_files = sorted(output_dir.glob("pattern_*_protein_aligned.pdb"))

    # 変換行列の計算（ベストパターンのatom mappingを使用）
    transform_rot = None
    transform_tran = None
    if args.from_probe and args.ligand_sdf:
        print(f"変換行列を計算中（match_index={args.match_index}, pattern_index={best_pattern_idx}）...")
        transform_rot, transform_tran = _compute_transform(
            args.ligand_sdf, args.from_probe, str(probe_pdb),
            match_index=args.match_index,
            pattern_index=best_pattern_idx,
        )
        print("  完了")
    elif args.from_probe or args.ligand_sdf:
        print("警告: --from-probe と --ligand-sdf は両方指定する必要があります。"
              "変換なしで描画します。", file=sys.stderr)

    if not pdb_files:
        print("エラー: パターンファイルが見つかりません", file=sys.stderr)
        return 1

    # Panel B: プローブ+マップ（プロファイルがある場合のみ、1回だけ）
    if profile_files:
        panel_b = str(output_dir / "probe_map.png")
        render_probe_with_maps(
            probe_pdb=str(probe_pdb),
            profile_files=profile_files,
            output_png=panel_b,
            view=None,
            isomesh_level=args.isomesh_level,
            ray_size=ray_size,
            dpi=args.dpi,
            save_pse=args.save_pse,
        )
        print(f"  probe_map.png")

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

        # ポケット視点（リガンド方向にカメラ）
        protein_view = compute_protein_view(str(pdb_path), ligand_sdf=str(sdf_path))

        # Panel A: 複合体
        render_complex(
            protein_pdb=str(pdb_path),
            ligand_sdf=str(sdf_path),
            output_png=str(output_dir / f"{pat}_complex.png"),
            view=protein_view,
            ray_size=ray_size,
            dpi=args.dpi,
            save_pse=args.save_pse,
        )

        # Panel C: 統合図（プロファイルがある場合のみ）
        if profile_files:
            render_combined(
                protein_pdb=str(pdb_path),
                ligand_sdf=str(sdf_path),
                probe_pdb=str(probe_pdb),
                profile_files=profile_files,
                output_png=str(output_dir / f"{pat}_combined.png"),
                view=protein_view,
                isomesh_level=args.isomesh_level,
                ray_size=ray_size,
                dpi=args.dpi,
                transform_rot=transform_rot,
                transform_tran=transform_tran,
                save_pse=args.save_pse,
            )

        rendered += 1
        print(f"  {pat} done")

    print(f"\n完了: {rendered} パターンを描画しました")
    return 0


if __name__ == "__main__":
    sys.exit(main())
