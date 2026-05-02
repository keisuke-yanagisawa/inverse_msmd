#!/usr/bin/env python3
"""旧 ratio 形式 (occupancy/bulk) の RIprofile dx を log 形式に変換する移行スクリプト。

新フォーマットでは ``combine_profiles`` が ``log(occupancy/bulk)`` を直接出力する。
それ以前 (ratio を保存していた版) で生成された dx ファイルが手元にある場合、
本スクリプトで ``np.log`` を当てて新フォーマットに揃える。

判定: dx の全ボクセル値が正なら旧形式とみなして変換する。
負値が混在していれば既に log 形式と判断してスキップする (冪等動作)。
"""
from __future__ import annotations

import argparse
import gzip
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
from gridData import Grid


def migrate_one(path: Path, dry_run: bool = False) -> str:
    g = Grid(str(path))
    if (g.grid <= 0).any():
        return f"skip (already log or non-positive present): {path.name}"
    if dry_run:
        return f"would convert: {path.name} (range [{g.grid.min():.4g}, {g.grid.max():.4g}])"

    g.grid = np.log(g.grid)
    is_gz = path.suffix == ".gz"
    with tempfile.TemporaryDirectory() as td:
        plain = Path(td) / "out.dx"
        g.export(str(plain))
        if is_gz:
            with plain.open("rb") as src, gzip.open(path, "wb", compresslevel=6) as dst:
                shutil.copyfileobj(src, dst)
        else:
            shutil.move(str(plain), str(path))
    return f"converted: {path.name} (log range [{g.grid.min():.3f}, {g.grid.max():.3f}])"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path, help="変換対象の dx / dx.gz ファイル (またはそれを含むディレクトリ)")
    parser.add_argument("--dry-run", action="store_true", help="書き換えずに対象を表示")
    args = parser.parse_args()

    targets: list[Path] = []
    for p in args.paths:
        if p.is_dir():
            targets.extend(sorted(p.glob("*.dx")))
            targets.extend(sorted(p.glob("*.dx.gz")))
        else:
            targets.append(p)

    if not targets:
        print("対象ファイルが見つかりません", file=sys.stderr)
        return 1

    for t in targets:
        print(migrate_one(t, dry_run=args.dry_run))
    return 0


if __name__ == "__main__":
    sys.exit(main())
