#!/usr/bin/env python
"""
単体部分構造置換+スコア計算CLIスクリプト

リガンド中の部分構造を別の部分構造で置換し、タンパク質構造の座標変換と
プロファイルマッチングスコア計算を行うコマンドラインツールです。

使用例:
    python scripts/run_single.py \\
        --ligand data/atom_matching/4hw3_A_lig.sdf \\
        --protein data/sample_proteins/4hw3_A.pdb \\
        --from-probe data/sample_probes/E23 \\
        --to-probe data/sample_probes/E24 \\
        --output output/results \\
        --profile-dir data/profiles \\
        --probe-id E24

    パッケージインストール後は inverse-msmd-run コマンドとしても利用可能です。
"""

import sys
from inverse_msmd.cli import run_single_main

if __name__ == "__main__":
    sys.exit(run_single_main())
