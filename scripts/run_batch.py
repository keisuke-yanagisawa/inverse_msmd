#!/usr/bin/env python
"""
バッチ処理CLIスクリプト

複数の置換パターン（from_probe, to_probe）に対して、
matching scoreを一括計算するコマンドラインツールです。

使用例:
    python scripts/run_batch.py \\
        --batch-csv experiments/batch_config.csv \\
        --ligand data/atom_matching/4hw3_A_lig.sdf \\
        --protein data/sample_proteins/4hw3_A.pdb \\
        --probe-dir data/sample_probes \\
        --profile-dir data/profiles \\
        --output output/batch_results

    パッケージインストール後は inverse-msmd-batch コマンドとしても利用可能です。
"""

import sys
from inverse_msmd.cli import run_batch_main

if __name__ == "__main__":
    sys.exit(run_batch_main())
