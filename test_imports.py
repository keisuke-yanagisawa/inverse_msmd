#!/usr/bin/env python3
"""基本インポートのテスト"""

# モジュールがインポートできることを確認
try:
    from inverse_msmd import substructure_replacement
    print("✓ モジュールのインポート成功")
except ImportError as e:
    print(f"✗ インポートエラー: {e}")
    exit(1)

# 必要な依存関係の確認
try:
    from rdkit import Chem
    from Bio.PDB import Structure
    import numpy as np
    print("✓ 全ての依存関係がインポート可能")
except ImportError as e:
    print(f"✗ 依存関係エラー: {e}")
    exit(1)

# 関数スタブが定義されているか確認
functions_to_check = [
    'find_substructure_in_ligand',
    'visualize_multiple_matches',
    'match_substructures',
    'calculate_transformation',
    'apply_transformation_to_protein',
    'replace_ligand_substructure',
    'integrated_substructure_replacement'
]

print("\n関数スタブの確認:")
for func_name in functions_to_check:
    if hasattr(substructure_replacement, func_name):
        print(f"  ✓ {func_name}")
    else:
        print(f"  ✗ {func_name} が見つかりません")
        exit(1)

print("\n✓ 全てのインポートチェック完了")