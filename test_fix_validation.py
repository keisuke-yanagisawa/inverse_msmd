#!/usr/bin/env python
"""
構造重ね合わせ修正の動作確認スクリプト
"""

import sys
from pathlib import Path

# 簡単な動作確認
try:
    from inverse_msmd.substructure_replacement import integrated_substructure_replacement
    print("✓ モジュールのインポート成功")
    
    # 最小限のテスト実行
    print("\n統合ワークフローをテスト中...")
    results = integrated_substructure_replacement(
        ligand_file="data/atom_matching/4hw3_A_lig.sdf",
        protein_file="data/sample_proteins/4hw3_A.pdb",
        from_file="data/sample_probes/E23",
        to_file="data/sample_probes/E24",
        output_dir="test_output/fix_validation/",
        match_index=0
    )
    
    print(f"\n✓ 処理完了: {len(results)} パターンの結果を生成")
    
    # 出力ファイルの確認
    for i, result in enumerate(results):
        ligand_file = Path(result['ligand_file'])
        protein_file = Path(result['protein_file'])
        
        if ligand_file.exists() and protein_file.exists():
            print(f"✓ パターン {i}: ファイル生成成功")
            print(f"  リガンド: {ligand_file}")
            print(f"  タンパク質: {protein_file}")
        else:
            print(f"✗ パターン {i}: ファイル生成失敗")
            sys.exit(1)
    
    print("\n✓ 全てのテスト成功")
    
except Exception as e:
    print(f"\n✗ エラー発生: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)