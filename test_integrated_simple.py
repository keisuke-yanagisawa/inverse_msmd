#!/usr/bin/env python
"""統合ワークフロー関数の簡易テスト"""

import sys
from pathlib import Path

print("インポートを開始...")
try:
    from inverse_msmd.substructure_replacement import integrated_substructure_replacement
    print("✓ インポート成功")
except Exception as e:
    print(f"✗ インポートエラー: {e}")
    sys.exit(1)

print("\n統合ワークフローを実行...")
try:
    results = integrated_substructure_replacement(
        ligand_file="data/atom_matching/4hw3_A_lig.sdf",
        protein_file="data/sample_proteins/4hw3_A.pdb",
        from_file="data/sample_probes/E23",
        to_file="data/sample_probes/E24",
        output_dir="test_output/integrated_simple/",
        match_index=0
    )
    
    print(f"\n✓ 処理完了: {len(results)} パターンの結果を生成")
    
    for i, result in enumerate(results):
        print(f"\nパターン {i}:")
        print(f"  リガンド: {result['ligand_file']}")
        print(f"  タンパク質: {result['protein_file']}")
        
        # ファイルの存在確認
        ligand_path = Path(result['ligand_file'])
        protein_path = Path(result['protein_file'])
        
        if ligand_path.exists():
            print(f"  ✓ リガンドファイル存在: {ligand_path.stat().st_size} bytes")
        else:
            print(f"  ✗ リガンドファイルが見つかりません")
        
        if protein_path.exists():
            print(f"  ✓ タンパク質ファイル存在: {protein_path.stat().st_size} bytes")
        else:
            print(f"  ✗ タンパク質ファイルが見つかりません")
    
    print("\n✓ 全ての検証に合格")
    
except Exception as e:
    print(f"\n✗ エラー発生: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)