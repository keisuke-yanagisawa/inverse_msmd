#!/usr/bin/env python
"""
バッチ処理の使用例

このスクリプトは、複数の置換パターンに対してmatching scoreを
一括計算するバッチ処理機能の使用方法を示します。
"""

from inverse_msmd.batch_processing import run_batch_processing
import sys
from pathlib import Path

def main():
    """バッチ処理を実行するメイン関数"""
    
    # プロジェクトのベースディレクトリ
    base_dir = Path(__file__).parent.parent
    
    # 入力ファイルのパス
    batch_csv = base_dir / "examples" / "batch_config_sample.csv"
    ligand_file = base_dir / "data" / "atom_matching" / "4hw3_A_lig.sdf"
    protein_file = base_dir / "data" / "sample_proteins" / "4hw3_A.pdb"
    probe_base_dir = base_dir / "data" / "sample_probes"
    profile_base_dir = base_dir / "data" / "profiles"
    output_base_dir = base_dir / "output" / "batch_results"
    
    # ファイルの存在確認
    required_files = {
        "バッチCSV": batch_csv,
        "リガンドファイル": ligand_file,
        "タンパク質ファイル": protein_file,
    }
    
    missing_files = []
    for name, filepath in required_files.items():
        if not filepath.exists():
            missing_files.append(f"{name}: {filepath}")
    
    if missing_files:
        print("エラー: 以下のファイルが見つかりません:")
        for item in missing_files:
            print(f"  - {item}")
        print("\nこのスクリプトを実行する前に、必要なデータファイルを準備してください。")
        sys.exit(1)
    
    print("=" * 70)
    print("バッチ処理の実行例")
    print("=" * 70)
    print(f"バッチCSV: {batch_csv}")
    print(f"リガンドファイル: {ligand_file}")
    print(f"タンパク質ファイル: {protein_file}")
    print(f"プローブディレクトリ: {probe_base_dir}")
    print(f"プロファイルディレクトリ: {profile_base_dir}")
    print(f"出力ディレクトリ: {output_base_dir}")
    print()
    
    # バッチ処理を実行
    try:
        result = run_batch_processing(
            batch_csv=str(batch_csv),
            ligand_file=str(ligand_file),
            protein_file=str(protein_file),
            probe_base_dir=str(probe_base_dir),
            profile_base_dir=str(profile_base_dir),
            output_base_dir=str(output_base_dir),
            parallel=False,  # Phase 2で並列処理を実装予定
            continue_on_error=True
        )
        
        # 結果のサマリーを表示
        print("\n" + "=" * 70)
        print("実行結果サマリー")
        print("=" * 70)
        print(f"総ジョブ数: {result.total_jobs}")
        print(f"成功: {result.num_success}")
        print(f"失敗: {result.num_failed}")
        print(f"スキップ: {result.num_skipped}")
        print(f"実行時間: {result.total_execution_time:.2f}秒")
        
        # 成功したジョブの詳細
        if result.num_success > 0:
            print("\n" + "=" * 70)
            print("成功したジョブの詳細")
            print("=" * 70)
            for job_result in result.job_results:
                if job_result.status == "success":
                    print(f"\n[{job_result.job_id}]")
                    print(f"  {job_result.from_probe} → {job_result.to_probe} (index={job_result.match_index})")
                    print(f"  パターン数: {job_result.num_patterns}")
                    if job_result.best_score is not None:
                        print(f"  最高スコア: {job_result.best_score:.2f} (パターン{job_result.best_pattern_index})")
                    print(f"  実行時間: {job_result.execution_time:.2f}秒")
        
        # 失敗したジョブの詳細
        if result.num_failed > 0:
            print("\n" + "=" * 70)
            print("失敗したジョブの詳細")
            print("=" * 70)
            for job_result in result.job_results:
                if job_result.status == "failed":
                    print(f"\n[{job_result.job_id}]")
                    print(f"  エラータイプ: {job_result.error_type}")
                    print(f"  エラーメッセージ: {job_result.error_message}")
        
        # 出力ファイルの案内
        print("\n" + "=" * 70)
        print("出力ファイル")
        print("=" * 70)
        print(f"サマリーCSV: {output_base_dir}/batch_summary.csv")
        print(f"サマリーJSON: {output_base_dir}/batch_summary.json")
        print(f"実行ログ: {output_base_dir}/batch_execution.log")
        print()
        print("各ジョブの詳細結果は以下のディレクトリに保存されています:")
        for job_result in result.job_results:
            if job_result.status == "success":
                print(f"  - {output_base_dir}/{job_result.job_id}/")
        
        return 0
        
    except Exception as e:
        print(f"\nエラー: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())