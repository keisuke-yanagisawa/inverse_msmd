#!/usr/bin/env python3
# coding: utf-8

"""
PDBファイルとSMIファイルからSDFファイルを作成

PDBファイルから3D座標情報を、SMIファイルから結合情報を取得し、
両方の情報を含むSDFファイルを生成します。

使用例:
    python create_sdf_from_pdb_smi.py input.pdb input.smi output.sdf
    python create_sdf_from_pdb_smi.py --help
"""

import argparse
import sys
from pathlib import Path
from inverse_msmd import read_mol_from_pdb_smi
from rdkit import Chem


def main():
    """メイン関数"""
    parser = argparse.ArgumentParser(
        description='PDBファイルとSMIファイルからSDFファイルを作成します。',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
使用例:
  %(prog)s data/sample_probes/A08.pdb data/sample_probes/A08.smi output/A08.sdf
  %(prog)s input.pdb input.smi output.sdf --name "My Molecule"
  
注意事項:
  - PDBファイルとSMIファイルの原子数が一致している必要があります
  - テンプレートマッチングにより、原子順序が異なっていても自動的にマッチングします
  - RDKitのインストールが必要です: conda install -c conda-forge rdkit
        '''
    )
    
    parser.add_argument('pdb_file', type=str,
                       help='入力PDBファイル（座標情報）')
    parser.add_argument('smi_file', type=str,
                       help='入力SMIファイル（結合情報）')
    parser.add_argument('output_file', type=str,
                       help='出力SDFファイル')
    parser.add_argument('--name', '-n', type=str, default=None,
                       help='分子名（省略時は出力ファイル名から自動設定）')
    
    args = parser.parse_args()
    
    # ファイルの存在確認
    if not Path(args.pdb_file).exists():
        print(f"エラー: PDBファイルが見つかりません: {args.pdb_file}", file=sys.stderr)
        sys.exit(1)
    
    if not Path(args.smi_file).exists():
        print(f"エラー: SMIファイルが見つかりません: {args.smi_file}", file=sys.stderr)
        sys.exit(1)
    
    # 出力ディレクトリを作成
    output_dir = Path(args.output_file).parent
    if output_dir != Path('.'):
        output_dir.mkdir(parents=True, exist_ok=True)
    
    # 分子を読み込み
    mol = read_mol_from_pdb_smi(
        pdb_file=args.pdb_file,
        smi_file=args.smi_file,
        verbose=True
    )
    
    # 分子名を設定
    if args.name:
        mol.SetProp("_Name", args.name)
    else:
        mol.SetProp("_Name", Path(args.output_file).stem)
    
    # SDFファイルとして保存
    print(f"SDFファイルを保存中: {args.output_file}")
    writer = Chem.SDWriter(args.output_file)
    writer.write(mol)
    writer.close()
    
    print("完了！")


if __name__ == '__main__':
    main()