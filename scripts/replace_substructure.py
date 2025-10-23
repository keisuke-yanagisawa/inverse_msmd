#!/usr/bin/env python3
"""
部分構造置換スクリプト

SDFファイル内の分子について、指定された部分構造を別の構造で置き換えます。

使用例:
    python scripts/replace_substructure.py \
        --input data/atom_matching/4hw3_A_lig.sdf \
        --from-smiles "Cc1cccc(C)c1Cl" \
        --to-smiles "c1ccc(-c2ccccc2)cc1" \
        --output output/4hw3_A_lig_E23toE24.sdf
"""

import argparse
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import matplotlib
matplotlib.use('Agg')  # GUIなしのバックエンドを使用（高速化）
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def find_and_replace_substructure(mol, from_smiles, to_smiles, verbose=False):
    """
    分子内の部分構造を検索して置き換えます。
    
    Parameters
    ----------
    mol : Chem.Mol
        対象の分子
    from_smiles : str
        検索する部分構造のSMILES
    to_smiles : str
        置き換える部分構造のSMILES
    verbose : bool
        詳細情報を出力するか
        
    Returns
    -------
    Chem.Mol
        部分構造が置き換えられた新しい分子、マッチしない場合はNone
    """
    # 検索する部分構造を作成
    from_mol = Chem.MolFromSmiles(from_smiles)
    if from_mol is None:
        raise ValueError(f"from_smilesが無効です: {from_smiles}")
    
    # 置き換える部分構造を作成
    to_mol = Chem.MolFromSmiles(to_smiles)
    if to_mol is None:
        raise ValueError(f"to_smilesが無効です: {to_smiles}")
    
    # 部分構造を検索
    matches = mol.GetSubstructMatches(from_mol)
    
    if not matches:
        raise ValueError(f"部分構造 '{from_smiles}' が見つかりませんでした")
    
    if verbose:
        print(f"部分構造 '{from_smiles}' を {len(matches)} 箇所で発見")
    
    # 最初のマッチを使用して置き換え
    match = matches[0]
    
    if verbose:
        print(f"マッチした原子インデックス: {match}")
        print(f"置き換え前の原子数: {mol.GetNumAtoms()}")
    
    try:
        # RDKitのReplaceSubstructsを使用
        # この関数は部分構造を置き換えた新しい分子を返します
        new_mol = Chem.ReplaceSubstructs(
            mol,
            from_mol,
            to_mol,
            replaceAll=False  # 最初のマッチのみ
        )[0]
        
        if verbose:
            print(f"置き換え後の原子数: {new_mol.GetNumAtoms()}")
            print("✓ 部分構造の置き換えに成功")
        
        return new_mol
        
    except Exception as e:
        print(f"エラー: 部分構造の置き換えに失敗しました: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(
        description="SDFファイル内の分子の部分構造を置き換えます"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="入力SDFファイル"
    )
    parser.add_argument(
        "--from-smiles",
        required=True,
        help="検索する部分構造のSMILES（例: E23の構造）"
    )
    parser.add_argument(
        "--to-smiles",
        required=True,
        help="置き換える部分構造のSMILES（例: E24の構造）"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="出力SDFファイル"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="詳細情報を表示"
    )
    parser.add_argument(
        "--draw", "-d",
        action="store_true",
        help="置換前後の構造を描画して保存"
    )
    parser.add_argument(
        "--draw-output",
        help="描画ファイルの出力先（デフォルト: 出力SDFファイル名.png）"
    )
    
    args = parser.parse_args()
    
    # 入力ファイルの確認
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"エラー: 入力ファイルが見つかりません: {args.input}")
        return 1
    
    # 出力ディレクトリの作成
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"入力ファイル: {args.input}")
    print(f"検索する部分構造: {args.from_smiles}")
    print(f"置き換える部分構造: {args.to_smiles}")
    print(f"出力ファイル: {args.output}")
    print()
    
    # SDFファイルを読み込み
    supplier = Chem.SDMolSupplier(str(input_path))
    
    # 出力用のSDFライターを作成
    writer = Chem.SDWriter(str(output_path))
    
    mol_count = 0
    replaced_count = 0
    original_mols = []
    replaced_mols = []
    
    for mol in supplier:
        if mol is None:
            continue
        
        mol_count += 1
        
        if args.verbose:
            print(f"\n--- 分子 {mol_count} ---")
            if mol.HasProp("_Name"):
                print(f"名前: {mol.GetProp('_Name')}")
        
        # 部分構造を置き換え
        try:
            new_mol = find_and_replace_substructure(
                mol,
                args.from_smiles,
                args.to_smiles,
                verbose=args.verbose
            )
            
            # 元のプロパティをコピー
            for prop in mol.GetPropNames():
                new_mol.SetProp(prop, mol.GetProp(prop))
            
            writer.write(new_mol)
            replaced_count += 1
            
            # 描画用に分子を保存
            if args.draw:
                original_mols.append(mol)
                replaced_mols.append(new_mol)
            
        except ValueError as e:
            print(f"エラー (分子 {mol_count}): {e}")
            if args.verbose:
                print("処理を中断します")
            writer.close()
            return 1
        except Exception as e:
            print(f"予期しないエラー (分子 {mol_count}): {e}")
            if args.verbose:
                print("処理を中断します")
            writer.close()
            return 1
    
    writer.close()
    
    print(f"\n処理完了:")
    print(f"  処理した分子数: {mol_count}")
    print(f"  置き換えた分子数: {replaced_count}")
    print(f"  出力ファイル: {args.output}")
    
    # 構造を描画
    if args.draw and replaced_count > 0:
        draw_output = args.draw_output
        if draw_output is None:
            # デフォルトの出力先を生成
            draw_output = str(output_path.with_suffix('.png'))
        
        print(f"\n構造を描画中...")
        draw_comparison(original_mols, replaced_mols, draw_output, args.from_smiles, args.to_smiles)
        print(f"  描画ファイル: {draw_output}")
    
    return 0


def draw_comparison(original_mols, replaced_mols, output_file, from_smiles, to_smiles):
    """
    置換前後の分子構造を並べて描画します。
    
    Parameters
    ----------
    original_mols : list of Chem.Mol
        元の分子のリスト
    replaced_mols : list of Chem.Mol
        置換後の分子のリスト
    output_file : str
        出力画像ファイルのパス
    from_smiles : str
        置換前の部分構造のSMILES
    to_smiles : str
        置換後の部分構造のSMILES
    """
    n_mols = len(original_mols)
    
    # 図の設定
    fig = plt.figure(figsize=(12, 4 * n_mols))
    gs = GridSpec(n_mols, 2, figure=fig, hspace=0.3, wspace=0.1)
    
    for i, (orig_mol, repl_mol) in enumerate(zip(original_mols, replaced_mols)):
        # 2D座標を生成（描画用）
        orig_mol_2d = Chem.Mol(orig_mol)
        AllChem.Compute2DCoords(orig_mol_2d)
        
        repl_mol_2d = Chem.Mol(repl_mol)
        AllChem.Compute2DCoords(repl_mol_2d)
        
        # 元の分子を描画
        ax1 = fig.add_subplot(gs[i, 0])
        img1 = Draw.MolToImage(orig_mol_2d, size=(400, 400))
        ax1.imshow(img1)
        ax1.axis('off')
        if i == 0:
            ax1.set_title(f'置換前\n(部分構造: {from_smiles})', fontsize=12, pad=10)
        else:
            ax1.set_title(f'分子 {i+1} - 置換前', fontsize=10)
        
        # 置換後の分子を描画
        ax2 = fig.add_subplot(gs[i, 1])
        img2 = Draw.MolToImage(repl_mol_2d, size=(400, 400))
        ax2.imshow(img2)
        ax2.axis('off')
        if i == 0:
            ax2.set_title(f'置換後\n(置換構造: {to_smiles})', fontsize=12, pad=10)
        else:
            ax2.set_title(f'分子 {i+1} - 置換後', fontsize=10)
    
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    exit(main())