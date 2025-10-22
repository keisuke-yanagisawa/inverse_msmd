"""
同位体ラベル付与スクリプト

SDFファイルの分子に対して、SMARTS記法で指定した部分構造に同位体番号を付与します。

複数の該当領域が存在する場合:
- デフォルト: 全ての該当領域に適用
- --interactive (-i): 対話的に選択
- --match-index N (-m N): N番目のマッチのみに適用（0始まり）

使用例:
    # 全マッチに適用（デフォルト）
    python add_isotope_labels.py input.sdf output.sdf "c1ccccc1" 13
    
    # 対話的に選択
    python add_isotope_labels.py input.sdf output.sdf "c1ccccc1" 13 --interactive
    
    # 2番目のマッチのみに適用
    python add_isotope_labels.py input.sdf output.sdf "c1ccccc1" 13 --match-index 1
"""

import sys
import argparse
import os
from rdkit import Chem
from rdkit.Chem import AllChem, Draw


def draw_molecule_with_isotopes(mol_before, mol_after, smarts_pattern, output_file,
                                  image_size=(800, 400)):
    """
    同位体ラベル付与前後の分子を並べて描画する
    
    Parameters
    ----------
    mol_before : rdkit.Chem.Mol
        元の分子
    mol_after : rdkit.Chem.Mol
        同位体番号が付与された分子
    smarts_pattern : str
        部分構造を指定するSMARTSパターン
    output_file : str
        出力画像ファイルのパス
    image_size : tuple, optional
        画像サイズ (width, height)
    """
    # SMARTSパターンでマッチする原子を取得
    pattern = Chem.MolFromSmarts(smarts_pattern)
    if pattern is None:
        print(f"警告: 描画用のSMARTSパターンが無効です: {smarts_pattern}")
        return
    
    matches = mol_after.GetSubstructMatches(pattern)
    if not matches:
        print(f"警告: 描画用にマッチする部分構造が見つかりませんでした")
        return
    
    # マッチした全ての原子のインデックスを収集
    highlight_atoms = set()
    for match in matches:
        highlight_atoms.update(match)
    
    # 2D座標を生成（視認性向上のため）
    mol_before_2d = Chem.Mol(mol_before)
    mol_after_2d = Chem.Mol(mol_after)
    
    # 2D座標を計算
    AllChem.Compute2DCoords(mol_before_2d)
    AllChem.Compute2DCoords(mol_after_2d)
    
    # 2つの分子を並べて描画
    # 元の分子と同位体ラベル付き分子を比較表示
    img = Draw.MolsToGridImage(
        [mol_before_2d, mol_after_2d],
        molsPerRow=2,
        subImgSize=(image_size[0]//2, image_size[1]),
        legends=["元の分子", f"同位体ラベル付与後 (ISO={mol_after_2d.GetAtomWithIdx(list(highlight_atoms)[0]).GetIsotope()})"],
        highlightAtomLists=[list(highlight_atoms), list(highlight_atoms)],
        returnPNG=False
    )
    
    # 画像を保存
    img.save(output_file)
    print(f"  → 描画画像を保存: {output_file}")


def add_isotope_labels(mol, smarts_pattern, isotope_number, match_index=None):
    """
    分子の指定した部分構造に同位体番号を付与する
    
    Parameters
    ----------
    mol : rdkit.Chem.Mol
        対象の分子
    smarts_pattern : str
        部分構造を指定するSMARTSパターン
    isotope_number : int
        付与する同位体番号
    match_index : int, optional
        使用するマッチのインデックス（0始まり）。Noneの場合は全マッチに適用
    
    Returns
    -------
    tuple
        (同位体番号が付与された分子, ラベル付与された原子のセット, マッチ数)
    """
    # 分子のコピーを作成（元の分子を変更しない）
    mol_copy = Chem.RWMol(mol)
    
    # SMARTSパターンから部分構造を作成
    pattern = Chem.MolFromSmarts(smarts_pattern)
    if pattern is None:
        raise ValueError(f"無効なSMARTSパターン: {smarts_pattern}")
    
    # 部分構造にマッチする原子のインデックスを取得
    matches = mol_copy.GetSubstructMatches(pattern)
    
    if not matches:
        print(f"警告: SMARTSパターン '{smarts_pattern}' にマッチする部分構造が見つかりませんでした")
        return mol_copy.GetMol(), set(), 0
    
    # match_indexが指定されている場合は、そのマッチのみを使用
    if match_index is not None:
        if match_index < 0 or match_index >= len(matches):
            print(f"警告: match_index {match_index} は範囲外です（0-{len(matches)-1}）")
            return mol_copy.GetMol(), set(), len(matches)
        matches_to_use = [matches[match_index]]
        print(f"マッチ {match_index} を使用（全マッチ数: {len(matches)}）")
    else:
        matches_to_use = matches
        print(f"マッチ数: {len(matches)}")
    
    # マッチした原子に同位体番号を設定
    labeled_atoms = set()
    for match in matches_to_use:
        for atom_idx in match:
            atom = mol_copy.GetAtomWithIdx(atom_idx)
            atom.SetIsotope(isotope_number)
            labeled_atoms.add(atom_idx)
    
    print(f"ラベル付与原子数: {len(labeled_atoms)}")
    
    return mol_copy.GetMol(), labeled_atoms, len(matches)


def process_sdf_file(input_file, output_file, smarts_pattern, isotope_number,
                     visualize=False, image_output_dir=None, interactive=False,
                     match_index=None):
    """
    SDFファイル内の分子に同位体ラベルを付与する
    
    Parameters
    ----------
    input_file : str
        入力SDFファイルのパス
    output_file : str
        出力SDFファイルのパス
    smarts_pattern : str
        部分構造を指定するSMARTSパターン
    isotope_number : int
        付与する同位体番号
    visualize : bool, optional
        Trueの場合、同位体ラベル付与前後の分子構造を描画（デフォルト: False）
    image_output_dir : str, optional
        画像の出力ディレクトリ（Noneの場合は入力ファイルと同じディレクトリ）
    interactive : bool, optional
        Trueの場合、複数マッチ時に対話的に選択（デフォルト: False）
    match_index : int, optional
        使用するマッチのインデックス（0始まり）。Noneの場合は全マッチに適用
    """
    # SDFファイルから分子を読み込む
    suppl = Chem.SDMolSupplier(input_file, removeHs=False)
    writer = Chem.SDWriter(output_file)
    
    # 画像出力ディレクトリの設定
    if visualize and image_output_dir is None:
        image_output_dir = os.path.dirname(output_file) or "."
    
    mol_count = 0
    processed_count = 0
    
    for mol in suppl:
        mol_count += 1
        if mol is None:
            print(f"警告: 分子 {mol_count} の読み込みに失敗しました")
            continue
        
        print(f"\n分子 {mol_count} を処理中...")
        
        # 対話的モードの場合、まずマッチ数を確認
        selected_match_index = match_index
        if interactive and match_index is None:
            # マッチ数を取得するため一度実行
            _, _, num_matches = add_isotope_labels(mol, smarts_pattern, isotope_number, match_index=None)
            
            if num_matches > 1:
                print(f"\n{num_matches} 個のマッチが見つかりました。")
                print("どのマッチを使用しますか？")
                print("  -1: 全てのマッチに適用")
                for i in range(num_matches):
                    print(f"  {i}: マッチ {i}")
                
                while True:
                    try:
                        choice = input("選択 (0-{}、または -1): ".format(num_matches - 1))
                        choice_int = int(choice)
                        if choice_int == -1:
                            selected_match_index = None
                            break
                        elif 0 <= choice_int < num_matches:
                            selected_match_index = choice_int
                            break
                        else:
                            print(f"無効な選択です。0-{num_matches-1} または -1 を入力してください。")
                    except ValueError:
                        print("数値を入力してください。")
                    except (EOFError, KeyboardInterrupt):
                        print("\n処理を中断しました。")
                        writer.close()
                        return
        
        # 同位体ラベルを付与
        labeled_mol, labeled_atoms, _ = add_isotope_labels(
            mol, smarts_pattern, isotope_number, match_index=selected_match_index
        )
        
        # 描画処理
        if visualize and labeled_atoms:
            # 出力ファイル名のベース名を取得
            base_name = os.path.splitext(os.path.basename(output_file))[0]
            image_file = os.path.join(image_output_dir, f"{base_name}_mol{mol_count}.png")
            draw_molecule_with_isotopes(mol, labeled_mol, smarts_pattern, image_file)
        
        # 結果を書き込む
        writer.write(labeled_mol)
        processed_count += 1
    
    writer.close()
    
    print(f"\n処理完了: {processed_count}/{mol_count} 個の分子を処理しました")
    print(f"出力ファイル: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="SDFファイルの分子にSMARTS記法で指定した部分構造に同位体番号を付与します",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用例:
  # ベンゼン環の炭素に同位体番号13を付与（全マッチ）
  python add_isotope_labels.py input.sdf output.sdf "c1ccccc1" 13
  
  # 対話的モードで複数マッチから選択
  python add_isotope_labels.py input.sdf output.sdf "c1ccccc1" 13 --interactive
  
  # 2番目のマッチのみに適用
  python add_isotope_labels.py input.sdf output.sdf "c1ccccc1" 13 --match-index 1
  
  # カルボニル基の炭素に同位体番号13、酸素に同位体番号18を付与（2段階で実行）
  python add_isotope_labels.py input.sdf temp.sdf "[C]=[O]" 13
  python add_isotope_labels.py temp.sdf output.sdf "[O]=[C]" 18
  
  # アミン基の窒素に同位体番号15を付与
  python add_isotope_labels.py input.sdf output.sdf "[NH2]" 15
        """
    )
    
    parser.add_argument("input", help="入力SDFファイル")
    parser.add_argument("output", help="出力SDFファイル")
    parser.add_argument("smarts", help="SMARTS記法による部分構造パターン")
    parser.add_argument("isotope", type=int, help="付与する同位体番号")
    parser.add_argument("--visualize", "-v", action="store_true",
                       help="同位体ラベル付与前後の分子構造を画像で出力")
    parser.add_argument("--image-dir", "-d", default=None,
                       help="画像ファイルの出力ディレクトリ（デフォルト: 出力SDFと同じディレクトリ）")
    
    parser.add_argument("--interactive", "-i", action="store_true",
                       help="複数マッチ時に対話的にマッチを選択する")
    parser.add_argument("--match-index", "-m", type=int, default=None,
                       help="使用するマッチのインデックス（0始まり）。指定しない場合は全マッチに適用")
    
    args = parser.parse_args()
    
    # 入力ファイルの存在確認
    try:
        with open(args.input, 'r') as f:
            pass
    except FileNotFoundError:
        print(f"エラー: 入力ファイル '{args.input}' が見つかりません")
        sys.exit(1)
    
    # 処理を実行
    try:
        process_sdf_file(args.input, args.output, args.smarts, args.isotope,
                        visualize=args.visualize, image_output_dir=args.image_dir,
                        interactive=args.interactive, match_index=args.match_index)
    except Exception as e:
        print(f"エラー: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()