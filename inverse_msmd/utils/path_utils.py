"""
パス処理ユーティリティモジュール

ファイルパスの展開と処理のためのユーティリティ関数を提供します。
"""
import os


def expandpath(path: str) -> str:
    """
    パス文字列を展開します。
    
    チルダ（~）とHOME環境変数、その他の環境変数を展開して、
    完全なパス文字列を返します。
    
    Parameters
    ----------
    path : str
        展開するパス文字列
        ~、$HOME、その他の環境変数（$VAR）を含むことができます
        
    Returns
    -------
    str
        展開されたパス文字列
        
    Examples
    --------
    >>> expandpath("~/data/protein.pdb")
    '/home/username/data/protein.pdb'
    
    >>> expandpath("$HOME/data/protein.pdb")
    '/home/username/data/protein.pdb'
    
    >>> expandpath("$WORK_DIR/results/output.pdb")
    '/path/to/work/results/output.pdb'
    
    Notes
    -----
    この関数は以下の順序で展開を行います：
    1. os.path.expanduser() によりチルダ（~）を展開
    2. os.path.expandvars() により環境変数を展開
    
    環境変数が定義されていない場合、そのまま（展開せずに）返されます。
    """
    path = os.path.expanduser(path)
    return os.path.expandvars(path)