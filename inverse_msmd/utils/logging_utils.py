"""
ロギングユーティリティ

このモジュールは、プロジェクト全体で統一されたロギング機能を提供します。
verbose出力の制御やログレベルの管理を行います。
"""

from typing import Literal

LogLevel = Literal["INFO", "WARNING", "ERROR", "DEBUG"]


def verbose_print(
    message: str,
    verbose: bool = False,
    level: LogLevel = "INFO",
    prefix: bool = True
) -> None:
    """
    Verbose出力ヘルパー関数。
    
    Parameters
    ----------
    message : str
        出力メッセージ
    verbose : bool, default=False
        verbose出力を行うか
    level : {"INFO", "WARNING", "ERROR", "DEBUG"}, default="INFO"
        ログレベル
    prefix : bool, default=True
        ログレベルプレフィックスを付けるか
    
    Examples
    --------
    >>> verbose_print("処理を開始します", verbose=True)
    INFO: 処理を開始します
    
    >>> verbose_print("警告メッセージ", verbose=True, level="WARNING")
    WARNING: 警告メッセージ
    
    >>> verbose_print("プレフィックスなし", verbose=True, prefix=False)
    プレフィックスなし
    
    >>> # verboseがFalseの場合は何も出力されない
    >>> verbose_print("このッセージは表示されません", verbose=False)
    
    Notes
    -----
    - verboseがFalseの場合、何も出力されません
    - prefixがTrueの場合、メッセージの前にログレベルが付加されます
    - 標準出力に直接出力されます（ファイルへのログ出力はサポートしていません）
    """
    if not verbose:
        return
    
    if prefix:
        prefix_str = f"{level}: "
    else:
        prefix_str = ""
    
    print(f"{prefix_str}{message}")


def format_progress_message(
    current: int,
    total: int,
    message: str = "",
    width: int = 50
) -> str:
    """
    進捗状況を表す文字列をフォーマットします。
    
    Parameters
    ----------
    current : int
        現在の進捗（1始まり）
    total : int
        全体の数
    message : str, default=""
        追加メッセージ
    width : int, default=50
        プログレスバーの幅（文字数）
    
    Returns
    -------
    str
        フォーマットされた進捗メッセージ
    
    Examples
    --------
    >>> msg = format_progress_message(5, 10, "処理中")
    >>> print(msg)
    [=========================.........................] 50.0% (5/10) 処理中
    
    >>> msg = format_progress_message(10, 10, "完了")
    >>> print(msg)
    [==================================================] 100.0% (10/10) 完了
    
    Notes
    -----
    - プログレスバーは'='で進捗を、'.'で残りを表示します
    - パーセンテージは小数点第1位まで表示されます
    """
    percentage = (current / total) * 100
    filled = int((current / total) * width)
    bar = "=" * filled + "." * (width - filled)
    
    progress_str = f"[{bar}] {percentage:.1f}% ({current}/{total})"
    if message:
        progress_str += f" {message}"
    
    return progress_str


def format_section_header(
    title: str,
    char: str = "=",
    width: int = 70
) -> str:
    """
    セクションヘッダーをフォーマットします。
    
    Parameters
    ----------
    title : str
        セクションタイトル
    char : str, default="="
        区切り文字
    width : int, default=70
        全体の幅
    
    Returns
    -------
    str
        フォーマットされたヘッダー文字列（複数行）
    
    Examples
    --------
    >>> header = format_section_header("データ処理")
    >>> print(header)
    ======================================================================
    データ処理
    ======================================================================
    
    >>> header = format_section_header("サブセクション", char="-", width=50)
    >>> print(header)
    --------------------------------------------------
    サブセクション
    --------------------------------------------------
    
    Notes
    -----
    - ヘッダーは3行（区切り線、タイトル、区切り線）で構成されます
    - タイトルが幅より長い場合、区切り線はタイトルの長さに合わせられます
    """
    separator = char * width
    if len(title) > width:
        separator = char * len(title)
    
    return f"{separator}\n{title}\n{separator}"


def format_key_value(
    key: str,
    value: str,
    key_width: int = 20,
    separator: str = ": "
) -> str:
    """
    キーと値のペアをフォーマットします。
    
    Parameters
    ----------
    key : str
        キー
    value : str
        値
    key_width : int, default=20
        キーの表示幅（左詰め）
    separator : str, default=": "
        キーと値の区切り文字
    
    Returns
    -------
    str
        フォーマットされた文字列
    
    Examples
    --------
    >>> print(format_key_value("ファイル名", "data.sdf"))
    ファイル名          : data.sdf
    
    >>> print(format_key_value("原子数", "42", key_width=15))
    原子数         : 42
    
    >>> print(format_key_value("スコア", "95.3", separator=" = "))
    スコア              = 95.3
    
    Notes
    -----
    - キーは指定された幅で左詰めされます
    - キーが指定幅より長い場合、そのまま表示されます
    """
    formatted_key = key.ljust(key_width)
    return f"{formatted_key}{separator}{value}"