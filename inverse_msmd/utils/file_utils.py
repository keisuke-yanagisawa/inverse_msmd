"""
ファイル処理ユーティリティ

このモジュールは、ファイルやディレクトリの操作に関する共通機能を提供します。
ファイルの存在確認、ディレクトリ作成、ファイルペアの検索などを含みます。
"""

from pathlib import Path
from typing import Tuple, Union


def ensure_output_dir(filepath: Union[str, Path]) -> Path:
    """
    出力ファイルのディレクトリを確実に作成します。
    
    Parameters
    ----------
    filepath : str or Path
        出力ファイルのパス
    
    Returns
    -------
    Path
        親ディレクトリのPathオブジェクト
    
    Raises
    ------
    PermissionError
        ディレクトリ作成権限がない場合
    
    Examples
    --------
    >>> output_dir = ensure_output_dir("output/results/data.csv")
    >>> print(output_dir)
    output/results
    
    >>> # ネストされたディレクトリも自動作成
    >>> output_dir = ensure_output_dir("a/b/c/d/file.txt")
    >>> output_dir.exists()
    True
    
    Notes
    -----
    - 親ディレクトリが存在しない場合、自動的に作成されます
    - カレントディレクトリ（'.'）の場合は何も作成しません
    - 既に存在するディレクトリの場合はエラーになりません
    """
    path_obj = Path(filepath)
    parent_dir = path_obj.parent
    
    if parent_dir != Path('.'):
        parent_dir.mkdir(parents=True, exist_ok=True)
    
    return parent_dir


def validate_file_exists(
    filepath: Union[str, Path],
    description: str = "File"
) -> Path:
    """
    ファイルの存在を検証します。
    
    Parameters
    ----------
    filepath : str or Path
        検証するファイルパス
    description : str, default="File"
        エラーメッセージ用の説明
    
    Returns
    -------
    Path
        検証済みのPathオブジェクト
    
    Raises
    ------
    FileNotFoundError
        ファイルが存在しない場合
    
    Examples
    --------
    >>> path = validate_file_exists("data/input.sdf", "Input SDF")
    >>> # ファイルが存在する場合、Pathオブジェクトが返される
    
    >>> # ファイルが存在しない場合
    >>> try:
    ...     path = validate_file_exists("missing.txt", "Missing file")
    ... except FileNotFoundError as e:
    ...     print(e)
    Missing file not found: missing.txt
    
    Notes
    -----
    - ファイルが存在しない場合、descriptionを含むエラーメッセージが表示されます
    - ディレクトリではなくファイルの存在のみを確認します
    """
    path_obj = Path(filepath)
    
    if not path_obj.exists():
        raise FileNotFoundError(
            f"{description} not found: {filepath}"
        )
    
    return path_obj


def find_pdb_smi_pair(
    base_path: Union[str, Path]
) -> Tuple[Path, Path]:
    """
    PDB+SMIファイルペアを検索・検証します。
    
    Parameters
    ----------
    base_path : str or Path
        拡張子なしのベースパス（例: "data/probes/E24"）
    
    Returns
    -------
    Tuple[Path, Path]
        (pdb_path, smi_path) のタプル
    
    Raises
    ------
    FileNotFoundError
        いずれかのファイルが存在しない場合
    
    Examples
    --------
    >>> pdb, smi = find_pdb_smi_pair("data/sample_probes/E24")
    >>> print(pdb, smi)
    data/sample_probes/E24.pdb data/sample_probes/E24.smi
    
    >>> # ファイルが見つからない場合
    >>> try:
    ...     pdb, smi = find_pdb_smi_pair("data/missing/probe")
    ... except FileNotFoundError as e:
    ...     print(e)
    PDB file not found: data/missing/probe.pdb
    SMI file not found: data/missing/probe.smi
    
    Notes
    -----
    - 両方のファイルが存在する場合のみ成功します
    - いずれかのファイルが存在しない場合、すべての欠落ファイルを列挙したエラーが発生します
    - ファイルの内容は検証されません（存在のみ確認）
    """
    base = Path(base_path)
    pdb_path = base.with_suffix('.pdb')
    smi_path = base.with_suffix('.smi')
    
    errors = []
    if not pdb_path.exists():
        errors.append(f"PDB file not found: {pdb_path}")
    if not smi_path.exists():
        errors.append(f"SMI file not found: {smi_path}")
    
    if errors:
        raise FileNotFoundError("\n".join(errors))
    
    return pdb_path, smi_path


def validate_output_path(
    filepath: Union[str, Path],
    must_not_exist: bool = False
) -> Path:
    """
    出力ファイルパスを検証します。
    
    Parameters
    ----------
    filepath : str or Path
        検証する出力ファイルパス
    must_not_exist : bool, default=False
        Trueの場合、ファイルが既に存在する場合にエラーを発生させます
    
    Returns
    -------
    Path
        検証済みのPathオブジェクト
    
    Raises
    ------
    FileExistsError
        must_not_exist=Trueでファイルが既に存在する場合
    
    Examples
    --------
    >>> # 既存ファイルの上書きを許可
    >>> path = validate_output_path("output/result.txt")
    
    >>> # 既存ファイルの上書きを禁止
    >>> try:
    ...     path = validate_output_path("existing.txt", must_not_exist=True)
    ... except FileExistsError as e:
    ...     print(e)
    Output file already exists: existing.txt
    
    Notes
    -----
    - 親ディレクトリが存在しない場合、自動的に作成されます
    - must_not_exist=Falseの場合、既存ファイルは上書きされます
    """
    path_obj = Path(filepath)
    
    # 親ディレクトリを作成
    ensure_output_dir(filepath)
    
    # 既存ファイルのチェック
    if must_not_exist and path_obj.exists():
        raise FileExistsError(f"Output file already exists: {filepath}")
    
    return path_obj


def get_file_extension(filepath: Union[str, Path]) -> str:
    """
    ファイルの拡張子を取得します（ドット付き）。
    
    Parameters
    ----------
    filepath : str or Path
        ファイルパス
    
    Returns
    -------
    str
        拡張子（ドット付き、例: ".txt", ".pdb"）
        拡張子がない場合は空文字列
    
    Examples
    --------
    >>> ext = get_file_extension("data/molecule.sdf")
    >>> print(ext)
    .sdf
    
    >>> ext = get_file_extension("file_without_extension")
    >>> print(ext)
    <empty string>
    
    >>> # 複数のドットを含むファイル名
    >>> ext = get_file_extension("archive.tar.gz")
    >>> print(ext)
    .gz
    
    Notes
    -----
    - 最後のドット以降が拡張子として扱われます
    - 隠しファイル（.bashrcなど）の場合、全体が拡張子として扱われます
    """
    return Path(filepath).suffix


def list_files_with_extension(
    directory: Union[str, Path],
    extension: str,
    recursive: bool = False
) -> list[Path]:
    """
    指定された拡張子を持つファイルをリストアップします。
    
    Parameters
    ----------
    directory : str or Path
        検索するディレクトリ
    extension : str
        ファイル拡張子（ドット付きまたはドットなし、例: ".sdf" または "sdf"）
    recursive : bool, default=False
        サブディレクトリも再帰的に検索するか
    
    Returns
    -------
    list[Path]
        見つかったファイルのPathオブジェクトのリスト
    
    Examples
    --------
    >>> # カレントディレクトリのSDFファイルを検索
    >>> files = list_files_with_extension(".", ".sdf")
    >>> for f in files:
    ...     print(f)
    
    >>> # 再帰的に検索
    >>> files = list_files_with_extension("data", "pdb", recursive=True)
    >>> print(f"Found {len(files)} PDB files")
    
    Notes
    -----
    - 拡張子の大文字小文字は区別されません
    - ディレクトリが存在しない場合、空のリストを返します
    """
    dir_path = Path(directory)
    
    if not dir_path.exists():
        return []
    
    # ドットを追加（必要な場合）
    if not extension.startswith('.'):
        extension = '.' + extension
    
    # 検索パターン
    if recursive:
        pattern = f"**/*{extension}"
    else:
        pattern = f"*{extension}"
    
    # 大文字小文字を区別しない検索
    files = []
    for file in dir_path.glob(pattern):
        if file.is_file() and file.suffix.lower() == extension.lower():
            files.append(file)
    
    return sorted(files)