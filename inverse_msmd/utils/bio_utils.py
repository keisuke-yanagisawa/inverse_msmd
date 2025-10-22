# coding: utf-8

"""
BioPythonユーティリティモジュール

構造操作と重ね合わせのためのBioPythonユーティリティを提供します。

このモジュールは、Bio/PDB.pyとBio/sklearn_interface.pyを統合したものです。
version: 1.0.0
Authors: Keisuke Yanagisawa
"""
import collections
import gzip
import os
from typing import Any, Callable, List, Literal, Optional, Union
import warnings
import numpy as np
from Bio import PDB as BioPDB
from Bio.PDB import PDBExceptions
from Bio.SVDSuperimposer import SVDSuperimposer
from sklearn.utils.validation import check_is_fitted
from sklearn.base import TransformerMixin, BaseEstimator
from collections.abc import Iterable
import tempfile
import io
from Bio.PDB.Atom import Atom
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
import numpy.typing as npt

from .spatial_utils import estimate_volume
from .path_utils import expandpath


# ============================================================================
# SuperImposer (from sklearn_interface.py)
# ============================================================================

class SuperImposer(TransformerMixin, BaseEstimator):
    """
    構造重ね合わせクラス
    
    BioPythonのSVDSuperimposerをscikit-learnのインターフェースでラップしたクラスです。
    fit/transform/inverse_transformメソッドを提供し、scikit-learnのパイプラインと
    互換性があります。
    
    Attributes
    ----------
    rot_ : np.ndarray
        回転行列（fit後に設定）
    tran_ : np.ndarray
        並進ベクトル（fit後に設定）
    
    Examples
    --------
    >>> si = SuperImposer()
    >>> si.fit(moving_coords, target_coords)
    >>> transformed = si.transform(coords_to_move)
    """

    rot_: npt.NDArray[np.float_]
    tran_: npt.NDArray[np.float_]

    def __init__(self):
        """SuperImposerを初期化します。"""
        pass

    def _reset(self):
        """内部状態をリセットします。"""
        if hasattr(self, "rot_"):
            del self.rot_
            del self.tran_

    def _superimpose(self, coords: npt.ArrayLike, reference_coords: npt.ArrayLike) -> None:
        """
        SVD法により重ね合わせパラメータを計算します。
        
        Parameters
        ----------
        coords : array-like
            移動させる座標群
        reference_coords : array-like
            参照座標群（重ね合わせ先）
        """
        sup = SVDSuperimposer()
        sup.set(reference_coords, coords)
        sup.run()
        self.rot_, self.tran_ = sup.get_rotran()  # type: ignore

    def fit(self, coords: npt.ArrayLike, reference_coords: npt.ArrayLike) -> "SuperImposer":
        """
        重ね合わせパラメータを計算します。
        
        与えられた2つの点群をなるべく重ねるような並進・回転移動を算出します。
        それぞれの点群は対応関係があることを仮定します（i番目の点同士が
        対応するように重ね合わせます）。

        Parameters
        ----------
        coords : array-like, shape (n_points, 3)
            重ね合わせのために移動させる点群
        reference_coords : array-like, shape (n_points, 3)
            重ね合わせ先の点群（参照座標）

        Returns
        -------
        SuperImposer
            fit済みのオブジェクト（メソッドチェーン用）
            
        Notes
        -----
        このメソッドは、最小二乗法により最適な回転行列と並進ベクトルを計算します。
        計算結果は rot_ および tran_ 属性に格納されます。
        """
        self._reset()
        self._superimpose(coords, reference_coords)
        return self

    def transform(self, coords: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """
        座標を変換します。
        
        fit()で計算された並進・回転に基づいて、与えられた点群を移動させます。

        Parameters
        ----------
        coords : np.ndarray, shape (n_points, 3)
            移動させる点群
            
        Returns
        -------
        np.ndarray, shape (n_points, 3)
            変換後の座標
            
        Raises
        ------
        sklearn.exceptions.NotFittedError
            fitメソッドが事前に呼ばれていない場合
        """
        check_is_fitted(self)
        coords = np.array(coords)
        return np.dot(coords, self.rot_) + self.tran_

    def inverse_transform(self, coords: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """
        逆変換を実行します。
        
        transform()の逆操作を行い、変換後の座標を元の座標系に戻します。

        Parameters
        ----------
        coords : np.ndarray, shape (n_points, 3)
            transform()した後の点群

        Returns
        -------
        np.ndarray, shape (n_points, 3)
            transform()する前の座標系での点群座標
            
        Raises
        ------
        sklearn.exceptions.NotFittedError
            fitメソッドが事前に呼ばれていない場合
        """
        coords = np.array(coords)
        check_is_fitted(self)
        return np.dot(coords - self.tran_, np.linalg.inv(self.rot_))


# ============================================================================
# PDB utilities (from PDB.py)
# ============================================================================

def get_structure(filepath: str, structname="") -> Structure:
    """
    PDBファイルを読み込みます。
    
    gzip圧縮されたファイル（.gz拡張子）にも対応しています。

    Parameters
    ----------
    filepath : str
        PDBファイルへのパス（.pdbまたは.pdb.gz）
    structname : str, optional
        構造に付ける名前（デフォルト: 空文字列）

    Returns
    -------
    Bio.PDB.Structure.Structure
        PDBファイルのStructureオブジェクト
        
    Examples
    --------
    >>> structure = get_structure("protein.pdb")
    >>> structure_gz = get_structure("protein.pdb.gz")
    """
    filepath = expandpath(filepath)
    if filepath.endswith(".gz"):
        fileobj = gzip.open(filepath, "rt")
    else:
        fileobj = open(filepath)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBExceptions.PDBConstructionWarning)
        return BioPDB.PDBParser(QUIET=True).get_structure(structname, fileobj)


def get_atom_attr(atom: Atom,
                  attr: Literal["resid", "resname", "coord", "element", "fullname"]
                  ) -> Union[int, str, npt.NDArray[np.float_], tuple]:
    """
    原子オブジェクトから属性を取得します。
    
    Bio.PDB.Atomオブジェクトから指定された属性を取得します。
    現在サポートされている属性: "resid", "resname", "coord", "element", "fullname"

    Parameters
    ----------
    atom : Bio.PDB.Atom.Atom
        原子オブジェクト
    attr : str
        取得する属性名
        - "resid": 残基番号
        - "resname": 残基名
        - "coord": 座標
        - "element": 元素記号
        - "fullname": 完全な原子名

    Returns
    -------
    int or str or np.ndarray or tuple
        指定された属性の値
        - resid: int
        - resname: str
        - coord: np.ndarray (shape: (3,))
        - element: str
        - fullname: tuple

    Raises
    ------
    NotImplementedError
        サポートされていない属性が指定された場合
    """

    if attr == "resid":
        return get_resi(atom)
    elif attr == "resname":
        return get_resname(atom)
    elif attr == "coord":
        return atom.get_coord()
    elif attr == "element":
        return atom.element
    elif attr == "fullname":
        return atom.fullname
    else:
        raise NotImplementedError(f"属性 {attr} はまだサポートされていません。")


def get_attr(model: Union[Structure, Model],
             attr: Literal["resid", "resname", "coord", "element", "fullname"],
             sele: Optional[Callable[[Atom], bool]] = None
             ) -> npt.NDArray[Any]:
    """
    モデルまたは構造から属性を取得します。
    
    Bio.PDB.ModelまたはStructureオブジェクトから、すべての原子（または選択された原子）の
    指定された属性を取得します。
    現在サポートされている属性: "resid", "resname", "coord", "element", "fullname"

    Parameters
    ----------
    model : Bio.PDB.Model.Model or Bio.PDB.Structure.Structure
        モデルまたは構造オブジェクト
    attr : str
        取得する属性名
    sele : function, optional
        原子選択関数。Atomオブジェクトを引数に取り、boolを返す。
        指定しない場合はすべての原子が選択されます。

    Returns
    -------
    np.ndarray
        選択されたすべての原子の属性値の配列

    Raises
    ------
    NotImplementedError
        サポートされていない属性が指定された場合
        
    Examples
    --------
    >>> protein = get_structure("protein.pdb")
    >>> coords = get_attr(protein, "coord")
    >>> ca_coords = get_attr(protein, "coord", sele=lambda a: a.get_name() == "CA")
    """

    data = []
    for atom in model.get_atoms():
        if sele is None or sele(atom):
            data.append(get_atom_attr(atom, attr))
    return np.array(data)


def get_resname(atom: Atom) -> str:
    """
    原子から残基名を取得します。

    Parameters
    ----------
    atom : Bio.PDB.Atom.Atom
        原子オブジェクト

    Returns
    -------
    str
        原子が属する残基の残基名（3文字コード）
        
    Examples
    --------
    >>> resname = get_resname(atom)  # "ALA", "GLY" など
    """
    return atom.get_parent().get_resname()  # type: ignore


def get_resi(atom: Atom) -> int:
    """
    原子から残基番号（残基ID）を取得します。

    Parameters
    ----------
    atom : Bio.PDB.Atom.Atom
        原子オブジェクト

    Returns
    -------
    int
        原子が属する残基の残基番号
        
    Examples
    --------
    >>> resi = get_resi(atom)  # 1, 2, 3, ...
    """
    return atom.get_full_id()[3][1]  # type: ignore


def set_attr(model: Model, attr: str, lst: npt.NDArray, sele=None) -> None:
    """
    モデルに属性を設定します。
    
    Bio.PDB.Modelオブジェクトの原子に対して属性を設定します。
    現在サポートされている属性: "coord"のみ

    Parameters
    ----------
    model : Bio.PDB.Model.Model
        モデルオブジェクト
    attr : str
        設定する属性名（現在は"coord"のみ対応）
    lst : np.ndarray
        設定する値の配列
    sele : function, optional
        原子選択関数。指定しない場合はすべての原子が対象

    Raises
    ------
    NotImplementedError
        "coord"以外の属性が指定された場合
        
    Notes
    -----
    lstの長さと選択された原子数が一致している必要があります。
    
    Examples
    --------
    >>> protein = get_structure("protein.pdb")
    >>> coords = get_attr(protein, "coord")
    >>> new_coords = coords + 10.0  # 10Å移動
    >>> set_attr(protein, "coord", new_coords)
    """

    # TODO: lstの長さと原子数の整合性チェック
    # 異なる場合は値を設定せずにエラーを出すべき

    lst_idx = 0
    for atom in model.get_atoms():
        if sele is None or sele(atom):
            if attr == "coord":
                atom.set_coord(lst[lst_idx])
            else:
                raise NotImplementedError(f"set_attr(attr={attr}) は実装されていません")
            lst_idx += 1


def save(structs, path) -> None:
    """
    構造をPDBファイルに保存します。
    
    1つまたは複数の構造をPDBファイルに保存します。
    複数の構造が与えられた場合は、MODELレコードを使用して保存されます。

    Parameters
    ----------
    structs : Bio.PDB.Structure.Structure or list of Bio.PDB.Structure.Structure
        保存する構造（単一または複数）
    path : str
        出力PDBファイルのパス
        
    Examples
    --------
    >>> # 単一構造の保存
    >>> save(protein, "output.pdb")
    >>> 
    >>> # 複数構造の保存
    >>> save([protein1, protein2], "multi_model.pdb")
    """
    path = expandpath(path)

    if not isinstance(structs, Iterable):
        structs = [structs]

    mod_structs = []
    for struct in structs:
        io = BioPDB.PDBIO()
        io.set_structure(struct)
        with tempfile.NamedTemporaryFile(suffix=".pdb") as fp:
            io.save(fp.name)
            mod_structs.append(get_structure(fp.name)[0])

    out_structure = Structure("")
    for struct in mod_structs:
        struct.id = len(out_structure)
        struct.serial_num = struct.id + 1
        out_structure.add(struct)

    io = BioPDB.PDBIO()
    io.set_structure(out_structure)
    io.save(path)


# PDB関数の名前空間オブジェクトを作成
class PDB:
    """
    PDBユーティリティ関数の名前空間クラス
    
    PDBファイルの読み書きと属性の取得・設定のための
    ユーティリティ関数を提供します。
    
    Examples
    --------
    >>> from inverse_msmd import PDB
    >>> 
    >>> # PDBファイルの読み込み
    >>> protein = PDB.get_structure("protein.pdb")
    >>> 
    >>> # 座標の取得
    >>> coords = PDB.get_attr(protein, "coord")
    >>> 
    >>> # 座標の設定
    >>> PDB.set_attr(protein, "coord", new_coords)
    >>> 
    >>> # 保存
    >>> PDB.save(protein, "output.pdb")
    """
    get_structure = staticmethod(get_structure)
    get_attr = staticmethod(get_attr)
    set_attr = staticmethod(set_attr)
    save = staticmethod(save)
    get_resname = staticmethod(get_resname)
    get_resi = staticmethod(get_resi)
    get_atom_attr = staticmethod(get_atom_attr)


__all__ = [
    "SuperImposer",
    "PDB",
    "get_structure",
    "get_attr",
    "set_attr",
    "save",
]