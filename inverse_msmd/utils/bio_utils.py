# coding: utf-8

"""
BioPython utilities for structure manipulation and superimposition.

Combined from Bio/PDB.py and Bio/sklearn_interface.py
version: 1.0.0
Authors: Keisuke Yanagisawa
"""
import collections
import gzip
import os
from typing import Any, Callable, List, Literal, Optional, Union
import warnings
import numpy as np
from Bio import PDB
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
    構造重ね合わせを行うBioPythonのクラスを
    scikit-learnのインターフェースでwrapしたクラス。
    """

    rot_: npt.NDArray[np.float_]
    tran_: npt.NDArray[np.float_]

    def __init__(self):
        pass

    def _reset(self):
        if hasattr(self, "rot_"):
            del self.rot_
            del self.tran_

    def _superimpose(self, coords: npt.ArrayLike, reference_coords: npt.ArrayLike) -> None:
        sup = SVDSuperimposer()
        sup.set(reference_coords, coords)
        sup.run()
        self.rot_, self.tran_ = sup.get_rotran()  # type: ignore

    def fit(self, coords: npt.ArrayLike, reference_coords: npt.ArrayLike) -> "SuperImposer":
        """
        与えられた2つの点群をなるべく重ねるような並行・回転移動を算出します。

        与えられた2つの点群はそれぞれ対応関係があることを仮定します。
        すなわち、それぞれの0番目の要素同士がなるべく重なるように、
        1番目の要素同士がなるべく重なるように…と重ね合わせを行います。

        Parameters
        ----------
        coords : list
            重ね合わせのために移動させる点群
        reference_coords : list
            重ね合わせ先の点群

        Returns
        -------
        SuperImposer
            fit済みのオブジェクト
        """
        self._reset()
        self._superimpose(coords, reference_coords)
        return self

    def transform(self, coords: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """
        fit()で計算された並進・回転に基づいて
        与えられた点群を移動させます。

        Parameters
        ----------
        coords : list
            移動させる点群
        """
        check_is_fitted(self)
        coords = np.array(coords)
        return np.dot(coords, self.rot_) + self.tran_

    def inverse_transform(self, coords: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """
        逆方向の移動を行います。

        Parameters
        ----------
        coords : list
            transform()した後の点群

        Returns
        -------
        np.array
            transform()する前の点群座標
        """
        coords = np.array(coords)
        check_is_fitted(self)
        return np.dot(coords - self.tran_, np.linalg.inv(self.rot_))


# ============================================================================
# PDB utilities (from PDB.py)
# ============================================================================

def get_structure(filepath: str, structname="") -> Structure:
    """
    Read PDB file.

    Parameters
    ----------
    filepath : str
        filepath to a PDB file

    Returns
    -------
    Bio.PDB.Structure
        Structure object of the PDB file
    """
    filepath = expandpath(filepath)
    if filepath.endswith(".gz"):
        fileobj = gzip.open(filepath, "rt")
    else:
        fileobj = open(filepath)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBExceptions.PDBConstructionWarning)
        return PDB.PDBParser(QUIET=True).get_structure(structname, fileobj)


def get_atom_attr(atom: Atom,
                  attr: Literal["resid", "resname", "coord", "element", "fullname"]
                  ) -> Union[int, str, npt.NDArray[np.float_], tuple]:
    """
    Get attribute from Bio.PDB.Atom object.
    {"resid", "resname", "coord", "element", "fullname"}
    are only acceptable as ``attr`` so far.
    Other attributes raises NotImplementedError.

    Parameters
    ----------
    atom : Atom
        An atom object.
    attr : str
        An attribute name which will be obtained.

    Returns
    -------
    int or str or np.array or tuple
        An attribute of the atom.

    Raises
    ------
    NotImplementedError
        If the ``attr`` is not "resid", "resname", coord", "element", nor "fullname".
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
        raise NotImplementedError(f"Attribute {attr} is not supported yet.")


def get_attr(model: Union[Structure, Model],
             attr: Literal["resid", "resname", "coord", "element", "fullname"],
             sele: Optional[Callable[[Atom], bool]] = None
             ) -> npt.NDArray[Any]:
    """
    Get attribute from Bio.PDB.Model object.
    {"resid", "resname", "coord", "element", "fullname"} 
    are only acceptable as ``attr`` so far.
    Other attributes raises NotImplementedError.

    Parameters
    ----------
    model : Model
        A model object.
    attr : str
        An attribute name which will be obtained.
    sele : function, optional
        Atom selector function. all atoms will be selected if ``sele`` is not
        provided.

    Returns
    -------
    np.array
        An array of attributes of all atoms selected by ``sele`` function.

    Raises
    ------
    NotImplementedError
        If the ``attr`` is not "resid", "resname", coord", "element", nor "fullname".
    """

    data = []
    for atom in model.get_atoms():
        if sele is None or sele(atom):
            data.append(get_atom_attr(atom, attr))
    return np.array(data)


def get_resname(atom: Atom) -> str:
    """
    Get residue name from Bio.PDB.Atom.

    Parameters
    ----------
    atom : Bio.PDB.Atom
        An atom object.

    Returns
    -------
    str
        A residue name of the atom of interest.
    """
    return atom.get_parent().get_resname()  # type: ignore


def get_resi(atom: Atom) -> int:
    """
    Get residue sequence number (residue ID) from Bio.PDB.Atom.

    Parameters
    ----------
    atom : Bio.PDB.Atom
        An atom object.

    Returns
    -------
    int
        A residue ID of the atom of interest.
    """
    return atom.get_full_id()[3][1]  # type: ignore


def set_attr(model: Model, attr: str, lst: npt.NDArray, sele=None) -> None:
    """
    Set attribute to Bio.PDB.Model object.
    attr == "coord" is only acceptable so far.
    Other attributes raises NotImplementedError.

    Parameters
    ----------
    model : Bio.PDB.Model
    attr : str
    lst : array_like
    sele : function, optional

    Raises
    ------
    NotImplementedError
        If the ``attr`` is not "coord".
    """

    # TODO check the length of lst and the number of atoms.
    # if they are different, set_attr() must not assign new values.

    lst_idx = 0
    for atom in model.get_atoms():
        if sele is None or sele(atom):
            if attr == "coord":
                atom.set_coord(lst[lst_idx])
            else:
                raise NotImplementedError(f"set_attr(attr={attr}) is not implemented")
            lst_idx += 1


def save(structs, path) -> None:
    """
    Save structure(s) to PDB file.

    Parameters
    ----------
    structs : Bio.PDB.Struct or list of Bio.PDB.Struct
    path : str
    """
    path = expandpath(path)

    if not isinstance(structs, Iterable):
        structs = [structs]

    mod_structs = []
    for struct in structs:
        io = PDB.PDBIO()
        io.set_structure(struct)
        with tempfile.NamedTemporaryFile(suffix=".pdb") as fp:
            io.save(fp.name)
            mod_structs.append(get_structure(fp.name)[0])

    out_structure = Structure("")
    for struct in mod_structs:
        struct.id = len(out_structure)
        struct.serial_num = struct.id + 1
        out_structure.add(struct)

    io = PDB.PDBIO()
    io.set_structure(out_structure)
    io.save(path)


# Create a namespace-like object for PDB functions
class PDB:
    """PDB utility functions namespace."""
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