"""
fragment_library.py

A lightweight, RDKit-backed container for small-molecule fragments.

Core features
-------------
*  Exact-member test using canonical SMILES   (`frag in lib`)
*  Substructure search with optional match maps
*  One-hot / count vectorisation for ML
*  Basic set algebra (union, intersection, difference)
*  Quick CSV|SDF I/O helpers
"""

from __future__ import annotations

import json
import logging
import pathlib
from dataclasses import dataclass, field
from typing import Iterable, List, Sequence, Tuple, Union

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.MolStandardize import rdMolStandardize as rdms
from rdkit.Chem.rdmolops import ReplaceSidechains
from rdkit.Chem.SaltRemover import SaltRemover

logger = logging.getLogger(__name__)

MolLike = Union[str, Chem.Mol]  # SMILES | Mol

MINIFRAGS_PATH = (
    pathlib.Path(__file__).parent.parent / "data" / "Enamine_MiniFrag_Library_80cmpds_20250123.sdf"
)
POISED_FRAGS_PATH = (
    pathlib.Path(__file__).parent.parent
    / "data"
    / "Enamine_DSi-Poised_Library_plated_860cmpds_20250309.sdf"
)


def _to_mol(obj: MolLike) -> Chem.Mol:
    if isinstance(obj, Chem.Mol):
        return obj
    m = Chem.MolFromSmiles(obj)
    if m is None:
        raise ValueError(f"Could not parse SMILES: {obj!r}")
    return m


def standardise(mol: Chem.Mol) -> Chem.Mol:
    """
    Standardise a molecule by:
      1. Cleanup (disconnect metals, remove isotopes, etc.)
      2. Stripping salts (using RDKit’s default salt list)
      3. Neutralising to a charge-balanced form
      4. Canonicalising to an isomeric SMILES
    """
    # 1. Cleanup
    clean_mol = rdms.Cleanup(mol)
    # 2. Strip salts
    remover = SaltRemover()
    clean_mol = remover.StripMol(clean_mol, dontRemoveEverything=False)
    # 3. Neutralise
    uncharger = rdms.Uncharger()
    clean_mol = uncharger.uncharge(clean_mol)
    # 4. SMILES-based canonicalisation & reparse
    can = Chem.MolToSmiles(clean_mol, isomericSmiles=True)

    # 5. Check if molecule has duplicate atoms
    if "." in can:
        split = can.split(".")
        if split[0] == split[1]:
            logger.warning(f"Duplicate fragment: {can}. Returning first half.")
            can = split[0]
        else:
            logger.warning(f"Duplicate fragment: {can}. Not returning anything.")
            return None

    return Chem.MolFromSmiles(can)


def _canonical_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, isomericSmiles=True)


@dataclass
class FragmentLibrary:
    """Container around a *set* of molecular fragments."""

    _mols: List[Chem.Mol] = field(default_factory=list, repr=False)
    _smiles: List[str] = field(default_factory=list, repr=False)
    _fps: List[DataStructs.cDataStructs.ExplicitBitVect] = field(default_factory=list, repr=False)

    # ---------------------- construction ----------------------

    @classmethod
    def from_iterable(cls, frags: Iterable[MolLike]) -> "FragmentLibrary":
        lib = cls()
        for f in frags:
            lib.add(f)
        return lib

    @classmethod
    def from_file(
        cls, path: Union[str, pathlib.Path], smiles_column: str | None = None
    ) -> "FragmentLibrary":
        path = pathlib.Path(path)
        if path.suffix.lower() == ".sdf":
            suppl = Chem.SDMolSupplier(str(path), removeHs=True)
            return cls.from_iterable([m for m in suppl if m])
        # assume CSV/TSV
        if smiles_column is None:
            raise ValueError("CSV/TSV input requires 'smiles_column' kwarg.")
        lines = path.read_text().splitlines()
        header = lines[0].split(",")
        idx = header.index(smiles_column)
        return cls.from_iterable(line.split(",")[idx] for line in lines[1:])

    def _from_file(self, path: Union[str, pathlib.Path], smiles_column: str | None = None) -> None:
        """
        Load fragments into this library in-place from an SDF/CSV file.
        """
        lib = FragmentLibrary.from_file(path, smiles_column)
        self._mols.clear()
        self._smiles.clear()
        self._fps.clear()
        self._mols.extend(lib._mols)
        self._smiles.extend(lib._smiles)
        self._fps.extend(lib._fps)

    # ---------------------- dunder bits ----------------------

    def __len__(self) -> int:
        return len(self._mols)

    def __iter__(self):
        return iter(self._mols)

    def __contains__(self, mol: MolLike) -> bool:
        smi = _canonical_smiles(_to_mol(mol))
        return smi in self._smiles

    # ---------------------- properties ----------------------

    @property
    def mols(self) -> List[Chem.Mol]:
        return self._mols

    @property
    def smiles(self) -> List[str]:
        return self._smiles

    @property
    def fps(self) -> List[DataStructs.cDataStructs.ExplicitBitVect]:
        return self._fps

    # ---------------------- mutators ----------------------

    def add(self, mol: MolLike) -> None:
        # parse & standardise the fragment
        raw = _to_mol(mol)
        m = standardise(raw)
        smi = _canonical_smiles(m)
        if smi in self._smiles:
            logger.warning(f"Duplicate fragment: {smi}")
            return
        if m is None or m.GetNumAtoms() == 0:
            logger.warning(f"Could not parse fragment: {mol}")
            return
        self._mols.append(m)
        self._smiles.append(smi)
        self._fps.append(AllChem.GetMorganFingerprintAsBitVect(m, radius=2, nBits=2048))

    # ---------------------- queries ----------------------

    def substructure_matches(
        self,
        mol: MolLike,
        return_atom_maps: bool = False,
    ) -> List[Tuple[str, Sequence[int]]]:
        """
        Return a list of (fragment_smiles, atom_map) pairs found in `mol`.

        If `return_atom_maps` is False, the atom_map tuple is empty.
        """
        m = _to_mol(mol)
        hits: List[Tuple[str, Chem.Mol, Sequence[int]]] = []
        for frag_smi, frag_mol in zip(self._smiles, self._mols):
            if m.HasSubstructMatch(frag_mol):
                amap = m.GetSubstructMatch(frag_mol) if return_atom_maps else ()
                # remove the matching fragment from the query, keeping its original 3D coordinates
                remaining = Chem.DeleteSubstructs(m, frag_mol)
                # keep the matching fragment
                # extract the fragment substructure using ReplaceSidechains
                frag = ReplaceSidechains(m, frag_mol)

                hits.append((frag_smi, remaining, frag, amap))
        return hits

    def contains_substructure(self, mol: MolLike) -> bool:
        """Boolean convenience wrapper around `substructure_matches`."""
        return bool(self.substructure_matches(mol))

    def get_similar(self, query_frag: MolLike, threshold: float = 0.5) -> List[Chem.Mol]:
        """
        Return a list of fragments from this library that are similar to `query_frag`.
        """
        query = _to_mol(query_frag)
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query, radius=2, nBits=2048)
        hits = []
        for frag_smi, frag_mol, frag_fp in zip(self._smiles, self._mols, self._fps):
            if DataStructs.TanimotoSimilarity(frag_fp, query_fp) >= threshold:
                hits.append(frag_mol)
        return hits

    # ---------------------- vectoriser ----------------------

    def featurise(self, mol: MolLike, mode: str = "binary") -> np.ndarray:
        """
        Convert a molecule to a presence/​count vector of fragments.

        `mode` ∈ {"binary", "count"}.
        """
        m = _to_mol(mol)
        v = np.zeros(len(self), dtype=int)
        for i, fmol in enumerate(self._mols):
            n = len(m.GetSubstructMatches(fmol))
            if n:
                v[i] = 1 if mode == "binary" else n
        return v

    # ---------------------- set algebra ----------------------

    def _operate(self, other: "FragmentLibrary", op) -> "FragmentLibrary":
        assert isinstance(other, FragmentLibrary)
        keep = [s for s in op(set(self._smiles), set(other._smiles))]
        return FragmentLibrary.from_iterable(keep)

    def union(self, other: "FragmentLibrary") -> "FragmentLibrary":
        return self._operate(other, set.union)

    def intersection(self, other: "FragmentLibrary") -> "FragmentLibrary":
        return self._operate(other, set.intersection)

    def difference(self, other: "FragmentLibrary") -> "FragmentLibrary":
        return self._operate(other, set.difference)

    # ---------------------- I/O ----------------------

    def to_smiles(self, path: Union[str, pathlib.Path]) -> None:
        pathlib.Path(path).write_text("\n".join(self._smiles))

    def to_json(self, path: Union[str, pathlib.Path]) -> None:
        meta = {
            "n_fragments": len(self),
            "canonical_smiles": self._smiles,
        }
        pathlib.Path(path).write_text(json.dumps(meta, indent=2))

    # ---------------------- pretty ----------------------

    def __repr__(self) -> str:  # noqa: D401
        return f"<FragmentLibrary n={len(self)}>"


class MiniFragsLib(FragmentLibrary):
    def __init__(self, path: Union[str, pathlib.Path] = MINIFRAGS_PATH):
        super().__init__()
        # correctly load the mini‐frag SDF into this instance
        self._from_file(path)

    def __repr__(self) -> str:
        return f"<MiniFragsLibrary n={len(self)}>"


class PoisedFragsLib(FragmentLibrary):
    def __init__(self, path: Union[str, pathlib.Path] = POISED_FRAGS_PATH):
        super().__init__()
        # correctly load the poised‐frag SDF into this instance
        self._from_file(path)

    def __repr__(self) -> str:
        return f"<PoisedFragsLibrary n={len(self)}>"


if __name__ == "__main__":
    # Two tiny benzamide fragments as a demo
    frags = ["c1ccccc1C(=O)N", "c1ccncc1"]
    lib = FragmentLibrary.from_iterable(frags)

    query = "c1ccc(cc1)C(=O)NC"  # acetanilide
    print(query in lib)  # False (exact)
    print(lib.contains_substructure(query))  # True
    print(lib.substructure_matches(query))  # [('c1ccccc1C(=O)N', tuple())]

    print(lib.featurise(query))  # [1 0]

    print(lib)
