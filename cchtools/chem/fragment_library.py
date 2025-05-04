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
import pickle
from dataclasses import dataclass, field
from typing import Iterable, List, Optional, Sequence, Tuple, Union, Dict, Any

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


def standardise(mol: Chem.Mol) -> Optional[Chem.Mol]:
    """
    Standardise a molecule by:
      1. Cleanup (disconnect metals, remove isotopes, etc.)
      2. Stripping salts (using RDKit's default salt list)
      3. Neutralising to a charge-balanced form
      4. Canonicalising to an isomeric SMILES
      
    Returns:
        The standardized molecule or None if standardization failed
    """
    if mol is None:
        logger.warning("Cannot standardize None molecule")
        return None
        
    try:
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
    except Exception as e:
        logger.error(f"Error standardizing molecule: {e}")
        return None


def _canonical_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, isomericSmiles=True)


@dataclass
class FragmentLibrary:
    """Container around a *set* of molecular fragments."""

    _mols: List[Chem.Mol] = field(default_factory=list, repr=False)
    _smiles: List[str] = field(default_factory=list, repr=False)
    _fps: List[DataStructs.cDataStructs.ExplicitBitVect] = field(default_factory=list, repr=False)
    _metadata: Dict[str, Any] = field(default_factory=dict, repr=False)

    # ---------------------- construction ----------------------

    @classmethod
    def from_iterable(cls, frags: Iterable[MolLike]) -> "FragmentLibrary":
        """
        Create a fragment library from an iterable of SMILES strings or RDKit Mol objects.
        
        Args:
            frags: An iterable of SMILES strings or RDKit Mol objects
            
        Returns:
            A new FragmentLibrary instance
        """
        lib = cls()
        for f in frags:
            lib.add(f)
        return lib

    @classmethod
    def from_file(
        cls, path: Union[str, pathlib.Path], smiles_column: str | None = None
    ) -> "FragmentLibrary":
        """
        Create a fragment library from an SDF or CSV/TSV file.
        
        Args:
            path: Path to the input file
            smiles_column: Name of the column containing SMILES strings (for CSV/TSV files)
            
        Returns:
            A new FragmentLibrary instance
            
        Raises:
            ValueError: If smiles_column is not provided for a CSV/TSV file
        """
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
        
    @classmethod
    def from_molecules(cls, mols: Iterable[Chem.Mol], extract_props: bool = False) -> "FragmentLibrary":
        """
        Create a fragment library from a list of RDKit molecules.
        
        Args:
            mols: An iterable of RDKit Mol objects
            extract_props: If True, extract properties from molecules and store in metadata
            
        Returns:
            A new FragmentLibrary instance
        """
        lib = cls()
        for mol in mols:
            if mol is None:
                continue
            
            # Extract properties if requested
            if extract_props:
                props = {}
                for prop_name in mol.GetPropNames():
                    props[prop_name] = mol.GetProp(prop_name)
                
                # Store properties indexed by SMILES
                smi = _canonical_smiles(mol)
                if "properties" not in lib._metadata:
                    lib._metadata["properties"] = {}
                lib._metadata["properties"][smi] = props
                
            lib.add(mol)
            
        return lib

    def _from_file(self, path: Union[str, pathlib.Path], smiles_column: str | None = None) -> None:
        """
        Load fragments into this library in-place from an SDF/CSV file.
        
        Args:
            path: Path to the input file
            smiles_column: Name of the column containing SMILES strings (for CSV/TSV files)
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
        
    def __add__(self, other: "FragmentLibrary") -> "FragmentLibrary":
        """Combine two fragment libraries using the + operator."""
        return self.union(other)

    # ---------------------- properties ----------------------

    @property
    def mols(self) -> List[Chem.Mol]:
        """Get the list of RDKit molecules in this library."""
        return self._mols

    @property
    def smiles(self) -> List[str]:
        """Get the list of canonical SMILES strings in this library."""
        return self._smiles

    @property
    def fps(self) -> List[DataStructs.cDataStructs.ExplicitBitVect]:
        """Get the list of molecular fingerprints in this library."""
        return self._fps
        
    @property
    def metadata(self) -> Dict[str, Any]:
        """Get the metadata dictionary for this library."""
        return self._metadata

    # ---------------------- mutators ----------------------

    def add(self, mol: MolLike) -> None:
        """
        Add a molecule to the fragment library.
        
        Args:
            mol: A SMILES string or RDKit Mol object
        """
        # parse & standardise the fragment
        try:
            raw = _to_mol(mol)
            m = standardise(raw)
            
            if m is None or m.GetNumAtoms() == 0:
                logger.warning(f"Could not parse or standardize fragment: {mol}")
                return
                
            smi = _canonical_smiles(m)
            if smi in self._smiles:
                logger.warning(f"Duplicate fragment: {smi}")
                return
                
            self._mols.append(m)
            self._smiles.append(smi)
            self._fps.append(AllChem.GetMorganFingerprintAsBitVect(m, radius=2, nBits=2048))
        except Exception as e:
            logger.error(f"Error adding molecule to library: {e}")

    # ---------------------- queries ----------------------

    def substructure_matches(
        self,
        mol: MolLike,
        return_atom_maps: bool = False,
    ) -> List[Tuple[str, Chem.Mol, Chem.Mol, Sequence[int]]]:
        """
        Return a list of (fragment_smiles, atom_map) pairs found in `mol`.

        If `return_atom_maps` is False, the atom_map tuple is empty.

        Args:
            mol: Query molecule as SMILES or RDKit Mol
            return_atom_maps: Whether to return atom maps for matches
            
        Returns:
            List[Tuple[str, Chem.Mol, Chem.Mol, Sequence[int]]]:
                A list of tuples containing:
                - fragment_smiles: The SMILES representation of the fragment
                - remaining_mol: The remaining molecule after removing the fragment
                - fragment_mol: The fragment molecule
                - atom_map: A tuple of atom indices in the fragment that match the atoms in the query molecule
        """
        m = _to_mol(mol)
        hits: List[Tuple[str, Chem.Mol, Chem.Mol, Sequence[int]]] = []
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
        """
        Boolean convenience wrapper around `substructure_matches`.
        
        Args:
            mol: Query molecule as SMILES or RDKit Mol
            
        Returns:
            True if any fragment in the library is a substructure of the query molecule
        """
        return bool(self.substructure_matches(mol))

    def get_similar(self, query_frag: MolLike, threshold: float = 0.5) -> List[Chem.Mol]:
        """
        Return a list of fragments from this library that are similar to `query_frag`.
        
        Args:
            query_frag: Query fragment as SMILES or RDKit Mol
            threshold: Minimum Tanimoto similarity (0-1) required for a match
            
        Returns:
            List of RDKit Mol objects that are similar to the query fragment
        """
        query = _to_mol(query_frag)
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query, radius=2, nBits=2048)
        hits = []
        for frag_smi, frag_mol, frag_fp in zip(self._smiles, self._mols, self._fps):
            if DataStructs.TanimotoSimilarity(frag_fp, query_fp) >= threshold:
                hits.append(frag_mol)
        return hits
        
    def get_similar_with_scores(self, query_frag: MolLike, threshold: float = 0.5) -> List[Tuple[Chem.Mol, float]]:
        """
        Return a list of fragments with similarity scores from this library.
        
        Args:
            query_frag: Query fragment as SMILES or RDKit Mol
            threshold: Minimum Tanimoto similarity (0-1) required for a match
            
        Returns:
            List of tuples containing (RDKit Mol, similarity score)
        """
        query = _to_mol(query_frag)
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query, radius=2, nBits=2048)
        hits = []
        for frag_smi, frag_mol, frag_fp in zip(self._smiles, self._mols, self._fps):
            score = DataStructs.TanimotoSimilarity(frag_fp, query_fp)
            if score >= threshold:
                hits.append((frag_mol, score))
        # Sort by descending similarity
        return sorted(hits, key=lambda x: x[1], reverse=True)

    # ---------------------- vectoriser ----------------------

    def featurise(self, mol: MolLike, mode: str = "binary") -> np.ndarray:
        """
        Convert a molecule to a presence/​count vector of fragments.

        Args:
            mol: Query molecule as SMILES or RDKit Mol
            mode: Either "binary" (0/1) or "count" (number of matches)
            
        Returns:
            NumPy array with one element per fragment in the library
            
        Raises:
            ValueError: If mode is not "binary" or "count"
        """
        if mode not in ["binary", "count"]:
            raise ValueError(f'Mode must be "binary" or "count", got {mode!r}')
            
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
        """
        Create a new library containing all fragments from both libraries.
        
        Args:
            other: Another FragmentLibrary
            
        Returns:
            A new FragmentLibrary with the union of fragments
        """
        return self._operate(other, set.union)

    def intersection(self, other: "FragmentLibrary") -> "FragmentLibrary":
        """
        Create a new library containing only fragments present in both libraries.
        
        Args:
            other: Another FragmentLibrary
            
        Returns:
            A new FragmentLibrary with the intersection of fragments
        """
        return self._operate(other, set.intersection)

    def difference(self, other: "FragmentLibrary") -> "FragmentLibrary":
        """
        Create a new library containing fragments in this library but not in the other.
        
        Args:
            other: Another FragmentLibrary
            
        Returns:
            A new FragmentLibrary with the difference of fragments
        """
        return self._operate(other, set.difference)

    # ---------------------- I/O ----------------------

    def to_smiles(self, path: Union[str, pathlib.Path]) -> None:
        """
        Save the fragment library as a simple SMILES file (one per line).
        
        Args:
            path: Output file path
        """
        pathlib.Path(path).write_text("\n".join(self._smiles))

    def to_json(self, path: Union[str, pathlib.Path]) -> None:
        """
        Save the fragment library as a JSON file.
        
        Args:
            path: Output file path
        """
        meta = {
            "n_fragments": len(self),
            "canonical_smiles": self._smiles,
        }
        
        # Add any custom metadata
        if self._metadata:
            meta["metadata"] = self._metadata
            
        pathlib.Path(path).write_text(json.dumps(meta, indent=2))

    def to_sdf(self, path: Union[str, pathlib.Path]) -> None:
        """
        Save the fragment library as an SDF file.
        
        Args:
            path: Output file path
        """
        writer = Chem.SDWriter(str(path))
        for mol in self._mols:
            writer.write(mol)
        writer.close()
        
    def to_pickle(self, path: Union[str, pathlib.Path]) -> None:
        """
        Save the fragment library to a pickle file for fast loading.
        
        Args:
            path: Output file path
        """
        with open(path, 'wb') as f:
            pickle.dump(self, f)
    
    @classmethod
    def from_pickle(cls, path: Union[str, pathlib.Path]) -> "FragmentLibrary":
        """
        Load a fragment library from a pickle file.
        
        Args:
            path: Path to the pickle file
            
        Returns:
            A FragmentLibrary instance
            
        Raises:
            ValueError: If the loaded object is not a FragmentLibrary
        """
        with open(path, 'rb') as f:
            lib = pickle.load(f)
            
        if not isinstance(lib, FragmentLibrary):
            raise ValueError(f"Loaded object is not a FragmentLibrary: {type(lib)}")
            
        return lib

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
