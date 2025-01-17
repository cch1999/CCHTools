from typing import Optional

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def get_coords(mol: Chem.Mol) -> np.ndarray:
    """Get coordinates of all atoms in an RDKit molecule.

    Args:
        mol: RDKit molecule object

    Returns:
        Array of shape (n_atoms, 3) containing XYZ coordinates
    """
    conf = mol.GetConformer()
    coords = conf.GetPositions()
    return np.array(coords)


def get_center_of_mass(mol: Chem.Mol) -> np.ndarray:
    """Calculate the center of mass of an RDKit molecule.

    Args:
        mol: RDKit molecule object

    Returns:
        Array of shape (3,) containing XYZ coordinates of center of mass
    """
    coords = get_coords(mol)
    center_of_mass = np.mean(coords, axis=0)
    return center_of_mass


def generate_conformer(mol: Chem.Mol, n_conformers: int = 1) -> Optional[Chem.Mol]:
    """Generate multiple conformers for an RDKit molecule.

    Args:
        mol: RDKit molecule object
        n_conformers: Number of conformers to generate (default: 1)

    Returns:
        List of RDKit molecule objects, each with a different conformer
    """
    try:
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.useSmallRingTorsions = True
        AllChem.EmbedMultipleConfs(mol, numConfs=n_conformers, params=params)
        AllChem.MMFFOptimizeMoleculeConfs(mol)
        mol = Chem.RemoveHs(mol)
    except ValueError:
        return None

    return mol
