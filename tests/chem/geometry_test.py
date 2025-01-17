import pytest
import numpy as np
from rdkit import Chem

import cchtools.chem.geometry as geometry


@pytest.fixture(scope="module")
def test_smiles():
    """Fixture providing test SMILES strings."""
    return ["C", "CC(=O)O", "c1ccccc1", "C1CCCCC1",
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C"]


@pytest.fixture(scope="module")
def test_mols(test_smiles):
    """Fixture providing test molecules with conformers."""
    mols = []
    for smile in test_smiles:
        mol = Chem.MolFromSmiles(smile)
        mol = geometry.generate_conformer(mol)
        mols.append(mol)
    return mols


def test_get_coords(test_mols):
    """Test getting coordinates from molecule."""
    for mol in test_mols:
        coords = geometry.get_coords(mol)
        assert isinstance(coords, np.ndarray)
        assert coords.shape[1] == 3  # XYZ coordinates
        assert coords.shape[0] == mol.GetNumAtoms()


def test_get_center_of_mass(test_mols):
    """Test calculating center of mass."""
    for mol in test_mols:
        com = geometry.get_center_of_mass(mol)
        assert isinstance(com, np.ndarray)
        assert com.shape == (3,)  # XYZ coordinates

        # Center of mass should be within the bounds of atomic coordinates
        coords = geometry.get_coords(mol)
        for i in range(3):
            assert coords[:, i].min() <= com[i] <= coords[:, i].max()


def test_generate_conformer():
    """Test conformer generation."""
    # Test with valid molecule
    mol = Chem.MolFromSmiles("C")
    conf_mol = geometry.generate_conformer(mol, n_conformers=1)
    assert isinstance(conf_mol, Chem.Mol)
    assert conf_mol.GetNumConformers() == 1

    # Test with multiple conformers
    conf_mol = geometry.generate_conformer(mol, n_conformers=3)
    assert conf_mol.GetNumConformers() == 3
