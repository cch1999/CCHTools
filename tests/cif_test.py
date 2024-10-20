import os
import tempfile

import pytest
from rdkit import Chem

from cchtools.rdkit_utils import fix_bond_orders_with_ccd
from cchtools.utils.cif import (
    cif_to_rdkit,
    download_ideal_ccd_structure,
    fetch_ideal_ccd_structure,
)


@pytest.fixture(scope="module")
def temp_dir():
    with tempfile.TemporaryDirectory() as tmpdirname:
        yield tmpdirname


@pytest.fixture(params=["HEM", "ATP", "NAD", "GTP", "FAD", "355"])
def test_ccd_code(request):
    return request.param


def test_download_ideal_ccd_structure(temp_dir, test_ccd_code):
    cif_path = download_ideal_ccd_structure(test_ccd_code, temp_dir)
    assert cif_path is not None
    assert os.path.exists(cif_path)
    assert cif_path.endswith(".cif")


def test_cif_to_rdkit(temp_dir, test_ccd_code):
    cif_path = download_ideal_ccd_structure(test_ccd_code, temp_dir)
    mol = cif_to_rdkit(cif_path)
    assert mol is not None
    assert isinstance(mol, Chem.Mol)


def test_fetch_ideal_ccd_structure(temp_dir, test_ccd_code):
    mol = fetch_ideal_ccd_structure(test_ccd_code, temp_dir)
    assert mol is not None
    assert isinstance(mol, Chem.Mol)


def test_fetch_ideal_ccd_structure_no_save(temp_dir, test_ccd_code):
    mol = fetch_ideal_ccd_structure(test_ccd_code, temp_dir, save_structures=False)
    assert mol is not None
    assert isinstance(mol, Chem.Mol)
    cif_path = os.path.join(temp_dir, f"{test_ccd_code}.cif")
    assert not os.path.exists(cif_path)


def test_invalid_ccd_code(temp_dir):
    invalid_code = "INVALID"
    cif_path = download_ideal_ccd_structure(invalid_code, temp_dir)
    assert cif_path is None
    mol = fetch_ideal_ccd_structure(invalid_code, temp_dir)
    assert mol is None


def test_fix_bond_orders_with_ccd(temp_dir, test_ccd_code):
    mol = fetch_ideal_ccd_structure(test_ccd_code, temp_dir)

    # shuffle atom ids
    mol = Chem.RenumberAtoms(mol, newOrder=list(range(mol.GetNumAtoms()))[::-1])

    # Set all bonds to single
    for bond in mol.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)

    fixed_mol = fix_bond_orders_with_ccd(mol, test_ccd_code, temp_dir)
    assert fixed_mol is not None
    assert isinstance(fixed_mol, Chem.Mol)
