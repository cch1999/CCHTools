import os

import datamol as dm
import pytest

from cchtools.constants import EXAMPLES_PDB, EXAMPLES_SDF
from cchtools.pymol_interface import PyMOLInterface


@pytest.fixture
def pymol_interface():
    return PyMOLInterface()


def test_init(pymol_interface):
    assert pymol_interface.background_color == "white"
    assert pymol_interface.cmd.get_viewport() == (400, 400)


def test_load_pdb(pymol_interface):
    pymol_interface.load_pdb(EXAMPLES_PDB)
    assert "protein" in pymol_interface.cmd.get_names()


def test_load_molecule(pymol_interface):
    ligand = dm.read_sdf(EXAMPLES_SDF)[0]
    pymol_interface.load_molecule(ligand)
    assert "ligand" in pymol_interface.cmd.get_names()


def test_set_protein_color(pymol_interface):
    pymol_interface.load_pdb(EXAMPLES_PDB)
    pymol_interface.set_protein_color("red")


def test_get_interacting_residues(pymol_interface):
    pymol_interface.load_pdb(EXAMPLES_PDB)
    ligand = dm.read_sdf(EXAMPLES_SDF)[0]
    pymol_interface.load_molecule(ligand)
    residues = pymol_interface.get_interacting_residues(distance=3.0)
    assert isinstance(residues, list)
    assert all(isinstance(res, int) for res in residues)


def test_save_image(pymol_interface, tmp_path):
    pymol_interface.load_pdb(EXAMPLES_PDB)
    image_path = tmp_path / "test_image.png"
    pymol_interface.save_image(str(image_path))
    assert os.path.exists(image_path)


def test_display_spheres(pymol_interface):
    positions = [(0, 0, 0), (1, 1, 1)]
    pymol_interface.display_spheres(positions)
    assert "sphere0" in pymol_interface.cmd.get_names()
    assert "sphere1" in pymol_interface.cmd.get_names()


def test_show_residue_sidechains(pymol_interface):
    pymol_interface.load_pdb(EXAMPLES_PDB)
    residue_ids = [1, 2, 3]
    pymol_interface.show_residue_sidechains(residue_ids)
    # This test is limited as we can't easily check the visual output
    # We're just ensuring the function runs without errors


def test_remove_all_hydrogens(pymol_interface):
    pymol_interface.load_pdb(EXAMPLES_PDB)
    initial_atom_count = pymol_interface.cmd.count_atoms()
    pymol_interface.remove_all_hydrogens()
    final_atom_count = pymol_interface.cmd.count_atoms()
    assert final_atom_count < initial_atom_count
