import pytest
import json
from rdkit import Chem
from cchtools.chem.fragment_library import FragmentLibrary, _to_mol, MiniFragsLib, PoisedFragsLib

def test_to_mol_accepts_smiles_and_mol_and_errors_on_invalid():
    # valid SMILES string
    m1 = _to_mol("CCO")
    assert isinstance(m1, Chem.Mol)
    # valid Mol instance
    m2 = Chem.MolFromSmiles("CC")
    assert _to_mol(m2) is m2
    # invalid SMILES raises
    with pytest.raises(ValueError):
        _to_mol("not_a_smiles")

def test_add_and_contains_and_len():
    lib = FragmentLibrary()
    assert len(lib) == 0
    lib.add("CCO")
    assert len(lib) == 1
    # exact and canonical membership
    assert "CCO" in lib
    assert "OCC" in lib  # canonical SMILES of CCO
    # adding duplicate does not increase length
    lib.add("OCC")
    assert len(lib) == 1

def test_from_iterable_and_iter_and_smiles_property():
    frags = ["CC", "CO", "N"]
    lib = FragmentLibrary.from_iterable(frags)
    # N is not a valid fragment, so it should be skipped
    assert len(lib) == 2
    # __iter__ yields RDKit Mol objects with matching smiles
    smiles_list = lib.smiles
    mols = list(lib)
    assert len(mols) == len(smiles_list)
    for mol, smi in zip(mols, smiles_list):
        assert isinstance(mol, Chem.Mol)
        assert Chem.MolToSmiles(mol, isomericSmiles=True) == smi

def test_substructure_matches_and_contains_substructure():
    lib = FragmentLibrary.from_iterable(["CC", "CO"])
    # default return_atom_maps=False
    hits = lib.substructure_matches("CCO")
    assert len(hits) == 2
    found = {h[0] for h in hits}
    assert found == {"CC", "CO"}
    # atom maps empty by default
    assert all(h[3] == () for h in hits)
    # return_atom_maps=True yields non-empty tuples
    hits_maps = lib.substructure_matches("CCO", return_atom_maps=True)
    assert all(isinstance(h[3], tuple) and len(h[3]) > 0 for h in hits_maps)
    # contains_substructure boolean wrapper
    assert lib.contains_substructure("CCO") is True
    assert lib.contains_substructure("NNN") is False

def test_get_similar_with_thresholds():
    lib = FragmentLibrary.from_iterable(["CCO", "CCC"])
    # threshold=1.0 should only return exact same fragment
    sims = lib.get_similar("CCO", threshold=1.0)
    sims_smiles = [Chem.MolToSmiles(m, isomericSmiles=True) for m in sims]
    assert sims_smiles == ["CCO"]
    # threshold=0 should return all fragments
    sims0 = lib.get_similar("CCO", threshold=0.0)
    sims0_smiles = set(Chem.MolToSmiles(m, isomericSmiles=True) for m in sims0)
    assert sims0_smiles == {"CCO", "CCC"}

def test_featurise_binary_and_count():
    lib = FragmentLibrary.from_iterable(["CC", "CO"])
    # binary mode
    v_bin = lib.featurise("CCO")
    assert list(v_bin) == [1, 1]
    # count mode: "CCCC" has three CC substructures, zero CO
    v_cnt = lib.featurise("CCCC", mode="count")
    assert list(v_cnt) == [3, 0]

def test_set_operations_and_repr():
    libA = FragmentLibrary.from_iterable(["CC", "CO"])
    libB = FragmentLibrary.from_iterable(["CO", "N"])
    u = libA.union(libB)
    # N is not a valid fragment, so it should be skipped
    assert set(u.smiles) == {"CC", "CO"}
    assert repr(u) == f"<FragmentLibrary n={len(u)}>"
    inter = libA.intersection(libB)
    assert set(inter.smiles) == {"CO"}
    diff = libA.difference(libB)
    assert set(diff.smiles) == {"CC"}

def test_to_smiles_and_to_json(tmp_path):
    lib = FragmentLibrary.from_iterable(["CC", "CO"])
    # to_smiles writes one-per-line
    smi_file = tmp_path / "frags.smi"
    lib.to_smiles(smi_file)
    lines = smi_file.read_text().splitlines()
    assert lines == lib.smiles
    # to_json writes metadata
    json_file = tmp_path / "frags.json"
    lib.to_json(json_file)
    data = json.loads(json_file.read_text())
    assert data["n_fragments"] == len(lib)
    assert data["canonical_smiles"] == lib.smiles

def test_from_file_csv_and_errors(tmp_path):
    # valid CSV with duplicates
    csv_path = tmp_path / "data.csv"
    csv_path.write_text("smiles,val\nCCO,1\nCCC,2\nCCO,3")
    lib = FragmentLibrary.from_file(csv_path, smiles_column="smiles")
    # duplicate CCO should only appear once
    assert set(lib.smiles) == {"CCO", "CCC"}
    assert len(lib) == 2
    # missing smiles_column raises
    bad_csv = tmp_path / "bad.csv"
    bad_csv.write_text("a,b\n1,2")
    with pytest.raises(ValueError):
        FragmentLibrary.from_file(bad_csv)

def test_minifrag_library():
    lib = MiniFragsLib()
    # Check basic properties
    assert len(lib) > 0
    assert all(isinstance(m, Chem.Mol) for m in lib.mols)
    assert all(isinstance(s, str) for s in lib.smiles)
    # Check repr
    assert repr(lib).startswith("<MiniFragsLibrary n=")


def test_poised_library():
    lib = PoisedFragsLib()
    # Check basic properties
    assert len(lib) > 0
    assert all(isinstance(m, Chem.Mol) for m in lib.mols)
    assert all(isinstance(s, str) for s in lib.smiles)
    # Check repr
    assert repr(lib).startswith("<PoisedFragsLibrary n=")
