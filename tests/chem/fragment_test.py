import json
import pickle
import os

import pytest
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from cchtools.chem.fragment_library import (
    FragmentLibrary,
    MiniFragsLib,
    PoisedFragsLib,
    _to_mol,
    standardise,
)


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


def test_standardise_handles_edge_cases():
    # Test None input
    assert standardise(None) is None
    
    # Test valid molecule
    mol = Chem.MolFromSmiles("CCO")
    result = standardise(mol)
    assert isinstance(result, Chem.Mol)
    
    # Test molecule with salt
    salt_mol = Chem.MolFromSmiles("CCO.Cl")
    result = standardise(salt_mol)
    assert isinstance(result, Chem.Mol)
    assert Chem.MolToSmiles(result) == "CCO"
    
    # Test duplicate fragments
    dup_mol = Chem.MolFromSmiles("CCO.CCO")
    result = standardise(dup_mol)
    assert isinstance(result, Chem.Mol)
    assert Chem.MolToSmiles(result) == "CCO"


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


def test_from_molecules():
    # Create some RDKit molecules with properties
    mols = []
    for smiles, name in [("CCO", "ethanol"), ("CCC", "propane")]:
        mol = Chem.MolFromSmiles(smiles)
        # Set properties as Chem.Props._ prefix is stripped during serialization
        mol.SetProp("Name", name)  # Changed from _Name to Name
        mol.SetProp("MW", str(round(rdMolDescriptors.CalcExactMolWt(mol), 2)))
        mols.append(mol)
        
    # Test without extracting properties
    lib = FragmentLibrary.from_molecules(mols, extract_props=False)
    assert len(lib) == 2
    assert set(lib.smiles) == {"CCO", "CCC"}
    assert not lib.metadata
    
    # Test with extracting properties
    lib_with_props = FragmentLibrary.from_molecules(mols, extract_props=True)
    assert len(lib_with_props) == 2
    assert "properties" in lib_with_props.metadata
    assert len(lib_with_props.metadata["properties"]) == 2
    
    # Check properties were extracted correctly
    ethanol_smi = Chem.MolToSmiles(Chem.MolFromSmiles("CCO"), isomericSmiles=True)
    assert lib_with_props.metadata["properties"][ethanol_smi]["Name"] == "ethanol"  # Changed from _Name to Name


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


def test_get_similar_with_scores():
    lib = FragmentLibrary.from_iterable(["CCO", "CCC", "c1ccccc1"])
    
    # Get similar fragments with scores
    sims = lib.get_similar_with_scores("CCO", threshold=0.0)
    
    # Check results format
    assert len(sims) == 3
    assert all(isinstance(mol, Chem.Mol) for mol, _ in sims)
    assert all(isinstance(score, float) for _, score in sims)
    
    # Check sorting (should be descending by score)
    scores = [score for _, score in sims]
    assert sorted(scores, reverse=True) == scores
    
    # The query itself should have highest similarity (1.0)
    top_mol, top_score = sims[0]
    assert Chem.MolToSmiles(top_mol) == "CCO"
    assert top_score == 1.0


def test_featurise_binary_and_count():
    lib = FragmentLibrary.from_iterable(["CC", "CO"])
    # binary mode
    v_bin = lib.featurise("CCO")
    assert list(v_bin) == [1, 1]
    # count mode: "CCCC" has three CC substructures, zero CO
    v_cnt = lib.featurise("CCCC", mode="count")
    assert list(v_cnt) == [3, 0]
    
    # Test invalid mode
    with pytest.raises(ValueError):
        lib.featurise("CCO", mode="invalid")


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


def test_operator_add():
    libA = FragmentLibrary.from_iterable(["CC", "CO"])
    libB = FragmentLibrary.from_iterable(["CO", "CCN"])
    
    # Test addition operator
    combined = libA + libB
    assert isinstance(combined, FragmentLibrary)
    assert len(combined) == 3
    assert set(combined.smiles) == {"CC", "CO", "CCN"}


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
    
    # Test with custom metadata
    lib._metadata["source"] = "test"
    lib.to_json(json_file)
    data = json.loads(json_file.read_text())
    assert "metadata" in data
    assert data["metadata"]["source"] == "test"


def test_to_pickle_and_from_pickle(tmp_path):
    # Create a library with some molecules
    lib = FragmentLibrary.from_iterable(["CC", "CO", "CCO"])
    lib._metadata["source"] = "test pickle"
    
    # Save to pickle
    pickle_path = tmp_path / "frags.pkl"
    lib.to_pickle(pickle_path)
    assert os.path.exists(pickle_path)
    
    # Load from pickle
    loaded_lib = FragmentLibrary.from_pickle(pickle_path)
    
    # Check if loaded library has the same content
    assert len(loaded_lib) == len(lib)
    assert loaded_lib.smiles == lib.smiles
    assert loaded_lib.metadata == lib.metadata
    
    # Test error handling for invalid pickle
    invalid_path = tmp_path / "invalid.pkl"
    with open(invalid_path, 'wb') as f:
        pickle.dump("not a library", f)
    
    with pytest.raises(ValueError):
        FragmentLibrary.from_pickle(invalid_path)


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


def test_error_handling():
    lib = FragmentLibrary()
    
    # Invalid SMILES should be handled gracefully
    lib.add("not_a_smiles")
    assert len(lib) == 0
    
    # Adding None should be handled
    lib.add(None)
    assert len(lib) == 0


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
