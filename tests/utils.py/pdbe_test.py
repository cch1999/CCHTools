import pytest
import cchtools.utils.pdbe as pdbe

def hetcodes():
    return ["ATP", "ADP", "355", "STI", "HEM"]

def hetcode_with_functional_annotation():
    # hetcodes that have functional annotations
    return ["COM", "HEM"]

@pytest.mark.parametrize("hetcode", hetcodes())
def test_get_pdbs_with_compound(hetcode):
    # Test with known compound
    pdbs = pdbe.get_pdbs_with_compound(hetcode)
    assert isinstance(pdbs, list)
    assert len(pdbs) > 0
    assert all(isinstance(pdb, str) for pdb in pdbs)


@pytest.mark.parametrize("hetcode", hetcodes())
def test_get_compound_atoms(hetcode):
    # Test with known compound
    atoms = pdbe.get_compound_atoms(hetcode)
    assert isinstance(atoms, list)
    assert len(atoms) > 0
    
    # Check required fields exist
    required_fields = ['stereo', 'leaving_atom', 'pdb_name', 'aromatic', 
                      'element', 'ideal_x', 'ideal_y', 'ideal_z',
                      'charge', 'atom_name']
    for atom in atoms:
        for field in required_fields:
            assert field in atom


@pytest.mark.parametrize("hetcode", hetcodes())
def test_get_compound_bonds(hetcode):
    # Test with known compound
    bonds = pdbe.get_compound_bonds(hetcode)
    assert isinstance(bonds, list)
    assert len(bonds) > 0
    
    # Check required fields exist
    required_fields = ['stereo', 'atom_1', 'atom_2', 'bond_type',
                      'bond_order', 'ideal_length', 'aromatic']
    for bond in bonds:
        for field in required_fields:
            assert field in bond

@pytest.mark.parametrize("hetcode", hetcode_with_functional_annotation())
def test_get_similar_hetcodes(hetcode):
    # Test with known compound
    similar = pdbe.get_similar_hetcodes(hetcode)
    assert isinstance(similar, list)
    assert len(similar) > 0
    
    # Check required fields exist
    required_fields = ['acts_as', 'class', 'chem_comp_ids']
    for entry in similar:
        for field in required_fields:
            assert field in entry

@pytest.mark.parametrize("hetcode", hetcodes())
def test_get_similar_ligands(hetcode):
    # Test with known compound
    similar = pdbe.get_similar_ligands(hetcode)
    assert isinstance(similar, list)
    assert len(similar) > 0 if hetcode != '355' else len(similar) == 0
    
    if len(similar) > 0:
        # Check required fields exist
        required_fields = ['stereoisomers', 'same_scaffold', 'similar_ligands']
        for entry in similar[0].keys():
            assert entry in required_fields

@pytest.mark.parametrize("hetcode", hetcodes())
def test_get_compound_substructures(hetcode):
    # Test with known compound
    substructures = pdbe.get_compound_substructures(hetcode)
    assert isinstance(substructures, list)
    assert len(substructures) > 0
    
    # Check required fields exist
    required_fields = ['fragments', 'scaffold']
    for entry in substructures:
        for field in required_fields:
            assert field in entry

@pytest.mark.parametrize("hetcode", hetcodes())
def test_get_compound_summary(hetcode):
    # Test with known compound
    summary = pdbe.get_compound_summary(hetcode)
    assert isinstance(summary, list)
    assert len(summary) > 0
    
    # Check required fields exist
    required_fields = ['name', 'formula', 'inchi', 'inchi_key', 'smiles',
                      'cross_links', 'synonyms', 'phys_chem_properties']
    for entry in summary:
        for field in required_fields:
            assert field in entry

def test_get_cofactor_summary():
    # Test basic functionality
    cofactors = pdbe.get_cofactor_summary()
    assert isinstance(cofactors, dict)
    assert len(cofactors) > 0
    
    # Check structure and content of returned data
    for cofactor_class, entries in cofactors.items():
        assert isinstance(cofactor_class, str)
        assert isinstance(entries, list)
        for entry in entries:
            assert 'EC' in entry
            assert 'cofactors' in entry
            assert isinstance(entry['EC'], list)
            assert isinstance(entry['cofactors'], list)
