import pytest
from rdkit import Chem

from cchtools.rdkit_utils import calculate_basic_metrics


@pytest.fixture
def sample_molecules():
    smiles = ["CC(=O)O", "C1=CC=CC=C1", "CCO"]  # Acetic acid  # Benzene  # Ethanol
    return [Chem.MolFromSmiles(smile) for smile in smiles]


@pytest.fixture
def reference_molecule():
    return Chem.MolFromSmiles("CC(=O)O")  # Acetic acid


def test_calculate_basic_metrics(sample_molecules):
    metrics = calculate_basic_metrics(sample_molecules)

    assert len(metrics) == len(sample_molecules)
    for metric in metrics:
        assert "QED" in metric
        assert "SA" in metric
        assert metric["Similarity"] is None
        assert 0 <= metric["QED"] <= 1
        assert metric["SA"] > 0


def test_calculate_basic_metrics_with_reference(sample_molecules, reference_molecule):
    metrics = calculate_basic_metrics(sample_molecules, reference_molecule)

    assert len(metrics) == len(sample_molecules)
    for metric in metrics:
        assert "QED" in metric
        assert "SA" in metric
        assert "Similarity" in metric
        assert 0 <= metric["QED"] <= 1
        assert metric["SA"] > 0
        assert 0 <= metric["Similarity"] <= 1


def test_calculate_basic_metrics_with_smiles_reference(sample_molecules):
    ref_smiles = "CC(=O)O"  # Acetic acid
    metrics = calculate_basic_metrics(sample_molecules, ref_smiles)

    assert len(metrics) == len(sample_molecules)
    for metric in metrics:
        assert "QED" in metric
        assert "SA" in metric
        assert "Similarity" in metric
        assert 0 <= metric["QED"] <= 1
        assert metric["SA"] > 0
        assert 0 <= metric["Similarity"] <= 1


def test_calculate_basic_metrics_empty_list():
    metrics = calculate_basic_metrics([])
    assert metrics == []


def test_calculate_basic_metrics_invalid_molecule():
    invalid_mol = Chem.MolFromSmiles("Invalid SMILES")
    with pytest.raises(Exception):
        calculate_basic_metrics([invalid_mol])
