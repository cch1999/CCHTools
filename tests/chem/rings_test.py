import pytest
from rdkit import Chem

from cchtools.chem.rings import generate_ring_systems
from cchtools.data.constants import RING_ATOMS, RING_SKELETONS


def test_generate_ring_systems():
    reduced_ring_skeletons = {
        "6": "c1ccccc1",
    }

    df = generate_ring_systems(ring_skeletons=reduced_ring_skeletons, atom_types=RING_ATOMS)

    # assert df has multiple examples of each skeleton
    assert df["skeleton"].nunique() == len(reduced_ring_skeletons)


if __name__ == "__main__":
    test_generate_ring_systems()
