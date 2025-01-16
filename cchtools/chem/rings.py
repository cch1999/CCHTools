from itertools import product
from typing import Dict

import pandas as pd
from rdkit import Chem
from tqdm import tqdm

from cchtools.data.constants import RING_ATOMS, RING_SKELETONS


def generate_ring_systems(
    ring_skeletons: Dict[str, str] = RING_SKELETONS, atom_types: Dict[str, str] = RING_ATOMS
) -> pd.DataFrame:
    """
    Generate ring systems from the given ring skeletons and atom types.

    This function takes in dictionaries of ring skeletons and atom types, and generates
    all possible ring systems by substituting the aromatic carbons in the skeletons with
    the provided atom types. The generated ring systems are then validated based on
    specific rules and returned as a pandas DataFrame.

    Parameters:
    ring_skeletons (Dict[str, str]): A dictionary where keys are skeleton names and values
                                     are SMILES strings representing the ring skeletons.
    atom_types (Dict[str, str]): A dictionary where keys are atom type names and values
                                 are their corresponding SMILES representations.

    Returns:
    pd.DataFrame: A DataFrame containing the valid generated ring systems.
    """

    def is_valid_scaffold(mol):
        if mol is None:
            return False

        # Check valence
        try:
            Chem.SanitizeMol(mol)
        except:
            print(f"Warning: Could not sanitize mol {mol}")
            return False

        # Rules from paper:
        # 1. Max 4 heteroatoms
        hetero_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() != 6)
        if hetero_count > 4:
            return False

        # 2. Max 2 carbonyl groups
        carbonyl_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C=O")))
        if carbonyl_count > 2:
            return False

        # 3. No unstable bonds (O-O, S-S, N-N)
        unstable_patterns = ["[O,S,N]-[O,S,N]"]
        for pattern in unstable_patterns:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                return False

        return True

    scaffolds = []

    # Generate scaffolds for each skeleton
    for skeleton_name, skeleton_smiles in tqdm(ring_skeletons.items()):
        base_mol = Chem.MolFromSmiles(skeleton_smiles)
        if base_mol is None:
            print(f"Warning: Could not parse SMILES for skeleton {skeleton_name}")
            continue

        # Add aromatic carbons to the mol
        base_mol = Chem.AddHs(base_mol)

        # Generate all possible combinations of heteroatom replacements
        num_aromatic_carbons = len([atom for atom in base_mol.GetAtoms() if atom.GetIsAromatic()])

        print(num_aromatic_carbons)

        # Generate all possible combinations of heteroatom replacements
        for positions in product(atom_types.values(), repeat=num_aromatic_carbons):
            new_smiles = skeleton_smiles
            for i, replacement in enumerate(positions):
                new_smiles = new_smiles.replace("c", replacement, 1)

            try:
                mol = Chem.MolFromSmiles(new_smiles)
                if mol and is_valid_scaffold(mol):
                    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
                    scaffolds.append(
                        {
                            "skeleton": skeleton_name,
                            "smiles": canonical_smiles,
                            "num_atoms": mol.GetNumAtoms(),
                            "num_hetero": sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() != 6),
                        }
                    )
                else:
                    print(f"Warning: Invalid scaffold {new_smiles}")
            except:
                print(f"Warning: Could not parse SMILES for skeleton {skeleton_name}")
                continue

    return pd.DataFrame(scaffolds)


if __name__ == "__main__":
    from rich import print

    reduced_ring_skeletons = {
        "5": "c1cccc1",
        "6": "c1ccccc1",
    }

    # Generate the scaffolds
    df_scaffolds = generate_ring_systems(ring_skeletons=reduced_ring_skeletons)

    # Basic analysis
    print(f"Total number of scaffolds generated: {len(df_scaffolds)}")
    print("\nScaffolds by skeleton type:")
    print(df_scaffolds["skeleton"].value_counts())
    print("\nAverage number of heteroatoms:", df_scaffolds["num_hetero"].mean())

    # Save results
    df_scaffolds.to_csv("ring_scaffolds.csv", index=False)
