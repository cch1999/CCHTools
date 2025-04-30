import logging
import os
import sys

from rdkit import Chem, DataStructs, RDConfig
from rdkit.Chem import QED, AllChem

sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer

from cchtools.utils.cif import fetch_ideal_ccd_structure

# Set up basic logging configuration
logging.basicConfig(level=logging.INFO)


def fix_bond_orders_with_ccd(mol: Chem.Mol, ccd_code: str, tmp_dir: str = "/tmp") -> Chem.Mol:
    """
    Fix the bond orders of an RDKit molecule using the ideal CCD structure.

    Args:
        mol (rdkit.Chem.Mol): The molecule to fix.
        ccd_code (str): The CCD code of the molecule.
        tmp_dir (str): The directory to save the ideal CCD structure.

    Returns:
        rdkit.Chem.Mol: The molecule with fixed bond orders.
    """

    ideal_mol = fetch_ideal_ccd_structure(ccd_code, tmp_dir)

    if ideal_mol is None:
        return mol  # Return original molecule if ideal structure couldn't be fetched

    # Transfer bond orders from ideal molecule to the main molecule
    fixed_mol = AllChem.AssignBondOrdersFromTemplate(ideal_mol, mol)

    if fixed_mol is None:
        logging.warning(f"Failed to fix bond orders for {ccd_code}")
        return mol  # Return original molecule if bond order transfer failed

    return fixed_mol


def calculate_basic_metrics(
    mols: list[Chem.Mol], ref_mol: Chem.Mol | str | None = None
) -> list[dict]:
    """
    Calculate basic metrics for a list of molecules, optionally comparing to a reference molecule.

    Args:
        mols (list[Chem.Mol]): List of RDKit molecules to analyze.
        ref_mol (Chem.Mol | str | None, optional): Reference molecule for similarity comparison.
            Can be an RDKit molecule, a SMILES string, or None. Defaults to None.

    Returns:
        list[dict]: List of dictionaries containing metrics for each molecule.
            Each dictionary contains 'QED', 'SA', and optionally 'Similarity' scores.
    """

    # Convert reference molecule from SMILES to RDKit Mol object if necessary
    if isinstance(ref_mol, str):
        ref_mol = Chem.MolFromSmiles(ref_mol, sanitize=False)

    metrics = []
    for mol in mols:
        # Calculate Quantitative Estimate of Drug-likeness (QED)
        qed = QED.qed(mol)

        # Calculate Synthetic Accessibility (SA) score
        sa = sascorer.calculateScore(mol)

        # Calculate similarity to reference molecule if provided
        if ref_mol is not None:
            sim = DataStructs.TanimotoSimilarity(
                Chem.RDKFingerprint(mol), Chem.RDKFingerprint(ref_mol)
            )
            sim = round(sim, 2)
        else:
            sim = None

        # Append calculated metrics to the list
        metrics.append({"QED": round(qed, 2), "SA": round(sa, 2), "Similarity": sim})

    return metrics


if __name__ == "__main__":
    # download the ideal structure for HEM
    from typing import Optional
    from cchtools.utils.cif import download_ideal_ccd_structure, cif_to_rdkit

    cif_path: Optional[str] = download_ideal_ccd_structure("HEM")
    if cif_path is not None:
        print(f"Downloaded CIF file to {cif_path}")

    # convert to rdkit
    if cif_path is not None:
        mol = cif_to_rdkit(cif_path)
        print(mol)

    # Fetch the ideal structure for HEM
    mol = fetch_ideal_ccd_structure("HEM")
    print(mol)

    # corrupt the bond order of HEM by setting all bonds to single
    if mol is not None:
        for bond in mol.GetBonds():
            bond.SetBondType(Chem.BondType.SINGLE)
        print(mol)
    else:
        print("Warning: mol is None, skipping bond processing")

    # fix the bond order of HEM
    mol = fix_bond_orders_with_ccd(mol, "HEM")
    print(mol)

    # Test the calculate_basic_metrics function
    smile = Chem.MolToSmiles(mol)
    metrics = calculate_basic_metrics([mol, mol], smile)
    print(metrics)
