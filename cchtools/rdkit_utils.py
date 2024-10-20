import logging

from rdkit import Chem
from rdkit.Chem import AllChem

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
