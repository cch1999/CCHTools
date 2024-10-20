import requests
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import os
from pathlib import Path
from typing import Optional, Union

from cchtools.utils.cif import download_ideal_ccd_structure, cif_to_rdkit, fetch_ideal_ccd_structure

import logging

# Set up basic logging configuration
logging.basicConfig(level=logging.INFO)


def fix_bond_orders_with_ccd(mol: Chem.Mol, ccd_code: str, tmp_dir: str = '/tmp') -> Chem.Mol:
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



if __name__ == '__main__':  

    # download the ideal structure for HEM
    cif_path = download_ideal_ccd_structure('HEM')
    print(f"Downloaded CIF file to {cif_path}")

    # convert to rdkit
    mol = cif_to_rdkit(cif_path)
    print(mol)

    # Fetch the ideal structure for HEM
    mol = fetch_ideal_ccd_structure('HEM')
    print(mol)

    # corrupt the bond order of HEM by setting all bonds to single
    for bond in mol.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
    print(mol)

    # fix the bond order of HEM
    mol = fix_bond_orders_with_ccd(mol, 'HEM')
    print(mol)
    
        
