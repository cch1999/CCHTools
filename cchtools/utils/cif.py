import requests
from rdkit import Chem
from rdkit.Chem import AllChem
import os
from pathlib import Path
from typing import Optional, Union

from pdbeccdutils.core import ccd_reader

import logging

# Set up basic logging configuration
logging.basicConfig(level=logging.INFO)


def download_ideal_ccd_structure(ccd_code: str, tmp_dir: str = '/tmp') -> Optional[str]:
    """
    Download the ideal CCD structure for a given component code.

    Args:
        ccd_code (str): The CCD component code.
        tmp_dir (str): Directory to save the downloaded file. Defaults to '/tmp'.

    Returns:
        Optional[str]: Path to the downloaded CIF file, or None if download failed.
    """
    url = f'https://files.rcsb.org/ligands/download/{ccd_code}.cif'
    cif_path = os.path.join(tmp_dir, f'{ccd_code}.cif')

    # Return existing file if already downloaded
    if os.path.exists(cif_path):
        return cif_path
    
    # Attempt to download the CIF file
    response = requests.get(url)
    if response.status_code != 200:
        logging.error(f"Failed to fetch CIF file for {ccd_code}")
        return None
    
    # Save the downloaded CIF file
    with open(cif_path, 'w') as f:
        f.write(response.text)
    return cif_path


def cif_to_rdkit(cif_path: Union[str, Path], sanitize: bool = True) -> Optional[Chem.Mol]:
    """
    Convert a CIF file to an RDKit molecule.

    Args:
        cif_path (Union[str, Path]): Path to the CIF file.
        sanitize (bool): Whether to sanitize the molecule. Defaults to True.

    Returns:
        Optional[Chem.Mol]: RDKit molecule object, or None if conversion failed.
    """
    ccd_component = ccd_reader.read_pdb_cif_file(cif_path, sanitize=sanitize).component
    return ccd_component.mol


def fetch_ideal_ccd_structure(ccd_code: str, tmp_dir: str = '/tmp', save_structures: bool = True) -> Optional[Chem.Mol]:
    """
    Fetch the ideal CCD structure and convert it to an RDKit molecule.

    This function will first attempt to download the CIF file for the given CCD code.
    If the file already exists in the specified directory, it will be used without re-downloading.
    The function will then convert the CIF file to an RDKit molecule and return the molecule object.
    If the download or conversion fails, the function will return None.

    Args:
        ccd_code (str): The CCD component code.
        tmp_dir (str): Directory to save the downloaded file. Defaults to '/tmp'.
        save_structures (bool): Whether to keep the downloaded CIF file. Defaults to True.

    Returns:
        Optional[Chem.Mol]: RDKit molecule object, or None if fetching or conversion failed.
    """
    cif_path = download_ideal_ccd_structure(ccd_code, tmp_dir)
    if not cif_path:
        return None

    rdkit_mol = cif_to_rdkit(cif_path)
    
    # Remove the CIF file if save_structures is False
    if not save_structures:
        os.remove(cif_path)

    return rdkit_mol