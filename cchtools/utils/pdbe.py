"""Functions to interacting with the PDBe Knowledge Graph API."""

from typing import Any, Dict, List

import requests


def _url_to_json(url: str) -> Dict[str, Any]:
    response = requests.get(url)
    # TODO: This is removed because the API returns 404 for real CCDs that have 0 annotations
    # This is a problem because it makes it difficult to check if a hetcode is in the CCD
    # Only solution is to check if the hetcode is in the CCD before querying the API?
    # if response.status_code == 404:
    #    raise ValueError(f"Query at {url} not found in the PDBe Knowledge Graph API.")
    return response.json()


def get_pdbs_with_compound(hetcode: str) -> List[str]:
    """Get PDB entries containing a specific compound from the PDBe Knowledge Graph API.

    Args:
        hetcode: The hetcode/compound ID from the PDB Chemical Component Dictionary

    Returns:
        List of PDB entry IDs that contain the specified compound

    Example:
        >>> get_pdbs_with_compound("3IP")
        ['1w7h', '3fty', '8qul']
    """
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/compound/in_pdb/{hetcode}"
    response = _url_to_json(url)
    return response.get(hetcode, [])


def get_compound_atoms(hetcode: str) -> List[Dict[str, Any]]:
    """Get atom information for a chemical component.

    Args:
        hetcode: The PDB chemical component code

    Returns:
        List of dictionaries containing atom information with fields:
        - stereo: Whether the atom is a stereocenter
        - leaving_atom: Whether the atom can leave during reactions
        - pdb_name: Alternate name of the atom
        - aromatic: Whether the atom is within an aromatic substructure
        - element: Element type of the atom
        - ideal_x/y/z: Calculated idealized coordinates
        - charge: Formal charge on the atom
        - atom_name: Atom identifier from chemical components dictionary

    Example:
        >>> atoms = get_compound_atoms("ADP")
        >>> sorted(atoms[0].keys())  # Show fields in a stable order
        ['aromatic', 'atom_name', 'charge', 'element', 'ideal_x', 'ideal_y', 'ideal_z', 'leaving_atom', 'pdb_name', 'stereo']
        >>> all(isinstance(x, bool) for x in [atoms[0]['stereo'], atoms[0]['leaving_atom'], atoms[0]['aromatic']])
        True
    """
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/compound/atoms/{hetcode}"
    response = _url_to_json(url)
    return response.get(hetcode, [])


def get_compound_bonds(hetcode: str) -> List[Dict[str, Any]]:
    """Get bond information for a specific compound from the PDBe Knowledge Graph API.

    Args:
        hetcode: The hetcode/compound ID from the PDB Chemical Component Dictionary

    Returns:
        List of dictionaries containing bond information including:
        - stereo: R/S stereochemistry for atoms (or E/Z for bonds) when applicable
        - atom_1: Name of first atom in the bond
        - atom_2: Name of second atom in the bond
        - bond_type: String describing bond type (e.g. 'sing', 'doub')
        - bond_order: Integer describing bond order
        - ideal_length: Target length of the bond
        - aromatic: Whether the bond is within an aromatic substructure
    """
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/compound/bonds/{hetcode}"
    response = _url_to_json(url)
    return response.get(hetcode, [])


def get_similar_hetcodes(hetcode: str) -> List[Dict[str, Any]]:
    """Get similar compounds to a specific hetcode from the PDBe Knowledge Graph API.

    NOTE: Returns hetcodes that share functional annotations with the query hetcode.


    Args:
        hetcode: The hetcode/compound ID from the PDB Chemical Component Dictionary

    Returns:
        List of dictionaries containing similar compound information including:
        - acts_as: Role of the compound (e.g. 'cofactor', 'reactant')
        - class: Classification of the compound group
        - chem_comp_ids: List of similar compounds with their hetcodes and names
    """
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/compound/cofactors/het/{hetcode}"
    response = _url_to_json(url)
    return response.get(hetcode, [])


def get_similar_ligands(hetcode: str) -> List[Dict[str, Any]]:
    """Get structurally similar ligands to a specific hetcode from the PDBe Knowledge Graph API.

    This returns ligands that are:
    - Stereoisomers of the query ligand
    - Share the same molecular scaffold
    - Have structural similarity >60% (using PARITY method)

    Args:
        hetcode: The hetcode/compound ID from the PDB Chemical Component Dictionary

    Returns:
        List of dictionaries containing similar compound information including:
        - stereoisomers: List of stereoisomer hetcode objects
        - same_scaffold: List of compounds sharing the same molecular scaffold
        - similar_ligands: List of structurally similar compounds
        Each compound entry includes:
        - chem_comp_id: The compound identifier
        - name: Chemical name
        - similarity_score: PARITY similarity score (0-1)
        - substructure_match: List of matching atom names
    """
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/compound/similarity/{hetcode}"
    response = _url_to_json(url)
    return response.get(hetcode, [])


def get_compound_substructures(hetcode: str) -> List[Dict[str, Any]]:
    """Get substructure information for a specific hetcode from the PDBe Knowledge Graph API.

    This returns information about molecular scaffolds and fragments found in the compound structure.

    Args:
        hetcode: The hetcode/compound ID from the PDB Chemical Component Dictionary

    Returns:
        List of dictionaries containing substructure information including:
        - fragments: Dictionary mapping fragment names to lists of atom names
        - scaffold: Dictionary mapping SMILES strings to lists of atom names

    Example:
        >>> substructures = get_compound_substructures("ATP")
        >>> substructures[0]["fragments"]["pyrimidine"][0][:3]
        ['C4', 'N3', 'C2']
        >>> list(substructures[0]["scaffold"].keys())[0]
        'c1ncc2ncn([C@H]3CCCO3)c2n1'
    """
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/compound/substructures/{hetcode}"
    response = _url_to_json(url)
    return response.get(hetcode, [])


def get_compound_summary(hetcode: str) -> List[Dict[str, Any]]:
    """Get summary information for a specific hetcode from the PDBe Knowledge Graph API.

    This returns comprehensive information about the chemical component including its
    chemical properties, cross references, and physical-chemical characteristics.

    Args:
        hetcode: The hetcode/compound ID from the PDB Chemical Component Dictionary

    Returns:
        List of dictionaries containing compound summary information including:
        - name: Chemical component name
        - formula: Chemical formula
        - inchi: InChI identifier
        - inchi_key: InChI key
        - smiles: SMILES representation
        - cross_links: External database references
        - synonyms: Alternative names
        - phys_chem_properties: Physical and chemical properties
    """
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/compound/summary/{hetcode}"
    response = _url_to_json(url)
    return response.get(hetcode, [])


def get_cofactor_summary() -> Dict[str, List[Dict[str, List[str]]]]:
    """Get summary information about all cofactor annotations in the PDB.

    This function retrieves comprehensive information about cofactor classes,
    including their associated enzyme classes (EC numbers) and member compounds.

    Returns:
        Dictionary mapping cofactor class names to lists of dictionaries containing:
        - EC: List of enzyme class numbers where these cofactors are active
        - cofactors: List of hetcodes that belong to this cofactor class

    Example:
        >>> cofactors = get_cofactor_summary()
        >>> cofactors["MIO"][0]["EC"][:2]
        ['4.3.1.23', '4.3.1.24']
        >>> cofactors["Coenzyme M"][0]["cofactors"]
        ['COM']
    """
    url = "https://www.ebi.ac.uk/pdbe/graph-api/compound/cofactors/"
    response = _url_to_json(url)
    return response


if __name__ == "__main__":
    from rich import print

    # #### COMPOUNDS
    # print(get_pdbs_with_compound("3IP"))
    # atoms = get_compound_atoms("3IP")
    # print(atoms[:2])
    # bonds = get_compound_bonds("3IP")
    # print(bonds[:5])
    # similar = get_similar_hetcodes("TDP")
    # print(similar)

    similar_ligands = get_similar_ligands("STI")
    print(similar_ligands)

    # substructures = get_compound_substructures("ATP")
    # print(substructures)

    # summary = get_compound_summary("ATP")
    # print(summary)

    # cofactors = get_cofactor_summary("COM")
    # print(cofactors)

    # cofactors = get_cofactor_summary("HEM")
    # print(cofactors)
