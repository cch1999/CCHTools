"""Functions to interacting with the PDBe Knowledge Graph API."""

import requests
from typing import Dict, Any, List

def _url_to_json(url: str) -> Dict[str, Any]:
    response = requests.get(url)
    if response.status_code == 404:
        raise ValueError(f"Query at {url} not found in the PDBe Knowledge Graph API.")
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
    """Get atom information for a specific compound from the PDBe Knowledge Graph API.
    
    Args:
        hetcode: The hetcode/compound ID from the PDB Chemical Component Dictionary
        
    Returns:
        List of dictionaries containing atom information including:
        - stereo: R/S stereochemistry for atoms (or E/Z for bonds) when applicable
        - leaving_atom: Flag indicating if this is a leaving atom
        - pdb_name: Alternate name of the atom
        - aromatic: Whether the atom is within an aromatic substructure
        - element: Element type of the atom
        - ideal_x/y/z: Calculated idealized coordinates
        - charge: Formal charge on the atom
        - atom_name: Atom identifier from chemical components dictionary

    Example:
        >>> atoms = get_compound_atoms("ADP")
        >>> atoms[0]
            {
        'stereo': False,
        'leaving_atom': False,
        'pdb_name': 'C13',
        'aromatic': True,
        'element': 'C',
        'ideal_y': 0.0,
        'ideal_x': 0.296,
        'charge': 0.0,
        'ideal_z': -4.67,
            'atom_name': 'C13'
        }
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

    Example:
        >>> bonds = get_compound_bonds("ADP")
        >>> bonds[0]
        {
            'stereo': False,
            'atom_1': "C1'",
            'atom_2': "H1'", 
            'bond_type': 'sing',
            'bond_order': 1,
            'ideal_length': 1.09,
            'aromatic': False
        }
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

    Example:
        >>> similar = get_similar_compounds("TDP")
        >>> similar[0]
        {
            'acts_as': 'cofactor',
            'class': 'Thiamine diphosphate',
            'chem_comp_ids': [
                {
                    'chem_comp_id': 'WWF',
                    'name': 'C2-1-HYDROXY-3-METHYL-BUTYL-THIAMIN'
                },
                ...
            ]
        }
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

    Example:
        >>> similar = get_similar_ligands("STI")
        >>> similar[0]["same_scaffold"][0]
        {
            'chem_comp_id': 'MPZ',
            'name': '4-[(4-METHYLPIPERAZIN-1-YL)METHYL]-N-{3-[(4-PYRIDIN-3-YLPYRIMIDIN-2-YL)AMINO]PHENYL}BENZAMIDE',
            'substructure_match': ['O16', 'C33', 'N34'],
            'similarity_score': 0.973
        }
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

    Example:
        >>> summary = get_compound_summary("ATP")
        >>> summary[0]["name"]
        'ADENOSINE-5\'-TRIPHOSPHATE'
        >>> summary[0]["formula"]
        'C10 H16 N5 O13 P3'
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
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/compound/cofactors/"
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
