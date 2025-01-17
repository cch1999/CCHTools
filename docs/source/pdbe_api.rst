PDBe Knowledge Graph API
=======================

The PDBe Knowledge Graph API module provides functions for interacting with the Protein Data Bank in Europe (PDBe) Knowledge Graph API.

Basic Usage
----------

.. code-block:: python

    from cchtools.utils.pdbe import get_pdbs_with_compound, get_compound_summary
    
    # Get PDB entries containing a specific compound
    pdbs = get_pdbs_with_compound("ATP")
    
    # Get detailed information about a compound
    summary = get_compound_summary("ATP")

API Reference
------------

.. automodule:: cchtools.utils.pdbe
   :members:
   :undoc-members:
   :show-inheritance:

Examples
--------

Getting Similar Compounds
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from cchtools.utils.pdbe import get_similar_ligands
    
    # Find structurally similar ligands
    similar = get_similar_ligands("ATP")
    
    # Print similar compounds and their similarity scores
    for compound in similar[0].get('similar_ligands', []):
        print(f"{compound['chem_comp_id']}: {compound['similarity_score']}")

Analyzing Substructures
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from cchtools.utils.pdbe import get_compound_substructures
    
    # Get substructure information
    substructures = get_compound_substructures("ATP")
    
    # Print fragment information
    for fragment_name, atoms in substructures[0]['fragments'].items():
        print(f"{fragment_name}: {len(atoms)} atoms")
