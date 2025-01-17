PyMOL Interface
==============

The ``PyMOLInterface`` class provides a convenient wrapper around PyMOL commands for molecular visualization.

Basic Usage
----------

.. code-block:: python

    from cchtools.pymol_interface import PyMOLInterface
    
    # Initialize the interface
    pymol = PyMOLInterface(background_color="white", view_size=(400, 400))
    
    # Load a protein structure
    pymol.fetch_protein("1abc", show_as="cartoon", color="cyan")
    
    # Load a molecule from an RDKit mol object
    pymol.load_molecule(mol, color="green")
    
    # Save the visualization
    pymol.save_image("output.png")

API Reference
------------

.. autoclass:: cchtools.pymol_interface.PyMOLInterface
   :members:
   :undoc-members:
   :show-inheritance:

Examples
--------

Loading a Protein and Ligand
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    pymol = PyMOLInterface()
    pymol.set_background_color()
    
    # Load protein
    pymol.load_pdb("protein.pdb", transparency=0.5, show_as="cartoon")
    
    # Load ligand
    ligand = dm.read_sdf("ligand.sdf")[0]
    pymol.load_molecule(ligand, color="green")
    
    # Show interacting residues
    residues = pymol.get_interacting_residues(distance=3.0)
    pymol.show_residue_sidechains(residues)
    
    # Save the result
    pymol.save_image("complex.png")

Visualization Settings
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    pymol = PyMOLInterface()
    
    # Customize view
    pymol.set_background_color()
    pymol.set_ray_trace_mode(3)
    
    # Set specific view coordinates
    view = [0.328213751, 0.371481866, -0.868490517, ...]
    pymol.set_view(view) 