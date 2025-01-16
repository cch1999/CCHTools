import os

# Get the project directory
PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Path to the examples directory
EXAMPLES_DIR = os.path.join(PROJECT_DIR, "examples")
EXAMPLES_PDB = os.path.join(EXAMPLES_DIR, "6w63.pdb")
EXAMPLES_SDF = os.path.join(EXAMPLES_DIR, "6w63_ref_ligand.sdf")
EXAMPLES_CIF = os.path.join(EXAMPLES_DIR, "HEM.cif")
EXAMPLES_PYMOL_IMAGE_1 = os.path.join(EXAMPLES_DIR, "seh_generated_0.png")
EXAMPLES_PYMOL_IMAGE_2 = os.path.join(EXAMPLES_DIR, "seh_ref_ligand.png")
