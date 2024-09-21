import os
from io import BytesIO

import datamol as dm
import pandas as pd
import rdkit.Chem as Chem
from PIL import Image
from pymol import cmd
import tempfile

from cchtools.constants import EXAMPLES_PDB


class PyMOLInterface:
    """
    A class to interface with PyMOL for molecular visualization.
    """

    def __init__(self, background_color: str = "white"):
        """
        Initializes the PyMOL interface.
        """
        self.cmd: cmd = cmd
        self.background_color = background_color

    def delete_all(self):
        self.cmd.delete("all")

    def fetch_protein(self, pdb_code: str):
        """
        Fetches a protein structure by its PDB code and displays it as a cartoon.

        Args:
            pdb_code (str): The PDB code of the protein to fetch.
        """
        print(f"Fetching protein {pdb_code}")

        self.cmd.fetch(pdb_code)
        self.cmd.hide("everything", "all")
        self.cmd.show("cartoon", "all")
        self.cmd.color("cyan", "all")
        self.cmd.set("cartoon_transparency", 0.5, "all")

    def load_pdb(self, pdb_path: str):
        """
        Loads a protein structure from a PDB file and displays it as a cartoon.

        Args:
            pdb_path (str): The path to the PDB file containing the protein structure.
        """
        # Load protein and call it 'protein'
        self.cmd.load(pdb_path, "protein")
        self.cmd.hide("everything", "protein")
        self.cmd.show("cartoon", "protein")
        self.cmd.color("cyan", "protein")
        self.cmd.set("cartoon_transparency", 0.5, "protein")

    def load_molecule(self, mol: Chem.Mol):
        """
        Loads a molecule from an SDF file and centers the view on it.

        Args:
            sdf_path (str): The path to the SDF file containing the molecule.
        """

        # save tmp sdf
        sdf_path = "ligand_pymol.sdf"
        dm.to_sdf(mol, sdf_path)

        self.cmd.load(sdf_path)
        self.cmd.center("all")
        self.cmd.zoom("ligand_pymol", 2)
        self.cmd.orient("ligand_pymol", 40)

        # Color only ligand by element
        # self.cmd.color('gray', 'elem C')
        # self.cmd.color('red', 'elem O')
        # self.cmd.color('blue', 'elem N')

        # self.cmd.select("ligand_pymol")
        # self.cmd.color("byelement", "ligand_pymol")

        # remove tmp sdf
        os.remove(sdf_path)

    def get_view(self):
        return self.cmd.get_view()

    def set_view(self, view: list):
        self.cmd.set_view(view)

    def set_background_color(self):
        self.cmd.bg_color(self.background_color)
        self.cmd.set("ray_opaque_background", 1)
        self.cmd.set("depth_cue", 0)

    def save_image(self, filename: str):
        """
        Saves the current view to an image file.

        Args:
            filename (str): The name of the file to save the image as.
        """
        self.cmd.png(filename, ray=1)

    def save_session(self, filename: str):
        """Saves the current PyMOL session to a file."""
        self.cmd.save(filename)

    def display_in_ipython(self, width: int = 400, height: int = 400):
        """Returns the image of the current PyMOL view for display in IPython."""
        return self.cmd.ipython_image(width=width, height=height)

    def to_matplotlib(self, width: int = 400, height: int = 400):
        """Returns the image of the current PyMOL view for display in Matplotlib."""
        img = self.display_in_ipython(width, height)
        img = Image.open(BytesIO(img.data))
        return img
    
    def set_ray_trace_mode(self, mode: str):
        """
        Sets the ray trace mode for the PyMOL view.

        Args:
            mode (str): The mode to set the ray trace mode to. Can be "off", "on", or "adaptive".
        """
        self.cmd.set("ray_trace_mode", mode)

    def set_ray_trace_image(self, filename: str):
        """
        Sets the ray trace image for the PyMOL view.

        Args:
            filename (str): The path to the file to save the ray trace image to.
        """
        self.cmd.ray_trace_image(filename)
        

    def display_spheres(self, positions, radius=1, color="red"):
        """
        Display spheres in PyMOL.

        Args:
            positions (list): A list of positions for the spheres.
            radius (int, optional): The radius of the spheres. Defaults to 1.
            color (str, optional): The color of the spheres. Defaults to 'red'.

        Returns:
            None
        """
        # break into list
        if isinstance(radius, int) or isinstance(radius, float):
            radius = [radius] * len(positions)
        positions = [list(pos) for pos in positions]
        for i, pos in enumerate(positions):
            self.cmd.pseudoatom(
                object="sphere" + str(i),
                resn="sphere" + str(i),
                resi=i,
                chain="P",
                elem="PS",
                label="sphere" + str(i),
                vdw=radius[i],
                hetatm=1,
                color=color,
                pos=pos,
            )
        self.cmd.show(representation="spheres")
        self.cmd.center(selection="all")

if __name__ == "__main__":
    
    pymol = PyMOLInterface()
    pymol.load_pdb(EXAMPLES_PDB)
    pymol.display_spheres([(0, 0, 0), (1, 1, 1)], radius=1, color="red")
    pymol.save_image("spheres.png")