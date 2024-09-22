import os
from io import BytesIO

import datamol as dm
import rdkit.Chem as Chem
from PIL import Image
from pymol import cmd


class PyMOLInterface:
    """
    A class to interface with PyMOL for molecular visualization.
    """

    def __init__(self, background_color: str = "white", view_size: tuple = (400, 400)):
        """
        Initializes the PyMOL interface.
        """
        self.cmd: cmd = cmd
        self.background_color = background_color

        # set view size
        self.cmd.viewport(*view_size)

    def delete_all(self):
        self.cmd.delete("all")

    def fetch_protein(
        self,
        pdb_code: str,
        show_as: str = "cartoon",  # cartoon, surface
        color: str = "cyan",
        transparency: float = 0.5,
    ):
        """
        Fetches a protein structure by its PDB code and displays it as a cartoon.

        Args:
            pdb_code (str): The PDB code of the protein to fetch.
        """
        print(f"Fetching protein {pdb_code}")

        self.cmd.fetch(pdb_code)
        self.cmd.hide("everything", "all")

        if show_as == "cartoon":
            self.cmd.show("cartoon", "all")
        elif show_as == "surface":
            self.cmd.show("surface", "all")
        else:
            raise ValueError(f"Invalid show_as value: {show_as}")

        self.cmd.color(color, "all")
        self.cmd.set("transparency", transparency, "all")

    def load_pdb(
        self,
        pdb_path: str,
        show_as: str = "cartoon",  # cartoon, surface
        color: str = "cyan",
        transparency: float = 0.5,
    ):
        """
        Loads a protein structure from a PDB file and displays it.

        Args:
            pdb_path (str): The path to the PDB file containing the protein structure.
            show_as (str): The representation style for the protein. Either "cartoon" or "surface".
            color (str): The color to apply to the protein.
            transparency (float): The transparency level of the protein representation.
        """
        # Load protein and call it 'protein'
        self.cmd.load(pdb_path, "protein")
        self.cmd.hide("everything", "protein")

        if show_as == "cartoon":
            self.cmd.show("cartoon", "protein")
            self.cmd.set("cartoon_transparency", transparency, "protein")
        elif show_as == "surface":
            self.cmd.show("surface", "protein")
            # self.cmd.set("surface_transparency", transparency, "protein")
        else:
            raise ValueError(f"Invalid show_as value: {show_as}")

        self.cmd.color(color, "protein")

    def set_protein_color(self, color: str):
        """
        Sets the color of the protein.

        Args:
            color (str): The color to set the protein to.
        """
        self.cmd.color(color, "protein")

    def load_molecule(self, mol: Chem.Mol, color: str = "green", by_element: bool = True):
        """
        Loads a molecule from an SDF file and centers the view on it.

        Args:
            sdf_path (str): The path to the SDF file containing the molecule.
        """

        # save tmp sdf
        sdf_path = "ligand.sdf"
        dm.to_sdf(mol, sdf_path)

        self.cmd.load(sdf_path)
        self.cmd.center("all")
        self.cmd.zoom("ligand", 2)
        self.cmd.orient("ligand", 40)

        # Color by element
        if by_element:
            self.cmd.color(color, "elem C")
            self.cmd.color("red", "elem O")
            self.cmd.color("blue", "elem N")

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

    def save_image(self, filename: str, dpi: int = 300, width: int = 400, height: int = 400):
        """
        Saves the current view to an image file.

        Args:
            filename (str): The name of the file to save the image as.
        """
        self.cmd.png(filename, ray=1, dpi=dpi, width=width, height=height)

    def save_image_as_svg(self, filename: str):
        """
        Saves the current view to an SVG file.

        Args:
            filename (str): The name of the file to save the image as.
        """
        self.cmd.svg(filename, ray=1)

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

    def set_ray_trace_image(self, width: int = 400, height: int = 400):
        """
        Sets the ray trace image for the PyMOL view.

        Args:
            filename (str): The path to the file to save the ray trace image to.
        """
        self.cmd.ray(width=width, height=height)

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

    def show_residue_sidechains(self, residue_ids: list[str], sidechain_helper: bool = True):
        """
        Shows only the sidechains coming out from the cartoon representation for the specified residues.

        Args:
            residue_ids (list[int]): A list of residue IDs to show the sidechains of.
        """
        residue_str = "+".join(str(id) for id in residue_ids)
        self.cmd.show("sticks", f"(resi {residue_str}) and not name c+n+o")
        self.cmd.hide("sticks", f"(resi {residue_str}) and name c+n+o")

        if sidechain_helper:
            self.cmd.set("cartoon_side_chain_helper", 1)

    def remove_all_hydrogens(self):
        """
        Removes all hydrogens from the protein.
        """
        self.cmd.remove("hydrogen")

    def get_interacting_residues(
        self, ligand_name: str = "ligand", distance: float = 4.0, include_details: bool = False
    ):
        """
        Get the protein residues interacting with the ligand within a specified distance.

        Args:
            ligand_name (str): The name of the ligand object in PyMOL. Defaults to "ligand".
            distance (float): The distance threshold for interaction. Defaults to 4.0 Angstroms.
            include_details (bool): If True, return detailed information about each interaction.

        Returns:
            list: A list of residue identifiers or detailed interaction information.
        """
        # Ensure the ligand object exists
        if ligand_name not in self.cmd.get_names("objects"):
            raise ValueError(f"No object named '{ligand_name}' found in the PyMOL session.")

        # Select interacting residues, excluding waters and focusing on protein
        selection_name = "interacting_residues"
        self.cmd.select(
            selection_name,
            f"(byres (protein and not hetatm) within {distance} of {ligand_name}) and not solvent",
        )

        # Get the number of selected atoms
        n_selected = self.cmd.count_atoms(selection_name)
        print(f"Number of atoms selected: {n_selected}")

        if include_details:
            # Get detailed information about the interactions
            interactions = []
            for atom in self.cmd.get_model(selection_name).atom:
                res_info = f"{atom.resn}{atom.resi}"
                dist = self.cmd.distance("tmp_dist", f"{ligand_name}", f"id {atom.id}")
                interactions.append((res_info, atom.name, dist))
            self.cmd.delete("tmp_dist")

            # Sort by distance and remove duplicates
            interactions.sort(key=lambda x: x[2])
            unique_interactions = []
            seen_residues = set()
            for res, atom, dist in interactions:
                if res not in seen_residues:
                    unique_interactions.append((res, atom, dist))
                    seen_residues.add(res)

            # Clean up the selection
            self.cmd.delete(selection_name)
            return unique_interactions
        else:
            residues = []
            # Get the residue identifiers
            self.cmd.iterate(
                selection_name, "residues.append(int(resi))", space={"residues": residues}
            )

            # Clean up the selection
            self.cmd.delete(selection_name)
            return residues


if __name__ == "__main__":
    from cchtools.constants import EXAMPLES_PDB, EXAMPLES_SDF

    view = (
        0.328213751,
        0.371481866,
        -0.868490517,
        0.564871073,
        0.659728706,
        0.495660394,
        0.757095039,
        -0.653266013,
        0.006699756,
        0.000289902,
        -0.000215106,
        -35.168968201,
        -21.841312408,
        18.973146439,
        -27.182886124,
        -304.886566162,
        375.222930908,
        -20.000000000,
    )

    pymol_interface = PyMOLInterface()
    pymol_interface.set_background_color()
    pymol_interface.load_pdb(EXAMPLES_PDB, transparency=0.5, show_as="cartoon")
    # pymol_interface.display_spheres([(0, 0, 0), (1, 1, 1)], radius=1, color="red")

    # load ligand
    ligand = dm.read_sdf(EXAMPLES_SDF)[0]
    pymol_interface.load_molecule(ligand, color="green")

    pymol_interface.set_protein_color("cyan")

    # Set view
    pymol_interface.set_view(view)

    # Set ray trace mode
    pymol_interface.set_ray_trace_mode(3)

    residues = pymol_interface.get_interacting_residues(distance=3.0)

    # show sidechains
    pymol_interface.show_residue_sidechains(residues)
    pymol_interface.remove_all_hydrogens()
    # pymol_interface.set_ray_trace_image(width=1000, height=1000)
    pymol_interface.save_image("spheres.png", dpi=20, width=500, height=500)
