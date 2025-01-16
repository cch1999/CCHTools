from typing import Tuple

import datamol as dm
import numpy as np
import rdkit.Chem as Chem
from PIL import Image


def open_image(path: str) -> Image.Image:
    """
    Open an image from a file path.
    """
    return Image.open(path)


def smile_to_image(
    smile: str, size: Tuple[int, int] = (400, 400), crop: bool = True
) -> Image.Image:
    """
    Convert a SMILES string to a cropped image of the molecule with transparent background.

    Args:
        smile (str): SMILES representation of the molecule
        size (tuple): Size of the image before cropping (width, height)
        crop (bool): Whether to crop the image to the molecule

    Returns:
        PIL.Image.Image: Cropped image of the molecule with transparent background
    """
    # Convert SMILES to RDKit mol object
    mol = dm.to_mol(smile)

    # Draw the molecule
    img = Chem.Draw.MolToImage(mol, size=size)

    # Convert to RGB mode if not already
    img = img.convert("RGB")

    # Convert image to numpy array
    img_array = np.array(img)

    # Find non-white pixels
    y, x = np.where(np.any(img_array < 255, axis=2))

    padding = 10

    # Add padding and crop
    if len(x) > 0 and len(y) > 0 and crop:
        img = img.crop((x.min() - padding, y.min() - padding, x.max() + padding, y.max() + padding))

    # Make background transparent
    img = img.convert("RGBA")
    datas = img.getdata()

    new_data = []
    for item in datas:
        if item[0] == 255 and item[1] == 255 and item[2] == 255:
            new_data.append((255, 255, 255, 0))
        else:
            new_data.append(item)

    img.putdata(new_data)

    return img


if __name__ == "__main__":
    # Example usage:
    smiles = [
        "CC(=O)NCCC1=CNc2c1cc(OC)cc2CC(=O)NCCc1c[nH]c2ccc(OC)cc12",
        "CC1=C(C(=O)C[C@@H]1OC(=O)[C@@H]2[C@H](C2(C)C)/C=C(\C)/C(=O)OC)C/C=C\C=C",
        "OCCc1c(C)[n+](cs1)Cc2cnc(C)nc2N",
    ]

    smile_imgs = [smile_to_image(smile) for smile in smiles]

    # Display an image (for testing purposes)
    smile_imgs[1].show()
