import matplotlib.pyplot as plt
import numpy as np
from typing import List, Tuple, Optional

from cchtools.utils.images import smile_to_image
from cchtools.utils.cif import cif_to_rdkit
from cchtools.rdkit_utils import calculate_basic_metrics

def create_comparison_plot(
    images: List[np.ndarray],
    smile_imgs: List[np.ndarray],
    texts: List[str],
    titles: Optional[List[str]] = None,
    text_correction: float = -0.55,
    header_correction: float = 1.17,
    mol_correction: float = 0.16,
    reference_img: bool = True,
    figsize: Optional[Tuple[float, float]] = None,
    line_width: float = 3.5
) -> plt.Figure:
    """
    Create a comparison plot with images, molecule representations, and text.
    
    Args:
        images (List[np.ndarray]): List of images to display in the top row.
        smile_imgs (List[np.ndarray]): List of molecule images to display in the bottom row.
        texts (List[str]): List of text strings to display below each column.
        titles (Optional[List[str]]): List of two strings for left and right titles. Defaults to None.
        text_correction (float): Correction factor for the text position. Defaults to -0.55.
        header_correction (float): Correction factor for the header position. Defaults to 1.17.
        mol_correction (float): Correction factor for the molecule position. Defaults to 0.16.
        reference_img (bool): Whether to include a reference image. Defaults to True.
        figsize (Optional[Tuple[float, float]]): Figure size. If None, calculated based on number of images.
        line_width (float): Width of the subplot borders. Defaults to 3.5.
    
    Returns:
        plt.Figure: The created figure.
    """
    
    num_images = len(images)
    if figsize is None:
        figsize = (num_images * 3.3333333, 8)
    
    # Create a figure and subplot
    fig, axs = plt.subplots(2, num_images, figsize=figsize, gridspec_kw={'wspace': 0, 'hspace': 0.00})

    for i in range(num_images):
        # Main plot
        axs[0, i].imshow(images[i])
        axs[0, i].set_xticks([])
        axs[0, i].set_yticks([])
        for spine in axs[0, i].spines.values():
            spine.set_linewidth(line_width)

        # Molecule image plot
        axs[1, i].imshow(smile_imgs[i])
        axs[1, i].set_xticks([])
        axs[1, i].set_yticks([])
        axs[1, i].axis('off')

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.0, hspace=0.00)

    # Move the last subplot to the right if reference image is included
    if reference_img:
        for row in [0, 1]:
            last_pos = axs[row, -1].get_position()
            axs[row, -1].set_position([last_pos.x0 + 0.03, last_pos.y0, last_pos.width, last_pos.height])

    # Adjust position of bottom row
    for ax in axs[1]:
        last_pos = ax.get_position()
        ax.set_position([last_pos.x0, last_pos.y0 + mol_correction, last_pos.width, last_pos.height])

    # Add text boxes below each subplot
    for i in range(num_images):
        axs[0, i].text(0.5, text_correction, texts[i], transform=axs[0, i].transAxes, 
                ha='center', va='top', fontsize=20, color='black', 
                bbox=dict(facecolor='white', edgecolor='none', alpha=1, pad=3))

    # Add titles if provided
    if titles and len(titles) == 2:
        axs[0, 0].text(0.03, header_correction, titles[0], fontsize=25, ha='left', va='top', transform=axs[0, 0].transAxes)
        if reference_img:
            axs[0, -1].text(0.03, header_correction, titles[1], fontsize=25, ha='left', va='top', transform=axs[0, -1].transAxes)

    return fig


if __name__ == "__main__":

    from cchtools.utils.images import open_image
    from cchtools.constants import EXAMPLES_PYMOL_IMAGE_1, EXAMPLES_PYMOL_IMAGE_2
    from cchtools.rdkit_utils import calculate_basic_metrics
    from rdkit import Chem
    img_1 = open_image(EXAMPLES_PYMOL_IMAGE_1)
    img_2 = open_image(EXAMPLES_PYMOL_IMAGE_2)

    smiles = ["OCCc1c(C)[n+](cs1)Cc2cnc(C)nc2N", "OCCc1c(C)[n+](cs1)Cc2cnc(C)nc2N", "CCO"]
    smile_imgs = [smile_to_image(smile) for smile in smiles]

    mols = [Chem.MolFromSmiles(smile) for smile in smiles]
    metrics = calculate_basic_metrics(mols, ref_mol=mols[-1])

    texts = [f"QED: {metric['QED']:.2f}\nSA: {metric['SA']:.2f}" for metric in metrics]

    fig = create_comparison_plot(images=[img_1, img_1, img_2],
                                 smile_imgs=smile_imgs, 
                                 texts=texts,
                                 reference_img=True,
                                 header_correction=1.1,
                                 titles=["Generated", "Reference"])
    
    plt.show()