import pytest
import numpy as np
from PIL import Image
from matplotlib.figure import Figure

from cchtools.plotting import create_comparison_plot
from cchtools.utils.images import smile_to_image, open_image
from cchtools.constants import EXAMPLES_PYMOL_IMAGE_1, EXAMPLES_PYMOL_IMAGE_2

@pytest.fixture
def sample_data():
    img_1 = open_image(EXAMPLES_PYMOL_IMAGE_1)
    img_2 = open_image(EXAMPLES_PYMOL_IMAGE_2)
    smiles = ["CCO", "CCN", "CCCl"]
    smile_imgs = [smile_to_image(smile) for smile in smiles]
    texts = ["Text 1", "Text 2", "Text 3"]
    return img_1, img_2, smile_imgs, texts

def test_create_comparison_plot(sample_data):
    img_1, img_2, smile_imgs, texts = sample_data
    fig = create_comparison_plot(
        images=[img_1, img_1, img_2],
        smile_imgs=smile_imgs,
        texts=texts,
        titles=["Generated", "Reference"]
    )
    
    assert isinstance(fig, Figure)
    assert len(fig.axes) == 6  # 2 rows * 3 columns

def test_create_comparison_plot_no_reference(sample_data):
    img_1, _, smile_imgs, texts = sample_data
    fig = create_comparison_plot(
        images=[img_1, img_1],
        smile_imgs=smile_imgs[:2],
        texts=texts[:2],
        reference_img=False
    )
    
    assert isinstance(fig, Figure)
    assert len(fig.axes) == 4  # 2 rows * 2 columns

def test_create_comparison_plot_custom_figsize(sample_data):
    img_1, img_2, smile_imgs, texts = sample_data
    custom_figsize = (12, 6)
    fig = create_comparison_plot(
        images=[img_1, img_1, img_2],
        smile_imgs=smile_imgs,
        texts=texts,
        figsize=custom_figsize
    )
    
    assert isinstance(fig, Figure)

def test_create_comparison_plot_custom_corrections(sample_data):
    img_1, img_2, smile_imgs, texts = sample_data
    fig = create_comparison_plot(
        images=[img_1, img_1, img_2],
        smile_imgs=smile_imgs,
        texts=texts,
        text_correction=-0.6,
        header_correction=1.2,
        mol_correction=0.2
    )
    
    assert isinstance(fig, Figure)

def test_create_comparison_plot_line_width(sample_data):
    img_1, img_2, smile_imgs, texts = sample_data
    custom_line_width = 5.0
    fig = create_comparison_plot(
        images=[img_1, img_1, img_2],
        smile_imgs=smile_imgs,
        texts=texts,
        line_width=custom_line_width
    )
    
    assert isinstance(fig, Figure)
    for ax in fig.axes[:3]:  # Check only top row axes
        for spine in ax.spines.values():
            assert spine.get_linewidth() == custom_line_width
