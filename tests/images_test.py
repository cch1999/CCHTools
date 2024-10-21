import pytest
from PIL import Image
import numpy as np
from cchtools.utils.images import smile_to_image

@pytest.fixture
def sample_smiles():
    return [
        'CC(=O)NCCC1=CNc2c1cc(OC)cc2CC(=O)NCCc1c[nH]c2ccc(OC)cc12',
        r'CC1=C(C(=O)C[C@@H]1OC(=O)[C@@H]2[C@H](C2(C)C)/C=C(\C)/C(=O)OC)C/C=C\C=C',
        'OCCc1c(C)[n+](cs1)Cc2cnc(C)nc2N'
    ]

def test_smile_to_image_output_type(sample_smiles):
    for smile in sample_smiles:
        img = smile_to_image(smile)
        assert isinstance(img, Image.Image)

def test_smile_to_image_size():
    smile = 'CC(=O)O'  # Acetate
    img = smile_to_image(smile, size=(300, 300), crop=False)
    assert img.size == (300, 300)

def test_smile_to_image_transparency():
    smile = 'CC(=O)O'  # Acetate
    img = smile_to_image(smile)
    assert img.mode == 'RGBA'
    
    # Check if there are transparent pixels
    img_array = np.array(img)
    assert np.any(img_array[:, :, 3] == 0)

def test_smile_to_image_cropping():
    smile = 'OCCc1c(C)[n+](cs1)Cc2cnc(C)nc2N'  # Pyrazine
    img_cropped = smile_to_image(smile, crop=True)
    img_uncropped = smile_to_image(smile, crop=False)
    
    assert img_cropped.size != img_uncropped.size
    assert img_cropped.size[0] <= img_uncropped.size[0]
    assert img_cropped.size[1] <= img_uncropped.size[1]
