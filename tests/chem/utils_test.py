import sys
import pytest
from cchtools.chem.utils import silence_rdkit, preserve_properties
from rdkit.Chem import MolFromSmiles, MolToSmiles

def test_silence_rdkit_suppresses_stderr(capsys):
    @silence_rdkit
    def noisy():
        # write something to stderr
        print("error message", file=sys.stderr)
        return 123
    result = noisy()
    captured = capsys.readouterr()
    assert result == 123
    # nothing should be emitted to stderr
    assert captured.err == ""

def test_without_silence_rdkit_emits_stderr(capsys):
    def noisy():
        print("visible error", file=sys.stderr)
        return 456
    result = noisy()
    captured = capsys.readouterr()
    assert result == 456
    assert captured.err.strip() == "visible error"

def test_preserve_properties_roundtrip_args_and_kwargs():
    original = MolFromSmiles("CCO")
    original.SetProp("a", "1")
    original.SetProp("b", "2")

    @preserve_properties
    def roundtrip(mol):
        # lose properties by converting to SMILES and back
        smi = MolToSmiles(mol)
        return MolFromSmiles(smi)

    # positional call
    new1 = roundtrip(original)
    assert new1.GetNumAtoms() == original.GetNumAtoms()
    assert new1.GetPropsAsDict() == original.GetPropsAsDict()

    # keyword call
    new2 = roundtrip(mol=original)
    assert new2.GetPropsAsDict() == original.GetPropsAsDict()

def test_preserve_properties_with_multiple_mols():
    mol1 = MolFromSmiles("CCC")
    mol1.SetProp("key", "val")
    mol2 = MolFromSmiles("CC")

    @preserve_properties
    def return_second(mol, other):
        # ignores the first and returns the second
        return other

    new = return_second(mol1, mol2)
    # new is mol2 but should have inherited props from mol1
    assert new is mol2
    assert new.GetPropsAsDict() == {"key": "val"}

def test_decorators_preserve_function_metadata():
    def custom_func(x):
        """custom doc"""
        return x

    # silence_rdkit should wrap but preserve metadata
    wrapped1 = silence_rdkit(custom_func)
    assert hasattr(wrapped1, "__wrapped__")
    assert wrapped1.__wrapped__ is custom_func
    assert wrapped1.__name__ == custom_func.__name__
    assert wrapped1.__doc__ == custom_func.__doc__

    # preserve_properties should do the same
    wrapped2 = preserve_properties(custom_func)
    assert hasattr(wrapped2, "__wrapped__")
    assert wrapped2.__wrapped__ is custom_func
    assert wrapped2.__name__ == custom_func.__name__
    assert wrapped2.__doc__ == custom_func.__doc__
