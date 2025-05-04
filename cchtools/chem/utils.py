from typing import Callable
from rdkit.Chem import Mol
from functools import wraps


def preserve_properties(func: Callable[[Mol, ...], Mol]) -> Callable[[Mol, ...], Mol]:
    """
    Decorator to copy RDKit molecule properties from the original molecule
    to the new molecule returned by the function.

    This ensures that any properties set via mol.SetProp(...) are
    preserved when the molecule is modified or converted.
    """
    @wraps(func)
    def wrapped(*args, **kwargs) -> Mol:
        # Find the first RDKit Mol in args or kwargs
        if "mol" in kwargs:
            mol = kwargs["mol"]
        else:
            mol = next(arg for arg in args if isinstance(arg, Mol))

        # Save all existing properties
        props = mol.GetPropsAsDict()

        # Call the original function
        new_mol = func(*args, **kwargs)

        # Restore saved properties
        for name, value in props.items():
            # ensure the value is a string for RDKit.SetProp
            new_mol.SetProp(name, str(value))

        return new_mol

    return wrapped