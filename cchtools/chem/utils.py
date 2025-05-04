import functools
import io
import logging
from contextlib import redirect_stderr
from rdkit import RDLogger

from typing import Callable
from rdkit.Chem import Mol
from functools import wraps

logger = logging.getLogger(__name__)


def silence_rdkit(func):
    """
    Decorator: run *func* while RDKit’s Python-level log messages and
    C++ stderr warnings are suppressed.

    Example
    -------
    >>> @silence_rdkit
    ... def make_mol(smiles):
    ...     from rdkit import Chem
    ...     return Chem.MolFromSmiles(smiles)
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # 1. Quiet RDKit’s Python logger
        logger = RDLogger.logger()
        RDLogger.DisableLog('rdApp.*')             # mute everything

        # 2. Catch C++ warnings that appear on stderr
        fake_err = io.StringIO()
        with redirect_stderr(fake_err):
            try:
                return func(*args, **kwargs)
            finally:
                # 3. Restore the logging state whatever happens
                RDLogger.EnableLog('rdApp.*')
    return wrapper


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