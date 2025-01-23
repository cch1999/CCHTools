from typing import Tuple, Union
import logging
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem, rdChemReactions

from itertools import chain

RDLogger.DisableLog("rdApp.*")


def clean_molecule(mol):
    """Clean up a molecule by removing Hs and sanitizing."""
    if mol is None:
        return None
    try:
        mol = Chem.RemoveHs(mol)
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return mol

def run_reaction_with_fallback(rxn, reactants):
    """Run reaction with fallback error handling."""
    try:
        products = rxn.RunReactants(reactants)
        if len(products) == 0:
            raise ValueError("Reaction did not yield any products")
        return products
    except Exception as e:
        logging.warning(f"Reaction failed: {e}")
        return []

class Reaction:
    def __init__(self, template: str):
        self.template = template
        self.rxn = self.__init_reaction()
        # Extract reactants, agents, products

    def __init_reaction(self) -> Chem.rdChemReactions.ChemicalReaction:
        """Initializes a reaction by converting the SMARTS-pattern to an `rdkit` object."""
        rxn = AllChem.ReactionFromSmarts(self.template)
        rdChemReactions.ChemicalReaction.Initialize(rxn)
        return rxn
    
    def __repr__(self):
        return f"Reaction(template={self.template})"

    def reverse_template(self):
        """Reverses a reaction template and returns an initialized, reversed reaction object."""
        rxn = AllChem.ChemicalReaction()
        for i in range(self.rxn.GetNumReactantTemplates()):
            rxn.AddProductTemplate(self.rxn.GetReactantTemplate(i))
        for i in range(self.rxn.GetNumProductTemplates()):
            rxn.AddReactantTemplate(self.rxn.GetProductTemplate(i))
        rdChemReactions.ChemicalReaction.Initialize(rxn)
        return rxn

    def is_reactant(self, mol: Chem.Mol, rxn: rdChemReactions = None) -> bool:
        """Checks if a molecule is a reactant for the reaction."""
        if rxn is None:
            rxn = self.rxn
        return self.rxn.IsMoleculeReactant(mol)

    def is_product(self, mol: Chem.Mol, rxn: rdChemReactions = None) -> bool:
        """Checks if a molecule is a reactant for the reaction."""
        if rxn is None:
            rxn = self.rxn
        return self.rxn.IsMoleculeProduct(mol)

    def run_reactants(
        self, reactants: Tuple[Union[Chem.Mol, str, None]], rxn: rdChemReactions = None,
    ) -> Union[Chem.Mol, None]:
        """Runs the reaction on a set of reactants and returns the products."""
        if rxn is None:
            rxn = self.rxn

        # Run the reaction
        products = run_reaction_with_fallback(rxn, reactants)

        if not products:
            return None
        
        # Clean the products
        prods_out = [
            clean_molecule(prod)
            for prod in chain.from_iterable(products)
            if clean_molecule(prod) is not None
        ]
        return prods_out

    def run_reverse_reactants(
        self, product: Tuple[Chem.Mol], rxn: rdChemReactions = None,
    ) -> Union[Chem.Mol, None]:
        """Runs the reverse reaction on a product to return the reactants."""
        if rxn is None:
            rxn = self.reverse_template()

        return self.run_reactants(product, rxn)
