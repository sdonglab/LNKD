### Author: Bradford Derby ###

class Chain:
    """Represents a chain/molecule in a PDB file or structure.

    Attributes:
        chain_id (int): The identifier of the chain.
        atoms (dict): A dictionary of atoms that belong to the chain where the key is the atom name and the value is the atom object.
        bonded_chains (list): A list of chains that are bonded to this chain.
        chain_type (str): The type of the chain.
    """

    def __init__(self, chain_id, chain_type):
        """Initializes a Chain object.

        Args:
            chain_id (int): The identifier of the chain.
            chain_type (str): The type of the chain.
        """
        self.chain_id = chain_id
        self.atoms = {}
        self.bonded_chains = []
        self.chain_type = chain_type

    def __hash__(self):
        return hash(self.chain_id)

    def __eq__(self, other):
        """Checks if two chains are equal based on chain_id"""
        return self.chain_id == other.chain_id

    def add_atom(self, atom):
        """Adds an atom to the atoms dict of this chain.

        Args:
            atom (Atom): The atom to be added to the chain.
        """
        self.atoms[atom.name] = atom

    def add_bonded_chain(self, chain):
        """Adds a chain to the list of bonded chains.

        Args:
            chain (Chain): The chain to be added to the list of bonded chains.
        """
        self.bonded_chains.append(chain)
