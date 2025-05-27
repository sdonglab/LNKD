### Author: Bradford Derby ###

class Pair:
    """Represents a potential pairing of atoms in a structure.

    Attributes:
        atom1 (Atom): The first atom in the pair.
        atom2 (Atom): The second atom in the pair.
        chain1 (Chain): The chain that the first atom belongs to.
        chain2 (Chain): The chain that the second atom belongs to.
        distance (float): The distance between the two atoms.
        probability (float): The probability of this pair being a valid pairing, initially set to 0.
    """

    def __init__(self, atom1, atom2, distance):
        """Initializes a Pair object.

        Args:
            atom1 (Atom): The first atom in the pair.
            atom2 (Atom): The second atom in the pair.
            distance (float): The distance between the two atoms.
        """
        self.atom1 = atom1
        self.atom2 = atom2
        self.chain1 = atom1.chain
        self.chain2 = atom2.chain
        self.distance = distance
        self.probability = 0

    def set_probability(self, probability):
        """Sets the probability attribute of the Pair object.

        Args:
            probability (float): The probability of this pair being a valid pairing.
        """
        self.probability = probability

    # override of less than for heap (but inverted for max heap)
    def __lt__(self, other):
        return self.probability > other.probability

    def __hash__(self):
        return hash(self.distance)

    def __eq__(self, other):
        return self.distance == other.distance

    # say you have pair A-B, A-C, D-C. Then if A-B is current pair and [A-C, D-C] is the potential_pair_list, then this function will return [A-C] sinc A is common in both A-B and A-C
    def get_chain_branching_pairs(self, potential_pairs):
        """Finds the pairs that share a chain in common with this pair.

        Args:
            potential_pairs (list[Pair]): A list of potential pairs to check for common chains.

        Returns:
            list[Pair]: A list of pairs that share a chain in common with this pair.
        """
        linking_pairs = []
        for pair in potential_pairs:
            if self.is_same_chain(pair):
                linking_pairs.append(pair)
        return linking_pairs

    def is_same_chain(self, other_pair):
        """Determines if this pair's atoms are part of the same chain as the given pair's atoms.

        Args:
            other_pair (Pair): The pair to compare with.

        Returns:
            bool: True if this pair's atoms are part of the same chain as the given pair's atoms, False otherwise.
        """
        chain_id_a = self.atom1.res_seq
        chain_id_b = self.atom2.res_seq
        chain_id_x = other_pair.atom1.res_seq
        chain_id_y = other_pair.atom2.res_seq

        return (chain_id_a in [chain_id_x, chain_id_y]) or (
            chain_id_b in [chain_id_x, chain_id_y]
        )
    
    def get_atom_branching_pairs(self, potential_pairs):
        """Finds the pairs that share an atom in common with this pair.

        Args:
            potential_pairs (list[Pair]): A list of potential pairs to check for common atoms.

        Returns:
            list[Pair]: A list of pairs that share an atom in common with this pair.
        """
        linking_pairs = []
        for pair in potential_pairs:
            if self.are_atoms_matching(pair):
                linking_pairs.append(pair)
        return linking_pairs

    def are_atoms_matching(self, other):
        """Determines if the atoms of this pair overlap (have a common atom) with the given pair.

        For instance, pair a1,a2 and a1,a3 would have a1 as an overlapping or common atom.

        Args:
            other (Pair): The pair to compare with.

        Returns:
            bool: True if this pair's atoms overlap with the given pair's atoms, False otherwise.
        """
        return self.atom1 in [other.atom1, other.atom2] or self.atom2 in [
            other.atom1,
            other.atom2,
        ]

    def bond_pair(self):
        """Updates the properties of the atoms in the pair to be bonded.

        This method updates the chains of the atoms to be bonded to each other and sets the atoms'
        'is_bonded_external' attribute to True.
        """
        self.chain1.add_bonded_chain(self.chain2)
        self.chain2.add_bonded_chain(self.chain1)

        self.atom1.is_bonded_external = True
        self.atom2.is_bonded_external = True
