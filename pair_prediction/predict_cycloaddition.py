###  Author: Bradford Derby   ###
### Modified by: Emma Stevens ###

from enums import SchemeOneBE, SchemeOneBN, SchemeTwoBN
from tools.pair import Pair
from tools.predict_bonds import PredictBonds
from tools.constraint_validation import is_valid_LN2_surface_pair


class PredictBondsSur(PredictBonds):
    """Implements the logic for prediction of surface cross-linking.

    This class extends the PredictBonds class and modifies the behavior for surface cross-linking prediction.

    Attributes:
        IDEAL_BOND_DIST (int): The equilibrium bonding distance in angstroms.
        BOND_DIST_STD (float): The standard deviation for the bond distance.
        query_radius (int): The radius within which to search for potential bonding partners.
    """

    IDEAL_BOND_DIST = 3.0  # angstroms
    BOND_DIST_STD = 2.5
    MAX_CONNECTIVITY = 8

    def __init__(self, pdb, reactive_input_file, query_radius, weight):
        """Initializes the PredictBondsSur with a PDB object and a reactive input file.

        Args:
            pdb (PDB): The PDB object to predict bonds for.
            reactive_input_file (str): The input file for the specific reactive atom names and residues for the structure.
            QR (float): The query radius for finding nearby reactive atoms.
            W (float): The weight for the degree of isolation.
        """
        super().__init__(pdb, reactive_input_file, query_radius, weight)

    def calculate_bond_potential(self, atom1, atom2, atoms_dist):

        if not is_valid_LN2_surface_pair(atom1, atom2):
            return 0
        
        if len(self.pdb.get_bonded_pairs_surf()) > 0:
            all_bonded_chains = []
            # gets all the bonded chains (with duplicates) in the structure
            for pair in self.pdb.get_bonded_pairs_surf():
                all_bonded_chains.append(pair.chain1)
                all_bonded_chains.append(pair.chain2)

            chain1_iso = all_bonded_chains.count(pair.chain1)
            chain2_iso = all_bonded_chains.count(pair.chain2)
            connectivity = chain1_iso + chain2_iso
            isolatedness = self.MAX_CONNECTIVITY - connectivity
        else:
            isolatedness = self.MAX_CONNECTIVITY

        bond_pot = self.bond_potential(
            atoms_dist, self.IDEAL_BOND_DIST, self.BOND_DIST_STD**2, self.weight, isolatedness, self.MAX_CONNECTIVITY
        )
        return bond_pot


    def add_pair_pdb(self, pair):
        pair1, pair2 = self.get_linking_pairs(pair)
        self.pdb.add_bonded_pair(pair1)
        self.pdb.add_bonded_pair(pair2)

        self.pdb.add_surface_bonded_pair(pair1)
        self.pdb.add_surface_bonded_pair(pair2)

    def find_radicals(self, reactive_atoms_dict):
        pass

    def get_linking_pairs(self, pair):
        atom_a = pair.atom1
        atom_b = pair.atom2

        atom_b_pair = self.find_linking_pair(atom_a)
        atom_a_pair = self.find_linking_pair(atom_b)

        pair1 = Pair(atom_a, atom_a_pair, 0)
        pair2 = Pair(atom_b, atom_b_pair, 0)

        return pair1, pair2

    def find_linking_pair(self, atom):

        if atom.element == SchemeOneBE.CARBON.value:
            atom_xyz = atom.name[1]
            corresponding_pair = self.lookup_corresponding_atom(
                atom, atom_xyz, SchemeOneBN.CARBON.value
            )
            return corresponding_pair
        elif atom.element == SchemeOneBE.NITROGEN.value:
            if atom.name[1] == "6":
                corresponding_pair = self.lookup_corresponding_atom(
                    atom, "", SchemeTwoBN.NITROGEN6.value
                )
            elif atom.name[1] == "3":
                corresponding_pair = self.lookup_corresponding_atom(
                    atom, "", SchemeTwoBN.NITROGEN3.value
                )
            return corresponding_pair

    def lookup_corresponding_atom(self, atom, atom_xyz, name_number_lookup):

        atom_string = f"{atom.element}{atom_xyz}{str(name_number_lookup)}"
        atom_chain = atom.chain
        atoms_in_chain = atom_chain.atoms

        # find the pairing atom for the current chain's atom so for CX3 it would return CX2s
        return atoms_in_chain[atom_string]

