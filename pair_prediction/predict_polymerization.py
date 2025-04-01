from pair_prediction.tools.predict_bonds import PredictBonds
from pair_prediction.tools.constraint_validation import is_valid_core_pair
from pair_prediction.tools.utils import calc_distance


class PredictBondsCore(PredictBonds):
    """Implements the logic for prediction of core cross-linking.

    This class extends the PredictBonds class and modifies the behavior for core cross-linking prediction.

    Attributes:
        IDEAL_BOND_DIST (int): The equilibrium bonding distance in angstroms.
        BOND_DIST_STD (float): The standard deviation for the bond distance.
        query_radius (int): The radius within which to search for potential bonding partners.
    """

    IDEAL_BOND_DIST = 3.0  # angstroms
    BOND_DIST_STD = 2.5
    MAX_CONNECTIVITY = 8

    def __init__(self, pdb, reactive_input_file, QR, W):
        """Initializes the PredictBondsSur with a PDB object and a reactive input file.

        Args:
            pdb (PDB): The PDB object to predict bonds for.
            reactive_input_file (str): The input file for the specific reactive atom names and residues for the structure..
        """
        super().__init__(pdb, reactive_input_file, QR, W)
        self.query_radius = self.QR

    def calculate_bond_potential(self, atom1, atom2, atoms_dist):

        if not is_valid_core_pair(atom1, atom2):
            return 0
        
        if len(self.pdb.get_bonded_pairs_core()) > 0:
            all_bonded_chains = []
            # gets all the bonded chains (with duplicates) in the structure
            for pair in self.pdb.get_bonded_pairs_core():
                all_bonded_chains.append(pair.chain1)
                all_bonded_chains.append(pair.chain2)

            chain1_iso = all_bonded_chains.count(pair.chain1)
            chain2_iso = all_bonded_chains.count(pair.chain2)
            connectivity = chain1_iso + chain2_iso
            isolatedness = self.MAX_CONNECTIVITY - connectivity
        else:
            isolatedness = self.MAX_CONNECTIVITY

        #print(isolatedness)

        bond_pot = self.bond_potential(
            atoms_dist, self.IDEAL_BOND_DIST, self.BOND_DIST_STD**2, self.W, isolatedness
        )

        return bond_pot

    def add_pair_pdb(self, pair):
        self.pdb.add_bonded_pair(pair)
        self.pdb.add_core_bonded_pair(pair)


    def find_radicals(self, reactive_atoms_dict):
        """Add the radicals of all bonded chains to the radicals field of the PDB object if there are any"""

        bonded_chains = set()
        # gets all the bonded chains (without duplicates) in the structure
        for pair in self.pdb.get_bonded_pairs():
            bonded_chains.add(pair.chain1)
            bonded_chains.add(pair.chain2)

        for chain in bonded_chains:
            self.assign_radical(chain, reactive_atoms_dict)

    def assign_radical(self, chain, reactive_atoms_dict):
        """Add the radicals of the given chain to the radicals field of the PDB object if there are any"""
        # get atom names of reactive atoms
        atom_names = reactive_atoms_dict[chain.chain_type]
        # get the atom objects of the atom names
        atoms = []
        for atom_name in atom_names:
            atom = chain.atoms[atom_name]
            atoms.append(atom)
        # for each atom find its pair based on the distance btwn the atoms
        for cur_atom in atoms:
            if cur_atom.is_paired:
                continue
            # calculate the distances btwn the current atom and the rest of the atoms in the list
            atom_distances = []
            for sel_atom in atoms:
                if sel_atom.is_paired:
                    continue
                distance = calc_distance(cur_atom, sel_atom)
                atom_distances.append((distance, sel_atom))

            atom_distances.sort(key=lambda x: x[0])

            pair_atom_tuple = atom_distances[1]

            pair_atom = pair_atom_tuple[1]
            # determine which atom is a radical (n/a to atoms that are both bonded)
            if cur_atom.is_bonded_external and not pair_atom.is_bonded_external:
                self.pdb.add_radical(pair_atom)

            elif not cur_atom.is_bonded_external and pair_atom.is_bonded_external:
                self.pdb.add_radical(cur_atom)
            # if we already paired that atoms we don't need to consider them for future distance calculations
            cur_atom.is_paired = True
            pair_atom.is_paired = True

