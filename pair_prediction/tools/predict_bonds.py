# Super class for bond prediction containing common methods and attributes
from sklearn.neighbors import KDTree
from abc import ABC, abstractmethod
from tools.pair import Pair
from tools.utils import atoms_to_coords, remove_1st_and_chain
import numpy as np
import heapq


class PredictBonds(ABC):
    """Functions as a super class for different implementations of crosslinking prediction.

    Attributes:
        pdb (str): The pdb file of the structure.
        atoms (list): The atoms of the structure.
        probability_heap (list): The heap structure used to store potential pairs.
        potential_pairs (list): A list to store potential pairs of atoms for bonding.
        reactive_atoms (list): The reactive atoms in the structure denoted by the user.
        query_radius (float): The float value used for the spatial query of reactive atoms.
        reactive_input_file (str): The input file for the specific reactive atom names and residues for the structure.
    """

    def __init__(self, pdb: str, reactive_input_file: str, query_radius: float, weight: float):
        """Inits PredictBonds with pdb and reactive_input_file.

        Args:
            pdb (str): The pdb file of the structure.
            reactive_input_file (str): The input file for the specific reactive atom names and residues for the structure.
        """
        self.pdb = pdb
        self.atoms = pdb.get_atoms()
        self.probability_heap = []
        self.potential_pairs = []
        self.reactive_atoms = []
        self.query_radius = query_radius
        self.weight = weight
        self.reactive_input_file = reactive_input_file

    @abstractmethod
    def calculate_bond_potential(
        self, atom1, atom2, atoms_dist: float, ideal_distance: float, weight: float, isolatedness: float, Cmax: int
    ) -> float:
        """Calculate the bond potential between two atoms.

        Args:
            atom1 (Atom): The first atom of the pair.
            atom2 (Atom): The second atom of the pair.
            atoms_dist (float): The distance between the two atoms.
            ideal_distance (float): The bonding equilibrium distance.

        Returns:
            float: The bond potential between the two atoms.
        """
        pass

    @abstractmethod
    def add_pair_pdb(self, pair):
        """Add the pair to the respective field in the PDB object.

        Args:
            pair (Pair): A selected pair from the selection sequence.
        """
        pass

    @abstractmethod
    def find_radicals(self, reactive_atoms_dict: dict[str, list[str]]):
        """Find radicals in the structure and adds them to the PDB object based on the implementation.

        Args:
            reactive_atoms_dict (dict[str, list[str]]): A dictionary with residue names as the keys and arrays of atom names as the values. These represent the reactive aspects of the structure.
        """
        pass

    def get_reactive_str_representation(
        self, reactive_input_file: str
    ) -> dict[str, list[str]]:
        """From the reactive atom/residue input file, extract the string representation.

        Args:
            reactive_input (str): The file containing the formatted atom names and residue names.

        Returns:
            dict[str, list[str]]: A dictionary where the keys are residue names and the values are an array of atom names corresponding to the residue name.
        """
        reactive_atom_dict = {}
        with open(reactive_input_file, "r") as f:
            for line in f:
                line = line.strip()
                columns = line.split()
                residue_name = columns[1]
                atom_name = columns[0]
                # If the residue name is not in the dictionary, add it with an empty list as the value
                if residue_name not in reactive_atom_dict:
                    reactive_atom_dict[residue_name] = []
                # Append the atom name to the list of atom names for this residue
                reactive_atom_dict[residue_name].append(atom_name)
        return reactive_atom_dict

    # the bond potential function to "rank" potential pairs
    def bond_potential(
        self, atom_dist: float, dist_equilibrium: float, dist_variance: float, iso_weight: float, isolatedness: float, Cmax: int
    ) -> float:
        """Calculate the bond potential between two atoms.

        Args:
            atom_dist (float): The distance between the two atoms.
            dist_equilibrium (float): The bonding distance equilibrium.
            dist_variance (float): The acceptable variance allowed for the bonding to occur.

        Returns:
            float: The bond potential between the two atoms.
        """
        return ((np.exp(-((atom_dist - dist_equilibrium) ** 2) / (2 * dist_variance)) + (iso_weight*isolatedness)) / (1 + iso_weight*Cmax)
        )
        

    def predict_bonding(self):
        """The predicting sequence for identifying preferred bonding sites for cross-linking.

        This method performs the following steps:
        1. Gets the reactive atoms from the input file.
        2. Initializes the reactive atoms.
        3. Converts the reactive atoms to coordinates.
        4. Creates a KDTree from the reactive atom coordinates.
        5. Performs a radius query on the KDTree to find nearby atoms.
        6. Filters out the first atom and same chain atoms from the query results.
        7. Initializes the potential pairs of atoms for bonding.
        8. Stores the potential pairs in a max heap for quick access to the pair with the highest potential.
        9. Performs the bond selection sequence.
        10. Finds the radicals in the structure that weren't bonded in the selection sequence.
        """
        reactive_atoms_dict = self.get_reactive_str_representation(
            self.reactive_input_file
        )

        self.init_reactive_atoms(reactive_atoms_dict)

        reactive_bonding_coords = atoms_to_coords(self.reactive_atoms)

        tree = KDTree(reactive_bonding_coords, leaf_size=30)
        # create a matrix of spatial query distances
        indices, distance = tree.query_radius(
            reactive_bonding_coords, r=self.query_radius, return_distance=True
        )
        # filter out the first atom and same chain atoms from the matrices
        indices_filtered, distances_filtered = remove_1st_and_chain(
            indices, distance, self.reactive_atoms
        )
        # initialize the pairs calculating their bond potential
        self.potential_pairs = self.initialize_potential_pairs(
            indices_filtered, distances_filtered
        )
        # store the potential pairs in a max heap for quick access to highest potential
        self.probability_heap = self.init_prob_heap()
        # the selection sequence for bond pairs
        self.bond_selection_loop()
        # find the radicals of the structure that weren't bonded in the selection sequence
        self.find_radicals(reactive_atoms_dict)

    # from the reactive atom and residue names return the atom objects from the structure
    def init_reactive_atoms(self, reactive_atoms_dict: dict[str, list[str]]):
        """Adds the reactive atoms to the self.reactive_atoms list.

        Args:
            reactive_atoms_dict (Dict[str, List[str]]): A dictionary with residue name as the key and an array of atom names as the values. These represent the reactive aspects of the structure.
        """
        for atom in self.atoms:
            if atom.name in reactive_atoms_dict.get(atom.res_name[:3], []):
                self.reactive_atoms.append(atom)

    def recal_probability_map(self, recalc_pairs: list):
        """Recalculates and sets the probability of the pairs in the recalculate pair list.

        Args:
            recalc_pairs (list[Pair]): The list of potential pairs that need to be recalculated based on the previous selected pair.
        """
        for pair in recalc_pairs:
            new_probability = self.calculate_bond_potential(
                pair.atom1, pair.atom2, pair.distance
            )
            pair.set_probability(new_probability)

    def select_highest_probability_pair(self):
        """Selects the root node of the heap and removes it from the potential pair set.

        This method also sets the properties of the pair to denote that they are bonded.

        Returns:
            Pair: The root node of the heap.
        """
        root_pair = heapq.heappop(self.probability_heap)
        self.potential_pairs.remove(root_pair)
        root_pair.bond_pair()
        return root_pair

    def init_prob_heap(self):
        """Initializes a heap of bonded pairs based on their bond potential.

        Returns:
            list[Pair]: A heap with the potential pairs.
        """
        probability_heap = []
        for pair in self.potential_pairs:
            heapq.heappush(probability_heap, pair)
        return probability_heap

    def initialize_potential_pairs(
        self, indices_filtered: np.ndarray, distances_filtered: np.ndarray
    ):
        """Initializes a list of potential pairs from the filtered indices and distances.

        Args:
            indices_filtered (np.ndarray): A numpy array of atom indices relating to the reactive_atoms instance variable list.
            distances_filtered (np.ndarray): A numpy array where each element is the distance between the atom at the corresponding index in the indices_filtered array and the atom at the index of that element.

        Returns:
            list[Pair]: A list of potential pairs.
        """
        potential_pairs = set()
        # iterate through the bond matrix and initialize the potential pairs
        for query_atom_idx, inner_array in enumerate(indices_filtered):
            current_atom = self.reactive_atoms[query_atom_idx]
            for idx_in_list, cur_nn_idx in enumerate(inner_array):

                cur_nn_atom = self.reactive_atoms[cur_nn_idx]

                pair_distance = distances_filtered[query_atom_idx][idx_in_list]

                probability_of_pair = self.calculate_bond_potential(
                    current_atom, cur_nn_atom, pair_distance
                )
                if probability_of_pair > 0:
                    bonded_pair = Pair(current_atom, cur_nn_atom, pair_distance)
                    bonded_pair.set_probability(probability_of_pair)
                    potential_pairs.add(bonded_pair)
        return list(potential_pairs)

    def bond_selection_loop(self):
        """Iteratively selects the potential pairs based on the highest bond potential and adds them to the respective field.

        This method performs the following steps:
        1. While the highest bond potential is greater than 0.001, continue the loop.
        2. Select the pair with the highest bond potential and add it to the PDB object.
        3. Find the pairs that need to be recalculated based on the selected pair.
        4. Recalculate the bond potential for the pairs that need to be recalculated.
        5. If there are pairs that were recalculated, reinitialize the heap of potential pairs.
        """
        while self.probability_heap[0].probability > 0.001:
            # reinitialize everytime since the last check (len(pair_to_recalc)) won't work correctly otherwise
            pairs_to_recalculate = []
            selected_pair = self.select_highest_probability_pair()
            self.add_pair_pdb(selected_pair)

            # finds the pairs that need to be recalculated
            pairs_to_recalculate = selected_pair.get_chain_branching_pairs(
                self.potential_pairs
            )
            # recalculates the probability of the pairs in the recalculate pair list by mutating the pair objects in the potential_pair_list
            self.recal_probability_map(pairs_to_recalculate)
            # reheapify the probability_heap with the updated probabilities only if you need, to optimize compute
            if len(pairs_to_recalculate) > 0:
                self.probability_heap = self.init_prob_heap()
