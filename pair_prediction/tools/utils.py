### Author: Bradford Derby ###

from tools.constraint_validation import are_atoms_same_chain
import numpy as np


def atoms_to_coords(atoms):
    """Converts a list of Atom objects to a numpy array of their coordinates.

    Args:
        atoms (List[Atom]): A list of Atom objects.

    Returns:
        np.ndarray: A numpy array where each row represents the coordinates of an atom.
    """
    coords = np.empty((len(atoms), 3))
    for i, atom in enumerate(atoms):
        coords[i] = [atom.x, atom.y, atom.z]
    return coords


def remove_1st_and_chain(
    nn_indices: np.ndarray, distances_matrix: np.ndarray, reactive_atoms
) -> tuple[np.ndarray, np.ndarray]:
    """Filters out the first column of the indices and distances, which is the atom itself, and those atoms that are in the same chain.

    Args:
        nn_indices (np.ndarray): An array of nearest neighbor indices.
        distances_matrix (np.ndarray): A matrix of distances between atoms.
        reactive_atoms (List[Atom]): A list of reactive atoms.

    Returns:
        Tuple[np.ndarray, np.ndarray]: A tuple containing the filtered indices and distances.
    """
    indices_filtered = []
    distances_filtered = []

    # Filter out the atoms that are in the same chain
    for query_atom_idx, inner_array in enumerate(nn_indices):
        current_atom = reactive_atoms[query_atom_idx]
        inner_indices = []
        inner_distances = []

        for idx_in_list, cur_nn_idx in enumerate(inner_array):
            cur_nn_atom = reactive_atoms[cur_nn_idx]

            if not are_atoms_same_chain(current_atom, cur_nn_atom):
                inner_indices.append(cur_nn_idx)
                inner_distances.append(distances_matrix[query_atom_idx][idx_in_list])

        indices_filtered.append(inner_indices)
        distances_filtered.append(inner_distances)

    return indices_filtered, distances_filtered


def calc_distance(atom1, atom2) -> float:
    """Calculates the Euclidean distance between two Atom objects.

    Args:
        atom1 (Atom): The first Atom object.
        atom2 (Atom): The second Atom object.

    Returns:
        float: The Euclidean distance between the two Atom objects.
    """
    distance = np.sqrt(
        (atom1.x - atom2.x) ** 2 + (atom1.y - atom2.y) ** 2 + (atom1.z - atom2.z) ** 2
    )
    return distance
