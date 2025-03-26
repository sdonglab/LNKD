def are_atoms_bonded(atom1, atom2):
    """Checks if either of the given atoms are bonded.

    Args:
        atom1 (Atom): The first atom to be checked.
        atom2 (Atom): The second atom to be checked.

    Returns:
        bool: True if either atom is bonded, False otherwise.
    """
    return atom1.is_bonded_external or atom2.is_bonded_external


def are_chains_bonded(chain1, chain2):
    """Checks if the first chain is in the list of chains bonded to the second chain.

    Args:
        chain1 (Chain): The chain to be checked.
        chain2 (Chain): The chain whose bonded chains will be checked.

    Returns:
        bool: True if the first chain is in the list of chains bonded to the second chain, False otherwise.
    """
    return chain1 in chain2.bonded_chains


def are_chains_same_type(chain1, chain2):
    """Checks if the two chains are of the same type.

    Args:
        chain1 (Chain): The first chain to be checked.
        chain2 (Chain): The second chain to be checked.

    Returns:
        bool: True if the chains are of the same type, False otherwise.
    """
    return chain1.chain_type == chain2.chain_type


def are_atoms_same_chain(atom1, atom2):
    """Checks if the two atoms are part of the same chain.

    Args:
        atom1 (Atom): The first atom to be checked.
        atom2 (Atom): The second atom to be checked.

    Returns:
        bool: True if the atoms are part of the same chain, False otherwise.
    """
    return atom1.chain == atom2.chain


def is_valid_surface_pair(atom1, atom2):
    """Checks if the given atoms form a valid surface pair.

    Args:
        atom1 (Atom): The first atom to be validated.
        atom2 (Atom): The second atom to be validated.

    Returns:
        bool: False if the atoms are already bonded or are in the same chain. True otherwise.
    """
    return (not are_atoms_bonded(atom1, atom2)) and (
        not are_chains_same_type(atom1.chain, atom2.chain)
    )

def is_valid_LN2_surface_pair(atom1, atom2):
    """Similar to is_valid_surface_pair but also checks to see if the chains have already formed a bond to avoid 2 alkynes on 1 surfactant bound to 2 azides on the same linker
    
    Args:
        atom1 (Atom): The first atom to be validated.
        atom2 (Atom): The second atom to be validated.

    Returns:
        bool: False if the atoms are already bonded or are in the same chain or if the chains are already bonded. True otherwise.
    """
    return (not are_atoms_bonded(atom1, atom2)) and (
        not are_chains_same_type(atom1.chain, atom2.chain)) and (
        not are_chains_bonded(atom1.chain, atom2.chain)
     )

def is_valid_core_pair(atom1, atom2):
    """Checks if the given atoms form a valid core pair.

    Args:
        atom1 (Atom): The first atom to be validated.
        atom2 (Atom): The second atom to be validated.

    Returns:
        bool: False if the atoms are bonded or if the chains that the atoms are a part of have already bonded. True otherwise.
    """
    return (not are_atoms_bonded(atom1, atom2)) and (
        not are_chains_bonded(atom1.chain, atom2.chain)
    )
