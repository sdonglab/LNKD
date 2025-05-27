### Author: Bradford Derby ###

class Atom:
    """Represents an atom in a PDB file or structure

    Attributes:
        recordName (str): The record name of the atom, defaults to "ATOM".
        serial (int): The serial number of the atom.
        name (str): The name of the atom.
        res_name (str): The name of the residue the atom belongs to.
        chain (Chain): The chain object the atom belongs to.
        res_seq (int): The sequence number of the residue the atom belongs to.
        i_code (str): The insertion code of the atom.
        x (float): The x-coordinate of the atom.
        y (float): The y-coordinate of the atom.
        z (float): The z-coordinate of the atom.
        occupancy (float): The occupancy of the atom.
        tempFactor (float): The temperature factor of the atom.
        element (str): The element symbol of the atom.
        charge (float): The charge of the atom.
        is_bonded_external (bool): A boolean value indicating whether the atom is externally bonded. Defaults to False.
        is_paired (bool): A boolean value indicating whether the atom is part of a radical pair. Defaults to False.
    """

    def __init__(
        self,
        serial,
        name,
        res_name,
        chain,
        res_seq,
        i_code,
        x,
        y,
        z,
        occupancy,
        tempFactor,
        element,
        charge,
        is_bonded_external=False,
        is_paired=False,
    ):
        self.recordName = "ATOM"
        self.serial = serial
        self.name = name
        self.res_name = res_name
        self.chain = chain  # chain object
        self.res_seq = res_seq
        self.i_code = i_code
        self.x = x
        self.y = y
        self.z = z
        self.occupancy = occupancy
        self.tempFactor = tempFactor
        self.element = element
        self.charge = charge
        self.is_bonded_external = is_bonded_external  # boolean value for bonding logic
        self.is_paired = is_paired  # boolean value for radical pair logic
        if self.occupancy:
            self.occupancy = float(self.occupancy)
        if self.tempFactor:
            self.tempFactor = float(self.tempFactor)

    def __eq__(self, other):
        return self.serial == other.serial
