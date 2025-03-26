from pair_prediction.tools.atom import Atom
from pair_prediction.tools.chain import Chain


class PDB:
    """Represents the PDB structure from the given pdb file.

    Attributes:
        file (str): The path to the pdb file.
        atoms (List[Atom]): The list of Atom objects in the structure.
        chains (Dict[int, Chain]): A dictionary where the key is the res_seq number and the value is the Chain object. This represents all the molecules/chains in the structure.
        bonded_pairs (List[Pair]): Represents all the bonded pairs in the structure.
        bonded_pairs_core (List[Pair]): Represents all the core bonded pairs in the structure.
        bonded_pairs_surface (List[Pair]): Represents all the surface bonded pairs in the structure.
        radicals (List[Atom]): Represents the radical atoms in the structure after the bond selection sequence.
    """

    def __init__(self, file):
        self.file = file
        self.atoms = []
        self.chains = {}
        self.bonded_pairs = []
        self.bonded_pairs_core = []
        self.bonded_pairs_surface = []
        self.parse()
        self.radicals = []

    def parse(self):
        """Parses the PDB file and stores the atoms in the structure.

        This method reads the PDB file line by line. If a line starts with "ATOM", it assigns the atom data from that line to the structure.
        """
        with open(self.file, "r") as f:
            for line in f:
                if line.startswith("ATOM"):
                    self.assign_atom(line)

    def write_pymol_pairs(self, file: str, pairs):
        """Formats the pairs for viewing in pymol easily.

        Args:
            file (str): The path to the file where the formatted pairs will be written.
            pairs (List[Pair]): A list of pairs to be formatted.

        This method writes each pair in the format "1242+3245" to the specified file.
        """
        for pair in pairs:
            file.write(f"+{pair.atom1.serial}+{pair.atom2.serial}")

    def write_columned_pairs(self, file: str, pairs):
        """Formats the pair output in space separated columns Atom1 Atom2.

        Args:
            file (str): The path to the file where the formatted pairs will be written.
            pairs (List[Pair]): A list of pairs to be formatted.

        This method writes each pair in the format "1242 3245" to the specified file,
        """
        for pair in pairs:
            file.write(f"{pair.atom1.serial} {pair.atom2.serial}\n")

    def write_pairs(self, file: str, pairs, output_pymol_pairs: bool):
        """Writes the pairs to the specified file.

        Args:
            file (str): The path to the file where the pairs will be written.
            pairs (List[Pair]): A list of pairs to be written to the file.
            output_pymol_pairs (bool): A boolean value to determine if pymol pair format should be written at the top of the file.

        This method writes each pair to the specified file. If output_pymol_pairs is True, it also writes the pairs in pymol format at the top of the file.
        """
        with open(file, "w") as f:
            if output_pymol_pairs:
                self.write_pymol_pairs(f, pairs)
                f.write("\n")
            self.write_columned_pairs(f, pairs)

    def write_radicals(self, file: str):
        """Writes the radicals to the specified file.

        Args:
            file (str): The path to the file where the radicals will be written.

        This method writes each radical's serial number to the specified file, each on a new line.
        """
        with open(file, "w") as f:
            for radical in self.radicals:
                f.write(f"{radical.serial}\n")

    def add_bonded_pair(self, pair):
        """Adds a bonded pair to the list of bonded pairs.

        Args:
            pair (Pair): The pair to be added to the list of bonded pairs.
        """
        self.bonded_pairs.append(pair)

    def add_core_bonded_pair(self, pair):
        """Adds a core bonded pair to the list of core bonded pairs.

        Args:
            pair (Pair): The pair to be added to the list of core bonded pairs.
        """
        self.bonded_pairs_core.append(pair)

    def add_surface_bonded_pair(self, pair):
        """Adds a surface bonded pair to the list of surface bonded pairs.

        Args:
            pair (Pair): The pair to be added to the list of surface bonded pairs.
        """
        self.bonded_pairs_surface.append(pair)

    def add_radical(self, radical):
        """Adds a radical atom to the list of radicals.

        Args:
            radical (Atom): The radical atom to be added to the list of radicals.
        """
        self.radicals.append(radical)

    def get_atoms(self):
        """Gets the list of atoms in the PDB structure.

        Returns:
            List[Atom]: The list of atoms in the PDB structure.
        """
        return self.atoms

    def get_bonded_pairs(self):
        """Gets the list of bonded pairs in the PDB structure.

        Returns:
            List[Pair]: The list of bonded pairs in the PDB structure.
        """
        return self.bonded_pairs
    
    def get_bonded_pairs_core(self):
        """Gets the list of bonded  corepairs in the PDB structure.
        """

        return self.bonded_pairs_core
    
    def get_bonded_pairs_surf(self):
        """Gets the list of bonded  corepairs in the PDB structure.
        """

        return self.bonded_pairs_surface

    def assign_chain(self, line: str):
        """Initializes a Chain object and adds it to the dictionary of chains.

        This method parses a line from a PDB file to extract chain information, creates a Chain object
        with this information if it doesn't already exist, and adds the chain to the dictionary of chains.

        Args:
            line (str): A line from a PDB file containing chain information.

        Returns:
            Chain: The Chain object corresponding to the given line.
        """
        res_name = line[17:20].strip()
        res_seq = int(line[22:26].strip())
        # if the chain is already initialized just return it as we don't want to reinitialize it
        if res_seq in self.chains:
            return self.chains[res_seq]
        else:
            self.chains[res_seq] = Chain(res_seq, res_name)
            return self.chains[res_seq]

    def assign_atom(self, line: str):
        """Initializes an Atom object and adds it to the chain and the list of atoms.

        This method parses a line from a PDB file to extract atom information, creates an Atom object
        with this information, adds the atom to its chain, and appends the atom to the list of atoms.

        Args:
            line (str): A line from a PDB file containing atom information.

        """
        serial = int(line[6:11].strip())
        name = line[12:17].strip()
        res_name = line[17:21].strip()
        chain = self.assign_chain(line)
        res_seq = int(line[22:26].strip())
        i_code = line[26:27].strip()
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        occupancy = line[54:60].strip()
        tempFactor = line[60:66].strip()
        element = line[76:78].strip()
        charge = line[78:80].strip()
        atom = Atom(
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
        )
        chain.add_atom(atom)
        self.atoms.append(atom)
