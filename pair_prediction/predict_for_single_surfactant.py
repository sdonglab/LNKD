import argparse
import os
from pair_prediction.tools.pdb import PDB
from pair_prediction.predict_polymerization import PredictBondsCore
from pair_prediction.predict_cycloaddition import PredictBondsSur


def controller():

    file_path_input, core_reactive_input, sur_reactive_input, output_directory = (
        handle_input()
    )
    pdb = PDB(file_path_input)

    predict_core_bonds = PredictBondsCore(pdb, core_reactive_input)
    predict_sur_bonds = PredictBondsSur(pdb, sur_reactive_input)
    predict_core_bonds.predict_bonding()
    predict_sur_bonds.predict_bonding()

    handle_output(file_path_input, output_directory, pdb)

    print("Finished!")

    return 1


# creates the parser object and handles the input from the user (via command line inputs)
def handle_input():
    parser = argparse.ArgumentParser(
        description="Prediction of cross linking bonds in micelle structure."
    )
    parser.add_argument("pdb_file_path", help="The path to the pdb input file")
    parser.add_argument(
        "core_reactive_input", help="The path to the core reactive input file"
    )
    parser.add_argument(
        "sur_reactive_input", help="The path to the surface reactive input file"
    )
    parser.add_argument(
        "-o", "--output_directory", help="The path to the output directory"
    )

    args = parser.parse_args()

    file_path_input = args.pdb_file_path
    core_reactive_input = args.core_reactive_input
    sur_reactive_input = args.sur_reactive_input

    if args.output_directory:
        output_dir = args.output_directory
    else:
        output_dir = os.path.dirname(file_path_input)

    return file_path_input, core_reactive_input, sur_reactive_input, output_dir


def handle_output(file_path_input, output_directory, pdb):

    base_path = os.path.splitext(os.path.basename(file_path_input))[0]

    pair_output_filename = base_path + "_pair_output.txt"
    sur_pair_output_filename = base_path + "_sur_pair_output.txt"
    core_pair_output_filename = base_path + "_core_pair_output.txt"
    radical_output_filename = base_path + "_radical_output.txt"

    pair_output_file = os.path.join(output_directory, pair_output_filename)
    sur_pair_output_file = os.path.join(output_directory, sur_pair_output_filename)
    core_pair_output_file = os.path.join(output_directory, core_pair_output_filename)
    radical_output_file = os.path.join(output_directory, radical_output_filename)

    pdb.write_pairs(pair_output_file, pdb.bonded_pairs, output_pymol_pairs=True)
    pdb.write_pairs(
        sur_pair_output_file, pdb.bonded_pairs_surface, output_pymol_pairs=False
    )
    pdb.write_pairs(
        core_pair_output_file, pdb.bonded_pairs_core, output_pymol_pairs=False
    )
    pdb.write_radicals(radical_output_file)


if __name__ == "__main__":
    controller()

