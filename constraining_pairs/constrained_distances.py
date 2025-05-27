### Author: Bradford Derby ###

import argparse
from md_utils import (
    write_distance_of_ranges,
    write_header,
    write_matheval,
    write_start_end_idx,
    write_upper_wall,
)


def controller():
    parser = argparse.ArgumentParser(
        prog="constrained_distances.py",
        description="Generates the eq file for the equilibrium protocol",
    )
    parser.add_argument("pdb_file", help="The pdb file to be parsed")
    parser.add_argument(
        "template", help="The template to search for - 1AR for instance"
    )
    parser.add_argument(
        "post_poly_bonds",
        help="The post polymerization bonds from the prediction sequence",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        help="The output file to write to... default is template.dat",
        default="template.dat",
    )
    args = parser.parse_args()
    generate_eq_add_file(
        args.pdb_file, args.template, args.post_poly_bonds, args.output_file
    )


def get_pair_data(post_poly_bonds):
    with open(post_poly_bonds, "r") as f:
        pairs_data = []
        lines = f.readlines()
        for i in range(0, len(lines), 2):
            line1 = lines[i].strip()
            line2 = lines[i + 1].strip()
            pair_data = process_lines_pair(line1, line2)
            pairs_data.append(pair_data)

    return pairs_data

def process_lines_pair(line1, line2):

    col_line1 = line1.split(" ")  # C3 N1
    col_line2 = line2.split(" ")  # N3 C2
    pair_data = col_line1 + col_line2

    return pair_data


def write_position_section(f, pairs_data):
    for i, pair_data in enumerate(pairs_data, 1):
        f.write(f"p{i:02}_c2: POSITION ATOM={pair_data[3]}\n")
    f.write("\n")

    for i, pair_data in enumerate(pairs_data, 1):
        f.write(f"p{i:02}_c3: POSITION ATOM={pair_data[0]}\n")
    f.write("\n")

    for i, pair_data in enumerate(pairs_data, 1):
        f.write(f"p{i:02}_n3: POSITION ATOM={pair_data[2]}\n")
    f.write("\n")

    for i, pair_data in enumerate(pairs_data, 1):
        f.write(f"p{i:02}_n1: POSITION ATOM={pair_data[1]}\n")
    f.write("\n")

def get_indices(pdb_file, template):
    prev_res_seq = ""
    start_idx = 0
    end_idx = 0
    dvb_templ_idx_ranges = []
    c6_indices = []
    nitrogen_indices = []
    with open(pdb_file, "r") as file:
        cur_serial_idx = 0
        for line in file:
            if line.startswith("ATOM"):
                cur_serial_idx += 1
                columns = line.split()
                res_name = columns[3][:3]
                res_seq = columns[4]
                atom_name = columns[2]
                if res_name in ["SUR"]:
                    if atom_name == "C6":
                        c6_indices.append(str(cur_serial_idx))
                    elif atom_name == "N":
                        nitrogen_indices.append(cur_serial_idx)

                elif res_name in ["DVO", "DVM", "DVP", template]:
                    if prev_res_seq != res_seq:
                        if prev_res_seq != "":
                            dvb_templ_idx_ranges.append((start_idx, end_idx))
                        prev_res_seq = res_seq
                        start_idx = cur_serial_idx
                    end_idx = cur_serial_idx

    if prev_res_seq != "":
        dvb_templ_idx_ranges.append((start_idx, end_idx))

    return c6_indices, nitrogen_indices, dvb_templ_idx_ranges

def generate_eq_add_file(pdb_file, template, post_poly_bonds, output_file):
    pairs_data = get_pair_data(post_poly_bonds)
    c6_indices, nitrogen_indices, dvb_templ_idx_ranges = get_indices(pdb_file, template)

    with open(output_file, "w") as f:

        write_header(f, pdb_file, dvb_templ_idx_ranges[-1][1])

        f.write("\n")

        write_position_section(f, pairs_data)

        for i, pair_data in enumerate(pairs_data, 1):
            f.write(f"p{i:02}_com: COM ATOMS={pair_data[3]},{pair_data[0]}\n")

        f.write("\n")

        for i, pair_data in enumerate(pairs_data, 1):
            f.write(
                f"p{i:02}_dist: DISTANCE ATOMS=p{i:02}_com,{int(pair_data[2]) - 1}\n"
            )

        f.write("\n")
        # Previous part -----------------------------------------------------------------------------
        f.write(f"centerTail: COM ATOMS={','.join(c6_indices)}\n")

        f.write("\n")

        write_start_end_idx(f, dvb_templ_idx_ranges)

        f.write(
            f"ALLDVB: COM ATOMS={dvb_templ_idx_ranges[1][0]}-{dvb_templ_idx_ranges[-1][1]}\n"
        )

        f.write("\n")

        write_distance_of_ranges(f, dvb_templ_idx_ranges)

        f.write("\n")

        for i, n_index in enumerate(nitrogen_indices, 1):
            f.write(f"distN{i:02}: DISTANCE ATOMS=ALLDVB,{n_index}\n")

        f.write("\n")
        # UPPER_WALLS ARG=dist01 AT=2.20 KAPPA=500000 EXP=2 OFFSET=0 LABEL=dist01-uwall
        for i in range(1, len(dvb_templ_idx_ranges) + 1):
            write_upper_wall(f, f"dist{i:02}", "2.20", f"dist{i:02}-uwall")

        f.write("\n")

        # UPPER_WALLS ARG=distN01 AT=2.50 KAPPA=500000 EXP=2 OFFSET=0 LABEL=distN01-uwall
        for i in range(1, len(nitrogen_indices) + 1):
            write_upper_wall(f, f"distN{i:02}", "2.50", f"distN{i:02}-uwall")

        f.write("\n")

        # UPPER_WALLS ARG=p01_dist AT=XXXX KAPPA=500000 EXP=2 OFFSET=0 LABEL=p01_dist-uwall
        for i in range(1, len(pairs_data) + 1):
            write_upper_wall(f, f"p{i:02}_dist", "XXXX", f"p{i:02}_dist-uwall")

        f.write("\n")

        f.write("ENDPLUMED\n")


if __name__ == "__main__":
    controller()
