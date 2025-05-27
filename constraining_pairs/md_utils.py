### Author: Bradford Derby ###

import os


def write_start_end_idx(output_file, ranges):

    for i, (start, end) in enumerate(ranges, 1):
        if i == 1:
            output_file.write(f"TEMPL: COM ATOMS={start}-{end}\n")
        else:
            output_file.write(f"DVB{i-1:02}: COM ATOMS={start}-{end}\n")


def write_distance_of_ranges(output_file, ranges):
    for i in range(1, len(ranges) + 1):
        if i == 1:
            output_file.write(f"dist{i:02}: DISTANCE ATOMS=centerTail,TEMPL\n")
        else:
            output_file.write(f"dist{i:02}: DISTANCE ATOMS=centerTail,DVB{i-1:02}\n")


def write_upper_wall(output_file, arg, at, label):
    output_file.write(
        f"UPPER_WALLS ARG={arg} AT={at} KAPPA=500000 EXP=2 OFFSET=0 LABEL={label}\n"
    )


def write_lower_wall(output_file, arg, at, label):
    output_file.write(
        f"LOWER_WALLS ARG={arg} AT={at} KAPPA=500000 EXP=2 OFFSET=0 LABEL={label}\n"
    )


def write_header(output_file, pdb_file, end_micelle_idx):
    output_file.write(
        f"""MOLINFO STRUCTURE={os.path.basename(pdb_file)}

WHOLEMOLECULES ...
        ENTITY0=1-{end_micelle_idx} # Micelle
        ADDREFERENCE 
...\n"""
    )


def write_matheval(output_file, leading_dec, arg, var, func):
    output_file.write(
        f"{leading_dec}: MATHEVAL ARG={arg} VAR={var} FUNC={func} PERIODIC=NO\n"
    )
