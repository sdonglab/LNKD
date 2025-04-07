# Preparing for Constrained MD Simulations

These directions indicate how to constrain the distances between the pairs predicted by LNKD in prepartion for topology updating.
NOTE: We perform our MD simulations with GROMACS 2022.5 compiled with PLUMED 2.9.0

## Preparing PLUMED File
The constrained_distances.py script takes the following arguments:
1. PDB file of your pre-crosslinked structure
2. Residue name of your template molecule (ex. TMP)
3. Text file with a list of atom pairs to constrain (output from LNKD)
4. Output file name (default = template.dat)
