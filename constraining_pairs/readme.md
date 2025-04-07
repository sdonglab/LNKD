# Preparing for Constrained MD Simulations

These directions indicate how to constrain the distances between the pairs predicted by LNKD in prepartion for topology updating.
NOTE: We perform our MD simulations with GROMACS 2022.5 compiled with PLUMED 2.9.0

## Preparing the PLUMED File
The constrained_distances.py script takes the following arguments:
1. PDB file of your pre-crosslinked structure
2. Residue name of your template molecule (ex. TMP)
3. Text file with a list of atom pairs to constrain (output from LNKD)
4. Output file name (default = template.dat)

For example, to generate a PLUMED file for the surface pair output in the example output directory, run the following command:
``` text
Python3 constraining_pairs/constrained_distances.py example/input/1-PSA_MINP_build1_round2_pair_opt.pdb PSA example/output/1-PSA_MINP_build1_round2_pair_opt_sur_pair_output.txt -o template.dat
```
