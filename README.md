# LNKD
A resource efficient algorithm for predicting polymer crosslinking by Linking Nodes in KD-trees (LNKD)

## Associated Manuscript
The following paper describes the methodology and benchmarking in detail. Please cite this paper should you use or reference this work:

Stevens, E.; Curtolo, F.; Derby, B.; Wallace, D.; Zhao, Y.; Dong, S. S. Modeling Molecularly Imprinted Nanoparticles with LNKD: A Resource Efficient Algorithm for Polymer Crosslinking. _ChemRxiv_ June 13, 2025. https://doi.org/10.26434/chemrxiv-2025-sfrfh.

## Motivation
Crosslinked polymer materials have extensive applications in industrial, pharmaceutical, energy, and health related fields. Computational efforts to model, characterize, and design these materials are limited in part due to challenges in solving their matrix-like and probabilistic structures experimentally. Despite work to predict polymer crosslinking, it remains a computationally expensive and underexplored task. We propose LNKD, a resource efficient algorithm for predicting pairs of reactive atoms in pre-crosslinked 3D structures to form crosslinked polymers. 

## Algorithm Overview
LNKD uses a dictionary of reactive atoms to parse a PDB file and generate a KD-tree with each node representing the coordinate of an atom able to participate in a crosslink. The KD-tree is then spatially queried to find all other reactive atoms within a query radius, resulting in a bonding matrix of all possible crosslinking pairs. A bond potential, or probability of forming a crosslink, is calculated for each pair (Equation 1 in associated manuscript) and stored in a max heap. The root node of the max heap is selected as a crosslinking pair and removed from the heap. Potential pairs whose probabilities of forming a bond have changed due to the selected pair have their bond potentials recalculated and the max heap is reordered. This bond selection loop is performed until the root node is sufficiently close to null. The output is a list of pairs of atoms predicted to form bonds in the original structure. In the case of radical polymerization, any remaining radicals are also output so the user can later add hydrogen atoms to correct the charges.

## Use Cases
LNKD can be modified by users to be broadly applied to crosslink prediction in any polymer network. We have used the method to simultaneously predict both of the following types of crosslinks in structures of pre-crosslinked micelles to model molecularly imprinted nanoparticles (MINPs) from their monomers:
1. Azide-alkyne cycloaddition to form triazoles on the surface of the micelles
2. radical polymerization of alkenes to form carbon-carbon bonds in the core of the micelles

## Using LNKD

### Requirements

Code runs on Python3 and uses the modules listed in the requirements.txt file - these can be installed using the following command:
```text
$ pip install -r requirements.txt
```

### Input
1. Path to PDB structure of pre-crosslinked system
2. Path to text file containing core reactive atoms (polymerization)
3. Query radius for finding nearby reactive atoms in the core (float, angstrom)
4. Weight for the degree of isolation when predicting core pairs (float)
5. Path to text file containing surface reactive atoms (cycloaddition)
6. Query radius for finding nearby reactive atoms in the surface (float, angstrom)
4. Weight for the degree of isolation when predicting surface pairs (float)
7. Path to output directory (optional, default = current directory)

The format of the text file of reactive atoms should be as follows:
```text
ATOMi ResidueI
ATOMj ResidueJ
ATOMk ResidueK
```
Where the first column has the atom name (stored in columns 13-16 of your PDB file) and the second column is the name of the residue type the atom belongs to (stored in columns 18-20 of your PDB file).

NOTE: Our surface pair prediction is set up for azide alkyne cycloaddition. The list of surface reactive atoms includes the terminal carbon of each of the three alkynes on a surfactant (SUR) and the terminal nitrogen of the two azides on our diazide linker (LN2). The code automatically uses those atoms to find the other carbon and nitrogen that would form a bond during the cycloaddition. This ensures the algorithm doesn't form bonds to two different azides from the same alkyne for example. This can be modified to work for any crosslink that includes the formation of more than one bond. Crosslinks involving the formation of only one bond, such as in radical polymerization, should use our core bond prediction scripts. In the provided example, core reactive atoms include carbons of all alkenes able to be polymerized, including methacrylate groups on SUR and alkenes on ortho, meta, or para substituted divinyl benzene (DVO, DVM, DVP).

### Running from Command Line

Below is the command used to generate the output files found in the 'example' directory. Edit paths, query radii, and weights, for your system.
```text
python3 ./pair_prediction/predict_for_single_surfactant.py ./example/input/1-PSA_MINP_build1_round2_pair_opt.pdb ./example/input/core_reactive.txt -core_QR 10 -core_W -0.05 ./example/input/sur_reactive.txt -sur_QR 15 -sur_W -0.05 -o ./example/output/
```

### Output
LNKD outputs four files. They are automatically named starting with "name_of_input_PDB_file_" and ending with the following four distinctions:
1. core_pair_output.txt
2. sur_pair_output.txt
3. radical_output.txt
4. pair_output.txt

The first two files contain the lists of predicted polymerization (core) and cycloaddition (surface) pairs respectively. The third includes a list of atoms from the core reactive atoms that became radicals due to their neighboring carbon being paired, but not being predicted to form a pair itself. The final output combines the list of surface and core pairs and can be used to visualize them in PyMOL or other molecular graphics software.

---
© 2025 Northeastern University. All rights reserved. Any commercial use of this work without explicit written permission from the copyright holder is strictly prohibited. 
