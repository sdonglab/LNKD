# LNKD
A resource efficient algorithm for predicting polymer crosslinking

## Associated Manuscript
Provide citation here:
Kindly cite the above paper should you use or reference this work.

## Motivation
Crosslinked polymer materials have extensive applications in industrial, pharmaceutical, energy, and health related fields. Computational efforts to model, characterize, and design these materials are limited in part due to challenges in solving their matrix-like and probabilistic structures experimentally. Despite work to predict polymer crosslinking, it remains a computationally expensive and underexplored task. We propose LNKD, a resource efficient algorithm for predicting pairs of reactive atoms in pre-crosslinked 3D structures to form crosslinked polymers. 

## Algorithm Overview
LNKD uses a dictionary of reactive atoms to parse a PDB file and generate a KD-tree with each node representing the coordinate of an atom able to participate in a crosslink. The KD-tree is then spatially queried to find all other reactive atoms within a query radius, resulting in a bonding matrix of all possible crosslinking pairs. A bond potential, or probability of forming a crosslink, is calculated for each pair as described in the associated paper and stored in a max heap. The root node of the max heap is selected as a crosslinking pair and removed from the heap. Potential pairs whose probabilities of forming a bond have changed due to the selected pair have their bond potentials recalculated and the max heap is reordered. This bond selection loop is performed until the root node is sufficiently close to null. The output is a list of pairs of atoms predicted to form bonds in the original structure. In the case of radical polymerization, any remaining radicals are also output so the user can later add hydrogen atoms to correct the charges.

## Use Cases
LNKD can modified by users to be broadly applied to crosslink prediction in any polymer network. We have used the method to simultaneously predict both of the following types of crosslinks in our structures:
1. Azide-alkyne cycloaddition to form triazoles
2. radical polymerization of alkenes to form carbon-carbon bonds

insert notes about what parts of our algorithm are specific to these cases/what they would have to change

## Using LNKD

### Input
1. Path to PDB structure of pre-crosslinked system
2. Path to text file containing core reactive atoms (polymerization)
3. Path to text file containing surface reactive atoms (cycloaddition)
4. Path to output directory (optional, default = current directory)

see if I can edit scripts to take weight and query radius for both the surface and the core as inputs as well so that users don't have to edit scripts each time

The format of the text file of reactive atoms should be as follows:
```text
ATOM1 ResidueName1
ATOM2 ResidueName2
ATOM3 ResidueName3
```
Where the first column has the atom name (stored in columns 13-16 of your PDB file) and the second column is the name of the residue type the atom belongs to (stored in columns 18-20 of your PDB file).

