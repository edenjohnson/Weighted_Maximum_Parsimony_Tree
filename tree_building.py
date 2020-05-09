########################################################################################################################
# tree_building.py
# Author: Eden Johnson
# Last Modified: May 9, 2020
# Program: tree_building.py opens up the aligned fasta sequences obtained from MSA.py and constructs a weighted maximum
#          parsimony tree on the results. File outputs an annotated tree that serves as input for closest_OTU.py to
#          return the closest OTU(s) for each OTU.
########################################################################################################################

# Import modules
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import *
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Phylo import PhyloXML


print("\n\nMSA completed. Computing weighted maximum parsimony tree now....\n")

# relabeling the species to their common names for input into the tree
with open("BRCA2_family_aligned.fasta", "rU") as input_handle, open("BRCA2_family_fixed.fasta", "w") as output_handle:
    sequences = SeqIO.parse(input_handle, "fasta")

    fixed_sequences = []

    for line in sequences:
        if "NM_000059.4" == line.id:
            record = SeqRecord(line.seq, "Human")
            fixed_sequences.append(record)
        elif "NM_001006653.4" == line.id:
            record = SeqRecord(line.seq, "Dog")
            fixed_sequences.append(record)
        elif "NM_001009858.1" == line.id:
            record = SeqRecord(line.seq, "Cat")
            fixed_sequences.append(record)
        elif "NM_001081001.2" == line.id:
            record = SeqRecord(line.seq, "Mouse")
            fixed_sequences.append(record)
        elif "NM_031542.2" == line.id:
            record = SeqRecord(line.seq, "Rat")
            fixed_sequences.append(record)
        elif "NM_204276.2" == line.id:
            record = SeqRecord(line.seq, "Chicken")
            fixed_sequences.append(record)

    SeqIO.write(fixed_sequences, output_handle, "fasta")

input_handle.close()
output_handle.close()


# convert the clustalW format to phylip for the program
from Bio import AlignIO
AlignIO.convert("BRCA2_family_fixed.fasta", "fasta", "BRCA2_family.phy", "phylip")

# Read the sequences and align
aln = AlignIO.read('BRCA2_family.phy', 'phylip')

# create a starting tree with NJ
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)
constructor = DistanceTreeConstructor(calculator, 'nj')
starting_tree = constructor.build_tree(aln)


# A substitution cost matrix, used from in-lecture excise (penalty of 2 for transversion and gap, penalty of 1 for
# # transition)
cost_matrix = [[0],
            [2,0],
            [1,2,0],
            [2,1,2,0],
            [2,2,2,2,0]]

# weighted cost matrix corresponds to
weight = DistanceMatrix(names=['A', 'C', 'G', 'T','-'], matrix=cost_matrix)

# ParsimonyScorer will use the Sankoff Algorithm when provided with a matrix argument
scorer = ParsimonyScorer(matrix=weight)

searcher = NNITreeSearcher(scorer)

# pars_constructor will create the information necessary to build the tree
pars_constructor = ParsimonyTreeConstructor(searcher, starting_tree)
pars_tree = pars_constructor.build_tree(aln)
pars_tree.rooted = False
print(pars_tree)


# Promote the basic tree to PhyloXML
pars_phy = pars_tree.as_phyloxml()
# Save the annotated phyloXML file
Phylo.write(pars_phy, 'parsimony_tree.xml', 'phyloxml')

# Draw the phylogenetic tree
Phylo.draw(pars_tree, branch_labels=lambda c: "%.3f" % (c.branch_length))

# Print the phylogenetic tree in the terminal
print('\nPhylogenetic Tree\n===================')
Phylo.draw_ascii(pars_tree)
