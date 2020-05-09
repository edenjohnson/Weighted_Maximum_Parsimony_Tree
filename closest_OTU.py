########################################################################################################################
# closest_OTU.py
# Author: Eden Johnson
# Last Modified: May 9, 2020
# Program: closest_OTU.py opens up the annotated tree file generated from tree_building.py and uses its information
#          to return the closest OTU(s) for each OTU.
########################################################################################################################

# import the required module
from Bio import Phylo

# open the tree file
pars_tree = Phylo.read('parsimony_tree.xml', 'phyloxml')


# for a given OTU, return its closest OTU
def closest_otu(otu):
    """ Function returns the closest OTU given an OTU. If OTU has no sister OTU, returns list of all closest OTUs. """
    non_terminals = pars_tree.get_nonterminals()
    leaf = pars_tree.get_terminals()
    min_leaves = len(leaf)
    smallest_clade = leaf
    for node in non_terminals:
        leaves = node.get_terminals()
        leaf_length = len(leaves)
        if leaf_length < min_leaves and otu in str(leaves):
            min_leaves = leaf_length
            smallest_clade = leaves
    leaf_names = []
    for leaf in smallest_clade:
        if leaf.name != otu:
            leaf_names.append(leaf.name)
    return leaf_names


# Results
print("Returning each OTU's closest OTU(s)")

print("\nDog's closest OTU(s):")
print(closest_otu("Dog"))

print("\nCat's closest OTU(s):")
print(closest_otu("Cat"))

print("\nHuman's closest OTU(s):")
print(closest_otu("Human"))

print("\nRat's closest OTU(s):")
print(closest_otu("Rat"))

print("\nMouse's closest OTU(s):")
print(closest_otu("Mouse"))

print("\nChicken's closest OTU(s):")
print(closest_otu("Chicken"))
