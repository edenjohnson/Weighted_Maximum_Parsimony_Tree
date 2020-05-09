########################################################################################################################
# wmpt_script.py
# Author: Eden Johnson
# Last Modified: May 9, 2020
# Program: Script imports the os module to call 5 other .py files to construct a weighted maximum parsimony tree for the
#          BRCA2 gene in Homo sapiens along with its orthologs. Output expected are updates showing script process,
#          along with the tree output. Additional output of returning each OTU's closest OTU(s) is provided. Note that
#          this script will take ~ 1 hour to execute since it will be performing a BLAST search over the internet and
#          performs a multiple sequence alignment locally.
# Dependencies: User MUST pip install biopython and pip install matplotlib and have the clustal_omega binary file in the
#               same directory as this script.
########################################################################################################################

# import the OS module
import os

# call the following files to the command line via the OS module
os.system("python3 BLAST.py")
os.system("python3 parse_sequences.py")
os.system("python3 MSA.py")
os.system("python3 tree_building.py")
os.system("python3 closest_OTU.py")
