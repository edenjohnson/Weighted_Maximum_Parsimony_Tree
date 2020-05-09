########################################################################################################################
# MSA.py
# Author: Eden Johnson
# Last Modified: May 9, 2020
# Program: MSA.py opens up the fasta sequences obtained from parse_sequences.py and performs a MSA on the sequences of
#          the file. This program uses command line tools of Clustal Omega. Aligned sequences are output to the
#          BRCA2_family_aligned.fasta file to serve as input for the tree_building.py file.
# Dependencies: User MUST have the clustal_omega binary file in the same directory as the wmpt_script.py file.
########################################################################################################################

# Import Clustal Omega wrapper
from Bio.Align.Applications import ClustalOmegaCommandline
import os

# create variables with paths to the input and output files
in_file = "BRCA2_family.fasta"
out_file = "BRCA2_family_aligned.fasta"

# Get the command for Clustal Omega
clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)

# adjust the executable command to run with the binary file in the directory
command_list = str(clustalomega_cline).split()
for n, i in enumerate(command_list):
    if i == "clustalo":
        command_list[n] = "./clustal-omega-1.2.3-macosx"

# create a string of the command to input into the command line
command = " ".join(command_list)

# call the command
os.system(command)

