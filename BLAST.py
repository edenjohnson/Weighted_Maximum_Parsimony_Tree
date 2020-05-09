########################################################################################################################
# BLAST.py
# Author: Eden Johnson
# Last Modified: May 9, 2020
# Program: BLAST.py performs a BLAST search for the BRCA2 gene in Homo sapiens via its accession number. An entrez query
#          is used to filter out unwanted predicted structures/models/unknown/experimental genes accession numbers
#          from the resulting subject ids generated from the search. Results are written to a file to be opened in
#          parse_sequences.py
########################################################################################################################

# Import modules
from Bio.Blast import NCBIWWW

# entrez query to filter out unwanted predicted structures/models/unknown/experimental genes
entrez_q = "(all [filter] NOT(environmental samples[organism] OR metagenomes[organism] OR Homo sapiens[Organism] OR " +\
           "Model[Text Word] OR Predicted[Text Word] OR unknown[Text Word] OR Experimental[Text Word]))"

print("Starting the query....\n")

# perform a BLAST query for the BRCA2 gene in Humans
result = NCBIWWW.qblast(program="blastn",database="refseq_rna",sequence="NM_000059.4",entrez_query=entrez_q,expect=0.01)

# Search for homologs of a protein sequence using BLAST.
with open("my_blast.xml", "w") as out_handle:
    out_handle.write(result.read())

print("\nWriting to the file....\n")

# close the result handle from the query search to only read from it once. We can now read from it from the xml file
result.close()
