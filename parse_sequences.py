########################################################################################################################
# parse_sequence.py
# Author: Eden Johnson
# Last Modified: May 9, 2020
# Program: parse_sequences.py opens up the results generated from BLAST.py and performs an Entrez query on those
#          top sequences that meet the thresholds established in this program and retrieves their RefSeqIDs and FASTA
#          sequences. The FASTA sequences are written to a separate file BRCA_family.fasta for input into MSA.py in
#          order to perform a MSA on the subject sequences with the query sequence of the BRCA2 gene in Homo sapiens.
########################################################################################################################


# Import modules
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Entrez


print("\nOpening the XML file to parse the BLAST results...\n")

# open the xml file, obtain the blast records.
results_handle = open("my_blast.xml")
blast_records = NCBIXML.read(results_handle)


def get_refseq_id(alignments, expect_threshold):
    """ Function isolates the refseq id of the top sequences that pass the threshold value parameter. """
    species_list = []                           # to avoid duplicate species
    refseq_list = ["NM_000059.4"]               # starting the list with the BRCA2 Human gene for MSA

    for aln in alignments:
        for hsp in aln.hsps:
            if hsp.expect == expect_threshold:             # previous knowledge that desired species have evalue = 0.0
                species = aln.title.split()[1:3]           # obtain the species name
                if species in species_list:
                    break
                else:
                    species_list.append(species)                   # avoid duplicates
                    refseq_id = aln.title.split()[0].split("|")[3] # obtain the refseq id from the title based on format
                    refseq_list.append(refseq_id)                  # add the refseq id to the list
                    break

    return refseq_list


# obtain the list of refseq ids
ref_seqs = get_refseq_id(blast_records.alignments, 0.0)


# get the fasta sequences of the ref seq IDs
def get_fasta(refseqids):
    """ Function takes in parameter of refseq IDs and performs an Entrez query to obtain their
        FASTA sequences. """
    for item in refseqids:
        Entrez.email = "eden.johnson@sjsu.edu"  # Always tell NCBI who you are
        handle = Entrez.efetch(db="nucleotide", id=item, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        sequence = str(record.seq)
        yield SeqRecord(Seq(sequence), id=record.id, description=record.description)


fasta_seqs = get_fasta(ref_seqs)
SeqIO.write(fasta_seqs, 'BRCA2_family.fasta', 'fasta')

print("\nFASTA sequences written to file. Beginning MSA now....\n\n")

# close the file
results_handle.close()
