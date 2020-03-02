#!/usr/bin/env python
import functions as fc
import wget
from Bio import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import pandas as pd
import numpy as np


IUPAC_dict = Seq.IUPAC.IUPACData.ambiguous_dna_values
fc.percentage(50, 100)
fc.iupac_to_dna('TNTM')

# Ensembl GRCh38 human genome. Chromosome 22 only (for faster testing)
wget.download('http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa')
# Download E. coli K12 gDNA
wget.download('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz')

input_file = 'GCF_000005845.2_ASM584v2_genomic.fna'

# SeqIO.parse creates a 'SeqRecord' object takes the record identifier '>' and returns
# information for that record, such as sequence length and any other information.
# Record.id is the first string after '>' and anything after a space is ignored.
# Record.seq is the sequence of that record

for record in SeqIO.parse(input_file, 'fasta'):
    print("Record " + record.id + ", length " + str(len(record.seq)))

# Returns record referring to 23 chromosomes, Y, and mitochondrial DNA
for record in SeqIO.parse('GRCh38_latest_genomic.fna', 'fasta'):
    if 'NC_' in record.id:
        print record.id

motif_data = {}
motifs = fc.iupac_to_dna('TNTC')
for record in SeqIO.parse(input_file, 'fasta'):
    sequence = record.seq
    rev_sequence = record.seq.reverse_complement()
    total_sites= 0
    total_sites_rev = 0
    motif_data['Record'] = record.id
    for i in range(0, len(motifs)):
        # Build up total number of sites over motif interations
        total_sites += sequence.count(motifs[i])
        total_sites_rev += rev_sequence.count(motifs[i])
        # Add individual counts for each motif to data dict
        motif_data[motifs[i] + " complement sites"] = sequence.count(motifs[i])
        motif_data[motifs[i] + " reverse complement sites"] = rev_sequence.count(motifs[i])

    # At end of motif search, build dictionary with total sites
    motif_data["Total complement sites"] = total_sites
    motif_data["Total reverse complement sites"] = total_sites_rev
    motif_data["Total number of sites"] = total_sites + total_sites_rev

# Data dict output
output = pd.DataFrame(data=motif_data, index=[0])
output.to_csv('output.csv')

t = {}
t = {'record' : ['bacterial gDNA', 'human gDNA'], motifs[0] : [1, 200]}
pd.DataFrame(data=t)






