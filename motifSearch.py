#!/usr/bin/env python
import functions as fc
from Bio import Seq
from Bio import SeqIO
import pandas as pd


IUPAC_dict = Seq.IUPAC.IUPACData.ambiguous_dna_values
fc.percentage(50, 100)
fc.iupac_to_dna('TNTM')


input_file = 'test-seq.fna'

fc.motifs_in_fasta('test-seq.fna', 'TNTC')
output = pd.DataFrame(data=fc.motifs_in_fasta(input_file, 'TNTC'))
output.to_csv('output.csv')

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
