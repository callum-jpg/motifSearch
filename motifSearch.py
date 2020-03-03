#!/usr/bin/env python
import functions as fc
import wget
from Bio import Seq
from Bio import SeqIO
import pandas as pd


IUPAC_dict = Seq.IUPAC.IUPACData.ambiguous_dna_values
fc.percentage(50, 100)
fc.iupac_to_dna('TNTM')
fc.motifs_in_fasta(input_file, 'TNTC')


# Ensembl GRCh38 human genome. Chromosome 22 only (for faster testing)
wget.download('http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa')
# Download E. coli K12 gDNA
wget.download('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz')

input_file = 'test-seq.fna'

fc.motifs_in_fasta(input_file, 'TNTC')
output = pd.DataFrame(data=motifs_in_fasta(input_file, 'TNTC'))
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
