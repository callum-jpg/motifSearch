#!/usr/bin/env python
import wget
from Bio import SeqIO

# Ensembl GRCh38 human genome. Chromosome 22 only (for faster testing)
wget.download('http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa')
# Download E. coli K12 gDNA
wget.download('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz')


# From human gDNA, move chr1:23 + Y into separate fasta
input_seq_iterator = SeqIO.parse('GRCh38_latest_genomic.fna', 'fasta')
chromosomes = [record for record in input_seq_iterator
               if 'NC_00' in record.id]
SeqIO.write(chromosomes, "chr1-23+Y.fna", "fasta")