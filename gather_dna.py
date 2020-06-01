#!/usr/bin/env python
import wget
from Bio import SeqIO
import functions as fc

# For realoding
import importlib
importlib.reload(fc)





# Bacterial gDNA of interest with accession number
bac_gdna_seqs = {
    'T.aquaticus': 'NZ_CP010822.1',
    'EPEC': 'LT827011.1',
    'E.coli(K12)': 'NC_000913.3',
    'M.tuberculosis': 'NC_000962.3',
    'P.aeruginosa': 'AE004091.2'
}
# Download a dictionary of dna
for k, v in bac_gdna_seqs.items():
    filename = k + '-gdna.fa'
    fc.download_dna(v, filename)


# Test for concatenate_fasta
fc.concatenate_fasta('downloaded_DNA', 'all-bac-seq.fa')



# Test for rename_fasta_id
bac_gdna_seqs = {
    'T.aquaticus': 'NZ_CP010822.1',
    'EPEC': 'LT827011.1',
    'E.coli(K12)': 'NC_000913.3',
    'M.tuberculosis': 'NC_000962.3',
    'P.aeruginosa': 'AE004091.2'
}
fc.rename_fasta_id(bac_gdna_seqs, 'downloaded_DNA/all-bac-seq.fa', 'downloaded_DNA/all-bac-renamed.fa')

fc.fasta_record_check('downloaded_DNA/all-bac-seq.fa')
fc.fasta_record_check('downloaded_DNA/all-bac-renamed.fa')





## OLD CODE ##

# Ensembl GRCh38 human genome. Chromosome 22 only (for faster testing)
wget.download('http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa')

# From human gDNA, move chr1:23 + Y into separate fasta
input_seq_iterator = SeqIO.parse('GRCh38_latest_genomic.fna', 'fasta')
chromosomes = [record for record in input_seq_iterator
               if 'NC_00' in record.id]
SeqIO.write(chromosomes, "chr1-23+Y.fna", "fasta")

input_seq_iterator = SeqIO.parse('chr1-23+Y.fna', 'fasta')
records = [record.id for record in input_seq_iterator]
chr_name = [records[i][7:9] for i, j in enumerate(records)]


# Extract just the gDNA based on record.id
def extract_gdna(input_file, gdna_record, output_file):
    """Takes an input FASTA file, extracts the desired gDNA containing
    record and outputs as a new FASTA file"""
    input_iter = SeqIO.parse(input_file, 'fasta')
    gdna = [record for record in input_iter
            if gdna_record in record.id]
    SeqIO.write(gdna, output_file, 'fasta')

extract_gdna('GCF_001399775.1_ASM139977v1_genomic.fna', 'NZ_CP010822', 'taq-gdna.fna')

# Download Mtb gDNA
extract_gdna('GCF_000195955.2_ASM19595v2_genomic.fna', 'NC_000962.3', 'mtb-gdna.fna')

