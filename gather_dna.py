#!/usr/bin/env python
import wget
from Bio import SeqIO
import functions as fc

# For realoding
import importlib
importlib.reload(fc)





# Bacterial gDNA of interest with RefSeq sequence
# Gather accession numbers from NCBI 'assembly' database
# For gDNA, select RefSeq that corresponds with the chromosome
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

bac_gdna_seq_extra = {
    **bac_gdna_seqs, # ** unpacks dict and adds key: value pairs and merges to dictionary
    'D.compactum': 'NZ_AP018316.1',
    'S.enterica': 'NC_003197.2',
    'G.thermoleovorans': 'NZ_CP014335.1',
    'B.naejangsanensis': 'NZ_CP015614.1',
    'B.cereus': 'NZ_CP034551.1',
    'C.abyssi': 'NZ_CP018099.1',
    'B.diminuta': 'NZ_CP035093.1',
    'V.cholerae': 'NZ_CP010812.1',
    'N.canis': 'NZ_LR134313.1',
    'T.scotoductus': 'NC_014974.1',
    'G.violaceus': 'NC_005125.1',
    'S.cyanosphaera': 'NC_019748.1',
    'L.trevisanii': 'NZ_AP019840.1',
    'A.giovannonii': 'NZ_CP042997.1',
    'P.bacterium': 'NZ_CP036264.1',
    'S.cyanobacteriorum': 'NZ_CP019633.1',
    'P.mikurensis': 'NC_017080.1',
    'R.baltica': 'NC_005027.1',
    'L.pneumophila': 'NZ_LR134380.1'
}
#len(bac_gdna_seq_extra)

# Pathogens
# 'A.baumannii': 'NZ_CP046654.1',
# 'K.pneumoniae': 'NC_016845.1',
# 'E.hormaechei': 'NZ_CP017179.1'

# Download a dictionary of dna
for k, v in bac_gdna_seq_extra.items():
    filename = k + '-gdna.fa'
    fc.download_dna(v, filename)


# Test for concatenate_fasta
fc.concatenate_fasta('downloaded_DNA', 'all-bac-seq.fa')

fc.concatenate_fasta('downloaded_DNA', 'all-bac-seq-extra.fa')



# Test for rename_fasta_id
bac_gdna_seqs = {
    'T.aquaticus': 'NZ_CP010822.1',
    'EPEC': 'LT827011.1',
    'E.coli(K12)': 'NC_000913.3',
    'M.tuberculosis': 'NC_000962.3',
    'P.aeruginosa': 'AE004091.2'
}
fc.rename_fasta_id(bac_gdna_seqs, 'downloaded_DNA/all-bac-seq.fa', 'downloaded_DNA/all-bac-renamed.fa')

fc.rename_fasta_id(bac_gdna_seq_extra, 'downloaded_DNA/all-bac-seq-extra.fa', 'downloaded_DNA/all-bac-extra-renamed.fa')

fc.fasta_record_check('downloaded_DNA/all-bac-seq.fa')
fc.fasta_record_check('downloaded_DNA/all-bac-renamed.fa')

fc.fasta_record_check('downloaded_DNA/all-bac-seq-extra.fa')
fc.fasta_record_check('downloaded_DNA/all-bac-extra-renamed.fa')

# Human gDNA
## OLD CODE ##

# Ensembl GRCh38 human genome. Chromosome 22 only (for faster testing)
wget.download('http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa')

# From human gDNA, move chr1:23 + Y into separate fasta
input_seq_iterator = SeqIO.parse('GRCh38_latest_genomic.fna', 'fasta')
chromosomes = [record for record in input_seq_iterator
               if 'NC_00' in record.id]
SeqIO.write(chromosomes, "chr1-23+Y.fna", "fasta")

fc.fasta_record_check('chr1-23+Y.fna')

human_chr = {
    "chr1": "NC_000001.11",
    "chr2": "NC_000002.12",
    "chr3": "NC_000003.12",
    "chr4": "NC_000004.12",
    "chr5": "NC_000005.10",
    "chr6": "NC_000006.12",
    "chr7": "NC_000007.14",
    "chr8": "NC_000008.11",
    "chr9": "NC_000009.12",
    "chr10": "NC_000010.11",
    "chr11": "NC_000011.10",
    "chr12": "NC_000012.12",
    "chr13": "NC_000013.11",
    "chr14": "NC_000014.9",
    "chr15": "NC_000015.10",
    "chr16": "NC_000016.10",
    "chr17": "NC_000017.11",
    "chr18": "NC_000018.10",
    "chr19": "NC_000019.10",
    "chr20": "NC_000020.11",
    "chr21": "NC_000021.9",
    "chr22": "NC_000022.11",
    "chrX": "NC_000023.11",
    "chrY": "NC_000024.10"
}

fc.rename_fasta_id(human_chr, 'chr1-23+Y.fna', "downloaded_DNA/human_gdna_renamed.fa")

fc.fasta_record_check('downloaded_DNA/human_gdna_renamed.fa')


# # Extract just the gDNA based on record.id
# def extract_gdna(input_file, gdna_record, output_file):
#     """Takes an input FASTA file, extracts the desired gDNA containing
#     record and outputs as a new FASTA file"""
#     input_iter = SeqIO.parse(input_file, 'fasta')
#     gdna = [record for record in input_iter
#             if gdna_record in record.id]
#     SeqIO.write(gdna, output_file, 'fasta')
#
# extract_gdna('GCF_001399775.1_ASM139977v1_genomic.fna', 'NZ_CP010822', 'taq-gdna.fna')
#

