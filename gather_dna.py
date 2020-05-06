#!/usr/bin/env python
import wget
from Bio import SeqIO

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


# Download E. coli K12 gDNA
wget.download('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz')

# Download T.aq
wget.download('ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Thermus_aquaticus/representative/GCF_001399775.1_ASM139977v1/GCF_001399775.1_ASM139977v1_genomic.fna.gz')

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






# Check sequence length for records in fasta
# Make into a function?
for record in SeqIO.parse('all-bac-corrected.fna', 'fasta'):
    print(
        #record.id, len(record.seq),
        record.description)





# Trying to find a better way of getting gDNA data
# To find FTP sequences, look at assembly_summary_genbank.txt
from ftplib import FTP

host = "ftp.ncbi.nlm.nih.gov"
wdir = "/genomes/refseq/"

# Download T.aq gDNA
ftp = FTP(host, 'anonymous', 'anonymous')
taqdir = '/genomes/refseq/bacteria/Thermus_aquaticus/representative/GCF_001399775.1_ASM139977v1/GCF_001399775.1_ASM139977v1_genomic.fna.gz'
list = ftp.nlst(wdir+taqdir)
ftp.quit()
for name in list1:
    #print(name)
    if 'Thermus' in name:
        print(name)
wget.download('ftp://'+host+wdir+taqdir)

# mtb
ftp = FTP(host, 'anonymous', 'anonymous')
mtbdir = 'bacteria/Mycobacterium_tuberculosis/reference/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz'
list = ftp.nlst(wdir+mtbdir)
ftp.quit()
for i in list:
    #print(i)
    if 'fna' in i:
       print(i)
wget.download('ftp://'+host+wdir+mtbdir)

# EPEC
ftp = FTP(host, 'anonymous', 'anonymous')
epecdir = '/genomes/all/GCA/900/149/915/GCA_900149915.1_EPEC-E2348_69-V2/GCA_900149915.1_EPEC-E2348_69-V2_genomic.fna.gz'
list = ftp.nlst('/genomes/all/GCA/900/149/915/GCA_900149915.1_EPEC-E2348_69-V2/GCA_900149915.1_EPEC-E2348_69-V2_genomic.fna.gz')
ftp.quit()
for i in list:
    print(i)
    if 'fna' in i:
       print(i)
wget.download('ftp://'+host+epecdir)
extract_gdna('GCA_900149915.1_EPEC-E2348_69-V2_genomic.fna', 'LT827011.1', 'epec-gdna.fna')

# Pseu
ftp = FTP(host, 'anonymous', 'anonymous')
pseudir = '/genomes/all/GCA/000/006/765/GCA_000006765.1_ASM676v1/GCA_000006765.1_ASM676v1_genomic.fna.gz'
list = ftp.nlst(pseudir)
ftp.quit()
for i in list:
    #print(i)
    if 'fna' in i:
       print(i)
wget.download('ftp://'+host+pseudir)
extract_gdna('GCA_000006765.1_ASM676v1_genomic.fna', 'AE004091.2', 'pseu-gdna.fna')

# assembly_summary_genbank.txt


# Combine gdna files
from os import listdir
dirlist = [i for i in listdir()
           if '-gdna.fna' in i]
print(dirlist)
sequences = []
for files in dirlist:
        sequences += [record for record in SeqIO.parse(files, 'fasta')]
SeqIO.write(sequences, 'all-bac-seq.fna', 'fasta')

# Renaming ID's to that of bacteria used
gb_to_sp = {'AE004091.2': 'P.aeruginosa',
            'LT827011.1': 'EPEC',
            'NC_000913.3': 'E.coli(K12)',
            'NC_000962.3': 'M.tuberculosis',
            'NZ_CP010822.1': 'T.aquaticus'}

with open('all-bac-seq.fna') as input, open('all-bac-corrected.fna', 'w') as output:
    for record in SeqIO.parse(input, 'fasta'):
        #print(record.description)
        if record.id in gb_to_sp.keys():
            #print(gb_to_sp[record.id])
            record.id = gb_to_sp[record.id]
            #record.decription = 'test'
        print('NEW', record.id)
        SeqIO.write(record, output, 'fasta')

# gb_to_sp['AE004091.2']
# 'P. aeruginosa'

