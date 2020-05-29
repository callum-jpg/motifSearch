#!/usr/bin/env python
import wget
import os
from Bio import SeqIO
from Bio import Entrez



Entrez.email = "john.smith@john.com"

def download_dna(id, filename):
    """
    Takes an NCBI accession number from the nucleotide database and downloads the DNA
    in FASTA format. Saves as the given filename. Does not overwrite files with the same
    name. Saves to downloaded_DNA directory in working directory.

    >>> download_dna('NC_000913.3', 'ecoli_k12_gdna.fa')
    :param id:
    :param filename:
    :return:
    """
    save_dir = 'downloaded_DNA'
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)
    path = os.path.join(os.getcwd(), save_dir, filename)
    # Check if the same filename exists
    if not os.path.isfile(path):
        # Fetch ID from Entrez NCBI in fasta format as a handle
        # Handle is a wrapper around the text information retrived from Entrez/NCBI
        with Entrez.efetch(db="nucleotide", id=id, rettype="fasta", retmode="text") as handle:
            with open(path, "w") as output_handle:
                # .read() reads the entire handle. .readline() would, obviously, read line by line
                output_handle.write(handle.read())
            print('\'{}\' saved in {}'.format(filename, path))
    else:
        print('Filename \'{}\' already exists in {}'.format(filename, path))
# Test for download_dna
id = 'NC_000913.3'
filename = 'test.fa'
download_dna(id, filename)

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
    download_dna(v, filename)



# Check record IDs and sequence length in a FASTA
# Good for checking things have worked
def fasta_record_check(filename):
    """
    Reads a FASTA file and gives a quick summary of the contents
    """
    num_records = [record.id for record in SeqIO.parse(filename, 'fasta')]
    print('{} contains {} records'.format(filename, len(num_records)))
    for record in SeqIO.parse(filename, 'fasta'):
        print('{} is {}bp long'.format(record.id, len(record.seq)))

# Test for fasta_record_check
fasta_record_check('downloaded_DNA/test.fa')


# Combine multiple FASTA files into one
def concatenate_fasta(dna_dir, output_filename):
    """
    Reads the contents of a directory for FASTA files ending in -gdna.fa and combines them into a
    new .fa file in the same directory.
    """
    dirList = [i for i in os.listdir(dna_dir) if '-gdna.fa' in i]
    save_path = dna_dir + os.path.sep + output_filename
    # Empty list that will store all sequences
    all_seq = []
    if not os.path.isfile(save_path):
        for files in dirList:
            all_seq += [record for record in SeqIO.parse(dna_dir+os.path.sep+files, "fasta")]
        SeqIO.write(all_seq, save_path, 'fasta')
        # Count records in the file
        num_records = [record.id for record in SeqIO.parse(dna_dir+os.path.sep+output_filename, 'fasta')]
        print('{} saved in {} and contains {} records'.format(output_filename, save_path, len(num_records)))
    else:
        print('{} already exists'.format(output_filename))

# Test for concatenate_fasta
concatenate_fasta('downloaded_DNA', 'all-bac-seq.fa')


def rename_fasta_id(reference_dict, input_filename, output_filename):
    """
    Reads over a fasta file and replaces record.id

    Takes a dict containing the accession number (value) and the associated species name (key)
    and replaces the accession number (record.id) with the species name.
    """
    if not os.path.isfile(output_filename):
        with open(input_filename) as input, open(output_filename, 'w') as output:
            for record in SeqIO.parse(input, 'fasta'):
                # Check if there's any records not covered by the dictionary
                if record.id not in reference_dict.values():
                    print('Record ID \'{}\' not found in reference dictionary. It has been excluded.'.format(record.id))
                for k, v in reference_dict.items():
                    if record.id in v:
                        record.id = k
                        SeqIO.write(record, output, 'fasta')
        num_records = [record.id for record in SeqIO.parse(output_filename, 'fasta')]
        print('\'{}\' created and contains {} records'.format(output_filename, len(num_records)))
    else:
        print('\'{}\' already exists'.format(output_filename))

# Test for rename_fasta_id
bac_gdna_seqs = {
    'T.aquaticus': 'NZ_CP010822.1',
    'EPEC': 'LT827011.1',
    'E.coli(K12)': 'NC_000913.3',
    'M.tuberculosis': 'NC_000962.3',
    'P.aeruginosa': 'AE004091.2'
}
rename_fasta_id(bac_gdna_seqs, 'downloaded_DNA/all-bac-seq.fa', 'downloaded_DNA/test.fa')


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

