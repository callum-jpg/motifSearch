#!/usr/bin/env python
import functions as fc
from script import IUPAC_dict
import wget
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import pandas as pd
import numpy as np

fc.num_combinations('TN')
fc.percentage(50, 100)
fc.max_iupac('N')

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
    print("Record " + record.id + ", length " + str(len(record.seq)) +
          "First/last 10 bases: " + record.seq[:10] + "..." + record.seq[-10:])

# Returns record referring to chromosomes, Y, and mitochondrial DNA
for record in SeqIO.parse('GRCh38_latest_genomic.fna', 'fasta'):
    if 'NC_' in record.id:
        print record.id


## Counting motifs
# This is inaccuate - only counts TTTC motifs
for record in SeqIO.parse(input_file, 'fasta', alphabet = IUPAC.unambiguous_dna):
    sequence = record.seq
    rev_sequence = record.seq.reverse_complement()
    ## Sequence.count must be converted to string for print
    print(record.seq.alphabet)
    print("Record " + record.id + " has a sequence length of " + str(len(record.seq)) +
          "bp. \nIt contains " +
          str(sequence.count('TTTC' or
          'TGTC' or
          'TCTC' or
          'TATC')) + " TNTC motifs and " +
          str(rev_sequence.count('TTTC' or
          'TGTC' or
          'TCTC' or
          'TATC')) + " TNTC motifs on the reverse complement.")

## Counting occurence of variations of TNTC motifs
for record in SeqIO.parse(input_file, 'fasta'):
    sequence = record.seq
    rev_sequence = record.seq.reverse_complement()
    modified = (rev_sequence.count('TTTC') +
                sequence.count('TGTC') +
                sequence.count('TCTC') +
                sequence.count('TATC') +
                rev_sequence.count('TTTC') +
                rev_sequence.count('TGTC') +
                rev_sequence.count('TCTC') +
                rev_sequence.count('TATC'))

    ## Sequence.count must be converted to string for print
    print("Record " + record.id + " is " + str(len(record.seq)) + "bp long.\n"
          " It contains " +
          str(sequence.count('TTTC')) + ' TTTC motifs, ' +
          str(sequence.count('TGTC')) + ' TGTC motifs, ' +
          str(sequence.count('TCTC')) + ' TCTC motifs, and ' +
          str(sequence.count('TATC')) + ' TATC motifs.\n' +
          str(rev_sequence.count('TTTC')) + ' TTTC motifs, ' +
          str(rev_sequence.count('TGTC')) + ' TGTC motifs, ' +
          str(rev_sequence.count('TCTC')) + ' TCTC motifs, and ' +
          str(rev_sequence.count('TATC')) + ' TATC motifs on the reverse complement')

    # 1.5% of E. coli K12 genome can be modified by Taq DarT
    ((modified / float(len(sequence) * 2)) * 100)





# Motif search test
test_seq = Seq("ATAGCTCTAGCTATGCTACGATACGTGTC")
test_seq2 = Seq("TGTCTGTCTGTC")



motif_data = {}
motifs = ['TTTC', 'TGTC', 'TCTC', 'TATC']
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

IUPAC = {'A': 'A',
         'T': 'T',
         'G': 'G',
         'C': 'C',
         'N': ['A', 'T', 'G', 'C'],
         'R': ['A', 'G'],
         'Y': ['C', 'T'],
         'S': ['G', 'C'],
         'W': ['A', 'T'],
         'K': ['G', 'T'],
         'M': ['A', 'C'],
         'B': ['C', 'G', 'T'],
         'D': ['A', 'G', 'T'],
         'H': ['A', 'C', 'T'],
         'V': ['A', 'C', 'G']
         }

def iupac_motif(motif_string):
    """
    Takes a motif in IUPAC format and generates a list of possible strings

    >>> IUPAC_motif("TNTC")
    >>> ['TTTC', 'TATC', 'TGTC', 'TCTC']
    """





# Find associated values with key
IUPAC['N']
# Find length of values for key 'N' - decides how many motifs to generate
len(IUPAC['N'])
# Print values for a specific key
for i in IUPAC['N']:
    print i
# Generate unique motifs

### Calculate  possible motifs for TN (4 total):

# First, calculate possible combinations:
fc.num_combinations('TBV')
fc.max_iupac('TBV')

tmotif = 'TN'
tlist = []
for i in range(0, fc.num_combinations(tmotif)):
    for j in range(1, fc.max_iupac(tmotif)):
        print j
        tlist.append(IUPAC_dict[tmotif[0]] + IUPAC_dict[tmotif[1]])


# If using append, generate the full motif within the append function
# This will require the IUPAC_dict.value position
# This will require any values below max(IUPAC_dict.values()) to be repeated
tlist.append("T" + "T")


# np.empty creates shape but doesn't use zeros - uses arbitrary data
# Scrapping the idea of saving motifs in an array - will instead store in a list
# For visualisation of all motifs, list could be coverted to array later.
np.empty((fc.num_combinations('TBV')/fc.max_iupac('TBV'), fc.max_iupac('TBV')))

np.empty((10, 10))




for i in 'TN':
    if len(IUPAC_dict[i]) > 1:
        print IUPAC_dict[i]

max([len(IUPAC_dict[i]) for i in 'TN'])




max([1, 2, 3, 4])

examp = (['TTTC', 'TATC', 'TCTC', 'TGTC'])
a = np.array([[1, 2, 3], [4, 5, 6]])
# i, j[start:end, start:end]
a[:2, 1:3]
a
# Calculate possible motifs for TNN (16 total):

a
# i, j[start:end, start:end]
a[:2, 1:3]
# Pull 0th column
a[:, 0]
# Pull index 1 row
a[1, :]







