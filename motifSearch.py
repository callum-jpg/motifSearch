#!/usr/bin/env python
import wget
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import pandas


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
          'TATC')) + " TNTC motifs on the reverse complement."
          )

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



motifs = ['TTTC', 'TGTC', 'TCTC', 'TATC']
for record in SeqIO.parse(input_file, 'fasta'):
    sequence = record.seq
    rev_sequence = record.seq.reverse_complement()
    total_sites= 0
    total_sites_rev = 0
    for i in range(0, len(motifs)):
        print(str(sequence.count(motifs[i])) + " " + str(motifs[i]))
        total_sites += sequence.count(motifs[i])
        print(str(rev_sequence.count(motifs[i])) + " " + str(motifs[i]))
        total_sites_rev += rev_sequence.count(motifs[i])

    print(str(total_sites) + " sites on complement")
    print(str(total_sites_rev) + " sites on reverse complement")
    print(str(total_sites + total_sites_rev) + " total sites for complement + reverse")
    print(str(percentage((total_sites + total_sites_rev), len(sequence + rev_sequence))) +
          "%% of genome can be modified")




