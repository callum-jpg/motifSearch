import wget
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Seq import generic_dna


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




## Counting motifs
for record in SeqIO.parse(input_file, 'fasta'):
    sequence = record.seq
    rev_sequence = record.seq.reverse_complement()
    ## Sequence.count must be converted to string for print
    print("Record " + record.id + " has a sequence length of " + str(len(record.seq)) +
          ". \nIt contains " +
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
    ## Sequence.count must be converted to string for print
    print("Record " + record.id + " contains " +
          str(sequence.count('TTTC')) + ' TTTC motifs, ' +
          str(sequence.count('TGTC')) + ' TGTC motifs, ' +
          str(sequence.count('TCTC')) + ' TCTC motifs, and ' +
          str(sequence.count('TATC')) + ' TATC motifs.')


# Motif search test
test_seq = Seq("ATAGCTCTAGCTATGCTACGATACGTGTC", generic_dna)

print(str(test_seq.count('TTTC')) + ' TTTC motifs, ' +
      str(test_seq.count('TGTC')) + ' TGTC motifs, ' +
      str(test_seq.count('TCTC')) + ' TCTC motifs, and ' +
      str(test_seq.count('TATC')) + ' TATC motifs.')



for index, record in enumerate(input_dna):
    print("index %i, ID = %s, length %i, with %i features"
          % (index, record.id, len(record.seq), len(record.features)))









