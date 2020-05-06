from itertools import product
from Bio import Seq
from Bio import SeqIO

IUPAC_dict = Seq.IUPAC.IUPACData.ambiguous_dna_values

def percentage(percent, whole):
    '''Calculates the percentage with float#
    >>> percentage(50, 100)
    >>> 50.0
    '''
    return ((percent/float(whole)) * 100.0)


def iupac_to_dna(sequence):
    """
    Takes an input DNA sequence in IUPAC format and returns all possible sequences.

    >>> iupac_to_dna('TN')
    >>> ['TA', 'TT', 'TG', 'TC']

    >>> iupac_to_dna('TNTM')
    >>> ['TATA', 'TATC', 'TTTA', 'TTTC', 'TGTA', 'TGTC', 'TCTA', 'TCTC']
    """
    # Cartesian product will generate all of the possible ordered motif strings from ambiguous IUPAC nucleotides
    # Product example: product(A, B) returns the same as ((x,y) for x in A for y in B)

    # [IUPAC_dict[j] for j in tmotif] creates an ordered 2d array (ordered, as in, the same order as the input motif)
    # of all possible bases that occur at each position of the sequence

    # The asterix (*) before IUPAC_dict is used to unpack the argument list into individual lists

    # "".join(i) takes each cartesian product ('x', 'y', 'z') from the cartesian output list
    # and collapses it into a single string ('xyz') on a item by item bases (hence the for i...)
    return ["".join(i) for i in product(*[IUPAC_dict[j] for j in sequence])]


def motifs_in_fasta(fasta, motif):
    """Returns a dictionary of the number of motifs found within each FASTA sequence identified.
    Parses through all available FASTA sequences. Identifies by name after >.
    Identifies motifs on both complement and reverse complement.

    fasta is a file in FASTA format
    motif is a short sequence in IUPAC format
    """
    motif_data = {} # Create empty dictionary to hold motif counts
    motif_variations = iupac_to_dna(motif)

    # Create empty lists in dict to store values
    output_cols = ["Record",
                   "Total complement sites",
                   "Total reverse complement sites",
                   "Total number of sites",
                   "Perc DNA modified (total)",
                   "motif seq"]
    # Populate dict with empty lists
    for i in output_cols:
        motif_data[i] = []

    for record in SeqIO.parse(fasta, 'fasta'):
        sequence = record.seq # Read the sequence associated with record
        rev_sequence = record.seq.reverse_complement() # Get the reverse complement
        total_sites = 0 # Counter for complement motifs
        total_sites_rev = 0 # counter for rev complement motifs

        for j in range(0, len(motif_variations)): # For every motif, count the number of occurrences
            total_sites += sequence.count(motif_variations[j])
            total_sites_rev += rev_sequence.count(motif_variations[j])

        motif_data[output_cols[0]].append(record.id) # Create a new dict entry to hold FASTA record name

        # Add counts to dict
        motif_data[output_cols[1]].append(total_sites)
        motif_data[output_cols[2]].append(total_sites_rev)
        motif_data[output_cols[3]].append((total_sites + total_sites_rev))
        motif_data[output_cols[4]].append(percentage((total_sites + total_sites_rev),
                                                        (len(sequence) + len(rev_sequence))))
        motif_data[output_cols[6]].append(str(motif))


    return motif_data
