from itertools import product

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