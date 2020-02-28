from script import IUPAC_dict

def percentage(percent, whole):
    '''Calculates the percentage with float#
    >>> percentage(50, 100)
    >>> 50.0
    '''
    return ((percent/float(whole)) * 100.0)


def num_combinations(motif):
    """Returns the number of possible combinations given a motif in
    IUPAC format

    >>> num_combinations('TNTC')
    >>> 4

    >>> num_combinations('TNNN')
    >>> 64
    """
    combi = 1
    for i in motif:
        combi = combi * len(IUPAC_dict[i])
    return combi

def max_iupac(motif):
    """
    Returns the maximum number of nucleotides that an IUPAC denotes

    >>> max_iupac('N')
    >>> 4
    """
    return max([len(IUPAC_dict[i]) for i in motif])