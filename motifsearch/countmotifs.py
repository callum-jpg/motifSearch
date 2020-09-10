from itertools import product
import os
from Bio import Seq
from Bio import SeqIO
from Bio import Entrez

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
                   "motif seq",
                   "record length"]
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
        motif_data[output_cols[5]].append(str(motif))
        motif_data[output_cols[6]].append(len(sequence))


    return motif_data

## motifs_in_fasta test
# mots = ['TNTC', 'TNTM', 'TYT', 'TTTH']
# output = pd.DataFrame()
# for i, _ in enumerate(mots):
# 	print("Counting sites for {0} motif".format(mots[i]))
# 	df = pd.DataFrame(data=fc.motifs_in_fasta('downloaded_DNA/all-bac-renamed.fa', str(mots[i])))
# 	output = output.append(df)
# output.to_csv('output.csv')



### Gathering DNA functions

Entrez.email = "john.smith@john.com"


def download_dna(dict, filename):
    '''
    Takes an NCBI accession number from the assembly database and downloads the gDNA
    in FASTA format. Input Sequences should be found in the values of a dictionary. 
    Multiple items in a dictionary will yield a file with multiple DNA sequences. Does
    not overwrite files. Saves to 'downloaded_DNA' in the current working directory.
    
    Saves as filename input. Remeber to include the desired path and .fa extention.
    
    >>> dna_dict = {
        'E.coli(K12)': 'NC_000913.3',                    
        }
    >>> download_dna(dna_dict, 'downloaded_DNA/ecoli_k12_gdna.fa')

    '''
    save_dir = 'downloaded_DNA'
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)
    path = os.path.join(os.getcwd(), save_dir, filename)
    # Check if the same filename exists
    if not os.path.isfile(path):
        # Create empty list to store DNA
        handle = []
        for k, v in dict.items():
            # Fetch DNA and append to list. Concatenates DNA sequences
            handle.append(Entrez.efetch(db="nucleotide", id=v, rettype="fasta", retmode="text"))
        
        # Save elements in the DNA list to a single file
        with open(path, 'w') as output:
            for record in handle:
                output.write(record.read())         
    else:
        print('File \'{}\' already exists in \'{}\''.format(filename, path))

    
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
#fasta_record_check('downloaded_DNA/test.fa')


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
                        # Strip spaces so biopython can identify by id
                        record.id = k.replace(' ', '')
                        SeqIO.write(record, output, 'fasta')
        num_records = [record.id for record in SeqIO.parse(output_filename, 'fasta')]
        print('\'{}\' created and contains {} records'.format(output_filename, len(num_records)))
    else:
        print('\'{}\' already exists'.format(output_filename))

# Test for rename_fasta_id
# bac_gdna_seqs = {
#     'T.aquaticus': 'NZ_CP010822.1',
#     'EPEC': 'LT827011.1',
#     'E.coli(K12)': 'NC_000913.3',
#     'M.tuberculosis': 'NC_000962.3',
#     'P.aeruginosa': 'AE004091.2'
# }
#rename_fasta_id(bac_gdna_seqs, 'downloaded_DNA/all-bac-seq.fa', 'downloaded_DNA/test.fa')
        
        
        
        
        
        
        
        