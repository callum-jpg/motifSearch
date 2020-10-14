#!/usr/bin/env python3
from motifsearch import countmotifs as ms

#%% Downloading bacterial gDNA fasta sequences

# Bacterial gDNA of interest with RefSeq sequence
# Gather accession numbers from NCBI 'assembly' database
# For gDNA, select RefSeq that corresponds with the chromosome
bac_gdna_seqs = {
    'E. coli': 'NC_000913.3',
    'T. aquaticus': 'NZ_CP010822.1',
    'M. tuberculosis': 'NC_000962.3',
    'P. aeruginosa': 'AE004091.2',
    'S. aureus': 'NC_007795.1',
    'S. enterica': 'NC_003197.2'
}

# Download dictionary of DNA sequences from NCBI
ms.download_dna(bac_gdna_seqs, 'bacterial-gDNA.fa')

# Rename the DNA sequence FASTA headers to match the key in the dictionary
ms.rename_fasta_id(bac_gdna_seqs, 'downloaded_DNA/bacterial-gDNA.fa', 'downloaded_DNA/bacterial-gDNA-renamed.fa')
    

#%% Checking the contents of the concatenated fasta file.

ms.fasta_record_check('downloaded_DNA/bacterial-gDNA.fa')
ms.fasta_record_check('downloaded_DNA/bacterial-gDNA-renamed.fa')


#%% Downloading human gDNA - autosomes only

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
}

# Download dictionary of human chromosome sequences from NCBI
ms.download_dna(human_chr, 'human-gDNA.fa')
    
ms.rename_fasta_id(human_chr, 'downloaded_DNA/human-gDNA.fa', 
                   'downloaded_DNA/human-gDNA-renamed.fa')

#%% Checking the contents of the concatenated fasta file.

ms.fasta_record_check('downloaded_DNA/human-gDNA.fa')

ms.fasta_record_check('downloaded_DNA/human-gDNA-renamed.fa')











