#!/usr/bin/env python
import functions as fc
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

#fc.percentage(50, 100)
#fc.iupac_to_dna('TNTM')

output = pd.DataFrame(data=fc.motifs_in_fasta('chr1-23+Y.fna', 'TNTC'))
output.to_csv('output.csv')


# Unreliable way to slice out chromosome numbers
input_seq_iterator = SeqIO.parse('chr1-23+Y.fna', 'fasta')
records = [record.id for record in input_seq_iterator]
chr_name = [records[i][7:9] for i, j in enumerate(records)]
output['Chromosome'] = chr_name # Append new column with chromosome numbers

plt.bar(output['Chromosome'].values, output['Perc DNA modified (total)'].values)