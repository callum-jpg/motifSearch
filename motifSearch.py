#!/usr/bin/env python
import matplotlib
matplotlib.use("TkAgg")
import functions as fc
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# For realoding
import importlib
importlib.reload(fc)

#fc.percentage(50, 100)
#fc.iupac_to_dna('TNTM')

# Human gDNA
output = pd.DataFrame(data=fc.motifs_in_fasta('chr1-23+Y.fna', 'TNTC'))
output.to_csv('output.csv')
# Unreliable way to slice out chromosome numbers
input_seq_iterator = SeqIO.parse('chr1-23+Y.fna', 'fasta')
records = [record.id for record in input_seq_iterator]
chr_name = [records[i][7:9] for i, j in enumerate(records)]
output['Chromosome'] = chr_name # Append new column with chromosome numbers
plt.bar(output['Chromosome'].values, output['Perc DNA modified (total)'].values)




### Bacterial gDNA
output = pd.DataFrame(data=fc.motifs_in_fasta('all-bac-corrected.fna', 'TNTC'))
output.to_csv('output.csv')

plt.bar(output['Record'].values, output['Perc DNA modified (total)'].values)
plt.suptitle(output['motif seq'][1], fontsize=25)


#werks
mots = ['TNTC', 'TNTM', 'TYT', 'TTTH']

fig, axes = plt.subplots(nrows=2, ncols=2, sharey=True)
fig.subplots_adjust(hspace=0.5)
#fig.setp(axes.xaxis.get_majorticklabels(), rotation=45)
for i, ax in enumerate(axes.flat):
	output = pd.DataFrame(data=fc.motifs_in_fasta('all-bac-corrected.fna', str(mots[i])))
	ax.bar(output['Record'].values, output['Perc DNA modified (total)'].values)
	ax.set(title=mots[i])

	print(i, ax)






# Calculating position of each DNA motif
test_seq = 'GATACGATTAAAAATCTCAAAAGCACCAGATCGA'
test_seq.find('TCTC') # only finds the first occurance

# List comprehension to find index of occurances
T = [i for i, j in enumerate('TAAT') if j == 'T']
A = [i for i, j in enumerate('TAAT') if j == 'A']


# Add sequence to object
seq = ''
for record in SeqIO.parse('k12-gdna.fna', 'fasta'):
	seq = record.seq
# Make a smaller, faster processing sequence
seq_sm = seq[:int(((len(seq))/1000))]
len(seq_sm)

# Count positions of a single letter
test_dict = {'A': [i for i, j in enumerate(seq_sm) if j == 'A'],
             'T': [i for i, j in enumerate(seq_sm) if j == 'T'],
             'G': [i for i, j in enumerate(seq_sm) if j == 'G'],
             'C': [i for i, j in enumerate(seq_sm) if j == 'C']
             }

# OK-ish plot. Unsure if at all useful
for k, v in test_dict.items():
	plt.xlim(0, len(seq_sm))
	sns.kdeplot(test_dict[k], bw=0.02, label=k)


# Make 4x subplots for each nucleotide
# fig, axes = plt.subplots(4, 1)
# for i, (k, v) in enumerate(test_dict.items()):
# 	sns.kdeplot(test_dict[k],
# 	            bw=0.02,
# 	            ax=i+1
# 	            )
# plt.show()


fig = plt.figure()
for i, (k, v) in enumerate(test_dict.items()):
	ax = fig.add_subplot(4, 1, i+1)
	plt.setp(ax, xlim=(0, len(seq_sm)), ylim=(0,0.0004))
	sns.kdeplot(test_dict[k], bw=0.02, ax=ax)
	ax.set_ylabel(k, fontsize=20)
#plt.xlim(0, len(seq_sm))
plt.xlabel('Position (kb)', fontsize=20)
plt.suptitle('Nucleotide distribution in K12 E. coli', fontsize=25)
plt.show