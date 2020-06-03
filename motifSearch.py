#!/usr/bin/env python
import matplotlib
matplotlib.use("TkAgg") # Backend to use for PyCharm interactive plots
matplotlib.use("GTK3Agg") # Backend to use for savefig with properly scaled DPI

import functions as fc
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

plt.plot([1,2,3])
plt.show()

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


# Bacterial gDNA
# Calculating data before loop
mots = ['TNTC', 'TNTM', 'TYT', 'TTTH']
output = pd.DataFrame()
for i, _ in enumerate(mots):
	print(i)
	df = pd.DataFrame(data=fc.motifs_in_fasta('downloaded_DNA/all-bac-renamed.fa', str(mots[i])))
	output = output.append(df)
output.to_csv('output.csv')


# Bacterial gDNA
mots = ['TNTC', 'TNTM', 'TYT', 'TTTH']
input_data = pd.read_csv('output.csv')
# Colours for plotting. From Solarized
colour_palette = [(211, 54, 130), (108, 113, 196), (38, 139, 210), (42, 161, 152)]
plot_colours = colour_palette
for i in range(len(colour_palette)):
	r, g, b = colour_palette[i]
	# Convert RGB (0, 255) to (0, 1) which matplotlib likes
	plot_colours[i] = (r / 255, g / 255, b / 255)


fig, axes = plt.subplots(nrows=2, ncols=2, sharey=True, sharex=True)
for i, ax in enumerate(axes.flatten()):
	subset = input_data[input_data['motif seq'] == str(mots[i])]
	ax.bar(subset['Record'], subset['Perc DNA modified (total)'], color=plot_colours[i])
	ax.set_xticklabels(labels=subset['Record'], rotation=45, ha="right", rotation_mode="anchor")
	ax.set_title('\'{0}\' motif'.format(mots[i]))
	# Remove plot borders
	ax.spines["top"].set_visible(False)
	#ax.spines["bottom"].set_visible(False)
	ax.spines["right"].set_visible(False)
	#ax.spines["left"].set_visible(False)
	# Tick marks
	ax.get_yaxis().tick_left()
	ax.tick_params(axis="both", which="both", bottom=False, left=False)
	# Display yticks with intervals of 1. Alternative oneliner: ax.set_yticks(ax.get_yticks()[::2])
	ax.set_ylim([0, 3])
	ymin, ymax = ax.get_ylim()
	ax.set_yticklabels(np.arange(ymin, ymax+1, 1, dtype=np.int))
	#ax.set_ylabel('% gDNA modified')
# fig.suptitle('This is the figure title')
fig.text(0.01, 0.5, '% gDNA modified', va='center', rotation='vertical') # Common Y axis label
#plt.setp(axes[:, 0], ylabel='% gDNA modified') # Apply yaxis label only to plots in column 1
fig.tight_layout(rect=[0, 0.03, 1, 0.9]) # Call tight_layout last


matplotlib.use("TkAgg") # Backend to use for PyCharm interactive plots
matplotlib.use("GTK3Agg") # Backend to use for savefig with properly scaled DPI
fig.savefig("plots/colours - common y label.png", dpi=300)






# ONE PLOT ONLY
input_data = pd.read_csv('output.csv')
subset = input_data[input_data['motif seq'] == 'TNTC']
#subset['Record']
#subset['Perc DNA modified (total)']
plt.bar(subset['Record'], subset['Perc DNA modified (total)'])



# Plot % gDNA modified vs gDNA length




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