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
output = pd.DataFrame(data=fc.motifs_in_fasta('downloaded_DNA/all-bac-renamed.fa', 'TNTC'))
output.to_csv('output.csv')

plt.bar(output['Record'].values, output['Perc DNA modified (total)'].values)
plt.suptitle(output['motif seq'][1], fontsize=25)


# Create 2x2 plot for 4 searched motifs
# Calculates motif data in each loop - calculate data before loop?
mots = ['TNTC', 'TNTM', 'TYT', 'TTTH']

fig, axes = plt.subplots(nrows=2, ncols=2, sharey=True)
fig.subplots_adjust(hspace=0.5)
#axes.setp(axes.xaxis.get_majorticklabels(), rotation=45)
#plt.xticks(rotation=45)
for i, ax in enumerate(axes.flat):
	output = pd.DataFrame(data=fc.motifs_in_fasta('downloaded_DNA/all-bac-renamed.fa', str(mots[i])))
	ax.bar(output['Record'].values, output['Perc DNA modified (total)'].values)
	ax.set(title=mots[i])
	print(i, ax)
plt.xticks(rotation=45) # only rotates labs on last plot
#axes.setp(plt.xticks()[0], rotation=45)
plt.show()






# Calculating data before loop
mots = ['TNTC', 'TNTM', 'TYT', 'TTTH']
output = pd.DataFrame()
for i, _ in enumerate(mots):
	print(i)
	df = pd.DataFrame(data=fc.motifs_in_fasta('downloaded_DNA/all-bac-renamed.fa', str(mots[i])))
	output = output.append(df)
output.to_csv('output.csv')



# Plot full dataset in one
# Doesn't work...
input_data = pd.read_csv('output.csv')


fig, axes = plt.subplots(nrows=2, ncols=2, sharey=True)

axes = axes.flatten()
for i, j in enumerate(input_data['motif seq'].unique()):
	#print(i, j)
	subset = input_data[input_data['motif seq'] == str(j)]
	t = subset[['Record', 'Perc DNA modified (total)']]
	axes.plot(subset['Record', 'Perc DNA modified (total)'], kind='bar')



# Another attempt, same as old.
input_data = pd.read_csv('output.csv')
fig, axes = plt.subplots(nrows=2, ncols=2)
fig.subplots_adjust(hspace=0.5)
keys = input_data['motif seq'].unique()
for ax, name in zip(axes.flatten(), keys):
	loop_data = input_data[input_data['motif seq'] == str(name)]
	#print(loop_data)
	ax.bar(loop_data['Record'].values, loop_data['Perc DNA modified (total)'].values)

# https://www.dataquest.io/blog/tutorial-advanced-for-loops-python-pandas/
input_data.values.T



# Trying again...
input_data = pd.read_csv('output.csv')

fig, axes = plt.subplots(nrows=2, ncols=2)
keys = input_data['motif seq'].unique()
for i, j in enumerate(keys):
	loop_data = input_data[input_data['motif seq'] == str(name)]
	loop_data_1 = input_data[['Record', 'Perc DNA modified (total)']]
	#print(loop_data_1)
	loop_data_1.plot(kind='bar', ax=axes[i])

t = input_data[['Record', 'motif seq']]

t = input_data[input_data['Perc DNA modified (total)'], input_data['Record'], input_data['motif seq' == 'TNTC']]

input_data[input_data['Record'] == 'P.aeruginosa']

pd.Series(t)
pd.Series(t, index=range(0, len(input_data)))




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