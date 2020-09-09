#!/usr/bin/env python3

#%% Import packages

from motifsearch import countmotifs as ms
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

#%% Create RGB colours for plotting

# Colours for plotting. From Solarized palette
colour_palette = [(211, 54, 130), (108, 113, 196), (38, 139, 210), (42, 161, 152)]
plot_colours = colour_palette
for i in range(len(colour_palette)):
	r, g, b = colour_palette[i]
	# Convert RGB (0, 255) to (0, 1) which matplotlib likes
	plot_colours[i] = (r / 255, g / 255, b / 255)

#%%


# Bacterial gDNA
# Calculating data before plot loop
# Assumes DNA has already been downloaded 
mots = ['TNTC', 'TNTM', 'TYT', 'TTTH']
output = pd.DataFrame()
for i, _ in enumerate(mots):
	print("Counting sites for {0} motif".format(mots[i]))
	df = pd.DataFrame(data=ms.motifs_in_fasta('tests/downloaded_DNA/all-bac-renamed.fa', str(mots[i])))
	output = output.append(df)
output.to_csv('output.csv')


#%% Plot gDNA lengths

input_data = pd.read_csv('output.csv')
# Subset to take TNTC only (motif doesn't matter, just don't want to repeatedly plot the same gDNA length)
subset = input_data[input_data['motif seq'] == 'TNTC']

fig, ax = plt.subplots()
ax.bar(subset['Record'], subset['record length'], color=plot_colours[1])
ax.set_title('gDNA length between species')
ax.set_xticklabels(labels=subset['Record'], rotation=45, ha="right", rotation_mode="anchor")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_yaxis().tick_left()
ax.tick_params(axis="both", which="both", bottom=False, left=False)
# Display yticks with intervals of 1. Alternative oneliner: ax.set_yticks(ax.get_yticks()[::2])
ax.set_ylim([0, 7e6])
ymin, ymax = ax.get_ylim()
ax.set_yticklabels(np.arange(ymin, ymax + 1e6, 1e6, dtype=np.int) / 1e6) #  Convert to Mb
ax.set_ylabel('gDNA length (Mb)')
fig.tight_layout(rect=[0, 0.03, 1, 0.9]) # Call tight_layout last

fig.savefig("tests/gDNA length between species.png", dpi=300)


#%% Plot output.csv for different motifs occurance in bacterial gDNA

input_data = pd.read_csv('output.csv')

fig, axes = plt.subplots(nrows=2, ncols=2, sharey=False, sharex=True)
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


fig.savefig("tests/colours - common x+y - all numbers presented.png", dpi=300)

#%% Plot gDNA length vs % DNA modified

linear_reg = LinearRegression() 
fig, axes = plt.subplots(nrows=2, ncols=2, sharey=False, sharex=True)
for i, ax in enumerate(axes.flatten()):
	subset = input_data[input_data['motif seq'] == str(mots[i])]
	# Reshape data for LinearRegression()
	x = np.reshape(subset['record length'].values, (len(subset['record length'].values), 1))
	y = np.reshape(subset['Perc DNA modified (total)'].values, (len(subset['Perc DNA modified (total)'].values), 1))
	# Linear regression predict
	linear_reg.fit(x, y)
	y_predict = linear_reg.predict(x) # Prediction holder
	# Plot
	ax.scatter(x, y)
	ax.plot(x, y_predict, color='red')  # Plot prediction
	# Format plot
	ax.set_title('\'{0}\' motif'.format(mots[i]))
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.get_yaxis().tick_left()
	ax.tick_params(axis="both", which="both", bottom=False, left=False)
	# Display yticks with intervals of 1. Alternative oneliner: ax.set_yticks(ax.get_yticks()[::2])
	ax.set_ylim([0, 3])
	ymin, ymax = ax.get_ylim()
	ax.set_yticklabels(np.arange(ymin, ymax+1, 1, dtype=np.int))
	# # x axis format
	ax.set_xlim([2e6, 6.5e6])
	xmin, xmax = ax.get_xlim()
	ax.set_xticklabels(np.arange(xmin, xmax + 1e6, 1e6, dtype=np.int) / 1e6)
fig.text(0.5, 0.04, 'gDNA length (Mb)', ha='center') # Common x axis label
fig.text(0.01, 0.5, '% gDNA modified', va='center', rotation='vertical')  # Common Y axis label
fig.tight_layout(rect=[0.02, 0.05, 1, 0.9]) # Call tight_layout last

fig.savefig("tests/gDNA length vs % gDNA modified.png", dpi=300)


#%% Human gDNA plotting




mots = ['TNTC', 'TNTM', 'TYT', 'TTTH']
output = pd.DataFrame()
for i, _ in enumerate(mots):
	print("Counting sites for {0} motif".format(mots[i]))
	df = pd.DataFrame(data=ms.motifs_in_fasta('downloaded_DNA/human-gdna/1-22-hs-gdna-renamed.fa', str(mots[i])))
	output = output.append(df)
output.to_csv('output_human_gdna.csv')

#%% Plot percentage of gDNA modified for all motifs

mots = ['TNTC', 'TNTM', 'TYT', 'TTTH']
input_data = pd.read_csv('output_human_gdna.csv')

fig, axes = plt.subplots(nrows=2, ncols=2, sharey=False, sharex=True, figsize=(12, 5))
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
	ax.set_ylim([0, 5])
	ymin, ymax = ax.get_ylim()
	ax.set_yticklabels(np.arange(ymin, ymax+1, 1, dtype=np.int))
	#ax.set_ylabel('% gDNA modified')
# fig.suptitle('This is the figure title')
fig.text(0.0001, 0.5, '% gDNA modified', va='center', rotation='vertical') # Common Y axis label
#plt.setp(axes[:, 0], ylabel='% gDNA modified') # Apply yaxis label only to plots in column 1
fig.tight_layout(rect=[0, 0.03, 1, 0.9]) # Call tight_layout last. (left, bottom, right, top)


fig.savefig("tests/human gDNA - all motifs", dpi=300)



#%% Plotting TNTC only for human gDNA

input_data = pd.read_csv('output_human_gdna.csv')

fig, axes = plt.subplots(figsize=(10, 2.5))
subset = input_data[input_data['motif seq'] == 'TNTC']
axes.bar(subset['Record'], subset['Perc DNA modified (total)'], color=plot_colours[i])
axes.set_xticklabels(labels=subset['Record'], rotation=45, ha="right", rotation_mode="anchor")
axes.set_title('\'{0}\' motif'.format('TNTC'))
# Remove plot borders
axes.spines["top"].set_visible(False)
axes.spines["right"].set_visible(False)
# Tick marks
axes.tick_params(axis="both", which="both", bottom=False, left=False)
# Display yticks with intervals of 1. Alternative oneliner: ax.set_yticks(ax.get_yticks()[::2])
axes.set_ylim([0, 3])
ymin, ymax = axes.get_ylim()
axes.set_yticks(np.arange(ymin, ymax+1, 1, dtype=np.int)) # Set the values of the ticks
axes.set_yticklabels(np.arange(ymin, ymax+1, 1, dtype=np.int)) # Set the displayed values of the ticks
#ax.set_ylabel('% gDNA modified')
# fig.suptitle('This is the figure title')
fig.text(0.0001, 0.5, '% gDNA modified', va='center', rotation='vertical') # Common Y axis label
#plt.setp(axes[:, 0], ylabel='% gDNA modified') # Apply yaxis label only to plots in column 1
fig.tight_layout(rect=[0, 0.03, 1, 0.9]) # Call tight_layout last. (left, bottom, right, top)


fig.savefig("tests/human gDNA - TNTC only", dpi=300)


#%% Modification of chromosomes vs length


np.arange(min(input_data['record length']), max(input_data['record length']), 10e6, dtype=np.int) / 1e6

#linear_reg = LinearRegression()
fig, axes = plt.subplots()
subset = input_data[input_data['motif seq'] == 'TNTC']
# Reshape data for LinearRegression()
x = np.reshape(subset['record length'].values, (len(subset['record length'].values), 1))
y = np.reshape(subset['Perc DNA modified (total)'].values, (len(subset['Perc DNA modified (total)'].values), 1))
# Linear regression predict
#linear_reg.fit(x, y)
#y_predict = linear_reg.predict(x) # Prediction holder
# Plot
axes.scatter(x, y)
#ax.plot(x, y_predict, color='red')  # Plot prediction
# Format plot
axes.set_title('\'{0}\' motif'.format('TNTC'))
axes.spines["top"].set_visible(False)
axes.spines["right"].set_visible(False)
axes.get_yaxis().tick_left()
axes.tick_params(axis="both", which="both", bottom=False, left=False)
# Display yticks with intervals of 1. Alternative oneliner: ax.set_yticks(ax.get_yticks()[::2])
axes.set_ylim([0, 3])
ymin, ymax = axes.get_ylim()
axes.set_yticks(np.arange(ymin, ymax+1, 1, dtype=np.int))
axes.set_yticklabels(np.arange(ymin, ymax+1, 1, dtype=np.int))
# # x axis format
axes.set_xlim([40e6, 260e6])
xmin, xmax = axes.get_xlim()
axes.set_xticks(np.arange(xmin, xmax, 10e6, dtype=np.int))
axes.set_xticklabels(np.arange(xmin, xmax, 10e6, dtype=np.int) / 1e6, rotation=45, ha="right", rotation_mode="anchor")
fig.text(0.5, 0.04, 'gDNA length (Mb)', ha='center') # Common x axis label
fig.text(0.01, 0.5, '% gDNA modified', va='center', rotation='vertical')  # Common Y axis label
fig.tight_layout(rect=[0.02, 0.05, 1, 0.9]) # Call tight_layout last



#%% Comparing bacterial gDNA with human gDNA - % modified

human_gdna = pd.read_csv("output_human_gdna.csv")
bac_gdna = pd.read_csv("output.csv")

human_gdna_tntc = human_gdna[(human_gdna['motif seq'] == 'TNTC')]
np.mean(human_gdna_tntc['Perc DNA modified (total)'])


# Subset the desired data for bacterial gDNA
taq_gdna = bac_gdna[(bac_gdna['Record'] == 'T.aquaticus') & (bac_gdna['motif seq'] == 'TNTC')]
# Conditionally subset two record types corresponding to one motif seq
taq_ecol_gdna = bac_gdna[((bac_gdna['Record'] == 'T.aquaticus') | (bac_gdna['Record'] == 'E.coli(K12)'))
                    & (bac_gdna['motif seq'] == 'TNTC')]

# Creating new dataset for just Taq and human gDNA
# Kinda hacky using values[0], but it works for now
gdna_dict = {'Record': [taq_ecol_gdna['Record'].values[0], taq_ecol_gdna['Record'].values[1], 'H.sapiens'],
             'Perc DNA modified (total)': [taq_ecol_gdna['Perc DNA modified (total)'].values[0], taq_ecol_gdna['Perc DNA modified (total)'].values[1], np.mean(human_gdna_tntc['Perc DNA modified (total)'])]
}

# Plotting
fig, axes = plt.subplots(figsize=(2, 3))
axes.bar(gdna_dict['Record'], gdna_dict['Perc DNA modified (total)'], color=plot_colours[i])
axes.set_xticklabels(labels=gdna_dict['Record'], rotation=45, ha="right", rotation_mode="anchor")
axes.set_title('\'{0}\' motif'.format('TNTC'))
# Remove plot borders
axes.spines["top"].set_visible(False)
axes.spines["right"].set_visible(False)
# Tick marks
axes.tick_params(axis="both", which="both", bottom=False, left=False)
# Display yticks with intervals of 1. Alternative oneliner: ax.set_yticks(ax.get_yticks()[::2])
axes.set_ylim([0, 2])
ymin, ymax = axes.get_ylim()
axes.set_yticks(np.arange(ymin, ymax+1, 1, dtype=np.int)) # Set the values of the ticks
axes.set_yticklabels(np.arange(ymin, ymax+1, 1, dtype=np.int)) # Set the displayed values of the ticks
#ax.set_ylabel('% gDNA modified')
# fig.suptitle('This is the figure title')
fig.text(0.0001, 0.5, '% gDNA modified', va='center', rotation='vertical') # Common Y axis label
#plt.setp(axes[:, 0], ylabel='% gDNA modified') # Apply yaxis label only to plots in column 1
fig.tight_layout(rect=[0, 0.03, 1, 0.9]) # Call tight_layout last. (left, bottom, right, top)


fig.savefig("tests/Taq, Ecol vs human gDNA - TNTC modification", dpi=300)































