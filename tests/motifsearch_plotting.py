#!/usr/bin/env python3
from motifsearch import motifSearch

#%%
# Counting motifs for the bacterial gDNA downloaded with downloadeding_dna.py

# Define motifs
mots = ['TNTC', 'GTCT', 'GGATC', 'CCTA']

# Instantiate motifsearch
motif = motifSearch()

# Count motifs in FASTA file
motif_df = motif.count('bacterial-gDNA-renamed.fa', mots)

# Save motif df to file
motif_df.to_csv('bac_counts.csv')

# Read saved df from file
# motif_df = pd.read_csv('bac_counts.csv')

# Plot motif modification burder bar plot for all motifs
motif.motif_bar(motif_df, mots)

# Plot genome length for species in Mb
motif.motif_bar_lengths(motif_df)

# Linear regression to compare relationship between modification burden and genome length
motif.count_vs_length(motif_df, mots)


#%% Saving example plots

motif.motif_bar(motif_df, mots).savefig('img/motifsearch-bar-plot', bbox_inches = 'tight', dpi=300)

motif.motif_bar_lengths(motif_df).savefig('img/motifsearch-bar-plot-lengths', bbox_inches = 'tight', dpi=300)

motif.count_vs_length(motif_df, mots).savefig('img/motifsearch-lin-reg', bbox_inches = 'tight', dpi=300)





#%%
# Counting motifs within human autosomes individually

from motifsearch import motifSearch
from natsort import index_natsorted
import numpy as np
import pandas as pd

mots = ['TNTC']

motif = motifSearch()

# Count motifs in autosomes. Slow to run
#human_counts = motif.count('downloaded_DNA/human-gdna/1-22-hs-gdna-renamed.fa', mots)

# Save human motif frequency df
# human_counts.to_csv('human_counts.csv')

# Read human motif frequency df
human_counts = pd.read_csv('human_counts.csv')

## Order 'chr1:22' into their natural order
# Key argument in sort_values allows for alternative sorting method to be applied
# index_natsorted first naturally order the 'Record' column
# np.argsort then returns an array for the sorted index_natsorted list, low to high, as a new array
# sort_values sorts 'Record' based on their returned, ordered indices, low to high 
# https://stackoverflow.com/questions/29580978/naturally-sorting-pandas-dataframe
# Understanding key:
# https://stackoverflow.com/questions/8966538/syntax-behind-sortedkey-lambda
human_counts_sorted = human_counts.sort_values(by='Record', key=lambda x: np.argsort(index_natsorted(human_counts['Record'])))

# Plot motif frequency for each chromosome
motif.motif_bar(human_counts_sorted, mots)

#%% Determine SEM of autosome modification
import numpy as np

human_std = np.std(human_counts['Perc DNA modified (total)'])

print(human_std)

#%% Comparing bacterial motif frequency with human


human_counts = pd.read_csv('human_counts.csv')

bac_counts = pd.read_csv('bac_counts.csv')

# Consolidate autosome data and find the average modification across chromosomes
hs_data = {
           "Record": ['H.sapiens'],
           "Total complement sites": [sum(human_counts['Total complement sites'])],
           "Total reverse complement sites": [sum(human_counts['Total reverse complement sites'])],
           "Total number of sites": [sum(human_counts['Total number of sites'])],
           "Perc DNA modified (total)": [np.mean(human_counts['Perc DNA modified (total)'])],
           "motif seq": ['TNTC'],
           "record length": [sum(human_counts['record length'])],
           "record length (Mb)": [sum(human_counts['record length (Mb)'])]
        }

# Save dict as dataframe
pd.DataFrame.from_dict(data=hs_data, orient='columns').to_csv('human_counts_summary.csv')

# Extract data for only E. coli and T. acquaticus for TNTC motif
taq_ecol_gdna = bac_counts[((bac_counts['Record'] == 'T.aquaticus') 
                            | (bac_counts['Record'] == 'E.coli(K12)'))
                    & (bac_counts['motif seq'] == 'TNTC')]

# Read consolidated human motif data
hs_gdna = pd.read_csv('human_counts_summary.csv')

# Concatenate bacterial with human data and save
taq_ecol_gdna.append(hs_gdna).to_csv('taq_ecol_human_gdna_counts.csv')

# Read bacterial-human motif frequency data
taq_ecol_human = pd.read_csv('taq_ecol_human_gdna_counts.csv')

# Plot for TNTC motif
motif = motifSearch()
motif.motif_bar(taq_ecol_human, 'TNTC')


#%% Adjust dimensions of the plot and save

plot = motif.motif_bar(taq_ecol_human, 'TNTC', width=2, height=3)

# Save the figure
plot.savefig("taq_ecol_human_counts.png", dpi=300)

#%% plotting gDNA length

plot = motif.motif_bar_lengths(taq_ecol_human, 2, 3)

plot.savefig("taq_ecol_human_lengths.png", dpi=300)