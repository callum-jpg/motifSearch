#!/usr/bin/env python3
from motifsearch.classmotifs import motifsearch

#%%
# Counting motifs with classes
from motifsearch.classmotifs import motifsearch

mots = ['TNTC', 'GTCT', 'GGATC']
#mots = ['GATC']

motif = motifsearch()

y = motif.count('downloaded_DNA/bacterial-gDNA-renamed.fa', mots)

y.to_csv('bac_counts.csv')

x = motif.motif_bar(y, mots)


#%%

# Save the figure
x.savefig("hello", dpi=300)


#%%
# Counting motifs within human autosomes individually

from motifsearch.classmotifs import motifsearch
from natsort import index_natsorted
import numpy as np
import pandas as pd

mots = ['TNTC']

motif = motifsearch()

#human_counts = motif.count('downloaded_DNA/human-gdna/1-22-hs-gdna-renamed.fa', mots)

human_counts.to_csv('human_counts.csv')

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

human_plot = motif.motif_bar(human_counts_sorted, mots)

#%%

human_counts = pd.read_csv('human_counts.csv')

#%% 


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


taq_ecol_gdna = bac_counts[((bac_counts['Record'] == 'T.aquaticus') 
                            | (bac_counts['Record'] == 'E.coli(K12)'))
                    & (bac_counts['motif seq'] == 'TNTC')]


hs_gdna = pd.read_csv('human_counts_summary.csv')

taq_ecol_gdna.append(hs_gdna).to_csv('taq_ecol_human_gdna_counts.csv')

taq_ecol_human = pd.read_csv('taq_ecol_human_gdna_counts.csv')


motif = motifsearch()
motif.motif_bar(taq_ecol_human, 'TNTC')

#%%

plot = motif.motif_bar(taq_ecol_human, 'TNTC', width=2, height=3)

# Save the figure
plot.savefig("taq_ecol_human_counts.png", dpi=300)

#%%

unique = ['T.acq', 'H.sap', 'M.tub', 'H.sap', 'M.tub']

# Return unique values from list
print(set(unique))

# Return number of rows, columns
print(y.shape)

print(len(set(y['Record'])))

# print(human_counts['Record'])
# print(sorted(human_counts['Record']))



    
#%%

print(y[['Record', 'Perc DNA modified (total)', 'motif seq']])


#%% Theoritical use of classes

# Count motifs
# Have count method return pd.df
motif_df = motifs.count(mots)

# Plot the motif df
motifs.plot_perc(motif_df)

# Plot gDNA lengths
motifs.plot_length(motif_df)

# Plot length vs perc modified
motifs.plot_length_mod(motif_df, regression = T)






