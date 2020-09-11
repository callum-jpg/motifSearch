#!/usr/bin/env python3
from motifsearch.classmotifs import motifsearch
import numpy as np

#%%
# Counting motifs with classes

mots = ['TNTC', 'TNTM', 'TYT', 'TTTH']


motif = motifsearch()

t = motif.anyone(mots)


#%%

y = motif.count('downloaded_DNA/bacterial-gDNA-renamed.fa', mots)

print(y)

# %%
from motifsearch.classmotifs import motifsearch

mots = ['GATC', 'AATT', 'TAATC', 'CATTT', 'TAATC', 'CATTT', 'TAATC', 'CATTT']

motif = motifsearch()

y = motif.count('downloaded_DNA/bacterial-gDNA-renamed.fa', mots)

x = motif.motif_bar1(mots)


#%%

# Save the figure
x.savefig("hello", dpi=300)


#%%

test = mots

if type(test) is str:
    test1 = [test]
    print('str')
elif type(test) is list:
    test1 = [x for x in test]
    print('list')
else:
    print('neither')
    
print(test)

#%%

test = ['tntnc']

print(len(test))
    
    
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






