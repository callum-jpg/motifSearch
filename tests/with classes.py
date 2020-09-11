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

mots = ['TNTC', 'TNTM', 'TYT', 'TTTH']

motif = motifsearch()

y = motif.count('downloaded_DNA/bacterial-gDNA-renamed.fa', ['TNTC'])

motif.motif_bar('TNTC')


#%%

test = '[1, 2, 3]'

if type(test) is str:
    print("hello")
else:
    print('waaa')
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






