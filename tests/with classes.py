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

y = motif.count('downloaded_DNA/bacterial-gDNA-renamed.fa', mots)

motif.motif_bar1(mots)


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
import matplotlib.pyplot as plt

fig, axes = plt.subplots(nrows=1, ncols=1, sharey=False, sharex=True)
y = np.array([axes])
for i, j in enumerate(y):
    print(i, j)
    
    
#%%

fig, ax = plt.subplots(2, 2, sharey=False, sharex=True)
axes = np.array(ax)
for i, axes1, in enumerate(axes.flatten()):
    print(i, axes1)


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






