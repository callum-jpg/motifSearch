#!/usr/bin/env python3
from motifsearch.classmotifs import motifsearch

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

mots = ['GATC', 'GTCT', 'GGATC']

motif = motifsearch()

y = motif.count('downloaded_DNA/bacterial-gDNA-renamed.fa', mots)

x = motif.motif_bar(mots)


#%%

# Save the figure
x.savefig("hello", dpi=300)

#%%

print([x for x in range(0, 3)])

#%%
import matplotlib.pyplot as plt

fig, ax_ = plt.subplots(2, 2, sharey=False, sharex=True)
# Make ax_ a 1/2D array, regardless of length
axes = np.array(ax_)

axes[0,1].remove()

# Smarter removing of excess subplots?
# position = np.array([[0, 1 ], [2, 3]])
# #print(np.where(position==3))
# np.where(position==3).remove()

#%%
import numpy as np

position = np.array([[0, 1 ], [2, 3]])

#print(np.where(position==3))

x = np.where(position==3)

print([i [0] for i in x])

#%%

mots = ['GATC', 'TNTC', 'CHC', 'TTTC']
ncol = 2
ncol = 2

for i, j in enumerate(mots):
    if i in range(0, len(mots) - 1):
        print("Plotted {}".format(i))
    else:
        print("I'll delete subplot {}".format(i))
        


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






