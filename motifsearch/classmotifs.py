#!/usr/bin/env python3
from motifsearch import countmotifs as ms
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

class motifsearch:
    
    def __init__(self):
        pass
        
    def __repr__(self):
        return "For finding motifs in a DNA sequence"
    
    def anyone(self, dna_motifs):
        return dna_motifs
    
    def count(self, input_file, motifs):
        output = pd.DataFrame()
        for i, _ in enumerate(motifs):
            print("Counting sites for {0} motif".format(motifs[i]))
            df = pd.DataFrame(data=ms.motifs_in_fasta(input_file, str(motifs[i])))
            output = output.append(df)
        #self.motif_data = output
        return output
        
    def motif_bar(self, input_data, motif, width=None, height=None):
        
        # Convert a str motif into a list element so the for loop works
        if type(motif) is str:
            motif = [motif]
        
        # Round up to calculate number of rows based on motifs
        self.rows = math.ceil(len(motif)/2)
        
        if len(motif) == 1:
            self.cols = 1
        else:
            self.cols = 2
            
        if width != None:
            figure_width = width
        elif self.cols == 1:
            figure_width = len(set(input_data['Record'])) * 0.5
        else:
            figure_width = (len(set(input_data['Record'])))
            
        if height != None:
            figure_height = height 
        # Adjust figure height based on row count
        elif self.rows == 1:
            figure_height = self.rows*0.5
        else:
            figure_height = self.rows*2
        
        # Colours for plotting. From Solarized palette
        colour_palette = [(42, 161, 152), (38, 139, 210), (108, 113, 196), (211, 54, 130)]
        plot_colours = colour_palette
        for i in range(len(colour_palette)):
        	r, g, b = colour_palette[i]
        	# Convert RGB (0, 255) to (0, 1) which matplotlib likes
        	plot_colours[i] = (r / 255, g / 255, b / 255)
        
        fig, ax_ = plt.subplots(self.rows, self.cols, sharey=False, sharex=True,
                                figsize=(figure_width, figure_height))
        # Make ax_ a 1/2D array, regardless of length
        axes = np.array(ax_)
        for i, ax in enumerate(axes.flatten()):
            # Plot bar only for subplots which has a corresponding motif
            # Delete unrequired subplot in the else statement
            if i in range(0, len(motif)):
                self.subset = input_data[input_data['motif seq'] == str(motif[i])]
                self.species = self.subset['Record']
                self.perc = self.subset['Perc DNA modified (total)']
                ax.bar(self.species, self.perc, color=plot_colours[i])
                # Customise plot
                ax.set_xticklabels(labels=self.species, rotation=45, ha="right", rotation_mode="anchor")
                ax.set_title('\'{0}\' motif'.format(motif[i]))
                # Remove plot borders    
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                ax.get_yaxis().tick_left()
                ax.tick_params(axis="both", which="both", bottom=False, left=False)
                # Display yticks with intervals of 1. Alternative oneliner: ax.set_yticks(ax.get_yticks()[::2])
                # Round up to the nearest int for the highest percentage
                ax.set_ylim([0, math.ceil(max(self.perc))])
                ymin, ymax = ax.get_ylim()
                ax.yaxis.set_ticks(np.arange(ymin, ymax+1, 1))  
            else:
                # I don't like this solution for removing a sub plot
                # But since the nrow/col is fixed at 2 you can assume that, assuming an odd number
                # of motifs is plotted, that the plot to delete is in the i = -1 and j = 1 of the 
                # np.array.
                # Assuming an even grid (2x2, 4x4 etc), this functionality will be ok. 
                axes[self.rows-1, self.cols-1].remove()
        # figure_width*1e-4 hopefully scales OK with varying dataset sizes
        fig.text(figure_width*1e-4, 0.55, '% gDNA modified', va='center', rotation='vertical') # Common Y axis label
        fig.tight_layout(rect=[0, 0.03, 1, 0.9]) # Call tight_layout last
        return fig
        
    
# for i, j in enumerate(mots):
#     if i in range(0, len(mots) - 1):
#         print("Plotted {}".format(i))
#     else:
#         print("I'll delete subplot {}".format(i))
        
    
# Useful code for getting Mb from bp ax.set_yticklabels([x/1000 for x in ax.get_yticks()])
    

	# Display yticks with intervals of 1. Alternative oneliner: ax.set_yticks(ax.get_yticks()[::2])
# 	ax.set_ylim([0, 3])
# 	ymin, ymax = ax.get_ylim()
# 	ax.set_yticklabels(np.arange(ymin, ymax+1, 1, dtype=np.int))