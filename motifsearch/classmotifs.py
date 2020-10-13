#!/usr/bin/env python3
from motifsearch import countmotifs as ms
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
from sklearn.linear_model import LinearRegression


class motifsearch:
    
    def __init__(self):
        # Colours for plotting. From Solarized palette
        colour_palette = [(42, 161, 152), (38, 139, 210), (108, 113, 196), (211, 54, 130)]
        self.plot_colours = colour_palette
        for i in range(len(colour_palette)):
        	r, g, b = colour_palette[i]
        	# Convert RGB (0, 255) to (0, 1) which matplotlib likes
        	self.plot_colours[i] = (r / 255, g / 255, b / 255)
        
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
        rows = math.ceil(len(motif)/2)
        
        if len(motif) == 1:
            cols = 1
        else:
            cols = 2
            
        if width != None:
            figure_width = width
        elif cols == 1:
            figure_width = len(set(input_data['Record'])) * 0.5
        else:
            figure_width = (len(set(input_data['Record'])))
            
        if height != None:
            figure_height = height 
        # Adjust figure height based on row count
        elif rows == 1:
            figure_height = rows*0.5
        else:
            figure_height = rows*2

        
        fig, ax_ = plt.subplots(rows, cols, sharey=False, sharex=True,
                                figsize=(figure_width, figure_height))
        # Make ax_ a 1/2D array, regardless of length
        axes = np.array(ax_)
        for i, ax in enumerate(axes.flatten()):
            # Plot bar only for subplots which has a corresponding motif
            # Delete unrequired subplot in the else statement
            if i in range(0, len(motif)):
                subset = input_data[input_data['motif seq'] == str(motif[i])]
                species = subset['Record']
                perc = subset['Perc DNA modified (total)']
                ax.bar(species, perc, color=self.plot_colours[i])
                # Customise plot
                ax.set_xticklabels(labels=species, rotation=45, ha="right", rotation_mode="anchor", fontstyle='italic')
                ax.set_title('\'{0}\' motif'.format(motif[i]))
                # Remove plot borders    
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                ax.get_yaxis().tick_left()
                ax.tick_params(axis="both", which="both", bottom=False, left=False)
                # Display yticks with intervals of 1. Alternative oneliner: ax.set_yticks(ax.get_yticks()[::2])
                # Round up to the nearest int for the highest percentage
                # ax.set_ylim([0, math.ceil(max(perc))])
                # ymin, ymax = ax.get_ylim()
                # ax.yaxis.set_ticks(np.arange(ymin, ymax+1, 1))  
                if max(perc) < 1:
                    # For modifications below 1, find rescale the y axis
                    upper_bound = np.round(max(perc)*1.5, 1)
                    ax.set_ylim([0, upper_bound])
                    ymin, ymax = ax.get_ylim()
                    np.arange(ymin, ymax, upper_bound)
                else:
                    ax.set_ylim([0, math.ceil(max(perc))])
                    ymin, ymax = ax.get_ylim()
                    ax.yaxis.set_ticks(np.arange(ymin, ymax+1, 1))  
            else:
                # I don't like this solution for removing a sub plot
                # But since the nrow/col is fixed at 2 you can assume that, assuming an odd number
                # of motifs is plotted, that the plot to delete is in the i = -1 and j = 1 of the 
                # np.array.
                # Assuming an even grid (2x2, 4x4 etc), this functionality will be ok. 
                axes[rows-1, cols-1].remove()
        # figure_width*1e-4 hopefully scales OK with varying dataset sizes
        fig.text(figure_width*1e-4, 0.55, '% gDNA modified', va='center', rotation='vertical') # Common Y axis label
        fig.tight_layout(rect=[0, 0.03, 1, 0.9]) # Call tight_layout last
        return fig
    
    def motif_bar_lengths(self, input_data, width=None, height=None):

        species = input_data['Record']
        dna_lengths = input_data['record length (Mb)']        

        if width != None:
            figure_width = width
        else:
            figure_width = len(set(input_data['Record']))*0.8
            
        if height != None:
            figure_height = height 
        else:
            figure_height = math.ceil(max(dna_lengths))*0.7
    
        # If there is a >50 Mb gap in recorded gDNA lengths, plot y as log
        if max(dna_lengths) - min(dna_lengths) > 50:
            log_y = True
        else:
            log_y = False
                
    
        fig, ax = plt.subplots(1, 1, sharey=False, sharex=True,
                                figsize=(figure_width, figure_height))
        ax.bar(species, dna_lengths, color=self.plot_colours[2])
        # Customise plot
        ax.set_xticklabels(labels=species, rotation=45, ha="right", rotation_mode="anchor", fontstyle='italic')
        # ax.set_title('DNA sequence lengths (Mb)')
        # Remove plot borders    
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.get_yaxis().tick_left()
        ax.tick_params(axis="both", which="both", bottom=False, left=False)
        # Display yticks with intervals of 1. Alternative oneliner: ax.set_yticks(ax.get_yticks()[::2])
        # Round up to the nearest int for the highest percentage
        if log_y == True:
            ax.set_yscale('log', basey=10)
            ax.set_ylim([1, math.ceil(max(dna_lengths))])            
        else:
            ax.set_ylim([0, math.ceil(max(dna_lengths))])
            ymin, ymax = ax.get_ylim()
            ax.yaxis.set_ticks(np.arange(ymin, ymax+1, 1))  
                
        # figure_width*1e-4 hopefully scales OK with varying dataset sizes
        fig.text(figure_width*1e-4, 0.55, 'DNA length (Mb)', va='center', rotation='vertical') # Common Y axis label
        fig.tight_layout(rect=[0, 0.03, 1, 0.9]) # Call tight_layout last
        return fig
        
    
    def count_vs_length(self, input_data, motif, width=None, height=None):
        linear_reg = LinearRegression()
        
        if width != None:
            figure_width = width
        else:
            figure_width = (len(set(input_data['Record'])))
            
        if height != None:
            figure_height = height 
        else:
            figure_height = math.ceil(max(input_data['record length (Mb)']  ))*0.8
            
            
        # Round up to calculate number of rows based on motifs
        rows = math.ceil(len(motif)/2)
        
        if len(motif) == 1:
            cols = 1
        else:
            cols = 2
        
        fig, ax_ = plt.subplots(rows, cols, sharey=False, sharex=False,
                                figsize=(figure_width, figure_height))
        # Make ax_ a 1/2D array, regardless of length
        axes = np.array(ax_)
        for i, ax in enumerate(axes.flatten()):
            if i in range(0, len(motif)): # If i isn't in range, else will remove the subplot
                subset = input_data[input_data['motif seq'] == str(motif[i])]
                
                # x, y for linear_reg
                x = np.reshape(subset['record length (Mb)'].values, (len(subset['record length (Mb)'].values), 1))
                y = np.reshape(subset['Perc DNA modified (total)'].values, (len(subset['Perc DNA modified (total)'].values), 1))
                	# Linear regression predict
                linear_reg.fit(x, y)
                y_pred = linear_reg.predict(x) # Prediction holder
                	# Plot
                ax.scatter(x, y, color=self.plot_colours[i])
                ax.plot(x, y_pred, color='red')

                # Customise plot
                ax.set_title('\'{0}\' motif'.format(motif[i]))
                # Remove plot borders    
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                ax.get_yaxis().tick_left()
                ax.tick_params(axis="both", which="both", bottom=False, left=False)
                # Display yticks with intervals of 1. Alternative oneliner: ax.set_yticks(ax.get_yticks()[::2])
                # Round up to the nearest int for the highest percentage
                if max(y) < 1:
                    # For modifications below 1, find rescale the y axis
                    upper_bound = np.round(max(y)*1.5, 1)
                    ax.set_ylim([0, upper_bound])
                    ymin, ymax = ax.get_ylim()
                    np.arange(ymin, ymax, upper_bound)
                else:
                    ax.set_ylim([0, math.ceil(max(y))])
                    ymin, ymax = ax.get_ylim()
                    ax.yaxis.set_ticks(np.arange(ymin, ymax+1, 1))  
            else:
                # I don't like this solution for removing a sub plot
                # But since the nrow/col is fixed at 2 you can assume that, assuming an odd number
                # of motifs is plotted, that the plot to delete is in the i = -1 and j = 1 of the 
                # np.array.
                # Assuming an even grid (2x2, 4x4 etc), this functionality will be ok.        
                axes[rows-1, cols-1].remove()
        # figure_width*1e-4 hopefully scales OK with varying dataset sizes
        fig.text(figure_width*1e-4, 0.5, '% gDNA modified', va='center', rotation='vertical') # Common Y axis label
        fig.text(0.5, 0.04, 'DNA length (Mb)', ha='center', va='center')
        fig.tight_layout(rect=[0, 0.03, 1, 0.9]) # Call tight_layout last
        return fig
        
        