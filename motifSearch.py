#!/usr/bin/env python
import functions as fc
import pandas as pd

#fc.percentage(50, 100)
#fc.iupac_to_dna('TNTM')

output = pd.DataFrame(data=fc.motifs_in_fasta('GRCh38_latest_genomic.fna', 'TNTC'))
output.to_csv('output.csv')
