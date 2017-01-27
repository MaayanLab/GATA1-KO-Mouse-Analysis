
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

library = 'ChEA'
term = 'STAT5_ChIP-Seq_MAMMARY-EPITHELIUM_Mouse'

# measure = 'p-value'
measure = 'CS'

signature = {}
frames = []

# def plot_term(measure, df):

for dirn in ['up', 'dn']:
	for tissue in ['marrow', 'spleen']:

		data = pd.read_table('Enrichment/'+'_'.join(['GATA1',tissue,library,dirn+'.txt']),
			index_col=1)

		signature[tissue] = data

		sig = data.loc[:,measure]
		sig.name = tissue+'_'+dirn
		print sig.shape
		# print data[measure].name(tissue+'_'+dirn)

		frames.append(sig)

print len(frames)

df = pd.concat(frames,ignore_index=True,axis=1)
CD
print df.shape



