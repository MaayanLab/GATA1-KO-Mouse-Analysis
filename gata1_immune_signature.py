
import pandas as pd
import numpy as np
from scipy.spatial.distance import jaccard

# library = 'ChEA'
library = 'MGI'

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
		# print data[measure].name(tissue+'_'+dirn)

		frames.append(sig)

df = pd.concat(frames,axis=1)

if measure == 'p-value':
	df = -np.log10(df)
	df = df.fillna(0.0)
	df = df[np.sum(df,axis=1) > 10]
	df = df[~(df<=5).all(axis=1)]
	df[df < 1.3] = 0.0

else:
	df = df.fillna(0.0)
	
	df = df[np.sum(df,axis=1) > 20]
	# df = df[~(df==0.0).all(axis=1)]
	# df = df[~(df<=2).all(axis=1)]
	df[df < 1] = 0.0

print df.shape

df.to_csv('GATA1_MGI_signature_'+measure+'.txt',sep='\t')

terms = df.index.values
codes = [ x.split(' ')[0].upper() for x in terms ]

gmt = pd.read_table('MGI_Mammalian_Phenotype_Level_4.txt',
	index_col=0, dtype='string', header=None)

grab_rows = [ x for x in gmt.index.values if x.split('_')[0] in codes ]

gmt = gmt.loc[grab_rows,:]

distmat = np.zeros((len(gmt),len(gmt)))

for i in range(len(gmt)):
	for j in range(len(gmt)):

		s1 = set(gmt.iloc[i,:].dropna().values)
		s2 = set(gmt.iloc[j,:].dropna().values)
		distmat[i,j] = len(s1&s2)/float(len(s1|s2))

pd.DataFrame(1-distmat, index=terms, columns=terms).to_csv('GATA1_MGI_distmat.txt',sep='\t')


