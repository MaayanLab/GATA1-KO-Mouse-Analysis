
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

pcter = 'STAT5 Genes/'

genelists = [	'STAT5_23275557_ChIP-Seq_MAMMARY-EPITHELIUM_Mouse.txt',
				# 'STAT5A_GM12878_hg19.txt',
				# 'STAT5A_K562_hg19.txt',
				'KEGG_JAK-STAT_genes.txt'
]

lists = []
for genelist in genelists:
	with open(pcter+genelist) as f: genes = f.read().splitlines()
	lists.append(genes)

consensus_genes = set.intersection(*map(set,lists))

print consensus_genes

# genelist = 'consensus_STAT5.txt'
# genelist = 'STAT5_23275557_ChIP-Seq_MAMMARY-EPITHELIUM_Mouse.txt'
# genelist = 'STAT5A_GM12878_hg19.txt'
# genelist = 'STAT5A_K562_hg19.txt' 
# genelist = 'KEGG_JAK-STAT_genes.txt'

method = 'fc'
# method = 'cd'

with open(pcter+'KEGG_JAK-STAT_genes.txt') as f: genes = f.read().splitlines()

print len(genes),"genes in the list"

spleen_data_fc = pd.read_table('GATA1_spleen_pctchange.txt', index_col=0)
marrow_data_fc = pd.read_table('GATA1_marrow_pctchange.txt', index_col=0)

spleen_data_cd = pd.read_table('GATA1_spleen_cdVector.txt',header=None, index_col=0, names=['SC'])
marrow_data_cd = pd.read_table('GATA1_marrow_cdVector.txt',header=None, index_col=0, names=['BMC'])

print spleen_data_fc.shape, marrow_data_fc.shape
print spleen_data_cd.shape, marrow_data_cd.shape

fc_mat = pd.concat([spleen_data_fc['logFC'],marrow_data_fc['logFC']],axis=1).dropna()
fc_mat.columns = ['SC','BMC']
fc_mat.to_csv('GATA1_fc.txt',sep='\t')

cd_mat = pd.concat([spleen_data_cd,marrow_data_cd],axis=1).dropna()
cd_mat.to_csv('GATA1_cd.txt',sep='\t')

print fc_mat.shape
print cd_mat.shape

if method == 'cd':
	df = cd_mat
	lim = 0.0015
	axis = [-0.02, 0.03, -0.02, 0.03]
else:
	df = fc_mat
	lim = 1.0
	axis = [-3.5, 3.5, -3.5, 3.5]
# genes = ['OSM', 'PTPN6', 'PIK3CG', 'MCL1']

overlap = [ x for x in genes if x in df.index.values ]
print len(overlap),"genes overlap with the expression matrix"

df_stat5 = df.loc[overlap,:]

colors = []
sizes = []
degs = []
markers = []
for gene in df_stat5.index.values:
	sc = df_stat5.loc[gene,'SC']
	bmc = df_stat5.loc[gene,'BMC']

	if abs(sc) < lim and abs(bmc) < lim:
		colors.append('lightgray')
		sizes.append(40)
		markers.append('.')
	else:
		if gene in consensus_genes:
			degs.append(gene)
			sizes.append(200)
			colors.append('red')
		else:
			sizes.append(40)
			colors.append('blue')

print df_stat5.shape

fig, ax = plt.subplots()
# plt.axis([-0.01, 0.01, -0.01, 0.01])
plt.axis(axis)
plt.plot([-1,-1],[-3.5,3.5],'w:',zorder=1)
plt.plot([1,1],[-3.5,3.5],'w:',zorder=1)
plt.plot([-3.5,3.5],[-1,-1],'w:',zorder=1)
plt.plot([-3.5,3.5],[1,1],'w:',zorder=1)
l1 = plt.plot([-3.5,3.5],[0,0],'w-',lw=1,zorder=1)
l2 = plt.plot([0,0],[-3.5,3.5],'w-',lw=1,zorder=1)
# plt.setp(l1,color='gray')
# plt.setp(l2,color='gray')

plt.scatter(df_stat5['SC'],df_stat5['BMC'],c=colors,s=sizes,zorder=2) #,marker=markers)
plt.xlabel('SC expression',fontsize=16)
plt.ylabel('BMC expression',fontsize=16)

for i, txt in enumerate(degs):
	# ax.annotate(txt, (df_stat5.loc[txt,'SC'],df_stat5.loc[txt,'BMC']))
	ax.annotate(txt, (df_stat5.loc[txt,:]+lim/10.0))

ax.set_axis_bgcolor('#CDCECB')

plt.show()

