
setwd("~/Documents/GATA1-Analysis")

library(gplots)
library(RColorBrewer)
library(heatmap3)

# lib = 'ChEA'
# direction = 'down'
# drug = 'SOR'
# drug_cutoff = 85

# data = read.csv('GATA1_MGI_signature_CS3.txt',sep='\t',header=TRUE,row.names=1)
data = read.csv('GATA1_MGI_signature_p-value3.txt',sep='\t',header=TRUE,row.names=1)

distmat = read.csv('GATA1_MGI_distmat.txt',sep='\t',header=TRUE,row.names=1)

dendro = as.dendrogram(hclust(as.dist(distmat)))

# colors = colorRampPalette(brewer.pal(20,'Reds'))
colors = colorRampPalette(brewer.pal(9,'YlOrRd'))
# colors = colorRampPalette(rev(brewer.pal(100,'RdYlBu')))
# colors = colorRampPalette(rev(brewer.pal(6,'PRGn')))(20:100)

display_data <- abs(as.matrix(data))



heatmap.2(display_data,
          trace='none',
          density.info='density', denscol='black', densadj = 1,
          margins =c(6,28),
          # Rowv = dendro,
          Colv = NULL,
          col=colors,
          colsep=c(1,2,3,4),
          offsetRow = 0.1,
          offsetCol = 0.1,
          cexCol=1,
          dendrogram='row',
          # ColSideColors = drug_labels,
          key.title='', key.xlab='Combined Score', key.ylab='',
          labCol = c('SC','BMC','SC','BMC'),
          adjCol = c(1,0.5),
          labRow = data$rownames)

