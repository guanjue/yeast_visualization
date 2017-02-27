library(gplots)
library(RColorBrewer)
#library(pheatmap)

args <- commandArgs(trailingOnly = TRUE)
data=read.table(args[1],header=T)

#dim(data)
data=as.matrix(as.matrix(data)[,3:dim(data)[2]])
data=apply(data, 2, as.numeric)
print(dim(data))

#class(data)
threh=as.numeric(args[5])
print(quantile(data,threh)[1])
my_palette <- colorRampPalette(c("white", args[4]))(n = 19)
if (quantile(data,threh)[1] >0){
		colors = seq(0,quantile(data,threh)[1],length=20)
	} else{
		colors = seq(0,1,length=20)
	}

my_palette1 = colorRampPalette(c('darkslategray1','tan','sienna1','thistle1','darkorchid1','lightcyan','plum2','palegreen1','gray'))

#colors = c(seq(0,quantile(data,0.95)[1],length=25),seq(quantile(data,0.95)[1]+0.1,(quantile(data,0.95)[1]+0.1)*2,length=25))
#colors = c(seq(0,8,length=25),seq(8+0.1,8*2,length=25))

height_used=dim(data)[1]
width_used=dim(data)[2]

if (height_used < 100){
	height_used=100
}
if (width_used < 100){
	width_used=100
}

png(paste(args[3],'png',sep='.'),res=300,width=500,height=800)
heatmap.2(data, Rowv=NULL, Colv=NULL, col=my_palette,breaks=colors, margins=c(0.01,0.01),symm=F,symkey=F,symbreaks=F, key=F, labRow=F, labCol=F, keysize=0.001, scale="none", density.info="none", trace="none", dendrogram='none')
#pheatmap(data, cluster_rows=F, cluster_cols=F, color=my_palette,breaks=colors,legend=F,show_rownames=F,show_colnames=F,labels_col=)
dev.off()

