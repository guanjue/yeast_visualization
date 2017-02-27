library(gplots)
library(RColorBrewer)
#library(pheatmap)

args <- commandArgs(trailingOnly = TRUE)
data=read.table(args[1],header=T)
data1=read.table(args[2],header=T)

#dim(data)
data=as.matrix(as.matrix(data)[,3:dim(data)[2]])
data=apply(data, 2, as.numeric)
print(dim(data))

gene_class_table=table(data1)
gene_class_num_bound=as.numeric(dimnames(gene_class_table)$data1)

print(gene_class_num_bound)

#class(data)
threh=as.numeric(args[5])
print(quantile(data,threh)[1])
my_palette <- colorRampPalette(c("white", args[4]))(n = 9)
if (quantile(data,threh)[1] >0){
		colors = seq(0,quantile(data,threh)[1],length=10)
	} else{
		colors = seq(0,1,length=10)
	}

my_palette1_colors = c('darkslategray1','tan','sienna1','thistle1','darkorchid1','lightcyan','plum2','palegreen1','gray')[gene_class_num_bound]
my_palette1 = colorRampPalette(my_palette1_colors)
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
jBrewColors <- brewer.pal(n = 8, name = "Pastel1")
color_bar=jBrewColors
color_bar=rainbow(8)
color_bar=my_palette1
data1=as.matrix(cbind(data1,data1))
png(paste(args[3],'color_split.png',sep='.'),res=300,height=800)
heatmap.2(data1, Rowv=NULL, Colv=NULL, col=color_bar, margins=c(0.01,0.01),symm=F,symkey=F,symbreaks=F, key=F, labRow=F, labCol=F, keysize=0.001, scale="none", density.info="none", trace="none", dendrogram='none')
#pheatmap(data, cluster_rows=F, cluster_cols=F, color=my_palette,breaks=colors,legend=F,show_rownames=F,show_colnames=F,labels_col=)
dev.off()
