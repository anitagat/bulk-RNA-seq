# bulk-RNA-seq

# LOAD LIBRARIES
library(ggplot2)
library(reshape2)
library(ggrepel)
library(tidyverse)
library(amap)
library(devEMF)
library(clusterProfiler) 
library(org.Hs.eg.db)
library(enrichplot)
library(GOSemSim)
library(DOSE)

# LOAD DATA
# expression matrix = em, e.g: 
# em = read.table("EM.csv", header=TRUE, row.names=1, sep="\t")

# differential expression statistics for each pair of samples generated with DeSeq2 = de
# de = read.table("DE_Senes_vs_Prolif.csv", header=TRUE, row.names=1, sep="\t")

# gene annotations downoalded from Ensembl: annotations 
# annotations = read.table("Human_Background_GRCh38.p13.csv", header=TRUE, sep="\t")

# sample group information: ss
# ss = read.table("sample_sheet.csv", header = TRUE, sep="\t")

# standardise column names 
names(annotations) = c("Gene.ID", "Gene.name", "Chromosome", "Start", "End", "Type")

# set Gene.ID as row names
row.names(annotations) = annotations[,"Gene.ID"]

# remove double column (Gene.id)
annotations = annotations[,-1]

# merge em and annotations 
em_annotated  = merge(em, annotations, by= 0)

# reset Gene.ID as row names 
row.names(em_annotated) = em_annotated[,"Row.names"]

# eliminate extra-column with gene names
em_annotated = em_annotated[,-1] 

# parse: creates a fully parsed dataset for further analysis
parse = function(em_annotated, de, p_threshold, fold_threshold){
  
  # merge annotated expression matrix to differential expression matrix
  new.data = merge(em_annotated, de, by = 0)
  
  # create column with information on significance (TRUE/FALSE)
  new.data$sig = as.factor(new.data$p.adj < p_threshold & abs(new.data$log2fold) > fold_threshold)
  
  # create new column with -log10(p-value)
  new.data$mlog10 = -log10(new.data$p.adj)
  
  # create a vector with ordered row names 
  new_order = order(new.data[,"p.adj"], decreasing=FALSE)
  
  # reorder the whole dataframe using the reordered index vector 
  new.data = new.data[new_order,]
  
  # assign new row name based on existing column
  row.names(new.data) = new.data[,"Gene.name"]
  
  # create new variable column and assigns default 
  new.data$direction = "No change"
  
  # assigns new value based on conditions (p-value and fold change thresholds)
  new.data[new.data$p.adj < p_threshold & new.data$log2fold > fold_threshold,]$direction = "Up-regulated"
  new.data[new.data$p.adj < p_threshold & new.data$log2fold < -fold_threshold,]$direction = "Down-regulated"
  
  # cast variable column into factor type
  new.data$direction = as.factor(new.data$direction)
  
  # reorder levels by custom order
  new.data$direction <- factor(new.data$direction, levels = c("No change", "Up-regulated", "Down-regulated"))
  
  # return dataset 
  return(new.data) 
}

# scale.df: performs all operations required for scaling a df 
scale.df = function(df.to.be.scaled){
  na.omit(data.frame(t(scale(t(df.to.be.scaled)))))
}

# make a vector of significant genes (to be used later for ordering)
sig_genes = subset(emde_parsed, p.adj < 0.01 & abs(log2fold) > 2.0)
sig_genes = sig_genes[,"Gene.name"]

# make table with scaled expression values for the  significant genes only
emde_scaled_sig = emde_scaled[sig_genes,]  #utilises prev defined gene name vector for selection
emde_scaled_sig = na.omit(emde_scaled_sig)  

# PLOTTING with ggplot2
# create custom theme (mainly for homogenising font sizes across plots)
my_theme = theme(
  plot.title = element_text(size=16, hjust = 0.5, face = "bold"),
  legend.title= element_text(face="bold", size = 16),
  legend.text = element_text(face="bold",size = 14),
  axis.title.x = element_text(face="bold",  size = 16),
  axis.title.y= element_text(face="bold",  size = 16), 
  axis.text = element_text(face = "bold", size = 12),
  panel.background = element_rect(fill = NA, color = "black")
)

# save.png: saves plots in PNG format 
save.png = function(plot, name = "plot.path.name.emf", height, width) {
  png(file = name, height = height, width = width)
  print(plot)
  dev.off()
}

# save.emf: saves plots in EMF format
save.emf = function(plot, name = "plot.path.name.emf", height, width) {
  emf(file = name, height=height, width=width)
  print(plot)
  dev.off()
}

# Volocano plot
plot_volcano = function(de_table,p_threshold,fold_threshold) {
  
  # make new table with sig genes
  de_table_sig=subset(de_table, sig == TRUE)
  
  # make new table with sig upreg genes
  de_table_sig_up=subset(de_table_sig, p < p_threshold & log2fold > 0)
  
  # make new table with sig downreg genes
  de_table_sig_down=subset(de_table_sig, p < p_threshold & log2fold < 0)
  
  # make tables with top 10 sig up & down reg genes (these genes' names will be shown on plot)
  de_table_sig_up_top10 = de_table_sig_up[1:10,]
  de_table_sig_down_top10 = de_table_sig_down[1:10,]
  
  # plot volcano
  ggp=ggplot(de_table,aes(x=log2fold, y=mlog10, colour=direction)) +
    geom_point(size=1) +
    labs(x="log2fold change", y="-log10 p-value", size=10) +
    geom_vline(xintercept = -fold_threshold, linetype="dashed", color="grey", linewidth=0.5) + 
    geom_vline(xintercept = fold_threshold, linetype="dashed", color="grey", linewidth=0.5) + 
    geom_hline(yintercept = -log10(0.01), linetype="dashed", color="grey", linewidth=0.5) +
    xlim(c(-20,20)) + 
    ylim(c(0,310)) +   
    geom_text_repel(data=de_table_sig_up_top10, aes(label=Gene.name), show.legend=FALSE) +
    geom_text_repel(data=de_table_sig_down_top10, aes(label=Gene.name), show.legend=FALSE) + 
    scale_color_manual(values = c("black", "red", "blue"), labels=c("No Change", "Up-reg", "Down-reg")) +
    theme_classic() +
    my_theme +
    theme(legend.position = 'bottom', 
          legend.title = element_blank(), 
          legend.text = element_text(size=11))
  
  # print plot to environment
  print(ggp)
}

# MA plot
ma_plot = function(de_table,fold_threshold, samples) {
  # calculate row mean expression
  de_table$MeanExpression = rowMeans(de_table[,as.vector(samples)])
  
  # make new table with significant genes
  de_table_sig=subset(de_table, sig == TRUE)
  
  # make new table with sig upregulated genes
  de_table_sig_up=subset(de_table_sig, log2fold > 0)
  
  # make new table with sign downregulated genes
  de_table_sig_down=subset(de_table_sig, log2fold < 0)
  
  # make tables with top 5 sig up & down reg genes 
  de_table_sig_up_top5 = de_table_sig_up[1:5,]
  de_table_sig_down_top5 = de_table_sig_down[1:5,]
  
  # plot MA plot
  ggp = ggplot(de_table, aes(x=log10(MeanExpression), y=log2fold, color=sig)) +
    geom_point(size=2)+
    labs(xlab= "log10(Mean Expression)", ylab="log2fold change", size=10) +
    geom_hline(yintercept = -log10(0.01), linetype="dashed", color="grey", size=0.5) +
    geom_hline(yintercept = log10(0.01), linetype="dashed", color="grey", size=0.5) +
    scale_color_manual(values=c("black","red"), labels=c("Not sig","Sig")) +
    my_theme +
    theme(legend.position='bottom', legend.title = element_blank()) 
  
  # print plot to env
  print(ggp)
}

# create vector with group names 
group.labels = unique(ss$SAMPLE_GROUP)
# note: change to scale_color_brewer() to adjust colors to sample numbers

# PCA function
pca = function(expression_matrix_scaled, samples, groups){
  # perform PCA
  pca = prcomp(t(as.matrix(sapply(expression_matrix_scaled, as.numeric))))  
  pca_coordinates = data.frame(pca$x)
  pca_coordinates$sample = samples
  
  # calculate %
  vars = apply(pca$x, 2, var)
  prop_x = round(vars["PC1"]/sum(vars),4)*100
  prop_y = round(vars["PC2"]/sum(vars),4)*100
  x_axis_label = paste("PC1","(", prop_x,"%)", sep = "")
  y_axis_label = paste("PC2","(", prop_y,"%)", sep = "")
  
  # plot PCA
  ggp = ggplot(data=pca_coordinates, aes(x=PC1, y=PC2), colour=groups) +
    geom_point(aes(colour=groups)) +
    geom_text_repel(aes(label=sample, colour=groups), show.legend=FALSE, size=5) +
    labs(x=x_axis_label, y=y_axis_label) +
    scale_color_manual(values = c("red", "blue"), labels=group.labels) +
    theme_classic() +
    my_theme +
    theme(legend.position = 'bottom', legend.title = element_blank())
  
  # print PCA
  print(ggp)
}

# Boxplot (faceted)
plot_boxplot = function(gene_vector, expression_matrix, samples, ncolumns){
  gene_data = data.frame(t(expression_matrix[gene_vector,]))
  gene_data$sample_group = as.factor(samples)
  gene_data.m = melt(gene_data, id.vars = "sample_group") 
  ggp = ggplot(gene_data.m, aes(x=sample_group, y=value, colour=sample_group)) + 
    # geom_boxplot() + 
    geom_jitter() +
    labs(y="Expression (Z-score)") +
    facet_wrap(~variable, ncol=ncolumns) +
    scale_color_manual(values = c("red","blue"), labels = group.labels) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_blank(),
      legend.text = element_text(face="bold",size=12),
      legend.title = element_blank(),
      legend.background = element_blank(),
      panel.grid.major = element_blank(),
      strip.text.x = element_text(size = 12, face="bold", vjust=1), 
      panel.spacing = unit(1, "lines"), 
      legend.position = "bottom",
      axis.ticks.x = element_blank(),
      strip.background = element_blank()
    )
  #print plot to env 
  print(ggp)
}

# Heatmap
plot_heatmap = function(em_scaled_sig, x.y.clust) { 
  # convert df into matrix to retain level information
  hm.matrix = as.matrix(em_scaled_sig) 
  
  # get distances between genes (method: spearman correlation)
  y.dist = Dist(hm.matrix, method="spearman") 
  
  # perform the clustering using the distances (build the dendrogram)
  y.cluster = hclust(y.dist, method="average")  
  
  # make dendrogram object
  y.dd = as.dendrogram(y.cluster)
  
  # disentangle the dendrogram (re-order)
  y.dd.reorder = reorder(y.dd,0,FUN="average")
  
  # get the new gene order (as a vector of row indexes) 
  y.order = order.dendrogram(y.dd.reorder)
  
  # cluster the X axis too: transpose to calculate distances by COL
  hm.matrix.t = t(hm.matrix)
  
  # calculate distance between values
  x.dist = Dist(hm.matrix.t, method = "spearman")
  
  # build dendrogram (method: averaging)
  x.cluster = hclust(x.dist, method = "average")
  
  # make dendrogram object
  x.dd = as.dendrogram(x.cluster)
  
  # reorder tree based on similarity nodes
  x.dd.reorder = reorder(x.dd,0,FUN="average")
  
  # get the new gene order (as a vector of row indexes) 
  x.order = order.dendrogram(x.dd.reorder)
  
  # reorder the heatmap matrix using new gene order, based on desired axes to reorder
  if (x.y.clust == "x") {
    hm.matrix_clustered = hm.matrix[,x.order]
  } 
  else if (x.y.clust == "y") {
    hm.matrix_clustered = hm.matrix[y.order,]
  }
  else if (x.y.clust == "x.y") {
    hm.matrix_clustered = hm.matrix[y.order,x.order]
  }
  else {
    print("Accepted arguments: x, y, or x.y")
  }
  
  # melt matrix
  hm.matrix_clustered= melt(hm.matrix_clustered)
  
  # set colours
  colours = c("blue","white","red")
  palette = colorRampPalette(colours)(100)
  
  # plot heatmap
  ggp = ggplot(data=hm.matrix_clustered, aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile() + 
    labs(title= "Proliferative VS Senescent", y="Genes",x="Samples") +
    scale_fill_gradientn(colours = colorRampPalette(colours)(100), name="Z-score") + 
    theme(axis.text.y = element_blank(), 
          axis.title.y = element_text(face="bold",size=14),
          axis.ticks = element_blank(),
          axis.text.x = element_text(face="plain",size=11),
          axis.title.x = element_blank(),
          plot.title = element_text(face = "bold",size=16, hjust = 0.5),
          legend.title = element_text(face="bold",size=12), 
          legend.spacing.x = unit(0.25, 'cm'))
  print(ggp)
}

# Pathway Analysis
pathway_analysis = function(sig_genes_list, fromType){
  
  # convert to ENTREZ
  if (fromType == "symbol"){
    sig_genes_entrez = bitr(sig_genes, fromType = "SYMBOL", toType="ENTREZID", OrgDb = org.Hs.eg.db)
  }
  else {
    sig_genes_entrez_ens = bitr(sig_genes_ensembl, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }
  # calculate gene ontology enrichment & save results in global env 
  pathway_results <<- enrichGO(gene=sig_genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, readable = T, ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
  
  # present results in multiple ways
  bar <<- barplot(pathway_results, showCategory = 5, font.size=16, label_format = 30) #showCategory specifies the number of gene-sets shown
  dot <<- dotplot(pathway_results, showCategory = 10)
  goplot <<- goplot(pathway_results, showCategory = 10)
  cnetplott <<- cnetplot(pathway_results, categorySize= "pvalue", circular = FALSE, node_label="gene", color.params = list(edge = TRUE))  
  
  # print plots
  print(bar)
  print(dot)
  print(goplot)
  print(cnetplott)
}

# plot emap to look at ontology clusters (edges represents overlapping genes)
catdist = pairwise_termsim(pathway_results)
emap = emapplot(catdist, category_node=1.5, category_label= 0.8, layout.params = list(layout = "kk"), showCategory=10, colorEdge=TRUE)

# get_pathway_genes: gest list of most enriched genes (n = nth most significant ontology) 
get_pathway_genes = function(pathway_results, n) {
  
  # extract gene ID, description, and p.adj from ORA analysis results
  gene_sets = pathway_results$geneID
  description = pathway_results$Description
  p.adj = pathway_results$p.adjust
  
  # create dataframe with extracted data from above
  ora_results = data.frame(cbind(gene_sets, description, p.adj))
  
  # select genes from ontologies (nth ontology ordered by enrichment p-value)
  enriched_gene_set = as.character(ora_results[n,1])
  
  # create vector of most enriched genes 
  enriched_genes = unlist(strsplit(enriched_gene_set, "/")) 
  return(enriched_genes)
}

# Boxplot (provide vector with most significantly enriched genes)
enriched_genes_boxplot = function(enriched_genes_list, em_scaled, sample_group, ngenes, title) {
  
  # get the n most enriched genes of most significant ontology
  enriched_genes_list = enriched_genes_list[1:ngenes]
  
  # select n most enriched genes from scaled expression matrix
  enriched_gene_data = em_scaled[enriched_genes_list,]
  enriched_gene_data = data.frame(t(enriched_gene_data))
  
  # add sample group col (factor)
  enriched_gene_data$sample_group = sample_group
  
  # melt dataframe
  enriched_gene_data.m = melt(enriched_gene_data, id.vars = "sample_group") #id.vars = "col name" determines what col to keep static
  enriched_gene_data.m$sample_group = as.factor(enriched_gene_data$sample_group)
  
  # make boxplot with n most enriched genes
  ggp = ggplot(enriched_gene_data.m, aes(x=variable, y=value)) + 
    geom_boxplot(aes(fill=sample_group)) + 
    #geom_jitter(aes(color=sample_group)) +
    labs(title = title, y="Expression (Z-score)") +
    #scale_fill_brewer(type="seq",palette = 1, name="Sample group") +
    scale_fill_manual(values=c("red","blue"), labels=group.labels) +
    theme_classic(base_size=12) +
    my_theme +
    theme(
      axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
      axis.text.y = element_text(face = "bold"),
      axis.title.y.left = element_text(face = "bold"),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size=12),
      plot.title = element_text(size=15, hjust = 0, face = "bold"),
      panel.grid = element_blank(),
      plot.tag = element_blank(),
      plot.subtitle = element_blank(),
      strip.text.x = element_text(size = 16, family="Arial", face="bold", vjust=1), panel.spacing = unit(1, "lines"), 
      legend.position = "right"
    )
  
  # print plot to env
  print(ggp)
}

