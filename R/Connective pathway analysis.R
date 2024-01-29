library(reticulate)
library(Matrix)
library(dplyr)
library(ggplot2)
library(tidyr)
library(nichenetr)
library(tidyverse)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(doParallel)
library(clusterProfiler)  
library(DOSE)
library(topGO) 
library(ReactomePA)


input.dir = paste0(wd,"/Your_pathway_folder")
outdir= paste0(wd,"/Outout_pathway_folder")
if (dir.exists(outdir)==F){dir.create(outdir,recursive = TRUE)
  print("outdir created")}

################################################################
# connective pathway analysis_ function####
run_cpa = function(geneset_list = c("IPA"),
                   output_dir = outdir,
                   n_cluster = 40,
                   width = 80,
                   height =6){
  if (dir.exists(output_dir)==F){dir.create(output_dir,recursive = TRUE)
    print("output_dir created")}
  ##  Merge pathways from different genesets####
  head(total_pathways)
  IPA = total_pathways[,c('Description','p.adjust' ,'geneID','celltype', 'geneset','z_score','Pathway_Type')]
  IPA$geneID = gsub(',','/',IPA$geneID)

  All_pathways_comb = data.frame()
  for (i in geneset_list){
    All_pathways_comb = rbind(All_pathways_comb,get(i))
  }
  table(All_pathways_comb$geneset)

  table(All_pathways_comb$geneset,All_pathways_comb$celltype)
  write.table(All_pathways_comb,file = paste0(output_dir,"/All_pathways_comb.txt"),row.names =F,quote = F,sep = "\t")

  ## Get unique genes of each pathways in combine of all celltypes ####
  pw_gene_list <- list()
  n = 1
  for (i in unique(All_pathways_comb$Description)) {
    a <- All_pathways_comb[All_pathways_comb$Description == i,]
    b <- strsplit(a$geneID,"/") %>% unlist() %>% unique()
    b = gsub("[ (includes others)]","",b)
    pw_gene_list[n] <- list(b)
    names(pw_gene_list)[n] <- i
    n = n+1
  }
  
  ## calculate Jaccard index and Jaccard Distance####
  JD_matrix <- matrix(ncol=3)
  colnames(JD_matrix) <- c("from","to","Jaccard_Distance")
  write.table(JD_matrix,file = paste0(output_dir,"/JD_matrix_full.txt"),sep = '\t',row.names =F)
  for (i in 1:length(pw_gene_list)){
    for (j in 1:length(pw_gene_list)){
      a = intersect(pw_gene_list[i] %>% unlist ,pw_gene_list[j] %>% unlist) %>% length #shared gene
      b = c(pw_gene_list[i] %>% unlist,pw_gene_list[j] %>% unlist) %>% unique  %>% length #all unique genes
      JI = a/b #jaccard index
      JD = 1-JI
      c = matrix(c(names(pw_gene_list)[i],names(pw_gene_list)[j],JD),ncol = 3)
      colnames(c) <- c("from","to","Jaccard_Distance")
      write.table(c,file = paste0(output_dir,"/JD_matrix_full.txt"),row.names =F,col.names = F,append=T,sep = '\t',quote=F)
      # JD_matrix = rbind(JD_matrix,c)
    }
    print(paste0("from pathway ",i," done!"))
  }

  JD_matrix = read.csv(paste0(output_dir,"/JD_matrix_full.txt"),sep = '\t')
  JD_matrix = JD_matrix %>% as.data.frame()
  JD_matrix = JD_matrix[2:dim(JD_matrix)[1],]
  head(JD_matrix)
  dim(JD_matrix)
  max(JD_matrix$Jaccard_Distance)
  
  ## MAKE Plot for Jaccard distance####
  library("reshape")
  xx = JD_matrix
  # xx = aggregate(Jaccard_Distance ~ ., data=xx, FUN=min)
  xx = reshape(xx,idvar="from", timevar="to", direction="wide")
  xx[1:5,1:5]
  rownames(xx) = xx$from
  xx = xx[,2:dim(xx)[2]]
  colnames(xx) = sapply(strsplit(colnames(xx),"Jaccard_Distance."),'[[',2)
  write.table(xx,file = paste0(output_dir,"/JD_matrix_full_wide.txt"),sep = '\t',row.names =T)
  # xx = read.table(paste0(output_dir,"/JD_matrix_full_wide.txt"),sep = '\t',header = T)
  
  ##get the cluster info####
  library(dendextend)
  hc.x = hclust(as.dist(xx),method='ward.D2') # xx contains Jaccard Index - symmetric matrix same ordering of pathways in columns and in rows
  clustersGlobal = cutree(hc.x,k=n_cluster,order_clusters_as_data=F)
  clustersGlobal = as.data.frame(clustersGlobal)
  clustersGlobal$key = rownames(clustersGlobal)
  hc.x$order
  rownames_after_cluster = rownames(xx[hc.x$order,])
  write.table(clustersGlobal,file = paste0(output_dir,"/Cluster_info_p",n_cluster,".txt"),sep = '\t',row.names =T)
  
  hc <- hc.x %>%
    color_branches(k = 40) %>%
    color_labels(k = 40) %>% 
    set("labels_cex", c(.5)) 
  pdf(file=paste0(output_dir,"/Dendro_circl_useJDclusters_",n_cluster,".pdf"), width = 15, height = 15)
  circlize_dendrogram(as.dendrogram(hc,split=40),
                            dend_track_height = 0.2)
  dev.off()
  
  ##get gene list of each program####
  gene_list_prog = data.frame()
  for (i in clustersGlobal$clustersGlobal){
    i_pathways = clustersGlobal[clustersGlobal$clustersGlobal ==i,"key"]
    i_genes = c()
    for (j in i_pathways){
      i_genes= c(i_genes,unlist(pw_gene_list[[j]])) %>% unique
    }
    i_genes = data.frame(genes = i_genes, cluster = i)
    gene_list_prog = rbind(gene_list_prog,i_genes)
  }
  gene_list_prog = unique(gene_list_prog)
  write.table(gene_list_prog,file = paste0(output_dir,"/Cluster_gene_info_p",n_cluster,".txt"),sep = '\t',row.names =F)
  
  ##  get the dendrogram and order based on Jaccard distance clustering####
  summary(dgp) #see the detail of dgp
  cluster <- dgp$tree_row #obtain the cluster information
  # plot(cluster,
  #      hang = -1, #Put the labels at the same height
  #      cex=0.3,axes=FALSE,ann=FALSE)
  dgp_roworder = rownames(xx[dgp$tree_row[["order"]],])
  identical(rownames_after_cluster,dgp_roworder)
  
  # Heatmaps############################################################################################################
  ## Heatmap using z-score and clustering from JD ####
  print("Start drawing Heatmap")
  #change format from long to wide
  All_pathways_comb <- read.csv(paste0(output_dir,"/All_pathways_comb.txt"),sep="\t",header = T)
  # colnames(All_pathways_comb) <- All_pathways_comb[1,]
  # All_pathways_comb <- All_pathways_comb[2:dim(All_pathways_comb)[1],]
  head(All_pathways_comb)
  
  df = All_pathways_comb[,c('Description','z_score','p.adjust','celltype','geneset','Pathway_Type')]
  df$ptw_geneset = paste0(df$Description,"_",df$geneset)
  head(df)
  dim(df)
  df = df[!duplicated(df),]
  dim(df)
  df = df[!duplicated(paste0(df$Description,df$celltype)),]
  
  library(reshape2)
  df$z_score = as.numeric(df$z_score)
  df.wide = dcast(df,Description ~ celltype, value.var = "z_score") %>% as.data.frame()
  rowname = df.wide$Description
  df.wide = df.wide[,2:dim(df.wide)[2]]
  df.wide[1:5,1:4]
  df.wide[is.na(df.wide)==1] <- -999
  # df.wide = apply(df.wide,2,as.numeric)  %>% as.data.frame()
  rownames(df.wide) = rowname
  # df.wide = as.matrix(df.wide)
  df.wide[1:5,1:4]
  df.wide = df.wide[colnames(xx),] #make the same order as xx
  df.wide[1:5,1:4]
  max(df.wide,na.rm=T)
  min(df.wide,na.rm=T)
  
  # df.wide = df.wide[which(rowSums(df.wide) != 0),]
  
  # add annotation
  Cell_type = colnames(df.wide)
  df_anno = df[,c("Description","Pathway_Type","geneset")] %>% unique
  df_anno = df_anno[!duplicated(df_anno$Description),]
  rownames(df_anno) = df_anno$Description
  pathway_type = df_anno[rownames(df.wide),'Pathway_Type'] %>% as.vector()
  geneset = df_anno[rownames(df.wide),'geneset'] %>% as.vector()
  
  anno_col = data.frame(pathway_type=pathway_type,
                        geneset=geneset)
  anno_col$pathway_type = factor(anno_col$pathway_type,levels=unique(anno_col$pathway_type))
  anno_col$geneset = factor(anno_col$geneset,levels=unique(anno_col$geneset))
  
  library(RColorBrewer)
  library(colorspace)
  
  ha <- HeatmapAnnotation(#pathway_type = anno_col$pathway_type,
                          geneset=geneset,
                          col = list(#pathway_type = c("Signaling" = "#74ae48", "Metabolic" = "#D7191C"), # rowAnnotation or HeatmapAnnotation
                                     geneset = c('IPA' = "#E16A86",'KEGG' = "#C18500",'GO' = "#00AB6E",'Reactome' = "#6C8EE6"))) #HeatmapAnnotation
  
  p1_zscore = Heatmap(t(df.wide),
                      # split
                      # row_split = 20,
                      column_split = n_cluster,  #split by dendrograms
                      name = "Z-score", #title of legend
                      column_title = "P_%s",
                      column_title_side = "top",   
                      column_title_gp = gpar(fontsize = 6, fontface = "bold"),  
                      column_names_gp = gpar(fontsize = 4),
                      show_column_names = T,
                      row_title = "Enriched Pathway",
                      row_title_side ="right",   
                      row_title_gp = gpar(fontsize = 8) ,
                      row_names_gp = gpar(fontsize = 7),  
                      # heatmap_legend_param = list(title="legend"),
                      top_annotation = ha, #put annotation
                      cluster_columns = color_branches(cluster, k = n_cluster),  
                      show_column_dend = T,
                      cluster_rows =F ,
                      col = colorRamp2(c(-999,-6,0,6,999), c("black","blue","white","red","grey")),
                      border = T
  ) 
  pdf(file=paste0(output_dir,"/Heatmap_zscore_useJDclusters_",n_cluster,".pdf"), width = width, height = height)
  print(p1_zscore)
  dev.off()

###########make a circlized heatmap###########
  library(circlize)
  library(dendextend)
  
  pdf(file=paste0(output_dir,"/Heatmap_zscore_useJDclusters_",n_cluster,"_circle.pdf"), width = 10, height = 10)
  circos.par("start.degree" = 90,cell.padding = c(0, 0, 0, 0), gap.after = c(rep(0.5,39), 10)) 
  col_fun1 = colorRamp2(c(-999,-6, 0, 6, 999), c("black","blue", "white", "red","grey"))
  split = clustersGlobal[colnames(xx),'clustersGlobal'] %>% factor
  circos.heatmap(df.wide, split = split,col = col_fun1, dend.side = "inside", rownames.side = "outside", 
                 rownames.cex = 0.28, track.height = 0.15,
                 cluster = F,show.sector.labels = F)
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    if(CELL_META$sector.numeric.index == 40) { # the last sector
      cn = colnames(df.wide)
      n = length(cn)
      circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                  n - 1:n + 0.5, cn, 
                  cex = 0.3, adj = c(0, 0.5), facing = "inside")
    }
  }, bg.border = NA)
  circos.clear()
  dev.off()
  
  pdf(file=paste0(output_dir,"/Heatmap_zscore_useJDclusters_",n_cluster,"_circle2.pdf"), width = 6, height = 6)
  circos.par("start.degree" = 90,cell.padding = c(0, 0, 0, 0), gap.after = c(rep(0.5,39), 10)) 
  circos.heatmap(df.wide, split = split,col = col_fun1, 
                 rownames.cex = 0.3, 
                 cluster = F,show.sector.labels = T)
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    if(CELL_META$sector.numeric.index == 40) { # the last sector
      cn = colnames(df.wide)
      n = length(cn)
      circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                  n - 1:n + 0.5, cn, 
                  cex = 0.3, adj = c(0, 0.5), facing = "inside")
    }
  }, bg.border = NA)
  circos.clear()
  dev.off()
  
  pdf(file=paste0(output_dir,"/Heatmap_zscore_useJDclusters_",n_cluster,"_circle3.pdf"), width = 15, height = 15)
  circos.par("start.degree" = 90,cell.padding = c(0, 0, 0, 0), gap.degree = 10) 
  hc <- hc.x %>%
    color_branches(k = 40) %>%
    color_labels(k = 40) %>% 
    set("labels_cex", c(.3)) 
  circlize_dendrogram(as.dendrogram(hc,split=40),
                      dend_track_height = 0.28,labels = F)
  circos.clear()
  dev.off()
  
}

run_cpa(geneset_list = c("IPA"),
        output_dir = paste0(outdir),
        n_cluster = 40,
        width = 40,
        height =6)





