setwd("D:/work/CQMU_CYP/project/CQMU_CQ_skin/pre-post-N2/dedoublets")

library(Seurat)
library(dplyr)
#library(Nebulosa)
library(ggplot2)
library(paletteer)
library(scCATCH)
#library(export)

source('D:/self/my_functions_R/my_palettes.R')
source('D:/self/my_functions_R/scRNA_plots.R')

pbmc <- readRDS('s8_raw.rds')
pbmc[["percent.rb"]] <- PercentageFeatureSet(pbmc,pattern = "^RP[SL]",assay = 'RNA')
a1 <- as.data.frame(table(pbmc$orig.ident))
colnames(a1) <- c("samples","beforeQC")
p1 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) & xlab(NULL)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & nCount_RNA < 50000 & percent.mt < 20)
a2 <- as.data.frame(table(pbmc$orig.ident))
colnames(a2) <- c("samples","afterQC")
a <- merge(a1,a2)
write.csv(a,'cell_num.csv',row.names = F)
p1 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) & xlab(NULL)
ggsave("after_qc.pdf", plot = p1, width = 16, height = 4.5)
pbmc <- readRDS('s8_qc.rds')
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 3000)

pbmc <- ScaleData(pbmc, verbose = FALSE,features = rownames(pbmc))
pbmc <- RunPCA(pbmc, npcs = 100, verbose = FALSE)
pbmc <- RunUMAP(pbmc, 
                reduction = "pca", 
                metric = "euclidean",
                n.neighbors = 40,
                #min.dist=0.1,
                n.epochs=500,
                negative.sample.rate = 10L,
                dims = 1:30)
pbmc <- RunTSNE(pbmc, reduction = "pca", dims = 1:30)
saveRDS(pbmc,'s8_dimenition.rds')

pbmc <- readRDS('s8_dimenition.rds')
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.6,random.seed = 1234) 

alpha.use <- 0.4
p1 <- DimPlot(pbmc,raster=FALSE) 
p1 <- DimPlot(pbmc,group.by = "orig.ident",raster=FALSE) + ggtitle(NULL) 
p1 <- DimPlot(pbmc,split.by = "orig.ident",ncol = 3,raster = F) + ggtitle(NULL) 
p1 <- DimPlot(pbmc,reduction = 'tsne',raster=FALSE) 
p1 <- DimPlot(pbmc,reduction = 'tsne',group.by = "orig.ident",raster=FALSE) + ggtitle(NULL) 
p1$layers[[1]]$mapping$alpha <- alpha.use
p1 <- p1 + scale_alpha_continuous(range=alpha.use,guide="none")
p1
ggsave("UMAP_3.pdf", plot = p1, width = 6, height = 4.5)
ggsave("TSNE_2.pdf", plot = p1, width = 6, height = 4.5)

DimPlot(pbmc,label = T,label.size = 6) + NoLegend()
sce.markers <- FindAllMarkers(pbmc, 
                              only.pos = TRUE, 
                              #min.diff.pct = 0.2,
                              #test.use="bimod",
                              #slot="counts",
                              min.pct = 0.1, 
                              logfc.threshold = .1,
                              return.thresh=0.05)
write.csv(sce.markers,'findallmarkers.csv')
#sce.markers <- read.csv('findallmarkers.csv',row.names = 1)
sce.markers <- sce.markers[order(sce.markers$cluster,sce.markers$avg_log2FC,decreasing = T),]
sce.markers1 <- sce.markers %>% group_by(cluster) %>% slice_head(n=30)
clu_ann <- scCATCH(sce.markers,
                   species = 'Human',
                   cancer = NULL,
                   tissue = 'Skin')
tops <- sce.markers %>% group_by(cluster) %>% slice_head(n=3)
DoHeatmap(object = pbmc, features = tops$gene,angle = 45) + NoLegend()
obj.aver <- AverageExpression(pbmc)
obj.aver <- data.frame(obj.aver$RNA)
obj.aver <- obj.aver[which(rowSums(obj.aver) > 0.005),]
obj.aver <- obj.aver[!grepl("\\.", rownames(obj.aver)),]
obj.aver <- obj.aver[!grepl("RPS", rownames(obj.aver)),]
write.csv(obj.aver,'allgene_avg.csv',quote = F)

markers <- c(
  "CLEC14A","ECSCR","PECAM1",#Endothelial
  "CDH5",#vascular endothelial cell(vEC) "PECAM1",
  "LYVE1","PROX1",#lymphatic endothelial cells(LECs)
  "TAGLN","ACTA2",#Mural
  #"RGS5","CD36","FABP4",#pericytes
  #"MYH11","CNN1","RERGL",#vascular Smooth muscle(vSMC)
  "COL1A1","DCN","SFRP2",#Fibroblast
  "NRXN1","SCN7A",#Neural cells
  "CPA3","VWA5A","CTSG",#mast
  "CD3D","TRAC", #T
  "CD68","C1QA",#Macrophage
  #"LYZ",#Monocyte
  #"HLA-DRA","CD1C","CLEC10A",#DC
  "DCT","TYRP1","MLANA",#Melanocyte????ϸ??
  "KRT14","KRT1","KRT5", #Keratinocyte????ϸ??"KRT19","KRT10",
  #"CCL21","LYVE1",#Langerhan cells
  #"THY1","CD44","ITGB1","ENG",#adipocyte
  #"CD79A","MS4A1",#B
  #"MKI67","STMN1","TOP2A",#Cycling
  "AQP5","MUCL1"#sweat gland cells
  #"RUNX1","CD34",#Hair follicle cell
  #"S100B","NCAM1",#Merkel
)

p1 <- DimPlot(pbmc,label = T,label.size = 6,raster = F) + NoLegend()
DefaultAssay(pbmc)<-"RNA"
p2 <- DotPlot(pbmc, features = markers) + theme_bw() + coord_flip() +  
  theme(panel.grid = element_blank(),  
        axis.text.x=element_text(hjust = 0.5,vjust=1)) + #angle=90,
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c('#330066','#336699','#66CC66','#FFCC33')) + 
  labs(x=NULL,y=NULL) 
DefaultAssay(pbmc)<-"integrated"
p1+p2

cluster_anno <- function(pbmc){
  C1 = c(28) 
  C2 = c(19,21) 
  #C3 = c(18) 
  C4 = c(2,4,6,7,11,12,23)
  C5 = c(8,13,14,15) 
  C6 = c(20) 
  C7 = c(0,1,3,5,9,24,27,31) 
  #C8 = c(11) 
  #C9 = c(2,23)
  C10 = c(18,22,25) 
  C11 = c(16,26,29) 
  C12 = c(17)
  C13 = c(10,30)
  current.cluster.ids <- c(C1,C2,C4,C5,C6,C7,C10,C11,C12,C13)
  new.cluster.ids <- c(rep("Melanocyte",length(C1)),
                       rep("Immune",length(C2)),
                       #rep("Myeloid",length(C3)),
                       rep("Fibroblast",length(C4)),
                       rep("VEC",length(C5)),
                       rep("LEC",length(C6)),
                       rep("Mural",length(C7)),
                       #rep("vSMC1",length(C8)),
                       #rep("vSMC2",length(C9)),
                       rep("Sweat gland cell",length(C10)),
                       rep("Keratinocyte",length(C11)),
                       rep("Mast",length(C12)),
                       rep("Neural",length(C13))
  )
  pbmc@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(pbmc@meta.data$seurat_clusters)), 
                                            from = current.cluster.ids, 
                                            to = new.cluster.ids)
  return(pbmc)
}

cluster_anno <- function(pbmc){
  C1 = c(28) 
  C2 = c(19,21) 
  #C3 = c(18) 
  C4 = c(2,4,6,7,11,12,23)
  C5 = c(8,13,14,15) 
  C6 = c(20) 
  C7 = c(0,1,3,5,9,24,27,31) 
  #C8 = c(11) 
  #C9 = c(2,23)
  C10 = c(18,22,25) 
  C11 = c(16,26,29) 
  C12 = c(17)
  C13 = c(10,30)
  current.cluster.ids <- c(C1,C2,C4,C5,C6,C7,C10,C11,C12,C13)
  new.cluster.ids <- c(rep("C06",length(C1)),
                       rep("C02",length(C2)),
                       #rep("Myeloid",length(C3)),
                       rep("C01",length(C4)),
                       rep("C10",length(C5)),
                       rep("C04",length(C6)),
                       rep("C07",length(C7)),
                       #rep("vSMC1",length(C8)),
                       #rep("vSMC2",length(C9)),
                       rep("C09",length(C10)),
                       rep("C03",length(C11)),
                       rep("C05",length(C12)),
                       rep("C08",length(C13))
  )
  pbmc@meta.data$celltype1 <- plyr::mapvalues(x = as.integer(as.character(pbmc@meta.data$seurat_clusters)), 
                                             from = current.cluster.ids, 
                                             to = new.cluster.ids)
  return(pbmc)
}

pbmc <- cluster_anno(pbmc)
saveRDS(pbmc,'s8_anno.rds')
pbmc <- readRDS('s8_anno.rds')
#mycolors <- c(paletteer_d("ggthemes::Superfishel_Stone")[1:8],c("#A39FC9FF","#F3A546FF","#59A14FFF","#BAB0ACFF"))
mycolors <- paletteer_d("ggthemes::Tableau_10")
source('D:/self/my_functions_R/scRNA_plots.R')
orders <- c("Normal","pre-treat","post-treat")
ann_plot(pbmc,
         mycolors = mycolors,
         alpha.use = 0.3,
         groups = "group",
         ncols = 2,
         celltypes = "celltype1",
         orders = orders)

p <- DimPlot(pbmc,group.by="celltype",raster = F) + ggtitle(NULL) + scale_color_manual(values=mycolors)
p <- DimPlot(pbmc,label = T,label.size = 5,group.by="celltype1",raster = F) +
  #NoLegend() +
  ggtitle(NULL) +
  scale_color_manual(values=mycolors)
ggsave("anno1.pdf", plot = p, width = 6, height = 4)
ggsave("anno.pdf", plot = p, width = 7, height = 4.5)
pbmc$orders <- factor(pbmc$group,levels = c("Normal","pre-treat","post-treat"))
p <- DimPlot(pbmc,group.by="celltype1",split.by = "orders",raster = F,label = T,label.size = 4) + 
  ggtitle(NULL) +
  scale_color_manual(values=mycolors) +
  NoLegend()
ggsave("umap_anno_split.pdf", plot = p, width = 10, height = 3.5)

relative_percentage_plot <- function(metadata){
  s1 <- as.data.frame(table(metadata$group,metadata$celltype))
  s1 <- reshape2::dcast(s1,Var1 ~ Var2)
  rownames(s1) <- s1[,1]
  s1 <- s1[-1]
  portion2 <- apply(s1, 1, function(x)  x/sum(x))
  #portion2 <- t(apply(portion2, 1, function(x)  x/sum(x)))
  portion <- portion2*100
  portion <- round(portion,2)
  portion <- reshape2::melt(portion)
  portion$Var1 <- factor(portion$Var1, levels=unique(as.character(portion$Var1)))
  orders <- c('Normal','pre-treat','post-treat')
  #orders <- c("N1","N2","case1_pre","case1_post","case2_pre","case2_post","case3_pre","case3_post")
  ggplot(portion,aes(x=Var2,y=value,fill=Var1)) + 
    geom_bar(stat="identity") +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values=mycolors) +
    #scale_fill_brewer(palette = "Paired") +
    guides(fill=guide_legend(title="Cluster")) +
    theme_classic(base_size = 14) +
    theme(axis.text = element_text(size = 12)) +
    scale_x_discrete(limits = orders) +
    ylab("The Percentage of Cluster(%)") +
    xlab(NULL)
  
}
#mycolor <- c("#FED439FF","#87CEFA","#8A9197FF","#FFA500","#D5E4A2FF","#197EC0FF","#F05C3BFF","#46732EFF")
#metas <- pbmc@meta.data
#metas <- metas[metas$celltype!="Keratinocyte",]
ppp <- relative_percentage_plot(pbmc@meta.data)
ggsave("ratio.pdf", plot = ppp, width = 6, height = 5)
ggsave("ratio_1.pdf", plot = ppp, width = 6, height = 4)

DefaultAssay(pbmc)<-"RNA"
p2 <- DotPlot(pbmc, features = markers,group.by = "celltype") + 
  theme_bw(base_size = 14) + #coord_flip() +  
  theme(panel.grid = element_blank(),  
        axis.text.x=element_text(hjust = 1,angle=90,vjust=0.5)) + #angle=90,
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c('#330066','#336699','#66CC66','#FFCC33')) + 
  scale_y_discrete(limits =c("Sweat_gland_cell","Keratinocyte","Melanocyte","Immune","Mast","Neural","Fibroblast","Mural","LEC","VEC")) +
  labs(x=NULL,y=NULL)
DefaultAssay(pbmc)<-"integrated"
p2
graph2ppt(file="plots.pptx",width=7.11, height=6.34,append=TRUE)
ggsave("marker_dotplot.pdf", plot = p2, width = 7.57, height = 3.44)

####################### anno version2
cluster_name <- unique(pbmc@meta.data[,c("seurat_clusters","celltype")])
rownames(cluster_name) <- NULL
colnames(cluster_name) <- c("cluster","labels")
cluster_name <- cluster_name[order(cluster_name$cluster),]
#cluster_name <- read.csv('anno.csv')
new.cluster.ids <- cluster_name$labels
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc,new.cluster.ids)
saveRDS(pbmc,'pre_post_normal_anno_V2.rds')
sce.markers <- FindAllMarkers(pbmc, 
                              only.pos = TRUE, 
                              test.use = 't',
                              min.diff.pct = 0.2,
                              #slot="counts",
                              min.pct = 0.1, 
                              logfc.threshold = .1,
                              return.thresh=0.05)


################################## consistence of gene expression 
GeneExpInGroup <- function(ct){
  t.cells <- subset(pbmc, celltype == ct)
  Idents(t.cells) <- "group"
  avg.t.cells <- log1p(AverageExpression(t.cells, verbose = FALSE)$RNA)
  avg.t.cells <- as.data.frame(avg.t.cells)
  avg.t.cells$gene <- rownames(avg.t.cells)
  colnames(avg.t.cells) <- c('Normal','Post','Pre','gene')
  outf <- paste(ct,'exp_avg.csv',sep='_')
  write.csv(avg.t.cells,outf,quote = F,row.names = F)
  
  avg.t.cells$N_Pre <- abs(avg.t.cells$Normal - avg.t.cells$Pre)
  avg.t.cells$N_Post <- abs(avg.t.cells$Normal - avg.t.cells$Post)
  avg.t.cells$Post_Pre <- abs(avg.t.cells$Post - avg.t.cells$Pre)
  genes.to.label1 = avg.t.cells[order(avg.t.cells$N_Pre,decreasing = T),]$gene[1:10]
  genes.to.label2 = avg.t.cells[order(avg.t.cells$Post_Pre,decreasing = T),]$gene[1:10]
  genes.to.label3 = avg.t.cells[order(avg.t.cells$N_Post,decreasing = T),]$gene[1:10]

  p1 <- ggplot(avg.t.cells, aes(Normal, Pre)) + 
    geom_point() + 
    ggtitle(ct) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) 
  p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE)
  p2 <- ggplot(avg.t.cells, aes(Post,Pre)) + 
    geom_point() + 
    ggtitle(ct) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) 
  p2 <- LabelPoints(plot = p2, points = genes.to.label2, repel = TRUE)
  p3 <- ggplot(avg.t.cells, aes(Normal,Post)) + 
    geom_point() + 
    ggtitle(ct) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) 
  p3 <- LabelPoints(plot = p3, points = genes.to.label3, repel = TRUE)
  p <- cowplot::plot_grid(p1, p2, p3)
  return(p)
}

p1 <- GeneExpInGroup("VEC")
ggsave("VEC_exp_group.pdf", plot = p1, width = 9, height = 8.5)

################################## ECs recluster
pbmc <- readRDS('s8_anno.rds')
ecs <- subset(pbmc,celltype=="VEC")
ecs <- FindVariableFeatures(ecs, selection.method = "vst", nfeatures = 3000)
ecs <- ScaleData(ecs, verbose = FALSE,features = rownames(ecs))
ecs <- RunPCA(ecs, npcs = 100, verbose = FALSE)
#ecs <- RunUMAP(ecs, reduction = "pca", dims = 1:30)
ecs <- RunUMAP(ecs, 
               reduction = "pca", 
               metric = "euclidean",
               n.epochs=500,
               negative.sample.rate = 10L,
               #min.dist=0.2,
               #n.neighbors = 20,
               dims = 1:30
               )
saveRDS(ecs,'s8_VEC_dimention.rds')
#ecs <- RunTSNE(ecs, reduction = "pca", dims = 1:30,n.epochs=500,negative.sample.rate = 10L,min.dist=0.3)
ecs <- readRDS('s8_VEC_dimention.rds')
ecs <- FindNeighbors(ecs, k.param = 20, prune.SNN =1/15, reduction = "pca", dims = 1:30)
ecs <- FindClusters(ecs, random.seed = 1234,resolution = 0.6,n.iter=10)
p1 <- DimPlot(ecs,label = T,label.size = 6) + NoLegend()
gavg <- AverageExpression(ecs)
gavg <- data.frame(gavg$RNA)
write.csv(gavg,'VEC_gexp_avg.csv')
markers <- c(
  "SEMA3G", #arteriole
  "MKI67","TOP2A",#Cycling
  "PLVAP", #capillaries
  "ACKR1","SELE",#post-capillary venules
  "FBLN2",#arteriole/venules
  #"LYVE1","PROX1",#Lymph EC
  "SLC2A1"
)

markers <- c(
  "SEMA3G","ELN","FBLN5", #artery
  "BTNL9",#artery-capillary
  "SELE",#post-capillary venules
  "SELP","ACKR1", #vein
  "RGCC","SLC2A1",#capillary
  "MKI67","TOP2A"
)

DefaultAssay(ecs)<-"RNA"
for(g in markers){
  pp <- FeaturePlot(ecs,features = g,min.cutoff = "q10")
  outfile <- paste(g,'FeaturePlot.pdf',sep = '_')
  ggsave(outfile, plot = pp, width = 6, height = 4.5)
}

plot_density(ecs, markers,reduction = "umap") + p1
DotPlot(ecs,features = markers) + theme(axis.text.x = element_text(angle = 90)) +p1
VlnPlot(ecs, features = markers,ncol = 3,log = TRUE) + p1
DefaultAssay(ecs)<-"integrated"
sce.markers <- FindAllMarkers(ecs, 
                              only.pos = TRUE, 
                              #test.use = 't',
                              assay = "RNA",
                              min.diff.pct = 0.05,
                              slot="counts",
                              min.pct = 0.05, 
                              logfc.threshold = .15,
                              return.thresh=0.05)
write.csv(sce.markers,'VEC_ann_findallmarker.csv')
top3 <- sce.markers %>% group_by(cluster) %>% top_n(n = 11, wt = avg_log2FC)
p2 <- DoHeatmap(ecs, features = top3$gene,angle = 90,size=4,draw.lines = F,group.bar.height = 0.01,group.colors = mycolor)
p2
ggsave("ECs_ann_pheatmap.pdf", plot = p2, width = 9, height = 10)


cluster_anno <- function(pbmc){
  C1 = c(0,11,12,13,15) 
  C2 = c(2,6,10) 
  C3 = c(8,5) 
  C4 = c(1,3,7,17,18)
  C5 = c(4,14,16,9)
  current.cluster.ids <- c(C1,C2,C3,C4,C5)
  new.cluster.ids <- c(rep("C1",length(C1)),
                       rep("C2",length(C2)),
                       rep("C3",length(C3)),
                       rep("C4",length(C4)),
                       rep("C5",length(C5))
  )
  pbmc@meta.data$subtype <- plyr::mapvalues(x = as.integer(as.character(pbmc@meta.data$seurat_clusters)), 
                                             from = current.cluster.ids, 
                                             to = new.cluster.ids)
  return(pbmc)
}
ecs <- cluster_anno(ecs)
saveRDS(ecs,'s8_VEC_recluster.rds')
ecs <- readRDS('s8_VEC_recluster.rds')
p1 <- DimPlot(ecs,group.by = "subtype",label = T,label.size = 6,cols = paletteer_d("ggthemes::Tableau_10")) + NoLegend() + ggtitle(NULL)
ggsave("VEC_recluster.pdf", plot = p1, width = 6, height = 4.5)
p2 <- FeaturePlot(ecs,features = markers,min.cutoff = "q10",ncol = 3) + p1 
ggsave("EC_marker1.pdf", plot = p2, width = 9, height = 7)
ecs$order <- factor(ecs$group,levels = c('Normal','pre-treat','post-treat'))
p3 <- DimPlot(ecs,label = T,label.size = 5,split.by = "order",cols = paletteer_d("ggthemes::Tableau_10"),group.by = "subtype") + NoLegend() + ggtitle(NULL)
ggsave("VEC_split.pdf", plot = p3, width = 11, height = 4)

HemEC <- c("PECAM1","VWF","SELE","KDR","TEK","CDH5")
HemPEC <- c("PROM1","CD34","KDR","PECAM1","MCAM","CDH5","VWF")
HemMSC <- c("ENG","SAMSN1","LYN","THY1","ITGB1","PROM1","ACTA2")
HemSC <- c("THY1","PECAM1","FLT1","KDR","MCAM","NRP1")
HemPericyte <- c("PDGFRB","CSPG4","DES","CNN1","THY1")
showgene <- c("ENG","THY1","ITGB1","RGS5","CD36","FABP4","MYH11","CNN1","RERGL",
              "GRB2","SRC","KDR","FLT1")
showgene <- c("EPAS1","EGLN1","PPARA")

ecs <- readRDS('s8_VEC_recluster.rds')
DefaultAssay(ecs)<-"RNA"
p3 <- FeaturePlot(ecs,features = c('ARRB1','ARRB2'),ncol = 2,min.cutoff = "q10")
FeaturePlot(ecs,features = c('ARRB1','ARRB2'),ncol = 2,min.cutoff = "q10",split.by = "group")
ggsave("VEC_gene.pdf", plot = p3, width = 9, height = 7)
#plot_density(smc, HemPericyte,reduction = "umap") + p1
DefaultAssay(ecs)<-"integrated"

cluster_name <- unique(ecs@meta.data[,c("seurat_clusters","subtype")])
rownames(cluster_name) <- NULL
colnames(cluster_name) <- c("cluster","labels")
cluster_name <- cluster_name[order(cluster_name$cluster),]
#cluster_name <- read.csv('anno.csv')
new.cluster.ids <- cluster_name$labels
names(new.cluster.ids) <- levels(ecs)
ecs <- RenameIdents(ecs,new.cluster.ids)
saveRDS(ecs,'VEC_recluster.V2.rds')
gavg <- AverageExpression(ecs)
gavg <- data.frame(gavg$RNA)
write.csv(gavg,'VEC_ann_gexp_avg.csv')


mycolor <- paletteer_d("ggthemes::Tableau_10")
relative_percentage_plot <- function(metadata){
  s1 <- as.data.frame(table(metadata$group,metadata$subtype))
  s1 <- reshape2::dcast(s1,Var1 ~ Var2)
  rownames(s1) <- s1[,1]
  s1 <- s1[-1]
  portion2 <- apply(s1, 1, function(x)  x/sum(x))
  #portion2 <- t(apply(portion2, 1, function(x)  x/sum(x)))
  portion <- portion2*100
  portion <- round(portion,2)
  portion <- reshape2::melt(portion)
  portion$Var1 <- factor(portion$Var1, levels=unique(as.character(portion$Var1)))
  orders <- c("Normal","pre-treat","post-treat")
  ggplot(portion,aes(x=Var2,y=value,fill=Var1)) + 
    geom_bar(stat="identity") +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values=mycolor) +
    #scale_fill_brewer(palette = "Paired") +
    guides(fill=guide_legend(title="Cluster")) +
    theme_classic() +
    theme(axis.text = element_text(size = 12)) +
    scale_x_discrete(limits = orders) +
    ylab("The Percentage of Cluster(%)") +
    xlab(NULL)
  
}
p1 <- relative_percentage_plot(ecs)
ggsave("VEC_ratio.pdf", plot = p1, width = 5, height = 4)

library(patchwork)
source('D:/self/MyscFuncPlots_doc/GeneExp_GroupCompare_plot.R')
tfs1 <- c('ADRB1','ADRB2')
mycolors <- c("#5773CC","#FFB900","#5CB85C")
plist2 <- GeneExpComparePlot(obj = ecs,
                             GroupBy = "group",
                             TargetGene = 'ARRB2',
                             CellType = "ALL",
                             Orders = c('Normal',"pre-treat","post-treat"),
                             mycolors = mycolors,
                             ctypecol = "subtype",
                             stat = "Y",
                             method = "wilcox.test",
                             hide = "Y")
pp <- wrap_plots(plist2,ncol = 2)


###################### smooth muscle recluster
pbmc <- readRDS('s8_anno.rds')
smc <- subset(pbmc,celltype=="Mural")
rm(pbmc)
gc()
smc <- FindVariableFeatures(smc, selection.method = "vst", nfeatures = 3000)
smc <- ScaleData(smc, verbose = FALSE,features = rownames(smc))
smc <- RunPCA(smc, npcs = 100, verbose = FALSE)
smc <- RunUMAP(smc, 
               reduction = "pca", 
               #metric = "euclidean",
               n.epochs=500,
               negative.sample.rate = 10L,
               #min.dist=0.1,
               #n.neighbors = 20,
               dims = 1:30
)
saveRDS(smc,'s8_mural_dimention.rds')
smc <- readRDS('s8_mural_dimention.rds')

#ecs <- RunTSNE(ecs, reduction = "pca", dims = 1:30,n.epochs=500,negative.sample.rate = 10L,min.dist=0.3)
smc <- FindNeighbors(smc, k.param = 20, prune.SNN =1/15, reduction = "pca", dims = 1:30)
smc <- FindClusters(smc, random.seed = 1234,resolution = 0.6,n.iter=10)
p1 <- DimPlot(smc,label = T,label.size = 6) + NoLegend()

HemPEC <- c("PROM1","CD34","KDR","PECAM1","MCAM","CDH5","VWF")
HemMSC <- c("ENG","GRB2","SRC","THY1","ITGB1","PROM1","ACTA2")
HemSC <- c("THY1","PECAM1","FLT1","KDR","MCAM","NRP1")
HemPericyte <- c("PDGFRB","CSPG4","DES","CNN1","THY1")
markers <- c("CSPG4","PDGFRB","ANPEP",
             "MCAM","DES","ALPL","ACTA2","3G5")

markers <- c(
  "RERGL",#vSMC
  "RGS5","FABP4","CD36","PDGFRB","CNN1",#Pericyte
  "MKI67","STMN1"
)

DefaultAssay(smc)<-"RNA"
for(g in markers){
  pp <- FeaturePlot(smc,features = g,min.cutoff = "q10")
  outfile <- paste(g,'FeaturePlot.pdf',sep = '_')
  ggsave(outfile, plot = pp, width = 6, height = 4.5)
}
p3 <- FeaturePlot(smc,features = markers,ncol = 3,min.cutoff = "q10") 
#plot_density(smc, HemPericyte,reduction = "umap") + p1
DefaultAssay(smc)<-"integrated"

sce.markers <- FindAllMarkers(smc, 
                              only.pos = TRUE, 
                              #test.use = 't',
                              assay = "RNA",
                              min.diff.pct = 0.05,
                              slot="counts",
                              min.pct = 0.05, 
                              logfc.threshold = .1,
                              return.thresh=0.05)
write.csv(sce.markers,'mural_ann_findallmarker.csv')
sce.markers <- read.csv('mural_findallmarker.csv',row.names = 1)
top3 <- sce.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
p2 <- DoHeatmap(smc, features = top3$gene,angle = 90,size=4,draw.lines = F,group.bar.height = 0.01,group.colors = mycolors)
ggsave("Mural_ann_pheatmap.pdf", plot = p2, width = 6.5, height = 8.5)


cluster_anno <- function(pbmc){
  C1 = c(12,13) 
  C2 = c(9) 
  C3 = c(0,1,5,14,15,2,3,4,6,8,10) 
  C4 = c(7) 
  C5 = c(11)
  current.cluster.ids <- c(C1,C2,C3,C4,C5)
  new.cluster.ids <- c(rep("C1",length(C1)),
                       rep("C2",length(C2)),
                       rep("C3",length(C3)),
                       rep("C4",length(C4)),
                       rep("C5",length(C5))
  )
  pbmc@meta.data$subtype <- plyr::mapvalues(x = as.integer(as.character(pbmc@meta.data$seurat_clusters)), 
                                            from = current.cluster.ids, 
                                            to = new.cluster.ids)
  return(pbmc)
}
smc <- cluster_anno(smc)
saveRDS(smc,'Mural_recluster.rds')
write.csv(smc@meta.data,'Mural_metadata.csv')
smc <- readRDS('Mural_recluster.rds')
mycolors <- as.character(paletteer_d("ggsci::default_jama")[c(2:3,5:7)])
p1 <- DimPlot(smc,group.by = "subtype",label = T,label.size = 6,cols = mycolors) + NoLegend() + ggtitle(NULL)
ggsave("Mural_recluster.pdf", plot = p1, width = 6, height = 4.5)
p2 <- FeaturePlot(smc,features = markers,min.cutoff = "q10",ncol = 3) 
ggsave("EC_marker1.pdf", plot = p2, width = 9, height = 7)
smc$order <- factor(smc$group,levels = c('Normal','pre-treat','post-treat'))
p3 <- DimPlot(smc,label = T,label.size = 5,split.by = "order",cols = mycolors,group.by = "subtype") + NoLegend() + ggtitle(NULL)
ggsave("Mural_split.pdf", plot = p3, width = 11, height = 4)

gavg <- AverageExpression(smc)
gavg <- data.frame(gavg$RNA)
write.csv(gavg,'mural_ann_gexp_avg.csv')

pp <- FeaturePlot(smc,features = "CENPF",min.cutoff = "q10")
ggsave("CENPF_umap.pdf", plot = pp, width = 5.5, height = 4)

################################## Fibroblast recluster
pbmc <- readRDS('s8_anno.rds')
ecs <- subset(pbmc,celltype=="Fibroblast")
ecs <- FindVariableFeatures(ecs, selection.method = "vst", nfeatures = 3000)
ecs <- ScaleData(ecs, verbose = FALSE,features = rownames(ecs))
ecs <- RunPCA(ecs, npcs = 100, verbose = FALSE)
#ecs <- RunUMAP(ecs, reduction = "pca", dims = 1:30)
ecs <- RunUMAP(ecs, 
               reduction = "pca", 
               metric = "euclidean",
               n.epochs=500,
               negative.sample.rate = 10L,
               #min.dist=0.2,
               #n.neighbors = 20,
               dims = 1:30
)
saveRDS(ecs,'s8_Fibro_dimention.rds')
#ecs <- RunTSNE(ecs, reduction = "pca", dims = 1:30,n.epochs=500,negative.sample.rate = 10L,min.dist=0.3)
ecs <- readRDS('s8_Fibro_dimention.rds')
ecs <- FindNeighbors(ecs, k.param = 20, prune.SNN =1/15, reduction = "pca", dims = 1:30)
ecs <- FindClusters(ecs, random.seed = 1234,resolution = 0.3,n.iter=10)
p1 <- DimPlot(ecs,label = T,label.size = 6) + NoLegend()
gavg <- AverageExpression(ecs)
gavg <- data.frame(gavg$RNA)
write.csv(gavg,'Fibro_gexp_avg.csv')
markers <- c(
  "CD74","KRT14",#C2  
  "FBLN1","PDGFRL",#C3
  "RGS5","MCAM","FABP4","CCL21",#C1
  "IGFBP2","SPON2","HOPX", #C4
  "CLDN1","FOXS1","SFRP4",#C6
  "TNN","ASPN",#C7
  "APOE","SFRP2" #C5
)
FeaturePlot(ecs,features = "GPC3",min.cutoff = "q10")
DefaultAssay(ecs)<-"RNA"
for(g in markers){
  pp <- FeaturePlot(ecs,features = g,min.cutoff = "q10")
  outfile <- paste(g,'FeaturePlot.pdf',sep = '_')
  ggsave(outfile, plot = pp, width = 6, height = 4.5)
}

plot_density(ecs, markers,reduction = "umap") + p1
DotPlot(ecs,features = markers) + theme(axis.text.x = element_text(angle = 90)) +p1
VlnPlot(ecs, features = markers,ncol = 3,log = TRUE) + p1
DefaultAssay(ecs)<-"integrated"
sce.markers <- FindAllMarkers(ecs, 
                              only.pos = TRUE, 
                              #test.use = 't',
                              assay = "RNA",
                              min.diff.pct = 0.05,
                              slot="counts",
                              min.pct = 0.05, 
                              logfc.threshold = .1,
                              return.thresh=0.05)
write.csv(sce.markers,'Fibro_ann_findallmarker.csv')
sce.markers <- sce.markers[sce.markers$p_val_adj<=0.05,]
top3 <- sce.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
p2 <- DoHeatmap(ecs, features = top3$gene,angle = 90,size=4,draw.lines = F,group.bar.height = 0.01,group.colors = mycolor)
p2
ggsave("Fibro_ann_pheatmap.pdf", plot = p2, width = 10, height = 10)


marker_cosg <- cosg(
  ecs,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  n_genes_user=10)

pdf('Fibro_res0.01_degs.pdf',width = 5, height = 4.2)
degs <- as.data.frame(marker_cosg[["names"]])
for(i in colnames(degs)){
  for(g in degs[,i]){
    pp <- FeaturePlot(ecs,features = g,min.cutoff = "q10")
    print(pp)
  }
}
dev.off()
write.csv(marker_cosg[["names"]],'Fibro_degs.csv',row.names = F)

cluster_anno <- function(pbmc){
  C1 = c(6) 
  C2 = c(7) 
  C3 = c(4) 
  C4 = c(5,0,1,2,11,9)
  C5 = c(3,8,10)
  current.cluster.ids <- c(C1,C2,C3,C4,C5)
  new.cluster.ids <- c(rep("C1",length(C1)),
                       rep("C2",length(C2)),
                       rep("C3",length(C3)),
                       rep("C4",length(C4)),
                       rep("C5",length(C5))
  )
  pbmc@meta.data$subtype <- plyr::mapvalues(x = as.integer(as.character(pbmc@meta.data$seurat_clusters)), 
                                            from = current.cluster.ids, 
                                            to = new.cluster.ids)
  return(pbmc)
}
ecs <- cluster_anno(ecs)
p1 <- DimPlot(ecs,group.by = "subtype",label = T,label.size = 6,cols = paletteer_d("ggthemes::Tableau_10")) + NoLegend() + ggtitle(NULL)
saveRDS(ecs,'s8_Fibro_recluster.rds')
ecs <- readRDS('s8_Fibro_recluster.rds')
ggsave("Fibro_recluster.pdf", plot = p1, width = 6, height = 4.5)
ecs$order <- factor(ecs$group,levels = c('Normal','pre-treat','post-treat'))
p3 <- DimPlot(ecs,label = T,label.size = 5,split.by = "order",cols = paletteer_d("ggthemes::Tableau_10"),group.by = "subtype") + NoLegend() + ggtitle(NULL)
ggsave("Fibro_split.pdf", plot = p3, width = 11, height = 4)


ecs <- readRDS('s8_VEC_recluster.rds')
DefaultAssay(ecs)<-"RNA"
p3 <- FeaturePlot(ecs,features = "EPAS1",ncol = 2,min.cutoff = "q10")
FeaturePlot(ecs,features = "EPAS1",ncol = 2,min.cutoff = "q10",split.by = "group")
ggsave("VEC_gene.pdf", plot = p3, width = 9, height = 7)
#plot_density(smc, HemPericyte,reduction = "umap") + p1
DefaultAssay(ecs)<-"integrated"

cluster_name <- unique(ecs@meta.data[,c("seurat_clusters","subtype")])
rownames(cluster_name) <- NULL
colnames(cluster_name) <- c("cluster","labels")
cluster_name <- cluster_name[order(cluster_name$cluster),]
#cluster_name <- read.csv('anno.csv')
new.cluster.ids <- cluster_name$labels
names(new.cluster.ids) <- levels(ecs)
ecs <- RenameIdents(ecs,new.cluster.ids)
gavg <- AverageExpression(ecs)
gavg <- data.frame(gavg$RNA)
write.csv(gavg,'Fibro_ann_gexp_avg.csv')


mycolor <- paletteer_d("ggthemes::Tableau_10")
relative_percentage_plot <- function(metadata){
  s1 <- as.data.frame(table(metadata$group,metadata$subtype))
  s1 <- reshape2::dcast(s1,Var1 ~ Var2)
  rownames(s1) <- s1[,1]
  s1 <- s1[-1]
  portion2 <- apply(s1, 1, function(x)  x/sum(x))
  #portion2 <- t(apply(portion2, 1, function(x)  x/sum(x)))
  portion <- portion2*100
  portion <- round(portion,2)
  portion <- reshape2::melt(portion)
  portion$Var1 <- factor(portion$Var1, levels=unique(as.character(portion$Var1)))
  orders <- c("Normal","pre-treat","post-treat")
  ggplot(portion,aes(x=Var2,y=value,fill=Var1)) + 
    geom_bar(stat="identity") +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values=mycolor) +
    #scale_fill_brewer(palette = "Paired") +
    guides(fill=guide_legend(title="Cluster")) +
    theme_classic() +
    theme(axis.text = element_text(size = 12)) +
    scale_x_discrete(limits = orders) +
    ylab("The Percentage of Cluster(%)") +
    xlab(NULL)
  
}
p1 <- relative_percentage_plot(ecs)
ggsave("Fibro_ratio.pdf", plot = p1, width = 5, height = 4)


################################## Immune recluster
pbmc <- readRDS('s8_anno.rds')
ecs <- subset(pbmc,celltype=="Immune")
ecs <- FindVariableFeatures(ecs, selection.method = "vst", nfeatures = 3000)
ecs <- ScaleData(ecs, verbose = FALSE,features = rownames(ecs))
ecs <- RunPCA(ecs, npcs = 100, verbose = FALSE)
#ecs <- RunUMAP(ecs, reduction = "pca", dims = 1:30)
ecs <- RunUMAP(ecs, 
               reduction = "pca", 
               metric = "euclidean",
               n.epochs=500,
               negative.sample.rate = 10L,
               #min.dist=0.2,
               #n.neighbors = 20,
               dims = 1:30
)
saveRDS(ecs,'s8_Immune_dimention.rds')
#ecs <- RunTSNE(ecs, reduction = "pca", dims = 1:30,n.epochs=500,negative.sample.rate = 10L,min.dist=0.3)
ecs <- readRDS('s8_Immune_dimention.rds')
ecs <- FindNeighbors(ecs, k.param = 20, prune.SNN =1/15, reduction = "pca", dims = 1:30)
ecs <- FindClusters(ecs, random.seed = 1234,resolution = 0.3,n.iter=10)
p1 <- DimPlot(ecs,label = T,label.size = 6) + NoLegend()
gavg <- AverageExpression(ecs)
gavg <- data.frame(gavg$RNA)
write.csv(gavg,'Immune_gexp_avg.csv')
markers <- c(
  "CD68", #Macro
  "MKI67",#Cycling
  "CD3D","CD3E", #T
  "GNLY","GZMB",#NK
  "CD79A","MS4A1",#B
  "CD14","S100A8",#Mono
  "CD1C","CLEC4C",#pDC
  "CLEC9A"#cDC
)

DefaultAssay(ecs)<-"RNA"
for(g in markers){
  pp <- FeaturePlot(ecs,features = g,min.cutoff = "q10")
  outfile <- paste(g,'FeaturePlot.pdf',sep = '_')
  ggsave(outfile, plot = pp, width = 6, height = 4.5)
}

plot_density(ecs, markers,reduction = "umap") + p1
DotPlot(ecs,features = markers) + theme(axis.text.x = element_text(angle = 90)) +p1
VlnPlot(ecs, features = markers,ncol = 3,log = TRUE) + p1
DefaultAssay(ecs)<-"integrated"
sce.markers <- FindAllMarkers(ecs, 
                              only.pos = TRUE, 
                              #test.use = 't',
                              assay = "RNA",
                              min.diff.pct = 0.05,
                              slot="counts",
                              min.pct = 0.05, 
                              logfc.threshold = .15,
                              return.thresh=0.05)
write.csv(sce.markers,'Immune_ann_findallmarker.csv')
top3 <- sce.markers %>% group_by(cluster) %>% top_n(n = 11, wt = avg_log2FC)
p2 <- DoHeatmap(ecs, features = top3$gene,angle = 90,size=4,draw.lines = F,group.bar.height = 0.01,group.colors = mycolor)
p2
ggsave("Immune_ann_pheatmap.pdf", plot = p2, width = 9, height = 10)


cluster_anno <- function(pbmc){
  C1 = c(0,5) 
  C2 = c(8) 
  C3 = c(1,2) 
  C4 = c(3)
  C5 = c(7)
  C6 = c(9)
  C7 = c(6)
  C8 = c(4)
  current.cluster.ids <- c(C1,C2,C3,C4,C5,C6,C7,C8)
  new.cluster.ids <- c(rep("T",length(C1)),
                       rep("NK",length(C2)),
                       rep("Macro",length(C3)),
                       rep("cDC",length(C4)),
                       rep("pDC",length(C5)),
                       rep("Mono",length(C6)),
                       rep("B",length(C7)),
                       rep("Others",length(C8))
  )
  pbmc@meta.data$subtype <- plyr::mapvalues(x = as.integer(as.character(pbmc@meta.data$seurat_clusters)), 
                                            from = current.cluster.ids, 
                                            to = new.cluster.ids)
  return(pbmc)
}
ecs <- cluster_anno(ecs)
p1 <- DimPlot(ecs,group.by = "subtype",label = T,label.size = 6,cols = paletteer_d("ggthemes::Tableau_10")) + NoLegend() + ggtitle(NULL)
saveRDS(ecs,'s8_Immune_recluster.rds')
ecs <- readRDS('s8_Immune_recluster.rds')
ggsave("Immune_recluster.pdf", plot = p1, width = 6, height = 4.5)
p2 <- FeaturePlot(ecs,features = markers,min.cutoff = "q10",ncol = 3) + p1 
ggsave("EC_marker1.pdf", plot = p2, width = 9, height = 7)
ecs$order <- factor(ecs$group,levels = c('Normal','pre-treat','post-treat'))
p3 <- DimPlot(ecs,label = T,label.size = 5,split.by = "order",cols = paletteer_d("ggthemes::Tableau_10"),group.by = "subtype") + NoLegend() + ggtitle(NULL)
ggsave("Immune_split.pdf", plot = p3, width = 11, height = 4)

showgene <- c("EPAS1","EGLN1","PPARA")

ecs <- readRDS('s8_Immune_recluster.rds')
DefaultAssay(ecs)<-"RNA"
p3 <- FeaturePlot(ecs,features = "PROM1",ncol = 2,min.cutoff = "q10")
FeaturePlot(ecs,features = "EPAS1",ncol = 2,min.cutoff = "q10",split.by = "group")
ggsave("VEC_gene.pdf", plot = p3, width = 9, height = 7)
#plot_density(smc, HemPericyte,reduction = "umap") + p1
DefaultAssay(ecs)<-"integrated"

cluster_name <- unique(ecs@meta.data[,c("seurat_clusters","subtype")])
rownames(cluster_name) <- NULL
colnames(cluster_name) <- c("cluster","labels")
cluster_name <- cluster_name[order(cluster_name$cluster),]
#cluster_name <- read.csv('anno.csv')
new.cluster.ids <- cluster_name$labels
names(new.cluster.ids) <- levels(ecs)
ecs <- RenameIdents(ecs,new.cluster.ids)
saveRDS(ecs,'VEC_recluster.V2.rds')
gavg <- AverageExpression(ecs)
gavg <- data.frame(gavg$RNA)
write.csv(gavg,'VEC_ann_gexp_avg.csv')


mycolor <- paletteer_d("ggthemes::Tableau_10")
relative_percentage_plot <- function(metadata){
  s1 <- as.data.frame(table(metadata$group,metadata$subtype))
  s1 <- reshape2::dcast(s1,Var1 ~ Var2)
  rownames(s1) <- s1[,1]
  s1 <- s1[-1]
  portion2 <- apply(s1, 1, function(x)  x/sum(x))
  #portion2 <- t(apply(portion2, 1, function(x)  x/sum(x)))
  portion <- portion2*100
  portion <- round(portion,2)
  portion <- reshape2::melt(portion)
  portion$Var1 <- factor(portion$Var1, levels=unique(as.character(portion$Var1)))
  orders <- c("Normal","pre-treat","post-treat")
  ggplot(portion,aes(x=Var2,y=value,fill=Var1)) + 
    geom_bar(stat="identity") +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values=mycolor) +
    #scale_fill_brewer(palette = "Paired") +
    guides(fill=guide_legend(title="Cluster")) +
    theme_classic() +
    theme(axis.text = element_text(size = 12)) +
    scale_x_discrete(limits = orders) +
    ylab("The Percentage of Cluster(%)") +
    xlab(NULL)
  
}
p1 <- relative_percentage_plot(ecs)
ggsave("Immune_ratio.pdf", plot = p1, width = 5, height = 4)

############################## copykat
obj <- readRDS('s8_anno_AddSubType_ck.rds')
ck <- obj@meta.data[,"copykat.pred",drop=F]
ck <- cbind(cell=rownames(ck),ck)
obj <- readRDS('s8_VEC_recluster.rds')
adds <- cbind(cell=rownames(obj@meta.data),obj@meta.data)
adds <- merge(adds,ck)
adds <- adds[,c("cell","copykat.pred")]
rownames(adds) <- adds[,1]
adds <- adds[-1]
obj <- AddMetaData(obj,adds)
saveRDS(obj,'s8_VEC_recluster_ck.rds')

alpha.use <- 0.3
p1 <- DimPlot(obj,group.by = "copykat.pred",raster = F) + ggtitle(NULL)
p1$layers[[1]]$mapping$alpha <- alpha.use
p1 <- p1 + scale_alpha_continuous(range=alpha.use,guide="none")
p1

obj$orders <- factor(obj$group,levels = c('Normal','pre-treat','post-treat'))
p1 <- DimPlot(obj,group.by = "copykat.pred",raster = F,split.by = 'orders') + ggtitle(NULL)
p1$layers[[1]]$mapping$alpha <- alpha.use
p1 <- p1 + scale_alpha_continuous(range=alpha.use,guide="none")
p1


relative_percentage_plot <- function(metadata,group){
  metadata <- metadata[metadata$group==group,]
  s1 <- as.data.frame(table(metadata$subtype,metadata$copykat.pred))
  s1 <- reshape2::dcast(s1,Var1 ~ Var2)
  rownames(s1) <- s1[,1]
  s1 <- s1[-1]
  portion2 <- apply(s1, 1, function(x)  x/sum(x))
  #portion2 <- t(apply(portion2, 1, function(x)  x/sum(x)))
  portion <- portion2*100
  portion <- round(portion,2)
  portion <- reshape2::melt(portion)
  portion$Var1 <- factor(portion$Var1, levels=unique(as.character(portion$Var1)))
  #orders <- c('Normal','pre-treat','post-treat')
  ggplot(portion,aes(x=Var2,y=value,fill=Var1)) + 
    geom_bar(stat="identity") +
    scale_y_continuous(expand = c(0,0)) +
    #scale_fill_manual(values=mycolor) +
    #scale_fill_brewer(palette = "Paired") +
    guides(fill=guide_legend(title="Cluster")) +
    theme_classic() +
    theme(axis.text = element_text(size = 12)) +
    #scale_x_discrete(limits = orders) +
    ylab("The Percentage of Cluster(%)") +
    xlab(NULL)
  
}
p1 <- relative_percentage_plot(obj@meta.data,"post-treat")

############################# pheatmap ourself
source("D:/self/my_functions_R/pheatmap_add_gene_name.R")

obj <- readRDS('Rat_S4_anno.rds')
Idents(obj) <- "celltype"
obj.aver <- AverageExpression(obj)
obj.aver <- data.frame(obj.aver$RNA)
obj.aver <- obj.aver[which(rowSums(obj.aver) > 0.001),]
obj.aver <- obj.aver[!grepl("\\.", rownames(obj.aver)),]
obj.aver <- obj.aver[!grepl("RPS", rownames(obj.aver)),]
write.csv(obj.aver,'allgene_ann_avg.csv',quote = F)
obj.aver <- read.csv('allgene_ann_avg.csv',row.names = 1)
shows <- sce.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
obj.aver <- obj.aver[,c("C5","C2","C4","C7",
                          "C1","C6","C3")]
p1 <- pheatmap(obj.aver,
               cluster_cols = F,
               color = paletteer_c("viridis::viridis", 30),
               clustering_distance_rows = "euclidean",
               clustering_method = "average",
               #show_rownames = F,
               scale = "row")
p2 <- add.flag(p1,
               kept.labels = shows$gene,
               repel.degree = 0.2)
p2 <- as.ggplot(p2)
ggsave("Mural_gene_pheatmap.pdf", plot = p2, width = 6, height = 7)

################################ add subtype to total metadata
pbmc <- readRDS('s8_anno.rds')
m1 <- read.csv('VEC_metadata.csv',row.names = 1)
m2 <- read.csv('Mural_metadata.csv',row.names = 1)
m3 <- read.csv('Fibro_metadata.csv',row.names = 1)
m1 <- m1[,"subtype_1",drop=F]
m2 <- m2[,"subtype_1",drop=F]
m3 <- m3[,"subtype_1",drop=F]
m <- Reduce(rbind, list(m1,m2,m3))
metadata <- pbmc@meta.data
metadata <- merge(metadata,m,by="row.names",all = T)
meta2 <- metadata[complete.cases(metadata),]
meta1 <- anti_join(metadata,meta2,by="Row.names")
meta1$subtype_1 <- meta1$celltype
metadata <- rbind(meta1,meta2)
metadata <- metadata[,c("Row.names","subtype_1")]
rownames(metadata) <- metadata[,1]
metadata <- metadata[,"subtype_1",drop=F]
pbmc <- AddMetaData(pbmc,metadata)
saveRDS(pbmc,'s8_anno_AddSubType.rds')

################################ cell propotion statistic
ecs <- readRDS('s8_VEC_recluster.rds')
s1 <- as.data.frame(table(ecs@meta.data$orig.ident,ecs@meta.data$subtype))
s1 <- reshape2::dcast(s1,Var1 ~ Var2)
rownames(s1) <- s1[,1]
s1 <- s1[-1]
portion2 <- apply(s1, 1, function(x)  x/sum(x))
portion <- portion2*100
portion <- round(portion,2)
portion <- reshape2::melt(portion)
s1 <- portion[portion$Var1=="C1",]
s1 <- s1[1:6,]
s1$group <- c("post-treat","pre-treat","post-treat","pre-treat","post-treat","pre-treat")
s1$group <- factor(s1$group,levels = c("pre-treat","post-treat"))
p <- ggplot(s1,aes(x=group,y=value,fill=group)) + 
  geom_bar(stat="summary",fun=mean) +
  coord_cartesian(ylim=c(0,65)) + 
  scale_y_continuous(expand = c(0,0)) +
  stat_compare_means(method="t.test",
                     paired = T,
                     comparison=list(c("post-treat","pre-treat")),
                     label.y = 55,
                     label="p.signif") +
  scale_fill_manual(values=c("#5773CC","#FFB900")) +
  theme_classic() +
  ylab("The Percentage of Cluster(%)") +
  theme(legend.position = "none") +
  theme(axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12)) +
  xlab(NULL) 
ggsave("VEC_C1_cellratio_compare.pdf", plot = p, width = 4, height = 4)  

#################################### target gene exp in groups
obj <- readRDS('s8_VEC_recluster.rds')
source("D:/self/my_functions_R/06_GeneExp_GroupCompare_plot.R")
tgenes <- c("IL6","STAT3")
APLN_TF <- c("ERG","ETS1","PRDM1","SOX18","TBX3","ZEB1")
#ecs <- subset(obj,group=="Normal",invert=T)
p1 <- GeneExp_Compare(obj = obj,
                      GroupBy = "group",
                      TargetGene = tgenes,
                      CellType = "ALL",
                      Orders = c('Normal',"pre-treat","post-treat"),
                      cols = c("#FF9E4A","#5773CCFF","#FFB900FF"))
p1
########################### monocle2
library(monocle)
library(ggpubr)
library(paletteer)
library(ggridges)

ctypes <- c("PC_1","PC_2","PC_3","PC_4","vSMC","Fibro_1","Fibro_2","Fibro_3","Fibro_4","Fibro_5")
cds <- readRDS("monocle2.rds")
pData(cds)$orders1 <- factor(pData(cds)$subtype_1,levels = rev(ctypes))
cols <- as.character(paletteer_d("ggthemes::Tableau_10"))
cols <- as.character(paletteer_d("ggthemes::Classic_Cyclic"))
p1 <- plot_cell_trajectory(cds,
                           show_branch_points=F,
                           cell_size = 0.5)
p3 <- plot_cell_trajectory(cds, 
                           color_by = "orders1",
                           cell_size = 0.65,
                           show_branch_points=F) + 
  scale_color_manual(values = cols)
ggsave("monocle2_3.pdf", plot = p3, width = 5, height = 4)
p7 <- plot_cell_trajectory(cds, 
                           color_by = "orders1",
                           cell_size = 0.65,
                           show_branch_points=F) + 
  scale_color_manual(values = cols) + 
  facet_wrap(~group,nrow=1)
ggsave("monocle2_7.pdf", plot = p7, width = 12, height = 4)

df <- pData(cds)
df$orders <- factor(df$group,levels = c("Normal","pre-treat","post-treat"))
p6 <- ggplot(df, aes(Pseudotime, colour = subtype_1, fill=subtype_1)) +
  geom_density(bw=0.5,size=1,alpha = 0.5) +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) +
  facet_wrap(~orders) +
  theme_classic2()
ggsave("monocle2_6.pdf", plot = p6, width = 14, height = 4)

df$colorder <- factor(df$subtype_1,levels = rev(ctypes))
p8 <- ggplot(df, aes(x = Pseudotime, y = colorder, fill = colorder)) +
  geom_density_ridges() +
  theme_ridges() + 
  facet_wrap(~orders) +
  scale_fill_manual(values=cols) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(title=element_text(size=14)) +
  xlab("Pseudotime")+
  ylab("celltype")
ggsave("monocle2_8.pdf", plot = p8, width = 10, height = 6)


s.genes <- c("ATF4","PROCR","APLN","APLNR","FBLN2")
pp1 <- plot_genes_in_pseudotime(cds[s.genes,], 
                         color_by = "State") 
ggsave("genes.pdf", plot = pp1, width = 4, height = 6)
pre <- subset(cds,group=="pre-treat")
pre <- cds[,rownames(pre)]
post <- subset(pData(cds),group=="post-treat")
post <- cds[,rownames(post)]
pp2 <- plot_genes_in_pseudotime(pre[s.genes,], 
                                color_by = "State") 
ggsave("genes_pre-treat.pdf", plot = pp2, width = 4, height = 6)
pp3 <- plot_genes_in_pseudotime(post[s.genes,], 
                                color_by = "State") 
ggsave("genes_post-treat.pdf", plot = pp3, width = 4, height = 6)

######################################## C1 pre-post DEGs and KEGG
library(clusterProfiler)

obj <- readRDS('s8_VEC_recluster.rds')
obj <- subset(obj,group=="Normal",invert=T)
Idents(obj) <- "group"
degs <- FindMarkers(obj, 
                    slot = "count",
                    assay = "RNA",
                    test.use = "t",
                    ident.1 = "post-treat", 
                    ident.2 = "pre-treat",
                    logfc.threshold = 0,
                    min.pct = 0.05)
degs <- cbind(genes=rownames(degs),degs)
logfc <- 0.2
degs$group <-  ifelse((degs$p_val > 0.05) | (abs(degs$avg_log2FC) < logfc) , "not-sig",
                      ifelse(degs$avg_log2FC > logfc, "up", "down"))
write.csv(degs,'VEC_DEGs.csv',row.names = F)
degs_sig <- degs[!degs$group %in% "not-sig",]
eg = bitr(degs_sig$genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
id = as.character(eg[,2])
kk <- enrichKEGG(gene = id,
                 organism = 'hsa', 
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 1)
kk_1 <- DOSE::setReadable(kk, OrgDb="org.Hs.eg.db", keyType='ENTREZID')
write.csv(kk_1@result,'VEC_KEGG.csv')
keggs <- dplyr::filter(kk_1@result, grepl('APLN', geneID))

############## plot target pathway
library(pathview)
eg_1 = bitr(degs$genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
degs_1 <- merge(degs,eg_1,by.x="genes",by.y="SYMBOL")
degs_1 <- degs_1[!duplicated(degs_1$genes),]
genelist <- as.numeric(degs_1$avg_log2FC) 
names(genelist) <- degs_1$ENTREZID
select_pathway <- "hsa04371"
pathview(gene.data     = genelist,
         pathway.id    = select_pathway,
         species       = 'hsa' ,  
         map.null      = FALSE ,
         kegg.native   = F,# TRUE输出完整pathway的png文件，F输出基因列表的pdf文件 
         new.signature = F, #pdf是否显示pathway标注
         limit         = list(gene=0.2, cpd=1) #图例color bar范围调整 
)

######################################## GSEA by Pi
# BiocManager::install("Pi")
# library(Pi)
library(Seurat)
library(GseaVis)
library(clusterProfiler)
#library(dplyr)

gmts <- read.gmt('D:/self/database/msigdb/msigdb.v7.5.1.symbols_1.gmt')
# obj <- readRDS('s8_VEC_recluster.rds')
# obj <- subset(obj,group=="Normal",invert=T)
degs <- read.csv('VEC_C1_DEGs.csv')
#eg <- bitr(degs$genes,fromType="SYMBOL",toType="ENTREZID",OrgDb = "org.Hs.eg.db")
#degs <- merge(degs,eg,by.x="genes",by.y="SYMBOL")
FCgenelist <- degs$avg_log2FC
names(FCgenelist) <- degs$genes
FCgenelist <- na.omit(FCgenelist)
FCgenelist <- sort(FCgenelist, decreasing = T)  
R.utils::setOption( "clusterProfiler.download.method",'auto' )
gsea <- GSEA(FCgenelist,
             minGSSize = 5,
             maxGSSize = 500,
             pvalueCutoff = 1,
             TERM2GENE = gmts)
write.csv(gsea@result,'VEC_C1_GSEA_Msigdb.csv',row.names = F)
gsea_1 <- gsea@result[gsea@result$pvalue<=0.05,]
sig_gmt <- gmts[gmts$term %in% unique(gsea_1$ID),]
sig_gmt <- sig_gmt[sig_gmt$gene %in% unique(degs$genes),]
sig_gmt_1 <- sig_gmt %>%
  group_by(term) %>%
  summarise(total_genes=paste(gene,collapse = '/'))

gsea_1 <- merge(gsea_1,sig_gmt_1,by.x = "ID",by.y = "term")
rownames(gsea_1) <- gsea_1[,1]
write.csv(gsea@result,'VEC_C1_GSEA_Msigdb_sig.csv',row.names = F)
gsea_1 <- gsea_1[-11]
colnames(gsea_1)[11] <- "core_enrichment"
gsea@result <- gsea_1

#write.csv(gsea_1,'VEC_C1_GSEA_Msigdb_sig.csv',row.names = F)
keggs <- dplyr::filter(gsea@result, grepl('/APLN/', core_enrichment))
write.csv(keggs,'VEC_C1_GSEA_APLN.csv')
save(gsea,file = 'VEC_C1_GSEA.rda')
gseaNb(object = gsea,
       geneSetID = "SANA_TNF_SIGNALING_DN",
       newGsea = T,
       addPoint = F,
       addPval = T,
       arrowType = 'open',
       geneCol = 'black',
       pCol = 'black',
       newHtCol = c("blue","white", "red"),
       addGene = c("APLN","APLNR"))

lapply(terms, function(x){
  gseaNb(object = gseaRes,
         geneSetID = x,
         addPval = T,
         pvalX = 0.75,pvalY = 0.75,
         pCol = 'black',
         pHjust = 0,
         subPlot = 2)
}) -> plotsall
cowplot::plot_grid(plotlist = gseaList1,ncol = 2,align = 'hv')

######################################## CellTypeDendrograms
library(ggtree)

obj <- readRDS('s8_VEC_recluster.rds')
#============= by sample
Idents(obj) <- "orig.ident"
celltypes <- unique(obj$subtype)
groups <- unique(obj@meta.data[,c("orig.ident","group")])
# write.csv(groups,'group.csv',row.names = F)
for (c in celltypes) {
  pbmc <- subset(obj,subtype == c)
  avgexp <- as.data.frame(AverageExpression(pbmc, verbose = FALSE)$RNA)
  avgexp <- avgexp[which(rowSums(avgexp) > 0.1),]
  #avgexp <- as.data.frame(t(avgexp))
  avgexp <- cor(avgexp, method = 'spearman')
  # outprex <- paste(c,'geneavg.csv',sep = '_')
  # write.csv(avgexp,outprex,row.names = F,quote = F)
  dists <- dist(avgexp,method = "euclidean") 
  
  hc <- hclust(dists, method = "ward.D2") #"ave"
  #dend1 <- as.dendrogram(hc)
  #### version1
  #labelColors = c("#CDB380", "#036564", "#EB6841", "#EDC951")
  #clusMember = cutree(hc, length(unique(groups$group2)))
  #plot(dend1, 
  #type = "rectangle", 
  #ylab="Height",
  #main=c)
  #### version2
  #ggtree(hc) + geom_tiplab() + layout_dendrogram()
  #clus <- cutree(hc,5)
  #g <- split(names(clus),clus)
  p <- ggtree(hc,linetype = 'dashed')
  #clades <- sapply(g, function(n) MRCA(p,n))
  #p <- groupClade(p,clades,group_name = 'subtree') #+ aes(color = subtree)
  rownames(groups) <- NULL
  colnames(groups) <- c("label","groups")
  p1 <- p %<+% groups +
    layout_dendrogram() +
    geom_tippoint(
      size = 3,
      shape = 21,
      aes(fill = groups), #,x = x + 0.5
      color = 'black'
    ) +
    #geom_tiplab(aes(label = groups),size = 3,hjust = .5,color = 'black') +
    geom_tiplab(angle = 25,hjust = 1,vjust = 1,show.legend = F,offset = -0.01,size=3) +
    ggtitle(c) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = as.character(paletteer_d("ggthemes::Superfishel_Stone"))) + 
    #scale_color_brewer(palette = "Set2",breaks = 1:4) +
    #theme_dendrogram(plot.margin=margin(6,6,80,6)) +
    theme(legend.position = c(.9,.7))
  ggsave(paste(c,'dendrogram.pdf',sep='_'), plot = p1, width = 6, height = 4)
}

#============= by celltype
Idents(obj) <- "subtype"
gs <- unique(obj$group)
plots <- list()
for(g in gs){
  pbmc <- subset(obj,group == g)
  avgexp <- as.data.frame(AverageExpression(pbmc, verbose = FALSE)$RNA)
  avgexp <- avgexp[which(rowSums(avgexp) > 0.1),]
  avgexp <- cor(avgexp, method = 'spearman')
  dists <- dist(avgexp,method = "euclidean") 
  hc <- hclust(dists, method = "ward.D2")
  p <- ggtree(hc,linetype = 'dashed')
  p1 <- p +
    layout_dendrogram() +
    geom_tippoint(
      size = 3,
      shape = 21,
      color = 'black') +
    geom_tiplab(angle = 90,hjust = 1,vjust = 1,show.legend = F,offset = -0.01,size=5) + 
    ggtitle(g) +
    theme(plot.title = element_text(hjust = 0.5)) 
  plots[[g]] <- p1
}
pp <- plots[[1]] + plots[[3]] + plots[[2]]
ggsave('dendrogram_subtype.pdf', plot = pp, width = 11.2, height = 3.5)


