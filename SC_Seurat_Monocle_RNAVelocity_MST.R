# Chen-Yu Liu cymaxliu@gmail.com 
# Feb  2022

# This script  is for the analysis of single cell RNAseq  with Seurat4 , Monocle2, using  a cell ranger output as well as alignment bam files.
library(Seurat)
library(dplyr)
library(Matrix)
library(Signac)
library(SeuratWrappers)
#library(monocle3)
library(monocle)
library(ggplot2)
library(patchwork)
library(velocyto.R)
library(clusterProfiler)
library(DOSE)
library("org.Mm.eg.db")
library("DEsingle")
#library(EnsDb.Mmusculus.v79)
#Original and revise from -filter_monocle2-1-22-fc0.5R0.05pca10_revise_1_25.R 
set.seed(1234)
setting="ProjectName_" #setup project name therefore all output will come with this prefix

setwd("working_folder_path") #define the working folder
KO <-Read10X(data.dir = "Cell_Ranger_output_folder/KO_out/outs/filtered_feature_bc_matrix") # load  cell ranger output and setup Seurat object  
KO<- CreateSeuratObject(counts  = KO, min.cells = 3, min.genes = 200, project = "KO")
WT <-Read10X(data.dir = "Cell_Ranger_output_folder/WT_out/outs/filtered_feature_bc_matrix")
WT<- CreateSeuratObject(counts  = WT, min.cells = 3, min.genes = 200, project = "WT")
#setup gene names for downstream differential expres sion analysis
genes <- read.table("Cell_Ranger_output_folder/KO_out/outs/filtered_feature_bc_matrix/features.tsv", stringsAsFactors = FALSE) 


mito.genes <- grep(pattern = "^mt-", x = rownames(x = KO), value = TRUE)
percent.mito <- Matrix::colSums(KO@assays$RNA[mito.genes, ])/Matrix::colSums(KO@assays$RNA)

# AddMetaData adds mito percentage columns to KO object
KO <- AddMetaData(object = KO, metadata = percent.mito, col.name = "percent.mito")

png(paste(setting,"violin_KO.png",sep="_"))
VlnPlot(object = KO, features= c("nCount_RNA", "nFeature_RNA", "percent.mito"))
dev.off()

png(paste(setting,"feature_scatter_KO.png",sep="_"))
par(mfrow = c(1, 2))
FeatureScatter(object = KO, feature1 = "nFeature_RNA", feature2 = "percent.mito")
FeatureScatter(object = KO, feature1 = "nFeature_RNA", feature2 = "nFeature_RNA")
dev.off()
# filtered out reads from KO object with criteria 
KO=subset(KO, subset = nFeature_RNA > 250 & nCount_RNA > 5000 & nCount_RNA <30000 & percent.mito < 15)


mito.genes <- grep(pattern = "^mt-", x = rownames(x = WT), value = TRUE)
percent.mito <- Matrix::colSums(WT@assays$RNA[mito.genes, ])/Matrix::colSums(WT@assays$RNA)
# AddMetaData adds mito percentage columns to WT object
WT <- AddMetaData(object = WT, metadata = percent.mito, col.name = "percent.mito")

png(paste(setting,"violin_WT.png",sep="_"))
VlnPlot(object = WT, features= c("nCount_RNA", "nFeature_RNA", "percent.mito"))
dev.off()

png(paste(setting,"feature_scatter_WT.png",sep="_"))
FeatureScatter(object = WT, feature1 = "nFeature_RNA", feature2 = "percent.mito")
FeatureScatter(object = WT, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")
dev.off()
# filtered out reads from WTobject with criteria 
WT=subset(WT, subset = nFeature_RNA >250  & nCount_RNA > 5000 & nCount_RNA <30000 & percent.mito < 15)

#Merge KO and WT objects
seurat_all <- merge(WT ,y = KO, add.cell.ids = c("WT", "KO"), project = "Some_Project")
seurat_all <- NormalizeData(seurat_all) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
seurat_all <-RunTSNE(seurat_all,reduction = "pca",verbose = FALSE)
seurat_all <- RunUMAP(seurat_all, reduction = "pca", dims = 1:10)

#check the number of cells after removing unwanted cells
KO_cell_num=length(colnames(KO))
WT_cell_num=length(colnames(WT))
 
#TSNE prior to Batch Effect Removal 
png(paste(setting,'Seurat_cluster_no_harmony_tsne.png',sep="_"))
DimPlot(seurat_all, group.by = c("orig.ident"),reduction = 'tsne')
dev.off()

##Batch Effect Removal
#seurat_all <- seurat_alld %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
#seurat_alld <- RunHarmony(seurat_all,group.by.vars="orig.ident",reduction = 'pca',assay = "RNA")
#seurat_alld <- RunUMAP(seurat_alld, reduction = "harmony", dims = 1:30)
#seurat_alld <-RunTSNE(seurat_alld,reduction = "pca",verbose = FALSE)
#png('Seurat_cluster_post_harmony_tsne.png')
#DimPlot(seurat_all, group.by = c("seurat_clusters"),reduction = 'tsne')
#dev.off()

seurat_alld <- FindNeighbors(seurat_all, reduction = "pca", dims = 1:10) %>% FindClusters(reduction.type = "umap",resolution = 0.09,) # define PCA dimension and clustering resolution

png(paste(setting,'Seurat_cluster_umap.png',sep="_"))
DimPlot(seurat_alld, group.by = c("seurat_clusters"),label = TRUE,)
dev.off()


seurat_alld.markers <- FindAllMarkers(seurat_alld, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5) # this criteria were used to identify the cell markers
save.image(file='Markers.RData')


#load('imonocle2-5000cut_low_resolutionRNAVelocity_fc0.5R0.09tsnepca10.RData')

# Define the marker genes
G_cell= c('Ddx4', 'Dazl', 'Sycp1', 'Uchl1', 'Crabp1', 'Stra8', 'Sohlh1')
S_cell= c('Sox9', 'Amhr2', 'Gstm6', 'Amh')
M_cell=  c('Acta2', 'Myh11', 'Fbxl22', 'Tagln')
E_cell=  c('Ly6c1', 'Pecam1', 'Vwf', 'Tek', 'Egfl7')
L_cell=  c('Hsd3b1', 'Cyp11a1', 'Cyp17a1', 'Star')
St_cell=  c('Igf1', 'Gng11', 'Inhba', 'Col1a1', 'Dlk1')
Ma_cell=  c('Adgre1', 'Csf1r', 'F13a1', 'Mrc1')
I_cell=  c('Cnmd','Igkc')

# Plot the features marker genei for each cell type  with default UMAP
png(paste(setting,'G_cell.png',sep="_"))
FeaturePlot(seurat_alld, features = c('Ddx4', 'Dazl', 'Sycp1', 'Uchl1', 'Crabp1', 'Stra8', 'Sohlh1'), min.cutoff = "q9")
dev.off()

png(paste(setting,'S_cell.png',sep="_"))
FeaturePlot(seurat_alld, features = c('Sox9', 'Amhr2', 'Gstm6', 'Amh'), min.cutoff = "q9")
dev.off()

png(paste(setting,'M_cell.png',sep="_"))
FeaturePlot(seurat_alld, features = c('Acta2', 'Myh11', 'Fbxl22', 'Tagln'), min.cutoff = "q9")
dev.off()

png(paste(setting,'E_cell.png',sep="_"))
FeaturePlot(seurat_alld, features = c('Ly6c1', 'Pecam1', 'Vwf', 'Tek', 'Egfl7'), min.cutoff = "q9")
dev.off()

png(paste(setting,'L_cell.png',sep="_"))
FeaturePlot(seurat_alld, features = c('Hsd3b1', 'Cyp11a1', 'Cyp17a1', 'Star'), min.cutoff = "q9")
dev.off()

png(paste(setting,'St_cell.png',sep="_"))
FeaturePlot(seurat_alld, features = c('Igf1', 'Gng11', 'Inhba', 'Col1a1', 'Dlk1'), min.cutoff = "q9")
dev.off()

png(paste(setting,'Ma_cell.png',sep="_"))
FeaturePlot(seurat_alld, features = c('Adgre1', 'Csf1r', 'F13a1', 'Mrc1'), min.cutoff = "q9")
dev.off()

png(paste(setting,'I_cell.png',sep="_"))
FeaturePlot(seurat_alld, features = c('Cnmd'), min.cutoff = "q9")
dev.off()

# Assign cell type for each cluster by the cell type marker genes listed above
# Usually a cluster will have marker gene for multiple cell types , only the cell type which has the largest number of marker genes will be assign to that cluster.
# germ_list is a list of clusters corresponding to G_ Cell  while "cluster" is a list for various cell types 
cluster=c()
germ_list=c()
for (i in 0:(length(levels(seurat_alld$seurat_clusters))-1))
{
   c1=seurat_alld.markers[which(seurat_alld.markers$cluster==i),]

   G_list=c1$gene
   
   name=""
   max=0
   if(length(which(G_list %in% G_cell))/length(G_cell)>max)
   {   max=length(which(G_list %in% G_cell))/length(G_cell)
    name="G_cell"
    germ_list=cbind(germ_list,i)
    print(i)
    print( "G_ Cell")
    print(G_list[which(G_list %in% G_cell)])
   }
   if(length(which(G_list %in% S_cell ))/length(S_cell)>max)
   { 
    max=length(which(G_list %in% S_cell ))/length(S_cell) 
    name="S_cell"
    print(i)
    print( "S_cell")
    print(G_list[which(G_list %in% S_cell)])
   }
   if(length(which(G_list %in% M_cell ))/length(M_cell)>max)
   {
    max=length(which(G_list %in% M_cell ))/length(M_cell)
    name="M_cell"
    print(i)
    print( "M_cell")
    print(G_list[which(G_list %in% M_cell)])
   }
   if(length(which(G_list %in% E_cell))/length(E_cell)>max)
   {
    max=length(which(G_list %in% E_cell))/length(E_cell)
    name="E_cell"
    print(i)
    print( "E_cell")
    print(G_list[which(G_list %in% E_cell)])
   }
   if(length(which(G_list %in% L_cell))/length(L_cell)>max)
   {
    max=length(which(G_list %in% L_cell))/length(L_cell)
    name="L_cell"
    print(i)
    print( "L_cell")
    print(G_list[which(G_list %in% L_cell)])
   }
   if(length(which(G_list %in% St_cell))/length(St_cell)>max)
   {
    max=length(which(G_list %in% St_cell))/length(St_cell)
    name="St_cell"
    print(i)
    print( "St_cell")
    print(G_list[which(G_list %in% St_cell)])
   }
   if(length(which(G_list%in% Ma_cell))/length(Ma_cell)>max)
   {
    max=length(which(G_list%in% Ma_cell))/length(Ma_cell)
    name="Ma_cell"
    print(i)
    print( "Ma_cell")
    print(G_list[which(G_list %in% Ma_cell)])
   }
   if(length(which(G_list %in% I_cell))/length(I_cell)>max)
   {
    max=length(which(G_list %in% I_cell))/length(I_cell)
    name="I_cell" 
    print(i)
    print( "I_cell")
    print(G_list[which(G_list %in% I_cell)])
   }
   cluster=cbind(cluster,name)
}


 # Rename the clusters based on the new cluster names
 current.cluster.ids=levels(seurat_alld$seurat_clusters)
 new.cluster.ids=as.vector(cluster)
 seurat_alld@meta.data$seurat_clusters <- plyr::mapvalues(x = seurat_alld@meta.data$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
  
 
cell_type_marker= c('Ddx4', 'Dazl', 'Sycp1', 'Uchl1', 'Crabp1', 'Stra8', 'Sohlh1','Sox9', 'Amhr2', 'Gstm6', 'Amh','Acta2', 'Myh11', 'Fbxl22', 'Tagln','Ly6c1', 'Pecam1', 'Vwf', 'Tek', 'Egfl7','Hsd3b1', 'Cyp11a1', 'Cyp17a1', 'Star','Igf1', 'Gng11', 'Inhba', 'Col1a1', 'Dlk1','Adgre1', 'Csf1r', 'F13a1', 'Mrc1')

 
 # UMAP for rename clusters
 png(paste(setting,'seurat_cluster_rename.png',sep="_"))
 DimPlot(seurat_alld, label = TRUE, group.by='seurat_clusters')
 dev.off()
 
 # Retrieve  the cells from germ cell clusters 
 G_=subset(x = seurat_alld, idents = (as.vector(as.character(germ_list))))
 G_.combined <- FindNeighbors(G_, reduction = "pca", dims = 1:30) %>% FindClusters(reduction.type = "umap",resolution = 0.05) #G_ cell as a whole
 
 # rename cells' identities
for(i in levels(seurat_alld)) {
  cells_to_reid <- WhichCells(seurat_alld, idents = i)
  newid <- names(sort(table(seurat_alld$seurat_clusters[cells_to_reid]),decreasing=TRUE))[1]
  Idents(seurat_alld, cells = cells_to_reid) <- newid
}

png(paste(setting,'dotplot_seurat_cluster_rename.png',sep="_"),width = 1800, height = 800, units = "px")
seurat_alld[['groups']] <- sample(x = c('WT', 'KO'), size = ncol(x = seurat_alld), replace = TRUE)
DotPlot(object = seurat_alld, features = cell_type_marker, split.by = 'groups',cols = c("blue", "red"), dot.scale = 30)+labs(y="Class", x = "Genes")+theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "gray",linetype = "dashed", size=0.35), panel.border = element_rect(colour = "black", fill=NA, size=2))+ facet_grid(~ unique(seurat_alld@meta.data$orig.ident))
dev.off()

# Print out various cell types' stat survey
cell_stat=c("WT","KO")
i=1
for ( icluster in levels(factor(cluster)) ) 
{

 for (condition in levels(factor(c("WT","KO"))))
 {
   item=subset(x = seurat_alld, orig.ident ==condition & seurat_clusters==icluster)
   if(condition=="WT")
   {
    wt=(length(colnames(item)))
   }
   else
   {
    ko=(length(colnames(item))) 
   }   
 }
 cell_stat=rbind(cell_stat,c(wt,ko)) 
  i=i+1
 rownames(cell_stat)[i]=icluster
}
write.csv(cell_stat, paste(setting,"cell_cluster_stat.csv",sep="_"))

germ_cell_stat=c("WT","KO")
i=1
for ( gcluster in levels(G_.combined)  )
{

 for (condition in levels(factor(c("WT","KO"))))
 {
   item=subset(x = G_.combined, orig.ident ==condition & seurat_clusters==gcluster)
   if(condition=="WT")
   {
    wt=(length(colnames(item)))
   }
   else
   {
    ko=(length(colnames(item)))
   }
 }
 germ_cell_stat=rbind(germ_cell_stat,c(wt,ko))
  i=i+1
 rownames(germ_cell_stat)[i]=gcluster
}
write.csv(germ_cell_stat, paste(setting,"cell_cluster_stat.csv",sep="_"))

#########################
####RNA Velocity Analysis 
########################
ldat <- read.loom.matrices("Cell_Ranger_output_folder/WT_out/velocyto/WT_out.loom")

##Rename cells for matching RNA Velocity Analysis
 
G_.combined=RenameCells(G_.combined, new.names = gsub("_", "_out:", colnames(G_.combined)))
G_.combined=RenameCells(G_.combined, new.names = gsub("-1", "x", colnames(G_.combined)))

##subset the cells into WT and KO groups, and output their embeddings/cluster/cell IDs
K=subset(x=G_.combined, subset = orig.ident == "KO") #extract cells from KO 
#write.csv(Embeddings(K, reduction = "umap"), file = "KO_cell_embeddings.csv")
#write.csv(Cells(K), file = "KO_cellID_obs.csv", row.names = FALSE)
#write.csv(K@meta.data$seurat_clusters, file ="KO_clusters.csv")

W=subset(x=G_.combined, subset = orig.ident == "WT") #extract cells from WT
#write.csv(Embeddings(W, reduction = "umap"), file = "WT_cell_embeddings.csv")
#write.csv(Cells(W), file = "WT_cellID_obs.csv", row.names = FALSE)
#write.csv(K@meta.data$seurat_clusters, file ="WT_clusters.csv")


emat <- ldat$spliced; 
nmat <- ldat$unspliced;

emat <- emat[,colnames(W)]
nmat <- nmat[,colnames(W)]

cluster.label <- W$seurat_clusters

current.cluster.ids=levels(cluster.label)
new.cluster.ids=as.vector(c("#b7ff98", "#ffc3e0", "#d7c5ff"))

cell.colors<- plyr::mapvalues(x =cluster.label, from = current.cluster.ids, to = new.cluster.ids) # change the cell label to colors



emb=Embeddings(object = G_.combined, reduction = "umap")# take umap embeddings
cell.dist <- as.dist(1-armaCor(t(emb)))
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)

fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile)

png(paste(setting,"WT_RNAVelocity.png",sep="_"))
show.velocity.on.embedding.cor(emb,rvel.cd,n=300,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.5),cex=0.8,arrow.scale=5,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
dev.off()

# repeat procedures above agin but for KO
ldat <- read.loom.matrices("Cell_Ranger_output_folder/KO_out/velocyto/KO_out.loom")

emat <- ldat$spliced;
nmat <- ldat$unspliced;

emat <- emat[,colnames(K)]
nmat <- nmat[,colnames(K)]

cluster.label <- K$seurat_clusters
current.cluster.ids=levels(cluster.label)
new.cluster.ids=as.vector(c("#b7ff98", "#ffc3e0", "#d7c5ff"))
cell.colors<- plyr::mapvalues(x =cluster.label, from = current.cluster.ids, to = new.cluster.ids)

emb=Embeddings(object = G_.combined, reduction = "umap")
cell.dist <- as.dist(1-armaCor(t(emb)))
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)

fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile)

png(paste(setting,"KO_RNAVelocity.png",sep="_"))
show.velocity.on.embedding.cor(emb,rvel.cd,n=300,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.5),cex=0.8,arrow.scale=5,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
dev.off()






######

png(paste(setting,'seurat_G_.png',sep="_")) # make plot of G_ cells
DimPlot(G_.combined , label = TRUE, group.by='seurat_clusters')
dev.off()


G_.combined <- RunUMAP(
  G_.combined,
  reduction = "pca",
  dims = 1:10,
  reduction.name = "UMAP"
)

G_.combined <- RunTSNE(
  G_.combined,
  reduction = "pca",
  dims = 1:10,
  reduction.name = "TSNE"
)

 current.cluster.ids=c("0","1","2")
 new.cluster.ids=c("SubclusterG2","SubclusterG1","SubclusterG3")
 G_.combined@meta.data$seurat_clusters <- plyr::mapvalues(x = G_.combined@meta.data$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids) # assign the new G_ cell clusters


# change the cell identity of G_ cells
#for(i in levels(G_.combined)) {
#  cells_to_reid <- WhichCells(G_.combined, idents = i)
#  newid <- names(sort(table(G_.combined$seurat_clusters[cells_to_reid]),decreasing=TRUE))[1]
#  Idents(G_.combined, cells = cells_to_reid) <- newid
#}
#G_.combined$assigned_celltype <- Idents(G_.combined)


#Plot G_ Cell with Markers

png(paste(setting,"cell_type.png",sep="_"))
DimPlot(G_.combined, label = TRUE)
dev.off()


png(paste(setting,'SubclusterG1_umap.png',sep="_"))
FeaturePlot(G_.combined, features = c('Id4', 'Etv5'), min.cutoff = 0,reduction= "umap")
dev.off()

png(paste(setting,'SubclusterG2_umap.png',sep="_"))
FeaturePlot(G_.combined, features = c('Gfra1', 'Sdc4'), min.cutoff = 0,reduction= "umap")
dev.off()

png(paste(setting,'SubclusterG3_umap.png',sep="_"))
FeaturePlot(G_.combined, features = c('Kit', 'Dmrtb1'), min.cutoff = 0,reduction= "umap")
dev.off()

png(paste(setting,'SubclusterG4_umap.png',sep="_"))
FeaturePlot(G_.combined, features = c('Prss50','Ly6k'), min.cutoff = 0,reduction= "umap")
dev.off()

png(paste(setting,'others_umap.png',sep="_"))
FeaturePlot(G_.combined, features = c('Chd1','Stra8','Ddx4','Pramef12'), min.cutoff = 0,reduction= "umap")
dev.off()

SubclusterG1_marker=c('Id4', 'Etv5', 'Gfra1', 'Lhx1', 'Ret')
SubclusterG2_marker=c('Ddit4', 'Utf1', 'Egr4', 'Neurog3', 'Sox3', 'Rarg', 'Nanos3')
SubclusterG3_marker=c('Kit', 'Stra8', 'Dnmt3b', 'Rhox10')
SubclusterG4_marker=c('Dmrtb1', 'Prss50', 'Ly6k', 'Prdm9')

# plot G_ Cell with Marker with TSNE 
png(paste(setting,'SubclusterG1_tsne.png',sep="_"))
FeaturePlot(G_.combined, features = c('Id4', 'Etv5'), min.cutoff = 0,reduction= "tsne")
dev.off()

png(paste(setting,'SubclusterG2_tsne.png',sep="_"))
FeaturePlot(G_.combined, features = c('Gfra1', 'Sdc4'), min.cutoff = 0,reduction= "tsne")
dev.off()        
                 
png(paste(setting,'SubclusterG3_tsne.png',sep="_"))
FeaturePlot(G_.combined, features = c('Kit', 'Dmrtb1'), min.cutoff = 0,reduction= "tsne")
dev.off()        
                 
png(paste(setting,'SubclusterG4_tsne.png',sep="_"))
FeaturePlot(G_.combined, features = c('Prss50','Ly6k'), min.cutoff = 0,reduction= "tsne")
dev.off()

png(paste(setting,'others_tsne.png',sep="_"))
FeaturePlot(G_.combined, features = c('Chd1','Stra8','Ddx4','Pramef12'), min.cutoff = 0,reduction= "tsne")
dev.off()

# plot WT and KO of G_ Cells 

G__WT=subset(x = G_.combined, subset = orig.ident== "WT")
png('G__WT_cluster_umap.png')
DimPlot(G__WT,label = TRUE,group.by = 'seurat_clusters', cols = c('SubclusterG1' = "#b7ff98", 'SubclusterG2' = "#ffc3e0", 'SubclusterG3' = "#d7c5ff"))
dev.off()

G__KO=subset(x = G_.combined, subset = orig.ident== "KO")
png('G__KO_cluster_umap.png')
DimPlot(G__KO,label = TRUE,group.by = 'seurat_clusters', cols = c('SubclusterG1' = "#b7ff98", 'SubclusterG2' = "#ffc3e0", 'SubclusterG3' = "#d7c5ff"))
dev.off()



# Construct CellDataSet object
pd <- data.frame(cell_id = colnames(G_.combined),
                 sample=G_.combined@ meta.data$orig.ident,
                 cell_type = G_.combined@ meta.data$seurat_clusters,
                 row.names = colnames(G_.combined))

pd <- new("AnnotatedDataFrame", data = pd)

fd <- data.frame(gene_id = rownames(G_.combined),
                 gene_short_name = rownames(G_.combined),
                 row.names = rownames(G_.combined))

fd <- new("AnnotatedDataFrame", data = fd)

GERMS <- newCellDataSet(G_.combined@assays$RNA@counts, phenoData = pd, featureData = fd,expressionFamily=negbinomial.size()) # create 

GERMS <- detectGenes(GERMS, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(GERMS), num_cells_expressed >= 50))

GERMS<- estimateSizeFactors(GERMS)
GERMS <- estimateDispersions(GERMS)
diff_test_res <- differentialGeneTest(GERMS[expressed_genes,], fullModelFormulaStr="~cell_type", cores=24) # Determine Differential Genes

ordering_genes <- row.names (subset(diff_test_res, qval < 0.1))


GERMS <- setOrderingFilter(GERMS, ordering_genes)
GERMS <- reduceDimension(GERMS, use_irlba=T) 
GERMS <- orderCells(GERMS, num_paths=1, reverse=T)

pdf(paste(setting,"monocle2_plot2_irlab.pdf",sep="_"))
plot_cell_trajectory(GERMS, color_by = "cell_type", cell_size = 1) + scale_color_manual(breaks = c("SubclusterG1", "SubclusterG2", "SubclusterG3"), values=c("#b7ff98", "#ffc3e0", "#d7c5ff"))
dev.off()


pdf(paste(setting,"monocle2_plot3_irlab.pdf",sep="_"))
plot_cell_trajectory(GERMS, color_by = "sample", cell_size = 1) 
dev.off()

SubclusterG2_marker=c('Ddit4', 'Utf1', 'Egr4', 'Neurog3', 'Sox3', 'Rarg', 'Nanos3')
SubclusterG3_marker=c('Kit', 'Stra8', 'Dnmt3b', 'Rhox10')
SubclusterG4_marker=c('Dmrtb1', 'Prss50', 'Ly6k', 'Prdm9')

pdf(paste(setting,"SubclusterG1_monocle2_plot3_irlab.pdf",sep="_"))
plot_cell_trajectory(GERMS, markers =c('Id4', 'Etv5', 'Gfra1', 'Lhx1', 'Ret'),use_color_gradient = TRUE)
dev.off()


pdf(paste(setting,"SubclusterG2_monocle2_plot3_irlab.pdf",sep="_"))
plot_cell_trajectory(GERMS, markers =c('Ddit4', 'Utf1', 'Egr4', 'Neurog3', 'Sox3', 'Rarg', 'Nanos3'),use_color_gradient = TRUE)
dev.off()

pdf(paste(setting,"SubclusterG3_monocle2_plot3_irlab.pdf",sep="_"))
plot_cell_trajectory(GERMS, markers =c('Kit', 'Stra8', 'Dnmt3b', 'Rhox10'),use_color_gradient = TRUE)
dev.off()

pdf(paste(setting,"SubclusterG4_monocle2_plot3_irlab.pdf",sep="_"))
plot_cell_trajectory(GERMS, markers =c('Dmrtb1', 'Prss50', 'Ly6k', 'Prdm9'),use_color_gradient = TRUE)
dev.off()

sig_genes <- subset(diff_test_res, qval < 1e-100)
sig_genes[,c("gene_short_name", "pval", "qval")]
sig_gene_names <- row.names(subset(diff_test_res, qval <1e-100 ))

png(paste(setting,"pseudotime_heatmap.png",sep="_"))
plot_pseudotime_heatmap(GERMS[sig_gene_names,], num_clusters = 3, cores = 1, show_rownames = T)
dev.off()

pdf(paste(setting,"monocle2_SubclusterG1.pdf",sep="_"))
GERMS_expressed_genes <-  row.names(subset(fData(GERMS),
num_cells_expressed >= 10))
GERMS_filtered <- GERMS[GERMS_expressed_genes,]
my_genes <- row.names(subset(fData(GERMS_filtered), gene_short_name %in% c('Id4', 'Etv5','Gfra1','Lhx1', 'Ret' )))
cds_subset <- GERMS_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "cell_type")
dev.off()

pdf(paste(setting,"monocle2_SubclusterG2.pdf",sep="_"))
GERMS_expressed_genes <-  row.names(subset(fData(GERMS),
num_cells_expressed >= 10))
GERMS_filtered <- GERMS[GERMS_expressed_genes,]
my_genes <- row.names(subset(fData(GERMS_filtered), gene_short_name %in% c( 'Ddit4', 'Utf1', 'Egr4', 'Neurog3', 'Sox3', 'Rarg', 'Nanos3'  )))
cds_subset <- GERMS_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "cell_type")
dev.off()


pdf(paste(setting,"monocle2_SubclusterG3.pdf",sep="_"))
GERMS_expressed_genes <-  row.names(subset(fData(GERMS),
num_cells_expressed >= 10))
GERMS_filtered <- GERMS[GERMS_expressed_genes,]
my_genes <- row.names(subset(fData(GERMS_filtered), gene_short_name %in% c('Kit','Dnmt3b', 'Rhox10','Stra8' )))
cds_subset <- GERMS_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "cell_type")
dev.off()

pdf(paste(setting,"monocle2_SubclusterG4.pdf",sep="_"))
GERMS_expressed_genes <-  row.names(subset(fData(GERMS),
num_cells_expressed >= 10))
GERMS_filtered <- GERMS[GERMS_expressed_genes,]
my_genes <- row.names(subset(fData(GERMS_filtered), gene_short_name %in% c('Dmrtb1', 'Prss50', 'Ly6k', 'Prdm9' )))
cds_subset <- GERMS_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "cell_type")
dev.off()


######################
# Go Term Enrichment 
######################

get.frac <- function(x){
    y <- as.numeric(strsplit(x, '/')[[1]])
    return(y[[1]] / y[[2]])
}
summarize.enrich <- function(ego, name, label){
lst <- list()
    d <- as.data.frame(ego)
    if (nrow(d) > 0){
      d$ontology <- name
      d$label <- label
      d$frac <- sapply(d$GeneRatio, get.frac)
    }
    lst[[name]] <- d
  df <- do.call(rbind, lst)
  return(df)
}

summarize.kegg <- function(ekg, labels){
  d <- as.data.frame(ekg)
  if (nrow(d) > 0){
    d$ontology <- 'kegg'
    for (label in names(labels)){
      d[label] <- labels[[label]]
    }
    d$frac <- sapply(d$GeneRatio, get.frac)
  }
  return(d)
}


sig_genes <- subset(diff_test_res, qval < 1e-100)
sig_genes[,c("gene_short_name", "pval", "qval")]
sig_gene_names <- row.names(subset(diff_test_res, qval <1e-100 ))

genename=sig_gene_names
ensembl_gene_ids=genes[which( genes[,2] %in% genename),1]
eg = bitr(toupper(ensembl_gene_ids), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

enrich.list <- list()

    # GO enrichment
    go.label <- paste(name,  'go', sep='.')
    sg <- summarize.enrich(enrichGO(eg[,2], OrgDb='org.Mm.eg.db',ont='BP',pvalueCutoff  = 0.01),'Enrichment','GO')
    if (!is.null(sg)){
      enrich.list[['GO']] <- sg
    }

    # KEGG enrichment
    sk <- summarize.enrich( enrichKEGG(gene= eg[,2], organism='mmu', pvalueCutoff = 0.05),'Enrichment','KEGG')
    if (!is.null(sk)){
      enrich.list[['KEGG']] <- sk
     }

    sd <- enrichDO(gene= eg[,2], ont= "DO",pvalueCutoff = 0.5,pAdjustMethod = "BH",qvalueCutoff  = 1)
    
    m <- do.call(rbind, enrich.list)

    m$phred <- -10 * log10(m$p.adjust)


    write.csv(sd, file=paste(setting, "_Disease_Mono_enrichment.csv",sep="_"))
    write.csv(m,file=paste(setting, "_Mono_enrichment.csv",sep="_"))

######################
#MAST APPROACH
######################
suppressPackageStartupMessages({
    library(ggplot2)
    library(GGally)
    library(GSEABase)
    library(limma)
    library(reshape2)
    library(data.table)
    library(knitr)
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    library(stringr)
    library(NMF)
    library(rsvd)
    library(RColorBrewer)
    library(MAST)
})
# Create a SingleCellExperiment  object from Seurat object 
pd <- data.frame(cell_id = colnames(G_.combined),
                 sample=G_.combined@ meta.data$orig.ident,
                 cell_type = G_.combined@ meta.data$seurat_clusters,
                 row.names = colnames(G_.combined))
fd <- data.frame(gene_id = rownames(G_.combined),
                 gene_short_name = rownames(G_.combined),
                 row.names = rownames(G_.combined))


newassay <- as.matrix((G_.combined@assays$RNA@counts))
sca <- FromMatrix(newassay,pd,fd, check_sanity = FALSE)


cond<-factor(colData(sca)$cell_type)
cond<-relevel(cond,"SubclusterG1")
colData(sca)$condition<-cond
zlmCond <- zlm(~condition , sca)

summaryDt<- summary(zlmCond, doLRT='(Intercept)')
#summaryDt[contrast=='cell_typeSubclusterG2' & component=='H']
Diffsummary=summaryDt$datatable[contrast=='(Intercept)' & component=='D']
Diffsummary=Diffsummary[order(Diffsummary[,4]),]
Diffsummary=Diffsummary[which(Diffsummary[,4]<1e-10)] 
sig_gene_names=Diffsummary$primerid
genename=sig_gene_names
ensembl_gene_ids=genes[which( genes[,2] %in% genename),1]
eg = bitr(toupper(ensembl_gene_ids), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

enrich.list <- list()

    # GO enrichment
    go.label <- paste(name,  'go', sep='.')
    sg <- summarize.enrich(enrichGO(eg[,2], OrgDb='org.Mm.eg.db',ont='BP',pvalueCutoff  = 0.01),'Enrichment','GO')
    if (!is.null(sg)){
      enrich.list[['GO']] <- sg
    }

    # KEGG enrichment
    sk <- summarize.enrich( enrichKEGG(gene= eg[,2], organism='mmu', pvalueCutoff = 0.05),'Enrichment','KEGG')
    if (!is.null(sk)){
      enrich.list[['KEGG']] <- sk
     }
    sd <- enrichDO(gene= eg[,2], ont= "DO",pvalueCutoff = 0.5,pAdjustMethod = "BH",qvalueCutoff  = 1)
  m <- do.call(rbind, enrich.list)

  m$phred <- -10 * log10(m$p.adjust)

 write.csv(sd, file=paste(setting, "_Disease_SubclusterG1_MAST_enrichment.csv",sep="_"))
 write.csv(m,file=paste(setting, "SubclusterG1_MAST_enrichment.csv",sep="_"))



summaryDt<- summary(zlmCond, doLRT='conditionSubclusterG2')
#summaryDt[contrast=='cell_typeSubclusterG2' & component=='H']
Diffsummary=summaryDt$datatable[contrast=='conditionSubclusterG2' & component=='D']
Diffsummary=Diffsummary[order(Diffsummary[,4]),]
Diffsummary=Diffsummary[which(Diffsummary[,4]<1e-10)]  
sig_gene_names=Diffsummary$primerid
genename=sig_gene_names
ensembl_gene_ids=genes[which( genes[,2] %in% genename),1]
eg = bitr(toupper(ensembl_gene_ids), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

enrich.list <- list()

    # GO enrichment
    go.label <- paste(name,  'go', sep='.')
    sg <- summarize.enrich(enrichGO(eg[,2], OrgDb='org.Mm.eg.db',ont='BP',pvalueCutoff  = 0.01),'Enrichment','GO')
    #sg <- summarize.go(enrichGO(eg[,2], OrgDb='org.Mm.eg.db'),'Enrichment','GO')#summarize.go(enrichGO(eg[,2], OrgDb='org.Mm.eg.db') , list( experiment='GO', name='test1'))
    if (!is.null(sg)){
      enrich.list[['GO']] <- sg
    }

    # KEGG enrichment
    sk <- summarize.enrich( enrichKEGG(gene= eg[,2], organism='mmu', pvalueCutoff = 0.05),'Enrichment','KEGG')
    if (!is.null(sk)){
      enrich.list[['KEGG']] <- sk
     }
    sd <- enrichDO(gene= eg[,2], ont= "DO",pvalueCutoff = 0.5,pAdjustMethod = "BH",qvalueCutoff  = 1)
  m <- do.call(rbind, enrich.list)

  m$phred <- -10 * log10(m$p.adjust)
 write.csv(sd, file=paste(setting, "_Disease_SubclusterG2_MAST_enrichment.csv",sep="_"))
 write.csv(m,file=paste(setting, "SubclusterG2_MAST_enrichment.csv",sep="_"))


summaryDt<- summary(zlmCond, doLRT='conditionSubclusterG3')
Diffsummary=summaryDt$datatable[contrast=='conditionSubclusterG2' & component=='D']
Diffsummary=Diffsummary[order(Diffsummary[,4]),]
Diffsummary=Diffsummary[which(Diffsummary[,4]<1e-10)]  
sig_gene_names=Diffsummary$primerid
genename=sig_gene_names
ensembl_gene_ids=genes[which( genes[,2] %in% genename),1]
eg = bitr(toupper(ensembl_gene_ids), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

enrich.list <- list()

    # GO enrichment
    go.label <- paste(name,  'go', sep='.')
    sg <- summarize.enrich(enrichGO(eg[,2], OrgDb='org.Mm.eg.db',ont='BP',pvalueCutoff  = 0.01),'Enrichment','GO')
    #sg <- summarize.go(enrichGO(eg[,2], OrgDb='org.Mm.eg.db'),'Enrichment','GO')#summarize.go(enrichGO(eg[,2], OrgDb='org.Mm.eg.db') , list( experiment='GO', name='test1'))
    if (!is.null(sg)){
      enrich.list[['GO']] <- sg
    }

    # KEGG enrichment
    sk <- summarize.enrich( enrichKEGG(gene= eg[,2], organism='mmu', pvalueCutoff = 0.05),'Enrichment','KEGG')
    if (!is.null(sk)){
      enrich.list[['KEGG']] <- sk
     }
    sd <- enrichDO(gene= eg[,2], ont= "DO",pvalueCutoff = 0.5,pAdjustMethod = "BH",qvalueCutoff  = 1)
  m <- do.call(rbind, enrich.list)

  m$phred <- -10 * log10(m$p.adjust)
 write.csv(sd, file=paste(setting, "_Disease_SubclusterG3_MAST_enrichment.csv",sep="_"))
 write.csv(m,file=paste(setting, "SubclusterG3_MAST_enrichment.csv",sep="_"))


######################
##   DEsingle APPROACH
######################

counts <- G_.combined@assays$RNA@counts 
metadata <- G_.combined@meta.data
metadata$cluster_id <- factor(G_.combined@active.ident)
# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)
results <- DEsingle(counts = counts(sce), group = factor(colData(sce)$orig.ident)) # Perform differential gene expression 



Diffsummary=results[which(results$chi2LR1<1),]
sig_gene_names=rownames(Diffsummary)
genename=sig_gene_names
ensembl_gene_ids=genes[which( genes[,2] %in% genename),1]
eg = bitr(toupper(ensembl_gene_ids), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

enrich.list <- list()

    # GO enrichment
    go.label <- paste(name,  'go', sep='.')
    sg <- summarize.enrich(enrichGO(eg[,2], OrgDb='org.Mm.eg.db',ont='BP',pvalueCutoff  = 0.01),'Enrichment','GO')
    if (!is.null(sg)){
      enrich.list[['GO']] <- sg
    }

    # KEGG enrichment
    sk <- summarize.enrich( enrichKEGG(gene= eg[,2], organism='mmu', pvalueCutoff = 0.05),'Enrichment','KEGG')
    if (!is.null(sk)){
      enrich.list[['KEGG']] <- sk
     }
    sd <- enrichDO(gene= eg[,2], ont= "DO",pvalueCutoff = 0.5,pAdjustMethod = "BH",qvalueCutoff  = 1)
  m <- do.call(rbind, enrich.list)

  m$phred <- -10 * log10(m$p.adjust)
  write.csv(sd, file=paste(setting, "_Disease_DEsingle_chi2LR1_enrichment.csv",sep="_"))
  write.csv(m,file=paste(setting,"DEsingle_chi2LR1_enrichment.csv",sep="_"))

Diffsummary=results[which(results$pvalue_LR2<1e-10),]
sig_gene_names=rownames(Diffsummary)
genename=sig_gene_names
ensembl_gene_ids=genes[which( genes[,2] %in% genename),1]
eg = bitr(toupper(ensembl_gene_ids), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

enrich.list <- list()

    # GO enrichment
    go.label <- paste(name,  'go', sep='.')
    sg <- summarize.enrich(enrichGO(eg[,2], OrgDb='org.Mm.eg.db',ont='BP',pvalueCutoff  = 0.01),'Enrichment','GO')
    if (!is.null(sg)){
      enrich.list[['GO']] <- sg
    }

    # KEGG enrichment
    sk <- summarize.enrich( enrichKEGG(gene= eg[,2], organism='mmu', pvalueCutoff = 0.05),'Enrichment','KEGG')
    if (!is.null(sk)){
      enrich.list[['KEGG']] <- sk
     }
    sd <- enrichDO(gene= eg[,2], ont= "DO",pvalueCutoff = 0.5,pAdjustMethod = "BH",qvalueCutoff  = 1)
  m <- do.call(rbind, enrich.list)

  m$phred <- -10 * log10(m$p.adjust)
  write.csv(sd, file=paste(setting, "_Disease_DEsingle_pvalue_LR2_enrichment.csv",sep="_"))
  write.csv(m,file=paste(setting,"DEsingle_pvalue_LR2_enrichment.csv",sep="_"))



Diffsummary=results[which(results$pvalue_LR3<1e-10),]
sig_gene_names=rownames(Diffsummary)
genename=sig_gene_names
ensembl_gene_ids=genes[which( genes[,2] %in% genename),1]
eg = bitr(toupper(ensembl_gene_ids), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

enrich.list <- list()

    # GO enrichment
    go.label <- paste(name,  'go', sep='.')
    sg <- summarize.enrich(enrichGO(eg[,2], OrgDb='org.Mm.eg.db',ont='BP',pvalueCutoff  = 0.01),'Enrichment','GO')
    if (!is.null(sg)){
      enrich.list[['GO']] <- sg
    }

    # KEGG enrichment
    sk <- summarize.enrich( enrichKEGG(gene= eg[,2], organism='mmu', pvalueCutoff = 0.05),'Enrichment','KEGG')
    if (!is.null(sk)){
      enrich.list[['KEGG']] <- sk
     }
    sd <- enrichDO(gene= eg[,2], ont= "DO",pvalueCutoff = 0.5,pAdjustMethod = "BH",qvalueCutoff  = 1)
  m <- do.call(rbind, enrich.list)

  m$phred <- -10 * log10(m$p.adjust)
  write.csv(sd, file=paste(setting, "_Disease_DEsingle_pvalue_LR3_enrichment.csv",sep="_"))
  write.csv(m,file=paste(setting,"DEsingle_pvalue_LR3_enrichment.csv",sep="_"))


