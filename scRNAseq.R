library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)

#Load datasets 

ds001.data <- Read10X(data.dir="O:/DC001")
ds002.data <- Read10X(data.dir="O:/DC002")
ds007.data <- Read10X(data.dir="O:/DC007")
ds008.data <- Read10X(data.dir="O:/DC008")
ds009.data <- Read10X(data.dir="O:/DC009")
ds012.data <- Read10X(data.dir="O:/DC012")
ds013.data <- Read10X(data.dir="O:/DC013")
ds014.data <- Read10X(data.dir="O:/DC014")

#Create Seurat objects

ds001 <- CreateSeuratObject(counts = ds001.data, project = "CFA", min.cells = 3, min.features = 200)
ds002 <- CreateSeuratObject(counts = ds002.data, project = "asyn", min.cells = 3, min.features = 200)
ds007 <- CreateSeuratObject(counts = ds007.data, project = "PBS", min.cells = 3, min.features = 200)
ds008 <- CreateSeuratObject(counts = ds008.data, project = "CFA", min.cells = 3, min.features = 200)
ds009 <- CreateSeuratObject(counts = ds009.data, project = "asyn", min.cells = 3, min.features = 200)
ds012 <- CreateSeuratObject(counts = ds012.data, project = "PBS", min.cells = 3, min.features = 200)
ds013 <- CreateSeuratObject(counts = ds013.data, project = "CFA", min.cells = 3, min.features = 200)
ds014 <- CreateSeuratObject(counts = ds014.data, project = "asyn", min.cells = 3, min.features = 200)

#Merge to one Seurat object

dsmerged <- merge(ds001, y = c(ds002, ds007, ds008, ds009, ds012, ds013, ds014), add.cell.ids = c("CFA1","asyn1","PBS1","CFA2","asyn2","PBS2","CFA3","asyn3"))

#QC and filtering

dsmerged[["percent.mt"]] <- PercentageFeatureSet(dsmerged, pattern = "^mt-")
VlnPlot(dsmerged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, ncol = 3)
dsmerged <- subset(dsmerged, subset = nCount_RNA > 1000 & nCount_RNA < 50000 & percent.mt < 15)
VlnPlot(dsmerged, features = c("percent.mt"), ncol = 3, split.by = "orig.ident")

#Remove ribo/mito genes

countstable = dsmerged@assays$RNA
metadata = dsmerged@meta.data
keepgenes=grep("^mt|^Rps|^Rpl",rownames(countstable), invert = T, val = T)
dsmerged <- CreateSeuratObject(counts = countstable[keepgenes,],meta.data = metadata)

#Normalize, find variable features

dsmerged <- NormalizeData(dsmerged)
dsmerged <- FindVariableFeatures(dsmerged, selection.method = "vst", nfeatures = 3000)

#Scale the data

all.genes <- rownames(dsmerged)
dsmerged <- ScaleData(dsmerged, features = all.genes)

dsmerged <- RunPCA(dsmerged, features = VariableFeatures(object = dsmerged))
print(dsmerged[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(dsmerged, dims = 1:5, reduction = "pca")

ElbowPlot(dsmerged,ndims=50)

dsmerged <- FindNeighbors(dsmerged, dims = 1:30)
dsmerged <- FindClusters(dsmerged, resolution = 1) 
dsmerged <- RunUMAP(dsmerged, dims = 1:30)
DimPlot(dsmerged, reduction = "umap",label=T) #displays the clusters
DimPlot(dsmerged, reduction = "umap",group.by="orig.ident") #shows the original identity of the cells

#Make batch column

batchvec <- rep("batch1", nrow(dsmerged@meta.data))
batchvec[grep("PBS1", rownames(dsmerged@meta.data))] <- "batch2"
batchvec[grep("CFA2", rownames(dsmerged@meta.data))] <- "batch2"
batchvec[grep("asyn2", rownames(dsmerged@meta.data))] <- "batch2"
batchvec[grep("PBS2", rownames(dsmerged@meta.data))] <- "batch3"
batchvec[grep("CFA3", rownames(dsmerged@meta.data))] <- "batch3"
batchvec[grep("asyn3", rownames(dsmerged@meta.data))] <- "batch3"

# Batch correction with Harmony

dsmerged_bc <- dsmerged %>%
  RunHarmony("batch", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(dsmerged_bc, 'harmony')
harmony_embeddings[1:5, 1:5]

dsmerged_bc <- dsmerged_bc %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 1) %>%
  identity()
DimPlot(dsmerged_bc, reduction = "umap", label=T)
DimPlot(dsmerged_bc, reduction = "umap", group.by="orig.ident")


# Cluster annotations

DimPlot(dsmerged_bc, reduction = "umap",label=T)
new.cluster.ids <- c("Plasma 1","IEL 1","CD8 T cell","Th1","Memory B cell","IEL 2","DC 1","Th17","ILC3","Granulocyte","NK/ILC1","Macrophage","Treg","NK T cell","DC 2","Induced IEL","Plasma 2","ILC2","GD T cell","Prolif T cell","Naive T cell","Monocyte","pDC","Plasma 3","Prolif DC","Unknown","Prolif B cell")
names(new.cluster.ids) <- levels(dsmerged_bc)
dsmerged_bc <- RenameIdents(dsmerged_bc, new.cluster.ids)
DimPlot(dsmerged_bc, reduction = "umap",label=T)+ NoLegend()
dsmerged_bc[["dsmerged_bc_ClusterNames"]] <- Idents(object = dsmerged_bc)


# Subset Th1 and Th17 separately, generate tissue resident memory signature heatmap

dsmerged_bc_Tcells <- subset(x=dsmerged_bc, idents = c("Th17"))

# replace Th17 with Th1 as needed

dsmerged_bc_Tcells <- AddMetaData(dsmerged_bc_Tcells, paste(dsmerged_bc_Tcells@meta.data$dsmerged_bc_ClusterNames,dsmerged_bc_Tcells@meta.data$orig.ident), col.name = "cell_cond")
head(dsmerged_bc_Tcells@meta.data)

Idents(dsmerged_bc_Tcells) <- dsmerged_bc_Tcells@meta.data$orig.ident

orig.levels <- levels(dsmerged_bc_Tcells)
Idents(dsmerged_bc_Tcells) <- gsub(pattern = " ", replacement = "_", x = Idents(dsmerged_bc_Tcells))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(dsmerged_bc_Tcells) <- orig.levels
cluster.averages <- AverageExpression(dsmerged_bc_Tcells, return.seurat = TRUE)
cluster.averages

ResidentGenes <- c("Tnfrsf4","Crem","Vps37b","Cstb","Nkg7","Hspa1a","Ccl5","Areg","H2afz","Sub1","Ctla2a","Ramp3","Fth1","Tnfrsf18","Mt1","H3f3b","Hif1a","Mmd","Fam107b","Hilpda")

DoHeatmap(cluster.averages, features = ResidentGenes, size = 3, disp.min = -1.5, disp.max = 1.5) + scale_fill_gradientn(colors = rev(brewer.pal(n=8, name = "RdBu")))

#Generate Resident Memory T cell feature plot for Th1 and Th17 subsets

dsmerged_bc_Tcells <- AddModuleScore(
  object = dsmerged_bc_Tcells,
  features = ResidentGenes,
  ctrl = 100,
  name = "ResGenes"
)
FeaturePlot(dsmerged_bc_Tcells, features = ResGenes1)


#Cluster annotations dotplot

cell.type.markers <- c("Cd3e","Cd4","Il17a","Il12rb2","Foxp3","Ccr7","Cd163l1","Cd8a","Cd8b1","Itgae","Trdc","Klrb1c","Gata3","Rorc","Zbtb16","Lcn2","Cx3cr1","Ccr2","Itgam","Clec9a","Cd209a","Siglech","Sdc1","Cd19","Cdc20")

Idents(dsmerged_bc) <- dsmerged_bc@meta.data$dsmerged_bc_ClusterNames
levels(dsmerged_bc) <- c("Unknown","Prolif DC","Prolif T cell","Prolif B cell","Memory B cell","Plasma 3","Plasma 2","Plasma 1","pDC","CD103+ CD11b+ DC","CD103+ CD11b- DC","Monocyte","Macrophage","Granulocyte","NK T cell","ILC3","ILC2","NK/ILC1","GD IEL","CD8aa IEL","CD8 Trm","CD8 T cell","GD T cell","Naive T cell","Treg","Th1","Th17")
DotPlot(dsmerged_bc, features = cell.type.markers, 
        cols = c("gray88","dodgerblue3"),
        dot.min = .05) + theme(axis.text.x=element_text(size=14,angle=45,hjust=1),
                               axis.text.y=element_text(size=13), axis.title = element_blank(), legend.text = element_text(size=11), legend.title=element_blank())


#Volcano plot (granulocytes)

dsmerged_grans <- subset(x=dsmerged_bc, idents = c("Granulocyte"))
Idents(dsmerged_grans) <- dsmerged_grans@meta.data$orig.ident

CFAvSYN <- FindMarkers(object = dsmerged_grans, ident.1 = ("asyn"), ident.2 = c("CFA"))

cells <- CFAvSYN[complete.cases(CFAvSYN), ]
ggplot(data = cells, aes(x=avg_log2FC, y=p_val_adj)) + geom_point()

cells <- cells %>%
  mutate(
    Expression = case_when(avg_log2FC >= 0.5 & p_val_adj <= 0.05 ~ "Up-regulated",
                           avg_log2FC <= -0.5 & p_val_adj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

cells_VP <- ggplot(data = cells, aes(x=avg_log2FC, y=-log10(p_val_adj))) + geom_point(aes(color = Expression), size = 1) +
  scale_color_manual(values = c("dodgerblue3","gray75","firebrick3")) +
  theme_minimal()

cells_VP
cells_VP2 <- cells_VP + geom_vline(xintercept=c(-0.5,0.5), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red")
cells_VP2

genes.to.label <- c("Chil3","Stfa2l1","Saa3","Mt1","Napsa","Ahnak","Wfdc21","Ngp","Prnp","Tgm2","Cstdc4","Crispld2","Lcn2","Mmp8","Fpr1","Il18rap","Wfdc17","Fpr2","S100a8","Klf13","Il1rap","S100a6","Cd14","Mt2","Cd33","Klf2","Dgat2","Traf1","Ifitm6","Batf","Trem1","Tnf","Fcer1g","Ifitm1","Ifitm3","Cxcr2","Ctla","H2afz","Ifitm2","Cd52","Adamdec1","Ccr3","Serpinb2","F5","Alox15","Fut4","Lta4h","P2ry14","Gsn","Adgre1","Ikzf3","Rgs1","Cxcl3","Adgre4","Mmp25","Dusp6","Ahr","Itgax","Arrdc3","Il1a","Pecam1","P2ry10","Socs2","Cd9","Il1rn","Vegfa","Fli1","Itga4","Rgs2","Ddx3x","Itgam","Ncf1","Cish")
LabelPoints(plot=cells_VP2, points = genes.to.label, repel = TRUE, xnudge = 0.7, ynudge = 0, max.overlaps = Inf)


#Th1 subset then re-cluster. Replace Th1 with Th17 or Granulocytes, etc as needed

dsmerged_Th1 <- subset(x=dsmerged_bc, idents = c("Th1"))
dsmerged_Th1 <- FindVariableFeatures(dsmerged_Th1, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(dsmerged_Th1)
dsmerged_Th1 <- ScaleData(dsmerged_Th1, features = all.genes)

dsmerged_Th1 <- RunPCA(dsmerged_Th1, features = VariableFeatures(object = dsmerged_Th1))

ElbowPlot(dsmerged_Th1,ndims=50)

dsmerged_Th1 <- FindNeighbors(dsmerged_Th1, dims = 1:18)
dsmerged_Th1 <- FindClusters(dsmerged_Th1, resolution = 0.5) #lower resolution, cluters more general
dsmerged_Th1 <- RunUMAP(dsmerged_Th1, dims = 1:18)
DimPlot(dsmerged_Th1, reduction = "umap")
DimPlot(dsmerged_Th1, reduction = "umap", group.by = "batch")

dsmerged_Th1 <- dsmerged_Th1 %>%
  RunHarmony("batch", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(dsmerged_Th1, 'harmony')
harmony_embeddings[1:5, 1:5]

dsmerged_Th1 <- dsmerged_Th1 %>%
  RunUMAP(reduction = "harmony", dims = 1:18) %>%
  FindNeighbors(reduction = "harmony", dims = 1:18) %>%
  FindClusters(resolution = 0.7) %>%
  identity()

DimPlot(dsmerged_Th1, reduction = "umap", label=T) #displays the clusters
DimPlot(dsmerged_Th1, reduction = "umap", group.by="batch") #shows the original identity of the cells
DimPlot(dsmerged_Th1, reduction = "umap", group.by = "orig.ident")

#Differential expression test

CFAvSYN <- FindMarkers(dsmerged_Th1, ident.1 = "asyn", ident.2 = "CFA", min.pct = 0, logfc.threshold = 0)



