setwd("C:/MSc_Dissertation")

# Load necessary libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(HGNChelper)
library(dplyr)
library(ggplot2)
library(Seurat)
remotes::install_github("satijalab/seurat-data")
library(SeuratData)
library(patchwork)
# load libraries
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

# Load dataset
sample1.data <- Read10X(data.dir = "C:/MSc_Dissertation/Dataset files/GSE229711_RAWNP/NPH51")

# Create Seurat objects and add sample identity
sample1 <- CreateSeuratObject(counts = sample1.data, project = "Sample1", min.cells = 3, min.features = 200)

sample1$sample <- "Sample1"

# Calculate percentage of mitochondrial genes
sample1[["percent.mt"]] <- PercentageFeatureSet(sample1, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(sample1, features = c("nFeature_RNA", "percent.mt"), ncol = 3)


#QC
sample1 <- subset(sample1,
                     subset = nFeature_RNA > 150 & 
                       nFeature_RNA < 2000 & 
                       percent.mt < 3)


VlnPlot(sample1, features = c("nFeature_RNA", "percent.mt"), ncol = 3)


#Normalization
sample1 <- NormalizeData(sample1, normalization.method = "LogNormalize", scale.factor = 10000)
sample1 <- FindVariableFeatures(sample1, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sample1), 50)

# scale and run PCA
all.genes <- rownames(sample1)
sample1 <- ScaleData(sample1, features = all.genes)
sample1 <- RunPCA(sample1, features = VariableFeatures(object = sample1))


# Examine and visualize PCA results a few different ways
print(sample1[["pca"]], dims = 1:5, nfeatures = 5)



DimPlot(sample1, reduction = "pca") + NoLegend()


# Check number of PC components 
ElbowPlot(sample1) #(selected 10 PCs for downstream analysis, based on Elbow plot)

# cluster and visualize
sample1 <- FindNeighbors(sample1, dims = 1:10)
sample1 <- FindClusters(sample1, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(sample1), 5)
sample1 <- RunUMAP(sample1, dims = 1:10)
DimPlot(sample1, reduction = "umap")


#scType Annotation
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

 # DB file
tissue <- "NucleusPulposus" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
library(readxl)

# Read the first sheet (or specify it by name)
db_ <- read_excel("C:/MSc_Dissertation/NP_CellTypeMarkers_ScType.xlsx", sheet = 1)

# Check structure
str(db_)
unique(db_$tissue)  # Now this should work
gene_sets_prepare <- function(db_, tissue) {
  db_filtered <- db_[db_$tissue == tissue, ]
  
  if (nrow(db_filtered) == 0) stop("No entries found for this tissue.")
  
  db_filtered$positive_markers <- strsplit(db_filtered$positive_markers, ",\\s*")
  db_filtered$negative_markers <- strsplit(db_filtered$negative_markers, ",\\s*")
  
  gs_positive <- setNames(db_filtered$positive_markers, db_filtered$cellName)
  gs_negative <- setNames(db_filtered$negative_markers, db_filtered$cellName)
  
  return(list(gs_positive = gs_positive, gs_negative = gs_negative))
}

gs_list <- gene_sets_prepare(db_, tissue)




# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(sample1[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# extract scaled scRNA-seq matrix
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(sample1[["RNA"]]$data) else as.matrix(sample1[["RNA"]]@data)

# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)


# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(sample1@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(sample1@meta.data[sample1@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sample1@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])

sample1@meta.data$sctype_classification = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  sample1@meta.data$sctype_classification[sample1@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(sample1, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification')        


#Assign cell labels
counts <- GetAssayData(sample1, assay = "RNA", layer = "counts")  # ✅
cell_types <- sample1@meta.data$sctype_classification
colnames(counts) <- as.character(cell_types)  # this causes R to make unique names

df <- as.data.frame(counts)
df <- cbind(Gene = rownames(df), df)
# If you already have your counts in a data frame with gene names as a column
write.table(df, file = "CIBERF.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Two way ANOVA
A <- c("Chondrocytes","Chondrocytes","Chondrocytes","Progenitors","Progenitors","Progenitors","Chondrocytes","Chondrocytes","Chondrocytes","Progenitors","Progenitors","Progenitors","Chondrocytes","Chondrocytes","Chondrocytes","Progenitors","Progenitors","Progenitors")
B <- c("FBS","FBS","FBS","FBS","FBS","FBS","Autoserum","Autoserum","Autoserum","Autoserum","Autoserum","Autoserum","platelet_lysate","platelet_lysate","platelet_lysate","platelet_lysate","platelet_lysate","platelet_lysate")
DV <- c(0.786,0.529,0.863,0.214,0.471,0.137,0.634,0.684,0.736,0.366,0.311,0.264,0.571,0.471,0.552,0.429,0.529,0.448)
ID <- c(0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5)
df1 <- data.frame(A,B,DV,ID)
Model1 <- aov(DV ~ A + B + A:B, data=df1)
summary(Model1)
res=residuals(object = Model1)


#Boxplot
ggplot(df1, aes(x = interaction(A, B), y = DV, fill = A)) +
      geom_boxplot() +
      labs(x = "A x B", y = "DV", title = "Two-way ANOVA Boxplot") +
      theme_minimal()
