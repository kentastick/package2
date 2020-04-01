library(Seurat)
library(tidyverse)
library(package2)


data <- readRDS("data/Seurat_object/Aizarani.rds")




# pbc case 1 --------------------------------------------------------------
setwd("C:/Users/ken/Documents/single_cell/single_cell_project/data/pbc_case1")

data <- Read10X(data.dir = ".")

object <- CreateSeuratObject(counts = data, project = "object", min.cells = 3, min.features = 200)
object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("plot1.jpg", device = "jpeg")

plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
ggsave("plot2.jpg", device = "jpeg")

#object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
ggplot2::ggsave("plot3.jpg")

all.genes <- rownames(object)

#normalization data mean =1 variaty =0
object <- ScaleData(object, features = all.genes)
object <- RunPCA(object, features = VariableFeatures(object = object))

print(object[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(object, dims = 1:2, reduction = "pca")
ggsave("plot4.jpg")

DimPlot(object, reduction = "pca")
ggsave("plot5.jpg")

DimHeatmap(object, dims = 15:30, cells = 500, balanced = TRUE)
ggsave("plot6.jpg")


#jackstraw analysis
object <- JackStraw(object, num.replicate = 100)
object <- ScoreJackStraw(object, dims = 1:20)

JackStrawPlot(object, dims = 1:20)
ggsave("plot7.jpg")



object <- FindNeighbors(object, dims = 1:27)
object <- FindClusters(object, resolution = 0.8)
object <- RunUMAP(object, dims = 1:30)
object <- RunTSNE(object, dims = 1:30)

DimPlot(object, reduction = "umap")
ggsave("plot8.jpg")

DimPlot(object, reduction = "tsne")
ggsave("plot9.jpg")

data <- object
rm(object)
tmap(c('ALB', 'KRT7', 'CD68', 'CD3D', 'CD79A', 'IGJ','CD34', 'GNLY', 'FGFBP2', 'CD14', "MZB1"))
ggsave("plot10.jpg")


saveRDS(data, "pbc_case1.rds")








# pbc case2 ---------------------------------------------------------------

setwd("C:/Users/ken/Documents/single_cell/single_cell_project/data/pbc_case2/")

data <- read_delim("pbc2_table_th1000andOver_list_7-Ecut-human_S4.txt","\t", escape_double = FALSE, trim_ws = TRUE)
typeof(data)
class(data)
data <- data %>% as.data.frame()
rownames(data) <- data[[1]]
data %>% dim
data <- data[-1]
data %>% dim


object <- CreateSeuratObject(counts = data, project = "object", min.cells = 3, min.features = 200)
object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("plot1.jpg", device = "jpeg")

plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
ggsave("plot2.jpg", device = "jpeg")

#object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
ggplot2::ggsave("plot3.jpg")

all.genes <- rownames(object)

#normalization data mean =1 variaty =0
object <- ScaleData(object, features = all.genes)
object <- RunPCA(object, features = VariableFeatures(object = object))

print(object[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(object, dims = 1:2, reduction = "pca")
ggsave("plot4.jpg")

DimPlot(object, reduction = "pca")
ggsave("plot5.jpg")

DimHeatmap(object, dims = 15:30, cells = 500, balanced = TRUE)
ggsave("plot6.jpg")


#jackstraw analysis
object <- JackStraw(object, num.replicate = 100)
object <- ScoreJackStraw(object, dims = 1:20)

JackStrawPlot(object, dims = 1:20)
ggsave("plot7.jpg")



object <- FindNeighbors(object, dims = 1:27)
object <- FindClusters(object, resolution = 0.8)
object <- RunUMAP(object, dims = 1:30)
object <- RunTSNE(object, dims = 1:30)

DimPlot(object, reduction = "umap")
ggsave("plot8.jpg")

DimPlot(object, reduction = "tsne")
ggsave("plot9.jpg")

data <- object
rm(object)
tmap(c('ALB', 'KRT7', 'CD68', 'CD3D', 'CD79A', 'IGJ','CD34', 'GNLY', 'FGFBP2', 'CD14', "MZB1"))
ggsave("plot10.jpg")


saveRDS(data, "pbc_case2.rds")


# GSE125449 carcinoma----------------------------------------------------------------
#trace(Read10X, edit=TRUE)

#setwd("C:/Users/ken/Documents/single_cell/single_cell_project/data/GSE125449/GSE125449_set1/")
setwd("C:/Users/ken/Documents/single_cell/single_cell_project/data/GSE125449/GSE125449_set2/")
data <- Read10X(data.dir = ".")

#add sample meta_data to cell_id from GSE125448_set1_samples.txt
#label <- read.delim("GSE125449_Set1_samples.txt/GSE125449_Set1_samples.txt")
label <- read.delim("GSE125449_Set2_samples.txt/GSE125449_Set2_samples.txt")
label %>% colnames
sample_name <- label$Sample %>% unique() %>% as.character()

#set1
label$Sample <- label$Sample %>% plyr::mapvalues(from = sample_name,
                                 to = c("HCC_1", "HCC_2", "HCC_3", "ICC_1", "ICC_2", "HCC_4", "ICC_3", "HCC_5", "ICC_4", "HCC_6","HCC_7", "ICC_5"))

sample_name
#set2
label$Sample <- label$Sample %>% plyr::mapvalues(from = sample_name,
                                 to = c("HCC_8", "ICC_6", "ICC_7", "ICC_8", "ICC_9", "HCC_9", "ICC_10"))

identical(colnames(ma_set2), as.character(label$Cell.Barcode))


#ma_set1$id <- label$Sample
ma_set2$id <- label$Sample

#sav(ma_set1)
sav(ma_set2)


#
object <- CreateSeuratObject(counts = data, project = "object", min.cells = 3, min.features = 200)
object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("plot1.jpg", device = "jpeg")

plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
ggsave("plot2.jpg", device = "jpeg")

#object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
ggplot2::ggsave("plot3.jpg")

all.genes <- rownames(object)

#normalization data mean =1 variaty =0
object <- ScaleData(object, features = all.genes)
object <- RunPCA(object, features = VariableFeatures(object = object))

print(object[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(object, dims = 1:2, reduction = "pca")
ggsave("plot4.jpg")

DimPlot(object, reduction = "pca")
ggsave("plot5.jpg")

DimHeatmap(object, dims = 15:30, cells = 500, balanced = TRUE)
ggsave("plot6.jpg")


#jackstraw analysis
object <- JackStraw(object, num.replicate = 100)
object <- ScoreJackStraw(object, dims = 1:20)

JackStrawPlot(object, dims = 1:20)
ggsave("plot7.jpg")



object <- FindNeighbors(object, dims = 1:27)
object <- FindClusters(object, resolution = 0.8)
object <- RunUMAP(object, dims = 1:30)
object <- RunTSNE(object, dims = 1:30)

DimPlot(object, reduction = "umap")
ggsave("plot8.jpg")

DimPlot(object, reduction = "tsne")
ggsave("plot9.jpg")

data <- object
rm(object)
tmap(c('ALB', 'KRT7', 'CD68', 'CD3D', 'CD79A', 'IGJ','CD34', 'GNLY', 'FGFBP2', 'CD14', "MZB1"))
ggsave("plot10.jpg")

saveRDS(data, "normal.rds")

# GSE124395 from Aizarani --------------------------------------------------------------
setwd("C:/Users/ken/Documents/single_cell/single_cell_project/data/GSE124395")
data <- readRDS("GSE124395_Normalhumanliverdata.RData")

object <- CreateSeuratObject(counts = data, project = "object", min.cells = 3, min.features = 200)

object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("plot1.jpg", device = "jpeg")

plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
ggsave("plot2.jpg", device = "jpeg")

#object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
ggplot2::ggsave("plot3.jpg")

all.genes <- rownames(object)

#normalization data mean =1 variaty =0
object <- ScaleData(object, features = all.genes)
object <- RunPCA(object, features = VariableFeatures(object = object))

print(object[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(object, dims = 1:2, reduction = "pca")
ggsave("plot4.jpg")

DimPlot(object, reduction = "pca")
ggsave("plot5.jpg")

DimHeatmap(object, dims = 15:30, cells = 500, balanced = TRUE)
ggsave("plot6.jpg")


#jackstraw analysis
object <- JackStraw(object, num.replicate = 100)
object <- ScoreJackStraw(object, dims = 1:20)

JackStrawPlot(object, dims = 1:20)
ggsave("plot7.jpg")



object <- FindNeighbors(object, dims = 1:27)
object <- FindClusters(object, resolution = 0.8)
object <- RunUMAP(object, dims = 1:30)
object <- RunTSNE(object, dims = 1:30)

DimPlot(object, reduction = "umap")
ggsave("plot8.jpg")

DimPlot(object, reduction = "tsne")
ggsave("plot9.jpg")

data <- object
rm(object)
tmap(c('ALB', 'KRT7', 'CD68', 'CD3D', 'CD79A', 'IGJ','CD34', 'GNLY', 'FGFBP2', 'CD14', "MZB1"))
ggsave("plot10.jpg")


saveRDS(data, "normal.rds")

marker <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write_csv(marker, path = "marker.csv")



# GSE115496Normal Liver from coudate Macpoland -------------------------------------
setwd("C:/Users/ken/Documents/single_cell/single_cell_project/data/GSE115469/count")

data <- read_csv(file = "GSE115469_Data.csv")
# data <- data %>% as.data.frame(data)
# rownames(data) <- data$X1
# data %>% dim
# data <- data[-1]
# data %>% dim
# data <- (2**data)-1


object <- CreateSeuratObject(counts = data, project = "object", min.cells = 3, min.features = 200)

object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("plot1.jpg", device = "jpeg")

plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
ggsave("plot2.jpg", device = "jpeg")

#object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
ggplot2::ggsave("plot3.jpg")

all.genes <- rownames(object)

#normalization data mean =1 variaty =0
object <- ScaleData(object, features = all.genes)
object <- RunPCA(object, features = VariableFeatures(object = object))

print(object[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(object, dims = 1:2, reduction = "pca")
ggsave("plot4.jpg")

DimPlot(object, reduction = "pca")
ggsave("plot5.jpg")

DimHeatmap(object, dims = 15:30, cells = 500, balanced = TRUE)
ggsave("plot6.jpg")


#jackstraw analysis
object <- JackStraw(object, num.replicate = 100)
object <- ScoreJackStraw(object, dims = 1:20)

JackStrawPlot(object, dims = 1:20)
ggsave("plot7.jpg")



object <- FindNeighbors(object, dims = 1:27)
object <- FindClusters(object, resolution = 0.8)
object <- RunUMAP(object, dims = 1:30)
object <- RunTSNE(object, dims = 1:30)

DimPlot(object, reduction = "umap")
ggsave("plot8.jpg")

DimPlot(object, reduction = "tsne")
ggsave("plot9.jpg")

data <- object
rm(object)
tmap(c('ALB', 'KRT19', "CD163","CD32",'CD68', 'CD3D', 'CD79A', 'IGJ','CD34', 'GNLY', 'FGFBP2', 'CD14', "MZB1"))
ggsave("plot10.jpg")


saveRDS(data, "macpoland.rds")

#hepatocyte extraction

data <- readRDS(file = "macpoland.rds")


tmap("ALB", min.cutoff = 3)
ggsave("alb.jpg")

marker <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write_csv(marker, path = "marker.csv")

top10 <- marker %>% group_by(cluster) %>% top_n(10, avg_logFC)

top10 %>% write_clip()

readRDS(file = "macpoland.rds")

signature_gene_list <- read_clip_tbl(header = F)
signature_gene_list <- signature_gene_list %>% column_to_rownames("V1") %>% df_to_list()

saveRDS(signature_gene_list, "../../signature_list.rds")


signature_plot(gene_list = gene_list)

#cluster 0, 1, 2, 8, 16, 19,21

hep_1 <- subset(data, idents = c(0, 1, 2, 8, 16, 19, 21 ))






# GSE136103 NASH liver data process ----------------------------------------------------------------
#data import
setwd("C:/Users/ken/Documents/single_cell/single_cell_project/data/GSE136103/cd45-ch/")


alldata <- list.files(pattern="matrix")
headname <- str_extract(alldata, ".*(?=_matrix.mtx.gz)")
headname <- headname[1:24]
temp_genes1 <- read.table(paste0(headname[1], "_genes.tsv.gz"),sep="\t",header = FALSE)
temp_cells <- NULL
nrow(temp_genes1)
eb_raw<-matrix(nrow = nrow(temp_genes1))



#combined each matrix
for(i in seq_along(headname)){
  temp_cells1 <- scan(paste0(headname[i],"_barcodes.tsv.gz"),what = character(), sep="\t")
  temp_cells1 <- paste0(headname[i],"_", temp_cells1)
  temp_cells1 <- gsub("-","_",temp_cells1)
  temp_cells <- c(temp_cells,temp_cells1)
  raw1<-readMM(paste0(headname[i],"_matrix.mtx.gz"))
  eb_raw <- cbind(eb_raw,raw1)

}

eb_raw<-eb_raw[,2:ncol(eb_raw)]
colnames(eb_raw) <- temp_cells
rownames(eb_raw) <- temp_genes1[,2]

#Remove duplicated gene names (a couple genes are in under their MGI and HGNC symbols)
temp_r <- rownames(eb_raw)[which(duplicated(toupper(rownames(eb_raw))))]
temp_r <- lapply(temp_r,function(X) grep(paste0("^",X,"$"),rownames(eb_raw),ignore.case=T))
temp_r <- which(rownames(eb_raw) %in%
                  names(sapply(temp_r,function(X) which.min(apply(eb_raw[X,],1,function(Y) sum(Y>0))))))
if (length(temp_r) > 0) { eb_raw <- eb_raw[-temp_r,] }
print(temp_r)

saveRDS(eb_raw, "gene_count_matrix.rds")

#pick up only cirrhosis and parencime cells

eb_raw <- readRDS("../gene_count_matrix.rds")
sub_raw <- str_detect(colnames(eb_raw), "cirrhotic\\d_cd45_") %>% eb_raw[, .]
sub_raw <- str_detect(colnames(eb_raw), "healthy\\d_cd45_") %>% eb_raw[, .]
sub_raw <- str_detect(colnames(eb_raw), "cirrhotic\\d_cd45+") %>% eb_raw[, .]
sub_raw <- str_detect(colnames(eb_raw), "healthy\\d_cd45+") %>% eb_raw[, .]
sub_raw <- str_detect(colnames(eb_raw), "blood\\d") %>% eb_raw[, .]


#seurat procedure

# procidure from ----------------------------------------------------------------------

setwd("C:/Users/ken/Documents/single_cell/single_cell_project/data/GSE136103/cd45+ht//")
sub_raw %>% dim
object <- CreateSeuratObject(counts = sub_raw, project = "object", min.cells = 3, min.features = 200)

object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("plot1.jpg", device = "jpeg")

plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
ggsave("plot2.jpg", device = "jpeg")

#object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
ggplot2::ggsave("plot3.jpg")

all.genes <- rownames(object)

#normalization data mean =1 variaty =0
object <- ScaleData(object, features = all.genes)
object <- RunPCA(object, features = VariableFeatures(object = object))

print(object[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(object, dims = 1:2, reduction = "pca")
ggsave("plot4.jpg")

DimPlot(object, reduction = "pca")
ggsave("plot5.jpg")

DimHeatmap(object, dims = 15:30, cells = 500, balanced = TRUE)
ggsave("plot6.jpg")


#jackstraw analysis
object <- JackStraw(object, num.replicate = 100)
object <- ScoreJackStraw(object, dims = 1:20)

JackStrawPlot(object, dims = 1:20)
ggsave("plot7.jpg")



object <- FindNeighbors(object, dims = 1:30)
object <- FindClusters(object, resolution = 0.8)
object <- RunUMAP(object, dims = 1:30)
object <- RunTSNE(object, dims = 1:30)

DimPlot(object, reduction = "umap")
ggsave("plot8.jpg")

DimPlot(object, reduction = "tsne")
ggsave("plot9.jpg")

data <- object
rm(object)
tmap(c('ALB', 'KRT7', 'CD68', 'CD3D', 'CD79A', 'IGJ','CD34', 'GNLY', 'FGFBP2', 'CD14', "MZB1"))
ggsave("plot10.jpg")


saveRDS(data, "data.rds")


# to ----------------------------------------------------------------------


#pickup hepatocytes
setwd("C:/Users/ken/Documents/single_cell/single_cell_analysis/data/GSE136103/cd45-ch/")

data <- readRDS("cd45-ch.rds")
signature_plot(gene_list = gene_list)
val <- sig_val(gene_list)
data@meta.data <- data@meta.data %>% cbind(val)
tmap("Cholangiocyte")

#read exel file
my_data <- readxl::read_excel("../../combined_0423/data/liver_marker.xlsx", sheet = 1, col_names = T, skip = 1)
my_data$cluster <- sheets_name[[1]]

for(i in 2:20){
  temp <- readxl::read_excel("../../combined_0423/data/liver_marker.xlsx", sheet = i, col_names = T, skip = 1)
  temp$cluster <- sheets_name[i]
  my_data <- my_data %>% rbind(temp)
}

#get sheets name
abpath("data/")
sheets_name <- readxl::excel_sheets(path = "../../combined_0423/data/liver_marker.xlsx")


a <- my_data %>% group_by(cluster) %>% do(head(., 10)) %>% select(cluster, Genes)
b <- a %>% group_by(cluster) %>% nest %>% mutate(data = (map(data, ~.$Genes))) %>% .$data
names(b) <- sheets_name

signature_plot(gene_list)


# GSE 130473 --------------------------------------------------------------

setwd("C:/Users/ken/Documents/single_cell/single_cell_project/data/GSE130473/")
data <- read.csv(file = "C:/Users/ken/Documents/single_cell/single_cell_project/data/GSE130473/GSE130473_Series_count_matrix.csv/GSE130473_Series_count_matrix.csv")
rownames(data) <- data[[1]]
data <- data[-1]

#ensembl id to symbol

#remove ercc
ensemblsIDS <- rownames(data)
ercc_nu <- ensemblsIDS %>% str_which("^ERCC-")
length(ercc_nu)
ensemblsIDS %>% length
ensemblsIDS <- ensemblsIDS[-ercc_nu]
ensemblsIDS %>% length
ensemblsIDS <- ensemblsIDS %>% str_extract(".*(?=\\.\\d{1,})")
ensemblsIDS %>% duplicated() %>% which()

#remove ercc_data from original matrix
data <- data[-ercc_nu,]
rownames(data) %>% length
ensemblsIDS %>% length
rownames(data) <- ensemblsIDS

#org.hg. method

library(org.Hs.eg.db)
symbols <- mapIds(org.Hs.eg.db, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")
gene_id <- as.character(names(unlist(symbols)))
length(gene_id)
unlist(gene_id) %>% is.na() %>% which
setdiff(ensemblsIDS, gene_id)
unlist(symbols) %>% length
identical(ensemblsIDS, gene_id)

gene_symbol <- unlist(symbols) %>% as.character()
rownames(data)[1:length(gene_symbol)] <- gene_symbol

#biomaRT method
library("biomaRt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

res2 <- getBM(attributes = c(
                            'external_gene_name',
                            'hgnc_symbol',
                            'ensembl_gene_id'),
             filters = 'ensembl_gene_id',
             values = ensemblsIDS,
             mart = ensembl)

res2$gene %>% duplicated() %>% which %>% res2$gene[.]
res2$id %>% duplicated() %>% which %>% res2$id[.]
dup_row <- res2$id %>% duplicated() %>% which
res2 <- res2[-dup_row,]
res2 <- res2 %>% dplyr::select(gene = external_gene_name, id = ensembl_gene_id)

# #detect non detected gene id
# setdiff(ensemblsIDS, res2$ensembl_gene_id) -> non_detected_gene
# non_detected_gene
# #convert at https://www.biotools.fr/mouse/ensembl_symbol_converter
# conv_non_detected_gene <- read_clip_tbl(header = F)
#
# #detct overwrap between two
# over_wrap_gene <- intersect(res2$gene, conv_non_detected_gene$gene)
#
# #filter only gene converted
# conv_non_detected_gene <- conv_non_detected_gene %>% filter(!is.na(V2))
# conv_non_detected_gene <- conv_non_tedetected_gene %>% dplyr::select(gene = V2, id = V1) %>%
#   filter(!gene %in% over_wrap_gene)
#
# conv_non_detected_gene$gene %>% duplicated() %>% which %>% conv_non_detected_gene$gene[.]
#
# res2 %>% nrow #56840
# #no shared gene for confirmation
# any(res2$id %in% non_detected_gene)
#
# #combine two detected data frame
#
# res <- res2 %>% bind_rows(conv_non_detected_gene)
#
# res$gene %>% duplicated() %>% which %>% res$gene[.]
#
#
#

 a <- tibble(id = rownames(data)) %>% left_join(as.tibble(res2), by = "id") %>% pull(gene)

detected_row <- which(!is.na(a))
data <- data[detected_row,]
 data<- as.matrix(data)
 rownames(data) <- a[detected_row]

#remove_duplicated row
temp_r <- rownames(data)[which(duplicated(rownames(data)))]
temp_r <- lapply(temp_r,function(X) grep(paste0("^",X,"$"),rownames(data),ignore.case=T))
temp_r <- which(rownames(data) %in%
                  names(sapply(temp_r,function(X) which.min(apply(data[X,],1,function(Y) sum(Y>0))))))

if (length(temp_r) > 0) { data <- data[-temp_r,] }

saveRDS(data, "filtered_matrix.rds")


object <- CreateSeuratObject(counts = data, project = "object", min.cells = 3, min.features = 200)

object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("plot1.jpg", device = "jpeg")

plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
ggsave("plot2.jpg", device = "jpeg")

#object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
ggplot2::ggsave("plot3.jpg")

all.genes <- rownames(object)

#normalization data mean =1 variaty =0
object <- ScaleData(object, features = all.genes)
object <- RunPCA(object, features = VariableFeatures(object = object))

print(object[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(object, dims = 1:2, reduction = "pca")
ggsave("plot4.jpg")

DimPlot(object, reduction = "pca")
ggsave("plot5.jpg")

DimHeatmap(object, dims = 15:30, cells = 500, balanced = TRUE)
ggsave("plot6.jpg")


#jackstraw analysis
object <- JackStraw(object, num.replicate = 100)
object <- ScoreJackStraw(object, dims = 1:20)

JackStrawPlot(object, dims = 1:20)
ggsave("plot7.jpg")



object <- FindNeighbors(object, dims = 1:27)
object <- FindClusters(object, resolution = 0.8)
object <- RunUMAP(object, dims = 1:30)
object <- RunTSNE(object, dims = 1:30)

DimPlot(object, reduction = "umap")
ggsave("plot8.jpg")

DimPlot(object, reduction = "tsne")
ggsave("plot9.jpg")

data <- object
rm(object)
tmap(c('ALB', 'KRT19', "CD163","CD32",'CD68', 'CD3D', 'CD79A', 'IGJ','CD34', 'GNLY', 'FGFBP2', 'CD14', "MZB1"))
ggsave("plot10.jpg")


saveRDS(data, "segal.rds")


# mouse atlas -------------------------------------------------------------


expr_mat <- read.table("Liver1_dge.txt", header = T, sep = " ")
sav(expr_mat)
mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse_gene <- rownames(expr_mat)
genesV2 = biomaRt::getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mouse_gene , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
sav(genesV2)

expr_mat %>% colnames
expr_mat %>% rownames

genesV2 <- genesV2 %>% group_by(MGI.symbol) %>% do((head(., 1)))

expr_mat_fil <- expr_mat[rownames(expr_mat) %in% genesV2$MGI.symbol, ]
genesV2$rowname <-  rownames(expr_mat_fil)

gene <- tibble(rowname = rownames(expr_mat_fil))
correspond_table <- gene %>% left_join(genesV2, by = c("rowname" = "MGI.symbol"))
sav(correspond_table)
expr_mat_fil <- as.matrix(expr_mat_fil)
rownames(expr_mat_fil) <- correspond_table$HGNC.symbol

temp_r <- rownames(expr_mat_fil)[which(duplicated(toupper(rownames(expr_mat_fil))))]
temp_r <- lapply(temp_r,function(X) grep(paste0("^",X,"$"),rownames(expr_mat_fil),ignore.case=T))
temp_r <- which(rownames(expr_mat_fil) %in%
                  names(sapply(temp_r,function(X) which.min(apply(expr_mat_fil[X,],1,function(Y) sum(Y>0))))))
if (length(temp_r) > 0) { expr_mat_fil <- expr_mat_fil[-temp_r,] }

data <- do_seurat(data = expr_mat_fil)


#use homolog table
a <- read.table("HOM_MouseHumanSequence.rpt", header = T, sep = "\t")
a <- a %>% mutate(orig = case_when(Common.Organism.Name == "mouse, laboratory"~"mouse",
                              TRUE~"human")) %>%
  dplyr::select(id = HomoloGene.ID, orig, symbol = Symbol)
a %>% group_by(id) %>% mutate(n = n()) %>% filter(n ==2) %>%
  dplyr::select(-n) %>%
  mutate(symbol = as.character(symbol)) %>%
  pivot_wider(names_from = orig,  values_from = symbol) ->b
 b <- b %>% mutate_all(as.character)


 expr_mat <- read.table("Liver2_dge.txt", header = T, sep = " ")
 sav(expr_mat)
expr_mat_fil2 <- expr_mat[rownames(expr_mat) %in% b$mouse, ]
expr_mat_fil2 %>% rownames %>% View

label <- tibble(rowname = rownames(expr_mat_fil2))
label <- label %>% left_join(b, by = c("rowname" = "mouse"))
correspond_table2 <- label
sav(correspond_table2)
rownames(expr_mat_fil2) <- correspond_table2$human
sav(expr_mat_fil2)
rownames(expr_mat_fil2) %>% duplicated %>% any()

data <- do_seurat(expr_mat_fil2)
