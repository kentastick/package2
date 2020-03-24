library(package2)
library(tidyverse)
library(purrr)
library(Seurat)
library(clipr)
library(gghighlight)


#add signature value to metadata


df <- sig_val(gene_list = gene_list, object = data)
data@meta.data <- data@meta.data %>% rownames_to_column(var = "rowname") %>%
  bind_cols(df) %>%
  column_to_rownames(var = "rowname")

use_meta_data <- data@meta.data %>% keep(is.numeric) %>% colnames
tmap(use_meta_data)



#calculate mean value of each gene within gene_list

gene_name <-  data@assays$RNA@data@Dimnames[[1]]


data@assays$RNA@data %>% as.matrix %>% rowSums()
data@assays$RNA@data[unlist(gene_list)[unlist(gene_list) %in% gene_name],] %>%
  apply(., 1,FUN = mean)



df %>% gather(-cluster, key = "signature", value = "value") %>%
  ggplot(aes(cluster, value)) + geom_jitter() + facet_wrap(~signature)





# extract specific cells from each seurat object ------------------------------


#list up of analysis data name
data_list <- list.files(pattern = ".rds")

#data_list <- data_list[!str_detect(data_list, "posi|blood")] #remove cd45posi(non-parenchyme cells include)

#dir_name <- data_list %>% str_extract("(?<=\\/)\\S{1,12}(?=.rds)")
data_name <- data_list %>% str_extract("\\S{1,20}(?=.rds)")
#subset_name <- paste0(dir_name, "_hepato_subset") #hepatocyte extract

dir_name <- file.path("HSC_combined", data_name)

#subset_name <- paste0(data_name, "_hepato_subset") # extract
subset_name <- paste0(data_name, "_mesencyme_subset") # extract


i = 1
for(i in seq_along(data_list)[3:4]){

data <- readRDS(data_list[i])

dir.create(dir_name[i])
#signature_plot

df <- sig_val(gene_list = gene_list, object = data, func = "me")

df2 <- sig_val2(score_mt = df, object = data, gene_list = gene_list)

signature_plot(df2)

ggsave(filename = paste0(dir_name[i], "/signature_plot.jpg"))

#calculate mean value of each signature in whole cells.
val_mean <- apply(df, 2, mean)

#select cluster hepatocyte val over 0.2
#hepato_cluster_no <- df2 %>% filter(signature == "Hepatocyte", fraction_of_cells>0.2) %>% pull(cluster)

use_cluster_no <- df2 %>% group_by(cluster) %>% mutate(no = row_number(-score)) %>% filter(signature == "Mesenchyme", no ==1) %>% pull(cluster)

#mesenchyme_cluster_no <- df2 %>% filter(signature == "Mesenchyme", fraction_of_cells>0.2) %>% pull(cluster)


# if(length(hepato_cluster_no) ==0){
#   hepato_cluster_no <- df2 %>% filter(signature == "Hepatocyte", fraction_of_cells>0) %>% pull(cluster)
# }
# if(length(mesenchyme_cluster_no) ==0){
#   hepato_cluster_no <- df2 %>% filter(signature == "Mesenchyme", fraction_of_cells>0) %>% pull(cluster)
# }


df$cluster <- data@meta.data$seurat_clusters


df %>% rownames_to_column("var" = "id") %>%
  filter(Hepatocyte> val_mean["Mesenchyme"], cluster %in% use_cluster_no) %>%
  pull(id) -> use_id

#select cells which have hepatocyte value more than 0 and belong to candidate hepatocyte cluster
# df %>% rownames_to_column("var" = "id") %>%
#   filter(Hepatocyte> 0, cluster %in% hepato_cluster_no) %>%
#   pull(id) -> hepato_id2


#make subset_object of hepato_id cells
sub_data <- subset(data, cells = use_id)
ts(sub_data)
ggsave(paste0(dir_name[i], "/subset_plot.jpg"))
tmap(object = sub_data, features =  gene_list[["Mesenchyme"]])
ggsave(paste0(dir_name[i], "/feature_plot.jpg"))
#save as a hepatocyte_subset object
saveRDS(sub_data, file = paste0(dir_name[i], "/",subset_name[i],"_Mesencyme", ".rds"))

}




# making combined cluster -------------------------------------------------

file_dir <- list.files(pattern = "hepato_subset", recursive = T)

file_name <- file_dir %>% str_split("/") %>% map(~.[1]) %>% unlist

for(i in seq_along(file_dir)){
  assign(x =file_name[i], readRDS(file_dir[i]))
}

#add meta info to each object
aizarani$batch <- "aizarani"
ramachandran_cd45nega_ch$batch <- "chandran_cd45nega_ch"
ramachandran_cd45nega_ht$batch <- "chandran_cd45nega_ht"
ramachandran_cd45posi_ch$batch <- "chandran_cd45posi_ch"
ramachandran_cd45posi_ht$batch <- "chandran_cd45posi_ht"


ramachandran_blood$batch <- "ramachandran_blood"
macpoland$batch <- "macpoland"
pbc_case1$batch <- "pbc_case1"
pbc_case2$batch <- "pbc_case2"





#normalization and gene detection

object.list <- list(pbc_case1, pbc_case2, aizarani, chandran_cd45nega_ch, chandran_cd45nega_ht, macpoland)
#remove other file for saving memory
rm(list = file_name)

object.list <- lapply(X = object.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

object.anchors <- FindIntegrationAnchors(object.list = object.list[c(2,1,3:6)], dims = 1:20)
#suceeded by changing object.list order but don't know the reason

object.combined <- IntegrateData(anchorset = object.anchors, dims = 1:20)

DefaultAssay(object.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
object.combined <- ScaleData(object.combined, verbose = FALSE)
object.combined <- RunPCA(object.combined, npcs = 30, verbose = FALSE)

# t-SNE and Clustering
object.combined <- RunUMAP(object.combined, reduction = "pca", dims = 1:20)
object.combined <- RunTSNE(object.combined, reduction = "pca", dims = 1:20)
object.combined <- FindNeighbors(object.combined, reduction = "pca", dims = 1:20)
object.combined <- FindClusters(object.combined, resolution = 0.5)

DimPlot(object.combined, reduction = "umap", label = TRUE)
DimPlot(object.combined, reduction = "tsne", label = TRUE)

DimPlot(object.combined, reduction = "umap", group.by = "batch")

tmap(gene_list$Mesenchyme, object.combined)
saveRDS(object.combined, file = "data_combined.rds")


# KRT7 analysis -----------------------------------------------------------
up()
ts()

ump(c("KRT7", "KRT19", "ALB", "CRP"))
tmap(c("KRT7", "KRT19", "ALB", "CRP"))
gene_list

df <- sig_val(gene_list = gene_list, func = "me")
df_gm <- sig_val(gene_list = gene_list,func = "gm_mean")

#add signature_value (mean, gm_mean) to seurat metadata
colnames(df_gm) <- paste0(colnames(df),"_gm")
data@meta.data <- data@meta.data %>%
  rownames_to_column(var = "name") %>% bind_cols(df) %>%
  column_to_rownames(var = "name")
data@meta.data <- data@meta.data %>%
  rownames_to_column(var = "name") %>% bind_cols(df_gm) %>%
  column_to_rownames(var = "name")

ump(c("Hepatocyte_gm", "Hepatocyte",  "Cholangiocyte_gm", "Cholangiocyte"))
markers <- FindAllMarkers(data, min.pct = 0.25, logfc.threshold = log(1.5), only.pos = T)

markers %>% group_by(cluster) %>% top_n(10, avg_logFC) %>% View()

saveRDS(markers, "markers.rds")
markers %>% group_by(cluster) %>% top_n(10, avg_logFC) %>% View

ts(group.by = "batch")
ts()




# HSC analysis ------------------------------------------------------------

#use make_subset function as a new tool for extracting specific cells from each seurat object

data_list = list.files(pattern = ".rds")


#execute this seurat_object folda

make_subset(data_list = data_list, "HSC_combined", cell_type = "Mesenchyme", func = "me")


#make list of subset_seurat object

object_list <- make_list("Mesenchyme")

k.filter <- min(200, min(sapply(object_list, ncol)))

#object_list_sub <- object_list[which(object_list %>% map(~.@meta.data %>% nrow()) %>% unlist() > 100)]

com <- combined(object.list = object_list, cell_type = "HSC", k.filter = k.filter)


#these procedure on combined function
# object.anchors <- FindIntegrationAnchors(object.list = object_list[c(9, 1:8, 10:11)], k.filter = 10,dims = 1:20)
# object.combined <- IntegrateData(anchorset = object.anchors, dims = 1:20)
# DefaultAssay(object.combined) <- "integrated"
# object.combined <- ScaleData(object.combined, verbose = FALSE)
# object.combined <- RunPCA(object.combined, npcs = 30, verbose = FALSE)
# object.combined <- RunUMAP(object.combined, reduction = "pca", dims = 1:20)
# object.combined <- RunTSNE(object.combined, reduction = "pca", dims = 1:20)
# object.combined <- FindNeighbors(object.combined, reduction = "pca", dims = 1:20)
# object.combined <- FindClusters(object.combined, resolution = 0.5)
# DimPlot(object.combined, reduction = "umap", label = TRUE)
# ggsave("umap.jpg", device = "jpeg")
# DimPlot(object.combined, reduction = "tsne", label = TRUE)
# DimPlot(object.combined, reduction = "umap", group.by = "batch")

data <- com

signature_plot(marker = "macpoland", n_liver_marker = 10)


library(clipr)
hsc_marker <- FindAllMarkers(data,min.pct = 0.25,logfc.threshold = log(2), only.pos = T)

hsc_marker %>% group_by(cluster) %>% top_n(20, avg_logFC) %>% write_clip()

saveRDS(hsc_marker, "hsc_marker.rds")




#select only mesenchyme cells

choose_cell <- up() %>% CellSelector()

data <- subset(data, ident.remove = c("3", "7", "8", "9", "10", "12" ))
saveRDS(data, "hsc_combined_filtered.rds")


#cell origin plot
data@meta.data %>% dplyr::select(seurat_clusters, batch) %>%
  pivot_longer(batch,values_to = "batch") %>%
  ggplot(aes(seurat_clusters, fill = batch)) + geom_bar(position = "fill")

pplot("batch")
pplot("batch_combined")

#filter by signature value
df <- sig_val(marker = "gene_list", use_func = "mean", filter = T)
use_id <- df %>% pivot_longer(cols = c(-id, -cluster), names_to = "signature", values_to = "score") %>%
  group_by(id) %>% mutate(rank = row_number(-score)) %>%
  filter(signature =="Mesenchyme", score!=0) %>% pull(id)

sub_data <-subset(data, cells = use_id)
DimPlot(sub_data)



hsc_marker <- FindAllMarkers(sub_data,min.pct = 0.15,  only.pos = T)
hsc_marker %>% group_by(cluster) %>% top_n(20, avg_logFC) %>% clipr::write_clip()


data <- sub_data

#chisq test for
data$seurat_clusters <- fct_drop(data$seurat_clusters)
data[[c("batch", "seurat_clusters")]] %>% table() #%>%
  chisq.test()

#combined batch factor

label %>% dplyr::select(sample, id = Cell.Barcode)

data$id <- colnames(data)

data@meta.data %>% filter(str_detect(batch, "^ma_set2"))  %>% pull(id)
  mutate(case = str_extract(id, "\\d{1,2}_\\d{1,2}$")) %>% pull(case) #%>% unique()
  sav(ma_set2)

  ma_set1$type <- ifelse(str_detect(ma_set1$id, "HCC"), "HCC", "ICC")
  ma_set2$type <- ifelse(str_detect(ma_set2$id, "HCC"), "HCC", "ICC")
  sa

  label <- label[c(1,2)]
  label$Cell.Barcode <- paste0(label$Cell.Barcode, "_6")

  any(data$id %in% label$Cell.Barcode)

  data$id <- data@meta.data %>% left_join(label, by = c("id" = "Cell.Barcode")) %>%
    mutate(id = if_else(batch =="ma_set1", paste0(Sample, "_", id), as.character(id))) %>%
    pull(id)



#
data$disease <- data$batch %>%
  fct_collapse(normal_1 = "macpoland",
               normal_2 = "aizarani",
               normal_3 = c("chandran_cd45nega_ht", "chandran_cd45posi_ht"),
               foetal = "segal",
               cirrhosis = c("chandran_cd45nega_ch", "chandran_cd45posi_ch"),
               carcinoma = c("ma_set1", "ma_set2"),
               pbc = c("pbc_case1", "pbc_case2"),
               blood = "ramachandran_blood")

data$disease <- data@meta.data %>% mutate(disease = case_when(str_detect(id, "^HCC") ~"HCC",
                                              str_detect(id, "^ICC") ~"ICC",
                                              TRUE ~ as.character(disease))) %>% pull(disease)


data$batch_combined <- data$batch %>%
  fct_collapse(Macpoland = "macpoland",
               Aizarani = "aizarani",
               Chandran = c("chandran_cd45nega_ht", "chandran_cd45posi_ht", "chandran_cd45nega_ch", "chandran_cd45posi_ch", "ramachandran_blood"),
               Segal = "segal",
               Ma = c("ma_set1", "ma_set2"),
               ours_case1 = "pbc_case1",
               ours_case2 = "pbc_case2")



#differential gene expression analysis
  batch_list <- unique(data$batch)
  marker<- tibble()
  for(i in seq_along(batch_list)){
    name <- substitute(batch_list[i])
    sub_temp <- subset(data, batch == batch_list[i])
    temp <- FindAllMarkers(object = sub_temp, min.pct = 0.15,  only.pos = T)
    temp$batch <- batch_list[i]
    marker <- marker %>% bind_rows(temp)
  }


  res_ref <- diff_test(batch_combined)

  saveRDS(marker, "batch_marker.rds")
  saveRDS(res_ref, "batch_reference_marker.rds")


  #search shared gene or whole gene whithin cluster

  res_ref %>% filter(p_val_adj < 0.05) %>% group_by(cluster, batch) %>%
    nest() %>% mutate(data = map(data, ~pull(., gene))) %>% group_by(cluster) %>% #dplyr::select(-batch) %>%
    nest() -> df_marker
  df_marker$data

  #show as a venn figure

  df_marker  %>% mutate(data = mIap2(data, cluster, ~zip(df = .x, name = .y) %>%
                                      VennDiagram::venn.diagram(filename = paste0("cluster_",.[2]), data = .[1], main = .[2])))

  #function converting df to list
    zip <- function(df, name) {
      list <- df$data
      names(list) <- df$batch
      res <- list(list, as.character(name))
      return(res)
    }










# hepatocyte and cholangiocyte --------------------------------------------
data_list <- list.files(pattern = ".rds")
data_list

make_subset_(data_list = data_list, cell_type = c("Hepatocyte", "Cholangiocyte"), save_folda = "hepato_cholangio")

setwd("hepato_cholangio/")
object_list <- make_list(cell_type = "Cholangiocyte")

k.filter <- min(200, min(sapply(object_list, ncol)))

data <- combined(object.list = object_list, cell_type = "Cholangiocyte", k.filter = k.filter)

data[[]]

data$reference <- data@meta.data %>% mutate(reference = fct_collapse(batch, Macparland = "macpoland",
               Aizarani = "aizarani",
               RamaChandran = c("chandran_cd45nega_ht", "chandran_cd45posi_ht", "chandran_cd45nega_ch", "chandran_cd45posi_ch", "ramachandran_blood"),
               ours_case1 = "pbc_case1",
               ours_case2 = "pbc_case2")) %>%
  mutate(reference = case_when(str_detect(batch, "fetal|adult")~"Segal",
                               str_detect(batch, "HCC|ICC")~"Ma",
                               TRUE~as.character(reference))) %>% pull(reference)

data$disease <- data@meta.data %>% mutate(disease = fct_collapse(batch,NL_1 = "macpoland",
               NL_2 = "aizarani",
               NL_3 = c("chandran_cd45nega_ht", "chandran_cd45posi_ht"),
               CH = c("chandran_cd45nega_ch", "chandran_cd45posi_ch"),
               PBC = c("pbc_case1", "pbc_case2"),
               BL = "ramachandran_blood")) %>%
  mutate(disease = case_when(str_detect(batch, "fetal")~"FL",
                             str_detect(batch, "adult")~"NL_4",
                             str_detect(batch, "HCC")~"HCC",
                             str_detect(batch, "ICC")~"ICC",
                             TRUE~as.character(disease))) %>% pull(disease)



up(group.by = "disease")
up(group.by = "reference")

saveRDS(data, file = "hepato_cholangio_combined.rds")

signature_plot()
ump("KRT7")


#filtering procedure

add_m <- function(df_list, add = "_m") {
  colnames(df_list) <- paste0(colnames(df_list), add)
  return(df_list)
}

df_gene_list_m <- sig_val() %>% add_m(add = "_m")
df_gene_list_gm <- sig_val(use_func = "gm_mean") %>% add_m(add = "_gm")
df_segal_list_m <- sig_val(marker = "segal_list") %>% add_m(add = "_m")
df_segal_list_gm <- sig_val(marker = "segal_list",use_func = "gm_mean") %>% add_m(add = "_gm")

df_com <- list(df_gene_list_m, df_gene_list_gm, df_segal_list_m, df_segal_list_gm) %>% reduce(cbind)
df_com <- df_com %>% keep(is.numeric)
data <- add_meta(df = df_com)
data <- add_info(data= data)
data[[]] %>% colnames

data <- sub(Hepatocyte_m >0&Cholangiocyte_m==0, seurat_clusters %in% c(1,3,4,5,11,14,15,16 ), !disease %in% c("BL", "HCC","ICC","unknown"))

ts(group.by = "reference")
ts(group.by = "disease")

sav()

bar_origin("disease")
bar_origin("reference")
tmap("KRT7")

Idents(data)
ts(group.by = "reference", pt.size = 2) + gghighlight(reference %>% str_detect("case"))
tmap("Hepatocyte_m")
feature_sig("zone_list")
sa(hepato_subset)

marker_hepato <- FindAllMarkers(object = data, min.pct = 0.25, min.diff.pct = 0.15, only.pos = T)
sav(marker_hepato)

marker %>% group_by(cluster) %>% top_n(20, avg_logFC)


marker_hepato %>% group_by(cluster) %>% nest()

DoHeatmap(data, features = top_10)



library(clusterProfiler)
library(org.Hs.eg.db)

data(geneList, package="DOSE")

gene <- names(geneList)[abs(geneList) > 2]

gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)

ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)


marker %>% filter(p_val_adj <0.05) %>% group_by(cluster) %>% arrange(cluster,-avg_logFC) %>% dplyr::select(cluster, gene) %>%
  nest() %>% mutate(data = map(data, ~pull(., gene)))

ts(data, group.by = "reference")  + gghighlight(reference == "Macparland")
ts(data, group.by = "reference")  + gghighlight(reference == "Our_case1")
ts(data, group.by = "reference")  + gghighlight(reference == "Our_case2")
ts(data, group.by = "disease")  + gghighlight(disease %in% c("HCC", "ICC", "FL"))
ts(data, group.by = "disease")  + gghighlight(disease %in% c("CH"))







# cholangio_hepato --------------------------------------------------------

up()
signature_plot_within()
bar_origin("disease")
marker_8 <- FindMarkers(data, ident.1 = 8, only.pos = T, min.diff.pct = 0.2, logfc.threshold = log(2))
marker_10 <- FindMarkers(data, ident.1 = 10, only.pos = T, min.diff.pct = 0.2, logfc.threshold = log(2))
marker_13 <- FindMarkers(data, ident.1 = 13, only.pos = T, min.diff.pct = 0.2, logfc.threshold = log(2))

marker_8$cluster <- 8
marker_8$gene <- rownames(marker_8)
marker_10$cluster <- 10
marker_10$gene <- rownames(marker_10)
marker_13$cluster <- 13
marker_13$gene <- rownames(marker_13)

marker_list <- list(marker, marker_8, marker_10, marker_13)
marker <- Reduce(x = marker_list, f = rbind)
marker <- marker %>% arrange(cluster)
saveRDS(marker, "hepato_cholangio_marker.rds")
marker %>% get_marker_table(p_val_adj<0.05)

#filter celll
data <- add_info(data)
data <- add_sig_val(data)
data <- sub_fil(data, !seurat_clusters %in% c(8, 10, 13), !disease %in% c("HCC", "ICC", "BL","unknown"))
use_id <- up() %>% CellSelector()
data <- sub_fil(data, id %in% use_id)
sa_data(hepato_chombined_filtered)




#hepato_cubset
data <- fil_cell("Hepatocyte", remove_cluster = c(0,2,6,7,9,12))
use_id <- up() %>% CellSelector()
data <- sub_fil(data, id %in% use_id)
up()
sa_data(hepato_subset)

combined()


#KRT7 cor

data@assays$RNA@data["KRT7", ] %>% hist(breaks = 200)
tmap(c("Cholangiocyte_m", "Hepatocyte_m"))
tmap(c("Cholangiocyte_gm", "Hepatocyte_gm"))

data_sub <- sub_fil(object = data, Cholangiocyte_gm>0)

use_id <- tmap(c("Cholangiocyte_m")) %>% CellSelector()
ts()
use_id <- up() %>% CellSelector()
data <- SetIdent(data, cells = use_id, value = "cholangio_like_cells")
up()
cholangio_like_marker <- FindMarkers(data, ident.1 = "cholangio_like_cells", min.pct = 0.2, only.pos = T)
cholangio_like_marker %>% View
defaultProto
tmap("KRT7")
tile("KRT7")
data_hibrid <- sub_fil(data, Cholangiocyte_gm>0, Hepatocyte_gm>0)
up(data_hibrid)
signature_plot("segal_list")
signature_plot("gene_list")



data2 <- sub_fil(data, Cholangiocyte_gm ==0)
tile(segal_list, object = data_hibrid)
use_id <- pick_id(seurat_clusters == "1", Cholangiocyte_gm >0)
data <-SetIdent(data, cells = use_id, value = "cholangio_like_hepato")
data$seurat_custers_2 <- data$seurat_clusters
data$seurat_clusters <- Idents(data)
tile(pick_gene("^MUC\\d{1,2}$"))


liver_marker_list_20_remove_dup <- remove_list_dup(liver_marker_list_20)
liver_marker_list_10_remove_dup <- remove_list_dup(liver_marker_list_10)
save_list(liver_marker_list_20_remove_dup)
save_list(liver_marker_list_10_remove_dup)

liver_marker <- liver_marker_list_10_remove_dup[names(liver_marker_list_10_remove_dup) %>% str_which("Hep|Chola")]
save_list(liver_marker)


tile(liver_marker)
tile_legend()

tmap(c("Cholangiocyte_gm", "Cholangiocyte_m", "Hepatocyte_m", "Hepatocyte_gm"))
ump(c("Cholangiocyte_gm", "Cholangiocyte_m", "Hepatocyte_m", "Hepatocyte_gm"))

krt7_cor$cor %>% hist(., main = "KRT7")

krt7_cor %>% filter(cor > 0.3) %>% nrow


#tile map of each cell
FetchData(data_sub, vars = unlist(segal_list)) %>% rownames_to_column("id") %>%
  mutate(id = row_number()) %>%
  pivot_longer(cols = -id, names_to = "gene", values_to = "log10CPM") %>%
  ggplot(aes(id, gene, fill = log10CPM)) + geom_tile() +
  scale_fill_gradientn(colours = c("red","yellow","white","lightblue","darkblue"),
                         values = c(1.0,0.7,0.6,0.4,0.3,0))

data$cholangio_m %>% summary
data$cholangio_m %>% mean()
data$Hepatocyte_m %>% mean()
FetchData(data, vars = "KRT7") %>% ggplot(aes(KRT7)) + geom_histogram(bins = 300)

FetchData(data, vars = "KRT7") %>% .$KRT7 %>% mean
FetchData(data, vars = "KRT7") %>% .$KRT7 %>%sd

vl(c("Cholangiocyte_gm", "Hepatocyte_gm"))
vl(c("KRT7"))



#before remove cholangiocyte_gm cells
df <- get_df(object = data)
krt7_cor <- do_cor(arg1 = df, arg2 = "krt7", gene = "KRT7")
krt7_cor %>% View

#after
df2 <- get_df(object = data2)
krt7_cor_2 <- do_cor(arg1 = df2, arg2 = "krt7_2", gene = "KRT7")
krt7_cor_2 %>% View
sav(krt7_cor)
sav(krt7_cor_2)

data$KRT7 <- if_else(df$KRT7 > 0.0843, "KRT7posi", "KRT7nega")
bar_origin("KRT7")
Idents(data) <- "KRT7"
up(data2)
ts(data2)
df2$KRT7 %>% mean

#only normal liver cell
data <- sub_fil(data, str_detect(disease, "NL"))
df <- get_df(data)

res_krt7 <- do_cor(arg1 = df, arg2 = "krt7_normal_cell", gene = "KRT7")


#
DefaultAssay(data) <- "RNA"
data <- sub_fil(data, Cholangiocyte_m ==0)



krt7_gene <- krt7_cor %>% filter(cor > 0.3) %>% pull(gene)
krt7_gene <- convert_gene(krt7_gene)
res <- ReactomePA::enrichPathway(krt7_gene$ENTREZID, pvalueCutoff = 0.5)
ump()

krt7_cor$cor %>% mean(na.rm = T)
krt7_cor$cor %>% sd(na.rm = T)

#remove

#pathway analysis
marker <- marker_list(marker)
marker <- do_geneano(marker)
marker <- marker %>% mutate(barplo = map2(res_enricher, cluster, bar))
marker <- marker %>% mutate(cnet = map2(res_enricher, cluster, ~cnet(.x, .y)))
anotation_marker <- marker
sav(anotation_marker)
sa_data(hepato_cholangio_combined_filter_3022)




marker_anotation <- marker %>% dplyr::select(cluster, gene_symbol) %>% unnest %>%
  left_join(a, by = c("gene_symbol" = "ID"))
marker_anotation %>% write_clip()

up(group.by = "disease") + gghighlight(str_detect(disease, "PBC"), label_key = T)

#cholangio vs hepatocyte
up()
data$seurat_clusters <- data$seurat_clusters2
Idents(data) <- "seurat_clusters"
up()
data$cell_type <- data@meta.data  %>% mutate(cell_type = fct_collapse(seurat_clusters, Hepatocyte = c("1","3","4","5","11","14","15","16"),
                                                   cholangiocyte = c("0","2","6","7","9","12"))) %>% pull(cell_type)

Idents(data) <- "cell_type"
up()

marker_hepato_vs_cholangio <- FindAllMarkers(data, only.pos = T)


#cholangio_vs cholangio_like_cells

DefaultAssay(data) <- "RNA"
cholangio_like_cell_id <- cell_id
cholangiocyte_id <- cell_id2

sav(cholangio_like_cell_id)
sav(cholangiocyte_id)

data <- SetIdent(data, cells = cell_id, value = "cholangiocyte_like_hep")
data <- SetIdent(data, cells = cell_id2, value = "cholangiocyte")
marker <- FindMarkers(data, ident.1 = "cholangiocyte_like_hep",ident.2 =  "cholangiocyte", min.pct = 0.2)

marker <-marker %>%  rownames_to_column(var = "gene") %>% mutate(cluster = if_else(avg_logFC>0, "cholangio_like_cell", "cholangiocyte"))
marker_cholangio_vs_cholangio_like <- marker
sav(marker_cholangio_vs_cholangio_like)
marker %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(50, avg_logFC) %>% 
  filter(cluster == "cholangio_like_cell") %>% pull(gene) %>% write_clip
marker %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(50, avg_logFC) %>% 
  filter(cluster == "cholangiocyte") %>% pull(gene) %>% write_clip

marker <- marker %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(50, avg_logFC)

marker_ano <- read_delim("cholangiocyte_marker_annotatio.txt", delim = "\t")
marker_ano2 <- read_delim("cholangio_like_cell_marker_annotation.txt", delim = "\t")

marker_ano %>% colnames
marker_ano2 %>% colnames
marker_ano2$COG_ONTOLOGY <- NULL
marker_com <- rbind(marker_ano, marker_ano2)

marker <- marker %>% left_join(marker_com, by = c("gene" = "ID"))
marker_cholangio_vs_cholangio_like_anotation <- marker
sav(marker_cholangio_vs_cholangio_like_anotation)

data$seurat_clusters2 <- data$seurat_clusters
data$seurat_clusters <- Idents(data)
signature_plot()





# macrophage -------------------------------------RNA-------------------------

data_list = list.files(pattern = ".rds")
data_list

make_subset(data_list, save_folda = "macrophage",cell_type = "MP", use_func = "mean")

object_list <- make_list("MP")

k.filter <- min(200, min(sapply(object_list, ncol)))

data <- combined(object.list = object_list,cell_type = "MP", k.filter = k.filter)

ts()
ts(group.by = "disease") + gghighlight(disease == "BL")
data <- add_info(data)
data <- add_sig_val(data)

mp_marker <- do_diff()
sav(mp_marker)
get_marker_table(mp_marker, p_val_adj < 0.05)

signature_plot_within()
data <- fil_cell(cell_type = "MP", remove_cluster = 6, remove_disease = "BL")

ts(group.by = "reference")

bar_origin("reference")
bar_cluster("disease")
bar_origin("reference")
bar_origin("disease")

marker <- do_diff(x = data)
sav(MP_marker)
data[[]]
marker_diff <- diff_test(x = "reference")
sav(marker_diff)
get_diff_test_marker(marker_diff)
make_venn(marker_diff)


# Tcell -------------------------------------------------------------------

data_list = list.files(pattern = ".rds")

load_list(data_list)

make_subset(data_list, save_folda = "Tcell",cell_type = "Tcell", use_func = "mean")

object_list <- make_list("Tcell")
sapply(object_list, ncol)
k.filter <- min(200, min(sapply(object_list, ncol)))

data <- combined(object.list = object_list,cell_type = "Tcell", k.filter = k.filter)

data <- add_info(data = data)
data[[]] %>% filter(is.na(batch))
bar_origin("reference")

ts()
data <- add_info(data)
data <- add_sig_val(data)

ts(group.by = "reference")
ts(group.by = "disease")
signature_plot()
signature_plot_within()
data <- fil_cell(cell_type = "Tcell", remove_cluster = 8, remove_disease = "BL")
bar_origin("reference")
bar_origin("disease")
bar_cluster("reference")
bar_cluster("disease")



Tcell_marker <-do_diff()
sav(Tcell_marker)
Tcell_marker_list <- marker_list(Tcell_marker)
Tcell_marker_list <- add_enrich_anotation(Tcell_marker_list)
sav(Tcell_marker_list)

ts(data, group.by = "disease", pt.size =2) + gghighlight(disease %>% str_detect("PBC"))
get_marker_table(Tcell_marker_list)

Tcell_diff_test_res <- diff_test("reference", object = sub_data)
diff_test_res <- diff_test(x ="reference", object = sub_data)

get_diff_test_marker(diff_test_res)

df <- make_venn(Tcell_diff_test_res)
df <- df %>% mutate(data = map(data, ~marker_list_map(.)))
df <- df %>%
  mutate(data = map(data, ~do_geneano(., use_func = geneano_enricher, res_name = "res_enricher", gene_type = gene_entrez)))

df <- df %>% mutate(data = map(data, ~mutate(.x, cluster = paste0("cluster_", cluster))))
df %>% mutate(data = map(data, ~mutate(., bar_plot = pmap(list(arg1 = res_enricher, arg2 = batch, arg3 = cluster), .f = bar))))
a %>% mutate(data = map(data, ~mutate(., bar_plot = pmap(list(arg1 = res_enricher, arg2 = batch, arg3 = cluster), .f = cnet))))

Tcell_diff_test_res %>% dplyr::select(cluster, data) %>% mutate(data = map(data, ~dplyr::select(.x,batch, gene_symbol))) %>%
  unnest() %>% mutate(gene_symbol = unlist(map(gene_symbol, ~paste0(.x, collapse = ", ")))) %>% write_clip()


# ILC ---------------------------------------------------------------------
data <- readRDS("~/single_cell/single_cell_project/data/Seurat_object/ILC/ILC_combined.rds")
ts()
signature_plot()
signature_plot_within()
data <- add_all(data)
data<- add_sig_val(data)
bar_origin("reference")
bar_origin("disease")

bar_cluster("reference")
data <- fil_cell("ILC", remove_cluster = "5", remove_disease = "BL")
sa_data(ILC_combined_filtered)

signature_plot_within()
bar_origin("disease")
bar_origin("reference")
bar_cluster("reference")
bar_cluster("disease")


tmap(c("FCGR3A", "GNLY", "CX3CR"))
tmap(c("CD3D", "PTPRC"), min.cutoff = 0)
ILC_marker <- do_diff()
ILC_marker_list <- marker_list(ILC_marker)
ILC_marker_list <- ILC_marker_list %>%
  do_geneano(use_func = geneano_enricher, res_name = "res_enricher", gene_type = gene_entrez)
ILC_marker_list <- ILC_marker_list %>% mutate(bar_plot = map2(res_enricher, cluster, ~bar(.x, .y)))
ILC_marker_list <- ILC_marker_list %>% mutate(cnet_plot = map2(res_enricher, cluster, ~cnet(.x, .y)))
sav(ILC_marker_list)
ILc_diff_test_res <- diff_test("reference")
sav(ILc_diff_test_res)

df <- make_venn(ILc_diff_test_res)
df <- df %>% mutate(data = map(data, ~marker_list_map(.)))

df <- df %>%
  mutate(data = map(data, ~do_geneano(., use_func = geneano_enricher, res_name = "res_enricher", gene_type = gene_entrez)))

df %>% mutate(data = map(data, ~mutate(bar_plot = map2(res_enricher, batch, ~bar(.x, .y)))))

Tcell_diff_test_res <- diff_test("reference")

Tcell_diff_test_res

# Bcell -------------------------------------------------------------------

data_list = list.files(pattern = ".rds")

load_list(data_list)

make_subset(data_list, save_folda = "Bcell",cell_type = "Bcell", use_func = "mean")

object_list <- make_list("Bcell")
object_list <- object_list[sapply(object_list, ncol)>50]
k.filter <- min(200, min(sapply(object_list, ncol)))

data <- combined(object.list = object_list,cell_type = "MP", k.filter = k.filter)
ts()

ts(group.by = "reference")
ts(group.by = "disease")
bar_origin
signature_plot()
signature_plot_within()

bcell_marker <- do_diff()
sav(bcell_marker)

bcell_marker %>% group_by(cluster) %>% top_n(20, avg_logFC) %>% write_clip()

data <- add_all(data)
data <- add_sig_val(data)

data <- fil_cell(cell_type = "Bcell", remove_cluster = c(4,9,10), remove_disease = "BL")
ts()
ts(group.by = "reference")
ts(group.by = "disease")
bar_origin("reference")
bar_origin("disease")
bar_cluster("reference")
bar_cluster("disease")
ts(group.by = "disease", pt.size =2) + gghighlight(disease %>% str_detect("PBC"), label_key = F)

tmap(bcell_list[[1]])
tmap(c("MZB1", "CD38", "CD27", "SDC1"), min.cutoff = 0)
sa_data(Bcell_combined_filtered)
bcell_marker <- do_diff()
bcell_marker_list <- marker_list(bcell_marker)
bcell_marker_list
bcell_marker_list <- bcell_marker_list %>%
  do_geneano(use_func = geneano_enricher, res_name = "res_enricher", gene_type = gene_entrez)
bcell_marker_list <- bcell_marker_list %>% mutate(bar_plot = map2(res_enricher, cluster, ~bar(.x, .y)))
bcell_marker_list <- bcell_marker_list %>% mutate(cnet_plot = map2(res_enricher, cluster, ~cnet(.x, .y)))
sav(bcell_marker_list)

# Mesenchyme --------------------------------------------------------------

data_list = list.files(pattern = ".rds")

load_list(data_list)

make_subset(data_list, save_folda = "HSC",cell_type = "Mesenchyme", use_func = "mean")

object_list <- make_list("Mesenchyme")

k.filter <- min(200, min(sapply(object_list, ncol)))

data <- combined(object.list = object_list,cell_type = "MP", k.filter = 25)
ts()

signature_plot()
signature_plot_within()
data <- add_all(data)
data <- fil_cell(cell_type = "Mesenchyme", remove_cluster = c(7,12), remove_disease =  c("BL","unknown"))
bar_origin("reference")
bar_origin("disease")
bar_cluster("reference")
bar_cluster("disease")
mesenchyme_marker <- do_diff()
mesenchyme_marker_list <- marker_list(mesenchyme_marker)
sav(mesenchyme_marker)
sa(Mesenchyme_combined_filtered)

tmap(c("PDGFRA", "PDGFRB","ACTA2A","VIM", "COL1A1","COL1A2"), min.cutoff= 0)
ts(group.by = "disease", pt.size =2) + gghighlight(disease %>% str_detect("PBC"))
ts(group.by = "disease", pt.size =2)

mesenchyme_marker_list <- mesenchyme_marker_list %>% do_geneano(use_func = geneano_msig_gsea, res_name = "res_gsea", gene_type = gene_list_entrez)
mesenchyme_marker_list <- mesenchyme_marker_list %>% do_geneano(use_func = geneano_enricher, res_name = "res_enricher", gene_type = gene_entrez)
mesenchyme_marker_list <- mesenchyme_marker_list %>% do_geneano(use_func = geneano_groupgo, res_name = "res_groupgo", gene_type = gene_entrez)

b <- b %>% mutate(bar = map2(res, cluster, ~bar(.x, .y)))
b <- b %>% mutate(bar2 = map2(res, cluster, ~bar(.x, .y)))
mesenchyme_marker_list <- mesenchyme_marker_list %>% mutate(cnet = map2(res_enricher, cluster, ~cnet(.x, .y)))


# nkt ---------------------------------------------------------------------

ts()
ts(group.by = "reference", pt.size = 2) + gghighlight(str_detect(reference, "ours"))
data <- add_info(data)
bar_origin("reference")
bar_origin("disease")
signature_plot_within()
ts(group.by = "disease") + gghighlight(str_detect(disease, "BL"))
marker <- do_diff()
marker_19 <- FindMarkers(data, ident.1 = "19", only.pos = T, min.diff.pct = 0.20, logfc.threshold = log(2))
marker_15 <- FindMarkers(data, ident.1 = "15", only.pos = T, min.diff.pct = 0.20, logfc.threshold = log(2))
marker_16 <- FindMarkers(data, ident.1 = "16", only.pos = T, min.diff.pct = 0.20, logfc.threshold = log(2))
sav(marker_19)
sav(marker_15)
sav(marker_16)

get_marker_table(nkt_marker, p_val_adj < 0.05)



# endothelia --------------------------------------------------------------
data <- readRDS("~/single_cell/single_cell_project/data/Seurat_object/Endothelia/Endothelia_combined.rds")
ts()
signature_plot()
signature_plot_within()
data <- add_all(data)
data
data <- fil_cell("Endothelia",remove_cluster = 9, remove_disease = "BL")
data
ts()
sa(Endothelia_combined_filtered)
bar_origin("reference")
bar_origin("disease")
bar_cluster("reference")
bar_cluster("disease")
up()
tmap(unlist(endothelial_list), min.cutoffc= 0)
endothelia_marker <- do_diff()
sav(endothelia_marker)

FeaturePlot()


# batch  -----------------------------------------------------------

