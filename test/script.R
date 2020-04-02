library(package2)
library(tidyverse)
library(purrr)
library(Seurat)
library(clipr)
library(gghighlight)


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

#filter celll
data <- add_info(data)
data <- add_sig_val(data)
data <- sub_fil(data, !seurat_clusters %in% c(8, 10, 13), !disease %in% c("HCC", "ICC", "BL","unknown"))
use_id <- up() %>% CellSelector()
data <- sub_fil(data, id %in% use_id)
sa_data(hepato_chombined_filtered)


#HCC ICC contain
data <- add_info(data)
data <- sub_fil(data, !seurat_clusters %in% c(8, 10, 13), !disease %in% c("BL","unknown"))
data <- sub_fil(data, !seurat_clusters %in% c(1,3,4,5, 11,14, 15,16))
data <- sub_fil(data, seurat_clusters != 6)


use_id <- up() %>% CellSelector()

data <-sub_fil(data, id %in% use_id)
sa_data(cholangio_ICC)

data <- SetIdent(data, cells = use_id, value = "add_1")
data$seurat_clusters2 <- data$seurat_clusters
data$seurat_clusters <- Idents(data)

upd("ICC")
cholangio_marker <- diff_test_batch(disease)

cholangio_marker
cholangio_marker %>% filter(batch == "ICC") %>% write_clip
data_icc<- sub_fil(data, disease == "ICC")
up(data_icc)
tile(gene = c("SFRP5", "ASGR1", "CFTR", "CD24", "NCAM1", "NES", "KRT19", "KRT7", "EPCAM","TACSTD2", "MUC6", "S100P", "CXCL8", "MUC13", "MUC17"),
     object = data)


#hepato_cubset
data <- fil_cell("Hepatocyte", remove_cluster = c(0,2,6,7,9,12))
use_id <- up() %>% CellSelector()
data <- sub_fil(data, id %in% use_id)
up()
sa_data(hepato_subset)


#compare with mouse cells

mouse_zonal_gene$Gene
sav(mouse_zonal_gene)

zonal_gene <- Seurat::AverageExpression(data, assays = "RNA", features = mouse_zonal_gene$Gene)

m <- apply(zonal_gene$RNA %>% as.data.frame(), MARGIN = 1, FUN = mean)
s <- apply(zonal_gene$RNA %>% as.data.frame(), MARGIN = 1, FUN = sd)

zonal_gene$m <- m
zonal_gene$s <- s

zonal_gene <- zonal_gene %>% mutate_all( .funs = ~(.-m)/s ) %>% dplyr::select(1:8) %>% add_column(gene = mouse_zonal_gene$Gene)
colnames(zonal_gene) <- colnames(zonal_gene) %>% paste(., "cluster", sep = "_")
res_cor <- zonal_gene %>% cbind(mouse_zonal_gene) %>% keep(is.numeric) %>% cor(method = "spearman")
res_cor[,1:8] %>% reshape2::melt() %>%
  filter(str_detect(Var1, "L")) %>% ggplot(aes(Var1, fct_reorder(Var2, value), fill = value)) +
           geom_tile() + scale_fill_gradient2(low = "blue", mid = "white", high = "red" )

zonal_gene

hemato_marker <- do_diff_test(data)
sav(hemato_marker)

hemato_marker <- hemato_marker %>% marker_list()
hemato_marker <- hemato_marker %>% do_geneano(use_func = geneano_enricher, res_name = "res_enricher", gene_type = gene_entrez)
hemato_marker <- hemato_marker %>% mutate(bar_plot = map2(res_enricher, cluster, ~bar(.x, .y)))
hemato_marker <- hemato_marker %>% mutate(cnet_plot = map2(res_enricher, cluster, ~cnet(.x, .y)))

hemato_marker <- hemato_marker %>% mutate(res_enricher_tab = map(res_enricher,~.@result)) %>%
  mutate(res_enricher_tab = map2(res_enricher_tab,cluster, ~mutate(.x, cluster = paste0("cluster_",.y))))

res_tab <- hemato_marker$res_enricher_tab %>% purrr::reduce(rbind)

bar_origin(disease, seurat_clusters)
bar_origin(seurat_clusters, disease)
data$disease %>% table

#cholangio_subset
data <- fil_cell("Cholangiocyte", remove_cluster = c(1,3,4,5, 11,14, 15,16))

use_id <- up() %>% CellSelector()

data <- SetIdent(data, cells = use_id, value = "cluster_far")

use_id <- up() %>% CellSelector()
data <- SetIdent(data, cells = use_id, value = "new_1")
up()
cholangio_marker <- FindAllMarkers(data, only.pos = T, logfc.threshold = 0.15)

cholangio_marker_fil <- cholangio_marker %>% group_by(cluster) %>% filter(p_val_adj<0.05)

cholangio_marker_anotatioin <- read_delim("cholangio_marker_annotation.txt", delim = "\t")

cholangio_marker_fil <- cholangio_marker_fil %>% left_join(cholangio_marker_anotatioin, by = c("gene" = "ID"))

sav(cholangio_marker_fil)

cholangio_marker_fil %>% write_clip()

package2::tile(c("SFRP5", "ASGR1", "CFTR", "CD24", "NCAM1", "NES", "KRT19", "KRT7", "EPCAM","TACSTD2", "MUC6", "S100P", "CXCL8", "MUC13", "MUC17"))

cholangio_marker <- cholangio_marker %>% marker_list
cholangio_marker <- cholangio_marker %>% do_geneano(.)
cholangio_marker <- cholangio_marker %>% mutate(barplo = map2(res_enricher, cluster, bar))
cholangio_marker <- cholangio_marker %>% mutate(cnet = map2(res_enricher, cluster, ~cnet(.x, .y)))

cholangio_marker <- cholangio_marker %>% mutate(barplo = map2(res_enricher, cluster, bar))

up(group.by = "seurat_clusters")
data$seurat_clusters





id_ch("seurat_clusters")
use_id <- up() %>% CellSelector()
data <- sub_fil(data, id %in% use_id)
data$seurat_clusters <- data$seurat_clusters %>% fct_recode(add_1 = "new_1")
up()
bar_origin(disease, seurat_clusters)
bar_origin(seurat_clusters, disease)
wrt_c_tab(data)

#vimentin analysis
package2::tile("VIM")
df <- get_df(data)
res_cor_vim <- do_cor(df, "vimentin", gene = "VIM")
res_cor_vim$cor %>% hist
res_cor_vim[1:100,] %>% .$gene %>% write_clip()
data@meta.data %>% colnames() %>% str_subset("_m$") %>% vl
sav(res_cor_vim)

res_cor_vim_ano <- read_delim("res_cor_vimentin_100_gene_ano.txt", delim = "\t")
res_cor_vim_ano_pathway <- read_delim("res_cor_vimentin_100_gene_pathway.txt", delim = "\t")
res_cor_vim_ano <- res_cor_vim[1:100,] %>% left_join(res_cor_vim_ano, by = c("gene" = "ID"))
sav(res_cor_vim_ano)
res_cor_vim_ano %>% write_clip()
res_cor_vim_ano_pathway %>% write_clip()
res_cor_vim %>% filter(str_detect(gene , "^FN1"))

vim_gene <- res_cor_vim_ano$gene %>% convert_gene(.)
res_enrich <- geneano_enricher(gene_entrez = vim_gene$ENTREZID)
res_enrich %>% barplot(showCategory = 30)
res_enrich %>% clusterProfiler::cnetplot(showCategory = 30)

geneano_enrichgo()
res_pathway <- read_delim("res_cor_vimentin_100_gene_pathway.txt", delim = "\t")

data <- add_meta_bi("VIM")
id_ch("VIM_bin")
vim_vs <- do_diff()
sav(vim_vs)
vim_vs %>% filter(pct.1-pct.2 > 0.2)
data$VIM %>% table
id_ch("seurat_clusters")
vl("VIM")

rna()
vl(c("KRT19", "CDH1", "TGFB1", "TGFB1R1", "ACTA2", "SNAI1", "SNAI2")) + scale_colour_manual(values = c("blue", "red"))


bar_origin(seurat_clusters, VIM) + scale_fill_discrete(c(red = "strong", blue = "ordinary"))
bar_origin(disease, VIM_bin)

vl("VIM", group.by = "disease")

res_cor_vim_ano %>% head(50) %>% ggplot(aes(fct_reorder(gene, cor),cor, fill = gene)) +
  geom_bar(stat = "identity") + coord_flip() + guides(fill = F)

#remove hepatocyte_like_cells
data <- sub_fil(data, !seurat_clusters %in% c("cluster_far", "6"))
sa_data(cholangio_subset_fil)





 #monocle
mono <- make_monocle3(data)
mono <- do_monocle(mono)
mono_marker <- do_diff(mono, group_cells_by = "seurat_clusters")
mop(mono, color_cells_by = "seurat_clusters")
mono_marker
mono_marker %>% arrange(cell_group, -mean_expression)
mono[]


pr_graph_test_res <- monocle3::graph_test(mono, neighbor_graph="knn", cores=8)

pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

gene_module_df <- monocle3::find_gene_modules(mono[pr_deg_ids,], resolution=1e-2)



#new_cell_id
new_id_1 <- up() %>% CellSelector()
data <- SetIdent(data, cells = new_id_1, value = "new_cell_1")
sav(new_id_1)




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
krt7_cor_2%>% head(50) %>%  ggplot(aes(fct_reorder(gene, cor), cor, fill = gene)) + geom_bar(stat = "identity") +
  guides(fill = F) + coord_flip()
krt7_cor %>% ggplot(aes(cor)) + geom_histogram(binwidth = 0.01)

#after removing cholangiocyte_gm positive cell
df2 <- get_df(object = data2)
krt7_cor_2 <- do_cor(arg1 = df2, arg2 = "krt7_2", gene = "KRT7")
krt7_cor_2 %>% View
sav(krt7_cor)
sav(krt7_cor_2)



#after removing cholangio_m positive cell
df3 <- get_df(data)
krt7_cor_3 <- do_cor(arg1 = df3, arg2 = "krt7_3", gene = "KRT7")
sav(krt7_cor_3)
krt7_cor_3%>% head(50) %>%  ggplot(aes(fct_reorder(gene, cor), cor, fill = gene)) + geom_bar(stat = "identity") +
  guides(fill = F) + coord_flip()
krt7_cor %>% ggplot(aes(cor)) + geom_histogram(binwidth = 0.01)

#compare each cor genes
krt7_gene <- krt7_cor %>% filter(cor>0.2) %>% .$gene
krt7_gene_3 <- krt7_cor_3 %>% filter(cor>0.2) %>% .$gene
sh_gene <- intersect(krt7_gene, krt7_gene_3)
krt7_gene <- convert_gene(krt7_gene)
krt7_gene_3 <- convert_gene(krt7_gene_3)
res_krt7_cor <- geneano_enricher(krt7_gene$ENTREZID)
res_krt7_cor_3 <- geneano_enricher(krt7_gene_3$ENTREZID, pvalueCutoff = 0.3, qvalueCutoff = 0.5)
krt7_cor_3$gene %>% write_clip()



data <- sub_fil(data, id %in% hepato_id)
up()

#cholangio and hepato_id
hepato_id <- pick_id(Hepatocyte_m -Cholangiocyte_m>0)
pick_id(Hepatocyte_m -Cholangiocyte_m  < 1) -> cholangio_id
data$cholangio_vs_hepato <- if_else(data$id %in% cholangio_id, "cholangio_like", "hepato_like")
id_ch("cholangio_vs_hepato")
hepato_vs_cholangio <- do_diff_test()
sav(hepato_vs_cholangio)
hepato_vs_cholangio %>% filter(p_val_adj <0.05) %>% write_clip()

data$hepato_vs_cholangio <- data@meta.data %>% mutate(hepato_vs_cholangio = case_when(id %in% cholangio_id~"cholangio_like",
                                                          id %in% hepato_id~"hepato",
                                                          TRUE~"cholangio")) %>% pull(hepato_vs_cholangio)
id_ch("hepato_vs_cholangio")
up()
tile()

hepa_vs_cholke_vs_cholangio <- do_diff_test()
hepa_vs_cholke_vs_cholangio %>% filter(p_val_adj <0.05) %>% write_clip()


marker_hepato_vs_cholangio %>% group_by(cluster) %>% top_n(30, avg_logFC) %>% dplyr::select(cluster, gene) %>%
  split(.$cluster) %>% map(~pull(., gene)) -> hepato_cholangio_list
sav(hepato_cholangio_list)
tile(hepato_cholangio_list, cluster_label = "hepato_vs_cholangio")

data$KRT7 <- if_else(df$KRT7 > 0.0843, "KRT7posi", "KRT7nega")
bar_origin("KRT7")
Idents(data) <- "KRT7"
up(data2)
ts(data2)
df2$KRT7 %>% mean





up()
data <- sub_fil(data, Cholangiocyte_m ==0)
vl("Cholangiocyte_m")
vl("Cholangiocyte_m", group.by = "disease")
vl("KRT7", group.by = "disease")
ump(c("Hepatocyte_m","Cholangiocyte_m"))

#only normal liver cell
data <- sub_fil(data, str_detect(disease, "NL"))
df <- get_df(data)

res_krt7 <- do_cor(arg1 = df, arg2 = "krt7_normal_cell", gene = "KRT7")


DefaultAssay(data) <- "RNA"
data <- sub_fil(data, Cholangiocyte_m ==0)



krt7_gene <- krt7_cor %>% filter(cor > 0.3) %>% pull(gene)
krt7_gene <- convert_gene(krt7_gene)
res <- ReactomePA::enrichPathway(krt7_gene$ENTREZID, pvalueCutoff = 0.5)
ump()





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

data$seurat_clusters2 <- data$seurat_clusters

Idents(data) <- "seurat_clusters"
up()
data$cell_type <- data@meta.data  %>% mutate(cell_type = fct_collapse(seurat_clusters, Hepatocyte = c("1","3","4","5","11","14","15","16"),
                                                   cholangiocyte = c("0","2","6","7","9","12"))) %>% pull(cell_type)

Idents(data) <- "cell_type"
up()

marker_hepato_vs_cholangio <- FindAllMarkers(data, only.pos = T)

sav(marker_hepato_vs_cholangio)
marker_hepato_vs_cholangio %>% group_by(cluster) %>% filter(p_val_adj <0.05) %>%
  filter(cluster == "Hepatocyte") %>% pull(gene) %>% write_clip()
signature_plot()

marker_hepato_vs_cholangio_all <- marker_hepato_vs_cholangio %>% group_by(cluster) %>%
  dplyr::select(cluster, gene)

marker_hepato_vs_cholangio_all_list <- marker_hepato_vs_cholangio_top20 %>% nest() %>%
  mutate(data = map(data, ~pull(.))) %>% pull(data)
names(marker_hepato_vs_cholangio_all_list) <- c("Cholangiocyte", "Hepatocyte")
marker_hepato_vs_cholangio_all_list <- marker_hepato_vs_cholangio_top20_list
marker_hepato_vs_cholangio_all_list %>% map(~head(., 20)) %>% tile()
save_list(marker_hepato_vs_cholangio_all_list)
get_list_name()

tile(marker_hepato_vs_cholangio_top20)

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






# MP ----------------------------------------------------------------------


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
#data making process
data_list = list.files(pattern = ".rds")

load_list(data_list)

make_subset(data_list, save_folda = "Tcell",cell_type = "Tcell", use_func = "mean")

object_list <- make_list("Tcell")
sapply(object_list, ncol)
k.filter <- min(200, min(sapply(object_list, ncol)))

data <- combined(object.list = object_list,cell_type = "Tcell", k.filter = k.filter)





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

data <- add_info(data)
add_sig_val(data, marker_list = "gene_list", use_func = "gm_mean")
tile(gene = gene_list, object = data)
up()
signature_plot()
signature_plot_within()
bar_origin(seurat_clusters, disease)

tile(Epithelial_list)

data <- fil_cell("Mesenchyme", remove_cluster = c(6,7,11,12), remove_disease = "BL")
sub_fil(object = data, disease != "HCC|ICC") %>% up

sa_data(Mesencyme_combined_filtered)

mesenchyme_marker %>% write_clip()
signature_plot_within()
bar_origin(seurat_clusters, disease)
annotation_table  <- read_delim(file = "mesenchyme_marker_annotation.txt", delim = "\t")

data$vs_ICC <- data@meta.data %>% mutate(vs_ICC = case_when(disease %in%c("HCC", "ICC")~"cancer",
                                                            TRUE~"non_cancer") ) %>% pull(vs_ICC)

vl(c("CD276","CD274", "PDCD1LG2", "NES" ), group.by = "vs_ICC")
marker <- do_diff()

marker <- marker %>% marker_list()
marker <- marker %>% do_geneano(use_func = geneano_enricher, res_name = "res_enricher", gene_type = gene_entrez)
marker <- marker %>% mutate(barplo = map2(res_enricher, cluster, ~bar(.x, .y)))
marker <- marker %>% mutate(cnet = map2(res_enricher, cluster, ~cnet(.x, .y)))

marker %>% dplyr::select(cluster, data) %>% unnest %>% left_join(annotation_table, by = c("gene" = "ID")) %>% write_clip

rna()


ump(unlist(mc_list))
ump(c("RBP4", "TGFB1"))
id <- up() %>% CellSelector()

data$two <-  if_else(data$id %in% id, "fbl_1", "fbl_2")
id_ch("two")
up(group.by = "disease")
upd(label = "PBC_2")

marker_two <- do_diff()
marker_two

diff_test()


ump(c("CCL2" ))

id_ch("seurat_clusters")

marker_2_vs_10 <- FindMarkers(data, ident.1 = 2, ident.2 = 10, only.pos = T)
marker_2_vs_10 <- marker_2_vs_10 %>% mutate(gene = rownames(marker_2_vs_10), cluster = if_else(pct.1-pct.2>0, "2", "10"))
sav(marker_2_vs_10)
marker_2_vs_10 %>% write_clip()

bar_origin(disease, seurat_clusters)
bar_origin(seurat_clusters, disease)
ump("S100A6")


# NKT ---------------------------------------------------------------------

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

up()
signature_plot()
signature_plot_within()


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




