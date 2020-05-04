library(Seurat)
library(package2)
library(tidyverse)
library(clipr)
library(gghighlight)


bar_origin(disease, label2)

data_n <- package2::sub_fil(data, !disease %in% c("HCC", "ICC"))

#marker_list

segal_list3 <-  list(Hepatocyte =c("AHSG", "ALB", "ASGR1", "FN1", "HNF4A", "HP","ORM1", "RBP4"),
                     Imature_Hepatocyte = c("CYP3A7", "GPC3", "AFP", "DLK1"),
                     BEC_common = c("KRT19", "KRT7", "SPP1", "CD133", "CD24", "SOX9", "STAT1", "EPCAM", "FGFR2"),
                     Mature_BEC = c("AGR2", "CLDN4", "CLDN10", "CXCL2", "LGALS2", "MMP7", "MUC1", "MUC5B", "NQO1", "ONECUT2", "TACSTD2", "TFF1", "TFF2", "TFF3", "TSPAN8"),
                     Imature_BEC = c("CAV1", "CDH6", "CLDN6", "CTNND2", "GPRC5B", "MCAM", "NCAM1", "SFRP5", "STAT4"),
                     Cycling = c("MKI67", "CCNA2", "CCNB2", "STMN1"))

save_list(segal_list3)

  list(hepatic = c("FGL1", "APOE", "APOH", "TF", "ORM1", "AHSG","APOA2", "APOC3"),

     Mature_BEC =c("MUC5B", "MUC3A", "ONECUT2", "TMEM45B","MUC1", "CD44", "TACSTD2", "CLDN4", "KRT7", "EPCAM", "CD24", "KRT19", "SPP1"),
     c("CYP3A7", "GPC3", "AFP", "DLK1"),
     c("NCAM1", "MUC6", "FGFR2", "GPRC5B", "SLC12A2", "PROM1", "SOX9", "CTNND2", "MCAM", "CDH6", "CXCL2", "ABCB9"))






#label4 is default
sa_data(hepato_cholangio_combined_filtered)
rna()

#cell number
data$disease %>% table ->a
b <- a %>% as.list %>% enframe %>% mutate(value = unlist(value), value = paste0(name, ": ", value)) %>% pull(value)
b %>% reduce(~paste(.x,.y, sep = ", ")) %>% write_clip()


data <- add_info(data)
#labeling
use_id <- up() %>% CellSelector()
prolife_id <- use_id; sav(prolife_id)
data <-SetIdent(data, use_id,value = 17)
data$seurat_clusters <- Idents(data)
data$disease <- data[[]] %>% mutate(disease = case_when(batch == "pbc_case1"~"PBC_1",
                                        batch == "pbc_case2"~"PBC_2",
                                        TRUE~disease)) %>% pull(disease)
data$disease <- fct_relevel(data$disease, c("PBC_1", "PBC_2", paste0("NL_", 1:4), "FL", "CH", "HCC", "ICC"))
data$label <- data$seurat_clusters %>% fct_collapse(Hep = c("1","3","4","5","11","14", "16"), im_Hep = "15",
                                                    HHyP = c("0","2","7","12"), BEC = c("9"),  Prolifer= "17" ) %>%
  fct_relevel(c("Hep", "im_Hep", "BEC", "HHyP", "Prolifer"))

data$seurat_clusters <- data$seurat_clusters %>% fct_reorder(as.numeric(data$label)) %>% fct_drop()
data$seurat_clusters
data$label2 <- data[[]] %>% mutate(n = plyr::mapvalues(seurat_clusters, from = levels(data$seurat_clusters),
                                                       to = c("", "", "", "", "", "", "","","", "(1)", "(2)", "(3)", "(4)", "")),
                                   name = paste0(label,n)) %>% mutate(name = fct_reorder(name, as.numeric(label))) %>%
  pull(name)
data$label3 <- paste0(data$label2, "_", data$disease)
data$label3 <- data[[]] %>% mutate(label3 = fct_reorder(label3, as.numeric(disease))) %>% pull(label3)
data$label3 <- data[[]] %>% mutate(label3 = fct_reorder(label3, as.numeric(label2))) %>% pull(label3)

id_ch("label2")
data <- sub_fil(data, !disease %in% c("HCC", "ICC"))

signature_tile("segal_list")
tile_plot("segal_list")

#remove low expresion cells
data <- sub_fil(data, !seurat_clusters %in% c("6", "8", "10","13"))

data <- sub_fil(data, !disease %in% c("BL","unknown"))


#remove dispersed cells
use_id2 <- up() %>% CellSelector()
remove_id <- use_id2
sav(remove_id)
data <- sub_fil(data, id %in% use_id2)

#zonation analysis
sub_nl <-sub_fil(data, str_detect(label2, "Hep"), !disease %in% c("CH", "HCC", "ICC", "FL"))

add_sig_val(sub_nl, "zone_list")
vl(object = sub_nl,c("zone1_", "zone2_"))
signature_tile("zone_list", object = sub_nl)
sub_nl$label4 <- paste0(sub_nl$label2, "_", sub_nl$disease) %>% fct_relevel(sort(unique(sub_nl$label4)))

id_ch("label4", object = sub_nl)


#meke clustring map
df <- sig_val(object = data, "segal_list")
df2 <- sig_val2(df)
df2 %>% pivot_wider(id_cols = c(cluster), names_from = signature, values_from = mean) -> df3
df3 %>% dist() %>% hclust() %>% plot

#note  contain cancer cells
sa_data(hepato_cholangio_combined_filtered)
id_ch("label2")
up()

#data filtering
data <- readRDS("~/single_cell/single_cell_project/data/Seurat_object/hepato_cholangio/hepato_cholangio_combined_filtered.rds")
data <- sub_fil(data, str_detect(label2, "Hep")&cholangio_gm==0|!str_detect(label2, "Hep")&cholangio_gm!=0)

data_hep <- sub_fil(data, str_detect(label2, "Hep"), !disease %in% c("HCC", "ICC"))

data_cho <- sub_fil(data, !str_detect(label2, "Hep"), !disease %in% c("HCC", "ICC"))


#bar_origin

bar_origin(disease, label2, data_hep)
bar_origin(label2, disease, data_hep)

bar_origin(disease, label2, data_cho)
bar_origin(label2, disease, data_cho)

bar_origin(disease, label2, data_)
rna(data_hep)

VlnPlot(data_hep, c("KRT7", "KRT19"), group.by = "disease")


FetchData(data_hep, c("KRT7", "KRT19"), slot = "data")


#Monocle


mono <- make_monocle3(data_cho)
mono <- do_monocle(mono)
mop(mono, color_cells_by = "label2")

mono_marker <- do_diff_mono(mono, group_cells_by = "seurat_clusters")
mop(mono, color_cells_by = "seurat_clusters")
mop(mono, color_cells_by = "label2")

mono_marker %>% arrange(cell_group, -mean_expression)
use_id1 <- mono %>% monocle3::choose_cells()
use_id2 <- mono %>% monocle3::choose_cells()



#correlation matrix
VlnPlot(data_hep, features = c("zone1_", "zone3_", "zone_ratio"))

up(data_hep)
id_ch("label2", object = data_hep)
av_df_hep <- AverageExpression(data_hep, assays = "RNA", features = data_hep@assays$RNA@var.features)
av_df_hep_batch <- batch_mat(av_df_hep[[1]], object = data_hep)
cor(av_df_hep_batch)
batch_cor_heatmap(av_df_hep_batch, method = "spearman")


VlnPlot(data_cho, features = "cholangio_gm", group.by = "disease")
data_cho <- sub_fil(data_cho, cholangio_gm !=0)
up(data_cho)
data_cho <- sub_fil(data, !str_detect(label2, "Hep"))
id_ch("label2", object = data_cho)
av_df_cho <- AverageExpression(data_cho, assays = "RNA", features = data_cho@assays$integrated@var.features)
av_df_cho_batch <- batch_mat(av_df_cho[[1]], object = data_cho)
batch_cor_heatmap(av_df_cho_batch, method = "spearman")

#sample

bar_origin(disease, label2)

#pathway plot
 #hepatocyte
  signature_tile("zone_list", object = data_hep)
package2::tile("zone_list",object = data_hep)
data_hep@assays$RNA
DefaultAssay(data_hep) <- "RNA"
marker_hep <- diff_test(data_hep)
sav(marker_hep)
marker_hep_rna %>% write_clip()
marker_hep <- marker_hep %>% left_join(marker_hep_ano, by = c("gene" = "ID"))
marker_hep %>% write_clip()
marker_hep_ano <- read_delim("marker_hep_anotation.txt", delim = "\t")
marker_hep <- marker_hep %>% diff_marker_convert()
marker_hep <- marker_hep %>% do_geneano()
marker_hep %>% mutate(bar = map2(res_enricher, cluster, ~barplot(.x, title = .y)))
marker_hep %>% mutate(a = map2(res_enricher, cluster,~.x@result %>% mutate(cluster = .y, n  = row_number()) %>% filter(p.adjust<0.05))) %>% pull(a) %>%
  purrr::reduce(rbind) %>% write_clip()
use_df <- read_clip_tbl() %>% filter(mark ==1)
use_df2 <- read_clip_tbl() %>% filter(mark ==1)
use_df3 <- rbind(use_df, use_df2)
use_colum <- use_df3 %>% group_by(cluster) %>% nest %>% pull(data)
marker_hep$res_enricher_fil <- use_colum
marker_hep <- marker_hep %>% mutate(res_enricher_fil = map2(res_enricher, res_enricher_fil, ~{a <- .x; a@result <- .y;return(a)}))
marker_hep <- marker_hep %>% mutate(bar = map2(res_enricher_fil, cluster, ~enrichplot:::barplot.enrichResult(.x, title = .y,font.size = 20)))
marker_hep$bar

 #cholangiocyte
marker_cho2 <- diff_test(data_cho)

saveRDS(object = data_cho, file = "data_cho.rds")
sav(marker_cho2)

marker_cho2 <- marker_cho2 %>% diff_marker_convert()
marker_cho2 <- marker_cho2 %>% do_geneano()
marker_cho2 <- marker_cho2 %>% mutate(bar = map2(res_enricher, cluster, ~barplot(.x, title = .y)))
marker_cho2 %>% mutate(a = map2(res_enricher, cluster,~.x@result %>% mutate(cluster = .y, n  = row_number()) %>% filter(p.adjust<0.05))) %>% pull(a) %>%
  purrr::reduce(rbind) %>% write_clip()

marker_cho_tabl <- read_clip_tbl()
marker_cho_tabl <- marker_cho_tabl %>% filter(X ==1)

b <- marker_cho_tabl %>% group_by(cluster) %>% nest %>% pull(data)
a <- marker_cho2
a$data2 <- b

a <- a %>% mutate(res_enricher_fil = map2(res_enricher, data2, ~{a <- .x; a@result <- .y; return(a)}))
a <- a %>% mutate(bar = map2(res_enricher_fil, cluster, ~barplot(.x, title = .y)))
marker_cho2_res_enrich <- a
sav(marker_cho2_res_enrich)


#KRT7 analaysis
VlnPlot(data_hep, features = "KRT7")
df_hep <- get_df(data_hep)
res_cor_krt7 <- do_cor(expr_df = df_hep, group_label = "krt", gene = "KRT7", method = "pearson")
res_cor_krt7 %>% head(30) %>%  ggplot(aes(fct_reorder(gene,cor), cor, fill = c("#2CD2DB"))) + geom_bar(stat = "identity") +
  coord_flip() + guides(fill = F)+ xlab("correlation")
use_gene <- res_cor_krt7 %>% head(100) %>% .$gene
df_hep_sub <- df_hep[c("KRT7", use_gene)]

df_hep %>% ggplot(aes(KRT7, FXYD2)) + geom_point()
df_hep %>% ggplot(aes(KRT7, TACSTD2)) + geom_point()
sav(res_cor_krt7)
add_meta_binval(object = data_hep ,gene = "KRT7")
id_ch("KRT7_bin", data_hep)

up(data_hep)

 #extract krt7 strong positive cells within hepatocyte
use_id <- data_hep[[]] %>% filter(KRT7_bin == "strong") %>% pull(id)
data_hep <- SetIdent(data_hep, use_id, value = "krt7Hep")
data_hep$label3 <- Idents(data_hep)
up(data_hep)
id_ch("label3")

data <- SetIdent(data, cells = use_id, value = "KRT7posi_Hep")
data$label3 <- Idents(data)
id_ch("label3")
up()

tile("segal_list2")
signature_tile("segal_list2")
pick_id(object = data_hep, str_detect(label2, "Hep")&KRT7_bin =="strong") ->use_id
marker_krt7_within_all <- diff_test(data)
sav(marker_krt7_within_all)
id_ch("label", object = data_hep)
data_hep <- SetIdent(data_hep, cells = use_id, value = "KRT7posi_Hep")
marker_krt7_within_all %>% write_clip()
up(data_hep)
marker_krt7_within_hep <- diff_test(data_hep)
sav(marker_krt7_within_hep)
marker_krt7_within_hep %>% write_clip()

marker_krt7_within_all %>%group_by(cluster) %>% top_n(5, wt = avg_logFC) %>%
  dplyr::select(cluster,gene) %>% nest() %>% mutate(data = map(data, ~pull(.,gene))) %>%
  .$data -> each_list
marker_krt7_within_all %>% distinct(cluster) %>% pull(cluster) -> use_name
names(each_list) <- use_name
each_list

data <- SetIdent(data, cells = use_id, value = "KRT7posi_Hep")
signature_tile(marker = "segal_list2")
id_ch("label4")


DoHeatmap(data, features = unlist(segal_list2), slot = "scale.data")



tile("marker_krt7_list")

tile(gene = use_gene[1:40])
data$label3 <- Idents(data)
data_sub <- sub_fil(data, str_detect(label3, "KRT7|HHyP"))
up(data_sub)
up()

marker_krt7_within_hhyp <- diff_test_vs(ident.1 = "KRT7posi_Hep", "HHyP")
marker_krt7_within_hhyp <- marker_krt7_within_hhyp %>% rownames_to_column(var = "gene") %>% mutate(cluster = if_else(avg_logFC>0, "KRT7posi_Hep", "HHyP"))
use_gene <- marker_krt7_within_hhyp %>% filter(avg_logFC>0) %>% head(50) %>% pull(gene)
use_gene <- marker_krt7_within_hhyp %>% filter(avg_logFC<0) %>% head(50) %>% pull(gene)
sav(marker_krt7_within_hhyp)
tile(use_gene)
marker_krt7_within_hhyp %>% arrange(-avg_logFC) %>% write_clip()


#KRT7+KRT19
df <- get_df(data_hep)
data_hep$KRT7_19 <- df %>% select(KRT7, KRT19) %>% transmute(a = (KRT7*KRT19)^1/2) %>% pull(a)
data_hep$KRT7_19
id_ch("KRT7_bin", data_hep)
id_ch("label2", data_hep)
up(data_hep)
vl(features = "KRT7_19", data_hep)
use_id <- pick_id(KRT7_19>0, object = data_hep)
data_hep <- SetIdent(data_hep, use_id, "krt7_19")
id_ch("label2", data_hep)
up(data_hep)

df_ <- data@assays$RNA@data[c("KRT7", "KRT19"),] %>% as.matrix() %>% t() %>% as.data.frame()
data$krt7_19 <- df_ %>% mutate(a = (KRT7*KRT19)^1/2) %>% pull(a)
id_ch("label2")
vl("krt7_19")
data<- SetIdent(object = data, cells = use_id, value = "krt7_19_")

tile("segal_list2")
up()
diff_test()
data$label3 <- Idents(data)
data_sub <- sub_fil(data, !str_detect(label3, "Hep"))
up(data_sub)
diff_test(data, logFC)
marker_ <- FindAllMarkers(data, logfc.threshold = log(2), min.pct = 0.25, only.pos = T)

#VIM
up(data_cho)
ump("VIM")
VlnPlot(data_cho, "VIM")
df_cho <- get_df(data_cho)
res_cor_vim <- do_cor(expr_df = df_cho, group_label = "vimentin", gene = "VIM", method = "pearson")
sav(res_cor_vim)
res_cor_vim %>% head(50) %>%
  ggplot(aes(fct_reorder(gene,cor), cor, fill = c("#D96464CF"))) + geom_bar(stat= "identity") + coord_flip() + guides(fill = F)
df_cho %>% ggplot(aes(VIM, TAGLN2)) + geom_point()
df_cho %>% ggplot(aes(VIM, MDK)) + geom_point()
res_cor_vim$gene %>% head(100) %>% write_clip()

use_id <- up(data_hep) %>% Seurat::CellSelector()
data_hep <- SetIdent(data_hep, use_id, "far")
add_meta_binval("KRT7", data_hep)
use_id_krt7 <- data_hep@meta.data %>% filter(KRT7_bin == "strong") %>% pull(id)
data_hep <- SetIdent(data_hep, use_id_krt7, "krt7")
up(data_hep)
id_ch("label2", data_hep)

marker_krt7_vs_far <- diff_test_vs("krt7", "far",object = data_hep)
marker_krt7_vs_far <- marker_krt7_vs_far %>% rownames_to_column(var = "gene") %>% mutate(cluster = if_else(avg_logFC>0, "krt7", "far"))
marker_krt7_vs_far %>% arrange(cluster, p_val_adj)
package2::tile("segal_list2", data_hep)


#zonation
use_gene <- df_hep[unlist(zone_list)] %>% apply(MARGIN = 2, mean) %>% sort(decreasing = T) %>% head(60) %>% names()
zone_list_fil <- zone_list %>% map(~.[.%in% use_gene])
save_list(zone_list_fil)

package2::tile("zone_list_fil", data_hep, order = F)
data_hep$label3 <- Idents(data_hep)
data_hep$label4 <- data_hep@meta.data %>%
  mutate(label4 = if_else(disease == "CH", paste0(label3,"_CH"), as.character(label3))) %>%
  pull(label4)

id_ch("label4", data_hep)
up(data_hep)
use_level <- data_hep$label4 %>% unique() %>% sort()
data_hep$label4 <- data_hep$label4 %>% fct_relevel(use_level)



up(data_hep)
signature_tile("zone_list_fil", object = data_hep)


#hormon

package2::tile("hormon_list")

#GS
add_sig_val_(data_hep, marker_list = "zone_list3", overwrite = T)
data_hep$zone_ratio3 <- data_hep@meta.data %>% mutate(zone_ratio = zone1_1 -zone3_1) %>% pull(zone_ratio)

VlnPlot(data_hep, c("zone1_1", "zone3_1", "zone_ratio3"))
VlnPlot(data_hep, c("zone1_1", "zone3_1", "zone_ratio3"), group.by = "disease")
res_cor_gs <- do_cor(expr_df = df_hep, group_label = "gs", gene = "GLUL", method = "pearson")
res_cor_gs %>% head(50) %>%
  ggplot(aes(fct_reorder(gene,cor), cor, fill = c("#D96464CF"))) + geom_bar(stat= "identity") + coord_flip() + guides(fill = F)
res_cor_gs %>% filter(cor>0.1) %>% pull(gene)-> gs_gene
signature_tile("zone_list2", data_hep)
res_enrich_gs <- gs_gene %>% convert_gene(.) %>% .$ENTREZID %>% geneano_enricher(gene_entrez = .)
res_enrich_gs %>% barplot(showCategory = 20)
res_enrich_gs@result %>% write_clip()
use_df <- read_clip_tbl() %>% filter(mark ==1)
res_enrich_gs@result <- use_df
res_enrich_gs_fil <- res_enrich_gs
sav(res_enrich_gs_fil)
res_enrich_gs %>% barplot
res_cor_gs <- res_cor_gs %>% mutate(color = if_else(gene %in% zone_list$zone3_, "red", "black"))
res_cor_gs %>% head(40) %>% ggplot(data = .,aes(fct_reorder(gene,cor), cor, fill = c("#D48C8C"))) + geom_bar(stat = "identity", show.legend = F) + coord_flip()+
  theme(axis.text.y = element_text(colour = res_cor_gs$color))

#HAL

add_sig_val(object = data_hep, "zone_list_fil", overwrite = T)
ump(c("zone1_", "zone3_"), object = data_hep)
data_hep$zone_ratio <- data_hep[[]] %>% mutate(zone1_ = scale(zone1_), zone3_ = scale(zone3_),zone_ratio = zone1_ - zone3_) %>% pull(zone_ratio)

ump(c("zone1_","zone3_", "zone_ratio2"),data_hep)

res_cor_hal <- do_cor(expr_df = df_hep, group_label = "HAL", gene = "HAL", method = "pearson")
res_cor_hal %>% head(50) %>%
  ggplot(aes(fct_reorder(gene,cor), cor, fill = c("#D96464CF"))) + geom_bar(stat= "identity") + coord_flip() + guides(fill = F)
res_cor_hal %>% filter(cor >0.1) %>% .$gene -> hal_gene
res_cor_hal_ <- res_cor_hal %>% filter(cor >0.1) %>% .$gene %>% convert_gene()
res_enrich_hal <- res_cor_hal_$ENTREZID %>% geneano_enricher(gene_entrez = .)
res_enrich_hal@result %>% write_clip()
res_cor_hal %>% filter(cor>0.1) %>% pull(gene)-> hal_gene
sav(res_enrich_hal)
use_df<- read_clip_tbl() %>% filter(mark == 1)
res_enrich_hal@result <- use_d
res_enrich_hal@pvalueCutoff <- 0.2
res_enrich_hal %>% enrichplot:::barplot.enrichResult(height = .,font.size = 20)
res_cor_hal <- res_cor_hal %>% mutate(color = if_else(gene %in% zone_list$zone1_, "red", "black"))
res_cor_hal %>% head(40) %>% ggplot(data = .,aes(fct_reorder(gene,cor), cor, fill = c("#D48C8C"))) + geom_bar(stat = "identity", show.legend = F) + coord_flip()+
  theme(axis.text.y = element_text(colour = res_cor_gs$color))

#cancer study
data$label_cancer <- data[[]] %>% mutate(a = if_else(disease %in% c("HCC", "ICC"), "C", "N" )) %>% pull(a)
id_ch("label_cancer")
marker_c_vs_n <- diff_test(data)

setwd("~/single_cell/single_cell_project/data/GSE125449/GSE125449_set1/GSE125449_Set1_samples.txt")
ma_id <- read_delim("GSE125449_Set1_samples.txt",delim =  "\t")
data[[]] %>% filter(str_detect(batch, "ma_set1")) %>% select(id)
ma_id$`Cell Barcode` <- paste0(ma_id$`Cell Barcode`, "_6")
ma_id$Type %>% table

malignant_id1 <- ma_id %>% filter(Type == "Malignant cell") %>% pull(`Cell Barcode`)


setwd("~/single_cell/single_cell_project/data/GSE125449/GSE125449_set2/GSE125449_Set2_samples.txt")
ma_id2 <- read_delim("GSE125449_Set2_samples.txt",delim =  "\t")
data[[]] %>% filter(str_detect(batch, "ma_set2")) %>% select(id)
ma_id2$`Cell Barcode` <- paste0(ma_id2$`Cell Barcode`, "_7")
malignant_id2 <-ma_id2 %>% filter(Type == "Malignant cell") %>% pull(`Cell Barcode`)
malignant_id <- c(malignant_id1, malignant_id2)

data$label_cancer <- data[[]] %>% mutate(a = if_else(id %in% malignant_id, "Malignant", "non_malignant")) %>% pull(a) %>%
  fct_relevel(c("Malignant", "non_malignant"))

data$label_cancer2 <- data[[]] %>% mutate(a = paste0(label2, "_", label_cancer), b = fct_reorder(a, as.numeric(label_cancer)),
                                          c = fct_reorder(b, as.numeric(label2))) %>% pull(c)


id_ch("label_cancer2")
up()
signature_tile("segal_list")
signature_tile("Epithelial_list2")
tile_plot("segal_list")
marker_c_vs_n <- diff_test(data)
add_sig_val(marker_list = "Epithelial_list", overwrite = T)
up() + scale_color_manual(values = c("#40C75B", "#CC5462"))
up() + scale_color_manual(values = gg_color_hue(2))

vl("epithelial")
tile_plot("Epithelial_list2")

marker_c_vs_n_batch <- diff_test_batch("label2")
marker_c_vs_n_batch %>% write_clip()
sub_c <- sub_fil(data, label_cancer == "Malignant")
id_ch(use_meta = "label2", sub_c)
marker_c <- diff_test(sub_c)
id_ch("disease", object = sub_c)
up(sub_c)


bar_origin(sample, label2, object = sub_c)
bar_origin(disease, label2, object = sub_c)
sub_c$sample <- sub_c$batch %>% str_extract("(?<=_)\\w{1,6}$") %>%  fct_relevel(c(paste0("HCC_",c(3,4,5,6,8,9)), paste0("ICC_",c(1,2,3,5,6,7,8,9,10))))
data$sample <- data$batch %>% str_extract("(?<=_)\\w{5,6}$") %>%  fct_relevel(c(paste0("HCC_",c(3,4,5,6,8,9)), paste0("ICC_",c(1,2,3,5,6,7,8,9,10))))

id_ch("label_cancer")
vl(c("cycle", "epithelial", colnames(data[[]])[3]), pt.size = 0)
sub_c$label_cancer3 <- paste0(sub_c$label2, "_", sub_c$sample)
sub_c$label_cancer3 <- sub_c[[]] %>% mutate(a = paste0(label2, "_", sample), b = fct_reorder(a, as.numeric(sample)),
                                          c = fct_reorder(b, as.numeric(label2))) %>% pull(c)

id_ch("label_cancer3", sub_c)
signature_tile("segal_list", object = sub_c)
marker_c
sav(sub_c)



sub_c <- sub_fil(data, disease %in% c("HCC", 'ICC'))
up(sub_c)
id_ch("label_cancer", sub_c)
res <- diff_test_batch("label2", object = sub_c)
up(sub_c)

mono <- make_monocle3(sub_c)
mono <- do_monocle(mono)
mop(mono, color_cells_by = "label_cancer")


id_ch("label_cancer2", sub_c)
signature_tile("segal_list", object = sub_c)
