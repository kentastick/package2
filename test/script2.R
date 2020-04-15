library(Seurat)
library(package2)
library(tidyverse)


#labeling
data$label <- data$seurat_clusters %>% fct_collapse(Hep = c("1","3","4","5","11","14", "16"), im_Hep = "15",
                                                    HHyP = c("0","2","7","12"), BEC = c("9"),  prolifer= "17" )
data$seurat_clusters <- data$seurat_clusters %>% fct_reorder(as.numeric(data$label))
data$label2 <- data[[]] %>% mutate(n = plyr::mapvalues(seurat_clusters, from = levels(data$seurat_clusters),
                                                       to = c("", "(1)", "(2)", "(3)", "(4)", "(1)", "(2)","(3)", "(4)", "(5)", "(6)", "(7)", "", "", "")),
                                   name = paste0(label,n)) %>% mutate(name = fct_reorder(name, as.numeric(label))) %>%
  pull(name)
sa_data(hepato_cholangio_combined_filtered)


#data filtering
data <- sub_fil(data, str_detect(label2, "Hep")&cholangio_gm==0|!str_detect(label2, "Hep")&cholangio_gm!=0)
data_hep <- sub_fil(data, str_detect(label2, "Hep"), cholangio_gm==0)
data_cho <- sub_fil(data, !str_detect(label2, "Hep"), seurat_clusters !=6, cholangio_gm!=0)

#Monocle
mono <- make_monocle3(data_cho)

mono <- do_monocle(mono)
mop(mono, color_cells_by = "disease")

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



up(data)
add_meta_binval(object = data_hep ,gene = "KRT7")
 #extract krt7 strong positive cells within hepatocyte
use_id <- data_hep[[]] %>% filter(KRT7_bin == "strong") %>% pull(id)
data_hep <- SetIdent(data_hep, use_id, value = "krt7Hep")
up(data_hep)
id_ch("label2")

data <- SetIdent(data, cells = use_id, value = "KRT7posi_Hep")
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
use_df <- read_clip_tbl() %>% filter(mark ==1)
res_enrich_gs@result <- use_df
res_enrich_gs %>% barplot
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

res_cor_hal %>% filter(cor>0.1) %>% pull(gene)-> hal_gene
sav(res_enrich_hal)
use_df<- read_clip_tbl() %>% filter(mark == 1)
res_enrich_hal@result <- use_df
res_enrich_hal %>% enrichplot:::barplot.enrichResult()
