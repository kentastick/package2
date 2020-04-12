

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

up(data_hep)
id_ch("label2", object = data_hep)
av_df_hep <- AverageExpression(data_hep, assays = "RNA", features = data_hep@assays$RNA@var.features)
av_df_hep_batch <- batch_mat(av_df_hep[[1]], object = data_hep)
cor(av_df_hep_batch, )
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
tile("zone_list",object = data_hep, )
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

signature

 #cholangiocyte
marker_cho2 <- diff_test(data_cho)

saveRDS(object = data_cho, file = "data_cho.rds")
sav(marker_cho)

marker_cho2 <- marker_cho2 %>% diff_marker_convert()
marker_cho2 <- marker_cho2 %>% do_geneano()
marker_cho2 <- marker_cho2 %>% mutate(bar = map2(res_enricher, cluster, ~barplot(.x, title = .y)))
marker_cho2 %>% mutate(a = map2(res_enricher, cluster,~.x@result %>% mutate(cluster = .y, n  = row_number()) %>% filter(p.adjust<0.05))) %>% pull(a) %>% 
  purrr::reduce(rbind) %>% write_clip()


marker_cho$bar


#KRT7 analaysis
VlnPlot(data_hep, features = "KRT7")
df_hep <- get_df(data_hep)
df_hep
res_cor_krt7 <- do_cor(expr_df = df_hep, group_label = "krt", gene = "KRT7", method = "pearson")
res_cor_krt7 %>% head(30) %>%  ggplot(aes(fct_reorder(gene,cor), cor, fill = c("#2CD2DB"))) + geom_bar(stat = "identity") +
  coord_flip() + guides(fill = F)+ xlab("correlation")
use_gene <- res_cor_krt7 %>% head(100) %>% .$gene 
df_hep_sub <- df_hep[c("KRT7", use_gene)]

df_hep_sub %>% ggplot(aes(KRT7, FXYD2)) + geom_point()
df_hep_sub %>% ggplot(aes(KRT7, TACSTD2)) + geom_point()
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


res_cor_vim <- do_cor(expr_df = df_hep, group_label = "krt", gene = "KRT7", method = "pearson")

