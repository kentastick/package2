library(package2)
library(Seurat)
library(tidyverse)
library(clipr)
library(gghighlight)

up()

id_ch("label4")
signature_tile("macrophage_list")

id_ch("label5")
signature_tile("macrophage_list")

bar_origin(disease, label4)

tile_plot("macrophage_list")
save_list(macrophage_list)

#marker list
macrophage_list <- list(Mono = c("FCN1", "MNDA", "LYZ", "VCAM", "S100A8", "S100A12"),
               KC =c("CD14", "CD163", "MARCO", "VCAM1", "TIMD4", "CD5L", "C1QA", "FCGR3A", "APOE"),
               DC = c("CD1C", "CD1A", "CLEC9A", "CLEC10A", "CD274", "FCER1A", "HLA-DQA1", "HLA-DQB1"),
               pDC = c("LILRB4", "IRF7", "IRF8", "TCF4", "GZMB", "TCL1A"),
               SAM = c("TREM2", "CD9"),
               Cycling = c("MKI67", "CCNA2", "CCNB2", "STMN1"))

immunocheck_list2 <- list(co_stimulatory = marker_list$co_stimulatory_apc, co_inhibitory = marker_list$co_inhibitoryr_apc)
save_list(immunocheck_list2)

#rename
id_ch("label5")

data$disease <- fct_relevel(data$disease, c("PBC_1", "PBC_2", paste0("NL_", 1:4), "FL", "CH", "HCC", "ICC"))

data$label <- data$seurat_clusters %>% fct_collapse(cMono = c("0","2"),
                                                    tMono =c("3", "4", "10", "11"),
                                                    DC = c("5"),
                                                    KC = c("1", "7"),
                                                    SAM = c("9"),
                                                    cycling ="13"
) %>% fct_relevel(c("cMono","tMono", "DC", "KC", "SAM", "cycling"))

data$seurat_clusters <- data$seurat_clusters %>% fct_reorder(as.numeric(data$label)) %>% fct_drop()

data$label2 <- data[[]] %>% mutate(n = plyr::mapvalues(seurat_clusters, from = levels(data$seurat_clusters),
                                                       to = c("(1)", "(2)", "(1)", "(2)", "(3)", "(4)", "","(1)", "(2)","","")), name = paste0(label, n)) %>%
  pull(name)




data$label4 <- data@meta.data %>% mutate(a = if_else(label2 %in% c("tMono(2)", "tMono(3)","tMono(4)"), "DC_like_Mono", as.character(label))) %>% pull(a)
data$label4 <- data$label4 %>% fct_recode(SAM = "SM")
data$label4 <- data$label4 %>% fct_relevel(c("tMono", "DC_like_Mono", "DC", "KC","SAM", "cycling"))

data$label3 <- paste0(data$label2, "_", data$disease) %>% fct_reorder(as.numeric(as.factor(data$label)))

data$label5 <- paste0(data$label4, "_", data$disease)
data$label5 <- data@meta.data %>% mutate(label5 = fct_reorder(label5, as.numeric(disease))) %>% pull(label5)
data$label5 <- data@meta.data %>% mutate(label5 = fct_reorder(label5, as.numeric(label4))) %>% pull(label5)
id_ch("label5")

up()

sa_data(MP_combined_last)


data <- sub_fil(data, disease != "BL", !label2 %in% c("cMono(1)", "cMono(2)"))



bar_origin(disease, label4)
bar_origin(label4, disease)

id_ch("label5")
signature_tile(marker = "macrophage_list")
tile_plot("macrophage_list")

#validation by sample correlation matrix
id_ch("label2")



signature_tile("NKT_list")
tile_plot("macrophage_list")
id_ch("label4")

av_df <- AverageExpression(data, assays = "RNA", features = data@assays$RNA@var.features)
av_df_batch <- batch_mat(av_df[[1]], object = data)
batch_cor_heatmap(av_df_batch)

#comparison object within cell
a <- tibble(cell_type = unique(data$label4))
a <- a %>% mutate(data = map(cell_type, ~sub_fil(data, label4 == .x)))
a$data <- a$data %>% map(~{a <- .;Idents(a) <- "disease";return(a)})

a <- a %>% mutate(signature = map2(data, cell_type,~{signature_tile("macrophage_list", object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_signature.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(signature = map2(data, cell_type,~{signature_tile("macrophage_list", object = .x, title = .y)->p}))

a <- a %>% mutate(cytokine = map2(data, cell_type,~{tile_plot(marker_list$cytokine, object = .x, title = .y)->p}))
a <- a %>% mutate(chemokine = map2(data, cell_type,~{tile_plot(marker_list$chemokine, object = .x, title = .y)->p}))
a %>% mutate(chemokine = map2(data, cell_type,~{tile_plot(marker_list$chemokine, object = .x, title = .y)->p; ggsave(plot = p, filename = paste0(.y, "_chemokine.jpg"), path = "comparison_cell")}))

a %>% mutate(TLR = map2(data, cell_type,~{tile_plot(marker_list$TLR, object = .x, title = .y)->p;
  ggsave(plot = p, filename = paste0(.y, "_TLR.jpg"), path = "comparison_cell", width = 3.32, height = 3.15)}))

a <- a %>% mutate(TLR = map2(data, cell_type,~{tile_plot(marker_list$TLR, object = .x, title = .y)->p}))
a <- a %>% mutate(co_stimu = map2(data, cell_type,~{tile_plot(marker_list$co_stimulatory_apc, object = .x, title = .y)->p}))
a <- a %>% mutate(co_inhibi = map2(data, cell_type,~{tile_plot(marker_list$co_inhibitoryr_apc, object = .x, title = .y)->p}))
a <- a %>% mutate(stimu_ibhibi = map2(data, cell_type,~{tile_plot("immunocheck_list2", object = .x, title = .y)->p}))
a <- a %>% mutate(stimu_ibhibi_sig = map2(data, cell_type,~{signature_tile("immunocheck_list2", object = .x, title = .y)->p}))




a <- a %>% mutate(signature2 = map2(data, cell_type,~{tile_plot(marker_list$cytokine, object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_cytokine.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(signature2 = map2(data, cell_type,~{tile_plot(marker_list$chemokine, object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_chemokine.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(signature2 = map2(data, cell_type,~{tile_plot(marker_list$Costimulator_apc, object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_costimulatory.jpg"), path = "comparison_cell")}))

a <- a %>% mutate(signature2 = map2(data, cell_type,~{tile_plot(marker_list$co_inhibitoryr_apc, object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_coinhibi.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(signature2 = map2(data, cell_type,~{up(.x)->p;ggsave(plot = p, filename = paste0(.y, "_up.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(vl = map2(data, cell_type,~{.x->a;DefaultAssay(a)<-"RNA"; vl(features = use_gene, object = a)->p}))

#comparison object within condition
b <- tibble(condition = unique(data$disease))
b <- b %>% mutate(data = map(condition, ~sub_fil(data, disease == .x)))
b$data <- b$data %>% map(~{a <- .;Idents(a) <- "disease";return(a)})
b$data <- b$data %>% map(~{a <- .; a$label2 <-a$label2 %>% fct_relevel(sort(unique(a$label2)));return(a) })

b <- b %>% mutate(signature = map2(data, condition,~{signature_tile("macrophage_list", object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_signature.jpg"), path = "comparison_fil")}))
b <- b %>% mutate(signature2 = map2(data, condition,~{tile_plot("macrophage_list", object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_signature2.jpg"), path = "comparison_fil")}))
b <- b %>% mutate(signature2 = map2(data, condition,~{tile_plot("gene_list", object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_norm.jpg"), path = "comparison_fil")}))
b <- b %>% mutate(cytokine = map2(data, condition,~{tile_plot(marker_list$cytokine, object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_cytokine.jpg"), path = "comparison_fil")}))
b <- b %>% mutate(chemokine = map2(data, condition,~{tile_plot(marker_list$chemokine, object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_chemokine.jpg"), path = "comparison_fil")}))

b <- b %>% mutate(ump = map2(data, condition,~{up(.x)->p;ggsave(plot = p, filename = paste0(.y, "_up.jpg"), path = "comparison_condition")}))


#removing low quolity and number cells

b %>% mutate(tmp = map2(data, condition, ~sig_val(.x) %>% sig_val2 %>%
                 group_by(cluster) %>%
                 mutate(rank = row_number(-mean)) %>%
                 arrange(cluster) %>% filter(rank == 1) %>%
                 filter(!signature %in% c("MP", "Cycling")) %>% pull(cluster) %>% paste0(., "_",.y))) %>% .$tmp %>%
  unlist-> rm_list


b %>%  mutate(a = map2(data, condition, ~.x@meta.data %>% group_by(label2) %>% summarise(n = n()) %>% filter(n<=3) %>% pull(label2) %>%
                  paste0(.y, "_",. ))) %>% .$a %>% unlist -> rm_list2

data <- sub_fil(data, !label3 %in% c(rm_list,rm_list2))


b$data %>% map(~.@meta.data$disease %>% table)


use_gene <- c("CCL5", "CCL8","CXCL10", "CXCL14" )



#monocle

mono <- make_monocle3(data)
mono <- do_monocle(mono)
mop(mono, color_cells_by = "disease")
mono <- monocle3::order_cells(mono)
mono_sub <- monocle3::choose_cells(mono)
monocle3::plot_cells(mono,
                     color_cells_by = "pseudotime",
                     label_cell_groups=FALSE,
                     label_leaves=FALSE,
                     label_branch_points=FALSE,
                     graph_label_size=1.5)

