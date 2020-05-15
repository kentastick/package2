library(package2)
library(Seurat)
library(tidyverse)
library(gghighlight)
library(clipr)

#marker
NKT_list2 <- list(Tcell = c("CD3D", "CD4", "CD8A", "CD8B", "TRAC", "TRBC2"),
                 NK = c("NCR1", "NCR2", "NCR3", "KLRF1", "NKG2D", "NLRC1", "NKG7", "GNLY", "FGFBP2"),
                 cNK = c("TBX21", "FCGR3A", "CX3CR1"),
                 lrNK = c("EOMES", "NCAM1", "CD69", "CXCR6", "CCR5"),
                 Naieve = c("IL7R", "CD4", "CCR7", "TCF7"),
                 CYtotoxic= c("GZMA", "GZMB", "GZMH", "GZMK", "PRF1"),
                 Dysfunction = c("PDCD1", "LAG3", "TIGIT", "CTLA4"),
                 Regurational = c("FOXP3", "IL2RA", "TNFRSF4"))
save_list(NKT_list2)
id_ch("label4")
signature_tile("NKT_list2")
bar_origin(disease, label2)

#data


#labeling
data$disease <- fct_relevel(data$disease, c("PBC_1", "PBC_2", paste0("NL_", 1:4), "FL", "CH", "HCC", "ICC"))

data$label2 <- data$seurat_clusters %>% fct_collapse(Naieve_like = c("0","17"), Dysfunctional = c("2", "9"),
                                                     Transitional = c("1","8", "14", "18"),
                                                     cNK = c("3", "13"), lrNK = c("4","5"), Cytotoxic = c("10"), Treg = c("12"))
data$label2 <- data$label2 %>% fct_relevel(c("Naieve_like", "Cytotoxic", "Transitional", "Dysfunctional", "Treg", "cNK","lrNK"))


data$seurat_clusters <- data$seurat_clusters %>% fct_reorder(as.numeric(data$label2)) %>% fct_drop()

data$label3 <- data[[]] %>% mutate(n = plyr::mapvalues(seurat_clusters, from = levels(data$seurat_clusters),
                                                       to = c("(1)", "(2)", "(3)", "(4)", "(5)","","(1)","(2)","(3)", "(4)","(1)","(2)","", "(1)", "(2)", "(1)","(2)")),
                                   name = paste0(label2,n)) %>% mutate(name = fct_reorder(name, as.numeric(label))) %>%
  pull(name)

data$label4 <- paste0(data$label2, "_", data$disease)

data$label4 <- data@meta.data %>% mutate(label4 = fct_reorder(label4, as.numeric(disease))) %>% pull(label4)
data$label4 <- data@meta.data %>% mutate(label4 = fct_reorder(label4, as.numeric(label2))) %>% pull(label4)


c("TNFSF15", "POU2AF1",
  "CD80", "IKZF3-ZPBP2-", "GSDMB-ORMDL3",
  "IL7R", "NFKB1",
  "STAT4", "TNFAIP2",
  "CXCR5", "MAP3K7IP1",
  "DENND1B")
data$label5 <- interaction(data$label2, data$disease, sep = "_", lex.order = T)
id_ch("label5")

id_ch("label2")
signature_tile("NKT_list2")
tile_plot("NKT_list2")
bar_origin(disease, label2,object = data_)

#removing peripheral cells
data <- sub_fil(data, disease != "BL", !seurat_clusters %in% c("6", "7", "11"))

sa_data(NKT_combined_filtered)



#validation by sample correlation matrix
id_ch("label2")
data$label4 <- paste0(data$label3, "_", data$disease) %
>% fct_reorder(as.numeric(as.factor(data$label2)))

signature_tile("NKT_list")
tile_plot("NKT_list")
id_ch("label4")

up()
av_df <- AverageExpression(data, assays = "RNA", features = data@assays$RNA@var.features)
av_df_batch <- batch_mat(av_df[[1]], object = data)
batch_cor_heatmap(av_df_batch, method = "spearman")


VlnPlot(data_cho, features = "cholangio_gm", group.by = "disease")
data_cho <- sub_fil(data_cho, cholangio_gm !=0)
up(data_cho)
data_cho <- sub_fil(data, !str_detect(label2, "Hep"))
id_ch("label2", object = data_cho)
av_df_cho <- AverageExpression(data_cho, assays = "RNA", features = data_cho@assays$integrated@var.features)
av_df_cho_batch <- batch_mat(av_df_cho[[1]], object = data_cho)
batch_cor_heatmap(av_df_cho_batch, method = "spearman")


#comparison by cells

a <- tibble(cell = unique(data$label2))
a <- a %>% mutate(data = map(cell, ~sub_fil(data, label2 == .x)))
a$data <- a$data %>% map(~{a <- .;Idents(a) <- "disease";return(a)})

a <- a %>% mutate(signature = map2(data, cell,~{up(object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_up.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(signature = map2(data, cell,~{signature_tile("NKT_list2", object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_signature.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(signature2 = map2(data, cell,~{signature_tile("NKT_list2", object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_signature2.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(cytokine = map2(data, cell,~{tile_plot(marker_list$cytokine, object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_cytokine.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(cytokine = map2(data, cell,~{tile_plot(marker_list$cytokine, object = .x, title = .y)->p}))
a <- a %>% mutate(chemokine = map2(data, cell,~{tile_plot(marker_list$chemokine, object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_chemokine.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(chemokine = map2(data, cell,~{tile_plot(marker_list$chemokine, object = .x, title = .y) ->p}))
a <- a %>% mutate(co_stimulatory = map2(data, cell,~{tile_plot(marker_list$co_stimulatory_tcell, object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_costimu.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(co_stimulatory = map2(data, cell,~{tile_plot(marker_list$co_stimulatory_tcell, object = .x, title = .y)->p}))
a <- a %>% mutate(co_inihibi = map2(data, cell,~{tile_plot(marker_list$co_inhibitory_tcell, object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_coinhibi.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(co_inihibi = map2(data, cell,~{tile_plot(marker_list$co_inhibitory_tcell, object = .x, title = .y)->p}))
a <- a %>% mutate(cytokine_sig = map2(data, cell,~{signature_tile("cytokine_list", object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_tile.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(cytokine_sig = map2(data, cell,~{signature_tile("cytokine_list", object = .x, title = .y)->p}))
a <- a %>% mutate(imche_sig = map2(data, cell,~{signature_tile("immunocheck_list", object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_immuncheck.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(imche_sig = map2(data, cell,~{signature_tile("immunocheck_list", object = .x,title = .y)->p}))
a <- a %>% mutate(pro_vs_anti_sig = map2(data, cell,~{signature_tile("inflamatory_list", object = .x,title = .y)->p}))
a <- a %>% mutate(pro_vs_anti = map2(data, cell,~{tile_plot("inflamatory_list", object = .x,title = .y)->p;ggsave(plot = p, filename = paste0(.y, "_anti_pro.jpg"), path = "comparison_cell", width = 3.18, height = 3.09 )}))
a <- a %>% mutate(pro_vs_anti = map2(data, cell,~{tile_plot("inflamatory_list", object = .x,title = .y)->p;ggsave(plot = p, filename = paste0(.y, "_anti_pro.jpg"), path = "comparison_cell", width = 3.18, height = 3.09 )}))
a <- a %>% mutate(stim_vs_inhi = map2(data, cell,~{tile_plot("immunocheck_list", object = .x,title = .y)->p;ggsave(plot = p, filename = paste0(.y, "_stim_inhibi.jpg"), path = "comparison_cell", width = 4.23, height = 3.98 )}))
a <- a %>% mutate(graph = map2(data,cell,~vl(c("IFNG","IL18","IL12B", "IL1B"), ncol =4, object = .x)))
a <- a %>% mutate(data = map(data,~{a <-.x; DefaultAssay(a) <-"RNA";return(a)}))


id_ch("label4")
signature_tile("inflamatory_list")
tile_plot(gene = marker_list$cytokine)


immunocheck_list <- list(co_stimulatory = marker_list$co_stimulatory_tcell,co_inhibitory = marker_list$co_inhibitory_tcell)
sav(a)


#Treg proportion
data@meta.data %>% select(label2, disease) %>% group_by(disease, label2) %>% summarise(n = n()) %>%
  group_by(disease) %>% mutate(all = sum(n), ratio = round(n/all*100, digits = 1)) %>%
  pivot_wider(id_cols = disease,names_from = label2, values_from = ratio) %>% write_clip()


id_ch("label4")
inflamatory_list <- list(pro = marker_list$pro_inflammatory, anti = marker_list$anti_inflammatory)

tile_plot(list(pro = marker_list$pro_inflammatory, anti = marker_list$anti_inflammatory),object = sub_fil(data, label2 == "Treg"))


signature_tile("inflamatory_list",sub_fil(data, label2 == "Treg"))

data[[]] %>% dplyr::select(disease, label2) %>% filter(!label2 %in% c("cNK", "lrNK")) %>%
  group_by(disease, label2) %>% summarise(n = n()) %>% group_by(disease) %>%
  mutate(sum = sum(n), prop = n/sum*100) %>%  filter(label2 %in%  c("Cytotoxic", "Treg")) %>% ggplot(aes(disease, prop, fill = label2)) +
  geom_bar(position = "dodge",stat = "identity") + facet_wrap(~label2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA)) +
  scale_fill_manual(values = gg_color_hue(8)[c(2,5)]) +
  labs(y = "Percent (%)")
  pivot_wider(id_cols = c(disease), names_from = label2, values_from = prop) %>% write_clip()
  mutate(Treg_ratio = Treg/sum*100)


res_enrich <- res_all %>% filter(cor>0.25) %>% pull(gene) %>%  convert_gene() %>% .$ENTREZID %>% geneano_enricher()
