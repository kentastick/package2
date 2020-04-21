library(package2)
library(Seurat)
library(tidyverse)

up()

id_ch("seurat_clusters")

bar_origin(disease, label3)
bar_origin(label2, disease)

id_ch("label2")

data$label2 <- data$seurat_clusters %>% fct_collapse(Naieve_like = c("0","6", "7", "11","17"), Dysfunctional = c("2", "9"),
                                                     Transitional = c("1","8", "14", "18"),
                                                     cNK = c("3", "13"), lrNK = c("4","5"), Cytotoxic = c("10"), Treg = c("12"))
data$label2 <- data$label2 %>% fct_relevel(c("Naieve_like", "Cytotoxic", "Transitional", "Dysfunctional", "Treg", "cNK","lrNK"))
data$seurat_clusters <- data$seurat_clusters %>% fct_reorder(as.numeric(data$label2)) %>% fct_drop()
data$label3 <- data[[]] %>% mutate(n = plyr::mapvalues(seurat_clusters, from = levels(data$seurat_clusters),
                                                       to = c("(1)", "(2)", "(3)", "(4)", "(5)","","(1)","(2)","(3)", "(4)","(1)","(2)","", "(1)", "(2)", "(1)","(2)")),
                                   name = paste0(label2,n)) %>% mutate(name = fct_reorder(name, as.numeric(label))) %>%
  pull(name)

id_ch("label3")
up()

data <- sub_fil(data, disease !="BL")

#validation by sample correlation matrix
id_ch("label2")
data$label4 <- paste0(data$label3, "_", data$disease) %>% fct_reorder(as.numeric(as.factor(data$label2)))

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

#reguratory Tcell
sub_reg <- sub_fil(data, label2 == "Treg")

sub_pbc <- sub_fil(data, str_detect(disease,"PBC"))
up(sub_pbc)
id_ch("label2", object = sub_pbc)
up(sub_pbc)
marker_pbc <- FindAllMarkers(sub_pbc, only.pos = T)

DefaultAssay(sub_pbc) <- "RNA"

marker_pbc_rna <- FindAllMarkers(sub_pbc, only.pos = T)

a <- data$label3 %>% unique() %>% map(~sub_fil(data, label3 == .))

names(a) <- data$label3 %>% unique()
enframe(a) -> a
a <- a %>% mutate(cytokine = map2(value, name,~{tile_plot(list(cytokine =marker_list$cytokine), object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_cytokine.jpg"), path = "plot")}))
a <- a %>% mutate(chemokine= map2(value, name,~{tile_plot(list(cytokine =marker_list$chemokine), object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_chemokine.jpg"), path = "plot")}))
a <- a %>% mutate(costimulatory= map2(value, name,~{tile_plot(list(cytokine =marker_list$co_stimulatory_tcell), object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_costimulatory.jpg"), path = "plot")}))
a <- a %>% mutate(inhibitory= map2(value, name,~{tile_plot(list(cytokine =marker_list$co_inhibitory_tcell), object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_inhibitory.jpg"), path = "plot")}))
a <- a %>% mutate(signature = map2(value, name,~{signature_tile("NKT_list", object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_signature.jpg"), path = "plot")}))

a$value <- a$value %>% map(~{a <- .;Idents(a) <- "label4";return(a)})
a$value %>% map()

unique(a$disease)

