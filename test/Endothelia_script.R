library(package2)
library(Seurat)
library(tidyverse)

up()
upd("PBC")
id_ch("label")
signature_tile("Endothelial_list2")
tile_plot("Endothelial_list2")

bar_origin(disease, label,object = me)
#endothelial list

Endothelial_list2 <- list(
     LSEC = c("CLEC4M", "CLEC4G", "STAB2", "LYVE1"),
     Lympha = c("PDPN", "CCL21"),
     non_LSEC = c("VWF", "PECAM1", "CD34","PVAP"),
     Central = c("RSPO3", "WNT2"),
     Aeterial = c("AIF1L", "PPA1"))
     #SAEC = c("VWA1", "ACKR1", "COL4A1", "PLVAP"))
save_list(Endothelial_list2)

data$disease <- fct_relevel(data$disease, c("PBC_1", "PBC_2", paste0("NL_", 1:4), "FL", "CH", "HCC", "ICC"))

data$label2 <- interaction(data$label, data$disease, sep = "_", lex.order = T)
data$label2 <- interaction(data$label3, data$disease, sep = "_", lex.order = T)
id_ch("label")

#cell number
data$disease %>% table ->a
b <- a %>% as.list %>% enframe %>% mutate(value = unlist(value), value = paste0(name, ": ", value)) %>% pull(value)
b %>% reduce(~paste(.x,.y, sep = ", ")) %>% write_clip()


#labelling
data$label <- data$seurat_clusters %>% fct_collapse(., non_LSECs = c("0","6","7", "5"),
                                                    LSEC = c("3", "10"),
                                                    Central = c("2"),
                                                    Arterial = c("4"),
                                                    Lympha = c("8")
) %>% fct_relevel("LSEC", "Central", "Arterial", "non_LSECs", "Lympha")
data$label3 <- data$label %>% fct_collapse(non_LSECs = c("non_LSECs", "Central", "Arterial"))
id_ch("label3")

data$seurat_clusters <- data[[]] %>% mutate(seurat_clusters = fct_reorder(seurat_clusters, as.numeric(label))) %>% pull(seurat_clusters) %>%
fct_drop()
data$seurat_clusters

data$label2 <- data[[]] %>% mutate(n = plyr::mapvalues(seurat_clusters,from = c(3, 10, 2,4,0,5,6,7,8),
                                                       to = c("(1)","(2)","","","(1)","(2)","(3)","(4)", "")),
                                   name = paste(label, n, sep = "")) %>%
  mutate(name = fct_reorder(name, as.numeric(label))) %>%
  pull(name)

data$label3 <- paste0(data$label, "_", data$disease)
data$label3 <- data[[]] %>% mutate(label3 = fct_reorder(label3, as.numeric(disease))) %>% pull(label3)
data$label3 <- data[[]] %>% mutate(label3 = fct_reorder(label3, as.numeric(label2))) %>% pull(label3)


data@meta.data %>% select(label2, disease) %>% group_by(disease, label2) %>% summarise(n = n()) %>%
  group_by(disease) %>% mutate(all = sum(n), ratio = round(n/all*100, digits = 1)) %>%
  pivot_wider(id_cols = disease,names_from = label2, values_from = ratio) %>% write_clip()



#comparison within condition
b <- tibble(condition = unique(data$disease))

b <- b %>% mutate(data = map(condition, ~sub_fil(data, disease == .x)))

b$data <- b$data %>% map(~{a <- .;Idents(a) <- "disease";return(a)})
b$data <- b$data %>% map(~{a <- .; a$label <-a$label %>% fct_relevel(sort(unique(a$label)));return(a) })

b %>% mutate(tmp = map2(data, condition, ~sig_val(.x) %>% sig_val2 %>%
                          group_by(cluster) %>%
                          mutate(rank = row_number(-mean)) %>%
                          arrange(cluster) %>% filter(rank == 1) %>%
                          filter(!signature %in% c("Endothelia")) %>% pull(cluster) %>% paste0(., "_",.y))) %>% .$tmp %>%
  unlist-> rm_list


b %>%  mutate(a = map2(data, condition, ~.x@meta.data %>% group_by(label) %>% summarise(n = n()) %>% filter(n<=3) %>% pull(label) %>%
                         paste0(.y, "_",. ))) %>% .$a %>% unlist -> rm_list2

b$data %>% map(~.@meta.data$label %>% table)
data <- sub_fil(data, !label3 %in% c(rm_list,rm_list2))

b <- b %>% mutate(umap = map(data, ~up(.)))
b %>% filter(condition == "NL_3") %>% .$umap

#comparison object within cell
a <- tibble(cell_type = unique(data$label2))
a$data <- a$data %>% map(~{a <- .;a$disease <- a$disease %>% fct_relevel(sort(unique(a$disease)));Idents(a) <- "disease";return(a)})
a <- a %>% mutate(data = map(cell_type, ~sub_fil(data, label2 == .x)))

a <- a %>% mutate(signature = map2(data, cell_type,~{signature_tile("Endothelial_list", object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_signature.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(signature2 = map2(data, cell_type,~{tile_plot("Endothelial_list", object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_sig_tile.jpg"), path = "comparison_cell")}))

a <- a %>% mutate(signature2 = map2(data, cell_type,~{tile_plot(marker_list$cytokine, object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_cytokine.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(signature2 = map2(data, cell_type,~{tile_plot(marker_list$chemokine, object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_chemokine.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(signature2 = map2(data, cell_type,~{tile_plot(marker_list$Costimulator_apc, object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_costimulatory.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(signature2 = map2(data, cell_type,~{tile_plot(marker_list$co_inhibitoryr_apc, object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_coinhibi.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(signature2 = map2(data, cell_type,~{up(.x)->p;ggsave(plot = p, filename = paste0(.y, "_cytokine.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(signture_cytokine = map2(data, cell_type, ~signature_tile("cytokine_list", .x, title = .y)))


#correlation
df <- get_df(data)
res_wnt2 <- do_cor(df,gene = "WNT2")
res_aif1l <- do_cor(df, gene = "AIF1L")
res_vwa <- do_cor(df, gene = "VWA1")
