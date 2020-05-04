library(package2)
library(tidyverse)
library(purrr)
library(Seurat)
library(clipr)
library(magrittr)
bar_origin(disease, label)
up()
sa_data(Mesenchyme_combined_filtered)

rna()
mesenchyme_marker <- diff_test(data)
marker_2_vs_10 <- diff_test_vs(2,10)
marker_2_vs_10 %>% arrange(-cluster)
mesenchyme_marker %>% filter(cluster == 10) %>% head(30) %>% pull(gene) %>% write_clip()
mesenchyme_marker %>% filter(cluster == 10) %>% head(10) %>% pull(gene) %>% vl()
mesenchyme_marker %>% filter(cluster == 2) %>% head(30) %>% pull(gene) %>% tile_plot()


#mesenchyme marker list

mesenchyme_list <- list(Myofibroblast = c("TIMP1", "COL1A1", "COL1A2", "COL3A1","ELN", "PDGFRA"),
                        Mesenchyme = c("ACTA2", "VIM", "CRYAB", "NGFR"))
save_list(mesenchyme_list)

tile_plot("mesenchyme_list")
signature_tile("mesenchyme_list")
bar_origin(seurat_clusters, disease)

data$disease <- fct_relevel(data$disease, c("PBC_1", "PBC_2", paste0("NL_", 1:4), "FL", "CH", "HCC", "ICC"))


data@meta.data %>% select(label, disease) %>% group_by(disease, label) %>% summarise(n = n()) %>%
  group_by(disease) %>% mutate(all = sum(n), ratio = round(n/all*100, digits = 1)) %>%
  pivot_wider(id_cols = disease,names_from = label, values_from = ratio) %>% write_clip()

data$disease %>% table ->a
b <- a %>% as.list %>% enframe %>% mutate(value = unlist(value), value = paste0(name, ": ", value)) %>% pull(value)
b %>% reduce(~paste(.x,.y, sep = ", ")) %>% write_clip()


#clustering
df <- sig_val(object = data, "MFB_list")
df2 <- sig_val2(df)
df2 %>% pivot_wider(id_cols = c(cluster), names_from = signature, values_from = score) -> df3
df3 %>% as.data.frame() %>% column_to_rownames(var = "cluster") %>%  dist() %>% hclust() %>% plot



#labeling
data$label <- data$seurat_clusters %>% fct_collapse(MFB = c("10","2"), Mesenchyme = c("0","1","3","4","5","8","9")) %>%
  fct_relevel(c("MFB", "Mesenchyme"))

data$seurat_clusters <- data$seurat_clusters %>% fct_reorder(as.numeric(data$label)) %>% fct_drop()

#data$label2 <- data$seurat_clusters %>% fct_collapse(MFB= c("2", "10"),Mesenchyme= c("0","1","3","4","5","8","9"))

data$label2 <- data[[]] %>% mutate(n = plyr::mapvalues(seurat_clusters, from = c("2","10","0","1","3","4","5","8","9"), to = c("(1)", "(2)", "(1)", "(2)","(3)", "(4)", "(5)", "(6)", "(7)")),
                                   name = paste0(label, n)) %>%
  mutate(a = fct_reorder(name, as.numeric(label))) %>%
  pull(a)


data$label3 <- paste0(data$label, "_", data$disease)
data$label3 <- data[[]] %>% mutate(label3 = fct_reorder(label3, as.numeric(disease))) %>% pull(label3)
data$label3 <- data[[]] %>% mutate(label3 = fct_reorder(label3, as.numeric(label))) %>% pull(label3)

data$label4 <- paste0(data$label, "_", data$disease)
data$label4 <- data[[]] %>% mutate(label4 = fct_reorder(label4, as.numeric(disease))) %>% pull(label4)
data$label4 <- data[[]] %>% mutate(label4 = fct_reorder(label4, as.numeric(label))) %>% pull(label4)

#tile plot
id_ch("label4")
signature_tile("mesenchyme_list")
tile_plot("mesenchyme_list")


sa_data(mesenchyme_combined_fil)

#ump
up()

upd("PBC")

#proportion
bar_origin(disease, label2)
bar_origin(label, disease)

#data
data$disease %>% table

signature_tile("mesenchyme_list")
tile_plot("MFB_list")

#
mesenchyme_marker %>% group_by(cluster) %>% top_n(10, wt = avg_logFC) %>% View
mesenchyme_marker_422 <- mesenchyme_marker
sav(mesenchyme_marker_422)



#comparison object within cell
a <- tibble(cell_type = unique(dlata$label))
a$data <- a$data %>% map(~{a <- .;a$disease <- a$disease %>% fct_relevel(sort(unique(a$disease)));Idents(a) <- "disease";return(a)})
a <- a %>% mutate(data = map(cell_type, ~sub_fil(data, label == .x)))
a <- a %>% mutate(signature = map2(data, cell_type,~{signature_tile("MFB_list", object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_signature.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(signature2 = map2(data, cell_type,~{tile_plot("MFB_list", object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_sig_tile.jpg"), path = "comparison_cell")}))

a <- a %>% mutate(signature2 = map2(data, cell_type,~{tile_plot(marker_list$cytokine, object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_cytokine.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(signature2 = map2(data, cell_type,~{tile_plot(marker_list$chemokine, object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_chemokine.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(signature2 = map2(data, cell_type,~{up(.x)->p;ggsave(plot = p, filename = paste0(.y, "up.jpg"), path = "comparison_cell")}))
a <- a %>% mutate(signture_cytokine = map2(data, cell_type, ~signature_tile("cytokine_list", .x, title = .y)))


#comparison object within condition
b <- tibble(condition = unique(data$disease))
b <- b %>% mutate(data = map(condition, ~sub_fil(data, disease == .x)))
b$data <- b$data %>% map(~{a <- .;Idents(a) <- "disease";return(a)})
b$data <- b$data %>% map(~{a <- .; a$label2 <-a$label2 %>% fct_relevel(sort(unique(a$label2)));return(a) })

b <- b %>% mutate(signature = map2(data, condition,~{signature_tile("MFB_list", object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_signature.jpg"), path = "comparison_condition")}))
b <- b %>% mutate(signature2 = map2(data, condition,~{tile_plot("MFB_list", object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_signature2.jpg"), path = "comparison_condition")}))
b <- b %>% mutate(signature2 = map2(data, condition,~{tile_plot("gene_list", object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_norm.jpg"), path = "comparison_condition")}))
b <- b %>% mutate(cytokine = map2(data, condition,~{tile_plot(marker_list$cytokine, object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_cytokine.jpg"), path = "comparison_condition")}))
b <- b %>% mutate(chemokine = map2(data, condition,~{tile_plot(marker_list$chemokine, object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_chemokine.jpg"), path = "comparison_condition")}))

b <- b %>% mutate(ump = map2(data, condition,~{up(.x)->p;ggsave(plot = p, filename = paste0(.y, "_up.jpg"), path = "comparison_condition")}))


#
data_sub <- sub_fil(data, disease != "FL", label2 %in% c("MFB(1)", "MFB(2)"))
mono <- make_monocle3(data_sub)
mono <- do_monocle(mono)
mop(mono, color_cells_by = "label2")

mono_marker <- do_diff_mono(mono, group_cells_by = "seurat_clusters")
mop(mono, color_cells_by = "seurat_clusters")
mop(mono, color_cells_by = "label2")




