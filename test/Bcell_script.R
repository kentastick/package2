library(package2)
library(Seurat)
library(tidyverse)
library(gghighlight)

up()

tile_plot("gene_list", Bcell)
tile_plot("Bcell_list2", Bcell)
tile_plot(marker_list$CXCR, Bcell)
save_list(Bcell_list)
Bcell_list <- edit(Bcell_list)
Bcell <- add_info(Bcell)

Bcell$disease <- fct_relevel(Bcell$disease, c("PBC_1", "PBC_2", paste0("NL_", 1:4), "FL", "CH", "HCC", "ICC"))
res_cor_mzb1 <- do_cor(expr_df = df, gene = "MZB1")
Bcell_list2 <- list(naieve = c("TCL1A", "IGHD", "CXCR4", "FCER2", "CD22", "HVCN1"),
                    Plasmablast = c("IGJ", "JCHAIN", "PPIB", "HSP90B1", "MZB1", "SRGN","SEC11C", "PDIA4", "TXNDC5", "XBP1", "FAM46C"),
                    Memory = c("IGHG1"))

data$disease %>% table
Bcell <- sub_fil(Bcell, !seurat_clusters %in% c("4","8", "9"))

Bcell$label <- Bcell$seurat_clusters %>% fct_collapse(Plasma_like = c("0", "10"),
                                                    Naieve = c("1"),
                                                    Naieve_lke = "7",
                                                    Memory = c("2", "3"),
                                                    Memory_like = "6",
                                                    Transitional = c("5") )

id_ch("label3", Bcell)
up(Bcell)
Bcell$label2 <- interaction(Bcell$label, Bcell$disease, sep = "_", lex.order = T)
Bcell$label3 <- interaction(Bcell$seurat_clusters, Bcell$disease, sep = "_", lex.order = T)
tile_plot(gene = "Bcell_list", object = Bcell)

data$label2 <- paste0(data$label, "_", data$disease)
data$label2 <- data[[]] %>% mutate(label2 = fct_reorder(label2, as.numeric(disease))) %>% pull(label2)
data$label2 <- data[[]] %>% mutate(label2 = fct_reorder(label2, as.numeric(label))) %>% pull(label2)

data$label3 <- paste0(data$seurat_clusters, "_", data$disease)
data$label3 <- data[[]] %>% mutate(label3 = fct_reorder(label3, as.numeric(disease))) %>% pull(label3)
data$label3 <- data[[]] %>% mutate(label3 = fct_reorder(label3, as.numeric(label))) %>% pull(label3)

id_ch("label2")

up(data)
marker_list$plasmacellsRA %>% head(30) %>% tile_plot
marker_list$plasmacellsRA %>% head(30) -> use_gene
pick_gene("^CD\\d{1,2}$") %>% tile_plot()

bar_origin(disease,label, object = Bcell)




#comparison object within condition
b <- tibble(condition = unique(data$disease))
b <- b %>% mutate(data = map(condition, ~sub_fil(data, disease == .x)))
b$data <- b$data %>% map(~{a <- .;Idents(a) <- "label";return(a)})
b$data <- b$data %>% map(~{a <- .; a$label2 <-a$label2 %>% fct_relevel(sort(unique(a$label2)));return(a) })

b <- b %>% mutate(signature = map2(data, condition,~{signature_tile("macrophage_list", object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_signature.jpg"), path = "comparison_condition")}))
b <- b %>% mutate(signature2 = map2(data, condition,~{tile_plot("macrophage_list", object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_signature2.jpg"), path = "comparison_condition")}))
b <- b %>% mutate(signature2 = map2(data, condition,~{tile_plot("gene_list", object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_norm.jpg"), path = "comparison_condition")}))
b <- b %>% mutate(cytokine = map2(data, condition,~{tile_plot(marker_list$cytokine, object = .x)->p;ggsave(plot = p, filename = paste0(.y, "cytokine.jpg"), path = "comparison_condition")}))
b <- b %>% mutate(chemokine = map2(data, condition,~{tile_plot(marker_list$chemokine, object = .x)->p;ggsave(plot = p, filename = paste0(.y, "_chemokine.jpg"), path = "comparison_condition")}))
b <- b %>% mutate(ump = map2(data, condition,~{up(.x)->p;ggsave(plot = p, filename = paste0(.y, "_up.jpg"), path = "comparison_condition")}))

#monocle analysis

mono <- make_monocle3(Bcell)
mono <- do_monocle(mono)
mop(mono, color_cells_by = "seurat_clusters")
