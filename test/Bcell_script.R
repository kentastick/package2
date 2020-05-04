library(package2)
library(Seurat)
library(tidyverse)
library(gghighlight)

up()
upd("PBC")

ump("MS4A1")

df <- get_df(data)
res_cor_cd9 <- do_cor(df, gene = "CD9")
res_cor_cd9 %>% head(200) %>% View


data$disease %>% table

data$label <- data$seurat_clusters %>% fct_collapse(Plasma_like = c("0", "8"),Bcell = c("1", "2", "3", "6","7"), MZBcell = c("5") )


id_ch("label")
up(data)
marker_list$plasmacellsRA %>% head(30) %>% tile_plot
marker_list$plasmacellsRA %>% head(30) -> use_gene
pick_gene("^CD\\d{1,2}$") %>% tile_plot()

bar_origin(disease, label)



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

