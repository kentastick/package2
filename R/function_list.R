# data import function ----------------------------------------------------

impdata <- function(name) {
  rds_list <- list.files(path = "C:/Users/ken/Documents/single_cell/single_cell_analysis/data", pattern = "rds", full.names = T, recursive = T)
  importdata <- rds_list %>% str_subset(name)
  if(length(importdata) != 1){
    print(importdata)
    stop("more than two duplicated pattern exist")
  }
  return(readRDS(importdata))
}


list_data <- function() {
  list.files(path = "C:/Users/ken/Documents/single_cell/single_cell_analysis/data", recursive = T, pattern = ".rds")
}

# do_seurat ---------------------------------------------------------------


do_seurat <- function(data) {
  object <- CreateSeuratObject(counts = data, project = "object", min.cells = 3, min.features = 200)

  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave("plot1.jpg", device = "jpeg")

  plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  CombinePlots(plots = list(plot1, plot2))
  ggsave("plot2.jpg", device = "jpeg")

  #object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)

  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(object), 10)

  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(object)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  CombinePlots(plots = list(plot1, plot2))
  ggplot2::ggsave("plot3.jpg")

  all.genes <- rownames(object)

  #normalization data mean =1 variaty =0
  object <- ScaleData(object, features = all.genes)
  object <- RunPCA(object, features = VariableFeatures(object = object))

  print(object[["pca"]], dims = 1:5, nfeatures = 5)

  VizDimLoadings(object, dims = 1:2, reduction = "pca")
  ggsave("plot4.jpg")

  DimPlot(object, reduction = "pca")
  ggsave("plot5.jpg")

  DimHeatmap(object, dims = 15:30, cells = 500, balanced = TRUE)
  ggsave("plot6.jpg")


  #jackstraw analysis
  object <- JackStraw(object, num.replicate = 100)
  object <- ScoreJackStraw(object, dims = 1:20)

  JackStrawPlot(object, dims = 1:20)
  ggsave("plot7.jpg")



  object <- FindNeighbors(object, dims = 1:27)
  object <- FindClusters(object, resolution = 0.8)
  object <- RunUMAP(object, dims = 1:30)
  object <- RunTSNE(object, dims = 1:30)

  DimPlot(object, reduction = "umap")
  ggsave("plot8.jpg")

  DimPlot(object, reduction = "tsne")
  ggsave("plot9.jpg")

  data <- object
  rm(object)
  tmap(c('ALB', 'KRT7', 'CD68', 'CD3D', 'CD79A', 'IGJ','CD34', 'GNLY', 'FGFBP2', 'CD14', "MZB1"))
  ggsave("plot10.jpg")


  saveRDS(data, "normal.rds")
  return(object)
}



# wrapper function --------------------------------------------------------

#create present time as name
make_time <- function() {
  Sys.time() %>% str_remove_all('[:punct:]|\\s')
}


#plot
up <- function(..., label= TRUE) {
  DimPlot(object = data, label = label,...)
}

ups <- function(..., label = TRUE,type = 'png') {
  DimPlot(object = data, label = label,...) %>% plot_save(type = type)
}

ts <- function(object = data, ...) {
  DimPlot(object = object, reduction = 'tsne', label = TRUE, ...)
}

tmap <- function(features, object = data, ...) {
  FeaturePlot(features = features,object = object, reduction = 'tsne',  cols = c('lightgray', 'red'),...)
}

tmps <- function(features, object = data, type = 'png') {
  tmap(features, object = object) %>% plot_save( name = features, type = type)
}

ump <- function(features, object = data,...) {
  FeaturePlot(features = features,object = object, cols = c('lightgray', 'red'), min.cutoff = 0,...)
}

umps <- function(features, type = 'png') {
  ump(features) %>% plot_save(type = type)
}

vl <- function(...) {
  VlnPlot(object = data,...)
}
vls <- function(..., type = 'png') {
  VlnPlot(object = data,...) %>% plot_save(type = type)
}




#file save

sav <- function(object, name, path) {
  saveRDS(object = object, file = paste(path, name,".rds"))
}

savc <- function(...) {
  sav(path = "data/combined/",...)
}



# scraping function -------------------------------------------------------

#page <- rvest::read_html("https://www.genecards.org/Search/Keyword?queryString=Chromogranin")
#page %>% rvest::html_nodes("a")



# ggsave wrapper ----------------------------------------------------------


plot_save <- function(object, name = make_time(), type = 'tiff') {
  ggsave(plot =  object, filename = paste0('data/combined/plot/', name,'.',type), device = type)
}



# get gene data -----------------------------------------------------------

get_gene <- function(object = data) {
  object@assays$RNA@data %>% rownames
}



# pick up cells by gene value ---------------------------------------------

pick_cell <- function(gene, value) {
  a <-com@assays$RNA@data[gene,] > value
  cell_label <- a[a] %>% names()
  return(cell_label)
}


# make ident class --------------------------------------------------------


make_class <- function(gene, value, object = data) {
  id <- pick_cell(gene, value)
  object[[]] %>% rownames_to_column(var = "ID") %>%
    mutate(new_class = if_else(id %in% krt7posi, "positive", "negative")) %>% pull(new_class)
}


# search gene in Seurat object --------------------------------------------


search_gene <- function() {
  data@assays$RNA@data@Dimnames[[1]] %>% View
}


# output gene vector searched from seurat object  -------------------------

pick_gene <- function(pattern) {
  data@assays$RNA@data@Dimnames[[1]] %>% str_subset(., pattern =pattern)
}



# return absolute path ----------------------------------------------------

abpath <- function(path = clipr::read_clip()) {
  normalizePath(path) %>% stringr::str_replace_all("\\\\", "/") %>% clipr::write_clip()
}


# signature value calcuration ---------------------------------------------

#calculate geometric_mean of each cells
sig_val <- function(gene_list= gene_list, object = data) {
  mt <- object@meta.data
  count_mt <- data@assays$RNA@data
  gene_name <- rownames(count_mt)
  gene_list <- purrr::map(gene_list, ~.[. %in% gene_name])
  for(i in seq_along(gene_list)){
    sub_mt <- count_mt[gene_list[[i]],]
    value <- apply(sub_mt,2, gm_mean1)
    mt[names(gene_list)[i]] <- value
  }
  mt <- mt[names(gene_list)]
  return(mt)
}




# geometric mean ----------------------------------------------------------

gm_mean1 = function(a){prod(a)^(1/length(a))}
gm_mean1(1:10)


gm_mean2 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


# signature plot ---------------------------------------------------------------------

sig_val2 <- function(score_mt, object = data, gene_list = gene_list, non_filter = F) {
  if(!non_filter) {val_mean <- apply(score_mt, 2, mean) # signature mean
  for(i in seq_along(gene_list)){
    temp <- score_mt[[names(gene_list)[i]]]
    score_mt[[names(gene_list)[i]]] <- if_else(temp> val_mean[[i]], temp, 0)
  }
}
  score_mt$cluster <- object@meta.data %>% .$seurat_clusters


  score_mt %>% gather(-cluster, key = "signature", value = "score") %>%
    group_by(signature, cluster) %>%
    summarise(fraction_of_cells = sum(score>0)/n(), score = mean(score)) %>%
    group_by(signature) %>%
    mutate(max = max(score)) %>%
    mutate(score = score/max)
}





signature_plot <- function(score_mt, use.color = c("#0099FF", "#FAF5F5", "#E32020")) {

    score_mt %>% ggplot(aes(cluster, signature, colour =score, size = fraction_of_cells)) + geom_point() +
    scale_colour_gradientn(colours = c("red","yellow","white","lightblue","darkblue"),
                           values = c(1.0,0.7,0.6,0.4,0.3,0))
}

# data.frame to list ------------------------------------------------------

df_to_list <- function(df) {
  df_list <- vector("list", nrow(df))
  df <- df %>% as.matrix()
  for (i in 1:nrow(df)) {
    sub_row <- df[i,]
    names(sub_row) <- NULL
    df_list[[i]] <- sub_row
    names(df_list)[i] <- rownames(df)[i]
  }
  df_list <- map(df_list, ~.[!is.na(.)])
  return(df_list)
}
