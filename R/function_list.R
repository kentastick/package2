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
sig_val <- function(object = data, marker = "gene_list", use_func = "mean", filter = F) {
  gene_list <- get_list(marker)
  mt <- object@meta.data
  use_func <- switch (use_func, "mean" = mean, "gm_mean" = gm_mean1)
  count_mt <- object@assays$RNA@data
  gene_name <- rownames(count_mt)
  gene_list <- purrr::map(gene_list, ~.[. %in% gene_name])
  for(i in seq_along(gene_list)){
    sub_mt <- count_mt[gene_list[[i]],]
    value <- apply(sub_mt,2, use_func)
    mt[names(gene_list)[i]] <- value
  }
  mt <- mt[names(gene_list)]

  if(filter){
    val_mean <- apply(mt, 2, mean)
    for(i in seq_along(gene_list))
      temp <- mt[[names(gene_list)[i]]]
    mt[[names(gene_list)[i]]] <- if_else(temp> val_mean[[i]], temp, 0)
  }
  mt$cluster <- object$seurat_clusters
  mt$cell_id <- rownames(mt)
  return(mt)
}

# #sig_val + filter
# filter_cell <- function(marker = "gene_list",  use_func = "mean", object = data) {
#   gene_list <- get_list(marker = marker)
#   df <- sig_val(marker = marker, use_func = use_func, object = object)
#   val_mean <- apply(df, 2, mean)
#   for(i in seq_along(gene_list))
#     temp <- df[[names(gene_list)[i]]]
#   df[[names(gene_list)[i]]] <- if_else(temp> val_mean[[i]], temp, 0)
#   return(df)
# }


# geometric mean ----------------------------------------------------------

gm_mean1 = function(a){prod(a)^(1/length(a))}
gm_mean1(1:10)


gm_mean2 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


# signature plot ---------------------------------------------------------------------

#make df of each cluster's signature score: each_value/max_value, fraction_rate: expressed_cell(>0)/total_cells

sig_val2 <- function(score_mt) {
#  gene_list <- get_list(marker = marker)
#  if(filter){
#   val_mean <- apply(score_mt, 2, mean) # signature mean
#   for(i in seq_along(gene_list)){
#     temp <- score_mt[[names(gene_list)[i]]]
#     score_mt[[names(gene_list)[i]]] <- if_else(temp> val_mean[[i]], temp, 0)
#   }
# }
#  score_mt$cluster <- object@meta.data %>% .$seurat_clusters


  score_mt %>% gather(-cluster, -cell_id, key = "signature", value = "score") %>%
    group_by(signature, cluster) %>%
    summarise(fraction_of_cells = sum(score>0)/n(), mean = mean(score)) %>%
    group_by(signature) %>%
    mutate(max = max(mean)) %>%
    mutate(score = mean/max)
}



#internal function in make_subset
signature_plot_ <- function(mat_value, use.color = c("#0099FF", "#FAF5F5", "#E32020")) {
  mat_value %>% ggplot(aes(cluster, signature, colour =score, size = fraction_of_cells)) + geom_point() +
    scale_colour_gradientn(colours = c("red","yellow","white","lightblue","darkblue"),
                           values = c(1.0,0.7,0.6,0.4,0.3,0))
}


#conduct through calculation in sig_val1, 2 to plot
signature_plot <- function(object = data, marker = "gene_list", use_func = "mean",filter = F, use.color = c("#0099FF", "#FAF5F5", "#E32020")) {
    df <- sig_val(marker = marker, use_func = use_func, object = object, filter =filter)
    df <- sig_val2(score_mt = df)
    df %>% ggplot(aes(cluster, signature, colour =score, size = fraction_of_cells)) + geom_point() +
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



# pick_up_specific cell type ----------------------------------------------

  #execute at seurat_object directory
#make_subset(data_list = data_list, "HSC_combined", signature = "Mesenchyme", func = "me)

  make_subset <- function(data_list, save_folda, cell_type, use_func = "mean", marker = "gene_list") {

    #data_list <- data_list[!str_detect(data_list, "posi|blood")] #remove cd45posi(non-parenchyme cells include)

    #dir_name <- data_list %>% str_extract("(?<=\\/)\\S{1,12}(?=.rds)")
    data_name <- data_list %>% str_extract("\\S{1,20}(?=.rds)")
    #subset_name <- paste0(dir_name, "_hepato_subset") #hepatocyte extract

    dir_name <- file.path(save_folda, data_name)

    #subset_name <- paste0(data_name, "_hepato_subset") # extract
    subset_name <- paste0(data_name,"_", cell_type) # extract

    gene_list <- get_list(marker = marker)

    cell_num <- vector()
    for(i in seq_along(data_list)){

      cat(data_name[i], "executing\n")

      if(!dir.exists(dir_name[i])){
        dir.create(dir_name[i])
      }else next

      if(!data_name[i] %in% ls(envir = globalenv())) data <- readRDS(data_list[i])

      if(class(data) != "Seurat") next

       #signature_plot

      df <- sig_val(marker = marker, object = data, use_func = use_func)

      df2 <- sig_val2(score_mt = df)

      signature_plot_(mat_value = df2)

      try(ggsave(filename = paste0(dir_name[i], "/signature_plot.jpg")))

      #calculate mean value of each signature in whole cells.
      val_mean <- apply(df[names(gene_list)], 2, mean)

       # use_cluster_no <- df2 %>% group_by(cluster) %>% mutate(n_clu = row_number(-score)) %>%
       #      group_by(signature) %>% mutate(n_sig = row_number(-mean)) %>%
       #      filter(signature == cell_type, n_clu %in% c(1, 2)|n_sig ==1) %>% pull(cluster)
       #
       # use_id <- df %>% filter(get(cell_type)> val_mean[cell_type], cluster %in% use_cluster_no) %>%
       #   pull(cell_id)

       use_cluster_no <- df2 %>% group_by(cluster) %>% mutate(n_clu = row_number(-score)) %>%
            group_by(signature) %>% mutate(n_sig = row_number(-mean)) %>%
            filter(signature == cell_type, n_clu ==1|n_sig ==1|score>0.5) %>% pull(cluster)

       use_id <- df %>% filter(get(cell_type)> val_mean[cell_type], cluster %in% use_cluster_no) %>%
         pull(cell_id)

      if(length(use_id)<100){
        use_cluster_no <- df2 %>% group_by(cluster) %>% mutate(n_clu = row_number(-score)) %>%
          group_by(signature) %>% mutate(n_sig = row_number(-mean)) %>%
          filter(signature == cell_type, n_clu %in% c(1,2)|n_sig %in% c(1)) %>% pull(cluster)
        use_id <- df %>% filter(get(cell_type)> val_mean[cell_type], cluster %in% use_cluster_no) %>%
          pull(cell_id)

      }


      cat("number of cell is ", length(use_id), "\n")
      cell_num[i] <- length(use_id)
      names(cell_num)[i] <- data_name[i]
      if(length(use_id)< 20){
        cat("cell number is too small. skip procedure\n")
        next
      }

      #make subset_object of hepato_id cells
      sub_data <- subset(data, cells = use_id)
      ts(object = sub_data)
      try(ggsave(paste0(dir_name[i], "/subset_plot.jpg")))

      ts(object = data)
      try(ggsave(paste0(dir_name[i], "/whole_plot.jpg")))

      tmap(object = data, features =  gene_list[[cell_type]])
      try(ggsave(paste0(dir_name[i], "/whole_feature_plot.jpg"), width = 20, height =20, units =  "cm" ))

      tmap(object = sub_data, features =  gene_list[[cell_type]])
      try(ggsave(paste0(dir_name[i], "/feature_plot.jpg"), width = 20, height =20, units =  "cm" ))

      #save as a hepatocyte_subset object

      saveRDS(sub_data, file = paste0(dir_name[i], "/",subset_name[i],".rds"))
      saveRDS(cell_num, file = paste0(dir_name[i], "/",subset_name[i],"_cell_num.rds"))
      assign("cell_num", cell_num,envir = globalenv() )

    }

  }


  make_subset_ <- function(data_list, save_folda, cell_type, use_func = "mean", marker = "gene_list") {

    #data_list <- data_list[!str_detect(data_list, "posi|blood")] #remove cd45posi(non-parenchyme cells include)

    #dir_name <- data_list %>% str_extract("(?<=\\/)\\S{1,12}(?=.rds)")
    data_name <- data_list %>% str_extract("\\S{1,20}(?=.rds)")
    #subset_name <- paste0(dir_name, "_hepato_subset") #hepatocyte extract

    dir_name <- file.path(save_folda, data_name)

    #subset_name <- paste0(data_name, "_hepato_subset") # extract
    subset_name <- paste0(data_name,"_", paste0(cell_type, collapse = "_")) # extract

    gene_list <- get_list(marker = marker)

    cell_num <- vector()
    for(i in seq_along(data_list)){

      cat(data_name[i], "executing\n")

      if(!dir.exists(dir_name[i])){
        dir.create(dir_name[i])
      }else next

      if(!data_name[i] %in% ls(envir = globalenv())) data <- readRDS(data_list[i])

      if(class(data) != "Seurat") next

       #signature_plot

      df <- sig_val(marker = marker, object = data, use_func = use_func)

      df2 <- sig_val2(score_mt = df)

      signature_plot_(mat_value = df2)

      try(ggsave(filename = paste0(dir_name[i], "/signature_plot.jpg")))

      #calculate mean value of each signature in whole cells.
      val_mean <- apply(df[names(gene_list)], 2, mean)

       # use_cluster_no <- df2 %>% group_by(cluster) %>% mutate(n_clu = row_number(-score)) %>%
       #      group_by(signature) %>% mutate(n_sig = row_number(-mean)) %>%
       #      filter(signature == cell_type, n_clu %in% c(1, 2)|n_sig ==1) %>% pull(cluster)
       #
       # use_id <- df %>% filter(get(cell_type)> val_mean[cell_type], cluster %in% use_cluster_no) %>%
       #   pull(cell_id)

       use_cluster_no <- df2 %>% group_by(cluster) %>% mutate(n_clu = row_number(-score)) %>%
            group_by(signature) %>% mutate(n_sig = row_number(-mean)) %>%
            filter(signature %in% cell_type, n_clu ==1|n_sig ==1|score>0.5) %>% pull(cluster)

       use_id <- df %>% filter(cluster %in% use_cluster_no) %>%
         pull(cell_id)

      if(length(use_id)<100){
        use_cluster_no <- df2 %>% group_by(cluster) %>% mutate(n_clu = row_number(-score)) %>%
          group_by(signature) %>% mutate(n_sig = row_number(-mean)) %>%
          filter(signature == cell_type, n_clu %in% c(1,2)|n_sig %in% c(1)) %>% pull(cluster)
        use_id <- df %>% filter(cluster %in% use_cluster_no) %>%
          pull(cell_id)

      }


      cat("number of cell is ", length(use_id), "\n")
      cell_num[i] <- length(use_id)
      names(cell_num)[i] <- data_name[i]
      if(length(use_id)< 20){
        cat("cell number is too small. skip procedure\n")
        next
      }

      #make subset_object of hepato_id cells
      sub_data <- subset(data, cells = use_id)
      ts(object = sub_data)
      try(ggsave(paste0(dir_name[i], "/subset_plot.jpg")))

      ts(object = data)
      try(ggsave(paste0(dir_name[i], "/whole_plot.jpg")))

     for(j in seq_along(cell_type)){
          tmap(object = data, features =  gene_list[[cell_type[j]]])
          try(ggsave(paste0(dir_name[i], "/whole_feature_plot_", cell_type[j],".jpg"), width = 20, height =20, units =  "cm" ))
          tmap(object = sub_data, features =  gene_list[[cell_type[j]]])
          try(ggsave(paste0(dir_name[i], "/feature_plot", cell_type[j],".jpg"), width = 20, height =20, units =  "cm" ))
          }



      #save as a hepatocyte_subset object

      saveRDS(sub_data, file = paste0(dir_name[i], "/",subset_name[i],".rds"))
      saveRDS(cell_num, file = paste0(dir_name[i], "/",subset_name[i],"_cell_num.rds"))
      assign("cell_num", cell_num,envir = globalenv() )

    }

  }




# combined method ---------------------------------------------------------
#execute at each folda for subset ex. HSC_subset

combined <- function(object.list,cell_type = "subset", k.filter = 200) {
  object.list <- lapply(X = object.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })

  object.anchors <- FindIntegrationAnchors(object.list = object.list,k.filter = k.filter, dims = 1:20)
  #suceeded by changing object.list order but don't know the reason

  object.combined <- IntegrateData(anchorset = object.anchors, dims = 1:20)

  DefaultAssay(object.combined) <- "integrated"

  # Run the standard workflow for visualization and clustering
  object.combined <- ScaleData(object.combined, verbose = FALSE)
  object.combined <- RunPCA(object.combined, npcs = 30, verbose = FALSE)

  # t-SNE and Clustering
  object.combined <- RunUMAP(object.combined, reduction = "pca", dims = 1:20)
  object.combined <- RunTSNE(object.combined, reduction = "pca", dims = 1:20)
  object.combined <- FindNeighbors(object.combined, reduction = "pca", dims = 1:20)
  object.combined <- FindClusters(object.combined, resolution = 0.5, )

  DimPlot(object.combined, reduction = "umap", label = TRUE)

  ggsave("umap.jpg", device = "jpeg")
  DimPlot(object.combined, reduction = "tsne", label = TRUE)

  DimPlot(object.combined, reduction = "umap", group.by = "batch")
  object.combined$id <- colnames(object.combined)
  saveRDS(object.combined, file = paste0(cell_type, "_combined.rds"))
  return(object.combined)
}



# read all subset data ----------------------------------------------------

make_list <- function(cell_type) {
  file_dir <- list.files(pattern = paste0(paste0(cell_type,collapse = "_"),".rds"), recursive = T)

  file_name <- file_dir %>% str_split("/") %>% map(~.[1]) %>% unlist
  object.list <- list()
  for(i in seq_along(file_dir)){
    assign(x =file_name[i], readRDS(file_dir[i]))
    object.list[i] <- get(file_name[i])
  }
  return(object.list)
}


# simple save function: just write object name on the same directo --------

# sav <- function(x) {
#   parse_arg <- substitute(x)
#   if(is.symbol(parse_arg)){
#     parse_arg <- deparse(parse_arg)
#     }
#   temp <- get(parse_arg)
#   saveRDS(object = temp, file = paste0(parse_arg, ".rds"))
# }

sa <- function(x) {
  parse_arg <- substitute(x)
  saveRDS(object = data, file = paste0(deparse(parse_arg), ".rds"))
}


sav <- function(x) {
  parse_arg <- substitute(x)
  saveRDS(object = x, file = paste0(deparse(parse_arg), ".rds"))
}

af <- function(x) {
  log(x +1)
}

# get gene_list -----------------------------------------------------------

get_liver_marker <- function(n = 20) {
  liver_marker_list <- readRDS("~/single_cell/single_cell_project/gene_list/liver_marker_list.rds.rds")
  liver_marker_list <- map(liver_marker_list, ~head(., n))
}


save_list <- function(marker) {
  parse_name <- deparse(substitute(marker))
  saveRDS(marker, paste0("~/single_cell/single_cell_project/gene_list/", parse_name , ".rds"))
}


get_list <- function(marker, output = T) {
  get_list <- readRDS(paste0("~/single_cell/single_cell_project/gene_list/", marker,".rds"))
  if(output)assign(x = marker, value  = get_list, envir =globalenv())
}

get_list_name <- function() {
  list.files(path = "~/single_cell/single_cell_project/gene_list/")
}



get_list_ <- function(marker, n_liver_marker) {
  gene_list <- switch(marker,
                      "gene_list" = readRDS("~/single_cell/single_cell_project/gene_list/gene_list.rds"),
                      "original_list" = readRDS("~/single_cell/single_cell_project/gene_list/original_gene_list.rds"),
                      "macpoland" = get_liver_marker(n = n_liver_marker))
  return(gene_list)
}


#list to df

# signature = data.frame(signature = names(gene_list))
# temp <- vector()
# for(i in seq_along(gene_list)){
#   temp[i] <- gene_list[i]
# }
#
# signature$gene <- temp






# version up package2 -----------------------------------------------------
###
restart <- function(remotes, install_github) {
  remotes::install_github("kentastick/package2")
  detach("package:package2", unload=TRUE)
  library("package2", lib.loc="~/R/win-library/3.6")
}



# proportion plot -------------------------------------------------------------

pplot <- function(x) {
  #a <- deparse(substitute(x))
  data@meta.data %>% dplyr::select(seurat_clusters, x) %>%
    pivot_longer(x,values_to = "batch") %>%
    ggplot(aes(seurat_clusters, fill = batch)) + geom_bar(position = "fill")
}



# differential gene expression analysis within same sample ----------------

diff_test <- function(x, min.pct = 0.15, min.diff.pct = 0.1, logfc.threshold = 0.25, ...) {
  x <- substitute(x)
  if(!is.character(x)){
    x <- deparse(x)
  }
  batch_list <- as.character(unique(pull(data@meta.data, x)))
  marker<- tibble()
  for(i in seq_along(batch_list)){
    cat("executing ", batch_list[i], "process\n")
    use_id <- rownames(data@meta.data[pull(data@meta.data, x) == batch_list[i],])
    sub_temp <- SubsetData(data, cells = use_id)
    temp <- FindAllMarkers(object = sub_temp, min.pct = min.pct, min.diff.pct = min.diff.pct,
                           logfc.threshold = logfc.threshold, ...)
    temp$batch <- batch_list[i]
    marker <- marker %>% bind_rows(temp)
  }
  return(marker)
}




# plot signature scpre ----------------------------------------------------


feature_sig <- function(marker, object = data, use_func = "mean", filter = FALSE) {
  df <- sig_val(marker = marker, use_func = use_func , filter = filter)
  df <- df %>% dplyr::select(-cluster, -cell_id)
  object@meta.data <- object@meta.data %>% rownames_to_column(var = "temp") %>%
    bind_cols(df[colnames(df)]) %>% column_to_rownames(var = "temp")
  ump(object = object, features = colnames(df))
}



# add meta.data -----------------------------------------------------------


add_meta <- function(df, object = data) {
  object@meta.data <- object@meta.data %>% rownames_to_column(var = "temp") %>%
    bind_cols(df[colname(df)]) %>% column_to_rownames(var = "temp")
  return(object)
}



# cell_origin bar plot ----------------------------------------------------

bar_origin <- function(meta_data, object= data) {
  object@meta.data %>% dplyr::select(seurat_clusters, meta_data) %>%
    pivot_longer(meta_data,values_to = "batch") %>%
    ggplot(aes(seurat_clusters, fill = batch)) + geom_bar(position = "fill")
}





# adding meta info --------------------------------------------------------

add_info <- function(data) {
  data$reference <- data@meta.data %>% mutate(reference = fct_collapse(batch, Macparland = "macpoland",
                                                                       Aizarani = "aizarani",
                                                                       RamaChandran = c("chandran_cd45nega_ht", "chandran_cd45posi_ht", "chandran_cd45nega_ch", "chandran_cd45posi_ch", "ramachandran_blood"),
                                                                       ours_case1 = "pbc_case1",
                                                                       ours_case2 = "pbc_case2")) %>%
    mutate(reference = case_when(str_detect(batch, "fetal|adult")~"Segal",
                                 str_detect(batch, "HCC|ICC")~"Ma",
                                 TRUE~as.character(reference))) %>% pull(reference)

  data$disease <- data@meta.data %>% mutate(disease = fct_collapse(batch,NL_1 = "macpoland",
                                                                   NL_2 = "aizarani",
                                                                   NL_3 = c("chandran_cd45nega_ht", "chandran_cd45posi_ht"),
                                                                   CH = c("chandran_cd45nega_ch", "chandran_cd45posi_ch"),
                                                                   PBC = c("pbc_case1", "pbc_case2"),
                                                                   BL = "ramachandran_blood")) %>%
    mutate(disease = case_when(str_detect(batch, "fetal")~"FL",
                               str_detect(batch, "adult")~"NL_4",
                               str_detect(batch, "HCC")~"HCC",
                               str_detect(batch, "ICC")~"ICC",
                               TRUE~as.character(disease))) %>% pull(disease)

  return(data)


}


