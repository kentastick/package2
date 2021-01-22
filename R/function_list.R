
# data to list ------------------------------------------------------------

#' transform dataframe into list
#' @import Seurat tidyverse ggplot2


df_to_list <- function(df) {
  df_list <- vector("list", nrow(df))
  df <- df %>% as.matrix()
  for (i in 1:nrow(df)) {
    sub_row <- df[i,]
    names(sub_row) <- NULL
    df_list[[i]] <- sub_row
    names(df_list)[i] <- rownames(df)[i]
  }
  df_list <- purrr::map(df_list, ~.[!is.na(.)])
  return(df_list)
}


# do_seurat ---------------------------------------------------------------

#' create seurat object from raw matrix data or seurat_data which have not been executed thorogh
#' @export

do_seurat <- function(data) {
  if(!class(data) == "Seurat"){
    object <- CreateSeuratObject(counts = data, project = "object", min.cells = 3, min.features = 200)
  }else{
    object <- data
  }
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")

  #object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)

  all.genes <- rownames(object)

  #normalization data mean =1 variaty =0
  object <- ScaleData(object, features = all.genes)
  object <- RunPCA(object, features = VariableFeatures(object = object))

  print(object[["pca"]], dims = 1:5, nfeatures = 5)

  #jackstraw analysis
  object <- JackStraw(object, num.replicate = 10)
  object <- ScoreJackStraw(object, dims = 1:20)

  object <- FindNeighbors(object, dims = 1:27)
  object <- FindClusters(object, resolution = 0.8)
  object <- RunUMAP(object, dims = 1:30)
  object <- RunTSNE(object, dims = 1:30)
  return(object)
}



# tips  function --------------------------------------------------------


#' return present time as name
#' @export

make_time <- function() {
  Sys.time() %>% str_remove_all('[:punct:]|\\s')
}





# plot function -----------------------------------------------------------


#' wrapper function of Dimplot return umap plot
#' @export

up <- function(object = data, label= TRUE,...) {
  DimPlot(object = object, label = label, ...)
}


#' wrapper function of Dimplot return tsne plot
#' @export

ts <- function(object = data, ...) {
  DimPlot(object = object, reduction = 'tsne', label = TRUE, ...)
}

#' wrapper function of FeaturePlot tsne
#' @export

tmap <- function(features, object = data, ...) {
  FeaturePlot(features = features,object = object, reduction = 'tsne',  cols = c('lightgray', 'red'), min.cutoff = 0, ...)
}

#' wrapper function of featureplot umap
#' @export

ump <- function(features, object = data,...) {
  FeaturePlot(features = features,object = object, cols = c('lightgray', 'red'), min.cutoff = 0,...)
}

#' violin plot
#' @export

vl <- function(features, object = data,...) {
  Seurat::VlnPlot(object = object, features = features, ...)
}


# get gene data -----------------------------------------------------------

#' search gene expressed in seurat_data
#' @export

search_gene <- function(data = data) {
  data@assays$RNA@data@Dimnames[[1]] %>% View
}

#' search gene expressed in seurat_data
#' @export

pick_gene <- function(pattern, data = data) {
  data@assays$RNA@data@Dimnames[[1]] %>% str_subset(., pattern =pattern)
}



# signature plot ---------------------------------------------------------------------

#' make df of each cluster's signature score: each_value/max_value, fraction_rate: expressed_cell(>0)/total_cells

#' @export

gm_mean1 = function(a){prod(a)^(1/length(a))}

gm_mean2 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


#' calculate module score
#' @export

sig_val <- function(object = data, marker = "gene_list", use_func = "mean", add_id_cluster = T,filter = F) {
  gene_list <- get_list(marker)
  mt <- object@meta.data
  use_func <- switch (use_func, "mean" = mean, "gm_mean" = gm_mean1)
  count_mt <- object@assays$RNA@data
  gene_name <- rownames(count_mt)
  gene_list <- purrr::map(gene_list, ~.[. %in% gene_name])
  gene_list <- gene_list[map(gene_list, length)>1]
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
  if(add_id_cluster){
    #mt$cluster <- object@meta.data[, label_name]
    mt$cluster <- object@active.ident
    mt$cell_id <- rownames(mt)
  }
  return(mt)
}

sig_val_ <- function(object = data, marker = "gene_list", use_func = "mean", add_id_cluster = T,filter = F) {
  gene_list <- get_list(marker)
  mt <- object@meta.data
  use_func <- switch (use_func, "mean" = mean, "gm_mean" = gm_mean1)
  count_mt <- object@assays$RNA@data
  gene_name <- rownames(count_mt)
  gene_list <- purrr::map(gene_list, ~.[. %in% gene_name])
  gene_list <- gene_list[map(gene_list, length)>1]
  for(i in seq_along(gene_list)){
    sub_mt <- count_mt[gene_list[[i]],]
    value <- apply(sub_mt,2, use_func)
    mt[names(gene_list)[i]] <- value
  }
  mt <- mt[names(gene_list)]
  mt <- mt %>% mutate_if(.predicate = ~is.numeric(.), .funs = ~{a <- max(.); ./a})

  if(filter){
    val_mean <- apply(mt, 2, mean)
    for(i in seq_along(gene_list))
      temp <- mt[[names(gene_list)[i]]]
    mt[[names(gene_list)[i]]] <- if_else(temp> val_mean[[i]], temp, 0)
  }
  if(add_id_cluster){
    #mt$cluster <- object@meta.data[, label_name]
    mt$cluster <- object@active.ident
    mt$cell_id <- rownames(mt)
  }
  return(mt)
}



sig_val2 <- function(score_mt) {
  score_mt %>% gather(-cluster, -cell_id, key = "signature", value = "score") %>%
    group_by(signature, cluster) %>%
    summarise(fraction_of_cells = sum(score>0)/n(), mean = mean(score)) %>%
    group_by(signature) %>%
    mutate(max = max(mean)) %>%
    mutate(score = mean/max)
}



#calculate whithin each gene score normalized by gene across the cluster
signature_plot <- function(marker = "gene_list", object = data, use_func = "mean", filter = F, use.color = c("#0099FF", "#FAF5F5", "#E32020")) {
  stopifnot(is.character(marker))
  stopifnot(class(object) == "Seurat")

  df <- sig_val(marker = marker, use_func = use_func, object = object, filter =filter)
    df <- sig_val2(score_mt = df)
    df %>% ggplot(aes(cluster, signature, colour =score, size = fraction_of_cells)) + geom_point() +
    scale_colour_gradientn(colours = c("red","yellow","white","lightblue","darkblue"),
                           values = c(1.0,0.7,0.6,0.4,0.3,0))+
      theme(axis.text.x = element_text(angle = 90))
}


signature_tile <- function(marker = "gene_list", object = data,title = "", use_func = "mean",filter = F, use.color = c("#0099FF", "#FAF5F5", "#E32020")) {
    df <- sig_val(marker = marker, use_func = use_func, object = object, filter =filter)
    df <- sig_val2(score_mt = df)
    df %>% ggplot(aes(cluster, signature, fill = score)) + geom_tile(color = "black") +
    #scale_fill_gradient2(low = "blue",  mid = "white", high = "red", midpoint = 0.5)
     scale_fill_gradientn(colours = c("red","yellow","white","lightblue","darkblue"),
                            values = c(1.0,0.7,0.6,0.4,0.3,0)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.2), axis.text.y = element_text(size = 12)) +
      labs(title = title)
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
    subset_name <- paste0(data_name,"_", paste0(cell_type,collapse = "")) # extract

    gene_list <- get_list(marker = marker)

    cell_num <- vector()
    for(i in seq_along(data_list)){

      cat(data_name[i], "executing\n")

      if(!dir.exists(dir_name[i])){
        dir.create(dir_name[i])
      }else next

      if(!data_name[i] %in% ls(envir = globalenv())){
        data <- readRDS(data_list[i])
      }else{
        data <- get(data_name[i])
      }

      if(class(data) != "Seurat") next

       #signature_plot

      df <- sig_val(marker = marker, object = data, use_func = use_func)

      df2 <- sig_val2(score_mt = df)

      signature_plot_(mat_value = df2)

      try(ggsave(filename = paste0(dir_name[i], "/signature_plot.jpg")))

      #calculate mean value of each signature in whole cells.
      #val_mean <- apply(df[names(gene_list)], 2, mean)


       # use_cluster_no <- df2 %>% group_by(cluster) %>% mutate(n_clu = row_number(-score)) %>%
       #      group_by(signature) %>% mutate(n_sig = row_number(-mean)) %>%
       #      filter(signature == cell_type, n_clu %in% c(1, 2)|n_sig ==1) %>% pull(cluster)
       #
       # use_id <- df %>% filter(get(cell_type)> val_mean[cell_type], cluster %in% use_cluster_no) %>%
       #   pull(cell_id)

       use_cluster_no <- df2 %>% group_by(cluster) %>% mutate(n_clu = row_number(-score)) %>%
            group_by(signature) %>% mutate(n_sig = row_number(-mean)) %>%
            filter(signature %in% cell_type, n_clu ==1|n_sig ==1,score>0.5) %>% pull(cluster)

       use_id <- df %>% filter(cluster %in% use_cluster_no) %>% pull(cell_id)
       # for(j in seq_along(cell_type)){
       #   df <- df %>% filter(get(cell_type[j])> val_mean[cell_type[j]])
       # }


      # if(length(use_id)<100){
      #   use_cluster_no <- df2 %>% group_by(cluster) %>% mutate(n_clu = row_number(-score)) %>%
      #     group_by(signature) %>% mutate(n_sig = row_number(-mean)) %>%
      #     filter(signature == cell_type, n_clu %in% c(1,2)|n_sig %in% c(1)) %>% pull(cluster)
      #   use_id <- pull(df, cell_id)
      # }

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
      try(ggsave(filename = paste0(dir_name[i], "/subset_plot.jpg")))

      ts(object = data)
      try(ggsave(filename = paste0(dir_name[i], "/whole_plot.jpg")))

      tmap(object = data, features =  unlist(gene_list[cell_type]))
      try(ggsave(filename = paste0(dir_name[i], "/whole_feature_plot.jpg"), width = 20, height =20, units =  "cm" ))

      tmap(object = sub_data, features =  unlist(gene_list[cell_type]))
      try(ggsave(filename = paste0(dir_name[i], "/feature_plot.jpg"), width = 20, height =20, units =  "cm" ))

      #save as a hepatocyte_subset object

      saveRDS(sub_data, file = paste0(dir_name[i], "/",subset_name[i],".rds"))
      saveRDS(cell_num, file = paste0(dir_name[i], "/",subset_name[i],"_cell_num.rds"))
      assign("cell_num", cell_num,envir = globalenv() )

    }

  }


make_subset_id <- function(data = data, pathname = "test") {
  file_list <- list.files(pattern = ".rds")
  id_list <- data@meta.data %>% dplyr::select(id, data_origin) %>% split(data$data_origin) %>%
    map(~.$id)
  id_list <- id_list %>% map(~str_remove_all(., "_\\d{1,2}$"))
  if(!dir.exists(pathname)) dir.create(pathname)
  name_list <- names(id_list)
  for(i in seq_along(id_list)){
    use_id <- id_list[[i]]
    temp <- readRDS(file = str_subset(file_list, name_list[i]))
    if(class(temp) != "Seurat") next
    temp <- subset(temp, cells = use_id)
    saveRDS(temp, file = file.path(pathname, paste0(name_list[i],"_subset.rds" ) ))
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

  # object.anchors <- try(FindIntegrationAnchors(object.list = object.list,k.filter = k.filter, dims = 1:20))
  # if(class(object.anchors)== "try-error")return(NULL)
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

make_list <- function() {
  file_dir <- list.files(pattern = ".rds", recursive = T)

  file_name <- file_dir %>% str_split("/") %>% map(~.[1]) %>% unlist
  object.list <- list()
  for(i in seq_along(file_dir)){
    assign(x =file_name[i], readRDS(file_dir[i]))
    object.list[i] <- get(file_name[i])
  }
  return(object.list)
}

load_list <- function(data_list) {
  file_name <-  data_list %>% str_extract(".*(?=.rds)")
  for(i in seq_along(data_list)){
    assign(x =file_name[i], value = readRDS(data_list[i]), envir = globalenv())
  }
}


# save function -----------------------------------------------------------


sav <- function(x) {
  parse_arg <- substitute(x)
  saveRDS(object = x, file = paste0(deparse(parse_arg), ".rds"))
}


# marker_list functions  -----------------------------------------------------------

save_list <- function(marker) {
  parse_name <- deparse(substitute(marker))
  gene_list_path = c(#"~/single_cell/single_cell_project/gene_list/",
                     #"E:/single_cell_project/gene_list/",
                     "~/single_cell/package2/test/gene_list/",
                     "E:/single_cell_project/package2/test/gene_list/",
                     "./test/gene_list/",
                     "~/package2/test/gene_list/"
                     )
  for (i in gene_list_path){
    try(saveRDS(marker, paste0(i, parse_name , ".rds")), silent = T)
  }

}


get_list <- function(marker, output = T) {
  gene_list_path = c(#"~/single_cell/single_cell_project/gene_list/",
                     #"E:/single_cell_project/gene_list/",
                     "~/single_cell/package2/test/gene_list/",
                     "E:/single_cell_project/package2/test/gene_list/",
                     "./test/gene_list/",
                     "~/package2/test/gene_list/"
                     )
  for (i in gene_list_path){
    get_list <- try(readRDS(paste0(i, marker,".rds")), silent = T)
  if(class(get_list) != "try-error")break
}
  if(output)assign(x = marker, value  = get_list, envir =globalenv())
}

get_list_name <- function() {
  gene_list_path = c(#"~/single_cell/single_cell_project/gene_list/",
                     #"E:/single_cell_project/gene_list/",
                     "~/single_cell/package2/test/gene_list/",
                     "E:/single_cell_project/package2/test/gene_list/",
                     "./test/gene_list/",
                     "~/package2/test/gene_list/"
                     )
  for (i in gene_list_path){
    res <- try(list.files(path = i))
    if(length(res) != 0) break
  }
  res
}





# differential gene expression analysis within same sample ----------------


#do findmarker per each selected condition
diff_test_batch <- function(x, object = data, min.pct = 0.1, min.diff.pct = 0.1, logfc.threshold = 0.25, only.pos =T, ...) {
  x <- substitute(x)
  if(!is.character(x)){
    x <- deparse(x)
  }
  batch_list <- as.character(unique(pull(object@meta.data, x)))
  marker<- tibble()
  for(i in seq_along(batch_list)){
    cat("executing ", batch_list[i], "process\n")
    use_id <- rownames(object@meta.data[pull(object@meta.data, x) == batch_list[i],])
    sub_temp <- SubsetData(object, cells = use_id)
    temp <- FindAllMarkers(object = sub_temp, min.pct = min.pct, min.diff.pct = min.diff.pct,
                           logfc.threshold = logfc.threshold, only.pos = only.pos,  ...)

    if(nrow(temp)==0)next

    temp$batch <- batch_list[i]
    marker <- marker %>% bind_rows(temp)

  }
  return(marker)
}



diff_test_batch_2 <- function (x, object = data, min.pct = 0.1, min.diff.pct = 0.1,
          logfc.threshold = 0.25, only.pos = T, ...)
{
  x <- substitute(x)
  if (!is.character(x)) {
    x <- deparse(x)
  }
  batch_list <- as.character(unique(pull(object@meta.data,
                                         x)))
  marker <- tibble()
  for (i in seq_along(batch_list)) {
    cat("executing ", batch_list[i], "process\n")
    use_id <- rownames(object@meta.data[pull(object@meta.data,
                                             x) == batch_list[i], ])
    sub_temp <- SubsetData(object, cells = use_id)
    temp <- FindAllMarkers(object = sub_temp, min.pct = min.pct,
                           min.diff.pct = min.diff.pct, logfc.threshold = logfc.threshold,
                           only.pos = only.pos, ...)
    if (nrow(temp) == 0)
      next
    temp$batch <- batch_list[i]
    marker <- marker %>% bind_rows(temp)
  }
  return(marker)
}

#diff_test


fil_marker <- function(marker_list, n = 30, filter = F, remove_gene = "^MT|^RPS|^RPL|^RBP|^LOC|^C1of|^RNA") {

  marker_list <- marker_list %>% group_by(cluster) %>% top_n(n, avg_logFC)

  if(filter){
    marker_list <- marker_list %>% filter(!str_detect(gene, remove_gene))
  }
  return(marker_list)
}






# cell_origin bar plot ----------------------------------------------------

bar_origin <- function(bar_x, bar_y,object= data,angle = 60, position = "fill", randam = T) {
  bar_x <- enquo(bar_x)
  bar_y <- enquo(bar_y)
  p <- object@meta.data %>% dplyr::select(!!bar_x, !!bar_y) %>%
    ggplot(aes(!!bar_x, fill = !!bar_y)) + geom_bar(position = position) +
    theme(axis.text.x = element_text(angle = angle, vjust = 0.5,size = rel(1.5))) + labs(fill = "Cell_type") +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) + labs(y = "ratio") +
    theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA, colour = "black"))
  if(randam){
    n_color <- object@meta.data %>% dplyr::select(!!bar_y) %>% n_distinct()
   p <- p + scale_fill_manual(values = gg_color_hue(n = n_color, randam = T))
  }
  print(p)
}

# meta data modifying function--------------------------------------------------------

#change current cell label
id_ch <- function(use_meta, object = data, x) {
  object_name <- as.character(substitute(object))
  Idents(object = object) <- use_meta
  assign(x = object_name, value = object, envir = .GlobalEnv)
}


#change reference sample_type infomation
add_info <- function(data) {
  data$reference <- data@meta.data %>% mutate(reference = fct_collapse(batch, Macparland = "macpoland",
                                                                       Aizarani = "aizarani",
                                                                       RamaChandran = c("chandran_cd45nega_ht", "chandran_cd45posi_ht", "chandran_cd45nega_ch", "chandran_cd45posi_ch", "ramachandran_blood"),
                                                                       ours_case1 = "pbc_case1",
                                                                       ours_case2 = "pbc_case2")) %>%
    mutate(reference = case_when(str_detect(batch, "fetal|adult")~"Segal",
                                 str_detect(batch, "HCC|ICC")~"Ma",
                                 TRUE~as.character(reference))) %>% pull(reference)
data$cancer


  data$disease <- data@meta.data %>% mutate(disease = fct_collapse(batch,NL_1 = "macpoland",
                                                                   NL_2 = "aizarani",
                                                                   NL_3 = c("chandran_cd45nega_ht", "chandran_cd45posi_ht"),
                                                                   CH = c("chandran_cd45nega_ch", "chandran_cd45posi_ch"),
                                                                   PBC_1 = "pbc_case1",
                                                                   PBC_2 = "pbc_case2",
                                                                   BL = "ramachandran_blood")) %>%
    mutate(disease = case_when(str_detect(batch, "fetal")~"FL",
                               str_detect(batch, "adult")~"NL_4",
                               str_detect(batch, "HCC")~"HCC",
                               str_detect(batch, "ICC")~"ICC",
                               TRUE~as.character(disease))) %>% pull(disease)

  data$condition  <- data@meta.data %>%
    mutate(condition = fct_collapse(disease, Healthy = c("NL_1", "NL2", "NL_3"),
                                    CH = "CH", PBC = c("PBC_1", "PBC_2")
                                    )) %>%
    pull(condition)

  data$disease2 <- data[[]] %>% mutate(disease = case_when(cancer=="non_tumor"~"NL_5",
                                                           cancer %in% c("suppliment", "tumor")~"ICC_2",
                                                           TRUE~as.character(disease))) %>% pull(disease) %>%
    fct_relevel(c("PBC_1", "PBC_2", paste0("NL_", 1:5), "CH", "ICC_1", "ICC_2", "HCC", "FL", "BL" ))



  data$data_origin <- data@meta.data %>% mutate(n = fct_collapse(batch, macparland = c("macpoland"), segal = c("adult", "fetal"))) %>%
    pull(n)

  data$id <- rownames(data@meta.data)
  return(data)

}


add_sig_val <- function(object = data, marker_list, use_func = "mean", label_name = "label", overwrite = T, add_signature_label = T){

  #object <- eval(as.name(object_name))
  object_name <- as.character(substitute(object))

  df_list <- vector("list", length(marker_list))

  for(i in seq_along(marker_list)){
    df_list[[i]] <- sig_val(object = object, marker = marker_list[i], use_func = use_func) %>%
      add_m(add = switch(use_func, "mean" = "", "gm_mean" = "_gm"))
  }

  df_com <- df_list %>% purrr::reduce(cbind)

  df_com <- df_com %>% keep(is.numeric)

  object <- add_meta(df = df_com, object = object)

  if(add_signature_label){

    for(i in seq_along(marker_list)){
    use_list <- get_list(marker = marker_list[[i]])
    if(use_func == "gm_mean"){
      use_name <- paste0(names(use_list), "_gm")
    }else {
      use_name <- names(use_list)

    }
    rank_table <- object@meta.data %>%  dplyr::select(seurat_clusters, all_of(use_name)) %>%
      group_by(seurat_clusters) %>% gather(-seurat_clusters, key = "signature", value = "value") %>%
      group_by(seurat_clusters, signature) %>%
      summarise(m = mean(value)) %>%
      group_by(signature) %>%
      mutate(max = max(m), score = m/max) %>%
      group_by(seurat_clusters) %>%
      mutate(rank = row_number(score)) %>%
      filter(rank == 1) %>%
      dplyr::select(seurat_clusters, signature)
    use_label <- object@meta.data %>% left_join(rank_table, by = "seurat_clusters") %>% pull(signature)
    object@meta.data[, marker_list[i]] <- use_label
    }
  }
  if(overwrite){
    assign(x = object_name, value = object, envir = .GlobalEnv)
  } else return(object)

}
add_sig_val_ <- function(object = data, marker_list, use_func = "mean", label_name = "label", overwrite = T, add_signature_label = T){

  #object <- eval(as.name(object_name))
  object_name <- as.character(substitute(object))

  df_list <- vector("list", length(marker_list))

  for(i in seq_along(marker_list)){
    df_list[[i]] <- sig_val_(object = object, marker = marker_list[i], use_func = use_func) %>%
      add_m(add = switch(use_func, "mean" = "", "gm_mean" = "_gm"))
  }

  df_com <- df_list %>% purrr::reduce(cbind)

  df_com <- df_com %>% keep(is.numeric)

  object <- add_meta(df = df_com, object = object)

  if(add_signature_label){

    for(i in seq_along(marker_list)){
    use_list <- get_list(marker = marker_list[[i]])
    if(use_func == "gm_mean"){
      use_name <- paste0(names(use_list), "_gm")
    }else {
      use_name <- names(use_list)

    }
    rank_table <- object@meta.data %>%  dplyr::select(seurat_clusters, all_of(use_name)) %>%
      group_by(seurat_clusters) %>% gather(-seurat_clusters, key = "signature", value = "value") %>%
      group_by(seurat_clusters, signature) %>%
      summarise(m = mean(value)) %>%
      group_by(signature) %>%
      mutate(max = max(m), score = m/max) %>%
      group_by(seurat_clusters) %>%
      mutate(rank = row_number(score)) %>%
      filter(rank == 1) %>%
      dplyr::select(seurat_clusters, signature)
    use_label <- object@meta.data %>% left_join(rank_table, by = "seurat_clusters") %>% pull(signature)
    object@meta.data[, marker_list[i]] <- use_label
    }
  }
  if(overwrite){
    assign(x = object_name, value = object, envir = .GlobalEnv)
  } else return(object)

}



#add strong/weak value of one gene by specific value
add_meta_binval <- function(gene, object = data) {
  object_name <- as.character(substitute(object))
  df <- FetchData(object = object, var = gene) %>% rownames_to_column(var = "id")
  gene_m <- df %>% pull(gene) %>%  mean()
  gene_sd <- df %>% pull(gene) %>%  sd()
  cut_off <- gene_m + 2*gene_sd
  colnames(df)[2] <- "temp"
  use_column <- df %>% dplyr::select(id, temp) %>% mutate(temp= if_else(temp> cut_off, "strong", "ordinary")) %>% dplyr::select(temp)
  colnames(use_column) <- paste0(gene,"_bin")
  object <- add_meta(df = use_column, object = object)
  assign(x = object_name, value = object, envir = .GlobalEnv)
}




# subset filter---------------------------------------------------------

sub_fil <- function(object = data, ...) {
  object$id <- colnames(object)
  use_id <- object@meta.data %>% filter(...) %>% pull(id)
  res <- Seurat:::subset.Seurat(x = object, cells = use_id)
  return(res)
}


# gene annotation analysis ------------------------------------------------

convert_gene <- function(x) {
  library(org.Hs.eg.db)
  res <- try(clusterProfiler::bitr(x, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db"))

  return(res)
}

convert_gene_to_ENTREZ <- function(x) {
  library(org.Hs.eg.db)
  res <-
    try(clusterProfiler::bitr(x,
                              fromType = "SYMBOL",
                              toType = c("ENTREZID"),
                              OrgDb = "org.Hs.eg.db"))
  return(res)
}
convert_gene_to_SYMBOL <- function(x) {
  library(org.Hs.eg.db)
  res <-
    try(clusterProfiler::bitr(x,
                              fromType = "ENTREZID",
                              toType = c("SYMBOL"),
                              OrgDb = "org.Hs.eg.db"))
  return(res)
}




convertMouseGene <- function(x){
  mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  genesV2 = biomaRt::getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
}




do_enrich_msig <- function(use_gene, name) {
  m_df <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, gene_symbol)
  em <- clusterProfiler::enricher(use_gene, TERM2GENE=m_df)
  res <- em@result
  return(res)
}


do_enrich_kegg <- function(use_gene) {
  tbl <- convert_gene_to_ENTREZ(use_gene)
  tbl_cnv <- tbl$SYMBOL
  names(tbl_cnv) <- tbl$ENTREZID
  res <- clusterProfiler::enrichKEGG(gene = tbl$ENTREZID,organism= 'hsa',pvalueCutoff = 0.05)
  res <- res@result %>% mutate(geneID= map(geneID, ~str_split(.x[[1]], pattern = "/") %>% map(~tbl_cnv[.] %>% as.character) %>% unlist %>% paste0(collapse = "/") %>% as.character) )
  res %>% mutate(geneID = map(geneID, ~.x[[1]]) %>% as.character()) %>% return()
}


do_enrich_reactome <- function(use_gene) {
  tbl <- convert_gene_to_ENTREZ(use_gene)
  tbl_cnv <- tbl$SYMBOL
  names(tbl_cnv) <- tbl$ENTREZID
  res <- ReactomePA::enrichPathway(gene = tbl$ENTREZID,pvalueCutoff = 0.05)
  res <- res@result %>% mutate(geneID= map(geneID, ~str_split(.x[[1]], pattern = "/") %>% map(~tbl_cnv[.] %>% as.character) %>% unlist %>% paste0(collapse = "/") %>% as.character) )
  res %>% mutate(geneID = map(geneID, ~.x[[1]]) %>% as.character()) %>% return()
}


do_enrich_go <- function(use_gene, ont = "BP") {
  tbl <- convert_gene_to_ENTREZ(use_gene)
  tbl_cnv <- tbl$SYMBOL
  names(tbl_cnv) <- tbl$ENTREZID
  res <- clusterProfiler::enrichGO(gene          = tbl$ENTREZID,
                                   OrgDb         = org.Hs.eg.db,
                                   ont           = ont,
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)
  res <- res@result #%>% mutate(geneID= map(geneID, ~str_split(.x[[1]], pattern = "/") %>% map(~tbl_cnv[.] %>% as.character) %>% unlist %>% paste0(collapse = "/") %>% as.character) )
  res %>% mutate(geneID = map(geneID, ~.x[[1]]) %>% as.character()) %>% return()
}


do_david <- function(use_gene, listName = "gene") {

  david<- RDAVIDWebService::DAVIDWebService(email="ktaka@stu.kanazawa-u.ac.jp", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

  gene_df <- convert_gene(x = use_gene)
  gene_tbl <- symbol_entrez_conv(object = gene_df)


  result<-RDAVIDWebService::addList(david, gene_df$ENTREZID,
                                    idType="ENTREZ_GENE_ID",
                                    listName= listName, listType="Gene")

  termCluster<- RDAVIDWebService::getClusterReport(david, type="Term")

  termCluster <- conv_david(termCluster, use_gene_tbl = gene_tbl)

  return(termCluster)
}

symbol_entrez_conv <- function(object) {
  temp <- object$SYMBOL
  names(temp) <- object$ENTREZID
  return(temp)
}

#' convert david result to tidy dataframe
#' @export

conv_david <- function(david_object, use_gene_tbl) {
  cluster_list <-map(david_object@cluster, ~.[[2]] %>% asS3)
  score_list <-map(david_object@cluster, ~.[[1]])
  length_res <- length(david_object@cluster)

  use_list <- list(x = cluster_list, y = 1:length_res,z = score_list)

  use_func <- function(x,y,z) mutate(.data = asS3(x), cluster = paste0("cluster_", y), score = z)
  result <- pmap(use_list, .f = use_func)

  result <- purrr::reduce(result, rbind)

  result <- result %>% mutate(Genes = str_split(Genes, pattern = ", "))

  result <- result %>% mutate(Genes = map(Genes, ~as.character(use_gene_tbl[.]) %>% paste0(collapse = ", ")) %>% unlist())

  return(result)
}

#' return david input gene information

david_gene_report <- function(use_gene) {


  david<- RDAVIDWebService::DAVIDWebService(email="ktaka@stu.kanazawa-u.ac.jp",
                                            url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

  gene_df <- convert_gene(x = use_gene)
  gene_tbl <- symbol_entrez_conv(object = gene_df)

  result<-RDAVIDWebService::addList(david, gene_df$ENTREZID,
                                    idType="ENTREZ_GENE_ID",
                                    listName= "gene", listType="Gene")

  gene_report <- RDAVIDWebService::getGeneListReport(david)
  rerun(gene_report)

}





#' do pathway analysis of several
#' @export

do_pathway <- function(data) {
  #stopifnot(all(class(data) == c("grouped_df", "tbl_df", "tbl", "data.frame")))
  cat("------------msig_procedure--------\n")
  data <- data %>% mutate(enrich_h = map(gene, ~try(do_enrich_msig(.x))))
  cat("------------kegg_procedure--------\n")
  data <- data %>% mutate(enrich_kegg = map(gene, ~try(do_enrich_kegg(.x))))
  cat("------------reactome_procedure--------\n")
  data <- data %>% mutate(enrich_reactome = map(gene, ~try(do_enrich_reactome(.x))))
  cat("------------david_procedure--------\n")
  data <- data %>% mutate(enrich_david = map(gene, ~try(do_david(.x))))
  cat("------------go_bp_procedure--------\n")
  data <- data %>% mutate(res_go_bp = map(gene, ~try(do_enrich_go(.x))))
  cat("------------go_cc_procedure--------\n")
  data <- data %>% mutate(res_go_cc = map(gene, ~try(do_enrich_go(.x, ont = "CC"))))
  return(data)
}


# venn  -------------------------------------------------------------------


do_venn <- function(arg1, arg2) {
  if(!dir.exists("batch_diff")) dir.create("batch_diff")
  jpeg(paste0("batch_diff/",as.character(arg2), "_venn.jpg"), width = 700, height = 700)
  venn::venn(x = arg1, zcolor = "style", ilcs = 2)
  dev.off()
  #ggsave( filename =  paste0("batch_diff/",as.character(arg2), "_venn.jpg"), device = "jpg")
}


make_venn <- function(df = dirr_test_res) {
  batch_diff_test <- df %>%
    group_by(batch, cluster) %>% nest() %>%
    mutate(data_filtered = map(data, ~filter(.) %>% pull(gene))) %>%
    group_by(cluster) %>% nest() %>%
    mutate(gene_list = map(data, ~zip(df = .))) %>%
    mutate(venn = map2(gene_list,cluster, ~do_venn(.x,.y) ))
  sav(batch_diff_test)
  return(batch_diff_test)

}



# tile_plot ---------------------------------------------------------------

fil_gene <- function(gene, object) {
  if(class(gene) =="list"){
    gene <- gene %>% map(., ~.[. %in% rownames(object)])
    gene <- gene[map(gene, length)>0]
  }else{
    gene <- gene[gene %in% rownames(object)]
  }
  return(gene)
}


remove_list_dup <- function(gene_list) {
  gene_list <- gene_list %>% unlist() #%>% split(., ngene_listmes(.))
  gene_list <- gene_list[!duplicated(gene_list)]
  names(gene_list) <- names(gene_list) %>% str_remove("\\d{1,2}$")
  gene_list %>% split(., names(.)) %>% map(as.vector)
}




tile_plot <- function(gene = "gene_list", object = data, title = "", order = F, plot_wrap = F, fil_val= NULL, color_label = T, ...) {
  stopifnot(class(gene) != "Seurat")
  if(str_detect(gene, "_list")){
    gene <- get_list(gene)
  }
  DefaultAssay(object = object) <- "RNA"
  if(class(gene) == "list" ){
    gene <- remove_list_dup(gene)
    gene <-   fil_gene(gene, object = object)
    feature <- unlist(gene)
    label_df <- enframe(gene, name = "label",value = "gene") %>% unnest

    for_tile_legend_df <<- label_df
  } else{
    feature <- fil_gene(gene, object = object)
  }
  use_id <- colnames(object)
  use_df <- object@assays$RNA@data[feature, use_id]


  cluster_label <- object@active.ident
  if(length(feature) ==1){
    use_df <- use_df %>% as.tibble()
  }else{
    use_df <- t(as.matrix(use_df)) %>% as.tibble()
  }
   use_df<- use_df %>% add_column(cluster = cluster_label)
  use_df <- use_df %>% tidyr::pivot_longer(-cluster, names_to = "gene", values_to = "logCPM") %>%
    group_by(cluster, gene) %>% summarise(avg_logCPM = mean(logCPM), pct = sum(logCPM>0)/n()) %>%
    group_by(gene) %>% mutate(score = avg_logCPM/max(avg_logCPM), m = mean(avg_logCPM))

  if(!is.null(fil_val)) {
    use_df <- use_df %>% filter(m > fil_val)
    use_gene <- unique(use_df$gene)
    label_df <- filter(gene %in% use_gene)

  }

  if(class(gene) == "list" ){
    use_df <- use_df %>% left_join(label_df, by = c("gene"))

    use_df <- use_df %>% ungroup() %>%  mutate(gene = fct_relevel(gene, feature))

  }
  if(order){
    use_df <- use_df %>% mutate( cluster = fct_reorder(cluster, avg_logCPM))
  }



  p <- use_df %>% ggplot(aes(cluster, gene, size = pct, colour = score)) + geom_point() +
    scale_colour_gradientn(colours = c("red","yellow","white","lightblue","darkblue"),
                           values = c(1.0,0.7,0.6,0.4,0.3,0))

  if(class(gene) == "list"){
    if(plot_wrap){
      p <- p + facet_grid(~label)
    }
    if(color_label){
      n <- length(gene)
      label_color <- gg_color_hue(n)
      label_color_use <- label_color[as.numeric(plyr::mapvalues(label_df$label, from = unique(label_df$label), to = 1:n))]

      use_df <- use_df %>% mutate(label_color = label_color[as.numeric(as.factor(label))])
      p <- use_df %>% ggplot(aes(cluster, gene, size = pct, fill = label, colour = score)) + geom_point() +
        scale_colour_gradientn(colours = c("red","yellow","white","lightblue","darkblue"),
                               values = c(1.0,0.7,0.6,0.4,0.3,0))
      p <- p + theme(axis.text.y = element_text(colour = label_color_use),
                     axis.text.x = element_text(angle = 90, vjust = 0.5)

                     ) +
        scale_fill_manual(values = label_color) + labs(title = title)
      return(p)
    }


  }else return(p +theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + labs(title = title))

}






tile_legend <- function(df = for_tile_legend_df) {
  n <- length(unique(df$label))
  df %>% mutate(gene_label = fct_relevel(label, unique(df$label)),color = as.numeric(as.factor(label))) %>%
    ggplot(aes(gene_label, color, fill = gene_label)) +
    geom_bar(stat = "identity")+ scale_fill_manual(values = gg_color_hue(n)) + guides(fill = guide_legend(reverse = T))
}



# monocle3 ----------------------------------------------------------------

make_monocle3 <- function(seurat_object, assay_type = "Integrated") {
  DefaultAssay(seurat_object) <- assay_type
  umi_matrix <- seurat_object@assays$RNA@data
  sample_info <- data.frame(seurat_object@meta.data,
                            stringsAsFactors = F)

  gene_annotation <- data.frame(gene_short_name = rownames(seurat_object))

  rownames(gene_annotation) <- rownames(umi_matrix)
  cds <- monocle3::new_cell_data_set(as(umi_matrix, "sparseMatrix"),
                           cell_metadata = sample_info,
                           gene_metadata = gene_annotation)
}

 do_monocle <- function(cds = mono) {
   cds <- monocle3::preprocess_cds(cds, num_dim = 30)
   cds <- monocle3::reduce_dimension(cds) #UMAP reduce dimension (defalt)
   cds = monocle3::cluster_cells(cds, k = 5, reduction_method = "UMAP")
   monocle3::plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster",
              group_label_size = 6, cell_size = 1.5) #plot by cluster
   cds <- monocle3::align_cds(cds)
   cds <- monocle3::learn_graph(cds)
   monocle3::plot_cells(cds, color_cells_by = "seurat_clusters", group_label_size = 5,  label_leaves = T)
   return(cds)
 }


 mop <- function(color_cells_by = "cluster", cds = mono, ...) {
   monocle3::plot_cells(cds = cds, color_cells_by = color_cells_by, group_label_size = 5, label_leaves = T, ...)
 }

 do_diff_mono <- function(group_cells_by = "cluster", cds = mono, reference_cells =NULL,...) {
   monocle3::top_markers(cds = cds, group_cells_by= group_cells_by, reference_cells= reference_cells, cores=8, reduction_method = "UMAP", ...)
 }



# cor_analysis ------------------------------------------------------------

get_df <- function(object) {
  df <- object@assays$RNA@data %>% as.matrix() %>% t() %>% as.data.frame()
  #df <- df %>% rownames_to_column(var = "id")
}

do_cor <- function(expr_df, gene, group_label = "subset", method = "pearson") {
  nu <- str_which(colnames(expr_df), paste0("^", gene, "$"))
  all_res <- cor(x = expr_df[gene], y = expr_df[-nu], method = method)
  all_res_df <- tibble(gene = colnames(all_res), cor = as.numeric(all_res[1,])) %>% arrange(-cor)
  all_res_df$batch <- group_label
  return(all_res_df)
  # all_res_df %>% head(50) %>%
  #   ggplot(aes(fct_reorder(gene,cor), cor, fill = gene)) + geom_bar(stat= "identity") + coord_flip() + guides(fill = F)
  # ggsave(filename = paste0(group_label, "_barplot.jpg"), device = "jpeg")
  # return(all_res_df)
}
do_cor2 <- function(expr_df, gene, group_label = "subset", method = "pearson") {
  nu <- str_which(colnames(expr_df), paste0("^", gene, "$"))
  all_res_p_value <- outer(expr_df[gene],expr_df[-nu], function(X, Y){
    mapply(function(...) cor.test(..., na.action = "na.exclude")$p.value,
           X, Y)
  })
  all_res_estimate <- outer(expr_df[gene],expr_df[-nu], function(X, Y){
    mapply(function(...) cor.test(..., na.action = "na.exclude")$estimate,
           X, Y)
  })
  all_res_df <- tibble(gene = colnames(all_res_estimate), cor = as.numeric(all_res_estimate[1,]), p_value = as.numeric(all_res_p_value[1,])) %>% arrange(-cor)
  all_res_df$batch <- group_label
  return(all_res_df)
  # all_res_df %>% head(50) %>%
  #   ggplot(aes(fct_reorder(gene,cor), cor, fill = gene)) + geom_bar(stat= "identity") + coord_flip() + guides(fill = F)
  # ggsave(filename = paste0(group_label, "_barplot.jpg"), device = "jpeg")
  # return(all_res_df)
}


do_cor_test <- function(expr_df) {
  outer(expr_df[, c(1)], expr_df[, c(2,4)], function(X, Y){
    mapply(function(...) cor.test(..., na.action = "na.exclude")$estimate,X, Y)
  })
}

# color function ----------------------------------------------------------

gg_color_hue <- function(n, l = 65, c = 100, randam =F) {
  hues = seq(15, 375, length = n + 1)
  color_vec <- hcl(h = hues, l = l, c = c)[1:n]

  if(randam){
   color_vec <- color_randam(color_vec = color_vec, n = n)
  }
  return(color_vec)
}

color_randam <- function(color_vec, n) {
  if(n <=2)return(color_vec)
  randam_vec <- vector(length = n)
  set.seed(100)
  randam_vec[seq(1,n,2)] <- color_vec[sample(seq(1,n,2), replace = F)]
  randam_vec[seq(2,n,2)] <- color_vec[sample(seq(2,n,2), replace = F)]
  color_vec <- randam_vec
}





# batch_test --------------------------------------------------------------

batch_cor_plot <- function(x, res_cor) {
  res_cor %>% reshape2::melt() %>% filter(str_detect(Var1, x), !str_detect(Var2, "_")) %>%
    ggplot(aes(Var1, Var2, fill =value)) + geom_tile(color = "black") +scale_fill_gradient(low = "white", high = "red")+
    theme(axis.text.x = element_text(angle = 90))
  ggsave(paste0(x,"_cor2.jpg" ))
}



batch_feature <- function(object = data) {
  disease_list <- unique(object$disease)
  for(i in seq_along(disease_list)){
    data_sub <- sub_fil(object = object, disease == disease_list[i])
    up(data_sub)
    ggsave(paste0(disease_list[i],"_ump.jpg"))
    tile(endothelial_list, object = data_sub)
    ggsave(paste0(disease_list[i],"_endo.jpg"))
    tile(gene_list,object = data_sub)
    ggsave(paste0(disease_list[i],"_gene.jpg"))
  }
}


batch_mat <- function(average_df = av_df, object = data) {
  disease_list <- unique(object$disease)
  use_features <- rownames(average_df)
  for(i in seq_along(disease_list)){
    cat("executing____", disease_list[i],"________\n")
    data_sub <- sub_fil(object = object, disease == disease_list[i])
    df_sub <- AverageExpression(data_sub, assays = "RNA",features = use_features) %>% .[[1]]
    colnames(df_sub) <- paste0(colnames(df_sub),"_", disease_list[i])

    average_df <- cbind(average_df, df_sub)
  }
  return(average_df)
}



batch_cor_heatmap <- function(av_df_batch, method = "pearson") {
  use_order <- av_df_batch %>% colnames() %>% sort()
  res_cor <- cor(av_df_batch, method = method)
  res_cor  %>% reshape2::melt() %>% mutate(Var1 = fct_relevel(Var1, use_order),
                                                Var2 = fct_relevel(Var2, use_order)) %>%
  ggplot(aes(Var1, Var2, fill = value)) +
    geom_tile()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5)
                       #axis.ticks.y = element_blank(),axis.text.y = element_blank(),
                       #axis.ticks.x = element_blank(),axis.text.x = element_blank()
                       )+
    #scale_fill_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 0.9)
    scale_fill_gradientn(colours = c("red","yellow","white","lightblue","darkblue"),
                         values = c(1.0,0.9,0.5,0.4,0.25,0))
}


# clustering --------------------------------------------------------------


aut_clust <- function(object = object, original_label_name = "seurat_clusters", marker) {
  Idents(object) <- original_label_name
  use_df <- sig_val(object = object, marker = marker) %>% sig_val2()

  sig_df <- use_df %>% group_by(cluster) %>% top_n(1, mean) %>% dplyr::arrange(cluster) %>% dplyr::select(signature, cluster)

  object[[]] %>%
    left_join(sig_df, by = c("seurat_clusters" = "cluster")) %>%
    pull(signature) -> label
  return(label)
  # object_name <- as.character(substitute(object))
  # assign(x = object_name, value = object, envir = .GlobalEnv)
}

make_av_df2 <- function (data, gene_list, use_label)
{

  mt <- FetchData(object = data, unique(unlist(gene_list)), slot = "data")
  mt$cluster <- data[[]] %>% pull(use_label)
  df <- mt %>% group_by(cluster) %>% summarise_all(mean)
  df <- column_to_rownames(df, var = "cluster")
}

clu_data <- function(data = data, gene_list = gene_list, use_label = "seurat_clusters", k = 6 ) {
  df <- make_av_df2(data = data, gene_list = gene_list, use_label = use_label)
  res <- dist(df)
  res <- hclust(res)
  clust <- cutree(res, k = k)
  tb <- data.frame(clu = as.character(clust), use_label = names(clust))
  res <- data[[use_label]] %>% left_join(tb, by = use_label) %>% pull(clu) %>% fct_relevel(as.character(1:k))
  return(res)
}





# network_analysis_singlecell_singal_R ------------------------------------

make_clust_list <- function(seurat_data, meta_name) {
  seurat_data$id <- colnames(seurat_data)
  cluster <-
    seurat_data[[]] %>% select(id, cluster_name = meta_name) %>%
    mutate(
      cluster_no = fct_drop(cluster_name) %>% as.numeric,
      id = paste(id, cluster_name, sep = ".")
    )

  use_name <-
    cluster %>% distinct(cluster_name, cluster_no) %>% arrange(cluster_no) %>% pull(cluster_name) %>% as.character()
  use_cluster <- cluster$cluster_no
  names(use_cluster) <- cluster$id

  cluster_list <- list(use_name = use_name,
                       use_cluster = use_cluster)

  return(cluster_list)
}


inter_signal <- function(seurat_data, cluster_list) {
  genes <- rownames(seurat_data@assays$RNA@counts )
  signal <-
    SingleCellSignalR::cell_signaling(
      data = as.matrix(seurat_data@assays$RNA@counts),
      genes = genes,
      cluster = cluster_list$use_cluster,
      c.names = cluster_list$use_name,
      write = FALSE
    )
  return(signal)
}


inter_net <- function(seurat_data, signal, cluster_list) {
  inter.net <-
    SingleCellSignalR::inter_network(
      data = as.matrix(seurat_data@assays$RNA@data),
      signal = signal,
      genes = rownames(seurat_data),
      cluster = cluster_list$use_cluster,
      c.names = cluster_list$use_name,
      write = FALSE
    )
  return(inter.net)

}

do_signal <- function(seurat_data, meta_name) {
  use_clust_list <-
    make_clust_list(seurat_data = seurat_data, meta_name = meta_name)
  signal <-
    inter_signal(seurat_data = seurat_data, cluster_list = use_clust_list)
}




# add_disease_label -------------------------------------------------------


com_label <- function(data = data, first_label = "first", second_label = "second", label_name = "com_label", add_num = T) {
  data[[label_name]] <- interaction(data[[first_label]][,], data[[second_label]][,])
  if(add_num){
    n_tab <- table(sub_hep[[label_name]][,])
    n_table <- as.vector(n_tab)
    names(n_table) <- names(n_tab)
    data[[label_name]][,] <- paste0(data[[label_name]][,], "(n=", n_table[data[[label_name]][,]], ")")
  }
  return(data)
}



# change_label ------------------------------------------------------------



label_change <- function(object, data, reference_name) {

  stopifnot(class(object)=="Seurat")
  stopifnot(class(data)=="Seurat")
  stopifnot(class(reference_name)== "character")

  reference_list <- unique(data$source_ref)
  stopifnot(reference_name %in% reference_list)

  use_tail <- data[[]] %>%
    filter(source_ref== reference_name) %>%
    pull(id) %>%
    str_extract("_\\d{1,2}$") %>%
    unique()

  object$id <- colnames(object) %>% paste0(use_tail)

  use_df <- data[[]] %>% select(id, use_label)

  compositional_list <- object[[]] %>% left_join(use_df, key = "id") %>%
    filter(!is.na(use_label)) %>%
    select(seurat_clusters, use_label) %>%
    group_by(seurat_clusters, use_label) %>%
    count() %>%
    group_by(seurat_clusters) %>%
    top_n(1, n) %>%
    select(-n)


  compositional_list

  object$use_label <- object[[]] %>%
    left_join(compositional_list, key = "seurat_cluters") %>%
    pull(use_label)

  return(object)
}



# gene_vs_analysis --------------------------------------------------------


gene_binom <- function(data, gene, cutoff_value) {

  stopifnot(class(data) == "Seurat")
  stopifnot(class(gene) == "character")
  stopifnot(gene %in% rownames(data))

  use_df <- FetchData(data, vars = gene)
  use_df <- use_df %>% mutate_all(
      ~if_else(.x > cutoff_value,"positive", "negative") %>% as.factor %>% fct_relevel("positive", "negative")
    )

  use_df <- use_df %>% dplyr::select(contains(gene))
  colnames(use_df) <- paste0(colnames(use_df), "_pos_nega")

  data@meta.data <- data@meta.data %>% cbind(use_df)
  return(data)
}

gene_binom_vs <- function(data, gene, cutoff_value = 0) {
  data <- gene_binom(data = data, gene = gene,cutoff_value = cutoff_value)
  Idents(data) <- paste0(gene,"_pos_nega")

  result <- FindAllMarkers(data, only.pos = T, min.pct = 0.2, min.diff.pct = 0.15)
  return(result)
}



sep_marker <- function(gene_name, cell_type, source_name) {
  Idents(data) <- paste0(gene_name,"_pos_nega")


  sub <- sub_fil(data, use_label == cell_type, source_dis == source_name)
  result <- FindAllMarkers(sub, only.pos = T, min.pct = 0.2, min.diff.pct = 0.15)
  result$cell_gene <- paste0(cell_type,"_", gene_name)
  return(result)
}


gene_vs <- function(data, gene, cutoff_value = 0) {

  stopifnot(class(data) == "Seurat")
  stopifnot(class(gene) == "character")
  stopifnot(gene %in% rownames(data))

  data <- gene_binom(data = data, gene = gene, cutoff_value = cutoff_value)

  use_df <- FetchData(data, c(gene, "source_dis", "use_label"))

  rank_df <- use_df %>%
    group_by(source_dis, use_label) %>%
    summarise_all(~(sum(.x > cutoff_value))/length(.x)) %>%
    pivot_longer(cols = -c(1,2), names_to = "gene", values_to = "value") %>%
    group_by(source_dis,gene) %>%
    mutate(rank = row_number(-value)) %>%
    filter(rank ==1)

  marker_df <- rank_df %>%
    mutate(def = pmap(list(gene, use_label, source_dis), ~try(sep_marker(..1,..2,..3))))

  marker_df <- marker_df %>% filter(class(def[[1]])!="try-error")

  return(marker_df)
}





a <- function() {
  use_df <- FetchData(data, c(gene, bb))
  enquo()
  rank_df <- use_df %>%
    group_by(bb) %>%
    summarise_all(~(sum(.x > cutoff_value))/length(.x)) %>%
    pivot_longer(cols = -c(1,2), names_to = "gene", values_to = "value") %>%
    group_by(source_dis,gene) %>%
    mutate(rank = row_number(-value)) %>%
    filter(rank ==1)
}



