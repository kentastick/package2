

# data to list ------------------------------------------------------------

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


# do_seurat ---------------------------------------------------------------

do_seurat <- function(data) {
  if(!class(data) == "Seurat"){
    object <- CreateSeuratObject(counts = data, project = "object", min.cells = 3, min.features = 200)
  }else{
    object <- data
  }
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
  object <- JackStraw(object, num.replicate = 10)
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

  saveRDS(object, "normal.rds")

  tmap(object = object, c('ALB', 'KRT7', 'CD68', 'CD3D', 'CD79A', 'IGJ','CD34', 'GNLY', 'FGFBP2', 'CD14', "MZB1"))
  ggsave("plot10.jpg")


  return(object)
}



#  function --------------------------------------------------------

#create present time as name
make_time <- function() {
  Sys.time() %>% str_remove_all('[:punct:]|\\s')
}


#plot
up <- function(object = data, label= TRUE,...) {
  DimPlot(object = object, label = label, ...)
}

upd <- function(label, object = data) {
  up(object = object, group.by = "disease", pt.size = 2) + gghighlight(str_detect(disease, label ), label_key = T)
}

upr <- function(label, object = data) {
  up(object = object, group.by = "reference", pt.size = 2) + gghighlight(str_detect(reference, label ), label_key = T)
}


ts <- function(object = data, ...) {
  DimPlot(object = object, reduction = 'tsne', label = TRUE, ...)
}

tmap <- function(features, object = data, ...) {
  FeaturePlot(features = features,object = object, reduction = 'tsne',  cols = c('lightgray', 'red'), min.cutoff = 0, ...)
}


ump <- function(features, object = data,...) {
  FeaturePlot(features = features,object = object, cols = c('lightgray', 'red'), min.cutoff = 0,...)
}


vl <- function(...) {
  VlnPlot(object = data,...)
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

search_gene <- function() {
  data@assays$RNA@data@Dimnames[[1]] %>% View
}

pick_gene <- function(pattern) {
  data@assays$RNA@data@Dimnames[[1]] %>% str_subset(., pattern =pattern)
}

# return absolute path ----------------------------------------------------

abpath <- function(path = clipr::read_clip()) {
  normalizePath(path) %>% stringr::str_replace_all("\\\\", "/") %>% clipr::write_clip()
}


# signature value calcuration ---------------------------------------------

#calculate mean expression of each marker list
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

# geometric mean ----------------------------------------------------------

gm_mean1 = function(a){prod(a)^(1/length(a))}
gm_mean1(1:10)


gm_mean2 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


# signature plot ---------------------------------------------------------------------

#make df of each cluster's signature score: each_value/max_value, fraction_rate: expressed_cell(>0)/total_cells

sig_val2 <- function(score_mt) {
  score_mt %>% gather(-cluster, -cell_id, key = "signature", value = "score") %>%
    group_by(signature, cluster) %>%
    summarise(fraction_of_cells = sum(score>0)/n(), mean = mean(score)) %>%
    group_by(signature) %>%
    mutate(max = max(mean)) %>%
    mutate(score = mean/max)
}


sig_val3 <- function(score_mt) {

  score_mt %>% gather(-cluster, -cell_id, key = "signature", value = "score") %>%
    group_by(cluster, signature) %>%
    summarise(fraction_of_cells = sum(score>0)/n(), mean = mean(score)) %>%
    group_by(cluster) %>%
    mutate(max = max(mean)) %>%
    mutate(score = mean/max)
}



#internal function in make_subset
signature_plot_ <- function(mat_value, use.color = c("#0099FF", "#FAF5F5", "#E32020")) {
  mat_value %>% ggplot(aes(cluster, signature, colour =score, size = fraction_of_cells)) + geom_point() +
    scale_colour_gradientn(colours = c("red","yellow","white","lightblue","darkblue"),
                           values = c(1.0,0.7,0.6,0.4,0.3,0))
}


#calculate whithin each gene score normalized by gene across the cluster
signature_plot <- function(marker = "gene_list", object = data, use_func = "mean", filter = F, use.color = c("#0099FF", "#FAF5F5", "#E32020")) {
    df <- sig_val(marker = marker, use_func = use_func, object = object, filter =filter)
    df <- sig_val2(score_mt = df)
    df %>% ggplot(aes(cluster, signature, colour =score, size = fraction_of_cells)) + geom_point() +
    scale_colour_gradientn(colours = c("red","yellow","white","lightblue","darkblue"),
                           values = c(1.0,0.7,0.6,0.4,0.3,0))
}
# signature_tile <- function(marker = "gene_list", object = data, use_func = "mean", label_name = "seurat_clusters", filter = F, use.color = c("#0099FF", "#FAF5F5", "#E32020")) {
#     df <- sig_val(marker = marker, use_func = use_func, object = object, filter =filter, label_name = label_name)
#     df <- sig_val2(score_mt = df)
#     df %>% ggplot(aes(cluster, signature, fill = score)) + geom_tile(color = "white") +
#     scale_fill_gradientn(colours = c("red","yellow","white","lightblue","darkblue"),
#                            values = c(1.0,0.7,0.6,0.4,0.3,0))
# }

signature_tile <- function(marker = "gene_list", object = data, use_func = "mean",filter = F, use.color = c("#0099FF", "#FAF5F5", "#E32020")) {
    df <- sig_val(marker = marker, use_func = use_func, object = object, filter =filter)
    df <- sig_val2(score_mt = df)
    df %>% ggplot(aes(cluster, signature, fill = score)) + geom_tile(color = "black") +
    #scale_fill_gradient2(low = "blue",  mid = "white", high = "red", midpoint = 0.5)
     scale_fill_gradientn(colours = c("red","yellow","white","lightblue","darkblue"),
                            values = c(1.0,0.7,0.6,0.4,0.3,0))
}


#calculate whithin each cluster score normalized by cluster value
signature_plot_within <- function(marker = "gene_list", object = data, use_func = "mean", filter = F, use.color = c("#0099FF", "#FAF5F5", "#E32020")) {
    df <- sig_val(marker = marker, use_func = use_func, object = object, filter =filter)
    df <- sig_val3(score_mt = df)
    df %>% ggplot(aes(signature, cluster, colour =score, size = fraction_of_cells)) + geom_point() +
    scale_colour_gradientn(colours = c("red","yellow","white","lightblue","darkblue"),
                           values = c(1.0,0.7,0.6,0.4,0.3,0))
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

load_list <- function(data_list) {
  file_name <-  data_list %>% str_extract(".*(?=.rds)")
  for(i in seq_along(data_list)){
    assign(x =file_name[i], value = readRDS(data_list[i]), envir = globalenv())
  }
}


# save function -----------------------------------------------------------


sa_data <- function(x) {
  parse_arg <- substitute(x)
  saveRDS(object = data, file = paste0(deparse(parse_arg), ".rds"))
}


sav <- function(x) {
  parse_arg <- substitute(x)
  saveRDS(object = x, file = paste0(deparse(parse_arg), ".rds"))
}


# marker_list functions  -----------------------------------------------------------


get_liver_marker <- function(n = 20, output =T) {
  gene_list_path = c(#"~/single_cell/single_cell_project/gene_list/",
    #"E:/single_cell_project/gene_list/",
    "~/single_cell/package2/test/gene_list/",
    "E:/single_cell_project/package2/test/gene_list/"
  )
  for (i in gene_list_path){
    liver_marker_list <- try(readRDS(file = paste0(i, "liver_marker_list.rds")), silent = T)
    if(class(liver_marker_list) != "try-error")break
  }
  liver_marker_list <- map(liver_marker_list, ~head(., n))
  if(output)assign(x = paste0("liver_marker_list_",n), value  = liver_marker_list, envir =globalenv())
}


save_list <- function(marker) {
  parse_name <- deparse(substitute(marker))
  gene_list_path = c(#"~/single_cell/single_cell_project/gene_list/",
                     #"E:/single_cell_project/gene_list/",
                     "~/single_cell/package2/test/gene_list/",
                     "E:/single_cell_project/package2/test/gene_list/"
                     )
  for (i in gene_list_path){
    saved_list <- try(saveRDS(marker, paste0(i, parse_name , ".rds")), silent = T)
    if(class(get_list) != "try-error")break
  }

}


get_list <- function(marker, output = T) {
  gene_list_path = c(#"~/single_cell/single_cell_project/gene_list/",
                     #"E:/single_cell_project/gene_list/",
                     "~/single_cell/package2/test/gene_list/",
                     "E:/single_cell_project/package2/test/gene_list/")
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
                     "E:/single_cell_project/package2/test/gene_list/")
  for (i in gene_list_path){
    res <- try(list.files(path = i))
    if(length(res) != 0) break
  }
  res
}


#remove duplicated gene in list

remove_list_dup <- function(gene_list) {
  gene_list <- gene_list %>% unlist() #%>% split(., ngene_listmes(.))
  gene_list <- gene_list[!duplicated(gene_list)]
  names(gene_list) <- names(gene_list) %>% str_remove("\\d{1,2}$")
  gene_list %>% split(., names(.)) %>% map(as.vector)
}


# version up package2 -----------------------------------------------------

restart <- function(remotes, install_github) {
  remotes::install_github("kentastick/package2")
  detach("package:package2", unload=TRUE)
  library("package2", lib.loc="~/R/win-library/3.6")
}




# differential gene expression analysis within same sample ----------------


#find marker wrapper function
diff_test <- function(object = data, ...) {
  FindAllMarkers(object = object, min.pct = 0.25, logfc.threshold = log(2), only.pos = T, ...)
}


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

#add entrez id list
diff_marker_convert <- function(marker_df) {
  parse_name <- deparse(substitute(marker_df))
  marker_df %>% group_by(cluster) %>% nest() %>%
    mutate(gene_symbol = map(data, ~filter(., p_val_adj<0.05) %>% pull(gene))) %>%
    mutate(gene_entrez_symbol = map(gene_symbol, ~convert_gene(.))) %>%
    mutate(gene_entrez = map(gene_entrez_symbol, ~pull(., ENTREZID))) %>%
    mutate(gene_list_entrez = map2(gene_entrez_symbol, data, ~make_gene_list(.x, .y) )) %>%
    mutate(gene_list_symbol = map2(data, gene_symbol, ~make_gene_list_2(.x, .y))) -> marker_df
  saveRDS(marker_df, paste0(parse_name, "_list.rds"))
  return(marker_df)
}

#sub lutins of diff_marker_convert

make_gene_list <- function(arg1, arg2) {
  df <-  arg1 %>% left_join(arg2,by = c("SYMBOL" = "gene") )
  val <- df$avg_logFC
  entrez_name <- df$ENTREZID
  names(val) <- entrez_name
  val <- sort(val, decreasing = T)
  return(val)
}

make_gene_list_2 <- function(arg1, arg2) {
  val <-  arg1 %>% filter(gene %in% arg2) %>% pull(avg_logFC)
  names(val) <- arg2
  val <- sort(val, decreasing = T)
  return(val)
}

#convert symbol to entrezid
convert_gene <- function(x) {
  library(org.Hs.eg.db)
  res <- try(clusterProfiler::bitr(x, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db"))

  return(res)
}

#
# marker_list_map <- function(df) {
#   df %>% mutate(gene_symbol = map(data, ~filter(., p_val_adj <  0.05) %>%
#                                     pull(gene))) %>%
#     mutate(gene_entrez_symbol = map(gene_symbol,  ~convert_gene(.))) %>%
#     mutate(gene_entrez = map(gene_entrez_symbol, ~try(pull(., ENTREZID)))) %>%
#     mutate(gene_list_entrez = map2(gene_entrez_symbol,data, ~try(make_gene_list(.x, .y)))) %>%
#     mutate(gene_list_symbol = map2(data, gene_symbol, ~try(make_gene_list_2(.x, .y))))
# }



# cell_origin bar plot ----------------------------------------------------

bar_origin <- function(bar_x, bar_y,object= data, position = "fill", randam = T) {
  bar_x <- enquo(bar_x)
  bar_y <- enquo(bar_y)
  p <- object@meta.data %>% dplyr::select(!!bar_x, !!bar_y) %>%
    ggplot(aes(!!bar_x, fill = !!bar_y)) + geom_bar(position = position)
  if(randam){
    n_color <- object@meta.data %>% select(!!bar_y) %>% n_distinct()
   p <- p +scale_fill_manual(values = gg_color_hue(n = n_color, randam = T))
  }
  print(p)
}

# meta data modifying function--------------------------------------------------------

#change current cell label
id_ch <- function(use_meta, object = data) {
  Idents(object = object) <- use_meta
  data <<- object
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
                                    CH = "CH", PBC = c("PBC_1", "PBC_2"),
                                    )) %>%
    pull(condition)
  data$id <- rownames(data@meta.data)

  return(data)


}


#add summarised marker_list average expression value

add_sig_val <- function(object = data, marker_list, use_func = "mean") {
  df_list <- vector("list", length(marker_list))
  for(i in seq_along(marker_list)){
    df_list[[i]] <- sig_val(object = object, marker = marker_list[i], use_func = use_func) %>%
      add_m(add = switch(use_func, "mean" = "_m", "gm_mean" = "_gm"))
  }

  df_com <- df_list %>% purrr::reduce(cbind)

  df_com <- df_com %>% keep(is.numeric)

  object <- add_meta(df = df_com, object = object)
  data <<- object
}


#add strong/weak value of one gene by specific value
add_meta_binval <- function(gene, object = data) {
  df <- FetchData(object = object, var = gene) %>% rownames_to_column(var = "id")
  gene_m <- df %>% pull(gene) %>%  mean()
  gene_sd <- df %>% pull(gene) %>%  sd()
  cut_off <- gene_m + 2*gene_sd
  colnames(df)[2] <- "temp"
  use_column <- df %>% dplyr::select(id, temp) %>% mutate(temp= if_else(temp> cut_off, "strong", "ordinary")) %>% dplyr::select(temp)
  colnames(use_column) <- paste0(gene,"_bin")
  data <- add_meta(df = use_column, object = data)
  data <<- data
}


#add meta data of data.frame to seurat object
add_meta <- function(df, object = data) {
  object@meta.data <- object@meta.data %>% rownames_to_column(var = "temp") %>%
    bind_cols(df) %>% column_to_rownames(var = "temp")
  object <<- object
}

add_m <- function(df_list, add = "_m") {
  colnames(df_list) <- paste0(colnames(df_list), add)
  return(df_list)
}




# subset filter---------------------------------------------------------

sub_fil <- function(object = data, ...) {
  use_id <- object@meta.data %>% filter(...) %>% pull(id)
  res <- subset(x = object, cells = use_id)
  return(res)
}

pick_id <- function(..., object =data) {
  object@meta.data %>% filter(...) %>% pull(id)
}


# gene annotation analysis ------------------------------------------------


convertMouseGene <- function(x){
  mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  genesV2 = biomaRt::getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
}


#reactomrPA
geneano_enricher <- function(gene_entrez, pvalueCutoff = 0.05, qvalueCutoff = 0.2) {
  res <- try(ReactomePA::enrichPathway(gene=gene_entrez, pvalueCutoff=pvalueCutoff, ,qvalueCutoff = qvalueCutoff, readable=T))
  #if(class(res)== "try-error")return(NULL)
}

#group go
geneano_groupgo <- function(gene, ont_type = c("MF","BP", "CC")) {
  library(org.Hs.eg.db)
  res_go <- list()
  for(i in seq_along(ont_type)){
  res_go[i] <- try(clusterProfiler::groupGO(gene     = gene,
            OrgDb    = org.Hs.eg.db,
            ont      = ont_type[i],
            level    = 3,
            readable = TRUE))
  }
  names(res_go) <- ont_type
  return(res_go)
}

#enrich go
geneano_enrichgo <- function(gene) {
        res <-  try(enrichGO(gene          = gene,
                  universe      = names(geneList),
                  OrgDb         = org.Hs.eg.db,
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE))
}


geneano_msig <- function(gene_symbol) {
  #gmtfile <- system.file("extdata", "c7.all.v7.0.symbols.gmt", package="clusterProfiler")
  gmtfile2 <- system.file("extdata", "c5.all.v7.0.symbols.gmt", package="clusterProfiler")
  #c7 <- read.gmt(gmtfile)
  c5 <- clusterProfiler::read.gmt(gmtfile2)
  res <- try(clusterProfiler::enricher(gene = gene_entrez, TERM2GENE=c5))

}

geneano_msig_gsea <- function(gene_list) {
  #gmtfile <- system.file("extdata", "c7.all.v7.0.symbols.gmt", package="clusterProfiler")
  gmtfile2 <- system.file("extdata", "c5.all.v7.0.symbols.gmt", package="clusterProfiler")
  #c7 <- read.gmt(gmtfile)
  c5 <- clusterProfiler::read.gmt(gmtfile2)
  res <- try(clusterProfiler::GSEA(gene = gene_list, TERM2GENE=c5))
}


do_geneano <- function(marker_df, res_name = "res_enricher",gene_type = gene_entrez, use_func = geneano_enricher) {

  gene_type <- enquo(gene_type)
  res_name <- enquo(res_name)
  marker_df %>% mutate(!!res_name := map(!!gene_type, ~use_func(gene = .)))
}


 add_enrich_anotation <- function(marker_list) {
   marker_list <- marker_list %>% do_geneano(use_func = geneano_enricher, res_name = "res_enricher", gene_type = gene_entrez)
   marker_list <- marker_list %>% mutate(bar_plot = map2(res_enricher, cluster, ~bar(.x, .y)))
   marker_list <- marker_list %>% mutate(cnet_plot = map2(res_enricher, cluster, ~cnet(.x, .y)))
 }





# filter  -----------------------------------------------------------------


add_all <- function(data) {
  data <- add_info(data = data)
  data <- add_sig_val(object = data)
}

fil_cell <- function(cell_type, remove_cluster = NULL, remove_disease = NULL) {

  data[[]] %>% dplyr::select(id,ends_with("_m")) %>% gather(-id,key = "type",value = "value" ) %>%
    group_by(id) %>%  mutate(rank = row_number(-value)) %>% arrange(id,rank) %>%
    filter(type == paste0(cell_type, "_m"), rank == 1|value>0) %>% pull(id)-> use_id
  data <- sub_fil(object = data, id %in%use_id, !seurat_clusters %in% remove_cluster, !disease %in% remove_disease)
}


fil_cell2 <- function(cell_type) {

  data[[]] %>% dplyr::select(id,ends_with("_m")) %>% gather(-id,key = "type",value = "value" ) %>%
    group_by(id) %>%  mutate(rank = row_number(-value)) %>% arrange(id,rank) %>%
    filter(type %in% paste0(cell_type, "_m"), rank == 1, value>0) %>% pull(id)-> use_id

  data <- SubsetData(data, cells = use_id)
}




# barplot -----------------------------------------------------------------



bar <- function(arg1, arg2, arg3 = "cluster", pathname = "pathway_plot") {

  if(!dir.exists(pathname)) dir.create(pathname)
  res <- try(barplot(arg1, title =as.character(arg2),showCategory = 20, supressResult = T))
  if(class(res)=="try-error"){
    return(NULL)
  }else res
    ggsave(filename = paste0(pathname, "/", arg3,"_", as.character(arg2), "_enrichplot.jpg"),
         device = "jpg", width = 30, height = 30, units = "cm")
  return(res)
}

cnet <- function(arg1, arg2, arg3 = "cluster", pathname = "pathway_plot") {
  if(!dir.exists(pathname)) dir.create(pathname)
  res <- try(clusterProfiler::cnetplot(x = arg1, title = as.character(arg2),showCategory = 20, supressResult = T))
   if(class(res)=="try-error"){
     return(NULL)
   }else res
  ggsave(filename = paste0(pathname, "/", arg3,"_", as.character(arg2), "_cnetplot.jpg"),
         device = "jpg", width = 30, height = 30, units = "cm")
  return(res)
}





# zip ---------------------------------------------------------------------

#function converting df to list
zip <- function(df) {
  list <- df$data_filtered
  names(list) <- df$batch
  return(list)
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


# get marker table --------------------------------------------------------

#result of findallmarker

get_marker_table <- function(marker_list, ...) {
   #filter_var <- enquo(filter_var)
  marker_list %>% filter(...) %>% dplyr::select(cluster, gene) %>%
    group_by(cluster) %>% nest %>%
    mutate(data = unlist(map(data, ~pull(.x, gene) %>% paste0(., collapse = ", ")))) %>%
    unnest() %>%
    write_clip()
}

get_diff_test_marker <- function(diff_test_res, ...) {
  diff_test_res %>% filter(...) %>% dplyr::select(cluster, batch, gene) %>%
    group_by(cluster, batch) %>% nest() %>%
    mutate(data = unlist(map(data, ~pull(.x, gene) %>% paste0(., collapse = ", ")))) %>%
    write_clip()
    #write_csv(., path = "diff_test_table.csv")
}


# tile_plot ---------------------------------------------------------------

tile <- function(gene, object = data, order = F, cluster_label = "seurat_clusters", plot_wrap = F, fil_val= NULL, color_label = T, ...) {
  #gene <- get_list(gene)
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
  use_id <- pick_id(object = object, ...)
  use_df <- object@assays$RNA@data[feature, use_id]
  cluster_label <- object@meta.data[,cluster_label]
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
    a<<- use_df

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
      #label_color_use <- label_color[as.numeric(as.factor(label_df$label))]
      use_df <- use_df %>% mutate(label_color = label_color[as.numeric(as.factor(label))])
      p <- use_df %>% ggplot(aes(cluster, gene, size = pct, fill = label, colour = score)) + geom_point() +
        scale_colour_gradientn(colours = c("red","yellow","white","lightblue","darkblue"),
                               values = c(1.0,0.7,0.6,0.4,0.3,0))
      p <- p + theme(axis.text.y = element_text(colour = label_color_use)) + scale_fill_manual(values = label_color)
      return(p)
    }


  }else return(p)

}



tile_legend <- function(df = for_tile_legend_df) {
  n <- length(unique(df$label))
  df %>% mutate(gene_label = fct_relevel(label, unique(df$label)),color = as.numeric(as.factor(label))) %>%
    ggplot(aes(gene_label, color, fill = gene_label)) +
    geom_bar(stat = "identity")+ scale_fill_manual(values = gg_color_hue(n)) + guides(fill = guide_legend(reverse = T))
}

tile2 <- function(gene, object = data, order =TRUE, ...) {
  if(class(gene) == "list" ){
    label_df <- enframe(gene, name = "label",value = "gene") %>% unnest
    feature <- unlist(gene)
  } else{
    feature <- gene
  }

  use_id <- pick_id(object = object, ...)
  use_df <-FetchData(object = object, feature, cells = use_id)
  cluster_label <- object@meta.data$seurat_clusters
  if(length(feature) ==1){
    use_df <- use_df %>% as.tibble()
  }else{
    use_df <- t(as.matrix(use_df)) %>% as.tibble()
  }
   use_df<- use_df %>% add_column(cluster = cluster_label)
  use_df <- use_df %>% tidyr::pivot_longer(-cluster, names_to = "gene", values_to = "logCPM") %>%
    group_by(cluster, gene) %>% summarise(avg_logCPM = mean(logCPM), pct = sum(logCPM>0)/n()) %>%
    group_by(gene) %>% mutate(score = avg_logCPM/max(avg_logCPM))
  if(class(gene) == "list" ){
    use_df <- use_df %>% left_join(label_df, by = c("gene"))
  }
  if(order){
    use_df <- use_df %>% mutate( cluster = fct_reorder(cluster, avg_logCPM))
  }

  p <- use_df %>% ggplot(aes(cluster, fct_relevel(gene,feature), size = pct, colour = score)) + geom_point() +
    scale_colour_gradientn(colours = c("red","yellow","white","lightblue","darkblue"),
                           values = c(1.0,0.7,0.6,0.4,0.3,0))

  if(class(gene) == "list" ){
   p + facet_grid(~label)
  }else return(p)



  }
#
# tile <- function(feature, object = data,...) {
#   use_id <- pick_id(object = object, ...)
#   use_df <- object@assays$RNA@data[feature, use_id]
#   cluster_label <- object@meta.data$seurat_clusters
#   if(length(feature) ==1){
#     use_df <- use_df %>% as.tibble()
#   }else{
#     use_df <- t(as.matrix(use_df)) %>% as.tibble()
#   }
#    use_df<- use_df %>% add_column(cluster = cluster_label)
#   use_df <- use_df %>% pivot_longer(-cluster, names_to = "gene", values_to = "logCPM") %>%
#     group_by(cluster, gene) %>% summarise(avg_logCPM = mean(logCPM))
#   a <<- use_df
#   use_df %>% ggplot(aes(cluster, gene,  fill = avg_logCPM)) + geom_tilet() +
#     scale_colour_gradientn(colours = c("red","yellow","white","lightblue","darkblue"),
#                            values = c(1.0,0.7,0.6,0.4,0.3,0))
# }



tile3 <- function(gene, object = data, order =TRUE, ...) {
  if(class(gene) == "list" ){
    label_df <- enframe(gene, name = "label",value = "gene") %>% unnest
    feature <- unlist(gene)
  } else{
    feature <- gene
  }
  use_df <-FetchData(object = object, feature, cells = use_id)

  use_df <- rownames_to_column("id") %>%
    mutate(id = row_number()) %>%
    pivot_longer(cols = -id, names_to = "gene", values_to = "log10CPM")

  if(class(gene) == "list" ){
    use_df <- use_df %>% left_join(label_df, by = c("gene"))
  }
  if(order){
    use_df <- use_df %>% mutate( cluster = fct_reorder(id, log10CPM))
  }

  p <- use_df %>% ggplot(aes(id, gene, fill = log10CPM)) + geom_tile() +
    scale_fill_gradientn(colours = c("red","yellow","white","lightblue","darkblue"),
                         values = c(1.0,0.7,0.6,0.4,0.3,0))
  if(color_label){
    n <- length(gene)
    label_color <- gg_color_hue(n)
    label_color_use <- label_color[as.numeric(plyr::mapvalues(label_df$label, from = unique(label_df$label), to = 1:n))]
    #label_color_use <- label_color[as.numeric(as.factor(label_df$label))]
    use_df <- use_df %>% mutate(label_color = label_color[as.numeric(as.factor(label))])
   p + theme(axis.text.y = element_text(colour = label_color_use)) + scale_fill_manual(values = label_color)
}

FetchData(data_sub, vars = unlist(segal_list)) %>% rownames_to_column("id") %>%
  mutate(id = row_number()) %>%
  pivot_longer(cols = -id, names_to = "gene", values_to = "log10CPM") %>%
  ggplot(aes(id, gene, fill = log10CPM)) + geom_tile() +
  scale_fill_gradientn(colours = c("red","yellow","white","lightblue","darkblue"),
                       values = c(1.0,0.7,0.6,0.4,0.3,0))
}

# monocle3 ----------------------------------------------------------------

make_monocle3 <- function(seurat_object) {
  umi_matrix <- seurat_object@assays$integrated@data
  sample_info <- data.frame(seurat_object@meta.data,
                            stringsAsFactors = F)

  gene_annotation <- data.frame(gene_short_name = rownames(seurat_object))

  rownames(gene_annotation) <- rownames(umi_matrix)
  cds <- monocle3::new_cell_data_set(as(umi_matrix, "sparseMatrix"),
                           cell_metadata = sample_info,
                           gene_metadata = gene_annotation)
}

 do_monocle <- function(cds = mono) {
   cds <- monocle3::preprocess_cds(cds, num_dim = 20)
   cds <- monocle3::reduce_dimension(cds) #UMAP reduce dimension (defalt)
   cds = monocle3::cluster_cells(cds, k = 7, reduction_method = "UMAP")
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

 get_earliest_principal_node <- function(cds, time_bin="130-170"){
   cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)

   closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
   closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
   root_pr_nodes <-
     igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]

   root_pr_nodes
 }



# Data2 <- monocle3::order_cells(Data1, reduction_method = "UMAP")
# monocle3::plot_cells(Data2, color_cells_by = "pseudotime" )




# filter_gene -------------------------------------------------------------


fil_gene <- function(gene, object) {
  if(class(gene) =="list"){
    gene <- gene %>% map(., ~.[. %in% rownames(object)])
    gene <- gene[map(gene, length)>0]
  }else{
    gene <- gene[gene %in% rownames(object)]
  }
  return(gene)
}


# cor_analysis ------------------------------------------------------------

get_df <- function(object) {
  df <- object@assays$RNA@data %>% as.matrix() %>% t() %>% as.data.frame()
  #df <- df %>% rownames_to_column(var = "id")
}

do_cor <- function(expr_df, gene, group_label = "subset", method = "pearson") {
  nu <- str_which(colnames(expr_df), paste0("^", gene, "$"))
  all_res <- cor(x = expr_d[gene], y = expr_df[-nu], method = method)
  all_res_df <- tibble(gene = colnames(all_res), cor = as.numeric(all_res[1,])) %>% arrange(-cor)
  all_res_df$batch <- group_label
  all_res_df %>% head(50) %>%
    ggplot(aes(fct_reorder(gene,cor), cor, fill = gene)) + geom_bar(stat= "identity") + coord_flip() + guides(fill = F)
  ggsave(filename = paste0(group_label, "_barplot.jpg"), device = "jpeg")
  return(all_res_df)
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
  randam_vec <- vector(length = n)
  set.seed(100)
  randam_vec[seq(1,n,2)] <- color_vec[sample(seq(1,n,2), replace = F)]
  randam_vec[seq(2,n,2)] <- color_vec[sample(seq(2,n,2), replace = F)]
  color_vec <- randam_vec
}



# assay type change -------------------------------------------------------


rna <- function() {
  DefaultAssay(object = data) <- "RNA"
  data <<- data
}
integ <- function() {
  DefaultAssay(object = data) <- "integrated"
  data <<- data
}



# write_count_table -------------------------------------------------------

wrt_c_tab <- function(data) {
  data$disease %>% table %>% as.data.frame() %>% t() %>% write_clip
}



