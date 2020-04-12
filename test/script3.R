ls() %>% map(., ~eval(parse(text = paste0("up(object =",.,")"))))
ls() %>% map(., ~try(eval(parse(text = paste0("ump(object =",.,", features = c(\"ACE2\", \"TMSRSS2\"))")))))
ls() %>% map(., ~try(eval(parse(text = paste0("ump(object =",.,", features = c(\"TMPRSS2\"))")))))
ls() %>% map(., ~try(eval(parse(text = paste0("VlnPlot(object =",.,", features = c(\"ACE2\"))")))))
ls() %>% map(., ~try(eval(parse(text = paste0("VlnPlot(object =",.,", features = c(\"ACE2\"))")))))
ls() %>% map(., ~try(add_sig_val(object = eval(as.name(.)), marker_list = "gene_list")))

setdiff(ls(), lsf.str()) %>% map(~try(add_sig_val(object_name = ., marker_list = "gene_list")))

names(gene_list)




do <- function() {
  up(object = eval(as.name("pbc_case1")))
}

as.character(substitute(eval(as.name("pbc_case1"))))


setdiff(ls(),lsf.str()) %>% map(~del_m(eval(as.name(.)))) -> object_list

object_list <- object_list %>% enframe()


object_list <- object_list %>% mutate(name = unlist(map(value,~.@meta.data$batch[[1]])))
object_list$name <- object_list$name %>% edit()
object_list <- object_list %>% mutate(value = map(value, ~add_sig_val(object = ., marker_list = "gene_list")))
object_list %>% mutate(value = map(value, ~add_sig_val(object = ., marker_list = "gene_list", add_signature_label = T)))

add_sig_val(pbc_case1, marker_list = "gene_list", label_name = "label")



data <- del_m(data)
