pkgs = c('RColorBrewer','plyr','ggplot2','GenomicRanges','wordcloud','reshape2','grid', 'optparse','grDevices')
lib = installed.packages()
installed=pkgs %in% rownames(lib)
not_installed = which(!installed)
if(length(not_installed)>0){
  for(i in not_installed) install.packages(pkgs[i])
}

color_clone_composition=c( 'M' = rgb(0,162,205, maxColorValue = 255), 
                           'B' = rgb(243,130,153, maxColorValue = 255), 
                           'P' = rgb(223,207,0, maxColorValue = 255))

color_density_plot=c( 'SNV' = rgb(172,152,199, maxColorValue = 255), 
                      'InDel' = rgb(214,147,114, maxColorValue = 255), 
                      'Gain' = rgb(53,154,95, maxColorValue = 255), 
                      'Loss' = rgb(227,171,108, maxColorValue = 255))

color_stats=c( 'SNV_diploid' = rgb(172,152,199, maxColorValue = 255), 'SNV_aneuploid'= 'grey20',
               'InDel_diploid' = rgb(214,147,114, maxColorValue = 255) , 'InDel_aneuploid'= 'grey40',
               'Gain_with_mutations' = rgb(53,154,95, maxColorValue = 255), 'Gain_no_mutations' = 'grey60',
               'Loss_with_mutations' = rgb(227,171,108, maxColorValue = 255), 'Loss_no_mutations'= 'grey80')

for(i in pkgs) suppressPackageStartupMessages(library(i, character.only =T))

mapply_pb = function(FUN, X, Y,  ...){
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)   
  
  # wrapper around FUN
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- mapply(wrapper, X, Y, ...)
  close(pb)
  res
}

lapply_pb = function(X, FUN, ...){
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)   
  
  # wrapper around FUN
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- lapply(X, wrapper, ...)
  close(pb)
  res
}

theme_cloneR=function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.background = element_blank(),
          panel.border     = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}

chr_levels = paste0('chr',c(as.character(1:22),'X','Y'))



write.bed=function(...){
  write.table(..., row.names=F, col.names=F, quote=F, sep="\t")
}

