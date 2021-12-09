library(pheatmap)
library(stringr)
library(matrixStats)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(ggplot2)
library(gridExtra)

#script for making the heatmaps for the pilot seq of tile 1.
#OBS: change paths to what is applicable for you. You can find all lines where paths are mentioned by searching for 'base_dir'
base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/results/enrich/'

#position of the first residue
start_pos <- 2
#last position
tile_end <- 39 
pos_WT <- seq(start_pos,tile_end,1)

#this is the WT aa seq, change for different tiles of dhfr
aa_WT <- c('Val','Gly','Ser','Leu','Asn','Cys','Ile','Val','Ala','Val','Ser','Gln','Asn','Met','Gly','Ile','Gly','Lys',
           'Asn','Gly','Asp','Leu','Pro','Trp','Pro','Pro','Leu','Arg','Asn','Glu','Phe','Arg','Tyr','Phe','Gln','Arg','Met','Thr')

#order of amino acids for the heatmap
order_aa <- c('His', 'Lys', 'Arg', 'Asp', 'Glu', 'Cys', 'Met', 'Asn', 'Gln', 'Ser', 'Thr', 'Ala', 
              'Ile', 'Leu', 'Val', 'Phe', 'Trp', 'Tyr', 'Gly', 'Pro', 'Ter', 'pos_median')

dat <- read.csv(paste0(base_dir, 'latest_re-run/ex123_all_wo_20_22_38/tsv/ex123_exp/main_synonymous_scores.tsv'), sep = '\t', header = F,
               stringsAsFactors = F, skip = 2)
colnames(dat) <- c('var', 'MTX_37_SE','MTX_37_eps','MTX_37_score','MTXCu_37_SE','MTXCu_37_eps','MTXCu_37_score',
                  'MTX_SE','MTX_eps','MTX_score','MTXCu_SE','MTXCu_eps','MTXCu_score')

make_heatmap <- function(condition, min_score=0,max_score=0) {
  
  if (condition == 'MTX') { 
    sy_score <- dat[1,10]
    wt_score <- dat[2,10]
    #now omit these two rows
    dat <- dat[-c(1,2),]
    d <- dat[!grepl(',',dat$var, fixed = T),c(1,10)]
  } else if (condition == 'MTXCu') {
    sy_score <- dat[1,13]
    wt_score <- dat[2,13]
    #now omit these two rows
    dat <- dat[-c(1,2),]
    d <- dat[!grepl(',',dat$var, fixed = T),c(1,13)]
  } else if  (condition == 'MTX_37') {
    sy_score <- dat[1,4]
    wt_score <- dat[2,4]
    #now omit these two rows
    dat <- dat[-c(1,2),]
    d <- dat[!grepl(',',dat$var, fixed = T),c(1,4)]
  } else if  (condition == 'MTXCu_37') {
    sy_score <- dat[1,7]
    wt_score <- dat[2,7]
    #now omit these two rows
    dat <- dat[-c(1,2),]
    d <- dat[!grepl(',',dat$var, fixed = T),c(1,7)]
  } else {
    print('Choose one of valid conditons: MTX, MTX_37, MTXCu, MTXCu_37')
    return()
  }
  
  colnames(d) <- c('var','score')
  #make cols for position and substiution
  d$pos <- as.numeric(str_extract(d$var, "([0-9]+)"))
  d$sub <- str_match(d$var, "[0-9]+([A-Za-z]+)")[,2] 
  #drop empty rows
  d <- d[complete.cases(d),]
  #order by position
  d <- arrange(d, pos)
  
  #cast into matrix format
  mat <- dcast(d, pos ~ sub, value.var = "score")
  rownames(mat) <- mat$pos
  #remove pos column
  mat <- mat[,-c(1)]
  mat$pos_median = rowMedians(as.matrix(mat), na.rm = T)
  
  #need to padd missing rows, i.e. positions where nothing has scored with NA
  all_rows <- seq(start_pos,tile_end)
  #missing rows:
  missing_rows <- all_rows[!all_rows %in% rownames(mat)]
  
  for (i in missing_rows) {
    mat[as.character(i),] <-  rep(NA,ncol(mat))
  }
  
  #arrange cols == amino acids according to properties and rows by row number
  #https://stackoverflow.com/questions/32593434/r-change-row-order
  mat <- mat[order(as.numeric(row.names(mat))),order_aa]
  
  #create dummy for marking the wild type in the plot
  dummy <- rep(wt_score,length(pos_WT))
  WT <- data.frame(pos_WT, aa_WT, dummy)
  
  if (min_score == 0)
  {min_score <- min(d$score, na.rm = T) - 0.1}
  if (max_score == 0)
  {max_score <- max(d$score, na.rm = T) + 0.1}
  
  print(paste0('min_score is: ', min_score))
  print(paste0('max_score is: ', max_score))
  
  p <- ggplot(melt(as.matrix(mat)), aes(x=Var1, y=Var2)) +
    geom_tile(aes(fill = value), color = "grey50") +
    #https://stackoverflow.com/questions/29476925/change-colors-of-raster-plot-in-ggplot2
    #scale_fill_gradientn(name = "MAVE Score", colors = c("#830823", "#c3393b", "#d6604d", "#f09b7a", "#fce2d3", 
    #                                                     "#f9f2ee", "white", "#e1edf3", "#c0ddec", "#a6cfe4", "#3683bb", 
    #                                                     "#1b5b9d", "#063264")) +
    
    #use this to find automatic limits to the color scale
    #scale_fill_gradient2(name = "MAVE Score", midpoint = 0, low = "#830823", high = "#063264", mid = 'white', 
    #                     limits = c(-max(abs(min(mat, na.rm = T)),abs(max(mat, na.rm = T))),max(abs(min(mat, na.rm = T)),abs(max(mat, na.rm = T))))) +
    #use this to use the same color scale limits for all heatmaps
    scale_fill_gradient2(name = "MAVE Score", midpoint = 0, low = "#830823", high = "#063264", mid = 'white', 
                         limits = c(min_score,max_score)) +
  

    #add WT markers. Lets try drawing black circles
    geom_tile(data = WT, aes(x = pos_WT, y=aa_WT, fill = dummy), color = 'grey50') + #scale_color_manual(values = c("white")) +
    geom_point(data = WT, aes(x = pos_WT, y=aa_WT), size = 3, color = 'black') +
    
    #which ticks to draw. Doesn't work properly if I start at 2 because it will just write 2 on the empty 
    #column it keeps printing in front
    scale_x_discrete(name ="pos", limits=as.factor(seq(1,tile_end,1))) + ylab('substitution') +
    #limit plot to from start of the tile to the end
    coord_cartesian(xlim = c(start_pos,tile_end)) +
    ggtitle(paste0("MAVE scores ", condition, ", 4 mismatch allowed, min base qual 20, average read qual 20")) + 
    theme(axis.text=element_text(size=16), axis.title=element_text(size=18), legend.text=element_text(size=18), 
          plot.title = element_text(size=22))
  
  return(p)
  
}

make_error_heatmap <- function(condition, max_err=0) {
  
  if (condition == 'MTX') { 
    sy_score <- dat[1,8]
    wt_score <- dat[2,8]
    #now omit these two rows
    dat <- dat[-c(1,2),]
    d <- dat[!grepl(',',dat$var, fixed = T),c(1,8)]
  } else if (condition == 'MTXCu') {
    sy_score <- dat[1,11]
    wt_score <- dat[2,11]
    #now omit these two rows
    dat <- dat[-c(1,2),]
    d <- dat[!grepl(',',dat$var, fixed = T),c(1,11)]
  } else if  (condition == 'MTX_37') {
    sy_score <- dat[1,2]
    wt_score <- dat[2,2]
    #now omit these two rows
    dat <- dat[-c(1,2),]
    d <- dat[!grepl(',',dat$var, fixed = T),c(1,2)]
  } else if  (condition == 'MTXCu_37') {
    sy_score <- dat[1,5]
    wt_score <- dat[2,5]
    #now omit these two rows
    dat <- dat[-c(1,2),]
    d <- dat[!grepl(',',dat$var, fixed = T),c(1,5)]
  } else {
    print('Choose one of valid conditons: MTX, MTX_37, MTXCu, MTXCu_37')
    return()
  }
  
  colnames(d) <- c('var','score')
  #make cols for position and substiution
  d$pos <- as.numeric(str_extract(d$var, "([0-9]+)"))
  d$sub <- str_match(d$var, "[0-9]+([A-Za-z]+)")[,2] 
  #drop empty rows
  d <- d[complete.cases(d),]
  #order by position
  d <- arrange(d, pos)
  
  #cast into matrix format
  mat <- dcast(d, pos ~ sub, value.var = "score")
  rownames(mat) <- mat$pos
  #remove pos column
  mat <- mat[,-c(1)]
  mat$pos_median = rowMedians(as.matrix(mat), na.rm = T)
  
  #need to padd missing rows, i.e. positions where nothing has scored with NA
  all_rows <- seq(start_pos,tile_end)
  #missing rows:
  missing_rows <- all_rows[!all_rows %in% rownames(mat)]
  
  for (i in missing_rows) {
    mat[as.character(i),] <-  rep(NA,ncol(mat))
  }
  
  #arrange cols == amino acids according to properties and rows by row number
  #https://stackoverflow.com/questions/32593434/r-change-row-order
  mat <- mat[order(as.numeric(row.names(mat))),order_aa]
  
  #create dummy for marking the wild type in the plot
  dummy <- rep(wt_score,length(pos_WT))
  WT <- data.frame(pos_WT, aa_WT, dummy)
  
  if (max_err == 0){
    max_err <- max(d$score) + 0.1  
  }
  print(paste0('max_err is: ', max_err))
  
  p <- ggplot(melt(as.matrix(mat)), aes(x=Var1, y=Var2)) +
    geom_tile(aes(fill = value), color = "grey50") +
    #https://stackoverflow.com/questions/29476925/change-colors-of-raster-plot-in-ggplot2
    #scale_fill_gradientn(name = "MAVE Score", colors = c("#830823", "#c3393b", "#d6604d", "#f09b7a", "#fce2d3", 
    #                                                     "#f9f2ee", "white", "#e1edf3", "#c0ddec", "#a6cfe4", "#3683bb", 
    #                                                     "#1b5b9d", "#063264")) +
    
    #use this to find automatic limits to the color scale
    #scale_fill_gradient2(name = "MAVE Score", midpoint = 0, low = "#830823", high = "#063264", mid = 'white', 
    #                     limits = c(-max(abs(min(mat, na.rm = T)),abs(max(mat, na.rm = T))),max(abs(min(mat, na.rm = T)),abs(max(mat, na.rm = T))))) +
    #use this to use the same color scale limits for all heatmaps
    scale_fill_gradient(name = "Standard Error", low = "white", high = "#830823", 
                         limits = c(0,max_err)) +
    
    
    #add WT markers. Lets try drawing black circles
    geom_tile(data = WT, aes(x = pos_WT, y=aa_WT, fill = dummy), color = 'grey50') + #scale_color_manual(values = c("white")) +
    geom_point(data = WT, aes(x = pos_WT, y=aa_WT), size = 3, color = 'black') +
    
    #which ticks to draw. Doesn't work properly if I start at 2 because it will just write 2 on the empty 
    #column it keeps printing in front
    scale_x_discrete(name ="pos", limits=as.factor(seq(1,tile_end,1))) + ylab('substitution') +
    #limit plot to from start of the tile to the end
    coord_cartesian(xlim = c(start_pos,tile_end)) +
    ggtitle(paste0("MAVE errors ", condition, ", 4 mismatch allowed, min base qual 20, average read qual 20")) + 
    theme(axis.text=element_text(size=16), axis.title=element_text(size=18), legend.text=element_text(size=18), 
          plot.title = element_text(size=22))
  
  return(p)
  
}

make_distance_heatmap <- function(cond1, cond2) {
  
  #d <- d[cond1] - d[cond2]
  wt_score <- dat[dat$var == '_wt', paste0(cond1,'_score')] - dat[dat$var == '_wt',paste0(cond2,'_score')]
  dat <- dat[-c(1,2),]
  d <- as.data.frame(dat$var[!grepl(',',dat$var, fixed = T)])
  d$diff <- dat[!grepl(',',dat$var, fixed = T), paste0(cond1,'_score')] - dat[!grepl(',',dat$var, fixed = T),paste0(cond2,'_score')]
  colnames(d) <- c('var','diff')


  #make cols for position and substiution
  d$pos <- as.numeric(str_extract(d$var, "([0-9]+)"))
  d$sub <- str_match(d$var, "[0-9]+([A-Za-z]+)")[,2] 
  #drop empty rows
  d <- d[complete.cases(d),]
  #order by position
  d <- arrange(d, pos)
  
  #cast into matrix format
  mat <- dcast(d, pos ~ sub, value.var = "diff")
  rownames(mat) <- mat$pos
  #remove pos column
  mat <- mat[,-c(1)]
  mat$pos_median = rowMedians(as.matrix(mat), na.rm = T)
  
  #need to padd missing rows, i.e. positions where nothing has scored with NA
  all_rows <- seq(start_pos,tile_end)
  #missing rows:
  missing_rows <- all_rows[!all_rows %in% rownames(mat)]
  
  for (i in missing_rows) {
    mat[as.character(i),] <-  rep(NA,ncol(mat))
  }
  
  #arrange cols == amino acids according to properties and rows by row number
  #https://stackoverflow.com/questions/32593434/r-change-row-order
  mat <- mat[order(as.numeric(row.names(mat))),order_aa]
  
  #create dummy for marking the wild type in the plot
  dummy <- rep(wt_score,length(pos_WT))
  WT <- data.frame(pos_WT, aa_WT, dummy)
  
  p <- ggplot(melt(as.matrix(mat)), aes(x=Var1, y=Var2)) +
    geom_tile(aes(fill = value), color = "grey50") +
    #https://stackoverflow.com/questions/29476925/change-colors-of-raster-plot-in-ggplot2
    #scale_fill_gradientn(name = "MAVE Score", colors = c("#830823", "#c3393b", "#d6604d", "#f09b7a", "#fce2d3", 
    #                                                     "#f9f2ee", "white", "#e1edf3", "#c0ddec", "#a6cfe4", "#3683bb", 
    #                                                     "#1b5b9d", "#063264")) +
    
    #use this to find automatic limits to the color scale
    #scale_fill_gradient2(name = "MAVE Score", midpoint = 0, low = "#830823", high = "#063264", mid = 'white', 
    #                     limits = c(-max(abs(min(mat, na.rm = T)),abs(max(mat, na.rm = T))),max(abs(min(mat, na.rm = T)),abs(max(mat, na.rm = T))))) +
    #use this to use the same color scale limits for all heatmaps
    scale_fill_gradient2(name = "Score difference", low="orange", mid="white",high="violet", midpoint = 0, 
                        limits = c(-2,2)) +
    
    
    #add WT markers. Lets try drawing black circles
    geom_tile(data = WT, aes(x = pos_WT, y=aa_WT, fill = dummy), color = 'grey50') + #scale_color_manual(values = c("white")) +
    geom_point(data = WT, aes(x = pos_WT, y=aa_WT), size = 3, color = 'black') +
    
    #which ticks to draw. Doesn't work properly if I start at 2 because it will just write 2 on the empty 
    #column it keeps printing in front
    scale_x_discrete(name ="pos", limits=as.factor(seq(1,tile_end,1))) + ylab('substitution') +
    #limit plot to from start of the tile to the end
    coord_cartesian(xlim = c(start_pos,tile_end)) +
    ggtitle(paste0("MAVE score difference: ", cond1, " - ", cond2)) + 
    #this size is good for the pdf of all plots
    #theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text=element_text(size=14))
    theme(axis.text=element_text(size=16), axis.title=element_text(size=18), legend.text=element_text(size=18), 
          plot.title = element_text(size=22))
  
  return(p)
  
}


#distance heatmaps
##################################

#all in one pdf
p1 <- make_distance_heatmap('MTX','MTX_37')
p2 <- make_distance_heatmap('MTX','MTXCu')
p3 <- make_distance_heatmap('MTX_37','MTXCu_37')
p4 <- make_distance_heatmap('MTXCu','MTXCu_37')

pdf(paste0(base_dir,'HZ_plots/distance_heatmaps.pdf'), width = 10, height = 16)
grid.arrange(p1,p2,p3,p4, nrow=4)
dev.off()

#make_distance_heatmap('MTX','MTX_37')

#make_error_heatmap(condition)
##################################

#all in one pdf
p1 <- make_error_heatmap('MTX')
p2 <- make_error_heatmap('MTXCu')
p3 <- make_error_heatmap('MTX_37')
p4 <- make_error_heatmap('MTXCu_37')

pdf(paste0(base_dir,'HZ_plots/all_flipped_heatmaps_stderr.pdf'), width = 10, height = 16)
grid.arrange(p1,p2,p3,p4, nrow=4)
dev.off()

for (condition in c('MTX', 'MTX_37', 'MTXCu', 'MTXCu_37')) {
  png(filename = paste0(base_dir,'HZ_plots/MAVE_',condition,'_flipped_heatmap_stderr.png'), width = 1400, height = 800)
  print(make_error_heatmap(condition = condition, max_err = 2.5))
  dev.off()
}

#'normal' heatmaps (showing the MAVE score) 
##############################

#all in one pdf
p1 <- make_heatmap('MTX', min_score = -2.5, max_score = 2.5)
p2 <- make_heatmap('MTXCu', min_score = -2.5, max_score = 2.5)
p3 <- make_heatmap('MTX37', min_score = -2.5, max_score = 2.5)
p4 <- make_heatmap('MTXCu37', min_score = -2.5, max_score = 2.5)

pdf(paste0(base_dir,'HZ_plots/all_flipped_heatmaps.pdf'), width = 10, height = 16)
grid.arrange(p1,p2,p3,p4, nrow=4)
dev.off()

for (condition in c('MTX', 'MTX_37', 'MTXCu', 'MTXCu_37')) {
  print(condition)
  png(filename = paste0(base_dir,'HZ_plots/MAVE_',condition,'_flipped_heatmap.png'), width = 1400, height = 800)
  print(make_heatmap(condition = condition, min_score = -2.7, max_score = 2.7))
  dev.off()
}



