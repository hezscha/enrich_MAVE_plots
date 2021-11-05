library(pheatmap)
library(stringr)
library(matrixStats)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(hash) #to store infos for the different tiles

#OBS: change paths to what is applicable for you. You can find all lines where paths are mentioned by searching for 'base_dir'
#based on heatmap_shorter.R . For the new full dhfr MAVE data that only has one condition but all 5 tiles (should be full DHFR)
#heatmap making functions from heatmap_shorter.R

base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/results/enrich_full_dhfr/'

tile <- "1"

#dict storing necessary infos on the 5 different tiles
############################
h <- hash() 
# set values
h[["1"]] <- hash()
h[['1']][['start_pos']] <- 2
h[['1']][['tile_end']] <- 39
h[['1']][['WT_seq']] <- c('Val','Gly','Ser','Leu','Asn','Cys','Ile','Val','Ala','Val','Ser','Gln','Asn','Met','Gly','Ile','Gly','Lys',
  'Asn','Gly','Asp','Leu','Pro','Trp','Pro','Pro','Leu','Arg','Asn','Glu','Phe','Arg','Tyr','Phe','Gln','Arg','Met','Thr')

h[["2"]] <- hash()
#the mutagenized region is 31 - 77 but we only have read overlap (from fw rv pair) from aa 37 to 72
h[['2']][['start_pos']] <- 37
h[['2']][['tile_end']] <- 72

# get values
h[["1"]]
h[['1']][['start_pos']]
h[[tile]][['start_pos']]

#################################

#position of the first residue
start_pos <- h[[tile]][['start_pos']]
#last position
tile_end <- h[[tile]][['tile_end']]
pos_WT <- seq(start_pos,tile_end,1)

#this is the WT aa seq, change for different tiles of dhfr
aa_WT <-  h[[tile]][['WT_seq']]

#order of amino acids for the heatmap
order_aa <- c('His', 'Lys', 'Arg', 'Asp', 'Glu', 'Cys', 'Met', 'Asn', 'Gln', 'Ser', 'Thr', 'Ala', 
              'Ile', 'Leu', 'Val', 'Phe', 'Trp', 'Tyr', 'Gly', 'Pro', 'Ter', 'pos_median')

dat <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps/tsv/tile',tile ,'_all_reps_exp/main_synonymous_scores.tsv'), 
                sep = '\t', header = F,
                stringsAsFactors = F, skip = 2,
                col.names = c('var', 'SE', 'epsilon', 'score'))


#for including different conditions see same function in script heatmaps_shorter.R
make_heatmap <- function(min_score=0,max_score=0) {
  
  #get the synonymous and WT scores and subset data to only scores for single vars 
  sy_score <- dat[1,'score']
  wt_score <- dat[2,'score']
  #now omit these two rows
  dat <- dat[-c(1,2),]
  #get only scores for single mutants
  d <- dat[!grepl(',',dat$var, fixed = T),c('var','score')]

  #colnames(d) <- c('var','score')
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
    ggtitle(paste0("MAVE scores tile ", tile, ", 0 mismatch allowed, min base qual 1, average read qual 1")) + 
    theme(axis.text=element_text(size=16), axis.title=element_text(size=18), legend.text=element_text(size=18), 
          plot.title = element_text(size=22))
  
  return(p)
  
}

make_error_heatmap <- function(max_err=0) {

  #get the synonymous and WT scores and subset data to only scores for single vars 
  sy_error <- dat[1,'SE']
  wt_error <- dat[2,'SE']
  #now omit these two rows
  dat <- dat[-c(1,2),]
  #get only errors for single mutants (this matrix should represent the same vars as make_heatmap but show the error instead of score)
  d <- dat[!grepl(',',dat$var, fixed = T),c('var','SE')]
  
  if (max_err == 0){
    max_err <- max(d$SE) + 0.1  
  }
  print(paste0('max_err is: ', max_err))
  
  #colnames(d) <- c('var','score')
  #make cols for position and substiution
  d$pos <- as.numeric(str_extract(d$var, "([0-9]+)"))
  d$sub <- str_match(d$var, "[0-9]+([A-Za-z]+)")[,2] 
  #drop empty rows
  d <- d[complete.cases(d),]
  #order by position
  d <- arrange(d, pos)
  
  #cast into matrix format
  mat <- dcast(d, pos ~ sub, value.var = "SE")
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
  dummy <- rep(wt_error,length(pos_WT))
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
    ggtitle(paste0("MAVE score errors tile ", tile, ", 0 mismatch allowed, min base qual 1, average read qual 1")) + 
    theme(axis.text=element_text(size=16), axis.title=element_text(size=18), legend.text=element_text(size=18), 
          plot.title = element_text(size=22))
  
  return(p)
  
}


png(filename = paste0(base_dir,'HZ_plots/heatmap_tile',tile,'.png'), width = 1400, height = 800)
print(make_heatmap(min_score = -2.7, max_score = 2.7))
dev.off()

png(filename = paste0(base_dir,'HZ_plots/heatmap_tile',tile,'.err_own_scale.png'), width = 1400, height = 800)
make_error_heatmap()
dev.off()

