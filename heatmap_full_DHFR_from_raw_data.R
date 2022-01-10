library(pheatmap)
library(stringr)
library(matrixStats)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(hash) #to store infos for the different tiles

# read in raw data as it comes from enrich. We use the file 'main_synonymous_scores.tsv' in the experiement folder of the tsv folder. 
# for the full seq of dhfr that is i.e. tile1_all_reps_minvarcount10_corr_sy/tsv/tile1_all_reps_exp/main_synonymous_scores.tsv
# removes double mutants
# does not do anything about vars that occur in two tiles since each tile is plotted separately
# tile borders have to specified inside the 'make_tile_data' function
# plots each tile from the specified start position to end position
# might have issues if positions are supplied for which there is no data. 
# Same with substitutions for which there is no data, which is why I take Ter out of the order of aa substitutions for tiles
# 3,4,5 since they appear to not have scores for early terminations.
# To make a heatmap encompassing all tiles, first the data needs to be compiled into a single dataframe. This is done in the script:
# prep_MAVE_data.R

#OBS: change paths to what is applicable for you. You can find all lines where paths are mentioned by searching for 'base_dir'
#based on heatmap_shorter.R . For the new full dhfr MAVE data that only has one condition but all 5 tiles (should be full DHFR)
#heatmap making functions from heatmap_shorter.R

base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/results/enrich_full_dhfr/'

corr_descr <- 'corrected with synonymous counts' #'corrected with counted reads' #'corrected with WT counts' # 'corrected with all reads' # 'corrected with counted reads'
corr <- '_corr_sy' #'_corr_complete' #'' #'_corr_full' #'_corr_complete'
#tile <- "5"

#dict storing necessary infos on the 5 different tiles
############################

make_tile_data <- function(){
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
  h[['2']][['WT_seq']] <- c('Arg', 'Met', 'Thr', 'Thr', 'Thr', 'Ser', 'Ser', 'Val', 'Glu', 'Gly', 'Lys', 'Gln', 'Asn', 'Leu',
                            'Val', 'Ile', 'Met', 'Gly', 'Lys', 'Lys', 'Thr', 'Trp', 'Phe', 'Ser', 'Ile', 'Pro', 'Glu', 'Lys', 
                            'Asn', 'Arg', 'Pro', 'Leu', 'Lys', 'Gly', 'Arg', 'Ile')
  
  h[["3"]] <- hash()
  #the mutagenized region is 75 - 121 but we only have read overlap (from fw rv pair) from aa 81 to 116
  h[['3']][['start_pos']] <- 81
  h[['3']][['tile_end']] <- 116
  h[['3']][['WT_seq']] <- c('Lys', 'Glu', 'Pro', 'Pro', 'Gln', 'Gly', 'Ala', 'His', 'Phe', 'Leu', 'Ser', 'Arg', 'Ser', 'Leu', 'Asp',
                            'Asp', 'Ala', 'Leu', 'Lys', 'Leu', 'Thr', 'Glu', 'Gln', 'Pro', 'Glu', 'Leu', 'Ala', 'Asn', 'Lys', 'Val',
                            'Asp', 'Met', 'Val', 'Trp', 'Ile', 'Val' )
  
  h[["4"]] <- hash()
  #the mutagenized region is 119 - 156 but we have read overlap (from fw rv pair) from aa 116 to 159 so I included counting those muts
  #there should be none if I understand correctly
  h[['4']][['start_pos']] <- 116
  h[['4']][['tile_end']] <- 159
  h[['4']][['WT_seq']] <- c('Val', 'Gly', 'Gly', 'Ser', 'Ser', 'Val', 'Tyr', 'Lys', 'Glu', 'Ala', 'Met', 'Asn', 'His', 'Pro', 'Gly',
                            'His', 'Leu', 'Lys', 'Leu', 'Phe', 'Val', 'Thr', 'Arg', 'Ile', 'Met', 'Gln', 'Asp', 'Phe', 'Glu', 'Ser',
                            'Asp', 'Thr', 'Phe', 'Phe', 'Pro', 'Glu', 'Ile', 'Asp', 'Leu', 'Glu', 'Lys', 'Tyr', 'Lys', 'Leu')
  
  h[["5"]] <- hash()
  #the mutagenized region is 153 - 187 but we have read overlap (from fw rv pair) from aa 146 til the stop codon (pos after 187)
  h[['5']][['start_pos']] <- 146
  h[['5']][['tile_end']] <- 187
  h[['5']][['WT_seq']] <- c('Asp', 'Thr', 'Phe', 'Phe', 'Pro', 'Glu', 'Ile', 'Asp', 'Leu', 'Glu', 'Lys', 'Tyr', 'Lys', 'Leu', 'Leu',
                            'Pro', 'Glu', 'Tyr', 'Pro', 'Gly', 'Val', 'Leu', 'Ser', 'Asp', 'Val', 'Gln', 'Glu', 'Glu', 'Lys', 'Gly',
                            'Ile', 'Lys', 'Tyr', 'Lys', 'Phe', 'Glu', 'Val', 'Tyr', 'Glu', 'Lys', 'Asn', 'Asp')
  
  return(h)
}
h <- make_tile_data()

# get values example
#h[["1"]]
#h[['1']][['start_pos']]
#h[[tile]][['start_pos']]

#################################

#for including different conditions see same function in script heatmaps_shorter.R
make_heatmap <- function(dat, start_pos, tile_end, pos_WT, aa_WT, order_aa, min_score=0,max_score=0, descr='', mid_point = 0) {
  
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
    
    #add WT markers. Lets try drawing black circles
    geom_tile(data = WT, aes(x = pos_WT, y=aa_WT, fill = dummy), color = 'grey50') + #scale_color_manual(values = c("white")) +
    geom_point(data = WT, aes(x = pos_WT, y=aa_WT), size = 3, color = 'black') +
    
    #this works but only displays x axis ticks every 10 positions
    #xlab('pos')+ ylab('substitution') +
    #which x-axis ticks to draw. It seems I need to use 1 - tile_end as limits, otherwise it just displays no x-axis scale.
    #the labels and breaks are me trying to get it to display only labels between tile start and tile end, not all the way from 1. 
    scale_x_discrete(name ="pos", breaks=seq(start_pos,tile_end,1), labels=as.character(seq(start_pos,tile_end,1)), 
                     limits=as.factor(seq(1,tile_end,1))) + ylab('substitution') +
    
    #limit plot to from start of the tile to the end so we cut off all the part from 1 to start_tile
    coord_cartesian(xlim = c(start_pos,tile_end)) +
    #ggtitle(paste0("MAVE scores tile ", tile, ", 0 mismatch allowed, min base qual 1, average read qual 1")) + 
    ggtitle(descr) +
    theme(axis.text=element_text(size=16), axis.title=element_text(size=18), legend.text=element_text(size=18), 
          plot.title = element_text(size=22))
  
  #add color scale later so we can determine what the mid point should be
  if(mid_point == 'wt'){
    #use this to use the same color scale limits for all heatmaps
    p <- p + scale_fill_gradient2(name = "MAVE Score", midpoint = wt_score, low = "#830823", high = "#063264", mid = 'white', 
                         limits = c(min_score,max_score))
  } else if(mid_point == 'sy'){
    p <- p + scale_fill_gradient2(name = "MAVE Score", midpoint = sy_score, low = "#830823", high = "#063264", mid = 'white', 
                         limits = c(min_score,max_score)) 
  } else {
    p <- p + scale_fill_gradient2(name = "MAVE Score", midpoint = mid_point, low = "#830823", high = "#063264", mid = 'white', 
                         limits = c(min_score,max_score)) 
  }
  
  #debugging
  #####################
  #png(filename = paste0(base_dir,'HZ_plots/heatmap_tile',tile,'.test.png'), width = 1400, height = 800)
  #print(p)
  #dev.off()
  
  return(p)
  
}

make_error_heatmap <- function(dat, start_pos, tile_end, pos_WT, aa_WT, order_aa, max_err=0, descr='') {

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
    
    #this works but only displays x axis ticks every 10 positions
    #xlab('pos')+ ylab('substitution') +
    #which x-axis ticks to draw. It seems I need to use 1 - tile_end as limits, otherwise it just displays no x-axis scale.
    #the labels and breaks are me trying to get it to display only labels between tile start and tile end, not all the way from 1. 
    scale_x_discrete(name ="pos", breaks=seq(start_pos,tile_end,1), labels=as.character(seq(start_pos,tile_end,1)), 
                     limits=as.factor(seq(1,tile_end,1))) + ylab('substitution') +
    #limit plot to from start of the tile to the end
    coord_cartesian(xlim = c(start_pos,tile_end)) +
    #ggtitle(paste0("MAVE score errors tile ", tile, ", 0 mismatch allowed, min base qual 1, average read qual 1")) + 
    ggtitle(descr) +
    theme(axis.text=element_text(size=16), axis.title=element_text(size=18), legend.text=element_text(size=18), 
          plot.title = element_text(size=22))
  
  return(p)
  
}

#go through all 5 tiles and make heatmaps
##########################
tile <- '1' #for debug
for (tile in c('1','2','3','4','5')) {
  print(c('tile', tile))
  
  #position of the first residue
  start_pos <- h[[tile]][['start_pos']]
  #last position
  tile_end <- h[[tile]][['tile_end']]
  pos_WT <- seq(start_pos,tile_end,1)
  
  #this is the WT aa seq, change for different tiles of dhfr
  aa_WT <-  h[[tile]][['WT_seq']]
  
  #order of amino acids for the heatmap
  #In the last tile there are apparently no early termination single vars so I took that out for tile 5. Same with tiles 3 and 4
  if(tile == '5' | tile == '3' | tile == '4'){
    order_aa <- c('His', 'Lys', 'Arg', 'Asp', 'Glu', 'Cys', 'Met', 'Asn', 'Gln', 'Ser', 'Thr', 'Ala', 
                  'Ile', 'Leu', 'Val', 'Phe', 'Trp', 'Tyr', 'Gly', 'Pro', 'pos_median')
  } else{
    order_aa <- c('His', 'Lys', 'Arg', 'Asp', 'Glu', 'Cys', 'Met', 'Asn', 'Gln', 'Ser', 'Thr', 'Ala', 
                  'Ile', 'Leu', 'Val', 'Phe', 'Trp', 'Tyr', 'Gly', 'Pro', 'Ter', 'pos_median')  
  }
  
  #minimum 10 read counts to accept a var as real
  ##########################
  dat <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps_minvarcount10',corr,'/tsv/tile',tile ,'_all_reps_exp/main_synonymous_scores.tsv'), 
                  sep = '\t', header = F,
                  stringsAsFactors = F, skip = 2,
                  col.names = c('var', 'SE', 'epsilon', 'score'))
  
  #find out max and min scores
  print('minimum 10 DNA var counts')
  
  mid_point <- 'wt'
  #run first without args to get the min and max scores
  make_heatmap(dat = dat, start_pos =  start_pos, tile_end = tile_end, pos_WT =  pos_WT, aa_WT = aa_WT, order_aa = order_aa)
  
  #now fill in min_score and max_score. I often put the same one for all tiles so they are on the same scale so that means
  #it needs to run once for each tile and then use the min of all min_scores and the max of all max_scores for making the 
  #plots to be saved as files. For the full DHFR MAVE I chose -3.5 and +4
  p1 <- make_heatmap(dat = dat, start_pos =  start_pos, tile_end = tile_end, pos_WT =  pos_WT, aa_WT = aa_WT, order_aa = order_aa
                     , min_score = -3.5, max_score = 4.0, mid_point = mid_point
                     , descr = paste0("MAVE scores tile ", tile, ", ", corr_descr, ", mid point: ", mid_point,
                                    "\n0 mismatch allowed, min base qual 1, average read qual 1, ",
                                    "only DNA vars with minimum 10 read counts"))
  
  p1
  
  png(filename = paste0(base_dir,'HZ_plots/heatmaps/heatmap_tile',tile,'.minvarcount10', corr,'.midpoint',mid_point,'.png'), 
      width = 1400, height = 800)
  print(p1)
  dev.off()
  
  #run first without args to get the max error
  make_error_heatmap(dat = dat, start_pos =  start_pos, tile_end = tile_end, pos_WT =  pos_WT, aa_WT = aa_WT, order_aa = order_aa)
  
  #I used 1.0 as max error for all tiles after looking at the different max for the different tiles  
  max_err = 1.0
  p2 <- make_error_heatmap(dat = dat, start_pos =  start_pos, tile_end = tile_end, pos_WT =  pos_WT, aa_WT = aa_WT, order_aa = order_aa
                     , max_err = max_err
                     , descr = paste0("MAVE score errors tile ", tile, ", ", corr_descr,
                                      "\n0 mismatch allowed, min base qual 1, average read qual 1, ",
                                      "only DNA vars with minimum 10 read counts"))
  p2
  png(filename = paste0(base_dir,'HZ_plots/heatmaps/heatmap_tile',tile,'.err.minvarcount10', corr,'.max_err',max_err,'.png'), 
      width = 1400, height = 800)
  print(p2)
  dev.off()
  
}
    







