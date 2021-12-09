library(pheatmap)
library(stringr)
library(matrixStats)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(hash) #to store infos for the different tiles

#heatmap of MAVE scores full DHFR seq all tiles together. data from prep_MAVE_data.R
base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/results/enrich_full_dhfr/'

corr_descr <- 'corrected with synonymous counts' 
corr <- '_corr_sy'

#make heatmap function
#########################
make_heatmap <- function(d, start_pos, tile_end, pos_WT, aa_WT, min_score=0,max_score=0, descr='') {
  
  order_aa <- c('His', 'Lys', 'Arg', 'Asp', 'Glu', 'Cys', 'Met', 'Asn', 'Gln', 'Ser', 'Thr', 'Ala', 
                'Ile', 'Leu', 'Val', 'Phe', 'Trp', 'Tyr', 'Gly', 'Pro', 'Ter', 'pos_median') 
  
  #I'm doing this prep manually because there are special cases for the combined map due to vars present in two tiles
  
  # #get the synonymous and WT scores and subset data to only scores for single vars 
  #sy_score <- dat[1,'score']
  #wt_score <- dat[2,'score']
  # #now omit these two rows
  # dat <- dat[-c(1,2),]
  # #get only scores for single mutants
  # d <- dat[!grepl(',',dat$var, fixed = T),c('var','score')]
  # 
  # #colnames(d) <- c('var','score')
  # #make cols for position and substiution
  # d$pos <- as.numeric(str_extract(d$var, "([0-9]+)"))
  # d$sub <- str_match(d$var, "[0-9]+([A-Za-z]+)")[,2] 
  
  #just hard code this since we're using this function only for plotting all tiles together which is anyway a special case 
  wt_score <- 0
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
    theme(axis.text=element_text(size=12), axis.title=element_text(size=18), legend.text=element_text(size=18), 
          plot.title = element_text(size=22), axis.text.x = element_text(angle = 90))
  
  #debugging
  #####################
  #png(filename = paste0(base_dir,'HZ_plots/heatmap_tile',tile,'.test.png'), width = 1400, height = 800)
  #print(p)
  #dev.off()
  
  return(p)
  
}


#put all scores into one long df
##########################
#all_df <- data.frame()
tile <- '1'

all_df <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps_minvarcount10',corr,'/tsv/tile',tile ,'_all_reps_exp/main_synonymous_scores.tsv'), 
                          sep = '\t', header = F,
                          stringsAsFactors = F, skip = 4,
                          col.names = c('var', 'SE', 'epsilon', 'score'))


for (tile in c('2','3','4','5')) {
  print(c('tile', tile))

  #we will use wt correction, except that we replaced wt counts with sy counts since the sy score is much closer to the 
  #peak of the score distr (most vars behave similar to it). Therefore, all scores will be on the 'same' scale where 0 is 
  #performance like wt/sy. Therefore, we do not need to read in more wt scores (they are all 0) or more sy scores (they are all NA because of our trick with replacing the wt counts with sy counts)
  #therefore we skip 4 lines instead of 2 when reading in the other tiles
  dat <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps_minvarcount10',corr,'/tsv/tile',tile ,'_all_reps_exp/main_synonymous_scores.tsv'), 
                  sep = '\t', header = F,
                  stringsAsFactors = F, skip = 4,
                  col.names = c('var', 'SE', 'epsilon', 'score'))
  all_df <- rbind(all_df, dat)
  
}

#necessary prep for making the all tiles map
##################################

#drop multi variants
all_df <- all_df[!grepl(',',all_df$var, fixed = T),c('var','score')]

#add position and sub
all_df$pos <- as.numeric(str_extract(all_df$var, "([0-9]+)"))
all_df$sub <- str_match(all_df$var, "[0-9]+([A-Za-z]+)")[,2]

#to avoid problems later on, every sub can only occur once per position. This means I need to remove the overlapping vars, or pick one
#or do the average
#example: all_df[all_df$var == 'p.Asp153Ala',]

#find combos of sub and position with two score values (becaues they were seen in multiple tiles)
combos <- table(interaction(all_df$pos, all_df$sub))
doubles <- combos[combos == 2]

#calc the average between the two scores and replace them with that average score
for (item in names(doubles)) {
  
  #print('\n\n')
  #print(item)
  
  #split the name into the position and the subs so we can adress the relevant lines in the df
  a <- strsplit(item, '.', fixed=T)
  #test if we get the correct lines. Ok
  #print(all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], ])
  var_name <- unique(all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], c('var')])
  #calc a mean over all scores observed for this sub and position (should only be 2 since a position is in max two different tiles)
  mean_score <- mean(all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], c('score')])
  #print(mean_score)
  
  #remove the rows and add a new one with the mean
  all_df <- all_df[!(all_df$pos == as.numeric(a[[1]][[1]]) & all_df$sub == a[[1]][[2]]), ]
  new_row <- data.frame(var=var_name,score=mean_score,pos=as.numeric(a[[1]][[1]]),sub= a[[1]][[2]]) 
  #print('new_row is:')
  #print(new_row)
  all_df <- rbind(all_df, new_row)
  #print(all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], ])
}


#position of the first residue
start_pos <- 2
#last position
tile_end <- 187
pos_WT <- seq(start_pos,tile_end,1)

#this is the complete WT aa seq
          #tile 1: 2 .. 36
aa_WT <- c('Val','Gly','Ser','Leu','Asn','Cys','Ile','Val','Ala','Val','Ser','Gln','Asn','Met','Gly','Ile','Gly','Lys',
           'Asn','Gly','Asp','Leu','Pro','Trp','Pro','Pro','Leu','Arg','Asn','Glu','Phe','Arg','Tyr','Phe','Gln',
           #tile 2: 37 .. 72
           'Arg', 'Met', 'Thr', 'Thr', 'Thr', 'Ser', 'Ser', 'Val', 'Glu', 'Gly', 'Lys', 'Gln', 'Asn', 'Leu',
           'Val', 'Ile', 'Met', 'Gly', 'Lys', 'Lys', 'Thr', 'Trp', 'Phe', 'Ser', 'Ile', 'Pro', 'Glu', 'Lys', 
           'Asn', 'Arg', 'Pro', 'Leu', 'Lys', 'Gly', 'Arg', 'Ile',
           #gap from 73 .. 80
           'Asn', 'Leu', 'Val', 'Leu', 'Ser', 'Arg', 'Glu', 'Leu',
           #tile 3: 81 .. 116
           'Lys', 'Glu', 'Pro', 'Pro', 'Gln', 'Gly', 'Ala', 'His', 'Phe', 'Leu', 'Ser', 'Arg', 'Ser', 'Leu', 'Asp',
           'Asp', 'Ala', 'Leu', 'Lys', 'Leu', 'Thr', 'Glu', 'Gln', 'Pro', 'Glu', 'Leu', 'Ala', 'Asn', 'Lys', 'Val',
           'Asp', 'Met', 'Val', 'Trp', 'Ile', 'Val',
           #gap from 117 .. 118
           'Gly', 'Gly', 
           #tile 4: 119 .. 152
           'Ser', 'Ser', 'Val', 'Tyr', 'Lys', 'Glu', 'Ala', 'Met', 'Asn', 'His', 'Pro', 'Gly',
           'His', 'Leu', 'Lys', 'Leu', 'Phe', 'Val', 'Thr', 'Arg', 'Ile', 'Met', 'Gln', 'Asp', 'Phe', 'Glu', 'Ser',
           'Asp', 'Thr', 'Phe', 'Phe', 'Pro', 'Glu', 'Ile', 'Asp', 'Leu', 'Glu', 'Lys',
           #tile 5: 153 .. 187
           'Tyr', 'Lys', 'Leu', 'Leu',
           'Pro', 'Glu', 'Tyr', 'Pro', 'Gly', 'Val', 'Leu', 'Ser', 'Asp', 'Val', 'Gln', 'Glu', 'Glu', 'Lys', 'Gly',
           'Ile', 'Lys', 'Tyr', 'Lys', 'Phe', 'Glu', 'Val', 'Tyr', 'Glu', 'Lys', 'Asn', 'Asp'
           )


p1 <- make_heatmap(d = all_df, start_pos =  start_pos, tile_end = tile_end, pos_WT =  pos_WT, aa_WT = aa_WT,
                   min_score = -2.5, max_score = 4.0,
                   descr = paste0("MAVE scores all tiles, 0 mismatch allowed, min base qual 1, average read qual 1",
                                  "\nonly DNA vars with minimum 10 read counts, ", corr_descr, ", mid point: sy score = 0"))

p1

png(filename = paste0(base_dir,'HZ_plots/heatmaps/heatmap_alltiles.minvarcount10',corr,'.png'), 
    width = 3000, height = 800, res = 100)
print(p1)
dev.off()

