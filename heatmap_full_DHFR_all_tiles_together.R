library(pheatmap)
library(stringr)
library(matrixStats)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(hash) #to store infos for the different tiles
library(egg)

#heatmap of MAVE scores full DHFR seq all tiles together
base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/results/enrich_full_dhfr/'
data_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/data/'

corr_descr <- 'corrected with synonymous counts' 
corr <- '_corr_sy'

#load prepared data
########################

all_df <- read.csv(paste0(base_dir, 'combined_R_dataframes/', 'all_tile_scores_long.csv'), 
                   sep = ',')

#all_df including the vars only scored in two selections instead of 3
all_df_new <- read.csv(paste0(base_dir, 'combined_R_dataframes/', 'all_tile_scores_long_plus2sels.csv'), 
                   sep = ',')

#load distance data
########################
mtx_dist <- read.csv(paste0(data_dir,'1u72_dist_to_MTX.dat'), 
                     sep = '\t', header = F,
                     col.names = c('pos_raw','chain', 'WTres', 'MTXdist'))

#put in correct position starting from 2
mtx_dist$pos <- mtx_dist$pos_raw+1

#remove distances of waters (empty residue field)
#mtx_dist <- mtx_dist[complete.cases(mtx_dist),]
mtx_dist <- mtx_dist[mtx_dist$WTres != '',]
#mtx_dist$molecule <- rep('MTX', nrow(mtx_dist))

nadph_dist <- read.csv(paste0(data_dir,'1u72_dist_to_NDP.dat'), 
                       sep = '\t', header = F,
                       col.names = c('pos_raw','chain', 'WTres', 'NADPHdist'))

#put in correct position starting from 2
nadph_dist$pos <- nadph_dist$pos_raw+1

#remove distances of waters (empty residue field)
#nadph_dist <- nadph_dist[complete.cases(nadph_dist),]
nadph_dist <- nadph_dist[nadph_dist$WTres != '',]
#nadph_dist$molecule <- rep('NADPH', nrow(nadph_dist))

#adding secondary structure elements
###############################

dssp_dat <- read.csv(paste0(data_dir,'prism_dssp_DHFR_P00374_1u72_parsed.txt'), 
                     sep = ' ', header = T, comment.char = '#')
dssp_dat$pos <- seq(start_pos,tile_end)
#do I need this in order to plot???
dssp_dat$y <- rep(1,nrow(dssp_dat))

#define the WT to mark in the heatmap
####################################
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


#make heatmap function
#########################
make_heatmap <- function(d, start_pos, tile_end, pos_WT, aa_WT, min_score=0,max_score=0, descr='', x_axis = T) {
  
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
    
    #which x-axis ticks to draw. It seems I need to use 1 - tile_end as limits, otherwise it just displays no x-axis scale.
    #the labels and breaks are me trying to get it to display only labels between tile start and tile end, not all the way from 1. 
    scale_x_discrete(name ="pos", breaks=seq(start_pos,tile_end,1), labels=as.character(seq(start_pos,tile_end,1)), 
                     limits=as.factor(seq(1,tile_end,1))) +
    
    ylab('substitution') +
    #add WT markers. Lets try drawing black circles
    geom_tile(data = WT, aes(x = pos_WT, y=aa_WT, fill = dummy), color = 'grey50') + #scale_color_manual(values = c("white")) +
    geom_point(data = WT, aes(x = pos_WT, y=aa_WT), size = 3, color = 'black') +
    
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
  
  if(!x_axis) {
    p <- p + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

  }
  
  
  
  return(p)
  
}


#orig heatmap of all tiles
p1 <- make_heatmap(d = all_df, start_pos =  start_pos, tile_end = tile_end, pos_WT =  pos_WT, aa_WT = aa_WT,
                   min_score = -2.5, max_score = 4.0,
                   descr = paste0("MAVE scores all tiles, 0 mismatch allowed, min base qual 1, average read qual 1",
                                  "\nonly DNA vars with minimum 10 read counts, ", corr_descr, ", mid point: sy score = 0"))

p1

#png(filename = paste0(base_dir,'HZ_plots/heatmaps/heatmap_alltiles.minvarcount10',corr,'.png'), 
#    width = 3000, height = 800, res = 100)
#print(p1)
#dev.off()

#heatmap of all tiles including vars scored only in two selection instead of 3. I got those by manually re-running enrich with only 2 
#sels at a time and adding the vars to a copy of the all_df, called all_df_new
p2 <- make_heatmap(d = all_df_new, start_pos =  start_pos, tile_end = tile_end, pos_WT =  pos_WT, aa_WT = aa_WT,
                   min_score = -2.5, max_score = 4.0,
                   descr = paste0("MAVE scores all tiles, including vars only scored in two replicates, 0 mismatch allowed, min base qual 1, average read qual 1",
                                  "\nonly DNA vars with minimum 10 read counts, ", corr_descr, ", mid point: sy score = 0"))

p2

#add mtx and nad distance plots and arrange
mtx_p <- ggplot(data = mtx_dist, aes(x=pos,y=MTXdist)) +
  geom_line() +
  scale_x_discrete(name ="pos", breaks=seq(start_pos,tile_end,1), labels=as.character(seq(start_pos,tile_end,1)), 
                   limits=as.factor(seq(1,tile_end,1))) +
  #limit the plot from 2 to 187 (omit the pos 1 in the plot but it still needs to be in the axis for scaling)
  coord_cartesian(xlim = c(start_pos,tile_end)) +
  theme_bw(base_size = 16) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

mtx_p


nad_p <- ggplot(data = nadph_dist, aes(x=pos,y=NADPHdist)) +
  geom_line() +
  scale_x_discrete(name ="pos", breaks=seq(start_pos,tile_end,1), labels=as.character(seq(start_pos,tile_end,1)), 
                   limits=as.factor(seq(1,tile_end,1))) +
  #limit the plot from 2 to 187 (omit the pos 1 in the plot but it still needs to be in the axis for scaling)
  coord_cartesian(xlim = c(start_pos,tile_end)) +
  theme_bw(base_size = 16) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=18), legend.text=element_text(size=18), 
        plot.title = element_text(size=22), axis.text.x = element_text(angle = 90))


nad_p

#adding plot of sse's
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
sse_p <- ggplot(dssp_dat, aes(x=pos,y=y,fill=SS)) +
  geom_tile() +
  scale_fill_manual(values=cbPalette) +
  scale_x_discrete(name ="pos", breaks=seq(start_pos,tile_end,1), labels=as.character(seq(start_pos,tile_end,1)), 
                   limits=as.factor(seq(1,tile_end,1))) +
  #limit the plot from 2 to 187 (omit the pos 1 in the plot but it still needs to be in the axis for scaling)
  coord_cartesian(xlim = c(start_pos,tile_end)) +
  theme_bw(base_size = 16) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=18), legend.text=element_text(size=18), 
        plot.title = element_text(size=22), axis.text.x = element_text(angle = 90),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

sse_p

#orig heatmap with dist plots
png(filename = paste0(base_dir,'HZ_plots/heatmaps/heatmap_alltiles.minvarcount10',corr,'_dist_sse.png'), 
    width = 4500, height = 1500, res = 120)
ggarrange(
  p1, sse_p, mtx_p, nad_p,
  nrow = 4, ncol = 1, heights = c(5,1,1,1)
)
dev.off()

#orig heatmap with dist plots
png(filename = paste0(base_dir,'HZ_plots/heatmaps/heatmap_alltiles.minvarcount10',corr,'_dist.png'), 
    width = 4500, height = 1500, res = 120)
ggarrange(
  p1, mtx_p, nad_p,
  nrow = 3, ncol = 1, heights = c(5,1,1)
)
dev.off()


#heatmap including the vars scored in only 2 selections, with dist plots
png(filename = paste0(base_dir,'HZ_plots/heatmaps/heatmap_alltiles.minvarcount10',corr,'_plus2sels_dist.png'), 
    width = 4500, height = 1500, res = 120)
ggarrange(
  p2, mtx_p, nad_p,
  nrow = 3, ncol = 1, heights = c(5, 1,1)
)
dev.off()

#heatmap including the vars scored in only 2 selections, with dist plots
png(filename = paste0(base_dir,'HZ_plots/heatmaps/heatmap_alltiles.minvarcount10',corr,'_plus2sels_dist_sse.png'), 
    width = 4500, height = 1500, res = 120)
ggarrange(
  p2, sse_p, mtx_p, nad_p,
  nrow = 4, ncol = 1, heights = c(5,1,1,1)
)
dev.off()






