library(ggplot2)
library(stringr)
library(reshape2)
library(matrixStats)
library(ggrepel)
library(dplyr) #arrange
library(hash)

#load in all MAVE scores and prepare and save two dfs:
#1. A long one listing the score and SE for each variant, plus pos, sub, tile. 
#Each comb of sub and pos can only occur once, so for vars scored in two tiles the score is the average and the error is found with error propagation between the two orig. errors
#2. A matrix of scores with the cols being subs and the rows being positions. Also has pos_median, position and tile info
#added some more, see below

base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/results/enrich_full_dhfr/'
data_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/data/'
corr_descr <- 'corrected with synonymous counts' # #'corrected with WT counts' # 'corrected with all reads' # 'corrected with counted reads'
corr <- '_corr_sy' #'_corr_complete' #'' #'_corr_full' #'_corr_complete'

order_aa <- c('His', 'Lys', 'Arg', 'Asp', 'Glu', 'Cys', 'Met', 'Asn', 'Gln', 'Ser', 'Thr', 'Ala', 
              'Ile', 'Leu', 'Val', 'Phe', 'Trp', 'Tyr', 'Gly', 'Pro', 'Ter', 'pos_median')

#Store infos on the location and WT seq of each tile in a hash (dictionary) 
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
  h[['4']][['start_pos']] <- 119
  h[['4']][['tile_end']] <- 156
  h[['4']][['WT_seq']] <- c('Val', 'Gly', 'Gly', 'Ser', 'Ser', 'Val', 'Tyr', 'Lys', 'Glu', 'Ala', 'Met', 'Asn', 'His', 'Pro', 'Gly',
                            'His', 'Leu', 'Lys', 'Leu', 'Phe', 'Val', 'Thr', 'Arg', 'Ile', 'Met', 'Gln', 'Asp', 'Phe', 'Glu', 'Ser',
                            'Asp', 'Thr', 'Phe', 'Phe', 'Pro', 'Glu', 'Ile', 'Asp', 'Leu', 'Glu', 'Lys', 'Tyr', 'Lys', 'Leu')
  
  h[["5"]] <- hash()
  #the mutagenized region is 153 - 187 but we have read overlap (from fw rv pair) from aa 146 til the stop codon (pos after 187)
  h[['5']][['start_pos']] <- 153
  h[['5']][['tile_end']] <- 187
  h[['5']][['WT_seq']] <- c('Asp', 'Thr', 'Phe', 'Phe', 'Pro', 'Glu', 'Ile', 'Asp', 'Leu', 'Glu', 'Lys', 'Tyr', 'Lys', 'Leu', 'Leu',
                            'Pro', 'Glu', 'Tyr', 'Pro', 'Gly', 'Val', 'Leu', 'Ser', 'Asp', 'Val', 'Gln', 'Glu', 'Glu', 'Lys', 'Gly',
                            'Ile', 'Lys', 'Tyr', 'Lys', 'Phe', 'Glu', 'Val', 'Tyr', 'Glu', 'Lys', 'Asn', 'Asp')
  
  return(h)
}
h <- make_tile_data()

#prep Mave data
#put all MAVE scores into one long df
##########################
#all_df <- data.frame()
tile <- '1'

all_df <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps_minvarcount10',corr,'/tsv/tile',tile ,'_all_reps_exp/main_synonymous_scores.tsv'), 
                   sep = '\t', header = F,
                   stringsAsFactors = F, skip = 4,
                   col.names = c('var', 'SE', 'epsilon', 'score'))
all_df$tile <- '1'


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
  dat$tile <- tile
  all_df <- rbind(all_df, dat)
  
}

#necessary prep for making the all tiles map

#drop multi variants
all_df <- all_df[!grepl(',',all_df$var, fixed = T),-which(names(all_df) %in% c('epsilon'))]

#add position and sub
all_df$pos <- as.numeric(str_extract(all_df$var, "([0-9]+)"))
all_df$sub <- str_match(all_df$var, "[0-9]+([A-Za-z]+)")[,2]

#to avoid problems later on, every sub can only occur once per position. This means I need to remove the overlapping vars, or pick one
#or do the average
#example: all_df[all_df$var == 'p.Asp153Ala',]

#find combos of sub and position with two score values (becaues they were seen in multiple tiles)
combos <- table(interaction(all_df$pos, all_df$sub))
doubles <- combos[combos == 2]

#debug
#a <- strsplit(item, '.', fixed=T)
#debug

#calc the average between the two scores and replace them with that average score
for (item in names(doubles)) {
  
  #print('\n\n')
  #print(item)
  #split the name into the position and the subs so we can adress the relevant lines in the df
  a <- strsplit(item, '.', fixed=T)
  #test if we get the correct lines. Ok
  #print(all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], ])
  
  #which tiles is this var in?
  double_pos = as.numeric(a[[1]][[1]])
  
  double_tiles <- c()
  for (tile in c('1','2','3','4','5')) {
    #print(c('tile', tile))
    if (double_pos >= h[[tile]][['start_pos']] & double_pos <= h[[tile]][['tile_end']]){
      double_tiles <- append(double_tiles, tile)
    }
  }
  
  var_name <- unique(all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], c('var')])
  #calc a mean over all scores observed for this sub and position (should only be 2 since a position is in max two different tiles)
  mean_score <- mean(all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], c('score')])
  #error progagation:
  new_error <- 1/2 * sqrt(all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], c('SE')][1]**2 + 
                            all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], c('SE')][2]**2)
  
  #print(mean_score)
  
  #remove the rows and add a new one with the mean
  all_df <- all_df[!(all_df$pos == as.numeric(a[[1]][[1]]) & all_df$sub == a[[1]][[2]]), ]
  new_row <- data.frame(var=var_name,score=mean_score,pos=as.numeric(a[[1]][[1]]),sub= a[[1]][[2]], SE = new_error, tile = paste(double_tiles, collapse = '+')) 
  print('new_row is:')
  print(new_row)
  all_df <- rbind(all_df, new_row)
  #print(all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], ])
}

all_df <- arrange(all_df, pos)

#MAVE scores across all tiles but per selection. Vars scored in two tiles averaged.
##########################################

rep <- 'A'
for (rep in c('A','B','C')) {
  #read in tile 1
  dat1 <- read.csv(paste0(base_dir,'tile1_all_reps_minvarcount10', corr,'/tsv/',rep,'_sel/main_synonymous_scores.tsv'), 
                   sep = '\t', header = F,
                   stringsAsFactors = F, skip = 3,
                   col.names = c('var', 'score',	'SE',	'SE_pctile',	'slope',	'intercept',	'SE_slope',	't', 'pvalue_raw'))
  dat1$tile <- '1'
  
  #read in all the other tiles
  for (tile in c('2','3','4','5')) {
    #print(c('tile', tile))
    
    #we will use wt correction, except that we replaced wt counts with sy counts since the sy score is much closer to the 
    #peak of the score distr (most vars behave similar to it). Therefore, all scores will be on the 'same' scale where 0 is 
    #performance like wt/sy. Therefore, we do not need to read in more wt scores (they are all 0) or more sy scores (they are all NA because of our trick with replacing the wt counts with sy counts)
    #therefore we skip 4 lines instead of 2 when reading in the other tiles
    dat <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps_minvarcount10', corr,'/tsv/',rep,'_sel/main_synonymous_scores.tsv'), 
                    sep = '\t', header = F,
                    stringsAsFactors = F, skip = 3,
                    col.names = c('var', 'score',	'SE',	'SE_pctile',	'slope',	'intercept',	'SE_slope',	't', 'pvalue_raw'))
    dat$tile <- tile
    dat1 <- rbind(dat1, dat)
    
  }
  
  #drop multi variants
  dat1 <- dat1[!grepl(',',dat1$var, fixed = T),-which(names(dat1) %in% c('SE_pctile', 'slope', 'intercept', 'SE_slope', 't', 'pvalue_raw'))]
  
  #add position and sub
  dat1$pos <- as.numeric(str_extract(dat1$var, "([0-9]+)"))
  dat1$sub <- str_match(dat1$var, "[0-9]+([A-Za-z]+)")[,2]
  
  #to avoid problems later on, every sub can only occur once per position. This means I need to remove the overlapping vars, or pick one
  #or do the average
  #example: dat1[dat1$var == 'p.Asp153Ala',]
  
  #find combos of sub and position with two score values (becaues they were seen in multiple tiles)
  combos <- table(interaction(dat1$pos, dat1$sub))
  doubles <- combos[combos == 2]
  
  #calc the average between the two scores and replace them with that average score
  for (item in names(doubles)) {
    
    #print('\n\n')
    #print(item)
    #split the name into the position and the subs so we can adress the relevant lines in the df
    a <- strsplit(item, '.', fixed=T)
    #test if we get the correct lines. Ok
    #print(dat1[dat1$pos == as.numeric(a[[1]][[1]]) &  dat1$sub == a[[1]][[2]], ])
    
    #which tiles is this var in?
    double_pos = as.numeric(a[[1]][[1]])
    
    double_tiles <- c()
    for (tile in c('1','2','3','4','5')) {
      #print(c('tile', tile))
      if (double_pos >= h[[tile]][['start_pos']] & double_pos <= h[[tile]][['tile_end']]){
        double_tiles <- append(double_tiles, tile)
      }
    }
    
    var_name <- unique(dat1[dat1$pos == as.numeric(a[[1]][[1]]) &  dat1$sub == a[[1]][[2]], c('var')])
    #calc a mean over all scores observed for this sub and position (should only be 2 since a position is in max two different tiles)
    mean_score <- mean(dat1[dat1$pos == as.numeric(a[[1]][[1]]) &  dat1$sub == a[[1]][[2]], c('score')])
    #error progagation:
    new_error <- 1/2 * sqrt(dat1[dat1$pos == as.numeric(a[[1]][[1]]) &  dat1$sub == a[[1]][[2]], c('SE')][1]**2 + 
                              dat1[dat1$pos == as.numeric(a[[1]][[1]]) &  dat1$sub == a[[1]][[2]], c('SE')][2]**2)
    
    #print(mean_score)
    
    #remove the rows and add a new one with the mean
    dat1 <- dat1[!(dat1$pos == as.numeric(a[[1]][[1]]) & dat1$sub == a[[1]][[2]]), ]
    new_row <- data.frame(var=var_name,score=mean_score,pos=as.numeric(a[[1]][[1]]),sub= a[[1]][[2]], SE = new_error, tile = paste(double_tiles, collapse = '+')) 
    print('new_row is:')
    print(new_row)
    dat1 <- rbind(dat1, new_row)
    #print(dat1[dat1$pos == as.numeric(a[[1]][[1]]) &  dat1$sub == a[[1]][[2]], ])
  }
  
  write.csv(dat1, file = paste0(base_dir, 'combined_R_dataframes/', 'all_tile_scores_long_rep', rep, '.csv'), quote = F, row.names = F)

}

#matrix with position based median MAVE scores
###########################################################

mat <- dcast(all_df, pos ~ sub, value.var = "score")
rownames(mat) <- mat$pos
#remove pos column
mat <- mat[,-c(1)]
mat$pos_median = rowMedians(as.matrix(mat), na.rm = T)

start_pos <- 2
tile_end <- 187
all_rows <- seq(start_pos,tile_end)
#missing rows:
missing_rows <- all_rows[!all_rows %in% rownames(mat)]

for (i in missing_rows) {
  mat[as.character(i),] <-  rep(NA,ncol(mat))
}

#order first, then add cols
mat <- mat[order(as.numeric(row.names(mat))),order_aa]
mat$pos <- as.numeric(rownames(mat))
mat$tile <- ifelse(mat$pos >= h[['1']][['start_pos']] &  mat$pos <= h[['1']][['tile_end']], 'tile1', 
                   ifelse(mat$pos >= h[['2']][['start_pos']] &  mat$pos <= h[['2']][['tile_end']], 'tile2', 
                          ifelse(mat$pos >= h[['3']][['start_pos']] &  mat$pos <= h[['3']][['tile_end']], 'tile3',
                                 ifelse(mat$pos >= h[['4']][['start_pos']] &  mat$pos <= h[['4']][['tile_end']], 'tile4', 'tile5'))))

#########################################

#save all_df and mat as csvs for other plots
write.csv(all_df, file = paste0(base_dir, 'combined_R_dataframes/', 'all_tile_scores_long.csv'), quote = F, row.names = F)
write.csv(mat, file = paste0(base_dir, 'combined_R_dataframes/', 'per_pos_matrix.csv'), quote = F, row.names = F)

