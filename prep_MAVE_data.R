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

#dict to translate 1 letter aa code to 3 letters
make_aa_hash <- function(){
  aa_1_to_3 <- hash()
  aa_1_to_3[["A"]] <- 'Ala'
  aa_1_to_3[["C"]] <- 'Cys'
  aa_1_to_3[["D"]] <- 'Asp'
  aa_1_to_3[["E"]] <- 'Glu'
  aa_1_to_3[["F"]] <- 'Phe'
  aa_1_to_3[["G"]] <- 'Gly'
  aa_1_to_3[["H"]] <- 'His'
  aa_1_to_3[["I"]] <- 'Ile'
  aa_1_to_3[["K"]] <- 'Lys'
  aa_1_to_3[["L"]] <- 'Leu'
  aa_1_to_3[["M"]] <- 'Met'
  aa_1_to_3[["N"]] <- 'Asn'
  aa_1_to_3[["P"]] <- 'Pro'
  aa_1_to_3[["Q"]] <- 'Gln'
  aa_1_to_3[["R"]] <- 'Arg'
  aa_1_to_3[["S"]] <- 'Ser'
  aa_1_to_3[["T"]] <- 'Thr'
  aa_1_to_3[["V"]] <- 'Val'
  aa_1_to_3[["W"]] <- 'Trp'
  aa_1_to_3[["Y"]] <- 'Tyr'
  return(aa_1_to_3)
}
aa_1_to_3 <- make_aa_hash()

#and reverse, translate 3 letter to 1 letter
make_aa_hash_rv <- function(){
  aa_3_to_1 <- hash()
  aa_3_to_1[["Ala"]] <- 'A'
  aa_3_to_1[["Cys"]] <- 'C'
  aa_3_to_1[["Asp"]] <- 'D'
  aa_3_to_1[["Glu"]] <- 'E'
  aa_3_to_1[["Phe"]] <- 'F'
  aa_3_to_1[["Gly"]] <- 'G'
  aa_3_to_1[["His"]] <- 'H'
  aa_3_to_1[["Ile"]] <- 'I'
  aa_3_to_1[["Lys"]] <- 'K'
  aa_3_to_1[["Leu"]] <- 'L'
  aa_3_to_1[["Met"]] <- 'M'
  aa_3_to_1[["Asn"]] <- 'N'
  aa_3_to_1[["Pro"]] <- 'P'
  aa_3_to_1[["Gln"]] <- 'Q'
  aa_3_to_1[["Arg"]] <- 'R'
  aa_3_to_1[["Ser"]] <- 'S'
  aa_3_to_1[["Thr"]] <- 'T'
  aa_3_to_1[["Val"]] <- 'V'
  aa_3_to_1[["Trp"]] <- 'W'
  aa_3_to_1[["Tyr"]] <- 'Y'
  aa_3_to_1[["Ter"]] <- '*'
  return(aa_3_to_1)
}
aa_3_to_1 <- make_aa_hash_rv()

##########################################
#prep Mave data
#1. put all MAVE scores into one long df
##########################
#all_df <- data.frame()
tile <- '1'

all_df <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps_minvarcount10',corr,'/tsv/tile',tile ,'_all_reps_exp/main_synonymous_scores.tsv'), 
                   sep = '\t', header = F,
                   stringsAsFactors = F, skip = 4,
                   col.names = c('var', 'SE', 'epsilon', 'MAVE_score'))
all_df$tile <- '1'

p_vals <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps_minvarcount10',corr,'/tsv/tile',tile ,'_all_reps_exp/main_synonymous_scores_pvalues_wt.tsv'), 
                   sep = '\t', header = F,
                   stringsAsFactors = F, skip = 3,
                   col.names = c('var', 'p_raw', 'z'))

all_df <- merge(all_df, p_vals, by = 'var')


for (tile in c('2','3','4','5')) {
  print(c('tile', tile))
  
  #we will use wt correction, except that we replaced wt counts with sy counts since the sy score is much closer to the 
  #peak of the score distr (most vars behave similar to it). Therefore, all scores will be on the 'same' scale where 0 is 
  #performance like wt/sy. Therefore, we do not need to read in more wt scores (they are all 0) or more sy scores (they are all NA because of our trick with replacing the wt counts with sy counts)
  #therefore we skip 4 lines instead of 2 when reading in the other tiles
  dat <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps_minvarcount10',corr,'/tsv/tile',tile ,'_all_reps_exp/main_synonymous_scores.tsv'), 
                  sep = '\t', header = F,
                  stringsAsFactors = F, skip = 4,
                  col.names = c('var', 'SE', 'epsilon', 'MAVE_score'))
  dat$tile <- tile
  
  p_vals <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps_minvarcount10',corr,'/tsv/tile',tile ,'_all_reps_exp/main_synonymous_scores_pvalues_wt.tsv'), 
                     sep = '\t', header = F,
                     stringsAsFactors = F, skip = 3,
                     col.names = c('var', 'p_raw', 'z'))
  
  dat <- merge(dat, p_vals, by = 'var')
  
  all_df <- rbind(all_df, dat)
  
}

#now correct the raw p-values with the number of tests. This number of tests includes multi var which actually we're not interested in
all_df$p_bonf <- p.adjust(all_df$p_raw, 'bonferroni')

#drop multi variants
all_df <- all_df[!grepl(',',all_df$var, fixed = T),-which(names(all_df) %in% c('epsilon'))]

#add position and sub
all_df$pos <- as.numeric(str_extract(all_df$var, "([0-9]+)"))
#all_df$sub <- trans_3_to_1(str_match(all_df$var, "[0-9]+([A-Za-z]+)")[,2])
#all_df$WT <- str_match(all_df$var, "([A-Za-z]+)[0-9]+")[,2]
#change to 1 letter code for sub and WT and create a prism style var name column
#for some reason all_df$sub became a weird nested list, so I'm casting as.character around it just to be sure the values become strings
all_df$sub <- as.character(sapply(str_match(all_df$var, "[0-9]+([A-Za-z]+)")[,2], function(i) aa_3_to_1[[i]]))
all_df$WT <- as.character(sapply(str_match(all_df$var, "([A-Za-z]+)[0-9]+")[,2], function(i) aa_3_to_1[[i]]))
all_df$prism_var <- paste0(all_df$WT, all_df$pos, all_df$sub)
all_df$sub_long <- str_match(all_df$var, "[0-9]+([A-Za-z]+)")[,2]

#drop orig var name column
#all_df <- subset(all_df, select = -c(var))

#to avoid problems later on, every sub can only occur once per position. This means I need to remove the overlapping vars, or pick one
#or do the average
#example: all_df[all_df$var == 'p.Asp153Ala',]

#find combos of sub and position with two score values (becaues they were seen in multiple tiles)
combos <- table(interaction(all_df$pos, all_df$sub))
doubles <- combos[combos == 2]

#debug: get just one of the vars with two scores 
a <- strsplit(names(doubles)[1], '.', fixed=T)
#debug

c <- 0
#calc the average between the two scores and replace them with that average score
for (item in names(doubles)) {
  
  c <- c+1
  print('')
  print(item)
  #split the name into the position and the subs so we can adress the relevant lines in the df
  a <- strsplit(item, '.', fixed=T)
  
  #test if we get the correct lines. Ok
  print(all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], ])
  
  #which tiles is this var in?
  double_pos = as.numeric(a[[1]][[1]])
  print(double_pos)
  double_tiles <- c()
  for (tile in c('1','2','3','4','5')) {
    #print(c('tile', tile))
    if (double_pos >= h[[tile]][['start_pos']] & double_pos <= h[[tile]][['tile_end']]){
      double_tiles <- append(double_tiles, tile)
    }
  }
  print(double_tiles)
  
  #the new line we want to replace the two orig lines with has the same values except those that need to be re-calc:
  #MAVE_score, SE (error propagation), tile, p_raw, z score, p_bonf
  
  #calc a mean over all scores observed for this sub and position (should only be 2 since a position is in max two different tiles)
  mean_score <- mean(all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], c('MAVE_score')])
  #error progagation:
  new_error <- 1/2 * sqrt(all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], c('SE')][1]**2 + 
                            all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], c('SE')][2]**2)
  #I should recalc the p-value with the new mean score and propagated error but I don't know how. I'll just make an average instead
  p_new =  mean(all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], c('p_raw')])
  bonf_new =  mean(all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], c('p_bonf')])
  z_new =  mean(all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], c('z')])
  
  #so the template for the new row will be the first of all rows that match this pos and sub (should only be two if the var is in 2 tiles)
  new_row <- all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], ][1,]
  #replace with the re-calc'ed values
  new_row$MAVE_score <- mean_score
  new_row$SE <- new_error
  new_row$p_raw <- p_new
  new_row$p_bonf <- bonf_new
  new_row$z <- z_new
  new_row$tile <- paste(double_tiles, collapse = '+')
  
  #remove the old rows
  all_df <- all_df[!(all_df$pos == as.numeric(a[[1]][[1]]) & all_df$sub == a[[1]][[2]]), ]
  #test that the rows have disappeared
  print(all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], ])
  
  #add the new row:
  print('new_row is:')
  print(new_row)
  all_df <- rbind(all_df, new_row)
  
  #verify the new row has been added
  print(all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], ])
  
}  

#re-sort the df rows so that the newly added rows go where they belong position-wise 
all_df <- arrange(all_df, pos)

#drop orig var name column and reorder so prism_var is the first
#all_df <- subset(all_df, select = -c(var))
#all_df <- all_df[, c('prism_var', 'SE', 'MAVE_score', 'tile', 'p_raw', 'pos', 'sub', 'WT')]

#save all_df and mat as csvs for other plots
#specify the order of columns and which ones to write with []
write.csv(all_df[, c('prism_var', 'var', 'SE', 'MAVE_score', 'tile', 'p_raw', 'p_bonf', 'pos', 'sub', 'WT')], 
          file = paste0(base_dir, 'combined_R_dataframes/', 'all_tile_scores_long.csv'), quote = F, row.names = F)

#I'll instead make a raw prism file of the data including vars scored in 2 selection
#for making a prism file, omit the pos, sub and WT cols (the prism parser will add these automatically) and 
#use space as the separator (write.table instead of write.csv)
#write.table(all_df[,c("prism_var", "SE", "MAVE_score", "tile", "p_raw","p_bonf")], 
#            file = paste0(base_dir, 'combined_R_dataframes/', 'raw_prism.csv'), quote = F, row.names = F, sep = ' ')


##########################################
#2. add vars that only have scores in two instead of all three selections/replicates (but they were seen in very time point of these two selections) 
##########################################

#subset all_df to only relevant columns for continuing
all_df <- all_df[, c('prism_var', 'var', 'SE', 'MAVE_score', 'tile', 'p_raw', 'p_bonf', 'pos', 'sub', 'WT')]
#remove correct pvals because we need to re-calc that with the additional p-values of the newly added vars
all_df_new <- subset(all_df, select = -c(p_bonf))
all_df_new$sels_scored_in <- 3

for (tile in c('1','2', '3','4','5')) {
  print(tile)
  for (comb in c('A_B', 'B_C', 'A_C')) {
    
    #open the re-runs I made with only two selections in the config file. 
    #ther are saved under i.e. enrich_full_dhfr/tile1_all_reps_minvarcount10_corr_sy_A_B
    data_path <- paste0(base_dir,'tile', tile ,'_all_reps_minvarcount10',corr,'_',comb,'/tsv/tile',
                        tile ,'_all_reps_exp/main_synonymous_scores.tsv')
    
    #I think there was some problem that one of these might not exist
    if(file.exists(data_path)){
      #two_sels is the main_synonymous_scores.tsv file of the current tile and selection combination, i.e. A_B
      two_sels <- read.csv(data_path, 
                           sep = '\t', header = F,
                           stringsAsFactors = F, skip = 4,
                           col.names = c('var', 'SE', 'epsilon', 'MAVE_score'))
      
      #add p-values
      p_val_path <- paste0(base_dir,'tile', tile ,'_all_reps_minvarcount10',corr,'_',comb,'/tsv/tile',
                           tile ,'_all_reps_exp/main_synonymous_scores_pvalues_wt.tsv')
      p_vals <- read.csv(p_val_path, 
                         sep = '\t', header = F,
                         stringsAsFactors = F, skip = 3,
                         col.names = c('var', 'p_raw', 'z'))
      
      two_sels <- merge(two_sels, p_vals[,c('var','p_raw')], by = 'var')
      
      #subset to only single vars
      two_sels <- two_sels[!grepl(',',two_sels$var, fixed = T),]
      #drop epsilon col and add tile col
      two_sels <- two_sels[,-which(names(two_sels) %in% c('epsilon'))]
      two_sels$tile <- tile
      two_sels$sels_scored_in <- 2
      
      #add pos, sub and WT cols to match the all_df cols
      two_sels$pos <- as.numeric(str_extract(two_sels$var, "([0-9]+)"))
      two_sels$sub <- as.character(sapply(str_match(two_sels$var, "[0-9]+([A-Za-z]+)")[,2], function(i) aa_3_to_1[[i]]))
      two_sels$WT <- as.character(sapply(str_match(two_sels$var, "([A-Za-z]+)[0-9]+")[,2], function(i) aa_3_to_1[[i]]))
      two_sels$prism_var <- paste0(two_sels$WT, two_sels$pos, two_sels$sub)
      
      #find rows/vars in two_sels that do not exist in all_df
      new_vars <- two_sels$var[!two_sels$var %in% all_df$var]
      
      print(comb)
      print(new_vars)
      
      all_df_new <- rbind(all_df_new, two_sels[two_sels$var %in% new_vars,]) 
      
    }
    
  }
  
}

#before saving the df, we need to take care of vars that appear multiple times. 
combos <- table(interaction(all_df_new$pos, all_df_new$sub))
doubles <- combos[combos == 2]

#checking:
#all_df_new[all_df_new$pos == 7 & all_df_new$sub == 'R',]

#re-add corrected p-values
all_df_new$p_bonf <- p.adjust(all_df_new$p_raw, 'bonferroni')

#writing out a df with only the new vars:
#the first rows of all_df_new are the already existing vars. The new ones have been added to the end, so I can get the rows from
#length(all_df) until len(all_df_new) and those will only be the newly added vars 
new_vars <- all_df_new[nrow(all_df)+1:nrow(all_df_new),]
write.csv(new_vars[, c('prism_var', 'var', 'SE', 'MAVE_score', 'tile', 'p_raw', 'p_bonf', 'pos', 'sub', 'WT')], 
          file = paste0(base_dir, 'combined_R_dataframes/', 'only2selvars.csv'), quote = F, row.names = F)

##########################################
#3. Categorize variants by their p-value against score = 0 and whether their MAVE is inside -2xSD(sy_vars) and +2xSD(sy_vars)
##########################################

#calc the standard deviation of MAVE scores of variants synonymous to the WT protein (but differnet on DNA level). 
#We're using main_variants_scores instead of main_synonymous_scores since we're looking at DNA level variants

#for debug
#tile <- 1

#dict to hold the 2xSD values for each tile
two_sigma <- hash()

for (tile in c('1', '2', '3', '4', '5')) {
  #load DNA level var data
  dna_dat <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps_minvarcount10',corr,'/tsv/tile',tile ,'_all_reps_exp/main_variants_scores.tsv'), 
                      sep = '\t', header = F,
                      stringsAsFactors = F, skip = 3,
                      col.names = c('var', 'SE', 'epsilon', 'MAVE_score'))
  
  #get the mutations into a list of only unique mutations
  #i.e. "c.6C>A (p.=), c.22A>T (p.Ile8Phe), c.24C>T (p.Ile8Phe)" should become "(p.=), (p.Ile8Phe)"
  #"c.6C>A (p.=), c.22A>T (p.=)" should be "(p.=)"
  dna_dat$uniq_muts <- sapply(str_extract_all(dna_dat$var, "[(]p[.].+?[)]"), function(i) unique(unlist(i)))
  #now, select only sy vars which are rows in which all muts are synonymous
  sy_vars <- dna_dat[dna_dat$uniq_muts == c('(p.=)'),]
  
  #get the standard dev
  two_sigma[[tile]] <- 2*sd(sy_vars$MAVE_score)
  
}

#check the values
two_sigma


#Add 2 sigma cutoff to dataframe
all_df_new$two_sigma <- ifelse(all_df_new$tile == '1', two_sigma[['1']], 
                           ifelse(all_df_new$tile == '2', two_sigma[['2']],
                                  ifelse(all_df_new$tile == '3', two_sigma[['3']],
                                         ifelse(all_df_new$tile == '4', two_sigma[['4']],
                                                ifelse(all_df_new$tile == '5', two_sigma[['5']],
                                                       ifelse(all_df_new$tile == '4+5', two_sigma[['4']],
                                                              NA))))))

#categorize variants according to pval and whether -2_sigma <= MAVE_score <= +2_sigma

#categories after raincloud_MAVE_categories_anno.ong
#cat 1: pval >= 0.05 and -2_sigma <= MAVE_score <= +2_sigma: WT-like
#cat 2: pval < 0.05 & MAVE_score < -2_sigma : del
#cat 3: pval < 0.05 & -2_sigma <= MAVE_score < 0 : WT-like del
#cat 4: pval < 0.05 & MAVE_score > 2_sigma : ben
#cat 5: pval < 0.05 & 0 < MAVE_score <= 2_sigma : WT-like ben
#cat 6: pval >= 0.05 and (MAVE_score < -2_sigma or MAVE_score > +2_sigma): non-WT but large error

all_df_new$classification <- ifelse(all_df_new$p_bonf >= 0.05 & -1*all_df_new$two_sigma <= all_df_new$MAVE_score & all_df_new$MAVE_score <= all_df_new$two_sigma, 'WT-like', 
                                ifelse(all_df_new$p_bonf < 0.05 & all_df_new$MAVE_score < -1*all_df_new$two_sigma, 'del',
                                       ifelse(all_df_new$p_bonf < 0.05 & -1*all_df_new$two_sigma <= all_df_new$MAVE_score & all_df_new$MAVE_score < 0, 'low_del',
                                              ifelse(all_df_new$p_bonf < 0.05 & all_df_new$MAVE_score > all_df_new$two_sigma, 'ben',
                                                     ifelse(all_df_new$p_bonf < 0.05 & 0 < all_df_new$MAVE_score & all_df_new$MAVE_score <= all_df_new$two_sigma, 'low_ben',
                                                            NA)))))


#print the whole df as csv
write.csv(all_df_new[, c('prism_var', 'var', 'SE', 'MAVE_score', 'tile', 'p_raw', 'p_bonf', 'pos', 'sub', 'WT', 'sels_scored_in', 'classification')], 
          file = paste0(base_dir, 'combined_R_dataframes/', 'all_tile_scores_long_plus2sels.csv'), quote = F, row.names = F)

#for making a prism file, omit the pos, sub and WT cols (the prism parser will add these automatically) and 
#use space as the separator (write.table instead of write.csv)
write.table(all_df_new[,c("prism_var", "SE", "MAVE_score", "tile", "p_raw","p_bonf", 'sels_scored_in', 'classification')], 
            file = paste0(base_dir, 'combined_R_dataframes/', 'raw_prism_plus2sels.csv'), quote = F, row.names = F, sep = ' ')




##########################################
#4. MAVE scores across all tiles but per selection. Vars scored in two tiles averaged.
##########################################

rep <- 'A'
for (rep in c('A','B','C')) {
  #read in tile 1
  dat1 <- read.csv(paste0(base_dir,'tile1_all_reps_minvarcount10', corr,'/tsv/',rep,'_sel/main_synonymous_scores.tsv'), 
                   sep = '\t', header = F,
                   stringsAsFactors = F, skip = 3,
                   col.names = c('var', 'MAVE_score',	'SE',	'SE_pctile',	'slope',	'intercept',	'SE_slope',	't', 'pvalue_raw'))
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
                    col.names = c('var', 'MAVE_score',	'SE',	'SE_pctile',	'slope',	'intercept',	'SE_slope',	't', 'pvalue_raw'))
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
    mean_score <- mean(dat1[dat1$pos == as.numeric(a[[1]][[1]]) &  dat1$sub == a[[1]][[2]], c('MAVE_score')])
    #error progagation:
    new_error <- 1/2 * sqrt(dat1[dat1$pos == as.numeric(a[[1]][[1]]) &  dat1$sub == a[[1]][[2]], c('SE')][1]**2 + 
                              dat1[dat1$pos == as.numeric(a[[1]][[1]]) &  dat1$sub == a[[1]][[2]], c('SE')][2]**2)
    
    #print(mean_score)
    
    #remove the rows and add a new one with the mean
    dat1 <- dat1[!(dat1$pos == as.numeric(a[[1]][[1]]) & dat1$sub == a[[1]][[2]]), ]
    new_row <- data.frame(var=var_name,MAVE_score=mean_score,pos=as.numeric(a[[1]][[1]]),sub= a[[1]][[2]], SE = new_error, tile = paste(double_tiles, collapse = '+')) 
    print('new_row is:')
    print(new_row)
    dat1 <- rbind(dat1, new_row)
    #print(dat1[dat1$pos == as.numeric(a[[1]][[1]]) &  dat1$sub == a[[1]][[2]], ])
  }
  
  write.csv(dat1, file = paste0(base_dir, 'combined_R_dataframes/', 'all_tile_scores_long_rep', rep, '.csv'), quote = F, row.names = F)

}

###########################################################
#5. matrix with position based median MAVE scores
###########################################################

mat <- dcast(all_df, pos ~ sub_long, value.var = "MAVE_score")
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

write.csv(mat, file = paste0(base_dir, 'combined_R_dataframes/', 'per_pos_matrix.csv'), quote = F, row.names = F)


##########################################
#6. prep and merge ecoli/human aligned dataframe with the other info.
#Add distance to ligands, ddE/G human DHFR, ddE/G ecoli DHFR 
##########################################

#load in the combined df I prepared in ipython notebook merge_dfs_human_ecoli.ipynb. 
#Positions that are aligned in the struct alignment are in the same row
#########################################
dat <- read.csv(paste0(data_dir,'merge_human_ecoli_fulldhfr.csv'), sep = '\t')
dat$pos <- as.numeric(str_extract(dat$human_variant, "([0-9]+)"))

dump_df <- dat[,c('human_variant', 'MAVE_score', 'SE', 'tile', 'pos', 'ecoli_variant', 'Lon_minus_score', 'Lon_minus_std_dev',
                  'Lon_minus_std_err', 'Lon_plus_score', 'Lon_plus_std_dev', 'Lon_plus_std_err')]

write.csv(dump_df, file = paste0(base_dir, 'combined_R_dataframes/', 'human_VS_ecoli_DHFR_data.csv'), quote = F, row.names = F)

#load in distance to NAPDH and MTX
############################

mtx_dist <- read.csv(paste0(data_dir,'1u72_dist_to_MTX.dat'), 
                     sep = '\t', header = F, col.names = c('pos_raw','chain', 'WTres', 'MTXdist'))

#put in correct position starting from 2
mtx_dist$pos <- mtx_dist$pos_raw+1

#remove distances of waters (empty residue field)
#mtx_dist <- mtx_dist[complete.cases(mtx_dist),]
mtx_dist <- mtx_dist[mtx_dist$WTres != '',]
mtx_dist$molecule <- rep('MTX', nrow(mtx_dist))

nadph_dist <- read.csv(paste0(data_dir,'1u72_dist_to_NDP.dat'), 
                       sep = '\t', header = F, col.names = c('pos_raw','chain', 'WTres', 'NAPDHdist'))

#put in correct position starting from 2
nadph_dist$pos <- nadph_dist$pos_raw+1

#remove distances of waters (empty residue field)
#nadph_dist <- nadph_dist[complete.cases(nadph_dist),]
nadph_dist <- nadph_dist[nadph_dist$WTres != '',]
nadph_dist$molecule <- rep('NADPH', nrow(nadph_dist))

#load merged human ddG/E data
###########################

dds <- read.csv(paste0(data_dir,'prism_merge_ddE_ddG_DHFR_P00374.txt'), comment.char = '#', sep = ' ')

#there are some rows at the end with ddGs for the non-vars (G16=). Remove those 
#for reasons I don't understand complete.cases doesn't want to operate on dds even though it is a dataframe, so I called it on the gemme score column instead
dds <- dds[complete.cases(dds$gemme_score_00),]
colnames(dds) <- c('variant', 'gemme_score_human', 'rosetta_ddg_score_human')

#load merged ecoli ddG/E data
###########################

dds_ecoli <- read.csv(paste0(data_dir,'prism_merge_ddE_ddG_ecoli_DHFR_P0ABQ4.txt'), comment.char = '#', sep = ' ')

#there are some rows at the end with ddGs for the non-vars (G16=). Remove those 
#for reasons I don't understand complete.cases doesn't want to operate on dds even though it is a dataframe, so I called it on the gemme score column instead
dds_ecoli <- dds_ecoli[complete.cases(dds_ecoli$gemme_score_00),]
colnames(dds_ecoli) <- c('variant', 'gemme_score_ecoli','gemme_coverage_ecoli', 'rosetta_ddg_score_ecoli', 'rosetta_ddg_std_ecoli')


#add distances and ddE/G to the human/ecoli aligned MAVE scores
############################

#distance to ligands is merged based on position in the human DHFR since this is what they were calculated on (in pymol on the human DHFR pdb struct)
all_merge <- merge(dump_df, mtx_dist[,c('MTXdist','pos')], by = 'pos', all.x = T)
all_merge <- merge(all_merge, nadph_dist[,c('NAPDHdist', 'pos')], by = 'pos', all.x = T)

#add some groupings for the distance to ligands (in human structure)
all_merge$dist_group_MTX <- ifelse(all_merge$MTXdist < 10,'MTX_close','MTX_far')
all_merge$dist_group_NADPH <- ifelse(all_merge$NAPDHdist < 10,'NADPH_close','NADPH_far')

#human ddG/E data is aligned to the human variant
all_merge <- merge(all_merge, dds, by.x = 'human_variant', by.y = 'variant', all.x = T)
#ecoli ddE/G data is aligned to the ecoli variant
all_merge <- merge(all_merge, dds_ecoli, by.x = 'ecoli_variant', by.y = 'variant', all.x = T)


# bla <- merge(dump_df[,c('human_variant', 'MAVE_score', 'SE', 'tile','ecoli_variant', 'Lon_minus_score', 'Lon_minus_std_dev',
#                         'Lon_minus_std_err', 'Lon_plus_score', 'Lon_plus_std_dev', 'Lon_plus_std_err')], 
#              dds, by.x = 'human_variant', by.y = 'variant', all = T)


write.csv(all_merge, file = paste0(base_dir, 'combined_R_dataframes/', 'human_VS_ecoli_DHFR_data_dds_dists.csv'), quote = F, row.names = F)

