library(ggplot2)
library(stringr)
library(hash)

#to plot the distributions of mave scores, ddE, ddG scores
#this is done with so called density plots

base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/results/enrich_full_dhfr/'
plot_dir <- paste0(base_dir, 'HZ_plots/score_densities/')

#which correction was employed and what should written in the plot headers
#the correction with the number of WT reads is default. It uses corr <- '' and corr_descr <- 'corrected with WT counts'.
#The corresponding enrich parameter is 'WT'
#correction with the counted reads showed better correlation in variants scored in both tiles 4 and 5
#it uses corr <- '_corr_complete' and corr_descr <- 'corrected with counted reads' . The corresponding enrich parameter is 'complete'
#There is also correction with the total number of reads. It gives very similar results to correcting with counted reads.

corr_descr <- 'corrected with synonymous counts' #'corrected with counted reads' #'corrected with WT counts' # 'corrected with all reads' # 'corrected with counted reads'
corr <- '_corr_sy' #'_corr_complete' #'' #'_corr_full' #'_corr_complete'

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

#plot distributions of standard errors per tile
############################

#separately
for (tile in c('1','2','3','4','5')) {
  print(c('tile', tile))

  dat10 <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps_minvarcount10/tsv/tile',tile ,'_all_reps_exp/main_synonymous_scores.tsv'), 
                    sep = '\t', header = F,
                    stringsAsFactors = F, skip = 2,
                    col.names = c('var', 'SE', 'epsilon', 'score'))
  
  p2 <- ggplot(dat10, aes(SE)) +
    geom_density(bw= 0.08) +
    xlim(0, 0.85) +
    ggtitle(paste0("SE tile", tile, ", min DNA var count 10")) +
    xlab('Bandwidth = 0.08') +
    theme_bw(base_size = 20) + theme(legend.position="bottom")
  p2
  
  png(filename = paste0(plot_dir,'error_tile',tile,'.png'),width = 800, height = 600)
  print(p2)
  dev.off()

}

#all tiles together
tile <- '1'
dat_all <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps_minvarcount10/tsv/tile',tile ,'_all_reps_exp/main_synonymous_scores.tsv'), 
                    sep = '\t', header = F,
                    stringsAsFactors = F, skip = 2,
                    col.names = c('var', 'SE', 'epsilon', 'score'))
dat_all$tile <- 'tile1'

for (tile in c('2','3','4','5')) {
  print(c('tile', tile))
  dat_temp <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps_minvarcount10/tsv/tile',tile ,'_all_reps_exp/main_synonymous_scores.tsv'), 
                       sep = '\t', header = F,
                       stringsAsFactors = F, skip = 2,
                       col.names = c('var', 'SE', 'epsilon', 'score'))
  dat_temp$tile <- paste0('tile',tile)
  dat_all <- rbind(dat_all, dat_temp)
}

p1 <- ggplot(dat_all, aes(SE,fill=as.factor(tile))) +
  geom_density(bw= 0.08, alpha = 0.5) +
  ggtitle(paste0("Errors per tile, min DNA var count 10")) +
  xlab('Bandwidth = 0.08') +
  theme_bw(base_size = 20) + theme(legend.position="bottom")

p1
png(filename = paste0(plot_dir,'alltiles.fill.error.png'),width = 800, height = 600)
print(p1)
dev.off()

p2 <- ggplot(dat_all, aes(SE,color=as.factor(tile))) +
  geom_density(bw= 0.08, alpha = 0.5) +
  ggtitle(paste0("Errors per tile, min DNA var count 10")) +
  xlab('Bandwidth = 0.08') +
  theme_bw(base_size = 20) + theme(legend.position="bottom")

p2
png(filename = paste0(plot_dir,'alltiles.error.png'),width = 800, height = 600)
print(p2)
dev.off()

#plot distribution of ddE and ddG scores
############################

data_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/data/'

ddG <- read.csv(paste0(data_dir,'prism_rosetta_xxx_DHFR_P00374_pdb4m6j_parsed.txt'), comment.char = '#', sep = ' ')
ddE <- read.csv(paste0(data_dir,'prism_gemme_999_DHFR_UNIPROT-kopi.txt'), comment.char = '#', sep = ' ')

ddG$aa_pos <- as.numeric(str_extract(ddG$variant, "(?<=[A-Z])([0-9]+)(?=[A-Z])"))
ddE$aa_pos <- as.numeric(str_extract(ddE$variant, "(?<=[A-Z])([0-9]+)(?=[A-Z])"))

#per tile: need to split up ddG/ddE data into the 5 regions  
for (tile in c('1','2','3','4','5')) {
  print(c('tile', tile))
  
  ddG_part <- ddG[ddG$aa_pos >= h[[tile]][['start_pos']] &  ddG$aa_pos <= h[[tile]][['tile_end']], ]
  
  p <- ggplot(ddG_part, aes(Rosetta_ddg_score)) +
    geom_density(bw= 0.08) +
    xlim(-3,20) +
    ggtitle(paste0("tile ", tile, ", ddG values")) +
    xlab('Bandwidth = 0.08') +
    theme_bw(base_size = 20)
  
  png(filename = paste0(plot_dir,'tile',tile,'.ddG.png'),width = 800, height = 550)
  print(p)
  dev.off()
  
  ddE_part <- ddE[ddE$aa_pos >= h[[tile]][['start_pos']] &  ddE$aa_pos <= h[[tile]][['tile_end']], ]
  
  p <- ggplot(ddE_part, aes(gemme_score)) +
    geom_density(bw= 0.08) +
    ggtitle(paste0("tile ", tile, ", ddE values")) +
    xlab('Bandwidth = 0.08') +
    theme_bw(base_size = 20)
  
  png(filename = paste0(plot_dir,'tile',tile,'.ddE.png'),width = 800, height = 550)
  print(p)
  dev.off()

}

#for all tiles together: ddG
#add tile column. The last else is 5 since the aa_pos wasn't in any of the other regions.
ddG$tile <- ifelse(ddG$aa_pos >= h[['1']][['start_pos']] &  ddG$aa_pos <= h[['1']][['tile_end']], 'tile1', 
                   ifelse(ddG$aa_pos >= h[['2']][['start_pos']] &  ddG$aa_pos <= h[['2']][['tile_end']], 'tile2', 
                          ifelse(ddG$aa_pos >= h[['3']][['start_pos']] &  ddG$aa_pos <= h[['3']][['tile_end']], 'tile3',
                                 ifelse(ddG$aa_pos >= h[['4']][['start_pos']] &  ddG$aa_pos <= h[['4']][['tile_end']], 'tile4', 'tile5'))))

ddG <- ddG[!is.na(ddG$tile),]

p <- ggplot(ddG, aes(Rosetta_ddg_score,color=as.factor(tile))) +
  geom_density(bw= 0.08) +
  xlim(-3,20) +
  ggtitle(paste0("ddG scores per tile")) +
  xlab('Bandwidth = 0.08') +
  theme_bw(base_size = 20) + theme(legend.position="bottom")

p
png(filename = paste0(plot_dir,'alltiles.ddG.png'),width = 800, height = 600)
print(p)
dev.off()

p2 <- ggplot(ddG, aes(Rosetta_ddg_score,fill=as.factor(tile))) +
  geom_density(bw= 0.08, alpha = 0.5) +
  xlim(-3,20) +
  ggtitle(paste0("ddG scores per tile")) +
  xlab('Bandwidth = 0.08') +
  theme_bw(base_size = 20) + theme(legend.position="bottom")

p2
png(filename = paste0(plot_dir,'alltiles.fill.ddG.png'),width = 800, height = 600)
print(p2)
dev.off()

#overlapping for all tiles: ddE
#add tile column. The last else is 5 since the aa_pos wasn't in any of the other regions.
ddE$tile <- ifelse(ddE$aa_pos >= h[['1']][['start_pos']] &  ddE$aa_pos <= h[['1']][['tile_end']], 'tile1', 
                   ifelse(ddE$aa_pos >= h[['2']][['start_pos']] &  ddE$aa_pos <= h[['2']][['tile_end']], 'tile2', 
                          ifelse(ddE$aa_pos >= h[['3']][['start_pos']] &  ddE$aa_pos <= h[['3']][['tile_end']], 'tile3',
                                 ifelse(ddE$aa_pos >= h[['4']][['start_pos']] &  ddE$aa_pos <= h[['4']][['tile_end']], 'tile4', 'tile5'))))
  
ddE <- ddE[!is.na(ddE$tile),]

p <- ggplot(ddE, aes(gemme_score,color=as.factor(tile))) +
  geom_density(bw= 0.08) +
  ggtitle(paste0("ddE scores per tile")) +
  xlab('Bandwidth = 0.08') +
  theme_bw(base_size = 20) + theme(legend.position="bottom")

p
png(filename = paste0(plot_dir,'alltiles.ddE.png'),width = 800, height = 600)
print(p)
dev.off()

p2 <- ggplot(ddE, aes(gemme_score,fill=as.factor(tile))) +
  geom_density(bw= 0.08, alpha = 0.5) +
  ggtitle(paste0("ddE scores per tile")) +
  xlab('Bandwidth = 0.08') +
  theme_bw(base_size = 20) + theme(legend.position="bottom")

p2
png(filename = paste0(plot_dir,'alltiles.fill.ddE.png'),width = 800, height = 600)
print(p2)
dev.off()

#number of scored vars with different enrich settings in the new data: 
#minimum 10 read counts to accept a DNA level var as real VS accepting all DNA level vars that pass quality control
##########################

tile <- '4'
for (tile in c('1','2','3','4','5')) {
  print(c('tile', tile))

  dat10 <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps_minvarcount10/tsv/tile',tile ,'_all_reps_exp/main_synonymous_scores.tsv'), 
                  sep = '\t', header = F,
                  stringsAsFactors = F, skip = 4,
                  col.names = c('var', 'SE', 'epsilon', 'score'))
  
  dat <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps/tsv/tile',tile ,'_all_reps_exp/main_synonymous_scores.tsv'), 
                  sep = '\t', header = F,
                  stringsAsFactors = F, skip = 4,
                  col.names = c('var', 'SE', 'epsilon', 'score'))
  
  n_single_10 <- nrow(dat10[!grepl(',',dat10$var, fixed = T),])
  n_single_all <- nrow(dat[!grepl(',',dat$var, fixed = T),])
  print(paste0('No. scored vars with only counting DNA vars counts >= 10: ',n_single_10))
  print(paste0('No. scored vars with counting all DNA vars passing qual check: ',n_single_all))

}

#debug
dat_10_single <- dat10[!grepl(',',dat10$var, fixed = T),]
dat_single <- dat[!grepl(',',dat$var, fixed = T),]
sum(grepl('Ter',dat_single$var, fixed = T))

#Plot MAVE score distributions
###############################


#plot distributions per nr mutations
###############################

#tile = '5' #for debugging
#base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/results/enrich_full_dhfr/'

for (tile in c('1','2','3','4','5')) {
  
  print(c('tile', tile))

  dat10 <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps_minvarcount10',corr,'/tsv/tile',tile ,
                           '_all_reps_exp/main_synonymous_scores.tsv'), 
                  sep = '\t', header = F,
                  stringsAsFactors = F, skip = 2,
                  col.names = c('var', 'SE', 'epsilon', 'score'))
  
  dat10$nr_subs <- str_count(dat10$var,",") + 1
  
  wt_score <- dat10[dat10$var == '_wt',c('score')]
  syn_score <- dat10[dat10$var == '_sy',c('score')]
  p <- ggplot(dat10, aes(score,color=as.factor(nr_subs))) +
    geom_density(bw= 0.08) +
    xlim(-3.5,3.5) + 
    #use aes to get a legend for the vertical line:
    geom_vline(aes(xintercept = wt_score, colour='WT score')) +
    geom_vline(aes(xintercept = syn_score, colour='Synonymous score'), linetype = 'dashed') +
    #according to https://stackoverflow.com/questions/39112735/using-colors-in-aes-function-in-ggplot2 one can use scale_color_manual to choose colors
    #need to specify a color for each aes mapping. 1-4 are the number of mutations, 'WT score' and "Synonymous score" have been named so by the two geom_vline calls above
    scale_colour_manual(values = c('1'='coral2', '2'= 'steelblue1', '3'='darkolivegreen3', 
                                   '4'='darkorchid2',"WT score" = "black", "Synonymous score" = "grey50")) +
  
    ggtitle(paste0("tile ", tile, " Reseq MAVE scores, min DNA var count 10,\n", corr_descr, 
                   "\n1 sub: ", table(dat10$nr_subs)[1], 
                   ' vars, 2 subs: ', table(dat10$nr_subs)[2], 
                   " vars, 3 subs: ", table(dat10$nr_subs)[3], ' vars, 4 subs: ', table(dat10$nr_subs)[4], " vars.")) +
    xlab('Bandwidth = 0.08') +
    theme_bw(base_size = 18)
  
  p #look at plot for debugging
  
  png(filename = paste0(plot_dir,'tile',tile,'.perNrMuts.minvarcount10',corr,'.png'),width = 800, height = 550)
  print(p)
  dev.off()
  
}

#plot distribution of all MAVE scores in pilot run (only seq'ing tile 1)
####################################
data_dir = "/home/henrike/Documents/PD_AS/projects/Sofie_Mave/"
ext = "results/enrich/latest_re-run/"

d1 <- read.csv(paste0(data_dir, ext, 'ex123_all_wo_20_22_38/tsv/ex123_exp/main_synonymous_scores.tsv'), 
               sep = '\t', header = F, 
               stringsAsFactors = F, skip = 4, 
               col.names = c('var','37_se',	'37_eps',	'37_score',	'37_Cu_se',	'37_Cu_eps',	'37_Cu_score',	
                             'MTX_se_old',	'MTX_eps_old',	'MTX_score_old',	
                             'MTX_Cu_se_old',	'MTX_Cu_eps_old',	'MTX_Cu_score_old'))

#select cols for MTX
MTX_scores_r <- d1[,c('var', 'MTX_se_old',	'MTX_eps_old',	'MTX_score_old')]
#drop empty rows: remaining 705
MTX_scores <- MTX_scores_r[complete.cases(MTX_scores_r),]
#compare with results of !is.na: 705
sum(!is.na(MTX_scores_r$MTX_score_old))
#only single variants
MTX_single <- MTX_scores[!grepl(',',MTX_scores$var, fixed = T),]

p <- ggplot(MTX_single, aes(x=MTX_score_old)) + 
  geom_density(bw= 0.08) +
  xlim(-3,3) + 
  ggtitle(paste0("MXT: single aa substitution MAVE scores.\nFirst sequencing. # vars: ",nrow(MTX_single))) +
  xlab('Bandwidth = 0.08') +
  theme_bw(base_size = 20)

png(filename = paste0(plot_dir,'pilot_seq_MTX_single_vars.png'),width = 800, height = 550)
p
dev.off()

#can also use density() function and then call plot o nthe result, but ggplot is more easily customizable
#MTX_single.density <- density(MTX_single$MTX_score_old, bw = 0.08)
#plot(MTX_single.density, main = paste0("MXT: single aa substitution MAVE scores. First sequencing. # vars: ",nrow(MTX_single)))

