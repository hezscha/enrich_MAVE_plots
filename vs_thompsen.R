library(ggplot2)
library(ggrepel)
library(stringr)
library(reshape2)
library(dplyr) #arrange
library(hash)

#Comparison between our MAVE on human dhfr and data from this study: DOI: https://doi.org/10.7554/eLife.53476 Altered expression of a quality control protease in E. coli reshapes the in vivo mutational landscape of a model enzyme, Thompsen et al 2019

base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/results/enrich_full_dhfr/'
data_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/data/'
plot_dir <- paste0(data_dir, '/HZ_plots/vs_thompsen/')

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

start_pos <- 2
tile_end <- 187


#load in the combined df I prepared in ipython notebook. Positions that are aligned in the struct alignment are in the same row
#########################################
dat <- read.csv(paste0(data_dir,'merge_human_ecoli_fulldhfr.csv'), sep = '\t')
dat$pos <- as.numeric(str_extract(dat$human_variant, "([0-9]+)"))

dump_df <- dat[,c('human_variant', 'score', 'SE', 'tile', 'pos', 'ecoli_variant', 'Lon_minus_score', 'Lon_minus_std_dev',
                  'Lon_minus_std_err', 'Lon_plus_score', 'Lon_plus_std_dev', 'Lon_plus_std_err')]

#write.csv(dump_df, file = paste0(base_dir, 'combined_R_dataframes/', 'human_VS_ecoli_DHFR_data.csv'), quote = F, row.names = F)


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

#load merged ddG/E data
###########################

dds <- read.csv(paste0(data_dir,'prism_merge_ddE_ddG_DHFR_P00374.txt'), comment.char = '#', sep = ' ')
#, col.names = c('var', 'gemme_score', 'ddG'))
#dds$pos <- as.numeric(str_extract(dds$variant, "(?<=[A-Z])([0-9]+)(?=[A-Z])"))
#get the (single letter) alt aa out of the variant name with str_match, then translate every single letter aa to the corresponding 3 letter code
#dds$sub <- sapply(str_match(dds$variant, "[0-9]+([A-Za-z])")[,2], function(i) aa_1_to_3[[i]])

#more explicit version: 
#dds$sub1 <- str_match(dds$variant, "[0-9]+([A-Za-z])")[,2]
#b <- sapply(dds$sub1, function(i) aa_1_to_3[[i]])

# dds$tile <- ifelse(dds$pos >= h[['1']][['start_pos']] &  dds$pos <= h[['1']][['tile_end']], 'tile1', 
#                    ifelse(dds$pos >= h[['2']][['start_pos']] &  dds$pos <= h[['2']][['tile_end']], 'tile2', 
#                           ifelse(dds$pos >= h[['3']][['start_pos']] &  dds$pos <= h[['3']][['tile_end']], 'tile3',
#                                  ifelse(dds$pos >= h[['4']][['start_pos']] &  dds$pos <= h[['4']][['tile_end']], 'tile4',
#                                         ifelse(dds$pos >= h[['5']][['start_pos']] &  dds$pos <= h[['5']][['tile_end']], 'tile5', NA)))))

#there are some rows at the end with ddGs for the non-vars (G16=). Remove those 
#for reasons I don't understand complete.cases doesn't want to operate on dds even though it is a dataframe, so I called it on the gemme score column instead
dds <- dds[complete.cases(dds$gemme_score_00),]

#add distances and ddE/G to the human/ecoli aligned MAVE scores
############################

all_merge <- merge(dump_df, mtx_dist[,c('MTXdist','pos')], by = 'pos', all.x = T)
all_merge <- merge(all_merge, nadph_dist[,c('NAPDHdist', 'pos')], by = 'pos', all.x = T)
all_merge <- merge(all_merge, dds, by.x = 'human_variant', by.y = 'variant', all.x = T)

all_merge$dist_group_MTX <- ifelse(all_merge$MTXdist < 10,'MTX_close','MTX_far')
all_merge$dist_group_NADPH <- ifelse(all_merge$NAPDHdist < 10,'NADPH_close','NADPH_far')

# bla <- merge(dump_df[,c('human_variant', 'score', 'SE', 'tile','ecoli_variant', 'Lon_minus_score', 'Lon_minus_std_dev',
#                         'Lon_minus_std_err', 'Lon_plus_score', 'Lon_plus_std_dev', 'Lon_plus_std_err')], 
#              dds, by.x = 'human_variant', by.y = 'variant', all = T)


#write.csv(all_merge, file = paste0(base_dir, 'combined_R_dataframes/', 'human_VS_ecoli_DHFR_data_dds_dists.csv'), quote = F, row.names = F)


#this was wrong! because the dds$long_sub referred to the pos in human DHFR and dump_df$ecoli_sub to the position in ecoli DHFR, so they didn't match up
# #combine human/ecoli aligned scores with ddE
# ##############################
# 
# dump_df$ecoli_sub <- str_match(dump_df$ecoli_variant, "([0-9]+[A-Za-z]+)")[,2]
# dds$long_sub <- str_match(dds$variant, "([0-9]+[A-Za-z]+)")[,2]
# 
# comb_dde <- merge(dump_df, dds, by.x = 'ecoli_sub', by.y = 'long_sub')


# #combine everything
# ########################################
# 
# m$ecoli_sub <- str_match(m$ecoli_variant, "([0-9]+[A-Za-z]+)")[,2]
# dds$long_sub <- str_match(dds$variant, "([0-9]+[A-Za-z]+)")[,2]
# 
# m_dds <-  merge(m, dds, by.x = 'ecoli_sub', by.y = 'long_sub')


#plot ddE/G VS thmopsen scores
####################################

min_score <- -2.5
max_score <- 4.0

#Rosetta_ddg_score_01
p <- ggplot(all_merge, aes(y=gemme_score_00, x=Lon_minus_score, color = score)) +
  geom_point() + 
  scale_color_gradient2(name = "Human DHFR\nMAVE Score", midpoint = 0, low = "#830823", high = "#063264", mid = 'white', 
                        limits = c(min_score,max_score)) +
  geom_text_repel(aes(label=ifelse((Lon_minus_score < -3 & !is.na(Lon_minus_score)), as.character(ecoli_variant),'') ),
                  color = 'black') +  
  theme(text = element_text(size = 20))
  #theme_bw(base_size = 20)

p

png(filename = paste0(plot_dir,'lon_minus_VS_human_ddE_color_humanMAVE.png'), 
    width = 2000, height = 1600, res = 200)
p
dev.off()

p <- ggplot(all_merge, aes(y=gemme_score_00, x=Lon_plus_score, color = score)) +
  geom_point() + 
  scale_color_gradient2(name = "Human DHFR\nMAVE Score", midpoint = 0, low = "#830823", high = "#063264", mid = 'white', 
                        limits = c(min_score,max_score)) +
  geom_text_repel(aes(label=ifelse((Lon_plus_score < -3 & !is.na(Lon_plus_score)), as.character(ecoli_variant),'') ),
                  color = 'black') +  
  theme(text = element_text(size = 20))
#theme_bw(base_size = 20)

p

png(filename = paste0(plot_dir,'lon_plus_VS_human_ddE_color_humanMAVE.png'), 
    width = 2000, height = 1600, res = 200)
p
dev.off()

#plot our MAVE scores VS Thompsen scores
###########################

min_score = min(all_merge$Lon_minus_score, all_merge$Lon_plus_score, all_merge$score, na.rm = T) - 
  max(abs(all_merge$Lon_minus_std_err), abs(all_merge$Lon_plus_std_err), abs(all_merge$SE),na.rm = T)
max_score = max(all_merge$Lon_minus_score, all_merge$Lon_plus_score, all_merge$score, na.rm = T) + 
  max(abs(all_merge$Lon_minus_std_err), abs(all_merge$Lon_plus_std_err), abs(all_merge$SE),na.rm = T)

#minus lon

#color by distance to MTX in human DHFR
p <- ggplot(all_merge, aes(x=score,y=Lon_minus_score, color = MTXdist)) +
  geom_point() +
  scale_color_gradient2(midpoint = 10, low = "#830823", high = "#063264", mid = 'white') +
  xlim(min_score, max_score) + ylim(min_score,max_score) +
  theme_bw(base_size = 20)

p
png(filename = paste0(plot_dir,'lon_minus_color_MTXdist.png'), 
    width = 800, height = 600)
p
dev.off()

#color by distance groupings
p <- ggplot(all_merge, aes(x=score,y=Lon_minus_score, color = interaction(dist_group_MTX, dist_group_NADPH))) +
  geom_point() +
  scale_color_discrete(name = 'dist', labels = c('both close', 'NADPH close', 'MTX close', 'both far')) +
  #scale_color_gradient2(midpoint = 10, low = "#830823", high = "#063264", mid = 'white') +
  #xlim(min_score, max_score) + ylim(min_score,max_score) +
  #geom_errorbar(aes(ymin=Lon_minus_score-Lon_minus_std_err, ymax=Lon_minus_score+Lon_minus_std_err), width=.2, linetype = 1, color = 'grey80', alpha = 0.5) +
  #geom_errorbar(aes(xmin=score-SE, xmax=score+SE), width=.2, linetype = 1, color = 'grey80', alpha = 0.5) +
  theme_bw(base_size = 20) + theme(legend.position = 'bottom')

p
png(filename = paste0(plot_dir,'lon_minus_color_dist_group_non_square.png'), 
    width = 2000, height = 1600, res = 200)
p
dev.off()

#plus lon
p2 <- ggplot(all_merge, aes(x=score,y=Lon_plus_score, color = tile)) +
  geom_point() +
  xlim(min_score, max_score) + ylim(min_score,max_score) +
  theme_bw(base_size = 20)
p2
png(filename = paste0(plot_dir,'lon_plus.png'), 
    width = 800, height = 600)
p2
dev.off()

#color by distance groupings
p <- ggplot(all_merge, aes(x=score,y=Lon_plus_score, color = interaction(dist_group_MTX, dist_group_NADPH))) +
  geom_point() +
  scale_color_discrete(name = 'dist', labels = c('both close', 'NADPH close', 'MTX close', 'both far')) +
  #scale_color_gradient2(midpoint = 10, low = "#830823", high = "#063264", mid = 'white') +
  #xlim(min_score, max_score) + ylim(min_score,max_score) +
  #geom_errorbar(aes(ymin=Lon_plus_score-Lon_plus_std_err, ymax=Lon_plus_score+Lon_plus_std_err), width=.2, linetype = 1, color = 'grey80', alpha = 0.5) +
  #geom_errorbar(aes(xmin=score-SE, xmax=score+SE), width=.2, linetype = 1, color = 'grey80', alpha = 0.5) +
  geom_text_repel(aes(label=ifelse((Lon_plus_score < -2.9 & !is.na(Lon_plus_score)), as.character(ecoli_variant),'') )) +
  theme_bw(base_size = 20) + theme(legend.position = 'bottom')

p
png(filename = paste0(plot_dir,'lon_plus_color_dist_group_non_square.png'), 
    width = 2000, height = 1600, res = 200)
p
dev.off()

#take apart the distance groups and plot them separately to see if we see something
################################

m_both_far <- all_merge[all_merge$dist_group_MTX == 'MTX_far' & all_merge$dist_group_NADPH == 'NADPH_far',]
m_both_close <- all_merge[all_merge$dist_group_MTX == 'MTX_close' & all_merge$dist_group_NADPH == 'NADPH_close',]
m_MTX_close <- all_merge[all_merge$dist_group_MTX == 'MTX_close' & all_merge$dist_group_NADPH == 'NADPH_far',]
m_NADPH_close <- all_merge[all_merge$dist_group_MTX == 'MTX_far' & all_merge$dist_group_NADPH == 'NADPH_close',]

#lon minus
p <- ggplot(m_both_far, aes(x=score,y=Lon_minus_score)) +
  geom_point() +
  ggtitle('Variant positions >10 A from both ligands') +
  #xlim(min_score, max_score) + ylim(min_score,max_score) +
  geom_errorbar(aes(ymin=Lon_minus_score-Lon_minus_std_err, ymax=Lon_minus_score+Lon_minus_std_err), width=.2, linetype = 1, color = 'grey80', alpha = 0.7) +
  geom_errorbar(aes(xmin=score-SE, xmax=score+SE), width=.2, linetype = 1, color = 'grey80', alpha = 0.7) +
  theme_bw(base_size = 16)

p
png(filename = paste0(plot_dir,'lon_minus_VS_MAVE_only_far_group_non_square.png'),
    width = 2000, height = 1600, res = 200)
p
dev.off()

p <- ggplot(m_both_close, aes(x=score,y=Lon_minus_score, color = gemme_score_00)) +
  geom_point() +
  ggtitle('Variant positions <=10 A from both ligands') + 
  #xlim(min_score, max_score) + ylim(min_score,max_score) +
  geom_text_repel(aes(label=ifelse((Lon_minus_score < -3 & !is.na(Lon_minus_score)), as.character(ecoli_variant),'') )) +
  geom_errorbar(aes(ymin=Lon_minus_score-Lon_minus_std_err, ymax=Lon_minus_score+Lon_minus_std_err), width=.2, linetype = 1, color = 'grey80', alpha = 0.7) +
  geom_errorbar(aes(xmin=score-SE, xmax=score+SE), width=.2, linetype = 1, color = 'grey80', alpha = 0.7) +
  theme_bw(base_size = 16) + theme(legend.position = 'bottom')

p
png(filename = paste0(plot_dir,'lon_minus_VS_MAVE_only_close_group_color_ddE_non_square.png'), 
    width = 2000, height = 1600, res = 200)
p
dev.off()

p <- ggplot(m_both_close, aes(x=score,y=Lon_minus_score, color = Rosetta_ddg_score_01)) +
  geom_point() +
  ggtitle('Variant positions <=10 A from both ligands') + 
  scale_color_gradient2(mid = 'darkorange', low = 'forestgreen', high = 'red', midpoint = 7, name = 'Rosetta ddG') +
  #scale_color_gradient(low = 'forestgreen', high = 'red',name = 'Rosetta ddG') +
  #xlim(min_score, max_score) + ylim(min_score,max_score) +
  geom_text_repel(aes(label=ifelse((Lon_minus_score < -3 & !is.na(Lon_minus_score)), as.character(ecoli_variant),'') )) +
  geom_errorbar(aes(ymin=Lon_minus_score-Lon_minus_std_err, ymax=Lon_minus_score+Lon_minus_std_err), width=.2, linetype = 1, color = 'grey80', alpha = 0.7) +
  geom_errorbar(aes(xmin=score-SE, xmax=score+SE), width=.2, linetype = 1, color = 'grey80', alpha = 0.7) +
  
  theme_bw(base_size = 16) + theme(legend.position = 'bottom')

p
png(filename = paste0(plot_dir,'lon_minus_VS_MAVE_only_close_group_color_ddG_non_square.png'),
    width = 2000, height = 1600, res = 200)
p
dev.off()

p <- ggplot(m_MTX_close, aes(x=score,y=Lon_minus_score)) +
  geom_point() +
  ggtitle('Variant positions <=10 A from MTX') + 
  #xlim(min_score, max_score) + ylim(min_score,max_score) +
  geom_errorbar(aes(ymin=Lon_minus_score-Lon_minus_std_err, ymax=Lon_minus_score+Lon_minus_std_err), width=.2, linetype = 1, color = 'grey80', alpha = 0.7) +
  geom_errorbar(aes(xmin=score-SE, xmax=score+SE), width=.2, linetype = 1, color = 'grey80', alpha = 0.7) +
  theme_bw(base_size = 16) 

p
png(filename = paste0(plot_dir,'lon_minus_VS_MAVE_only_MTX_close_non_square.png'),
    width = 2000, height = 1600, res = 200)
p
dev.off()

p <- ggplot(m_NADPH_close, aes(x=score,y=Lon_minus_score)) +
  geom_point() +
  ggtitle('Variant positions <=10 A from NADPH') + 
  #xlim(min_score, max_score) + ylim(min_score,max_score) +
  geom_errorbar(aes(ymin=Lon_minus_score-Lon_minus_std_err, ymax=Lon_minus_score+Lon_minus_std_err), width=.2, linetype = 1, color = 'grey80', alpha = 0.7) +
  geom_errorbar(aes(xmin=score-SE, xmax=score+SE), width=.2, linetype = 1, color = 'grey80', alpha = 0.7) +
  theme_bw(base_size = 16) 

p
png(filename = paste0(plot_dir,'lon_minus_VS_MAVE_only_NADPH_close_non_square.png'),
    width = 2000, height = 1600, res = 200)
p
dev.off()
 
#lon plus
################################
p <- ggplot(m_both_far, aes(x=score,y=Lon_minus_score)) +
  geom_point() +
  ggtitle('Variant positions >10 A from both ligands') +
  #xlim(min_score, max_score) + ylim(min_score,max_score) +
  geom_errorbar(aes(ymin=Lon_plus_score-Lon_plus_std_err, ymax=Lon_plus_score+Lon_plus_std_err), width=.2, linetype = 1, color = 'grey80', alpha = 0.7) +
  geom_errorbar(aes(xmin=score-SE, xmax=score+SE), width=.2, linetype = 1, color = 'grey80', alpha = 0.7) +
  theme_bw(base_size = 16)

p
png(filename = paste0(plot_dir,'Lon_plus_VS_MAVE_only_far_group_non_square.png'),
    width = 2000, height = 1600, res = 200)
p
dev.off()

p <- ggplot(m_both_close, aes(x=score,y=Lon_plus_score, color = gemme_score_00)) +
  geom_point() +
  ggtitle('Variant positions <=10 A from both ligands') + 
  #xlim(min_score, max_score) + ylim(min_score,max_score) +
  geom_text_repel(aes(label=ifelse((Lon_plus_score < -3 & !is.na(Lon_plus_score)), as.character(ecoli_variant),'') )) +
  geom_errorbar(aes(ymin=Lon_plus_score-Lon_plus_std_err, ymax=Lon_plus_score+Lon_plus_std_err), width=.2, linetype = 1, color = 'grey80', alpha = 0.7) +
  geom_errorbar(aes(xmin=score-SE, xmax=score+SE), width=.2, linetype = 1, color = 'grey80', alpha = 0.7) +
  theme_bw(base_size = 16) + theme(legend.position = 'bottom')

p
png(filename = paste0(plot_dir,'Lon_plus_VS_MAVE_only_close_group_color_ddE_non_square.png'), 
    width = 2000, height = 1600, res = 200)
p
dev.off()

p <- ggplot(m_both_close, aes(x=score,y=Lon_plus_score, color = Rosetta_ddg_score_01)) +
  geom_point() +
  ggtitle('Variant positions <=10 A from both ligands') + 
  scale_color_gradient2(mid = 'darkorange', low = 'forestgreen', high = 'red', midpoint = 7, name = 'Rosetta ddG') +
  #scale_color_gradient(low = 'forestgreen', high = 'red',name = 'Rosetta ddG') +
  #xlim(min_score, max_score) + ylim(min_score,max_score) +
  geom_text_repel(aes(label=ifelse((Lon_plus_score < -3 & !is.na(Lon_plus_score)), as.character(ecoli_variant),'') )) +
  geom_errorbar(aes(ymin=Lon_plus_score-Lon_plus_std_err, ymax=Lon_plus_score+Lon_plus_std_err), width=.2, linetype = 1, color = 'grey80', alpha = 0.7) +
  geom_errorbar(aes(xmin=score-SE, xmax=score+SE), width=.2, linetype = 1, color = 'grey80', alpha = 0.7) +
  
  theme_bw(base_size = 16) + theme(legend.position = 'bottom')

p
png(filename = paste0(plot_dir,'Lon_plus_VS_MAVE_only_close_group_color_ddG_non_square.png'),
    width = 2000, height = 1600, res = 200)
p
dev.off()

p <- ggplot(m_MTX_close, aes(x=score,y=Lon_plus_score)) +
  geom_point() +
  ggtitle('Variant positions <=10 A from MTX') + 
  #xlim(min_score, max_score) + ylim(min_score,max_score) +
  geom_errorbar(aes(ymin=Lon_plus_score-Lon_plus_std_err, ymax=Lon_plus_score+Lon_plus_std_err), width=.2, linetype = 1, color = 'grey80', alpha = 0.7) +
  geom_errorbar(aes(xmin=score-SE, xmax=score+SE), width=.2, linetype = 1, color = 'grey80', alpha = 0.7) +
  theme_bw(base_size = 16) 

p
png(filename = paste0(plot_dir,'Lon_plus_VS_MAVE_only_MTX_close_non_square.png'),
    width = 2000, height = 1600, res = 200)
p
dev.off()

p <- ggplot(m_NADPH_close, aes(x=score,y=Lon_plus_score)) +
  geom_point() +
  ggtitle('Variant positions <=10 A from NADPH') + 
  #xlim(min_score, max_score) + ylim(min_score,max_score) +
  geom_errorbar(aes(ymin=Lon_plus_score-Lon_plus_std_err, ymax=Lon_plus_score+Lon_plus_std_err), width=.2, linetype = 1, color = 'grey80', alpha = 0.7) +
  geom_errorbar(aes(xmin=score-SE, xmax=score+SE), width=.2, linetype = 1, color = 'grey80', alpha = 0.7) +
  theme_bw(base_size = 16) 

p
png(filename = paste0(plot_dir,'Lon_plus_VS_MAVE_only_NADPH_close_non_square.png'),
    width = 2000, height = 1600, res = 200)
p
dev.off()
 
