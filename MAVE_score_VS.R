library(ggplot2)
library(stringr)
library(reshape2)
library(matrixStats)
library(ggrepel)
library(dplyr) #arrange
library(hash)
#add marginal distr to dde vs ddg plot
library(gridExtra)
library(cowplot)
library(egg)

#plot average or indiv MAVE score by position against various things including ddsp accessible surface area

base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/results/enrich_full_dhfr/'
data_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/data/'
plot_dir <- paste0(data_dir, '/HZ_plots/vs_plots/')

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

#load RSA data
###########################

rsa <- read.csv(paste0(data_dir,'prism_dssp_DHFR_P00374_1u72_parsed.txt'), 
                sep = ' ', header = T, comment.char = '#')
rsa$pos <- seq(start_pos,tile_end)

#load in distance to NAPDH and MTX
############################

mtx_dist <- read.csv(paste0(data_dir,'1u72_dist_to_MTX.dat'), 
                       sep = '\t', header = F,
                     col.names = c('pos_raw','chain', 'WTres', 'MTXdist'))

#put in correct position starting from 2
mtx_dist$pos <- mtx_dist$pos_raw+1

#remove distances of waters (empty residue field)
#mtx_dist <- mtx_dist[complete.cases(mtx_dist),]
mtx_dist <- mtx_dist[mtx_dist$WTres != '',]
mtx_dist$molecule <- rep('MTX', nrow(mtx_dist))

nadph_dist <- read.csv(paste0(data_dir,'1u72_dist_to_NDP.dat'), 
                     sep = '\t', header = F,
                     col.names = c('pos_raw','chain', 'WTres', 'NADPHdist'))

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
dds$pos <- as.numeric(str_extract(dds$variant, "(?<=[A-Z])([0-9]+)(?=[A-Z])"))
#get the (single letter) alt aa out of the variant name with str_match, then translate every single letter aa to the corresponding 3 letter code
dds$sub <- sapply(str_match(dds$variant, "[0-9]+([A-Za-z])")[,2], function(i) aa_1_to_3[[i]])

#more explicit version: 
#dds$sub1 <- str_match(dds$variant, "[0-9]+([A-Za-z])")[,2]
#b <- sapply(dds$sub1, function(i) aa_1_to_3[[i]])

# dds$tile <- ifelse(dds$pos >= h[['1']][['start_pos']] &  dds$pos <= h[['1']][['tile_end']], 'tile1', 
#                    ifelse(dds$pos >= h[['2']][['start_pos']] &  dds$pos <= h[['2']][['tile_end']], 'tile2', 
#                           ifelse(dds$pos >= h[['3']][['start_pos']] &  dds$pos <= h[['3']][['tile_end']], 'tile3',
#                                  ifelse(dds$pos >= h[['4']][['start_pos']] &  dds$pos <= h[['4']][['tile_end']], 'tile4',
#                                         ifelse(dds$pos >= h[['5']][['start_pos']] &  dds$pos <= h[['5']][['tile_end']], 'tile5', NA)))))



############################

#load in prepped MAVE data
############################
all_df <- read.csv(paste0(base_dir, 'combined_R_dataframes/', 'all_tile_scores_long.csv'), 
                          sep = ',')
mat <- read.csv(paste0(base_dir, 'combined_R_dataframes/', 'per_pos_matrix.csv'), 
                sep = ',')

#merge everything
############################
#long table: one line per var
comb <- merge(dds, all_df, by = c('pos', 'sub'))
comb <- merge(comb, mtx_dist[,c('MTXdist', 'pos')], by = 'pos')
comb <- merge(comb, nadph_dist[,c('NADPHdist', 'pos')], by = 'pos')
comb <- merge(comb, rsa[,c('ASA','pos')], by = 'pos')
comb$dist_group <- ifelse((comb$MTXdist <= 10 & comb$NADPHdist <= 10), 'both_close', 
                          ifelse((comb$MTXdist <= 10 & comb$NADPHdist > 10), 'MTX_close', 
                          ifelse((comb$MTXdist > 10 & comb$NADPHdist <= 10), 'NADPH_close', 'both_far')))
  
table(comb$dist_group)

#matrix: one line per position, incl all var scores plus median pos
dist_mave <- merge(mat[,c('pos_median','pos')], mtx_dist[,c('pos','MTXdist')], by = 'pos')
dist_mave <- merge(dist_mave, nadph_dist[,c('pos','NADPHdist')], by = 'pos')
dist_mave <- merge(dist_mave, rsa[,c('ASA','pos')], by = 'pos')

dist_mave$tile <- ifelse(dist_mave$pos >= h[['1']][['start_pos']] &  dist_mave$pos <= h[['1']][['tile_end']], 'tile1',
                   ifelse(dist_mave$pos >= h[['2']][['start_pos']] &  dist_mave$pos <= h[['2']][['tile_end']], 'tile2',
                          ifelse(dist_mave$pos >= h[['3']][['start_pos']] &  dist_mave$pos <= h[['3']][['tile_end']], 'tile3',
                                 ifelse(dist_mave$pos >= h[['4']][['start_pos']] &  dist_mave$pos <= h[['4']][['tile_end']], 'tile4',
                                        ifelse(dist_mave$pos >= h[['5']][['start_pos']] &  dist_mave$pos <= h[['5']][['tile_end']], 'tile5', NA)))))


range(mat$pos_median, na.rm = T)

######################################
#plots


#avrg mave score per position VS distance to NADPH and MTX
########################################

p <- ggplot(dist_mave, aes(x=pos_median, y=MTXdist)) +
  geom_point() +
  ggtitle('Distance to the substrate analog MTX') +
  xlab('Median MAVE score') + ylab('Distance to MTX') +
  geom_text_repel(aes(label=ifelse(abs(pos_median)>=1, as.character(pos),''))) +
  theme_bw(base_size = 20)
p
png(filename = paste0(plot_dir,'median_mave_VS_dist_MTX.png'),width = 800, height = 550)
print(p)
dev.off()

p <- ggplot(dist_mave, aes(x=pos_median, y=NADPHdist)) +
  geom_point() +
  ggtitle('Distance to the ligand NADPH') +
  xlab('Median MAVE score') + ylab('Distance to NADPH') +
  theme_bw(base_size = 20)

png(filename = paste0(plot_dir,'median_mave_VS_dist_NADPH.png'),width = 800, height = 550)
print(p)
dev.off()

both_dist <- melt(dist_mave, id.vars = c('pos','pos_median'), variable.name = 'ligand', value.name = 'dist')
#both_dist <- rbind(dist_mave, dist_mave_nadph)

p2 <- ggplot(both_dist, aes(x=pos_median, y=dist, color=ligand)) +
  geom_point() +
  ggtitle('Distance to both ligands') +
  xlab('Median MAVE score') + ylab('Distance') +
  theme_bw(base_size = 20)
p2
png(filename = paste0(plot_dir,'median_mave_VS_dist_both.png'),width = 800, height = 550)
print(p2)
dev.off()


#RSA VS MAVE
##########################################

#add median mave score per position to rsa for plotting
# rsa$mave_score <- mat$pos_median
# rsa$pos <- seq(start_pos,tile_end)
# rsa$tile <- mat$tile

p1 <- ggplot(dist_mave, aes(x=pos_median, y=ASA, color = tile)) +
  geom_point(alpha = 0.7) +
  ggtitle('RSA VS MAVE scores, all tiles') +
  xlab('Median MAVE score') + ylab('Relative solvent accessible area (RSA)') +
  #geom_text_repel(aes(label=pos)) +
  theme_bw(base_size = 20)
p1
png(filename = paste0(plot_dir,'median_mave_VS_RSA_per_tile.png'),width = 800, height = 550)
print(p1)
dev.off()
##########################################





#ddg VS ddE
############################

#mave score scale limits from heatmaps
min_score <- -2.5
max_score <- 4.0

#introduce max ddG cutoff
cutoff_ddG <- 15
comb$Rosetta_ddg_score_01 <- ifelse(comb$Rosetta_ddg_score_01 > cutoff_ddG, cutoff_ddG, comb$Rosetta_ddg_score_01)

p2 <- ggplot(comb, aes(x=gemme_score_00,y=Rosetta_ddg_score_01, color = score)) +
  geom_point(size = 3) +
  scale_color_gradient2(name = "MAVE Score", midpoint = 0, low = "#830823", high = "#063264", mid = 'white', 
                       limits = c(min_score,max_score)) +
  #ggtitle(paste0('One point per variant: ddG VS ddE scores colored by MAVE score\nddG values capped at ', cutoff_ddG)) +
  xlab('Gemme evolutionary conservation (ddE)') + ylab('Rosetta ddG') +
  theme(text = element_text(size = 20), legend.position = 'bottom') 

p2

png(filename = paste0(plot_dir,'ddG_VS_ddE_indiv_vars_cap',cutoff_ddG,'.png'),width = 1000, height = 800)
print(p2)
dev.off()

#add distr of dde and ddg values to plot
distr_y <- ggplot(comb, aes(Rosetta_ddg_score_01)) +
  geom_density(bw= 0.08) +
  #xlim(-3,20) +
  #ggtitle(paste0("tile ", tile, ", ddG values")) +
  #xlab('Bandwidth = 0.08') +
  coord_flip() +
  theme_bw(base_size = 20)  +
  theme(axis.title.y=element_blank())

distr_y

distr_x <- ggplot(comb, aes(gemme_score_00)) +
  geom_density(bw= 0.08) +
  ggtitle(paste0('One point per variant: ddG VS ddE scores colored by MAVE score\nddG values capped at ', cutoff_ddG)) +
  #xlim(-3,20) +
  #ggtitle(paste0("tile ", tile, ", ddG values")) +
  #xlab('Bandwidth = 0.08') +
  theme_bw(base_size = 20) +
  theme(axis.title.x=element_blank())

distr_x

#little blank square in the top right corner
gg_empty <- ggplot(comb, aes(x=gemme_score_00,y=Rosetta_ddg_score_01))+
  geom_blank() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        panel.background = element_blank())

png(filename = paste0(plot_dir,'ddG_VS_ddE_indiv_vars_cap',cutoff_ddG,'_plus_marginals.png'),width = 1000, height = 900)
ggarrange(
  distr_x, gg_empty, p2, distr_y,
  nrow = 2, ncol = 2, widths = c(5, 1), heights = c(1, 5)
)
dev.off()


#should do another one with one dot per position
# p3 <- ggplot(comb, aes(x=gemme_score_00,y=Rosetta_ddg_score_01, color = score)) +
#   geom_point() +
#   scale_color_gradient2(name = "MAVE Score", midpoint = 0, low = "#830823", high = "#063264", mid = 'white', 
#                         limits = c(min_score,max_score)) +
#   ggtitle('One point per position: ddG VS ddE scores colored by average MAVE score') +
#   xlab('Gemme evolutionary conservation (ddE)') + ylab('Rosetta ddG')
#   theme_bw(base_size = 20)
# 
# p3

#ddg VS MAVE
##########################################
min_score <- -2.5
max_score <- 4.0

#introduce max ddG cutoff
cutoff_ddG <- 10
comb$Rosetta_ddg_score_01 <- ifelse(comb$Rosetta_ddg_score_01 > cutoff_ddG, cutoff_ddG, comb$Rosetta_ddg_score_01)

p2 <- ggplot(comb, aes(x=score,y=Rosetta_ddg_score_01, color =score)) +
  geom_point() +
  scale_color_gradient2(name = "MAVE Score", midpoint = 0, low = "#830823", high = "#063264", mid = 'white', 
                        limits = c(min_score,max_score)) +
  ggtitle(paste0('One point per variant: ddG VS MAVE score\nddG values capped at ', cutoff_ddG)) +
  xlab('Mave score') + ylab('Rosetta ddG') +
  theme_bw(base_size = 20)

p2

png(filename = paste0(plot_dir,'indiv_mave_VS_ddG_cap',cutoff_ddG,'.png'),width = 1000, height = 800)
print(p2)
dev.off()

p2 <- ggplot(comb, aes(x=score,y=Rosetta_ddg_score_01, color = MTXdist)) +
  geom_point() +
  #scale_color_gradient2(name = "Distance to MTX", midpoint = 0, low = "#830823", high = "#063264", mid = 'white', 
  #                      limits = c(min_score,max_score)) +
  ggtitle(paste0('One point per variant: ddG VS MAVE score\nddG values capped at ', cutoff_ddG)) +
  xlab('Mave score') + ylab('Rosetta ddG') +
  theme(text = element_text(size = 20))

p2

png(filename = paste0(plot_dir,'indiv_mave_VS_ddG_cap',cutoff_ddG,'_color_mtxdist.png'),width = 1000, height = 800)
print(p2)
dev.off()

#ddE VS MAVE
##########################################

p2 <- ggplot(comb, aes(x=score,y=gemme_score_00, color = score)) +
  geom_point() +
  scale_color_gradient2(name = "MAVE Score", midpoint = 0, low = "#830823", high = "#063264", mid = 'white', 
                        limits = c(min_score,max_score)) +
  ggtitle(paste0('One point per variant: ddE VS MAVE score')) +
  xlab('Mave score') + ylab('Gemme evolutionary conservation (ddE)') +
  theme_bw(base_size = 20)

p2

png(filename = paste0(plot_dir,'indiv_mave_VS_ddE.png'),width = 1000, height = 800)
print(p2)
dev.off()

p2 <- ggplot(comb, aes(x=score,y=gemme_score_00, color = MTXdist)) +
  geom_point() +
  scale_color_gradient(name = "Distance to MTX", low = "darkred", high = "gold") +
  ggtitle(paste0('One point per variant: ddE VS MAVE score')) +
  xlab('Mave score') + ylab('Gemme evolutionary conservation (ddE)') +
  theme(text = element_text(size = 20), legend.position = 'bottom')

p2

png(filename = paste0(plot_dir,'indiv_mave_VS_ddE_color_mtxdist.png'),width = 1000, height = 800)
print(p2)
dev.off()

p2 <- ggplot(comb, aes(x=score,y=gemme_score_00, color = NADPHdist)) +
  geom_point() +
  scale_color_gradient(name = "Distance to NADPH", low = "chartreuse4", high = "darkolivegreen1") +
  ggtitle(paste0('One point per variant: ddE VS MAVE score')) +
  xlab('Mave score') + ylab('Gemme evolutionary conservation (ddE)') +
  theme(text = element_text(size = 20), legend.position = 'bottom')

p2

png(filename = paste0(plot_dir,'indiv_mave_VS_ddE_color_nadphdist.png'),width = 1000, height = 800)
print(p2)
dev.off()

p2 <- ggplot(comb, aes(x=score,y=gemme_score_00, color = dist_group)) +
  geom_point() +
  #scale_color_gradient(name = "Distance to NADPH", low = "chartreuse4", high = "darkolivegreen1") +
  ggtitle(paste0('One point per variant: ddE VS MAVE score')) +
  xlab('Mave score') + ylab('Gemme evolutionary conservation (ddE)') +
  theme_bw(base_size = 20) + theme(legend.position = 'bottom')

p2

png(filename = paste0(plot_dir,'indiv_mave_VS_ddE_color_distgroup.png'),width = 1000, height = 800)
print(p2)
dev.off()

#dist MTX VS dist NAD
############################

dist_vs <- merge(mtx_dist[,c('pos','dist', 'molecule')], nadph_dist[,c('pos','dist', 'molecule')], by = 'pos')
dist_vs$pos_median <- mat$pos_median

range(dist_vs$dist.x)
range(dist_vs$dist.y)

p <- ggplot(dist_vs, aes(x=dist.x,y=dist.y,color=pos_median)) +
  geom_point(size = 3) +
  scale_color_gradient2(name = "median MAVE Score", midpoint = 0, low = "#830823", high = "#063264", mid = 'white', 
                        limits = c(min_score,max_score)) +
  ggtitle(paste0('Comparing distances to ligands')) +
  xlab('Distance MTX') + ylab('Distance NADPH') +
  theme(text = element_text(size = 20), panel.background = element_rect(fill = 'grey70', color = 'grey70'))

png(filename = paste0(plot_dir,'dist_MTX_vs_dist_NADPH_by_median_MAVE.png'),width = 1000, height = 800)
print(p)
dev.off()

cutoff_dist <- 15
dist_vs$dist.x <- ifelse(dist_vs$dist.x > cutoff_dist, cutoff_dist, dist_vs$dist.x)
dist_vs$dist.y <- ifelse(dist_vs$dist.y > cutoff_dist, cutoff_dist, dist_vs$dist.y)

p <- ggplot(dist_vs, aes(x=dist.x,y=dist.y,color=pos_median)) +
  geom_point(size = 3) +
  scale_color_gradient2(name = "median MAVE Score", midpoint = 0, low = "#830823", high = "#063264", mid = 'white', 
                        limits = c(min_score,max_score)) +
  ggtitle(paste0('Comparing distances to ligands\nDistance cut off at ',cutoff_dist,' Angstrom')) +
  xlab('Distance MTX') + ylab('Distance NADPH') +
  theme(text = element_text(size = 20), panel.background = element_rect(fill = 'grey70', color = 'grey70'))

p
png(filename = paste0(plot_dir,'dist_MTX_vs_dist_NADPH_by_median_MAVE_cutoff',cutoff_dist,'.png'),width = 1000, height = 800)
print(p)
dev.off()


#distr of median MAVE score in positions close to the ligand (dist < 10) and far from the ligand
############################

dist_vs$dist_group_MTX <- ifelse(dist_vs$dist.x < 10,'close','far')
dist_vs$dist_group_NADPH <- ifelse(dist_vs$dist.y < 10,'close','far')

p <- ggplot(dist_vs, aes(pos_median,color=dist_group_MTX)) +
  #geom_density(bw= 0.08) +
  geom_boxplot() +
  coord_flip() +
  ggtitle(paste0("Distribution of median MAVE scores by distance from MTX")) +
  xlab('Median MAVE score per position') +
  theme_bw(base_size = 20) + theme(legend.position="bottom")

p
png(filename = paste0(base_dir,'HZ_plots/score_densities/by_dist_MTX.png'),width = 800, height = 550)
print(p)
dev.off()

p <- ggplot(dist_vs, aes(pos_median,color=dist_group_NADPH)) +
  geom_density(bw= 0.08) +
  ggtitle(paste0("Distribution of median MAVE scores by distance from NADPH")) +
  xlab('Median MAVE score per position\nBandwidth = 0.08') +
  theme_bw(base_size = 20) + theme(legend.position="bottom")

p
png(filename = paste0(base_dir,'HZ_plots/score_densities/by_dist_NADPH.png'),width = 800, height = 550)
print(p)
dev.off()

#both in the same plot
#######################################

dist_mave$group <- ifelse(dist_mave$dist < 10,'MTX_close','MTX_far')
dist_mave_nadph$group <- ifelse(dist_mave_nadph$dist < 10,'NADPH_close','NADPH_far')
temp_df <- rbind(dist_mave,dist_mave_nadph)
temp_df$group <- factor(temp_df$group, levels=c('MTX_close','NADPH_close','MTX_far','NADPH_far'))

p <- ggplot(temp_df, aes(x=group,y=pos_median,fill=group)) +
  #geom_boxplot() +
  geom_violin(scale = "count") +
  #geom_jitter(height = 0, width = 0.1) +
  ggtitle(paste0("Distribution of median MAVE scores by distance from the ligands")) +
  ylab('Median MAVE score per position') +
  theme_bw(base_size = 20) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), legend.position = 'bottom')


p
png(filename = paste0(base_dir,'HZ_plots/score_densities/violine_by_dist_NADPH_MTX.png'),width = 800, height = 550)
print(p)
dev.off()

p2 <- ggplot(rbind(dist_mave,dist_mave_nadph), aes(pos_median,color=group)) +
  geom_density(bw= 0.08) +
  ggtitle(paste0("Distribution of median MAVE scores by distance from the ligands")) +
  xlab('Median MAVE score per position\nBandwidth = 0.08') +
  theme_bw(base_size = 20) #+ theme(legend.position="bottom")

p2
png(filename = paste0(base_dir,'HZ_plots/score_densities/by_dist_NADPH_MTX.png'),width = 800, height = 550)
print(p2)
dev.off()

