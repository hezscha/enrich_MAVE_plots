library(ggplot2)
library(stringr)

#compare MAVE scores for vars that are caputred in two tiles

base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/results/enrich_full_dhfr/'
corr_descr <- 'corrected with counted reads'#'corrected with WT counts' # 'corrected with all reads' # 'corrected with counted reads'
corr <-'_corr_complete' #'' #'_corr_full' #'_corr_complete'

order_aa <- c('His', 'Lys', 'Arg', 'Asp', 'Glu', 'Cys', 'Met', 'Asn', 'Gln', 'Ser', 'Thr', 'Ala', 
              'Ile', 'Leu', 'Val', 'Phe', 'Trp', 'Tyr', 'Gly', 'Pro', 'Ter', 'pos_median')

#tiles 1 vs 2: overlap is aa's 37, 38 , 39
#-> there are 3 vars in common between tiles 1 and 2
##################################

tile <- '1'
dat_t1 <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps/tsv/tile', tile ,'_all_reps_exp/main_synonymous_scores.tsv'), 
                sep = '\t', header = F,
                stringsAsFactors = F, skip = 2,
                col.names = c('var', paste0('SE_', tile), paste0('epsilon_', tile),  paste0('score_', tile)))
#add nr of subs col and aa_pos col
dat_t1$aa_pos <- as.numeric(str_extract(dat_t1$var, "(?<=p\\.[A-Z][a-z]{2})([0-9]{2})(?=[A-Z][a-z]{2})"))
dat_t1$nr_subs <- str_count(dat_t1$var,",") + 1

tile <- '2'
dat_t2 <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps/tsv/tile', tile,'_all_reps_exp/main_synonymous_scores.tsv'), 
                sep = '\t', header = F,
                stringsAsFactors = F, skip = 2,
                col.names = c('var', paste0('SE_', tile), paste0('epsilon_', tile),  paste0('score_', tile)))
#add nr of subs col and aa_pos col
dat_t2$aa_pos <- as.numeric(str_extract(dat_t2$var, "(?<=p\\.[A-Z][a-z]{2})([0-9]{2})(?=[A-Z][a-z]{2})"))
dat_t2$nr_subs <- str_count(dat_t2$var,",") + 1

#subset to the overlap region
dat_t1_overlap <- dat_t1[dat_t1$aa_pos >= 37 & !(is.na(dat_t1$aa_pos)),] 
dat_t2_overlap <- dat_t2[dat_t2$aa_pos <= 39 & !(is.na(dat_t2$aa_pos)),] 

m_match <- merge(dat_t1_overlap[, c('var','SE_1','score_1')], dat_t2_overlap[, c('var','SE_2', 'score_2')], by ='var')
m_all <- merge(dat_t1_overlap[, c('var', 'score_1')], dat_t2_overlap[, c('var', 'score_2')], by ='var', all = T)

max_score <- max(m_match$SE_1,m_match$score_2) + 0.25
min_score <- min(m_match$SE_1,m_match$score_2) - 0.25

p <- ggplot(m_match, aes(x=score_1, y=score_2)) +
  geom_point() +
  geom_errorbar(aes(ymin=score_2-SE_2, ymax=score_2+SE_2), width=.2) +
  geom_errorbar(aes(xmin=score_1-SE_1, xmax=score_1+SE_1), width=.2) +
  xlim(min_score, max_score) + ylim(min_score, max_score) +
  ggtitle('Vars scored in both tile 1 and tile 2') +
  geom_text(aes(label=var),hjust=1,vjust=-1)+
  theme_bw(base_size = 20) 

p
png(filename = paste0(base_dir,'HZ_plots/quality_control/overlap_tiles1_2.png'),width = 800, height = 550)
print(p)
dev.off()

#tiles 2 vs 3: no overlap. Tile 2 has overlap data until S72 and tile 3 has overlap data from K81
##################################

#tiles 3 vs 4: no overlap. Tile 3 has overlap data until V116 and tile 4's mutagenized region starts from S119 (overlap data from 116, but those should have no mutations since that region is not mutagenize)
##################################

#tiles 4 vs 5: overlap is aa's 153, 154, 155, and 156. 
##################################

tile1 <- '4'
#dat_t1 <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps/tsv/tile', tile ,'_all_reps_exp/main_synonymous_scores.tsv'), 
dat_t1 <- read.csv(paste0(base_dir,'tile', tile1 ,'_all_reps_minvarcount10',corr,'/tsv/tile',tile1 ,'_all_reps_exp/main_synonymous_scores.tsv'), 
                   sep = '\t', header = F,
                   stringsAsFactors = F, skip = 2,
                   col.names = c('var', paste0('SE_', tile1), paste0('epsilon_', tile1),  paste0('score_', tile1)))
#add nr of subs col and aa_pos col
dat_t1$aa_pos <- as.numeric(str_extract(dat_t1$var, "(?<=p\\.[A-Z][a-z]{2})([0-9]+)(?=[A-Z][a-z]{2})"))
dat_t1$nr_subs <- str_count(dat_t1$var,",") + 1

tile2 <- '5'
#dat_t2 <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps/tsv/tile', tile,'_all_reps_exp/main_synonymous_scores.tsv'), 
dat_t2 <- read.csv(paste0(base_dir,'tile', tile2 ,'_all_reps_minvarcount10',corr,'/tsv/tile',tile2 ,'_all_reps_exp/main_synonymous_scores.tsv'), 
                   sep = '\t', header = F,
                   stringsAsFactors = F, skip = 2,
                   col.names = c('var', paste0('SE_', tile2), paste0('epsilon_', tile2),  paste0('score_', tile2)))
#add nr of subs col and aa_pos col
dat_t2$aa_pos <- as.numeric(str_extract(dat_t2$var, "(?<=p\\.[A-Z][a-z]{2})([0-9]+)(?=[A-Z][a-z]{2})"))
dat_t2$nr_subs <- str_count(dat_t2$var,",") + 1


#subset to the overlap region
dat_t1_overlap <- dat_t1[dat_t1$aa_pos >= 153 & !(is.na(dat_t1$aa_pos)),] 
dat_t2_overlap <- dat_t2[dat_t2$aa_pos <= 156 & !(is.na(dat_t2$aa_pos)),] 

m_match <- merge(dat_t1_overlap[, c('var','SE_4','score_4')], dat_t2_overlap[, c('var','SE_5', 'score_5')], by ='var')
m_all <- merge(dat_t1_overlap[, c('var','SE_4','score_4')], dat_t2_overlap[, c('var','SE_5', 'score_5')], by ='var', all = T)

max_score <- max(m_match$SE_4,m_match$score_5) + 0.25
min_score <- min(m_match$SE_4,m_match$score_5) - 0.25

cor_4_5 <- cor(m_match$score_4, m_match$score_5)

p <- ggplot(m_match, aes(x=score_4, y=score_5)) +
  geom_point() +
  geom_errorbar(aes(ymin=score_5-SE_5, ymax=score_5+SE_5), width=.2, linetype = 3, color = 'grey50') +
  geom_errorbar(aes(xmin=score_4-SE_4, xmax=score_4+SE_4), width=.2, linetype = 3, color = 'grey50') +
  xlim(min_score, max_score) + ylim(min_score, max_score) +
  ggtitle(paste0('Vars scored in both tile ', tile1, ' and tile ', tile2, ',\nminimum DNA var count 10 and scores ', corr_descr)) +
  #geom_text(aes(label=var),hjust=1,vjust=-1)+
  annotate("text", x=min_score+0.1, y=min_score+0.1, label= paste0("cor: ", round(cor_4_5, 2)), size=7) + 
  theme_bw(base_size = 20) 

p
png(filename = paste0(base_dir,'HZ_plots/quality_control/overlap_tiles4_5minvarcount10', corr, '.png'),width = 800, height = 550)
print(p)
dev.off()




