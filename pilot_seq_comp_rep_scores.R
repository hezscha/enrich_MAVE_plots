library(ggplot2)

#compare replicate scores of the pilot sequencing
#first DHFR tile seq, 4 different conditions: compare score between replicates
#OBS: change paths to what is applicable for you. You can find all lines where paths are mentioned by searching for 'base_dir'


tile <- '1'
base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/results/enrich/'

rep1 <- 'ex2'
rep2 <- 'ex3'
cond <- 'MTX'

make_versus_plot <- function(rep1, rep2, cond, NA_score) {

  dat1 <- read.csv(paste0(base_dir,'latest_re-run/ex123_all_wo_20_22_38/tsv/',rep1,'_',cond,'_sel/main_synonymous_scores.tsv'), 
                   sep = '\t', header = F,
                   stringsAsFactors = F, skip = 2,
                   col.names = c('var', 'score_1',	'SE_1',	'SE_pctile_1',	'slope_1',	'intercept_1',	'SE_slope_1',	't_1', 'pvalue_raw_1'))
  
  dat2 <- read.csv(paste0(base_dir,'latest_re-run/ex123_all_wo_20_22_38/tsv/',rep2,'_',cond,'_sel/main_synonymous_scores.tsv'), 
                   sep = '\t', header = F,
                   stringsAsFactors = F, skip = 2,
                   col.names = c('var', 'score_2',	'SE_2',	'SE_pctile_2',	'slope_2',	'intercept_2',	'SE_slope_2',	't_2', 'pvalue_raw_2'))
  
  shared <- nrow(merge(dat1[,c('var', 'score_1')], dat2[,c('var','score_2')], by = 'var'))
  m <- merge(dat1[,c('var', 'score_1')], dat2[,c('var','score_2')], by = 'var', all = T)
  m_single <- m[!grepl(',',m$var, fixed = T),]
  shared_single <- nrow(m_single[(!is.na(m_single$score_1) & !is.na(m_single$score_2)),])
  
  min_score <- min(m[,c('score_1','score_2')], na.rm = T) - 0.1
  
  #put an arbitrary but out of the range score instead of NA to include those dots
  m[is.na(m)] <- NA_score
  m_single[is.na(m_single)] <- NA_score
  
  p1 <- ggplot(m, aes(x=score_1,y=score_2))+
    geom_point() +
    xlab(paste0('scores ', rep1)) + ylab(paste0('scores ', rep2)) +
    xlim(min_score,NA_score+0.5) + ylim(min_score,NA_score+0.5) +
    ggtitle(paste0("Comparing replicate scores on condition ", cond, ". All vars.\nNA set to ", NA_score, ".\n",
                   "Shared vars: ", shared, 
                   ', only ', rep1, ': ', nrow(m[m$score_2 == NA_score,]), 
                   ', only ', rep2, ': ', nrow(m[m$score_1 == NA_score,]))) +
    theme_bw(base_size = 20)
  
  png(filename = paste0(base_dir,'HZ_plots/scores_', cond, '_' , rep1, '_VS_', rep2,'.png'), width = 800, height = 600)
  print(p1)
  dev.off()
  
  p2 <- ggplot(m_single, aes(x=score_1,y=score_2))+
    geom_point() +
    xlab(paste0('scores ', rep1)) + ylab(paste0('scores ', rep2)) +
    xlim(min_score,NA_score+0.5) + ylim(min_score,NA_score+0.5) +
    ggtitle(paste0("Comparing replicate scores on condition ", cond, ". Single vars.\nNA set to ", NA_score, ".\n",
                   "Shared vars: ", shared_single, 
                   ', only ', rep1, ': ', nrow(m_single[m_single$score_2 == NA_score,]), 
                   ', only ', rep2, ': ', nrow(m_single[m_single$score_1 == NA_score,]))) +
    theme_bw(base_size = 20)
  
  png(filename = paste0(base_dir,'HZ_plots/scores_', cond, '_' , rep1, '_VS_', rep2,'.single.png'), width = 800, height = 600)
  print(p2)
  dev.off()
  
  
}

p1 <- make_versus_plot(rep1 = rep1,rep2 = rep2, cond = cond, 7)
p1
