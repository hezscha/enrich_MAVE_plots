library(ggplot2)

#OBS: change paths to what is applicable for you. You can find all lines where paths are mentioned by searching for 'base_dir'
#full DHFR seq: compare score between replicates to investigate the high error

tile <- '1'
base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/results/enrich_full_dhfr/'

rep1 <- 'A'
rep2 <- 'B'


make_versus_plot <- function(rep1, rep2) {
  
  dat1 <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps/tsv/',rep1,'_sel/main_synonymous_scores.tsv'), 
                   sep = '\t', header = F,
                   stringsAsFactors = F, skip = 2,
                   col.names = c('var', 'score_1',	'SE_1',	'SE_pctile_1',	'slope_1',	'intercept_1',	'SE_slope_1',	't_1', 'pvalue_raw_1'))
  
  dat2 <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps/tsv/',rep2,'_sel/main_synonymous_scores.tsv'), 
                   sep = '\t', header = F,
                   stringsAsFactors = F, skip = 2,
                   col.names = c('var', 'score_2',	'SE_2',	'SE_pctile_2',	'slope_2',	'intercept_2',	'SE_slope_2',	't_2', 'pvalue_raw_2'))
  
  shared <- nrow(merge(dat1[,c('var', 'score_1')], dat2[,c('var','score_2')], by = 'var'))
  m <- merge(dat1[,c('var', 'score_1')], dat2[,c('var','score_2')], by = 'var', all = T)
  m_single <- m[!grepl(',',m$var, fixed = T),]
  shared_single <- nrow(m_single[(!is.na(m_single$score_1) & !is.na(m_single$score_2)),])
  
  min_score <- min(m[,c('score_1','score_2')], na.rm = T) - 0.1
  
  #put an arbitrary but out of the range score instead of NA to include those dots
  m[is.na(m)] <- 6
  m_single[is.na(m_single)] <- 6
  
  p1 <- ggplot(m, aes(x=score_1,y=score_2))+
    geom_point() +
    xlab(paste0('scores ', rep1)) + ylab(paste0('scores ', rep2)) +
    xlim(min_score,6.5) + ylim(min_score,6.5) +
    ggtitle(paste0("Comparing replicate scores on run 'tile1_all_reps'. All vars.\nNA set to 6.\n",
                   "Shared vars: ", shared, 
                   ', only ', rep1, ': ', nrow(m[m$score_2 == 6,]), 
                   ', only ', rep2, ': ', nrow(m[m$score_1 == 6,]))) +
    theme_bw(base_size = 20)
  
  png(filename = paste0(base_dir,'HZ_plots/scores_', rep1, '_VS_', rep2,'.png'), width = 800, height = 600)
  print(p1)
  dev.off()
  
  p2 <- ggplot(m_single, aes(x=score_1,y=score_2))+
    geom_point() +
    xlab(paste0('scores ', rep1)) + ylab(paste0('scores ', rep2)) +
    xlim(min_score,6.5) + ylim(min_score,6.5) +
    ggtitle(paste0("Comparing replicate scores on run 'tile1_all_reps'. Single vars.\nNA set to 6.\n",
                   "Shared vars: ", shared_single, 
                   ', only ', rep1, ': ', nrow(m_single[m_single$score_2 == 6,]), 
                   ', only ', rep2, ': ', nrow(m_single[m_single$score_1 == 6,]))) +
    theme_bw(base_size = 20)
  
  png(filename = paste0(base_dir,'HZ_plots/scores_', rep1, '_VS_', rep2,'.single.png'), width = 800, height = 600)
  print(p2)
  dev.off()
  
  
}

p1 <- make_versus_plot('B','C')
p1


# A_vs_B <- merge(datA[,c('var', 'score_A')], datB[,c('var','score_B')], by = 'var')
# #need all because I also want vars that only have score in either A or B
# A_vs_B_all <- merge(datA[,c('var', 'score_A')], datB[,c('var','score_B')], by = 'var', all = T)
# #put an arbitrary but out of the range score instead of NA to include those dots
# A_vs_B_all[is.na(A_vs_B_all)] <- 6
# 
# #the below confirms there are no vars that are NA in both A and B
# #sum(is.na(A_vs_B_all$score_A) && is.na(A_vs_B_all$score_B))
# #bla <- A_vs_B_all[is.na(A_vs_B_all$score_B),]
# 
# #max_score <- max(A_vs_B_all[,c('score_A','score_B')], na.rm = T) + 0.1
# min_score <- min(A_vs_B_all[,c('score_A','score_B')], na.rm = T) - 0.1
# 
# p1 <- ggplot(A_vs_B_all, aes(x=score_A,y=score_B))+
#   geom_point() +
#   #xlim(min_score,max_score) + ylim(min_score,max_score) +
#   xlim(min_score,6.5) + ylim(min_score,6.5) +
#   ggtitle(paste0("Comparing replicate scores on run 'tile1_all_reps'. All vars.\nNA set to 6.\n",
#                  "Shared vars: ", nrow(A_vs_B), 
#                  ', only A: ', nrow(A_vs_B_all[A_vs_B_all$score_B == 6,]), 
#                  ', only B: ', nrow(A_vs_B_all[A_vs_B_all$score_A == 6,]))) +
#   theme_bw(base_size = 20)
# 
# png(filename = paste0(base_dir,'HZ_plots/scores_A_VS_B.png'), width = 800, height = 600)
# print(p1)
# dev.off()
