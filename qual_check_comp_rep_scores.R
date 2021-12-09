library(ggplot2)

#full DHFR seq: compare score between replicates (called selections in Enrich) to investigate the error/how sure we are of the experiment-level scores

base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/results/enrich_full_dhfr/'
corr_descr <- 'corrected with syn counts' # #'corrected with WT counts' # 'corrected with all reads' # 'corrected with counted reads'
corr <- '_corr_sy' #'_corr_complete' #'' #'_corr_full' #'_corr_complete'

rep1 <- 'A'
rep2 <- 'C'

make_versus_plot_all_tiles <- function(rep1, rep2){
  #read in prepped data: one long df with one score for each var across all tiles, but per selection. 
  #vars scored in several tiles were averaged and their errors propagated
  dat1 <- read.csv(paste0(base_dir, '/combined_R_dataframes/all_tile_scores_long_rep',rep1,'.csv'), sep = ',',
                   col.names = c('var', 'score_1',	'SE_1', 'tile', 'pos', 'sub'))
  #colnames(dat1) <- paste0(colnames(dat1), '_1')

  dat2 <- read.csv(paste0(base_dir, '/combined_R_dataframes/all_tile_scores_long_rep',rep2,'.csv'), sep = ',',
                   col.names = c('var', 'score_2',	'SE_2', 'tile', 'pos', 'sub'))
  #colnames(dat2) <- paste0(colnames(dat2), '_2')
  
  #combine data
  shared <- nrow(merge(dat1[,c('var', 'score_1')], dat2[,c('var','score_2')], by = 'var'))
  m <- merge(dat1[,c('var', 'score_1')], dat2[,c('var','score_2')], by = 'var', all = T)
  #m_single <- m[!grepl(',',m$var, fixed = T),]
  #shared_single <- nrow(m_single[(!is.na(m_single$score_1) & !is.na(m_single$score_2)),])
  
  min_score <- min(m[,c('score_1','score_2')], na.rm = T) - 0.1
  
  #cor_scores_single <- cor(m_single$score_1, m_single$score_2, use = 'complete.obs')
  cor_scores_all <- cor(m$score_1, m$score_2, use = 'complete.obs')
  
  #put an arbitrary but out of the range score instead of NA to include those dots
  m[is.na(m)] <- 6
  #m_single[is.na(m_single)] <- 6
  
  #make plot
  #I have taken out double/multi variants when I compiled the long tables for all tiles, therefore there is only a single var plot
  p1 <- ggplot(m, aes(x=score_1,y=score_2))+
    geom_point() +
    xlab(paste0('scores ', rep1)) + ylab(paste0('scores ', rep2)) +
    xlim(min_score,6.5) + ylim(min_score,6.5) +
    ggtitle(paste0("Comparing replicate scores on all tiles ", corr_descr, ".\nSingle vars. NA set to 6.",
                   "Shared vars: ", shared, 
                   ', only ', rep1, ': ', nrow(m[m$score_2 == 6,]), 
                   ', only ', rep2, ': ', nrow(m[m$score_1 == 6,]))) +
    annotate("text", x=6, y=6, label= paste0("pcc: ", round(cor_scores_all, 3)), size=7) +
    theme_bw(base_size = 20)
  
  png(filename = paste0(base_dir,'HZ_plots/quality_control/scores_all_tiles',corr, '_', rep1, '_VS_', rep2,'.single.png'), 
      width = 800, height = 600)
  print(p1)
  dev.off()
  
  # p2 <- ggplot(m_single, aes(x=score_1,y=score_2))+
  #   geom_point() +
  #   xlab(paste0('scores ', rep1)) + ylab(paste0('scores ', rep2)) +
  #   xlim(min_score,6.5) + ylim(min_score,6.5) +
  #   ggtitle(paste0("Comparing replicate scores on all tiles ", corr_descr, ".\nSingle vars. NA set to 6.",
  #                  "Shared vars: ", shared_single, 
  #                  ', only ', rep1, ': ', nrow(m_single[m_single$score_2 == 6,]), 
  #                  ', only ', rep2, ': ', nrow(m_single[m_single$score_1 == 6,]))) +
  #   annotate("text", x=6, y=6, label= paste0("cor: ", round(cor_scores_single, 2)), size=7) +
  #   theme_bw(base_size = 20)
  # 
  # png(filename = paste0(base_dir,'HZ_plots/quality_control/scores_all_tiles', corr, '_', rep1, '_VS_', rep2,'.single.png'), 
  #     width = 800, height = 600)
  # print(p2)
  # dev.off()
  
}

make_versus_plot <- function(rep1, rep2) {
  
  for (tile in c('1','2','3','4','5')) {
    dat1 <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps_minvarcount10', corr,'/tsv/',rep1,'_sel/main_synonymous_scores.tsv'), 
                     sep = '\t', header = F,
                     stringsAsFactors = F, skip = 3,
                     col.names = c('var', 'score_1',	'SE_1',	'SE_pctile_1',	'slope_1',	'intercept_1',	'SE_slope_1',	't_1', 'pvalue_raw_1'))
    
    dat2 <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps_minvarcount10',corr,'/tsv/',rep2,'_sel/main_synonymous_scores.tsv'), 
                     sep = '\t', header = F,
                     stringsAsFactors = F, skip = 3,
                     col.names = c('var', 'score_2',	'SE_2',	'SE_pctile_2',	'slope_2',	'intercept_2',	'SE_slope_2',	't_2', 'pvalue_raw_2'))
    
    shared <- nrow(merge(dat1[,c('var', 'score_1')], dat2[,c('var','score_2')], by = 'var'))
    m <- merge(dat1[,c('var', 'score_1')], dat2[,c('var','score_2')], by = 'var', all = T)
    m_single <- m[!grepl(',',m$var, fixed = T),]
    shared_single <- nrow(m_single[(!is.na(m_single$score_1) & !is.na(m_single$score_2)),])
    
    min_score <- min(m[,c('score_1','score_2')], na.rm = T) - 0.1
    
    cor_scores_single <- cor(m_single$score_1, m_single$score_2, use = 'complete.obs')
    cor_scores_all <- cor(m$score_1, m$score_2, use = 'complete.obs')
    
    #put an arbitrary but out of the range score instead of NA to include those dots
    m[is.na(m)] <- 6
    m_single[is.na(m_single)] <- 6
    
    p1 <- ggplot(m, aes(x=score_1,y=score_2))+
      geom_point() +
      xlab(paste0('scores ', rep1)) + ylab(paste0('scores ', rep2)) +
      xlim(min_score,6.5) + ylim(min_score,6.5) +
      ggtitle(paste0("Comparing replicate scores on tile ", tile, " ", corr_descr, ".\nAll vars. NA set to 6.",
                     "Shared vars: ", shared, 
                     ', only ', rep1, ': ', nrow(m[m$score_2 == 6,]), 
                     ', only ', rep2, ': ', nrow(m[m$score_1 == 6,]))) +
      annotate("text", x=6, y=6, label= paste0("cor: ", round(cor_scores_all, 2)), size=7) +
      theme_bw(base_size = 20)
    
    png(filename = paste0(base_dir,'HZ_plots/quality_control/scores_tile', tile, corr, '_', rep1, '_VS_', rep2,'.png'), width = 800, height = 600)
    print(p1)
    dev.off()
    
    p2 <- ggplot(m_single, aes(x=score_1,y=score_2))+
      geom_point() +
      xlab(paste0('scores ', rep1)) + ylab(paste0('scores ', rep2)) +
      xlim(min_score,6.5) + ylim(min_score,6.5) +
      ggtitle(paste0("Comparing replicate scores on tile ", tile, " ", corr_descr, ".\nSingle vars. NA set to 6.",
                     "Shared vars: ", shared_single, 
                     ', only ', rep1, ': ', nrow(m_single[m_single$score_2 == 6,]), 
                     ', only ', rep2, ': ', nrow(m_single[m_single$score_1 == 6,]))) +
      annotate("text", x=6, y=6, label= paste0("cor: ", round(cor_scores_single, 2)), size=7) +
      theme_bw(base_size = 20)
    
    png(filename = paste0(base_dir,'HZ_plots/quality_control/scores_tile', tile, corr, '_', rep1, '_VS_', rep2,'.single.png'), width = 800, height = 600)
    print(p2)
    dev.off()
  
  }
  
}

make_versus_plot('A','C')
make_versus_plot_all_tiles('B','C')



