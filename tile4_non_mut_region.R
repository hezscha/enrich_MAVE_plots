library(ggplot2)
library(stringr)
library(dplyr)
library(reshape2)

#plots counts on DNA variant level, i.e.counts per ALT base per DNA position.
#the read overlap of region 4 extends before and after the mutagenized region. I want to look at the distr of DNA level var counts 
#per aa pos. So for V116A, each DNA var that leads to that gives one count in the distr and same for all other aa vars at that pos
#I want to see how many counts vars at pos 116, 117 and 118 have because these should all not be there so they are likely seq errors/miscalls


#try out to get counts per alt base and for all selections
###########################

tile <- '4'
#rep <- 'C'
base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/results/enrich_full_dhfr/'

# for (tile in c('1','2','3','4','5')) {
#   
#   print(c('tile', tile))
all_sels <- data.frame()
for (sel in c('A','B','C')) {
  #DNA var level counts
  DNA_counts <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps/tsv/',sel,'_sel/main_variants_counts_unfiltered.tsv'), 
                               sep = '\t', header = F,
                               stringsAsFactors = F, skip = 2,
                               col.names = c('var', 'count_t0',	'count_t1',	'count_t2',	'count_t3'))
  
  #add sel
  DNA_counts$sel <- sel
  #add nr of subs to each var. OBS: This is number of subs on DNA level!
  DNA_counts$nr_subs <- str_count(DNA_counts$var,",") + 1
  #single vars
  single <- DNA_counts[DNA_counts$nr_subs == 1,]
  
  #position of DNA change
  single$dna_pos <- as.numeric(str_extract(single$var, "(?<=c\\.)[0-9]+(?=[A-Z])"))
  #alt base:
  single$alt_base <- str_extract(single$var, "(?<=c\\.[0-9]{1,3}[A-Z]\\>)[A-Z]")
  all_sels <- rbind(all_sels, single)
}

m_single <- melt(all_sels, id.vars = c('dna_pos', 'var', 'nr_subs', 'alt_base', 'sel'))

#try facet wrap instead
###############################
#part 1 of the tile 4 region 
part1 <- m_single[m_single$dna_pos <= 391, ]
part1$dna_pos <- as.factor(part1$dna_pos)

p4 <- ggplot(part1, aes(y=log10(value), x=alt_base, color=alt_base)) +
  geom_boxplot(position = position_dodge(0.5), width = 0.25) +
  ylab('log10 counts') +
  ggtitle(paste0('Tile ', tile, ': Counts of DNA level vars per DNA position for 346 - 391.',
                 '\nMutagenized region of tile 4 goes from 355 - 468, but we have read overlap coverage from 346 - 477.')) +
  #scale_x_continuous(breaks = seq(min(m_single$dna_pos), 412,by=1),
  #                   labels = seq(min(m_single$dna_pos), 412,by=1))+
  theme_bw(base_size = 20) +
  facet_wrap(~dna_pos, nrow = 1)

#p4

png(filename = paste0(base_dir,'HZ_plots/quality_control/DNA_var_counts_by_var_part1_tile',tile,'_allsels.facet.log.byDNApos.png'), 
    width = 5000, height = 1100, res = 200, units = 'px')
#png(filename = paste0(base_dir,'HZ_plots/quality_control/DNA_var_counts_by_var_part1_tile',tile,'_allsels.facet.log.byDNApos.png'), 
#   res = 300, units = 'px')
print(p4)
dev.off()

#part 2 of the whole tile 4 region
########

part2 <- m_single[m_single$dna_pos > 391 & m_single$dna_pos <= 436, ]
part2$dna_pos <- as.factor(part2$dna_pos)

p4 <- ggplot(part2, aes(y=log10(value), x=alt_base, color=alt_base)) +
  geom_boxplot(position = position_dodge(0.5), width = 0.25) +
  ylab('log10 counts') +
  ggtitle(paste0('Tile ', tile, ': Counts of DNA level vars per DNA position for 391 - 436.',
                 '\nMutagenized region of tile 4 goes from 355 - 468, but we have read overlap coverage from 346 - 477.')) +
  #scale_x_continuous(breaks = seq(min(m_single$dna_pos), 412,by=1),
  #                   labels = seq(min(m_single$dna_pos), 412,by=1))+
  theme_bw(base_size = 20) +
  facet_wrap(~dna_pos, nrow = 1)

#p4

png(filename = paste0(base_dir,'HZ_plots/quality_control/DNA_var_counts_by_var_part2_tile',tile,'_allsels.facet.log.byDNApos.png'), 
    width = 5000, height = 1100, res = 200, units = 'px')
print(p4)
dev.off()

#part 3 of the whole tile 4 region
##########
part3 <- m_single[m_single$dna_pos > 436, ]
part3$dna_pos <- as.factor(part3$dna_pos)

p4 <- ggplot(part3, aes(y=log10(value), x=alt_base, color=alt_base)) +
  geom_boxplot(position = position_dodge(0.5), width = 0.25) +
  ylab('log10 counts') +
  ggtitle(paste0('Tile ', tile, ': Counts of DNA level vars per DNA position for 436 - 477.',
                 '\nMutagenized region of tile 4 goes from 355 - 468, but we have read overlap coverage from 346 - 477.')) +
  #scale_x_continuous(breaks = seq(min(m_single$dna_pos), 412,by=1),
  #                   labels = seq(min(m_single$dna_pos), 412,by=1))+
  theme_bw(base_size = 20) +
  facet_wrap(~dna_pos, nrow = 1)

#p4

png(filename = paste0(base_dir,'HZ_plots/quality_control/DNA_var_counts_by_var_part3_tile',tile,'_allsels.facet.log.byDNApos.png'), 
    width = 5000, height = 1100, res = 200, units = 'px')
print(p4)
dev.off()


#bar plot
#################################
p <- ggplot(m_single[m_single$dna_pos <= 412, ], aes(y=value, x=dna_pos, fill=alt_base)) +
  geom_bar(position = 'dodge', stat= 'identity') +
  ylab('counts') +
  ggtitle(paste0('Tile4: Counts of DNA level vars per DNA position\nMutagenized region of tile 4 goes from 355 - 468, ',
                 'but we have read overlap coverage from 346 - 477.')) +
  scale_x_continuous(breaks = seq(min(m_single$dna_pos), 412,by=1),
                     labels = seq(min(m_single$dna_pos), 412,by=1))+
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90))

png(filename = paste0(base_dir,'HZ_plots/quality_control/DNA_var_counts_by_var_tile',tile,'_sel',sel,'.bar.byDNApos.png'), width = 2000, height = 550)
print(p)
dev.off()


##############################################################

part1 <- m_single[m_single$dna_pos <= 412, ]
part1$dna_pos <- as.factor(part1$dna_pos)

p2 <- ggplot(part1, aes(y=value, x=dna_pos, color=alt_base)) +
  geom_boxplot(position = position_dodge(0.5), width = 0.25) +
  ylab('counts') +
  ggtitle(paste0('Tile4: Counts of DNA level vars per DNA position\nMutagenized region of tile 4 goes from 355 - 468, ',
                 'but we have read overlap coverage from 346 - 477.')) +
  #scale_x_continuous(breaks = seq(min(m_single$dna_pos), 412,by=1),
  #                   labels = seq(min(m_single$dna_pos), 412,by=1))+
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90))

p2
png(filename = paste0(base_dir,'HZ_plots/quality_control/DNA_var_counts_by_var_tile',tile,'_sel',sel,'.box.byDNApos.png'), 
    width = 2500, height = 550)
print(p2)
dev.off()

p3 <- ggplot(part1, aes(y=log10(value), x=dna_pos, color=alt_base)) +
  geom_boxplot(position = position_dodge(0.5), width = 0.25) +
  ylab('log10 counts') +
  ggtitle(paste0('Tile4: Counts of DNA level vars per DNA position\nMutagenized region of tile 4 goes from 355 - 468, ',
                 'but we have read overlap coverage from 346 - 477.')) +
  #scale_x_continuous(breaks = seq(min(m_single$dna_pos), 412,by=1),
  #                   labels = seq(min(m_single$dna_pos), 412,by=1))+
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90))

png(filename = paste0(base_dir,'HZ_plots/quality_control/DNA_var_counts_by_var_tile',tile,'_sel',sel,'.box.log.byDNApos.png'), 
    width = 2500, height = 550)
print(p3)
dev.off()


# make DNA_var_counts_tile4.byDNApos.log.png 
###############################

tile <- '4'
rep <- 'C'
base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/results/enrich_full_dhfr/'

#DNA var level counts
A_sel_DNA_counts <- read.csv(paste0(base_dir,'tile', tile ,'_all_reps/tsv/',rep,'_sel/main_variants_counts_unfiltered.tsv'), 
                 sep = '\t', header = F,
                 stringsAsFactors = F, skip = 2,
                 col.names = c('var', 'count_t0',	'count_t1',	'count_t2',	'count_t3'))

#add nr of subs to each var. OBS: This is number of subs on DNA level!
A_sel_DNA_counts$nr_subs <- str_count(A_sel_DNA_counts$var,",") + 1

#single vars
single <- A_sel_DNA_counts[A_sel_DNA_counts$nr_subs == 1,]

#could also do a data structure that has two cols, count and position. So for every var, one line is created for each time point
# ct  pos
#  5   116 -> DNA_var_1_t0 : c.346G>A (p.Val116Ile) t0
#  1   116 -> DNA_var_1_t1 : c.346G>A (p.Val116Ile) t1
#  1   116 -> DNA_var_2_t0 : c.346G>C (p.Val116Leu) t0
#  7   117

#need to melt single so we get one line for each count and then add mutation aa pos by extracting it from the var
m_single <- melt(single)
m_single$aa_pos <- as.numeric(str_extract(m_single$var, "(?<=p\\.[A-Z][a-z]{2})([0-9]+)(?=[A-Z][a-z]{2})"))

p1 <- ggplot(m_single, aes(x=aa_pos, y=value)) +
  geom_boxplot() + ylab('Read counts') +ggtitle('Tile4: Counts of DNA level vars per protein level var\nMutagenized region of tile 4 goes from 119 - 156 but we have read overlap coverage from 116 - 159.\nNA pos are vars that lead to same protein as WT so no position could be extracted from the var name.') +
  theme_bw(base_size = 20)

p1

png(filename = paste0(base_dir,'HZ_plots/DNA_var_counts_tile',tile,'_sel',rep,'.png'), width = 1500, height = 550)
print(p1)
dev.off()

p2 <- ggplot(m_single, aes(x=aa_pos, y=log10(value))) +
  geom_boxplot() + ylab('Log10 read counts') +ggtitle('Tile4: Counts of DNA level vars per protein level var\nMutagenized region of tile 4 goes from 119 - 156 but we have read overlap coverage from 116 - 159.\nNA pos are vars that lead to same protein as WT so no position could be extracted from the var name.') +
  theme_bw(base_size = 20)

p2

png(filename = paste0(base_dir,'HZ_plots/DNA_var_counts_tile',tile,'_sel',rep,'.log.png'), width = 1500, height = 550)
print(p2)
dev.off()

#index on DNa level instead:
m_single$dna_pos <- str_extract(m_single$var, "(?<=c\\.)[0-9]+(?=[A-Z])")

p3 <- ggplot(m_single, aes(x=dna_pos, y=log10(value))) +
  geom_boxplot() + ylab('Log10 read counts') +
  ggtitle('Tile4: Counts of DNA level vars per DNA position\nMutagenized region of tile 4 goes from 355 - 466 but we have read overlap coverage from 346 - 477.') +
  theme_bw(base_size = 20) + theme(axis.text.x = element_text(angle = 90))

p3

png(filename = paste0(base_dir,'HZ_plots/DNA_var_counts_tile',tile,'_sel',rep,'.byDNApos.log.png'), width = 2000, height = 550)
print(p3)
dev.off()



#detail to str_extract and trying out on pos 116
#######################
#first step: get all DNA level vars that lead to V116X
pos_116 <- single[grepl('116',single$var, fixed = T),]
m_pos_116 <- melt(pos_116)

#extract DNA pos
str_extract(m_pos_116$var, "(?<=c\\.)[0-9]+(?=[A-Z])")
str_extract(m_pos_116$var, "(?<=p\\.[A-Z][a-z]{2})([0-9]{3})(?=[A-Z][a-z]{2})")

#trying to assemble df where each col is list of counts with 1 sub, 2 subs, 3 subs, ect but went off track
##################
#var is the first col, nr_subs is 6th. Can't figure out how to deselect column by name
single_cts <- as.numeric(as.list(unlist(single[,-c(1,6)])))
single_cts_log <- log10(single_cts)
hist(single_cts, labels = T)

#double vars
double <- A_sel_DNA_counts[A_sel_DNA_counts$nr_subs == 2,]
double_cts_log <- log10(as.numeric(as.list(unlist(double[,-c(1,6)]))))

#triple
three <- A_sel_DNA_counts[A_sel_DNA_counts$nr_subs == 3,]
three_cts_log <- log10(as.numeric(as.list(unlist(three[,-c(1,6)]))))

#four
four <- A_sel_DNA_counts[A_sel_DNA_counts$nr_subs == 4,]
four_cts_log <- log10(as.numeric(as.list(unlist(four[,-c(1,6)]))))

#assemble a df
lists = list(single_DNA_mut = single_cts_log,
             double_DNA_muts = double_cts_log,
             three_DNA_muts = three_cts_log,
             four_DNA_muts = four_cts_log)
read_counts <- as.data.frame(do.call(cbind, lists))






#trying stuff
###################

#time point 0 for now
tp0 <- A_sel_DNA_counts[,c('var','count_t0')]

#add nr of subs to each var. OBS: This is number of subs on DNA level!
tp0$nr_subs <- str_count(tp0$var,",") + 1
table(tp0$nr_subs)

#compare count distr between single and multivars
A_sel_DNA_counts_single <- A_sel_DNA_counts[!grepl(',',A_sel_DNA_counts$var, fixed = T),]
#table(A_sel_DNA_counts_single$nr_subs) #check that we have only single vars

A_sel_DNA_counts_multi <- A_sel_DNA_counts[!grepl(',',A_sel_DNA_counts$var, fixed = T),]


#lets look only at single vars



