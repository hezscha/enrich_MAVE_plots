library(ggplot2)
library(stringr)
library(reshape2)
library(ggrepel)
library(dplyr) #arrange
library(hash)
#add marginal distr to dde vs ddg plot
library(gridExtra)
library(egg)
#library(gghalves) #more customizable jitter for rain clouds
#library(ggdist) #halfeye for rain clouds

#plot about the way vars are categorized into beneficial, del, WT-like, NA

###################################################
#1. functions and global vars
##############################################
base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/results/enrich_full_dhfr/'
data_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/data/'
plot_dir <- paste0(base_dir, '/HZ_plots/category_plots/')

# class_palette_cb <- c('del' = "#D55E00", #orange-red 
#                       'ben' = "#0072B2", #a dark blue
#                       'WT-like' = 'lightgoldenrodyellow') #cream/off-white

class_palette_nona <- c('del' = "#D55E00", #orange-red 
                       'ben' = "#0072B2", #a dark blue
                       'WT-like' = 'lightgoldenrod2')

class_palette_red <- c('del' = "#D55E00", #orange-red 
                       'ben' = "#0072B2", #a dark blue
                       'WT-like' = 'lightgoldenrod2',
                       'unclassified' = 'azure4')

class_palette_full <- c('del' = "#D55E00", #orange-red 
                       'ben' = "#0072B2", #a dark blue
                       'low_del' = 'indianred2',
                       'low_ben' = 'darkolivegreen2',
                       'WT-like' = 'lightgoldenrod2',
                       'unclassified' = 'azure4')

###################################################
#2. load data
###################################################

all_df_full <- read.csv(paste0(data_dir, 'full_MAVE_data.csv'))
#all_df_full$class_no_na <- ifelse(is.na(all_df_full$classification), 'unclassified', all_df_full$classification)
#all_df_full$tile_no_shared <- ifelse(all_df_full$tile == '4+5', '4', all_df_full$tile)

#subset to only rows that contain MAVE scores. There should be 3106 rows left 
all_df_full <- all_df_full[!is.na(all_df_full$MAVE_score),]

#make another cloumn where the NA class has a real name so it doesn't drop and is displayed in plots
all_df_full$class_no_na <- ifelse(is.na(all_df_full$classification), 'unclassified', all_df_full$classification)

#for some plots, re-assign shared vars to tile 4. 
#This is sort of like saying we're here assigning the overlap region
#to be tile 4, since this is region-wise way of looking at variants.
all_df_full$tile_no_shared <- ifelse(all_df_full$tile == '4+5', '4', all_df_full$tile)

#subset to omit vars in the low_ben and low_del categories. These we expect to functionally behave like WT
#this also omits vars where the classification is NA (this is a characteristic of subset. It doesn't happen when using [] to subset the df)
#but when using [], it sets the entire row to NA???
#all_df_sub1 <- all_df[(!all_df$classification == 'low_ben' & !all_df$classification == 'low_del'), ]
all_df <- subset(all_df_full, (!class_no_na == 'low_ben' & !class_no_na == 'low_del'))

#re-assign NA class vars to WT-like
all_df$na_is_wt <- ifelse(is.na(all_df$classification), 'WT-like', all_df$classification)

#add distance to ligands grouping
all_df$dist_group <- ifelse((all_df$MTXdist <= 10 & all_df$NADPHdist <= 10), 'both_close', 
                            ifelse((all_df$MTXdist <= 10 & all_df$NADPHdist > 10), 'MTX_close', 
                                   ifelse((all_df$MTXdist > 10 & all_df$NADPHdist <= 10), 'NADPH_close',
                                          'both_far')))


all_df$tile_no_shared <- ifelse(all_df$tile == '4+5', '4', all_df$tile)

###################################################
#3. plot
###################################################

#3.1. MAVE score per classification group
########################################

p <- ggplot(all_df, aes(x=MAVE_score,y=classification)) +
  geom_point(alpha = 0.6) +
  theme_bw(base_size = 16)

p




#how many benefical vars come from tile 1, 2 ect
#####################################################################



p3 <- ggplot(all_df, aes(x=classification,fill = tile_no_shared)) +
  geom_bar(position="dodge", alpha = 0.7, color = 'grey50') +
  scale_fill_discrete(name = 'tile') +
  theme_bw(base_size = 16) +
  ggtitle('Variants classed per tile')

p3

  png(filename = paste0(plot_dir, 'MAVE_tiles_per_category_all.png'),
    width = 800, height = 550, res = 100)
print(p3)
dev.off()  

p3 <- ggplot(all_df, aes(x=na_is_wt,fill = tile_no_shared)) +
  geom_bar(position="dodge", alpha = 0.7, color = 'grey50') +
  scale_fill_discrete(name = 'tile') +
  theme_bw(base_size = 16) +
  ggtitle('Variants classed per tile') + xlab('classification')

p3

png(filename = paste0(plot_dir, 'MAVE_tiles_per_category.png'),
    width = 800, height = 550, res = 100)
print(p3)
dev.off()  



#var classes per tile
########################
#to get labels of the total number of scored variants in each tile into the plot below
tot_vars <- as.data.frame(table(all_df$tile_no_shared))
colnames(tot_vars) <- c('tile','total')

p4 <- ggplot(all_df, aes(x=tile_no_shared,fill = class_no_na)) +
  geom_bar(position="fill", alpha = 0.7, color = 'grey50') +
  scale_fill_manual(values = class_palette_red, name = 'classification') +
  ggtitle('Fraction of variants sorted into') +
  ylab('fraction') + xlab('tile') +
  theme_bw(base_size = 16) +
  annotate("text", x = tot_vars$tile, y = 1.05,
           label = tot_vars$total)
p4

png(filename = paste0(plot_dir, 'MAVE_categories_per_tile_fraction.png'),
    width = 800, height = 550, res = 100)
print(p4)
dev.off()                                       

#to get labels of the total number of scored variants in each tile into the plot below
#tot_vars <- as.data.frame(table(all_df$tile))
#colnames(tot_vars) <- c('tile','total')

p5 <- ggplot(all_df, aes(x=tile_no_shared,fill = class_no_na)) +
  geom_bar(position="dodge", alpha = 0.7, color = 'grey50') +
  scale_fill_manual(values = class_palette_red, name = 'classification') +
  ggtitle('Variant classifications per tile') +
  ylab('No. variants') + xlab('tile') +
  theme_bw(base_size = 16) 

p5

png(filename = paste0(plot_dir, 'MAVE_categories_per_tile.png'),
    width = 800, height = 550, res = 100)
print(p5)
dev.off()  

#na as WT
#to get labels of the total number of scored variants in each tile into the plot below
tot_vars <- as.data.frame(table(all_df$tile_no_shared))
colnames(tot_vars) <- c('tile','total')

p4 <- ggplot(all_df, aes(x=tile_no_shared,fill = na_is_wt)) +
  geom_bar(position="fill", alpha = 0.7, color = 'grey50') +
  scale_fill_manual(values = class_palette_nona, name = 'classification') +
  ggtitle('Fraction of variants sorted into') +
  ylab('fraction') + xlab('tile') +
  theme_bw(base_size = 16) +
  annotate("text", x = tot_vars$tile, y = 1.05,
           label = tot_vars$total)
p4

png(filename = paste0(plot_dir, 'MAVE_categories_per_tile_fraction_naiswt.png'),
    width = 800, height = 550, res = 100)
print(p4)
dev.off()                                       

#to get labels of the total number of scored variants in each tile into the plot below
#tot_vars <- as.data.frame(table(all_df$tile))
#colnames(tot_vars) <- c('tile','total')

p5 <- ggplot(all_df, aes(x=tile_no_shared,fill = na_is_wt)) +
  geom_bar(position="dodge", alpha = 0.7, color = 'grey50') +
  scale_fill_manual(values = class_palette_nona, name = 'classification') +
  ggtitle('Variant classifications per tile') +
  ylab('No. variants') + xlab('tile') +
  theme_bw(base_size = 16) 

p5

png(filename = paste0(plot_dir, 'MAVE_categories_per_tile_naiswt.png'),
    width = 800, height = 550, res = 100)
print(p5)
dev.off()

#all categories (including NA and low_ben, low_del)
p5 <- ggplot(all_df_full, aes(x=tile_no_shared,fill = class_no_na)) +
  geom_bar(position="dodge", alpha = 0.7, color = 'grey50') +
  scale_fill_manual(values = class_palette_full, name = 'classification') +
  ggtitle('Variant classifications per tile') +
  ylab('No. variants') + xlab('tile') +
  theme_bw(base_size = 16) 

p5

png(filename = paste0(plot_dir, 'MAVE_categories_all_per_tile.png'),
    width = 800, height = 550, res = 100)
print(p5)
dev.off()  


#var classes per dist group
#################################

#to get labels of the total number of scored variants in each tile into the plot below
tot_vars <- as.data.frame(table(all_df$dist_group))
colnames(tot_vars) <- c('dist','total')

p6 <- ggplot(all_df, aes(x=dist_group,fill = class_no_na)) +
  geom_bar(position="fill", alpha = 0.7, color = 'grey50') +
  scale_fill_manual(values = class_palette_red, name = 'classification') +
  ggtitle('Fraction of variants sorted into') +
  ylab('fraction') +
  theme_bw(base_size = 16) +
  annotate("text", x = tot_vars$dist, y = 1.05,
           label = tot_vars$total)
p6

png(filename = paste0(plot_dir,'MAVE_categories_per_dist_fraction.png'),
    width = 800, height = 550, res = 100)
print(p6)
dev.off()                                  

#to get labels of the total number of scored variants in each tile into the plot below
#tot_vars <- as.data.frame(table(all_df$tile))
#colnames(tot_vars) <- c('tile','total')

p5 <- ggplot(all_df, aes(x=dist_group,fill = class_no_na)) +
  geom_bar(position="dodge", alpha = 0.7, color = 'grey50') +
  scale_fill_manual(values = class_palette_red, name = 'classification') +
  ggtitle('Variant classifications per distance group') +
  ylab('No. variants') +
  theme_bw(base_size = 16) 

p5

png(filename = paste0(plot_dir,'MAVE_categories_per_dist.png'),
    width = 800, height = 550, res = 100)
print(p5)
dev.off()  

#na as WT

#to get labels of the total number of scored variants in each tile into the plot below
tot_vars <- as.data.frame(table(all_df$dist_group))
colnames(tot_vars) <- c('dist','total')

p6 <- ggplot(all_df, aes(x=dist_group,fill = na_is_wt)) +
  geom_bar(position="fill", alpha = 0.7, color = 'grey50') +
  scale_fill_manual(values = class_palette_nona, name = 'classification') +
  ggtitle('Fraction of variants sorted into') +
  ylab('fraction') +
  theme_bw(base_size = 16) +
  annotate("text", x = tot_vars$dist, y = 1.05,
           label = tot_vars$total)
p6

png(filename = paste0(plot_dir,'MAVE_categories_per_dist_fraction_naiswt.png'),
    width = 800, height = 550, res = 100)
print(p6)
dev.off()                                  

#to get labels of the total number of scored variants in each tile into the plot below
#tot_vars <- as.data.frame(table(all_df$tile))
#colnames(tot_vars) <- c('tile','total')

p5 <- ggplot(all_df, aes(x=dist_group,fill = na_is_wt)) +
  geom_bar(position="dodge", alpha = 0.7, color = 'grey50') +
  scale_fill_manual(values = class_palette_nona, name = 'classification') +
  ggtitle('Variant classifications per distance group') +
  ylab('No. variants') +
  theme_bw(base_size = 16) 

p5

png(filename = paste0(plot_dir,'MAVE_categories_per_dist_naiswt.png'),
    width = 800, height = 550, res = 100)
print(p5)
dev.off()  



#3.2. ddG per classification group
########################################

#with unclassified
p <- ggplot(all_df, aes(ddG, fill = class_no_na, color = class_no_na)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = class_palette_red, name = 'classification') +
  scale_color_manual(values = class_palette_red, name = 'classification') +
  ggtitle('ddG per variant class') +
  theme_bw(base_size = 16)

p

png(filename = paste0(plot_dir,'ddG_per_cat.png'),
    width = 1200, height = 800, res = 100)
print(p)
dev.off() 

#unclassified is WT
p <- ggplot(all_df, aes(ddG, fill = na_is_wt, color = na_is_wt)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = class_palette_nona, name = 'classification') +
  scale_color_manual(values = class_palette_nano, name = 'classification') +
  ggtitle('ddG per variant class') +
  theme_bw(base_size = 16)

p

png(filename = paste0(plot_dir,'ddG_per_cat_naiswt.png'),
    width = 1200, height = 800, res = 100)
print(p)
dev.off() 

#3.3. ddE per classification group
########################################

#with unclassified
p <- ggplot(all_df, aes(ddE, fill = class_no_na, color = class_no_na)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = class_palette_red, name = 'classification') +
  scale_color_manual(values = class_palette_red, name = 'classification') +
  ggtitle('ddE per variant class') +
  theme_bw(base_size = 16)

p

png(filename = paste0(plot_dir,'ddE_per_cat.png'),
    width = 1200, height = 800, res = 100)
print(p)
dev.off() 

#unclassified is WT
p <- ggplot(all_df, aes(ddE, fill = na_is_wt, color = na_is_wt)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = class_palette_nona, name = 'classification') +
  scale_color_manual(values = class_palette_nona, name = 'classification') +
  ggtitle('ddE per variant class') +
  theme_bw(base_size = 16)

p

png(filename = paste0(plot_dir,'ddE_per_cat_naiswt.png'),
    width = 1200, height = 800, res = 100)
print(p)
dev.off() 

#3.4. RSA per classification group
########################################

#with unclassified
p <- ggplot(all_df, aes(RSA, fill = class_no_na, color = class_no_na)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = class_palette_red, name = 'classification') +
  scale_color_manual(values = class_palette_red, name = 'classification') +
  ggtitle('Accessible surface area per variant class') +
  theme_bw(base_size = 16)

p

png(filename = paste0(plot_dir,'RSA_per_cat.png'),
    width = 1200, height = 800, res = 100)
print(p)
dev.off() 

#unclassified is WT
p <- ggplot(all_df, aes(RSA, fill = na_is_wt, color = na_is_wt)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = class_palette_nona, name = 'classification') +
  scale_color_manual(values = class_palette_nona, name = 'classification') +
  ggtitle('Accessible surface area per variant class') +
  theme_bw(base_size = 16)

p

png(filename = paste0(plot_dir,'RSA_per_cat_naiswt.png'),
    width = 1200, height = 800, res = 100)
print(p)
dev.off() 
