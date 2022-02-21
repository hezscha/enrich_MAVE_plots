library(ggplot2)
library(stringr)
library(reshape2)
library(ggrepel)
library(dplyr) #arrange
library(hash)
#add marginal distr to dde vs ddg plot
library(gridExtra)
library(egg)
library(gghalves) #more customizable jitter for rain clouds
library(ggdist) #halfeye for rain clouds

#plot average or indiv MAVE score by position against various things including ddsp accessible surface area

###################################################
#1. functions and global vars
##############################################
base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/results/enrich_full_dhfr/'
data_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/data/'
plot_dir <- paste0(base_dir, '/HZ_plots/vs_plots/')

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

start_prot <- 1
end_prot <- 187

#color palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

dist_palette <- c("both_close" = "#009E73", #green
                  "MTX_close" = "#E69F00",  #egg-yellow
                  "NADPH_close" = "#56B4E9", #light blue
                  "both_far" =  "#CC79A7") #pink

#class_palette <- c('del' = 'darkred', 'ben' = 'forestgreen', 'WT' = 'grey50')
class_palette <- c('del' = 'darkred', 
                   'low_del' = 'indianred2',
                   'low_ben' = 'darkolivegreen2',
                   'ben' = 'chartreuse4', 
                   'WT-like' = 'lightgoldenrodyellow')

class_palette_cb <- c('del' = "#D55E00", #orange-red 
                      'ben' = "#0072B2", #a dark blue
                      'WT-like' = 'lightgoldenrodyellow') #cream/off-white

class_palette_uncl <- c('del' = "#D55E00", #orange-red 
                      'ben' = "#0072B2", #a dark blue
                      'WT-like' = 'lightgoldenrodyellow', #cream/off-white
                      'unclassified' = 'azure4') 

#this scatter plot can be colored with a continous color scale (scale_color_gradient)
plot_scatter_color_cont <- function(dat, x, y, outpath, xlab, ylab, title,
                                    pt_size = '',
                                    color='', color_name='', col_low = '', col_high = '',
                                    print_labs = F, lab_string='', lab_score='',
                                    width = 1100, height = 800){
  
  #generate the aesthetic mapping for the plot
  #do we want a color mapping? If yes, pass the column to color by, if no pass nothing (default value is '')
  if(color == ''){
    p <- ggplot(dat, aes_string(x=x, y=y))
  } else {
    p <- ggplot(dat, aes_string(x=x, y=y, color = color))
    
    #do we want to define the colors?
    if(col_low != '') {
      
      #do we want to name the color scale?
      if(color_name != '') {
        p <- p + scale_color_gradient(name = color_name, low = col_low, high = col_high)
      } else {
        p <- p + scale_color_gradient(low = col_low, high = col_high)
      }
      
    } else{
      #if we do not want to change the low and high cols we can still name the scale
      if(color_name != '') {
        p <- p + scale_color_gradient(name = color_name)
      }
    }
    
  }
  
  #add a geom_point with the desired point size (or default)
  if(pt_size == ''){
    p <-  p + geom_point()
  } else {
    p <- p + geom_point(size = pt_size)
  }
  
  #add x and y labs
  p <- p + xlab(xlab) + ylab(ylab)
  #add title
  p <- p + ggtitle(title)
  #add a larger sized text
  p <- p  + theme_bw(base_size = 16)
  
  
  #do we want labels?
  #lab_string is which column to use for labels
  #lab_score is which column to use to decided whether we put on a label, i.e. the MAVE_score
  if(print_labs){
    #print(lab_string)
    #print('\n')
    labs <- ifelse(abs(dat[,c(lab_score)])>=1, lab_string,'')
    #print(labs)
    p <- p + geom_text_repel(aes(label=labs), color = 'black')
  }
  
  png(filename = outpath,width = width, height = height)
  print(p)
  dev.off()
  
}

#divergent scale isntead: use scale_color_gradient2. We can also set the mid point if desired
plot_scatter_color_div <- function(dat, x, y, outpath, xlab, ylab, title,
                                    pt_size = '',
                                    color='', color_name='', col_low = '', col_high = '', midpoint = '',
                                    print_labs = F, lab_string='', lab_score='',
                                    width = 1100, height = 800){
  
  #calc midpoint unless supplied
  if(midpoint == '' & color != '') {
    midpoint <- mean(dat[[color]], na.rm = T)
    print(paste0('midpoint:', midpoint))
  }
  
  #generate the aesthetic mapping for the plot
  #do we want a color mapping? If yes, pass the column to color by, if no pass nothing (default value is '')
  if(color == ''){
    p <- ggplot(dat, aes_string(x=x, y=y))
  } else {
    p <- ggplot(dat, aes_string(x=x, y=y, color = color))
    
    #do we want to define the colors?
    if(col_low != '') {
      
      #do we want to name the color scale?
      if(color_name != '') {
        p <- p + scale_color_gradient2(name = color_name, low = col_low, high = col_high, midpoint = midpoint)
      } else {
        p <- p + scale_color_gradient2(low = col_low, high = col_high, midpoint = midpoint)
      }
      
    } else{
      #if we do not want to change the low and high cols we can still name the scale
      if(color_name != '') {
        p <- p + scale_color_gradient2(name = color_name, midpoint = midpoint)
      }
    }
    
  }
  
  #add a geom_point with the desired point size (or default)
  if(pt_size == ''){
    p <-  p + geom_point()
  } else {
    p <- p + geom_point(size = pt_size)
  }
  
  #add x and y labs
  p <- p + xlab(xlab) + ylab(ylab)
  #add title
  p <- p + ggtitle(title)
  #add a larger sized text
  p <- p  + theme_bw(base_size = 16)
  
  
  #do we want labels?
  #lab_string is which column to use for labels
  #lab_score is which column to use to decided whether we put on a label, i.e. the MAVE_score
  if(print_labs){
    #print(lab_string)
    #print('\n')
    labs <- ifelse(abs(dat[,c(lab_score)])>=1, lab_string,'')
    #print(labs)
    p <- p + geom_text_repel(aes(label=labs), color = 'black')
  }
  
  png(filename = outpath,width = width, height = height)
  print(p)
  dev.off()
  
}

#This scatter plot is colored with categorical data
plot_scatter_color_cat <- function(dat, x, y, outpath, xlab, ylab, title, 
                                   color, color_palette, color_name='', pt_size = '', shape = '',
                                   print_labs = F, lab_string='', lab_score='', legend_bot = F,
                                   width = 1100, height = 800){
  
  #generate the aesthetic mapping for the plot
  #do we want a shape mapping?
  if(shape == '') {
    p <- ggplot(dat, aes_string(x=x, y=y, color = color))
  } else {
    #add the shape
    p <- ggplot(dat, aes_string(x=x, y=y, color = color, shape = shape))
  }
  
  #do we want to name the color scale?
  if(color_name != '') {
    #add the palette
    p <- p + scale_color_manual(values = color_palette, name = color_name)
  } else {
    #add the palette
    p <- p + scale_color_manual(values = color_palette)
  }
  
  #add a geom_point with the desired point size (or default)
  if(pt_size == ''){
    p <-  p + geom_point()
  } else {
    p <- p + geom_point(size = pt_size)
  }
  
  #add x and y labs
  p <- p + xlab(xlab) + ylab(ylab)
  #add title
  p <- p + ggtitle(title)
  #add a larger sized text
  p <- p  + theme_bw(base_size = 16)
  
  
  #do we want labels?
  #lab_string is which column to use for labels
  #lab_score is which column to use to decided whether we put on a label, i.e. the MAVE_score
  if(print_labs){
    #print(lab_string)
    #print('\n')
    labs <- ifelse(abs(dat[,c(lab_score)])>=1, lab_string,'')
    #print(labs)
    p <- p + geom_text_repel(aes(label=labs), color = 'black')
  }
  
  if(legend_bot){
    p <- p + theme(legend.position = 'bottom')
  }
  
  png(filename = outpath,width = width, height = height)
  print(p)
  dev.off()
  
  return(p)
  
}

# #usage examples
# #trying default color scale but named
# plot_scatter_color_cont(dat = avg_df, x = 'med_MAVE', y = 'MTXdist', xlab = 'Median MAVE score', ylab ='Distance to MTX',
#                         title = 'Distance to the substrate analog MTX',
#                         color = 'med_ddE', color_name = 'Median\nddE',
#                         print_labs = T, lab_string = paste0(avg_df[,c('WT')],avg_df[,c('pos')]), lab_score = 'med_MAVE',
#                         outpath = paste0(plot_dir, folder_path, 'pos_VS_MTXdist_color_dde_scatter.png'))
# 
# #custom color scale
# plot_scatter_color_cont(dat = avg_df, x = 'med_MAVE', y = 'MTXdist', xlab = 'Median MAVE score', ylab ='Distance to MTX',
#                         title = 'Distance to the substrate analog MTX', pt_size = 3,
#                         color = 'med_ddE', color_name = 'Median\nddE', col_low = 'darkred', col_high = 'gold',
#                         print_labs = T, lab_string = paste0(avg_df[,c('WT')],avg_df[,c('pos')]), lab_score = 'med_MAVE',
#                         outpath = paste0(plot_dir, folder_path, 'pos_VS_MTXdist_color_dde_scatter2.png'))

#####################################################################
#the plot I deconstructed into the plot_scatter_color_cont() function

# p <- ggplot(avg_df, aes(x=med_MAVE, y=MTXdist, color = med_ddE)) +
#   geom_point(size=3) +
#   scale_color_gradient(name = 'Median\nddE', low = "darkred", high = "gold") +
#   ggtitle('Distance to the substrate analog MTX') +
#   xlab('Median MAVE score') + ylab('Distance to MTX') +
#   geom_text_repel(aes(label=ifelse(abs(med_MAVE)>=1, paste0(WT,pos),'')), color = 'black') +
#   theme_bw(base_size = 16)
# 
# png(filename = paste0(plot_dir, folder_path, 'pos_VS_MTXdist_color_dde.png'),width = 1200, height = 800)
# print(p)
# dev.off()

###################################################
#2. load data
###################################################

all_df_full <- read.csv(paste0(data_dir, 'full_MAVE_data.csv'))

#subset to only rows that contain MAVE scores. There should be 3106 rows left 
all_df <- all_df_full[!is.na(all_df_full$MAVE_score),]

#make another cloumn where the NA class has a real name so it doesn't drop and is displayed in plots
all_df$class_no_na <- ifelse(is.na(all_df$classification), 'unclassified', all_df$classification)

#subset to omit vars in the low_ben and low_del categories. These we expect to functionally behave like WT
#this also omits vars where the classification is NA (this is a characteristic of subset. It doesn't happen when using [] to subset the df)
#but when using [], it sets the entire row to NA???
#all_df_sub1 <- all_df[(!all_df$classification == 'low_ben' & !all_df$classification == 'low_del'), ]
all_df <- subset(all_df, (!class_no_na == 'low_ben' & !class_no_na == 'low_del'))

#add distance to ligands grouping
all_df$dist_group <- ifelse((all_df$MTXdist <= 10 & all_df$NADPHdist <= 10), 'both_close', 
                            ifelse((all_df$MTXdist <= 10 & all_df$NADPHdist > 10), 'MTX_close', 
                                   ifelse((all_df$MTXdist > 10 & all_df$NADPHdist <= 10), 'NADPH_close',
                                   'both_far')))

#position medians, min, max
#################################

#this is on the full data, so all measured variants are included into the median of their position
#regardless of their i.e. variant class (benefical, unclassified, low benefical, ect).

avg_df <- read.csv(paste0(base_dir, 'combined_R_dataframes/', 'position_averages.csv'))
#add distance to ligands grouping
avg_df$dist_group <- ifelse((avg_df$MTXdist <= 10 & avg_df$NADPHdist <= 10), 'both_close', 
                            ifelse((avg_df$MTXdist <= 10 & avg_df$NADPHdist > 10), 'MTX_close', 
                                   ifelse((avg_df$MTXdist > 10 & avg_df$NADPHdist <= 10), 'NADPH_close',
                                          'both_far')))

###################################################
#3. scatter plots
###################################################

###################################################
#3.1 position specific plots
###################################################
#the dataframe is the averaged data and I will usually use the WT + pos as labels
dat <- avg_df
lab_string = paste0(avg_df[,c('WT')],avg_df[,c('pos')])
set <- 'pos'

#3.1.1 avrg mave score per position VS distance to ligands
########################################
folder_path = '/dist/by_pos/'
width = 800
height = 550

#set up for plots in this section: They are all MAVE score VS distance
x = 'med_MAVE'
y = 'MTXdist'
xlab = 'Median MAVE score'
ylab ='Distance to MTX'
title = 'Distance to the substrate analog MTX, per position'
lab_score = 'med_MAVE'

#3.1.1.1 color by median ddE
########################################
plot_scatter_color_cont(dat, x, y, xlab, ylab, title, pt_size = 3,
                        color = 'med_ddE', color_name = 'Median\nddE', col_low = 'darkred', col_high = 'gold',
                        print_labs = T, lab_string, lab_score, width = width, height = height,
                        outpath = paste0(plot_dir, folder_path, set, '_VS_MTXdist_color_dde.png'))

#3.1.1.2 color by min ddE
########################################
#color by min ddE, i.e. the most extreme exclusion of a variant (lower ddE means the pos is more conserved, i.e. high penalty to change it)
plot_scatter_color_cont(dat, x, y, xlab, ylab, title, pt_size = 3,
                        color = 'min_ddE', color_name = 'Min\nddE', col_low = 'darkred', col_high = 'gold',
                        print_labs = T, lab_string, lab_score, width = width, height = height,
                        outpath = paste0(plot_dir, folder_path, set, '_VS_MTXdist_color_dde_min.png'))


#3.1.1.3 color by max ddG
########################################
plot_scatter_color_cont(dat, x, y, xlab, ylab, title, pt_size = 3,
                        color = 'max_ddG', color_name = 'Max\nddG', col_low = 'gold', col_high = 'darkred',
                        print_labs = T, lab_string, lab_score, width = width, height = height,
                        outpath = paste0(plot_dir, folder_path, set, '_VS_MTXdist_color_ddg_max.png'))

#3.1.1.3 color by median ddG
########################################
plot_scatter_color_cont(dat, x, y, xlab, ylab, title, pt_size = 3,
                        color = 'med_ddG', color_name = 'Median\nddG', col_low = 'gold', col_high = 'darkred',
                        print_labs = T, lab_string, lab_score, width = width, height = height,
                        outpath = paste0(plot_dir, folder_path, set, '_VS_MTXdist_color_ddg_med.png'))


#3.1.5 ddG Vs MAVE
########################################
folder_path = '/ddG_VS_MAVE/by_pos/'

x = 'med_MAVE'
y = 'med_ddG'
xlab = 'median MAVE score'
ylab ='median ddG score'
title = 'ddG Vs med MAVE, per pos'


#3.2.4.2 colored by distance to ligands group
########################################
color <- 'dist_group'
plot_scatter_color_cat(dat, x, y, xlab, ylab, title = title, pt_size = 3,
                       color = color, color_palette = dist_palette, legend_bot = T,
                       #print_labs = T, lab_string = lab_string, lab_score = lab_score,
                       outpath = paste0(plot_dir, folder_path, set, '_ddG_VS_MAVE_color_dist.png'))


#3.1.6 RSA VS med MAVE
########################################
folder_path = '/RSA_VS_MAVE/by_pos/'

x = 'med_MAVE'
y = 'ASA'
xlab = 'median MAVE score'
ylab ='relative accessible surface area'
title = 'RSA VS Median MAVE'

set <- 'pos'
#color = 'classification'

plot_scatter_color_cont(dat, x, y, xlab, ylab, title = title, pt_size = 3,
                       #color = color, color_palette = dist_palette, legend_bot = T,
                       #print_labs = T, lab_string = lab_string, lab_score = lab_score,
                       outpath = paste0(plot_dir, folder_path, set, '_RSA_VS_MAVE.png'))


###################################################
#3.2 variant specific plots
###################################################
dat <- all_df
lab_string = all_df$var

#3.2.1 mave score VS distance to ligands
########################################
folder_path = '/dist/by_var/'

x = 'MAVE_score'
y = 'MTXdist'
xlab = 'MAVE score'
ylab ='Distance to MTX'
lab_score = 'MAVE_score'



#3.2.1.1 colored by ddE
########################################

#all vars
set <- 'all'
#color by ddE
plot_scatter_color_cont(dat, x, y, xlab, ylab, title = 'Distance to the substrate analog MTX, all vars', pt_size = 3,
                        color = 'ddE', col_low = 'darkred', col_high = 'gold',
                        #print_labs = T, lab_string = lab_string, lab_score = lab_score,
                        outpath = paste0(plot_dir, folder_path, set, '_VS_MTXdist_color_dde.png'))

set <- 'all'
#color by ddE
plot_scatter_color_cont(dat, x, y, xlab, ylab, title = 'Distance to the substrate analog MTX, all vars', pt_size = 3,
                        color = 'ddG', col_low = 'gold', col_high = 'darkred',
                        #print_labs = T, lab_string = lab_string, lab_score = lab_score,
                        outpath = paste0(plot_dir, folder_path, set, '_VS_MTXdist_color_ddG.png'))


#neutral vars (WT like)
#set <- 'neutral'
#color by ddE
#plot_scatter_color_cont(neutral, x, y, xlab, ylab, title = 'Distance to the substrate analog MTX, WT-like vars', pt_size = 3,
#                        color = 'ddE', col_low = 'darkred', col_high = 'gold',
#                        outpath = paste0(plot_dir, folder_path, set, '_VS_MTXdist_color_dde.png'))
# 
# set <- 'beneficial'
# #color by ddE
# plot_scatter_color_cont(good, x, y, xlab, ylab, title = 'Distance to the substrate analog MTX, beneficial vars', pt_size = 2,
#                         color = 'ddE', col_low = 'darkred', col_high = 'gold',
#                         outpath = paste0(plot_dir, folder_path, set, '_VS_MTXdist_color_dde.png'))
# 
# set <- 'deleterious'
# #color by ddE
# plot_scatter_color_cont(bad, x, y, xlab, ylab, title = 'Distance to the substrate analog MTX, deleterious vars', pt_size = 2,
#                         color = 'ddE', col_low = 'darkred', col_high = 'gold',
#                         outpath = paste0(plot_dir, folder_path, set, '_VS_MTXdist_color_dde.png'))


#3.2.2 MTXdist Vs NADPHdist
########################################
folder_path = '/dist/by_var/'

x = 'NADPHdist'
y = 'MTXdist'
xlab = 'Distance to NADPH'
ylab ='Distance to MTX'
color = 'classification'

p2 <- ggplot(all_df, aes_string(x=x,y=y, color = color)) +
  geom_point(size = 1) +
  xlab(xlab) + ylab(ylab)+
  scale_color_manual(values = class_palette_cb)

p2

#3.2.3 ddE Vs ddG
########################################
folder_path = '/ddE_VS_ddG/by_var/'

x = 'ddE'
y = 'ddG'
xlab = 'gemme score'
ylab ='ddG score'
title = 'ddE Vs ddG, all vars'


#3.2.3.1 colored by var classification
########################################
color = 'classification'

set <- 'all'
p1 <- plot_scatter_color_cat(all_df, x, y, xlab, ylab, title = title, pt_size = 3,
                        color = 'classification', color_palette = class_palette_cb,
                        legend_bot = T,
                        #print_labs = T, lab_string = lab_string, lab_score = lab_score,
                        outpath = paste0(plot_dir, folder_path, set, 'ddE_VS_ddG_color_class.png'))
p1



#3.2.3.1.1 plus marginals
########################################

#add distributions per class to both ddE and ddG
distr_y <- ggplot(all_df, aes(ddG, color = classification)) +
  geom_density(bw= 0.08) +
  scale_color_manual(values = class_palette_cb) + 
  #xlim(-3,20) +
  #ggtitle(paste0("tile ", tile, ", ddG values")) +
  #xlab('Bandwidth = 0.08') +
  coord_flip() +
  theme_bw(base_size = 20)  +
  theme(axis.title.y=element_blank(), legend.position = 'none')

distr_y

distr_x <- ggplot(all_df, aes(ddE,color = classification)) +
  geom_density(bw= 0.08) +
  #ggtitle(paste0('One point per variant: ddG VS ddE scores colored by MAVE score\nddG values capped at ', cutoff_ddG)) +
  scale_color_manual(values = class_palette_cb) + 
  #xlim(-3,20) +
  #ggtitle(paste0("tile ", tile, ", ddG values")) +
  #xlab('Bandwidth = 0.08') +
  theme_bw(base_size = 20) +
  theme(axis.title.x=element_blank(), legend.position = 'none')

distr_x

#little blank square in the top right corner
gg_empty <- ggplot(all_df, aes(x=ddE,y=ddG))+
  geom_blank() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        panel.background = element_blank())

png(filename = paste0(plot_dir, folder_path, set, 'ddE_VS_ddG_color_class_marginals.png'),width = 1200, height = 900)
ggarrange(
  distr_x, gg_empty, p1, distr_y,
  nrow = 2, ncol = 2, widths = c(5, 1), heights = c(1, 5)
)
dev.off()
#########################################


#3.2.3.2 colored by distance to ligands group
########################################
set <- 'all'
color <- 'dist_group'
plot_scatter_color_cat(all_df, x, y, xlab, ylab, title = title, pt_size = 3,
                       color = color, color_palette = dist_palette,
                       #print_labs = T, lab_string = lab_string, lab_score = lab_score,
                       outpath = paste0(plot_dir, folder_path, set, 'ddE_VS_ddG_color_dist.png'))

set <- 'all'
plot_scatter_color_cat(all_df, x, y, xlab, ylab, title = title, pt_size = 3,
                       shape = 'classification',
                       color = 'dist_group', color_palette = dist_palette,
                       #print_labs = T, lab_string = lab_string, lab_score = lab_score,
                       outpath = paste0(plot_dir, folder_path, set, 'ddE_VS_ddG_color_dist_shape_class.png'))


#3.2.4 ddE Vs MAVE
########################################
folder_path = '/ddE_VS_MAVE/by_var/'

x = 'MAVE_score'
y = 'ddE'
xlab = 'MAVE score'
ylab ='gemme score'
title = 'ddE Vs MAVE, all vars'
set <- 'all'

#3.2.4.2 colored by distance to ligands group
########################################
color <- 'dist_group'
plot_scatter_color_cat(all_df, x, y, xlab, ylab, title = title, pt_size = 3,
                       color = color, color_palette = dist_palette, legend_bot = T,
                       #print_labs = T, lab_string = lab_string, lab_score = lab_score,
                       outpath = paste0(plot_dir, folder_path, set, '_ddE_VS_MAVE_color_dist.png'))


#3.2.5 ddG Vs MAVE
########################################
folder_path = '/ddG_VS_MAVE/by_var/'

x = 'MAVE_score'
y = 'ddG'
xlab = 'MAVE score'
ylab ='ddG score'
title = 'ddG Vs MAVE, all vars'
set <- 'all'

#3.2.5.2 colored by distance to ligands group
########################################
color <- 'dist_group'
plot_scatter_color_cat(all_df, x, y, xlab, ylab, title = title, pt_size = 2,
                       color = color, color_palette = dist_palette, legend_bot = T,
                       #print_labs = T, lab_string = lab_string, lab_score = lab_score,
                       outpath = paste0(plot_dir, folder_path, set, '_ddG_VS_MAVE_color_dist.png'))


#not a great idea because all the vars at the same pos are gonna be on top of each other so you can't see them properly
#could work with median mave score instead but then how do we subset to only del vars? 
# #3.2.6 NADdist vs MTXdist
# ########################################
# folder_path = '/dist/by_pos/'
# 
# x = 'MTXdist'
# y = 'NADPHdist'
# xlab = 'Distance to MTX'
# ylab ='Distance to NADPH'
# title = 'MTXdist VS NADPHdist, only del vars'
# 
# set <- 'del'
# del_set <- subset(all_df, classification == 'del')
# 
# color <- 'MAVE_score'
# plot_scatter_color_div(del_set, x, y, xlab, ylab, title = title, pt_size = 3,
#                        color = color, color_name = 'MAVE_score', midpoint = 0,
#                        #print_labs = T, lab_string = lab_string, lab_score = lab_score,
#                        outpath = paste0(plot_dir, folder_path, set, '_MTXdist_VS_NADdist_color_MAVE.png'))


######################################
#4. box plots
######################################

#4.1 distance to ligands per variant class
######################################
folder_path = '/dist/boxplot/'
set <- 'all'
width = 800
height = 550


#MTXdist
p <- ggplot(all_df,aes(x= class_no_na, y=MTXdist, fill = class_no_na)) +
  geom_violin(alpha = 0.7) +
  geom_jitter(width = 0.2, color = 'grey50', alpha = 0.9) +
  scale_fill_manual(values = class_palette_uncl, name = 'classification') +
  theme_bw(base_size = 16) +
  ggtitle('Distance to MTX by variant class')

png(filename = paste0(plot_dir, folder_path, set, '_MTXdist_by_var_class.png'),width = width, height = height)
print(p)
dev.off()

#NADPHdist
p <- ggplot(all_df,aes(x= class_no_na, y=NADPHdist, fill = class_no_na)) +
  geom_violin(alpha = 0.7) +
  geom_jitter(width = 0.2, color = 'grey50', alpha = 0.9) +
  scale_fill_manual(values = class_palette_uncl, name = 'classification') +
  theme_bw(base_size = 16) +
  ggtitle('Distance to NADPH by variant class')
p

png(filename = paste0(plot_dir, folder_path, set, '_NADdist_by_var_class.png'),width = width, height = height)
print(p)
dev.off()

#4.2. ddG of only del vars per distance group: are del vars far from the binding site more destabilizing?
######################################
set <- 'del'
del_set <- subset(all_df, classification == 'del')

p <- ggplot(del_set, aes(x= dist_group, y=ddG, fill = dist_group)) +
  geom_violin(alpha = 0.7) +
  geom_jitter(width = 0.2, color = 'grey50', alpha = 0.7) +
  scale_fill_manual(values = dist_palette) +
  theme_bw(base_size = 16) +
  ggtitle('ddG of deleterious vars by distance from the ligands')

p
png(filename = paste0(plot_dir, folder_path ,set, '_ddG_by_dist.png'),
    width = 800, height = 550, res = 100)
print(p)
dev.off()

#4.2. ddE of only del vars per distance group
######################################
set <- 'del'
del_set <- subset(all_df, classification == 'del')

p <- ggplot(del_set, aes(x= dist_group, y=ddE, fill = dist_group)) +
  geom_violin(alpha = 0.7) +
  geom_jitter(width = 0.2, color = 'grey50', alpha = 0.7) +
  scale_fill_manual(values = dist_palette) +
  theme_bw(base_size = 16) +
  ggtitle('ddE of deleterious vars by distance from the ligands')

p
png(filename = paste0(plot_dir, folder_path, set,'_ddE_by_dist.png'),
    width = 800, height = 550, res = 100)
print(p)
dev.off() 

######################################
#5. bar plots 
######################################

#5.1 
#how do the classification classes interact with the distance groups, i.e. both_close, only MTX close, ect
#may a stacked bar of classification per dist group

#5.1.1. distance group per var classification
#################################
#to get labels of the total number of scored variants in each tile into the plot below
tot_vars <- as.data.frame(table(all_df$dist_group))
colnames(tot_vars) <- c('dist_group','total')

p4 <- ggplot(all_df, aes(x=dist_group,fill = classification)) +
  geom_bar(position="fill", alpha = 0.7, color = 'grey50') +
  scale_fill_manual(values = class_palette_cb) +
  ggtitle('Fraction of variants sorted into') +
  ylab('fraction') +
  theme_bw(base_size = 16) +
  annotate("text", x = tot_vars$dist_group, y = 1.05, label = tot_vars$total)

p4

png(filename = paste0(plot_dir, '/dist/by_var/','MAVE_categories_per_dist_group.png'),
    width = 800, height = 550, res = 100)
print(p4)
dev.off()  

#5.1.2. var classification per distance group  
#################################
#to get labels of the total number of scored variants in each tile into the plot below
tot_vars <- as.data.frame(table(all_df$classification))
colnames(tot_vars) <- c('classification','total')

p5 <- ggplot(all_df, aes(x=classification,fill = dist_group)) +
  geom_bar(position="fill", alpha = 0.7, color = 'grey50') +
  scale_fill_manual(values = dist_palette) +
  ggtitle('Fraction of variants sorted into') +
  ylab('fraction') +
  theme_bw(base_size = 16) +
  annotate("text", x = tot_vars$classification, y = 1.05, label = tot_vars$total)

p5

png(filename = paste0(plot_dir, '/dist/by_var/','MAVE_dists_per_cat.png'),
    width = 800, height = 550, res = 100)
print(p5)
dev.off()  

###########################
 



############################################################################
###########################
#6. plots per variant classification
###########################
