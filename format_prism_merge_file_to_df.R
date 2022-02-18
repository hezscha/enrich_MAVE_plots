library(ggplot2)
library(stringr)
library(reshape2)
library(matrixStats)
library(ggrepel)
library(dplyr) #arrange
library(hash)

#load in the big merged prism file, re-add pos and sub cols and change the column names, write out a csv

base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/results/enrich_full_dhfr/'
data_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/data/'

all_df <- read.csv(paste0(data_dir, 'prism_merge_all_DHFR-human_P00374_04_02_2022.txt'), 
                   sep = ' ', comment.char = '#', 
                   col.names = c('var', 'ddE', 'sse', 'RSA', 'chainID', 'ddG', 'ddG_SE',
                                 'MTXdist', 'NADPHdist', 'MAVE_SE', 'MAVE_score', 'tile',
                                 'p_raw', 'p_bonf', 'scored_in','classification'))

#add position and sub
all_df$pos <- as.numeric(str_extract(all_df$var, "([0-9]+)"))
all_df$sub <- str_match(all_df$var, "[0-9]+([A-Z*~=])")[,2]
all_df$WT <- str_match(all_df$var, "([A-Za-z]+)[0-9]+")[,2]

write.csv(all_df, paste0(data_dir, 'full_MAVE_data.csv'),
          quote = F, row.names = F)
