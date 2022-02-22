library(stringr) #str_extract
library(dplyr) #arrange

base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/'

thompsen1 <- read.csv(paste0(base_dir, 'data/elife-53476-supp1-v2.csv'), sep = ',', header = F, comment.char = '#')
colnames(thompsen1) <- c('var', 'Lon_minus_score', 'Lon_minus_std_dev', 'Lon_minus_std_err')
thompsen2 <- read.csv(paste0(base_dir, 'data/elife-53476-supp2-v2.csv'), sep = ',', header = F, comment.char = '#')
colnames(thompsen2) <- c('var', 'Lon_plus_score', 'Lon_plus_std_dev', 'Lon_plus_std_err')

thompsen <- merge(thompsen1, thompsen2, by = 'var', all = T)

#add pos col + sort
thompsen$pos <- as.numeric(str_extract(thompsen$var, "([0-9]+)"))
thompsen$WT <- str_match(thompsen$var, "([A-Z])[0-9]+")[,2]
#https://stackoverflow.com/questions/21003311/how-to-combine-multiple-character-columns-into-a-single-column-in-an-r-data-fram
thompsen$WT_pos <- with(thompsen, paste0(WT, pos))
#thompsen1$sub <- str_match(d$var, "[0-9]+([A-Za-z]+)")[,2] 
#order by position
thompsen <- arrange(thompsen, pos)
unique(thompsen$WT_pos)

write.table(thompsen[,c('var','Lon_minus_score', 'Lon_minus_std_dev', 'Lon_minus_std_err','Lon_plus_score', 'Lon_plus_std_dev', 'Lon_plus_std_err')], 
            file = paste0(base_dir, 'data/prism_mave_thompsen_minus_Lon.raw'),
          quote = F, sep = ' ', row.names = F)
