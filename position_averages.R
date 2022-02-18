library(ggplot2)
library(stringr)
library(reshape2)
library(ggrepel)
library(dplyr) #arrange
library(hash)
library(matrixStats) #rowMedians

#load in the merged prism file and create averages per position for numerical columns
#perhaps once with only vars scored in all 3 reps and once with all vars (that have mave scores)

base_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/results/enrich_full_dhfr/'
data_dir <- '/home/henrike/Documents/PD_AS/projects/Sofie_Mave/data/'
plot_dir <- paste0(base_dir, '/HZ_plots/vs_plots/')

make_aa_hash_rv <- function(){
  aa_3_to_1 <- hash()
  aa_3_to_1[["Ala"]] <- 'A'
  aa_3_to_1[["Cys"]] <- 'C'
  aa_3_to_1[["Asp"]] <- 'D'
  aa_3_to_1[["Glu"]] <- 'E'
  aa_3_to_1[["Phe"]] <- 'F'
  aa_3_to_1[["Gly"]] <- 'G'
  aa_3_to_1[["His"]] <- 'H'
  aa_3_to_1[["Ile"]] <- 'I'
  aa_3_to_1[["Lys"]] <- 'K'
  aa_3_to_1[["Leu"]] <- 'L'
  aa_3_to_1[["Met"]] <- 'M'
  aa_3_to_1[["Asn"]] <- 'N'
  aa_3_to_1[["Pro"]] <- 'P'
  aa_3_to_1[["Gln"]] <- 'Q'
  aa_3_to_1[["Arg"]] <- 'R'
  aa_3_to_1[["Ser"]] <- 'S'
  aa_3_to_1[["Thr"]] <- 'T'
  aa_3_to_1[["Val"]] <- 'V'
  aa_3_to_1[["Trp"]] <- 'W'
  aa_3_to_1[["Tyr"]] <- 'Y'
  aa_3_to_1[["Ter"]] <- '*'
  return(aa_3_to_1)
}
aa_3_to_1 <- make_aa_hash_rv()

start_prot <- 1
end_prot <- 187

#load in and prep data
##########################################

all_df <- read.csv(paste0(data_dir, 'prism_merge_all_DHFR-human_P00374_28_01_2022.txt'), 
                   sep = ' ', comment.char = '#')
#add position and sub
all_df$pos <- as.numeric(str_extract(all_df$var, "([0-9]+)"))
all_df$sub <- str_match(all_df$var, "[0-9]+([A-Z*~=])")[,2]
all_df$WT <- str_match(all_df$var, "([A-Za-z]+)[0-9]+")[,2]


########################################
#two kinds of cols:
#numeric columns we need to average: gemme_score_00, norm_ddG_02, MAVE_Score_04, SE_04  

#cols who's original values are per position: MTXdist_03, NADPHdist_03, ASA_01
#how do I get these out in a clever way? Perhaps I can cast unique on the position, or just take the first instance like when I define the new row for vars that are scored in two tiles
#new_row <- all_df[all_df$pos == as.numeric(a[[1]][[1]]) &  all_df$sub == a[[1]][[2]], ][1,]

##########################################

#pos specific stats: we need all rows of the df
#########################################

#debug/trying out
#pos <- 4
#get all rows where pos == 4 and then subset to only the first one
#all_df[all_df$pos == pos, ][1,]

#init df: https://www.statology.org/create-empty-data-frame-in-r/
#we could actually omit the cols where we're gonna calc the average since we add these later when the df alsready has the right shape
avg <- data.frame(pos=integer(),
                  WT = character(),
                  MTXdist=double(),
                  NADPHdist=double(),
                  ASA=double(),
                  stringsAsFactors=FALSE)

#fill the df with the correct data for the position specific stats
for (pos in seq(start_prot,end_prot)) {
  #get all rows where pos == current_pos and then subset to only the first one
  pos_dat <- all_df[all_df$pos == pos, ][1,]
  #put the values for the position specific fields into the new df and leave the rest on NA, we'll calc them later
  new_row <- data.frame(pos=pos_dat$pos, WT=pos_dat$WT, MTXdist = pos_dat$MTXdist_03, NADPHdist = pos_dat$NADPHdist_03, 
                        ASA = pos_dat$ASA_01)
  avg <- rbind(avg, new_row)
}

#averaging var specific dat: actually we can also do that with the full df, just need to use rm.na = T when making the average
###########################################


#I may need to reomve the syn and ter vars 


##############################################
#MAVE score

mat <- dcast(all_df, pos ~ sub, value.var = "MAVE_score_04")
rownames(mat) <- mat$pos
#remove pos column
mat <- mat[,-c(1)]

pos_median = rowMedians(as.matrix(mat), na.rm = T)
pos_mean = rowMeans(as.matrix(mat), na.rm = T)
min_val = rowMins(as.matrix(mat), na.rm = T)
max_val = rowMaxs(as.matrix(mat), na.rm = T)

avg$med_MAVE <- pos_median
avg$min_MAVE <- min_val
avg$max_MAVE <- max_val

#MAVE error
#################
mat <- dcast(all_df, pos ~ sub, value.var = "SE_04")
rownames(mat) <- mat$pos
#remove pos column
mat <- mat[,-c(1)]

pos_median = rowMedians(as.matrix(mat), na.rm = T)
min_val = rowMins(as.matrix(mat), na.rm = T)
max_val = rowMaxs(as.matrix(mat), na.rm = T)

avg$med_MAVE_SE <- pos_median
avg$min_MAVE_SE <- min_val
avg$max_MAVE_SE <- max_val

######################################
#ddE

mat <- dcast(all_df, pos ~ sub, value.var = "gemme_score_00")
rownames(mat) <- mat$pos
#remove pos column
mat <- mat[,-c(1)]

pos_median = rowMedians(as.matrix(mat), na.rm = T)
pos_mean = rowMeans(as.matrix(mat), na.rm = T)
min_val = rowMins(as.matrix(mat), na.rm = T)
max_val = rowMaxs(as.matrix(mat), na.rm = T)

#avg$avg_ddE <- pos_mean
avg$med_ddE <- pos_median
avg$min_ddE <- min_val
avg$max_ddE <- max_val

#check how different median and mean are
# p <- ggplot(avg, aes(x=avg_ddE,y=med_ddE)) + 
#   geom_point()
# 
# p

######################################
#ddG

mat <- dcast(all_df, pos ~ sub, value.var = "norm_ddG_02")
rownames(mat) <- mat$pos
#remove pos column
mat <- mat[,-c(1)]

pos_median = rowMedians(as.matrix(mat), na.rm = T)
pos_mean = rowMeans(as.matrix(mat), na.rm = T)
min_val = rowMins(as.matrix(mat), na.rm = T)
max_val = rowMaxs(as.matrix(mat), na.rm = T)

#avg$avg_ddG <- pos_mean
avg$med_ddG <- pos_median
avg$min_ddG <- min_val
avg$max_ddG <- max_val

#check how different median and mean are
#p <- ggplot(avg, aes(x=avg_ddG,y=med_ddG)) + 
#   geom_point()
# 
# p

#ddg error
#################
mat <- dcast(all_df, pos ~ sub, value.var = "std_ddG_02")
rownames(mat) <- mat$pos
#remove pos column
mat <- mat[,-c(1)]

pos_median = rowMedians(as.matrix(mat), na.rm = T)
min_val = rowMins(as.matrix(mat), na.rm = T)
max_val = rowMaxs(as.matrix(mat), na.rm = T)

avg$med_ddG_std <- pos_median
avg$min_ddG_std <- min_val
avg$max_ddG_std <- max_val

#replace the inf and -inf caused by calc'ing mins and maxs with NA's
#################################
#https://statisticsglobe.com/replace-inf-with-na-in-vector-and-data-frame-in-r
avg <- do.call(data.frame, lapply(avg, function(x) replace(x, is.infinite(x), NA)))

#write out file
########################################

write.csv(avg, file = paste0(base_dir, 'combined_R_dataframes/', 'position_averages.csv'), quote = F, row.names = F)





