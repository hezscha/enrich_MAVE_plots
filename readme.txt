
Scripts for the analysis of DHRF MAVE data.

Each R script contains below the imports a short description of what it is for.

Short desciptions:
format_prism_merge_file_to_df.R - read in the merged prism file and prep a dataframe with the same info but additional columns: sub, pos, WT
heatmap_full_DHFR_all_tiles_together.R - heatmap of all tiles together, plus distance to ligands and secondary structure elements
heatmap_full_DHFR_from_raw_data.R - plotting heatmaps directly from the csv files that come out of enrich2
heatmap_pilot_seq_tile1.R - heatmap of the pilot sequencing on tile 1 in 4 different conditions
MAVE_category_plots.R - some violine and density plots of i.e. ddE/ddG per variant category
MAVE_score_VS_new.R - plots of the MAVE score VS various measures, using the combined data and 3 plotting wrappers for continous, divergent and category coloring
MAVE_score_VS.R - the original MAVE score VS various measures for reference. Some plots I did not redo in the new version because I didn't think they made a lot of sense
pilot_seq_comp_rep_scores.R - comparing replicate scores from the pilot seq
position_averages.R - load in the merged prism file and create median, min and max per position for numerical columns
prep_MAVE_data.R - Merging MAVE data from all 5 tiles, averaging scores for vars scored in two tiles, adding vars only scored in 2 replicates, correcting p-values, calculating standard deviation of synonymous vars and categorizing vars according to p-value and scores with prep_MAVE_data.R. Produces a raw prism file
qual_check_cmp_overlap_var_scores.R - compare scores for vars that are scored in two tiles
score_density_plots.R - plotting the distributions of MAVE scores, both from the merged file and directly from the enrich output
tile4_non_mut_region.R - Tile 4 has a region sequenced that should not be mutagenized so all variants there should be false-positives. I made some analysis plots on that
vs_thompsen.R - comparing our MAVE scores to those in the Thompsen paper which did a DHFR MAVE in Ecoli
thompsen_data.R - make prism file from Ecoli csv data

R_tips.txt is a collection of commands and links I find useful for plotting.

The jupyter notebook replace_wt_w_sy_count_python2.ipynb shows how to open and manipulate hdf5 files, which is the format enrich uses to store counts/scores and performs calculations on. The tsv files are just for the user.

merge_dfs_human_ecoli.ipynb: read in tcoffee alignment and create merge dataframe where each row has the respective human and ecoli residues at that aligned position and their scores

fish_out_tile3_reads.sh is a workflow I used to subset reads from tile 3 to only the mutagenized region of interest (we had two regions sequenced for some reason) before running enrich2.

Instructions on how to run Erich2 and MAVE how to can be found in this shared google drive:
https://drive.google.com/drive/folders/1asocXbUFKnvw2Zy4_a5gIVR9do0h_T8b
Amelie, Sofie, Aleks and Henrike currently have access.

