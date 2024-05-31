

# library(data.table)
library(readxl)
# library(stringr)
# library(vegan)
library(picante)
library(igraph)
library(psych)
library(adespatial)
# library(betapart)
library(tidyverse)
# library(tidypaleo)

inputpath  <- "//smb.isipd.dmawi.de/projects/biodiv/user/metabarcoding_plants/cores/12_100id_filtered/"

source(paste0(inputpath,"box.cox.chord.R"))



save_output    <- "no"
corr_threshold <- 0.3

# dir_out <- paste0(inputpath,"Community_investigation/corrthres_",corr_threshold)
# if(!file.exists(dir_out)) {
#   dir.create(dir_out)
#   dir.create(paste0(dir_out,"/Figures"))
#   fs::dir_tree(path = dir_out, recurse = TRUE)
# }


dir_out  <- "C:/Users/roschw001/Documents/R/SDM/SDM" #C:/Data/Temp_CommunityAnalysis/"


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

### load data:

gh_datasets       <- read_xlsx("//smb.isipd.dmawi.de/projects/biodiv/user/metabarcoding_plants/cores/overview_cores_gh datasets.xlsx", sheet="overview", col_names=TRUE, trim_ws=TRUE, col_types = c("numeric", "text", "text", "numeric", "numeric", "text", "text"))

tab_seq_taxanames <- read.csv(paste0(inputpath,"Bryophytes/combined_gh_datasets_table_sequences_taxanames_eco_revised.csv"), header=TRUE)
timeslicetable    <- read.csv(paste0(inputpath,"combined_metabarcoding_dataset_timeslice_table_2000years.csv"), header=TRUE) 

df_timeslice_raw  <- read.csv(paste0(inputpath,"combined_metabarcoding_dataset_sequences_filtered_3samples_2cores_100reads_timeslices_2000years_longformat.csv"), header=TRUE) 

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

### search for aquatic taxa and sqrt the counts; terrestrial taxa stay as original counts

{
  
  aquatics <- tab_seq_taxanames[tab_seq_taxanames$Aquatic == 1, ]
  
  df_timeslice_aquatic <- df_timeslice_raw %>% 
    filter(ID %in% unique(aquatics$Seq_number)) %>% 
    mutate(counts = sqrt(sqrt(counts)))
  
  df_timeslice_terrestrial <- df_timeslice_raw %>% 
    filter(!(ID %in% unique(aquatics$Seq_number)))
  
  df_timeslice_aqua_terra <- rbind(df_timeslice_terrestrial, df_timeslice_aquatic)
  
}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

### Network analysis

{
  
  overview_taxanumbers_age <- df_timeslice_aqua_terra %>%
    mutate(short_label = paste(ID, scientific_name, sep = "_")) %>%
    select(core_time, core, lake_name, Longitude, Latitude, country_region, short_label, counts, age, time_slice) %>%
    group_by(core, age) %>% 
    filter(counts > 0) %>% 
    mutate(taxanum = length(unique(short_label))) %>% 
    ungroup() %>% 
    select(core_time, taxanum) %>% 
    distinct()
  
  overview_taxanumbers_age <- overview_taxanumbers_age %>% filter(taxanum > 20)
  
  
  df_taxa_age <- df_timeslice_aqua_terra %>% 
    filter(core_time %in% overview_taxanumbers_age$core_time)
  
  colnames(df_taxa_age) <- c("core_time", "core", "lake_name", "Longitude", "Latitude", "country_region", "NUC_SEQ", "identity", "scientific_name", "family", "genus", "species", "sample_name", "Extraction_number", "counts", "depth", "age", "SeqID", "time_slice")
  
  
  df_raw_info <- df_taxa_age %>% select(NUC_SEQ, SeqID, identity, family, genus, species, scientific_name) %>% distinct()
  
  # -------------------------------------------------------------------------------------------------
  
  # Measure some metrics: total reads count per raxa, number of sample occurence, number of cores occurence and total number of reads per sample
  df_raw1 <- df_taxa_age %>% 
    mutate(sample_id = paste(core, age, sep = "_")) %>% 
    group_by(SeqID) %>% mutate(reads_seq = sum(counts), nsamples = n_distinct(sample_id), ncores = n_distinct(core)) %>% ungroup()
  
  df_raw1_clean <- df_raw1 %>% arrange(core, age) %>% select(SeqID, sample_id, counts, nsamples, ncores, reads_seq) %>% distinct() 
  
  df_raw_wide <- df_raw1_clean %>% pivot_wider(names_from = sample_id, values_from = counts) #set as wide dataset
  df_raw_wide[is.na(df_raw_wide)] <- 0 # replace NAs per 0
  df_raw_wide <- df_raw_wide %>% mutate(id = 1:nrow(df_raw_wide)) # set an id for each row of the data
  
  df_raw_wide1 <- left_join(df_raw_info, df_raw_wide, by = "SeqID")
  df_raw_wide1 <- df_raw_wide1 %>% mutate(uniq_id = paste(SeqID, scientific_name, identity, sep = "_")) # set a uniq taxid for each ASV
  
  # get the total number of reads we handle
  sum(df_raw_wide1$reads_seq)
  
  # add uniq - short labels
  df <- df_raw_wide1 %>% mutate(short_label = paste(SeqID, scientific_name, sep = "_"))
  
  ############################################################################
  # 1 - Set parameters and save general information
  ############################################################################
  
  # create metadata taxa
  metadata_taxa <- df %>% select(uniq_id, short_label, SeqID, NUC_SEQ, identity, family, genus, species, scientific_name, nsamples, ncores, reads_seq, id)
  
  # subset for minimum 100 reads
  df_100 <- df %>% arrange(reads_seq) %>%  ##%>% filter(reads_seq >= 100)
    select(-id, -uniq_id, -SeqID, -NUC_SEQ, -identity, -family, -genus, -species, -scientific_name, -nsamples, -ncores, -reads_seq) 
  
  # Put taxa_id as row_names
  df_rdy_com <- as.data.frame(df_100)
  rownames(df_rdy_com) <- df_rdy_com$short_label
  df_rdy_com$short_label <- NULL
  
  ############################################################################
  # 2 - Co-occurence and community building -> TCorr.test and iGraph
  ############################################################################
  
  # box.cox.chord transform the counts:
  df_trans <- box.cox.chord(df_rdy_com, bc.exp=0)
  
  # -------------------------------------------------------------------------------------------------
  # save output table:
  # if(save_output == "yes"){save(df_trans, file = paste0(dir_out,"/df_for_co-occurence_test.rda"))}
  # -------------------------------------------------------------------------------------------------
  
  # Do correlation of occurence matrix - very time consuming
  corr_df <- corr.test(t(df_trans),       # a matrix or dataframe
                       use = "pairwise",  # pairwise" is the default value and will do pairwise deletion of cases. use="complete" will select just complete cases
                       method="spearman", # default is Pearson. spearman is slower especially for large datasets
                       adjust="holm",     # What adjustment for multiple tests should be used? ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none").
                       alpha=.05,         # alpha level of confidence intervals
                       ci=FALSE,          # default is TRUE but set false. For very large matrices (>200x200), noticeable speed improvement if confidence interval are not found
                       minlength=5,       # What is the minimum length for abbreviations. Defaults to 5.
                       normal=F)          # normal = FALSE -> we do not have normal distribution. But it is much slower to process
  
  # -------------------------------------------------------------------------------------------------
  # save output table:
  if(save_output == "yes"){save(corr_df, file = paste0(dir_out,"/Co-occurence_raw_results.rda"))}
  # -------------------------------------------------------------------------------------------------
  
  ############################################################################
  # 3 - Prepare for network analysis
  ############################################################################
  
  # Make p-values as long
  corr_p      <- as.data.frame(corr_df$p)
  corr_p$sp1  <- rownames(corr_p)
  corr_p_long <- gather(corr_p, sp2, prob, all_of(colnames(as.data.frame(corr_df$p)))) 
  
  # Make uniq_ids
  corr_p_long$uniqid <- paste(corr_p_long$sp1, corr_p_long$sp2, sep = "_")
  corr_p_long        <- corr_p_long %>% select(uniqid, prob)
  
  # Re-format correlation matrix and merge p-value
  corr_r      <- as.data.frame(corr_df$r)
  corr_r$sp1  <- rownames(corr_r)
  corr_r_long <- gather(corr_r, sp2, corr, all_of(colnames(as.data.frame(corr_df$r))))
  
  # Make uniq_ids
  corr_r_long$uniqid <- paste(corr_r_long$sp1, corr_r_long$sp2, sep = "_")
  
  # join both dataframes:
  corr_net <- left_join(corr_r_long, corr_p_long, by = "uniqid") %>% 
    as_tibble() %>% 
    select(sp1, sp2, corr, prob)
  
  # remove corr = 1 (correlation to same asv)
  corr_net <- corr_net %>% subset(!corr == 1) 
  
  # save nodes information
  nodes <- data.frame(label = rownames(df_trans))
  nodes <- left_join(nodes, metadata_taxa, by = c("label" = "short_label")) %>% as_tibble()
  
  ############################################################################
  # 4 - Keep only the positive and sup to corr_score correlation scores
  ############################################################################
  
  corr_score   <- corr_threshold # minimum correlation score we use in the manuscript
  corr_net_pos <- corr_net %>% subset(corr > 0) %>% subset(corr >= corr_score)
  
  #Save the correlation results (just in case and time saving)
  links_pos <- data.frame(from = corr_net_pos$sp1, 
                          to = corr_net_pos$sp2, 
                          corr = corr_net_pos$corr)
  
  # -------------------------------------------------------------------------------------------------
  # save output table:
  if(save_output == "yes"){write_delim(links_pos, paste0(dir_out,"/links_pos_new.csv"), delim = ",", col_names = TRUE)}
  if(save_output == "yes"){write_delim(nodes,     paste0(dir_out,"/nodes_new.csv"), delim = ",", col_names = TRUE)}
  # -------------------------------------------------------------------------------------------------
  
  ############################################################################
  # 5 - Make network analysis using  graph_from_data_frame 
  ############################################################################
  
  net_pos <- igraph::graph_from_data_frame(d=links_pos, vertices=nodes, directed=F)
  #plot(net_pos)
  
  # change nodes info
  V(net_pos)$label <- paste(nodes$family, strtrim(nodes$identity, 4), sep = ".")
  
  ############################################################################
  # 6 - Detect communities with the cluster louvain method to optimise modularity
  ############################################################################
  
  set.seed(123)
  cl <- cluster_louvain(net_pos, resolution = 1.1) # alternative used in the manuscript
  
  # # check number of communities
  # length(cl) # 87 communities detected
  # 
  # # check membership
  # membership(cl)
  
  ############################################################################
  # 7 - Create a dataframe with membership info and Add number of ASV per community as metadata
  ############################################################################
  
  membercl <- data.frame(group = cl$membership, label = cl$names)
  membercl <- membercl %>% group_by(group) %>% mutate(nb_in_group = n_distinct(label)) %>% 
    left_join(metadata_taxa, by = c("label" = "short_label"))
  
  # -------------------------------------------------------------------------------------------------
  # save output table:
  if(save_output == "yes"){write_delim(membercl, paste0(dir_out,"/Communities_louvain_corrsup_",corr_score,"_more_100_reads_new.csv"), delim = ",", col_names = TRUE)}
  if(save_output == "yes"){write_delim(df, paste0(dir_out,"/DF_Communities_louvain_corrsup_",corr_score,"_more_100_reads_new.csv"), delim = ",", col_names = TRUE)}
  #--------------------------------------------------------------------------------------------------
  
  ############################################################################
  # 8 - Filter the ASVs into final potential taxa
  ############################################################################
  
  # keep reads  > 100
  df_join <- df %>% filter(reads_seq >= 100) %>% arrange(reads_seq) %>% 
    select(-id, -SeqID, -NUC_SEQ, -identity, -family, -genus, -species, -scientific_name, -nsamples, -ncores, -reads_seq)
  
  
  # Get only ASVs part of communities with minimum 5 ASVs
  member                <- left_join(membercl, df_join, by = "uniq_id") %>% ungroup()
  member_sup5_group     <- member %>% filter(nb_in_group >= 5) %>% arrange(group)
  member_sup5_not_group <- member %>% filter(nb_in_group < 5) %>% arrange(group)
  
  # -------------------------------------------------------------------------------------------------
  
  unique(member_sup5_group$group) 
  length(unique(member_sup5_group$group))
  dim(member_sup5_group)
  dim(member_sup5_not_group)
  
  commun_groups <- unique(member_sup5_group$group)
  
  # -------------------------------------------------------------------------------------------------
  
  sequence_taxanames_groups <- member_sup5_group %>% select(group, SeqID, NUC_SEQ, short_label, identity, family, genus, species, scientific_name) %>% distinct()
  
  member_sequence_sup5 <- member_sup5_group %>% select(group, uniq_id, NUC_SEQ, identity, family, genus, species, scientific_name) %>% mutate(status = ifelse(identity == 1, 1, "cand")) %>% mutate(uniq_label = paste(group, status, family, scientific_name, sep = "_"))
  
  suptable1 <- member_sup5_group %>% mutate(type = ifelse(identity == 1, "dbASV", "candidate ASV")) %>%
    select(NUC_SEQ, type, identity, family, genus, species, scientific_name, nsamples, ncores, reads_seq) %>% arrange(family) %>% arrange(desc(type))
  
  samples_name <- member %>% 
    select(-group, -label, -nb_in_group, -uniq_id, -SeqID, -NUC_SEQ, -identity, -family, -genus, 
           -species, -scientific_name, -nsamples, -ncores, -reads_seq, -id, -short_label) %>%
    names()
  
  # -------------------------------------------------------------------------------------------------
  # save output table:
  if(save_output == "yes"){write_delim(member_sequence_sup5, paste0(dir_out,"/ASVs_taxa_after_louvain_community.csv"), delim = ",", col_names = TRUE)}
  if(save_output == "yes"){write_delim(suptable1, paste0(dir_out,"/Supplementary_table1.csv"), delim = ",", col_names = TRUE)}
  if(save_output == "yes"){write_delim(member_sup5_group, paste0(dir_out,"/Supplementary_info_after_louvain_community_detection.csv"), delim = ",", col_names = TRUE)}
  # -------------------------------------------------------------------------------------------------
  
  ############################################################################
  # 8.1 - Merge ASVs with same taxa names - genus level (unused in manuscript)
  ############################################################################
  
  merge_assignments_prep <- member_sup5_group %>% select(group, nb_in_group, identity, family, genus, species, scientific_name, all_of(samples_name)) %>%
    pivot_longer(cols = all_of(samples_name), names_to = "sample_id", values_to = "reads_seq")
  
  merge_assignments_100percent <- merge_assignments_prep %>% filter(identity == 1 & reads_seq > 0) %>% arrange(sample_id) %>%
    group_by(group, scientific_name, sample_id) %>% mutate(reads_seq = sum(reads_seq)) %>% ungroup() %>% distinct()
  
  ############################################################################
  # 8.2 - Merge ASVs with same taxa names - family level (used in manuscript)
  ############################################################################
  
  merge_assignments1 <- merge_assignments_100percent
  
  merge_assignments1 <- merge_assignments1 %>% group_by(group, nb_in_group, family) %>% mutate(minid = min(identity), maxid = max(identity)) %>% 
    arrange(group) %>% ungroup() 
  
  merge_assignments_clean1 <- merge_assignments1 %>% select(-minid, -maxid) %>% distinct() %>% 
    group_by(group, identity, scientific_name) %>% mutate(tot_reads = sum(reads_seq), nsamples = n_distinct(sample_id)) %>% ungroup() %>% 
    filter(tot_reads >= 100) %>% filter(nsamples >= 10)
  
  merge_assignments_wide1 <- merge_assignments_clean1 %>% pivot_wider(names_from = sample_id, values_from = reads_seq) %>% arrange(tot_reads) %>% 
    group_by(group) %>% mutate(nb_in_group = n_distinct(scientific_name)) %>% arrange(group) %>% ungroup()
  
  long_to_wide1 <- merge_assignments_wide1 %>% pivot_longer(cols = -c(1:9), names_to = "sample_id", values_to = "reads_seq")
  
  long_to_wide1$reads_seq[is.na(long_to_wide1$reads_seq)] <- 0
  merge_assignments_wide_clean1 <- long_to_wide1 %>% pivot_wider(names_from = sample_id, values_from = reads_seq)
  
  # -------------------------------------------------------------------------------------------------
  # save output table:
  if(save_output == "yes"){write_delim(merge_assignments_wide_clean1, paste0(dir_out,"/louvain_merged_assignment_rdy_for_resampling_merge_family_level.csv"), delim = ",", col_names = TRUE)}
  # -------------------------------------------------------------------------------------------------
  
}

