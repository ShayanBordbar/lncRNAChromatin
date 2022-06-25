args <- commandArgs(trailingOnly = TRUE)


my_lnc_name <- args[1]
my_partition_type_cur <- args[2] # either "R" or "C"
my_output_address <- args[3] # "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_lowres1MG/Gviz"
my_model_names <- args[4:length(args)] #c("dist", "kmerFq", "motSc_pairs_triplx", "sequ_dist", "sequ_dist_trnsp","sequ_dist_Chrom_Meth_AX_ChIP_trnsp")

source("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/RF_evaluation_functions.R")
load("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/distance_feature_partition_6_lowRes_dfGRs.RData")
load("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_random_chunk_cv_df_updated.RData")
load("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_dfs_GR.RData")
load("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/distance_feature_partition_6_lowRes_df.RData")

aaw <- which(Partition_6_random_chunk_cv_df$owner %in% my_lnc_name)
my_partition <- Partition_6_random_chunk_cv_df[aaw,]
distance_feature_partition_6_lowRes_dfGR$score <- distance_feature_partition_6_lowRes_dfGR$lowerres1Mg

aammmfl <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_lowres1MG/Learned_models/", my_lnc_name)
aamy_feature_Granges_list <- list()
my_distance_feature_partition_6_lowRes_df <- distance_feature_partition_6_lowRes_dfGR[distance_feature_partition_6_lowRes_df$owner == my_lnc_name]
#my_new_dist <- my_distance_feature_partition_6_lowRes_df$lowerres1Mg[match(aanammmsp, my_distance_feature_partition_6_lowRes_df$tile_name)]
aamy_feature_Granges_list[[1]] <- my_distance_feature_partition_6_lowRes_df
names(aamy_feature_Granges_list) <- "Genomic_distance"
plot_gviz_lncRNA(my_tile_GR=Partition_6_dfs_GR,
                 my_partition_type= my_partition_type_cur,
                 my_partition_dataset = my_partition,
                 my_model_folder = aammmfl,
                 my_model_name = my_model_names,
                 my_nu_partitions=5,
                 my_feature_Granges_list= aamy_feature_Granges_list,
                 my_Gnome = "mm9",
                 output_address = my_output_address,
                 file_pre_name = "Gvisual",
                 my_coord_grange = numeric(0),
                 lnc_name = my_lnc_name)