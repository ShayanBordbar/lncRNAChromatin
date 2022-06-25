args <- commandArgs(trailingOnly = TRUE)

my_lnc_name <- args[1]
my_partition_type_cur <- args[2] # either "R" or "C"
my_output_address <- args[3] # "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_200feat/Gviz"
nu_feat_show <- as.numeric(args[4])
my_model_names <- args[5] #c("dist")
model_input_address <- args[6] # "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_CVfirst_dataset_2410003L11Rik.RData"
to_be_filtered <- args[7]
source("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/RF_evaluation_functions.R")
source("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/RBP_filter_func.R")

load("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/distance_feature_partition_6_lowRes_dfGRs.RData")
load("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_random_chunk_cv_df_updated.RData")
load("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_dfs_GR.RData")
load("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/distance_feature_partition_6_lowRes_df.RData")
load(model_input_address)

stopifnot(length(my_model_names) == 1)

if(to_be_filtered != "no_filter"){
  my_filter_index <- as.numeric(read.table(to_be_filtered)$V1)
  print("entering Filter_RBP")
  my_Dataset <- Filter_RBP(inp_dataset = my_Dataset,
                           exp_feature_name = c("GROseq", "PROseq"),
                           inp_nameDic = my_name_dic,
                           feature_index=my_filter_index,
                           partition_Exp_features_address = "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_expression_features.RData")
}

aaw <- which(Partition_6_random_chunk_cv_df$owner %in% my_lnc_name)
my_partition <- Partition_6_random_chunk_cv_df[aaw,]
distance_feature_partition_6_lowRes_dfGR$score <- distance_feature_partition_6_lowRes_dfGR$lowerres1Mg
aammmfl <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_200feat/Learned_models/", my_lnc_name)
aamy_feature_Granges_list <- list()
my_distance_feature_partition_6_lowRes_df <- distance_feature_partition_6_lowRes_dfGR[distance_feature_partition_6_lowRes_df$owner == my_lnc_name]

aamy_feature_Granges_list[[1]] <- my_distance_feature_partition_6_lowRes_df
names(aamy_feature_Granges_list) <- "Genomic_distance"


aa_files <- list.files(aammmfl, pattern = "*_perf_pred_varIMP.RData")
model_file_list <- aa_files[grep(pattern = paste0("___",my_model_names,"___"), x = aa_files)]
model_file_list <- model_file_list[grep(pattern = paste0("___",my_partition_type_cur,"CV"), x = model_file_list)]
model_part_list <- unlist(lapply(strsplit(model_file_list, "___"), "[[", 2))


imporance_list <- list()
aa_allnn <- character(0)
for(cur_pp in 1:length(model_part_list)){
  load(paste0(aammmfl, "/", model_file_list[cur_pp]))
  imporance_list[[cur_pp]] <- variable_importance
  aa_allnn <- union(aa_allnn, names(variable_importance))
}
#aa_allnn <- sort(aa_allnn)
imporance_mat <- matrix(nrow = length(model_file_list), ncol = length(aa_allnn))
colnames(imporance_mat) <- aa_allnn
#sort(variable_importance, decreasing = T)[1:min(50, length(variable_importance))]
for(cur_pp in 1:length(model_part_list)){
  #load(paste0(aammmfl, "/", model_file_list[cur_pp]))
  #aawwww <- which(colnames(imporance_mat) %in% names(imporance_list[[cur_pp]]))
  imporance_mat[cur_pp,match(names(imporance_list[[cur_pp]]), colnames(imporance_mat))] <- imporance_list[[cur_pp]]
}
#print(imporance_mat)
imporance_mat_sum <- colSums(imporance_mat, na.rm = T)
important_features <- colnames(imporance_mat)[sort(imporance_mat_sum, decreasing = T, index.return = T)$ix[1:(min(nu_feat_show, ncol(imporance_mat)))]]
#print(important_features)
aatransformed <- grep(pattern = "_TR_", important_features)
if(length(aatransformed) > 0){
  print("removing transformed features")
  important_features <- important_features[setdiff(c(1:length(important_features)), aatransformed)]
}
aaamymatch <- match(important_features, my_name_dic[,1])
stopifnot(sum(is.na(aaamymatch)) == 0)
for(i in 1:length(important_features)){
  
  aamy_feature_Granges_list[[i + 1]] <- Partition_6_dfs_GR[match(my_partition$tile_name , Partition_6_dfs_GR$tile)]
  aamy_feature_Granges_list[[i + 1]]$score <- my_Dataset[, aaamymatch[i]]
}
names(aamy_feature_Granges_list)[2:length(aamy_feature_Granges_list)] <- important_features
print("names(aamy_feature_Granges_list)")
print(names(aamy_feature_Granges_list))

aattsst <- evaluate_byWindow(my_tile_GR = Partition_6_dfs_GR,
                             my_partition_dataset = my_partition,
                             my_model_folder = aammmfl,
                             my_model_name = my_model_names)
aamy_range <- aattsst[[1]][sort(aattsst[[1]]$auprc, decreasing = T, index.return=T, na.last=T)$ix[1:min(length(aattsst[[1]]), 20)]]
aamy_range2 <- aattsst[[1]][sort(aattsst[[1]]$auroc, decreasing = T, index.return=T, na.last=T)$ix[1:min(length(aattsst[[1]]), 20)]]
aamy_range <- c(aamy_range, aamy_range2)

plot_gviz_lncRNA(my_tile_GR=Partition_6_dfs_GR,
                 my_partition_type= my_partition_type_cur,
                 my_partition_dataset = my_partition,
                 my_model_folder = aammmfl,
                 my_model_name = my_model_names,
                 my_nu_partitions=5,
                 my_feature_Granges_list= aamy_feature_Granges_list,
                 my_Gnome = "mm9",
                 output_address = my_output_address,
                 file_pre_name = "TopFeatures",
                 my_coord_grange = aamy_range,
                 lnc_name = my_lnc_name)

print("Top feature plot is now drawn")

