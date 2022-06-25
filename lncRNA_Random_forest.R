# running random forest on lncRNA dataset

setwd("~/Documents/Shayan/BioInf/lncRNA/")

library(aod)
library(ROCR)
library(caret)
library(DMwR)
library(ROSE)
library(VGAM)
library(LiblineaR)
library(MLmetrics)
library(PRROC)
library(purrr)
library(dplyr)
library(scales)
library(ggplot2)
library(ranger)
library(gplots)
####################################################################################################################################################################################
####################################################################################################################################################################################
########################################################################      FUNCTIONS      #######################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
source("~/Documents/Shayan/BioInf/lncRNA/RF_evaluation_functions.R")


################################################################################################


marginal_perf_impro <- function(perf_list,
                                default_mat, 
                                return_comp = F){
  # perf_list : is a list where each entry is performance for one lncRNA, something like aa_all_auroc_list_cv2TX_NoD -->
  #  each entry is a matrix with each row corresponding to a model and each column corresponding to a condition
  # return_comp : if True it won't average over comparisons and the larger dataframe which has a column for comparisons will be reported, if False, all comparisons for each family are averaged.
  my_family_mat_list_all2TX_NoD <- list()
  
  for(i in 1:length(perf_list)){
    print(names(perf_list)[i])
    family_mat_list <- list()
    aana <- rownames(perf_list[[i]])
    aana_sp <- strsplit(aana, "_")
    aana_sp_sort <- lapply(aana_sp, sort)
    aana_allfam <- sort(unique(unlist(aana_sp)))
    if("sequ" %in% aana_allfam){
      aana_allfam <- setdiff(aana_allfam,c("sequ"))
    }
    
    for(j in 1:length(aana_allfam)){
      aatmp1 <- which(unlist(lapply(aana_sp, function(x) aana_allfam[j] %in% x)))
      if(!(aana_allfam[j] %in% c("TXRB", "TXRE"))){
        aatmp111 <- which(unlist(lapply(aana_sp[aatmp1], function(x) "TXRB" %in% x)))
        aatmp1112 <- which(unlist(lapply(aana_sp[aatmp1], function(x) "TXRE" %in% x)))
        aatmp1 <- aatmp1[setdiff(c(1:length(aatmp1)), union(aatmp111, aatmp1112))]
      }else{
        aatmp111 <- which(unlist(lapply(aana_sp[aatmp1], function(x) "sequ" %in% x)))
        #aatmp1 <- aatmp1[aatmp111] # evaluating the importance of transcriptional filtering only in models that already use the transcription information separately --> so the improvement is not due to introduction of new information
        
      }
      #print(aana_allfam[j])
      #print(aana[aatmp1])
      aatmp2 <- numeric(length = length(aatmp1))
      for(k in 1:length(aatmp1)){
        if(length(aana_sp[[aatmp1[k]]]) > 1){
          aanew_name <- setdiff(aana_sp_sort[[aatmp1[k]]], aana_allfam[j])
          #if((length(aanew_name) == 2) & identical(aanew_name, c("dist", "sequ"))){
          #aatmp2[k] <- 0
          #}else{
          aastm <- which(unlist(lapply(aana_sp_sort, function(x) identical(aanew_name, x))))
          if(length(aastm) == 0){
            print(paste("nomatch", names(perf_list)[i], aana_allfam[j],"###", paste(aana_sp_sort[[aatmp1[k]]], collapse = "_"),"###",
                        paste(aanew_name, collapse = "||")))
            aatmp2[k] <- 0
          }else{
            aatmp2[k] <- aastm 
          }
        }else{
          aatmp2[k] <- 0
        }
        
      }
      family_mat_list[[j]] <- cbind(aatmp1, aatmp2)
    }
    names(family_mat_list) <- aana_allfam
    my_family_mat_list_all2TX_NoD[[i]] <- family_mat_list
  }
  
  names(my_family_mat_list_all2TX_NoD) <- names(perf_list)
  
  
  # calcuating the difference in perf
  
  
  aatmp_lnc_df_list_cv2TX_NoD <- list()
  for(i in 1:length(my_family_mat_list_all2TX_NoD)){
    aatmp_lnc_fam_df_list <- list()
    print(names(my_family_mat_list_all2TX_NoD)[i])
    for(j in 1:length(my_family_mat_list_all2TX_NoD[[i]])){#going over families
      print(names(my_family_mat_list_all2TX_NoD[[i]])[j])
      aatmptroc <- matrix(nrow = nrow(my_family_mat_list_all2TX_NoD[[i]][[j]]), 
                          ncol = ncol(perf_list[[i]]))
      # aatmptprc <- matrix(nrow = nrow(my_family_mat_list_all2TX_NoD[[i]][[j]]), 
      #                     ncol = ncol(aa_all_auprc_list_cv2TX_NoD[[i]]))
      comp1_name <- character(nrow(my_family_mat_list_all2TX_NoD[[i]][[j]]))
      comp2_name <- character(nrow(my_family_mat_list_all2TX_NoD[[i]][[j]]))
      
      for(k in 1:nrow(my_family_mat_list_all2TX_NoD[[i]][[j]])){ #going over comparisons
        if(my_family_mat_list_all2TX_NoD[[i]][[j]][k,2] == 0){
          aatmpcomproc <- default_mat[i, ]
          comp2_name[k] <- "baseline"
          #aatmpcompprc <- aa_auprc_base[i, ]
        }else{
          aatmpcomproc <- perf_list[[i]][my_family_mat_list_all2TX_NoD[[i]][[j]][k,2],]
          comp2_name[k] <- rownames(perf_list[[i]])[my_family_mat_list_all2TX_NoD[[i]][[j]][k,2]]
          #aatmpcompprc <- aa_all_auprc_list_cv2TX_NoD[[i]][my_family_mat_list_all2TX_NoD[[i]][[j]][k,2],]
        }
        comp1_name[k] <- rownames(perf_list[[i]])[my_family_mat_list_all2TX_NoD[[i]][[j]][k,1]]
        aatmptroc[k,] <- perf_list[[i]][my_family_mat_list_all2TX_NoD[[i]][[j]][k,1],] - aatmpcomproc
        #aatmptprc[k,] <- aa_all_auprc_list_cv2TX_NoD[[i]][my_family_mat_list_all2TX_NoD[[i]][[j]][k,1],] - aatmpcompprc
      }
      aatmptroc_num <- colSums(!is.na(aatmptroc))
      #aatmptprc_num <- colSums(!is.na(aatmptprc))
      
      aatmptroc_mean <- colMeans(aatmptroc, na.rm = T)
      #aatmptprc_mean <- colMeans(aatmptprc, na.rm = T)
      aatmptroc_sd <- apply(aatmptroc, MARGIN = 2, FUN = sd, na.rm = T)
      #aatmptprc_sd <- apply(aatmptprc, MARGIN = 2, FUN = sd, na.rm = T)
      aatmptroc_max <- apply(aatmptroc, MARGIN = 2, FUN = max, na.rm = T)
      #aatmptprc_max <- apply(aatmptprc, MARGIN = 2, FUN = max, na.rm = T)
      aa_my_pt <- c(paste0("CCV", c(1:5)),paste0("RCV", c(1:5)))
      aa_my_pt_type <- c(rep("chunk", 5), rep("random", 5))
      if(return_comp){
        aatmp_lnc_fam_df_list[[j]] <- data.frame(lncRNA = rep(names(my_family_mat_list_all2TX_NoD)[i], length(aatmptroc)), 
                                                 feature_family = rep(names(my_family_mat_list_all2TX_NoD[[i]])[j], length(aatmptroc)), 
                                                 partition=   rep(aa_my_pt, each = nrow(aatmptroc)),
                                                 partition_type = rep(aa_my_pt_type, each = nrow(aatmptroc)),
                                                 comp1 = rep(comp1_name, ncol(aatmptroc)), 
                                                 comp2 = rep(comp2_name, ncol(aatmptroc)),
                                                 perf_delta = as.numeric(aatmptroc))
        
        
      }else{
        aatmp_lnc_fam_df_list[[j]] <- data.frame(lncRNA = rep(names(my_family_mat_list_all2TX_NoD)[i], length(aatmptroc_mean)), 
                                                 feature_family = rep(names(my_family_mat_list_all2TX_NoD[[i]])[j], length(aatmptroc_mean)), 
                                                 partition=c(paste0("CCV", c(1:5)),paste0("RCV", c(1:5))),
                                                 partition_type = c(rep("chunk", 5), rep("random", 5)),
                                                 nu_comparison = aatmptroc_num, 
                                                 perf_mean_delta = aatmptroc_mean,
                                                 perf_sd_delta = aatmptroc_sd,
                                                 #auprc_mean_delta = aatmptprc_mean,
                                                 #auprc_sd_delta = aatmptprc_sd,
                                                 perf_max_delta = aatmptroc_max)#,
        #auprc_max_delta = aatmptprc_max)
      }

      
    } 
    aatmp_lnc_df_list_cv2TX_NoD[[i]] <- do.call(rbind, aatmp_lnc_fam_df_list)
    
  }
  
  aatmp_lnc_df_cv2TX_NoD <- do.call(rbind, aatmp_lnc_df_list_cv2TX_NoD)
  return(list(marginal_perf_df = aatmp_lnc_df_cv2TX_NoD))
  
}
################################################################################################
# example
aadefmat <- aa_auprc_base


aatstss <- marginal_perf_impro(perf_list = aa_all_auprc_list_cv2TX_NoD,
                               default_mat = aadefmat)

aatstss2 <- marginal_perf_impro(perf_list = aa_all_auprc_list_cv2TX_NoD,
                               default_mat = aadefmat,return_comp = T)

################################################################################################
# function to get marginal improvement of pairs of features

marginal_perf_impro_pair <- function(perf_list,
                                     my_pairs = numeric(0),
                                default_mat, 
                                return_comp = F){
  my_family_mat_list_all2TX_NoD <- list()
  
  if(length(my_pairs) == 0){
    my_pairs <- rbind(c("kmerFq", "motSc"),
                      c("kmerFq", "pairs"),
                      c("kmerFq", "Rep"),
                      c("kmerFq", "triplx"),
                      c("motSc", "pairs"),
                      c("motSc", "Rep"),
                      c("motSc", "triplx"),
                      c("pairs", "Rep"),
                      c("pairs", "triplx"),
                      c("Rep", "triplx"),
                      c("AX", "ChIP"),
                      c("AX", "Chrom"),
                      c("AX", "Meth"),
                      c("AX", "trnsp"),
                      c("ChIP", "Chrom"),
                      c("ChIP", "Meth"),
                      c("ChIP", "trnsp"),
                      c("Chrom", "Meth"),
                      c("Chrom", "trnsp"),
                      c("trnsp", "Meth")
                      )
  }
  for(i in 1:length(perf_list)){
    print(names(perf_list)[i])
    family_mat_list <- list()
    aana <- rownames(perf_list[[i]])
    aana_sp <- strsplit(aana, "_")
    aana_sp_sort <- lapply(aana_sp, sort)
    aana_allfam <- sort(unique(unlist(aana_sp)))
    if("sequ" %in% aana_allfam){
      aana_allfam <- setdiff(aana_allfam,c("sequ"))
    }
    aa_pair_name <- character(nrow(my_pairs))
    for(j in 1:nrow(my_pairs)){
      aa_pair_name[j] <- paste0( my_pairs[j, 1], "_AND_" ,my_pairs[j, 2])
      aatmp1p <- which(unlist(lapply(aana_sp, function(x) my_pairs[j, 1] %in% x)))
      aatmp2p <- which(unlist(lapply(aana_sp, function(x) my_pairs[j, 2] %in% x)))
      aatmp1both <- intersect(aatmp1p, aatmp2p)
      aatmp111 <- which(unlist(lapply(aana_sp[aatmp1both], function(x) "TXRB" %in% x)))
      aatmp1112 <- which(unlist(lapply(aana_sp[aatmp1both], function(x) "TXRE" %in% x)))
      aatmp1both <- aatmp1both[setdiff(c(1:length(aatmp1both)), union(aatmp111, aatmp1112))]
  
      #print(aana_allfam[j])
      #print(aana[aatmp1])
      aatmp2 <- numeric(length = length(aatmp1both))
      for(k in 1:length(aatmp1both)){
        if(length(aana_sp[[aatmp1both[k]]]) > 2){
          aanew_name <- setdiff(aana_sp_sort[[aatmp1both[k]]], my_pairs[j, ])

          aastm <- which(unlist(lapply(aana_sp_sort, function(x) identical(aanew_name, x))))
          if(length(aastm) == 0){
            print(paste("nomatch", names(perf_list)[i], my_pairs[j, 1],my_pairs[j, 2],"###", 
                        paste(aana_sp_sort[[aatmp1both[k]]], collapse = "_"),"###",
                        paste(aanew_name, collapse = "||")))
            aatmp2[k] <- 0
          }else{
            aatmp2[k] <- aastm 
          }
        }else{
          aatmp2[k] <- 0
        }
        
      }
      family_mat_list[[j]] <- cbind(aatmp1both, aatmp2)
    }
    names(family_mat_list) <- aa_pair_name
    my_family_mat_list_all2TX_NoD[[i]] <- family_mat_list
  }
  
  names(my_family_mat_list_all2TX_NoD) <- names(perf_list)
  
  
  # calcuating the difference in perf
  
  
  aatmp_lnc_df_list_cv2TX_NoD <- list()
  for(i in 1:length(my_family_mat_list_all2TX_NoD)){ # going over lncRNAs
    aatmp_lnc_fam_df_list <- list()
    print(names(my_family_mat_list_all2TX_NoD)[i])
    for(j in 1:length(my_family_mat_list_all2TX_NoD[[i]])){#going over family pairs
      print(names(my_family_mat_list_all2TX_NoD[[i]])[j])
      aatmptroc <- matrix(nrow = nrow(my_family_mat_list_all2TX_NoD[[i]][[j]]), 
                          ncol = ncol(perf_list[[i]]))

      comp1_name <- character(nrow(my_family_mat_list_all2TX_NoD[[i]][[j]]))
      comp2_name <- character(nrow(my_family_mat_list_all2TX_NoD[[i]][[j]]))
      
      for(k in 1:nrow(my_family_mat_list_all2TX_NoD[[i]][[j]])){ #going over comparisons
        if(my_family_mat_list_all2TX_NoD[[i]][[j]][k,2] == 0){
          aatmpcomproc <- default_mat[i, ]
          comp2_name[k] <- "baseline"
          
        }else{
          aatmpcomproc <- perf_list[[i]][my_family_mat_list_all2TX_NoD[[i]][[j]][k,2],]
          comp2_name[k] <- rownames(perf_list[[i]])[my_family_mat_list_all2TX_NoD[[i]][[j]][k,2]]
          
        }
        comp1_name[k] <- rownames(perf_list[[i]])[my_family_mat_list_all2TX_NoD[[i]][[j]][k,1]]
        aatmptroc[k,] <- perf_list[[i]][my_family_mat_list_all2TX_NoD[[i]][[j]][k,1],] - aatmpcomproc
        
      }
      aatmptroc_num <- colSums(!is.na(aatmptroc))
      aatmptroc_mean <- colMeans(aatmptroc, na.rm = T)
      aatmptroc_sd <- apply(aatmptroc, MARGIN = 2, FUN = sd, na.rm = T)
      aatmptroc_max <- apply(aatmptroc, MARGIN = 2, FUN = max, na.rm = T)
      aa_my_pt <- c(paste0("CCV", c(1:5)),paste0("RCV", c(1:5)))
      aa_my_pt_type <- c(rep("chunk", 5), rep("random", 5))
      if(return_comp){
        aatmp_lnc_fam_df_list[[j]] <- data.frame(lncRNA = rep(names(my_family_mat_list_all2TX_NoD)[i], length(aatmptroc)), 
                                                 feature_family = rep(names(my_family_mat_list_all2TX_NoD[[i]])[j], length(aatmptroc)), 
                                                 partition=   rep(aa_my_pt, each = nrow(aatmptroc)),
                                                 partition_type = rep(aa_my_pt_type, each = nrow(aatmptroc)),
                                                 comp1 = rep(comp1_name, ncol(aatmptroc)), 
                                                 comp2 = rep(comp2_name, ncol(aatmptroc)),
                                                 perf_delta = as.numeric(aatmptroc))
        
        
      }else{
        aatmp_lnc_fam_df_list[[j]] <- data.frame(lncRNA = rep(names(my_family_mat_list_all2TX_NoD)[i], length(aatmptroc_mean)), 
                                                 feature_family = rep(names(my_family_mat_list_all2TX_NoD[[i]])[j], length(aatmptroc_mean)), 
                                                 partition=c(paste0("CCV", c(1:5)),paste0("RCV", c(1:5))),
                                                 partition_type = c(rep("chunk", 5), rep("random", 5)),
                                                 nu_comparison = aatmptroc_num, 
                                                 perf_mean_delta = aatmptroc_mean,
                                                 perf_sd_delta = aatmptroc_sd,
                                                 #auprc_mean_delta = aatmptprc_mean,
                                                 #auprc_sd_delta = aatmptprc_sd,
                                                 perf_max_delta = aatmptroc_max)#,
        #auprc_max_delta = aatmptprc_max)
      }
      
      
    } 
    aatmp_lnc_df_list_cv2TX_NoD[[i]] <- do.call(rbind, aatmp_lnc_fam_df_list)
    
  }
  
  aatmp_lnc_df_cv2TX_NoD <- do.call(rbind, aatmp_lnc_df_list_cv2TX_NoD)
  return(list(marginal_perf_df = aatmp_lnc_df_cv2TX_NoD))
  
}
####################################################################################################################################################################################
#example

aatyty3 <- marginal_perf_impro_pair(perf_list = aa_all_auprc_list_cv2TX_NoD,
                                    default_mat = aadefmat)
aatyty4 <- marginal_perf_impro_pair(perf_list = aa_all_auprc_list_cv2TX_NoD,
                                    default_mat = aadefmat,
                                    return_comp = T)
####################################################################################################################################################################################
####################################################################################################################################################################################
# compare pair with sum of individuals
pair_indiv_compare <- function(pair_df, indiv_df){
  aa_all_pairs <- unique(pair_df$feature_family)
  aa_all_pairs_sp <- strsplit(aa_all_pairs, split = "_AND_")
  aa_all_lnc <-  unique(pair_df$lncRNA)
  aa_lnc_list <- list()
  for(i in 1:length(aa_all_lnc)){ # going over lncRNAs
    cur_df_pair <- pair_df[pair_df$lncRNA %in% aa_all_lnc[i],]
    cur_df_indi <- indiv_df[indiv_df$lncRNA %in% aa_all_lnc[i],]
    aa_pair_list <- list()
    for(j in 1:length(aa_all_pairs)){
      cur_df_pair_1 <- cur_df_pair[cur_df_pair$feature_family %in% aa_all_pairs[j],]
      cur_df_indi_1 <- cur_df_indi[cur_df_indi$feature_family %in% aa_all_pairs_sp[[j]][1], ]
      cur_df_indi_2 <- cur_df_indi[cur_df_indi$feature_family %in% aa_all_pairs_sp[[j]][2], ]
      colnames(cur_df_pair_1)[ncol(cur_df_pair_1)] <- "perf_pair"
      colnames(cur_df_indi_1)[ncol(cur_df_indi_1)] <- "perf_1"
      colnames(cur_df_indi_2)[ncol(cur_df_indi_2)] <- "perf_2"
      cur_df_merge1 <- merge(cur_df_indi_1, cur_df_indi_2, by = c("lncRNA","partition_type","partition", "comp1")) 
      cur_df_merge2 <- merge(cur_df_pair_1, cur_df_merge1, by = c("lncRNA","partition_type","partition", "comp1")) 
      aa_pair_list[[j]] <- cur_df_merge2
    }
    aa_lnc_list[[i]] <- do.call(rbind, aa_pair_list)
  }
  aa_all_df <- do.call(rbind, aa_lnc_list)
  return(aa_all_df)
}
####################################
# example
aamytss <- pair_indiv_compare(pair_df = aatstss31$marginal_perf_df, indiv_df = aatstss3$marginal_perf_df)
####################################################################################################################################################################################
####################################################################################################################################################################################
########################################################################      Partition1      #######################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################

Partition_1_feature_mat_tile
Partition_1_feature_mat_owner
Partition_1_feature_mat_pair


save(list = c("Partition_1_feature_mat_tile", "Partition_1_feature_mat_owner", "Partition_1_feature_mat_pair", "Partition_1_dfs", "Partition_1_chunk_dfs"), 
     file = "~/Documents/Shayan/BioInf/lncRNA/Random_forest_features.RData")

##################################################################################################################################################################################
#small example for testing the rf
aa_samp <- sample(x = c(1:nrow(ChIPATLAS_features_partition1)), size = 1000, replace = F)
aa_minidataset <- cbind(ChIPATLAS_features_partition1[aa_samp, 1:10],
                        #Partition_1_feature_mat_owner[1:1000, 1:10],
                        Partition_1_feature_mat_pair[aa_samp, ],
                        Partition_1_dfs$label[aa_samp])
colnames(aa_minidataset)[ncol(aa_minidataset)] <- "label"

#aa_minidataset[, ncol(aa_minidataset)] <- as.factor(aa_minidataset[, ncol(aa_minidataset)])
aa_minidataset <- as.data.frame(aa_minidataset)
aa_minidataset$label[aa_minidataset$label == 0] <- "neg"
aa_minidataset$label[aa_minidataset$label == 1] <- "pos"

aa_minidataset$label <- factor(aa_minidataset$label , levels = c("neg", "pos"))
table(Partition_1_dfs$dataset[aa_samp])
aa_minidataset[is.na(aa_minidataset)] <- 0
aa_minidataset_train <- aa_minidataset[Partition_1_dfs$dataset[aa_samp] %in% c("train", "valid"),]
aa_minidataset_test <- aa_minidataset[Partition_1_dfs$dataset[aa_samp] %in% c("test"),]

aactrl <- trainControl(method = "repeatedcv",
                       number = 6,
                       repeats = 1,
                       classProbs = TRUE,
                       summaryFunction = defaultSummary,
                       savePredictions = TRUE,
                       ## new option here:
                       sampling = "down", 
                       search = 'random')

aaRF_down <- train(label ~ ., data = aa_minidataset_train,
                   method = "rf",
                   #preProc = c("center", "scale"),
                   ntree = 100,
                   metric = "Kappa",
                   trControl = aactrl
                   #,tuneLength  = aagridLen[1]
)

aarfdown <- predict(aaRF_down, 
                    aa_minidataset_test,
                    type = "prob")
cbind(aarfdown, aa_minidataset_test$label)
varImp(aaRF_down,scale = F)

#########
aa_all_feat <- cbind(Partition_1_feature_mat_tile, Partition_1_feature_mat_owner, Partition_1_feature_mat_pair, Partition_1_dfs$label)
aa_zerovar <- apply(aa_all_feat, MARGIN = 2, FUN = var)
Zero_variance_columns <- colnames(aa_all_feat)[aa_zerovar == 0]
aa_all_feat <- aa_all_feat[, -which(aa_zerovar == 0)]
aadescrCor <- cor(aa_all_feat[, 1:(ncol(aa_all_feat) - 1)])
aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .90)
Highly_correlating_columns <- colnames(aa_all_feat)[aahighlyCorDescr]

# 
# aanzv <- nearZeroVar(aa_all_feat, saveMetrics= F)
# if(length(aanzv) > 0){
#   aa_all_data <- aa_all_feat[,-aanzv]
# }else{
#   aa_all_data <-aa_all_feat
# }
aa_all_data<-aa_all_feat
#  Identifying Correlated Predictors
if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
  aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
  aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .90)
  if(length(aahighlyCorDescr) > 0){
    aa_all_data <- aa_all_data[,-aahighlyCorDescr]
  }
}
print("#########################################################################################")
print(paste0("number of features in original feature set number:"))
print(ncol(aa_all_feat))
print("number of features after removing near zero variance and highly correlated features ( > 0.9)")
print(ncol(aa_all_data))
print("#########################################################################################")

aatrainind <- which(Partition_1_dfs$dataset %in% c("train", "valid"))
aagridLen <- c(5)

aax <- colnames(Partition_1_modified_dataset)
aax2 <- gsub(pattern = " ", replacement = "", x = aax)
colnames(Partition_1_modified_dataset) <- aax2
name_dic <- cbind(colnames(Partition_1_modified_dataset), paste0("feature_", c(1:ncol(Partition_1_modified_dataset))))
name_dic[ncol(Partition_1_modified_dataset), 2] <- name_dic[ncol(Partition_1_modified_dataset), 1] 
colnames(Partition_1_modified_dataset) <- name_dic[, 2]

aa_train_data <- aa_all_data[aatrainind, ]
aa_test_data <-  aa_all_data[-aatrainind, ]
aactrl <- trainControl(method = "repeatedcv",
                       number = 6,
                       repeats = 3,
                       classProbs = TRUE,
                       summaryFunction = defaultSummary,
                       savePredictions = TRUE,
                       ## new option here:
                       sampling = "down", 
                       search = 'grid')
aaRF_down <- train(label ~ ., data = aa_train_data,
                   method = "rf",
                   #preProc = c("center", "scale"),
                   ntree = 1000,
                   metric = "Kappa",
                   trControl = aactrl
                   #,tuneLength  = aagridLen[1]
)
my_RF_results_1 <- aaRF_down


save(list = c("my_RF_results_1"),
     file = "RF_results_random_partitioning.RData")

# prediction
my_RF_results_1_prediction_list <- predict(my_RF_results_1, newdata=aa_test_data[, c(1:(ncol(aa_cur_test) - 1))], type="response")
aapr <- prediction(my_RF_results_1_prediction_list, aa_test_data$label)
my_RF_results_1_prediction_list_perf_list[[1]][[i]] <- performance(aapr, measure = "tpr", x.measure = "fpr")




aatrainind <- which(Partition_1_chunk_dfs$dataset %in% c("train", "valid"))
# aagridLen <- c(5, 30, 50)
aa_train_data <- aa_all_data[aatrainind, ]
aa_test_data <-  aa_all_data[-aatrainind, ]
aactrl <- trainControl(method = "repeatedcv",
                       number = 6,
                       repeats = 3,
                       classProbs = TRUE,
                       summaryFunction = defaultSummary,
                       savePredictions = TRUE,
                       ## new option here:
                       sampling = "down", 
                       search = 'grid')
aaRF_down2 <- train(label ~ ., data = aa_train_data,
                    method = "rf",
                    #preProc = c("center", "scale"),
                    ntree = 1000,
                    metric = "Kappa",
                    trControl = aactrl
                    #,tuneLength  = aagridLen[1]
)
my_RF_results_2 <- aaRF_down2


save(list = c("my_RF_results_2"),
     file = "RF_results_chunk_partitioning.RData")
######################################################################
# evaluating the first random and chunk partitions
# random:
load("RF_results_random_partitioning.RData")
my_RF_results_1
aapred1 <- prediction(my_logistic_prediction_list[[1]][[i]], aa_cur_test$label)

# chunk
load("RF_results_chunk_partitioning.RData")
my_RF_results_2

save(list = c("MM9_1kb_tiled_owner_3Mfilter"), file = "~/Documents/Shayan/BioInf/lncRNA/MM9_1kb_tiled_owner_3Mfilter.RData")

aa_train_1 <- perf_eval_ROC_PRC(model_fit = fit,
                                my_dataset = aa_train_data,
                                my_label = aa_train_data[, ncol(aa_train_data)], 
                                file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_train_perf_partition1_random.png")
aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_test_data), names(MM9_1kb_tiled_owner_3Mfilter))]

aa_tst_1_perclass <- perf_eval_ROC_PRC_perclass(model_fit = rf_model_partition1_random_varimp,
                                                my_dataset = aa_test_data,
                                                my_label = aa_test_data[, ncol(aa_test_data)], 
                                                file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_test_perf_partition1_random_perclass.png", 
                                                class_name = aaclass_name)
aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_train_data), names(MM9_1kb_tiled_owner_3Mfilter))]
aa_train_1_perclass <- perf_eval_ROC_PRC_perclass(model_fit = rf_model_partition1_random_varimp,
                                                  my_dataset = aa_train_data,
                                                  my_label = aa_train_data[, ncol(aa_train_data)], 
                                                  file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_train_perf_partition1_random_perclass.png",
                                                  class_name = aaclass_name, 
                                                  train_perf = T)


aa_tst_2 <- perf_eval_ROC_PRC(model_fit = fit2,
                              my_dataset = aa_test_data,
                              my_label = aa_test_data[, ncol(aa_test_data)], 
                              file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_test_perf_partition1_random_classWeighted.png")
aa_train_2 <- perf_eval_ROC_PRC(model_fit = fit2,
                                my_dataset = aa_train_data,
                                my_label = aa_train_data[, ncol(aa_train_data)], 
                                file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_train_perf_partition1_random_classWeighted.png")
aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_test_data), names(MM9_1kb_tiled_owner_3Mfilter))]
aa_tst_2_perclass <- perf_eval_ROC_PRC_perclass(model_fit = rf_model_partition1_random_varimp_purity_classWeight,
                                                my_dataset = aa_test_data,
                                                my_label = aa_test_data[, ncol(aa_test_data)], 
                                                file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_test_perf_partition1_random_classWeighted_perclass.png", 
                                                class_name = aaclass_name)
aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_train_data), names(MM9_1kb_tiled_owner_3Mfilter))]
aa_train_2_perclass <- perf_eval_ROC_PRC_perclass(model_fit = rf_model_partition1_random_varimp_purity_classWeight,
                                                  my_dataset = aa_train_data,
                                                  my_label = aa_train_data[, ncol(aa_train_data)], 
                                                  file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_train_perf_partition1_random_classWeighted_perclass.png",
                                                  class_name = aaclass_name, 
                                                  train_perf = T)

# for chunk partition eval
aa_tst_1_chunk <- perf_eval_ROC_PRC(model_fit = fit,
                                    my_dataset = aa_test_data,
                                    my_label = aa_test_data[, ncol(aa_test_data)], 
                                    file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_test_perf_partition1_chunk.png")
aa_train_1_chunk <- perf_eval_ROC_PRC(model_fit = fit,
                                      my_dataset = aa_train_data,
                                      my_label = aa_train_data[, ncol(aa_train_data)], 
                                      file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_train_perf_partition1_chunk.png")

aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_test_data), names(MM9_1kb_tiled_owner_3Mfilter))]
aa_tst_1_chunk_perclass <- perf_eval_ROC_PRC_perclass(model_fit = rf_model_partition1_chunk_varimp_purity,
                                                      my_dataset = aa_test_data,
                                                      my_label = aa_test_data[, ncol(aa_test_data)], 
                                                      file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_test_perf_partition1_chunk_perclass.png", 
                                                      class_name = aaclass_name)
aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_train_data), names(MM9_1kb_tiled_owner_3Mfilter))]
aa_train_1_chunk_perclass <- perf_eval_ROC_PRC_perclass(model_fit = rf_model_partition1_chunk_varimp_purity,
                                                        my_dataset = aa_train_data,
                                                        my_label = aa_train_data[, ncol(aa_train_data)], 
                                                        file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_train_perf_partition1_chunk_perclass.png",
                                                        class_name = aaclass_name, 
                                                        train_perf = T)
aa_tst_2_chunk <- perf_eval_ROC_PRC(model_fit = fit2,
                                    my_dataset = aa_test_data,
                                    my_label = aa_test_data[, ncol(aa_test_data)], 
                                    file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_test_perf_partition1_chunk_classWeighted.png")
aa_train_2_chunk <- perf_eval_ROC_PRC(model_fit = fit2,
                                      my_dataset = aa_train_data,
                                      my_label = aa_train_data[, ncol(aa_train_data)], 
                                      file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_train_perf_partition1_chunk_classWeighted.png")
aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_test_data), names(MM9_1kb_tiled_owner_3Mfilter))]
aa_tst_2_chunk_perclass <- perf_eval_ROC_PRC_perclass(model_fit = rf_model_partition1_chunk_varimp_purity_classWeight,
                                                      my_dataset = aa_test_data,
                                                      my_label = aa_test_data[, ncol(aa_test_data)], 
                                                      file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_test_perf_partition1_chunk_classWeighted_perclass.png", 
                                                      class_name = aaclass_name)
aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_train_data), names(MM9_1kb_tiled_owner_3Mfilter))]
aa_train_2_chunk_perclass <- perf_eval_ROC_PRC_perclass(model_fit = rf_model_partition1_chunk_varimp_purity_classWeight,
                                                        my_dataset = aa_train_data,
                                                        my_label = aa_train_data[, ncol(aa_train_data)], 
                                                        file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_train_perf_partition1_chunk_classWeighted_perclass.png",
                                                        class_name = aaclass_name, 
                                                        train_perf = T)

aaimp <- importance(fit)
#names(aaimp) <- name_dic[1:length(aaimp),1]
#names(aaimp)[sort(aaimp, decreasing = T, index.return = T)$ix[1:30]]
png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_varimp_partition1_random.png",    # create PNG for the heat map        
    width = 6*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,1), mar = c(4,9,4,4))

barplot(sort(aaimp, decreasing = T)[1:50], horiz = T, las =2)
dev.off()




aaimp <- importance(fit)
#names(aaimp) <- name_dic[1:length(aaimp),1]
#names(aaimp)[sort(aaimp, decreasing = T, index.return = T)$ix[1:30]]
png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_varimp_partition1_random.png",    # create PNG for the heat map        
    width = 6*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,1), mar = c(4,9,4,4))

barplot(sort(aaimp, decreasing = T)[1:50], horiz = T, las =2)
dev.off()



aaimp <- importance(fit)
names(aaimp) <- name_dic[1:length(aaimp),1]
#names(aaimp)[sort(aaimp, decreasing = T, index.return = T)$ix[1:30]]
png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_varimp_partition1_chunk.png",    # create PNG for the heat map        
    width = 6*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,1), mar = c(4,9,4,4))

barplot(sort(aaimp, decreasing = T)[1:50], horiz = T, las =2)
dev.off()


aaimp2 <- importance(fit2)
names(aaimp2) <- name_dic[1:length(aaimp2),1]
#names(aaimp2)[sort(aaimp2, decreasing = T, index.return = T)$ix[1:30]]
png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_varimp_partition1_random_weighted.png",    # create PNG for the heat map        
    width = 6*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,1), mar = c(4,9,4,4))

barplot(sort(aaimp2, decreasing = T)[1:50], horiz = T, las =2)
dev.off()


aaimp2 <- importance(fit2)
names(aaimp2) <- name_dic[1:length(aaimp2),1]
#names(aaimp2)[sort(aaimp2, decreasing = T, index.return = T)$ix[1:30]]
png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_varimp_partition1_chunk_weighted.png",    # create PNG for the heat map        
    width = 6*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,1), mar = c(4,9,4,4))

barplot(sort(aaimp2, decreasing = T)[1:50], horiz = T, las =2)
dev.off()



save(list = c("MM9_1kb_tiled_owner_labels_binary_3Mfilter"), file = "~/Documents/Shayan/BioInf/lncRNA/MM9_1kb_tiled_owner_labels_binary_3Mfilter.RData")

save(list = c("lncRNA_chosen_gt1k_uniqTiles"), file = "~/Documents/Shayan/BioInf/lncRNA/lncRNA_chosen_gt1k_uniqTiles.RData")


aamy_pred <- predict(rf_model_partition1_chunk_varimp_purity, data=aa_test_data)
aapreddd <-  aamy_pred$predictions[, 2]
names(aapreddd) <- rownames(aa_test_data)
aatst <- performance_viz_chromosome(all_tile_label=MM9_1kb_tiled_owner_labels_binary_3Mfilter,
                                    show_points = c("all"), 
                                    subset_points= c( "train", "test", "valid"),
                                    partition_dataset = Partition_1_chunk_dfs,
                                    predicted_pos_probilbilty = aapreddd,
                                    positive_thersh = 0.2,
                                    specific_chr = character(0),
                                    my_file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/performance_per_chr_chunk_partition1_thresh2.png",
                                    lncRNA_coor_df= lncRNA_chosen_gt1k_uniqTiles,
                                    chr_sizes_file="~/Documents/Shayan/BioInf/lncRNA/mm9_chr_sizes.txt")


aapreddd2 <-  rf_model_partition1_chunk_varimp_purity$predictions[, 2]
names(aapreddd2) <- rownames(aa_train_data)
aapreddd <- c(aapreddd, aapreddd2)

aatst_2_2 <- performance_viz_chromosome(all_tile_label=MM9_1kb_tiled_owner_labels_binary_3Mfilter,
                                        show_points = c("used"), 
                                        subset_points= c( "train", "test", "valid"),
                                        partition_dataset = Partition_1_chunk_dfs,
                                        predicted_pos_probilbilty = aapreddd,
                                        positive_thersh = 0.2,
                                        specific_chr = character(0),
                                        my_file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/performance_per_chr_chunk_partition1_thresh2_used.png",
                                        lncRNA_coor_df= lncRNA_chosen_gt1k_uniqTiles,
                                        chr_sizes_file="~/Documents/Shayan/BioInf/lncRNA/mm9_chr_sizes.txt")

aamy_pred <- predict(rf_model_partition1_chunk_varimp_purity_classWeight, data=aa_test_data)
aapreddd <-  aamy_pred$predictions[, 2]
names(aapreddd) <- rownames(aa_test_data)
aatst2 <- performance_viz_chromosome(all_tile_label=MM9_1kb_tiled_owner_labels_binary_3Mfilter,
                                     show_points = c("all"), 
                                     subset_points= c( "train", "test", "valid"),
                                     partition_dataset = Partition_1_chunk_dfs,
                                     predicted_pos_probilbilty = aapreddd,
                                     positive_thersh = 0.2,
                                     specific_chr = character(0),
                                     my_file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/performance_per_chr_chunk_partition1_class_weighted_thresh2.png",
                                     lncRNA_coor_df= lncRNA_chosen_gt1k_uniqTiles,
                                     chr_sizes_file="~/Documents/Shayan/BioInf/lncRNA/mm9_chr_sizes.txt")
#aamy_pred <- predict(rf_model_partition1_chunk_varimp_purity_classWeight, data=aa_test_data)
aapreddd2 <-  rf_model_partition1_chunk_varimp_purity_classWeight$predictions[, 2]
names(aapreddd2) <- rownames(aa_train_data)
aapreddd <- c(aapreddd, aapreddd2)
aatst3 <- performance_viz_chromosome(all_tile_label=MM9_1kb_tiled_owner_labels_binary_3Mfilter,
                                     show_points = c("used"), 
                                     subset_points= c( "train", "test", "valid"),
                                     partition_dataset = Partition_1_chunk_dfs,
                                     predicted_pos_probilbilty = aapreddd,
                                     positive_thersh = 0.2,
                                     specific_chr = character(0),
                                     my_file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/performance_per_chr_chunk_partition1_class_weighted_thresh2_used.png",
                                     lncRNA_coor_df= lncRNA_chosen_gt1k_uniqTiles,
                                     chr_sizes_file="~/Documents/Shayan/BioInf/lncRNA/mm9_chr_sizes.txt")
aamy_pred <- predict(rf_model_partition1_chunk_varimp_purity_classWeight, data=aa_test_data)
aapreddd <-  aamy_pred$predictions[, 2]
names(aapreddd) <- rownames(aa_test_data)
aapreddd2 <-  rf_model_partition1_chunk_varimp_purity_classWeight$predictions[, 2]
names(aapreddd2) <- rownames(aa_train_data)
aapreddd <- c(aapreddd, aapreddd2)
aatst3 <- performance_viz_chromosome(all_tile_label=MM9_1kb_tiled_owner_labels_binary_3Mfilter,
                                     show_points = c("used"), 
                                     subset_points= c( "train", "test", "valid"),
                                     partition_dataset = Partition_1_chunk_dfs,
                                     predicted_pos_probilbilty = aapreddd,
                                     positive_thersh = 0.2,
                                     specific_chr = character(0),
                                     my_file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/performance_per_chr_chunk_partition1_class_weighted_thresh2_used.png",
                                     lncRNA_coor_df= lncRNA_chosen_gt1k_uniqTiles,
                                     chr_sizes_file="~/Documents/Shayan/BioInf/lncRNA/mm9_chr_sizes.txt")


aamy_pred <- predict(rf_model_partition1_random_varimp, data=aa_test_data)
aapreddd <-  aamy_pred$predictions[, 2]
names(aapreddd) <- rownames(aa_test_data)
aapreddd2 <-  rf_model_partition1_random_varimp$predictions[, 2]
names(aapreddd2) <- rownames(aa_train_data)
aapreddd <- c(aapreddd, aapreddd2)
aatst4 <- performance_viz_chromosome(all_tile_label=MM9_1kb_tiled_owner_labels_binary_3Mfilter,
                                     show_points = c("used"), 
                                     subset_points= c( "train", "test", "valid"),
                                     partition_dataset = Partition_1_dfs,
                                     predicted_pos_probilbilty = aapreddd,
                                     positive_thersh = 0.1335,
                                     specific_chr = character(0),
                                     my_file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/performance_per_chr_random_partition1_thresh1335_used.png",
                                     lncRNA_coor_df= lncRNA_chosen_gt1k_uniqTiles,
                                     chr_sizes_file="~/Documents/Shayan/BioInf/lncRNA/mm9_chr_sizes.txt")

aamy_pred <- predict(rf_model_partition1_random_varimp_purity_classWeight, data=aa_test_data)
aapreddd <-  aamy_pred$predictions[, 2]
names(aapreddd) <- rownames(aa_test_data)
aapreddd2 <-  rf_model_partition1_random_varimp_purity_classWeight$predictions[, 2]
names(aapreddd2) <- rownames(aa_train_data)
aapreddd <- c(aapreddd, aapreddd2)
aatst4_wght <- performance_viz_chromosome(all_tile_label=MM9_1kb_tiled_owner_labels_binary_3Mfilter,
                                          show_points = c("used"), 
                                          subset_points= c( "train", "test", "valid"),
                                          partition_dataset = Partition_1_dfs,
                                          predicted_pos_probilbilty = aapreddd,
                                          positive_thersh = 0.1375,
                                          specific_chr = character(0),
                                          my_file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/performance_per_chr_random_partition1_thresh1375_used_weighted.png",
                                          lncRNA_coor_df= lncRNA_chosen_gt1k_uniqTiles,
                                          chr_sizes_file="~/Documents/Shayan/BioInf/lncRNA/mm9_chr_sizes.txt")
####################################################################################################################################################################################
####################################################################################################################################################################################
########################################################################      Partition2      #######################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################

aa_all_feat <- cbind(Partition_2_feature_mat_tile, Partition_2_feature_mat_owner, Partition_2_feature_mat_pair, Partition_2_dfs$label)
aa_all_feat[is.na(aa_all_feat)] <- 0
aavar_all <- apply(aa_all_feat, MARGIN = 2, var)


colnames(aa_all_feat)[ncol(aa_all_feat)] <- "label"
Zero_var_cols_partition2 <- colnames(aa_all_feat)[aavar_all == 0]
aa_all_feat <- aa_all_feat[, -which(aavar_all == 0)]

aa_all_feat <- as.data.frame(aa_all_feat)
aa_all_feat$label[aa_all_feat$label == 0] <- "neg"
aa_all_feat$label[aa_all_feat$label == 1] <- "pos"

aa_all_feat$label <- factor(aa_all_feat$label , levels = c("neg", "pos"))



#  Identifying Correlated Predictors
aadescrCor <- cor(aa_all_feat[, 1:(ncol(aa_all_feat) - 1)])
aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .90)
  # if(length(aahighlyCorDescr) > 0){
  #   aa_all_feat <- aa_all_feat[,-aahighlyCorDescr]
  # }
Highly_correlating_features_partition2 <- colnames(aa_all_feat)[aahighlyCorDescr]
aa_all_feat <- aa_all_feat[,-aahighlyCorDescr]


Partition_2_modified_dataset <- aa_all_feat
name_dic_p2 <- cbind(colnames(Partition_2_modified_dataset), paste0("feature_", c(1:ncol(Partition_2_modified_dataset))))
name_dic_p2[ncol(Partition_2_modified_dataset), 2] <- name_dic_p2[ncol(Partition_2_modified_dataset), 1] 
colnames(Partition_2_modified_dataset) <- name_dic_p2[, 2]
#save.image("Partition_2_modified_dataset_and_removedColumnNames.RData")

#random-class weighted

aatrainind <- which(Partition_2_dfs$dataset %in% c("train", "valid"))
aa_train_data <- Partition_2_modified_dataset[aatrainind, ]
aa_test_data <-  Partition_2_modified_dataset[-aatrainind, ]
rm(Partition_2_modified_dataset)

aa_table <- table(aa_train_data$label)
aa_class_weight <- c(aa_table[2]/nrow(aa_train_data), aa_table[1]/nrow(aa_train_data))
names(aa_class_weight) <- c("neg", "pos")
aa_case_weight <- numeric(length = nrow(aa_train_data))
aa_case_weight[aa_train_data$label == "pos"] <- aa_class_weight[2]
aa_case_weight[aa_train_data$label == "neg"] <- aa_class_weight[1]
set.seed(1233)
RF_partition2_random_classWght <- ranger(label ~ ., 
                                         data = aa_train_data, 
                                         num.trees = 1000,
                                         max.depth = 8,
                                         probability = TRUE, 
                                         importance  = "impurity",
                                         case.weights = aa_case_weight)


# chunk class-weighted

aatrainind <- which(Partition_2_chunk_dfs$dataset %in% c("train", "valid"))
aa_train_data <- Partition_2_modified_dataset[aatrainind, ]
aa_test_data <-  Partition_2_modified_dataset[-aatrainind, ]
rm(Partition_2_modified_dataset)

aa_table <- table(aa_train_data$label)
aa_class_weight <- c(aa_table[2]/nrow(aa_train_data), aa_table[1]/nrow(aa_train_data))
names(aa_class_weight) <- c("neg", "pos")
aa_case_weight <- numeric(length = nrow(aa_train_data))
aa_case_weight[aa_train_data$label == "pos"] <- aa_class_weight[2]
aa_case_weight[aa_train_data$label == "neg"] <- aa_class_weight[1]
set.seed(1234)
RF_partition2_chunk_classWght <- ranger(label ~ ., 
                                         data = aa_train_data, 
                                         num.trees = 1000,
                                         max.depth = 8,
                                         probability = TRUE, 
                                         importance  = "impurity",
                                        case.weights = aa_case_weight)








############################################################################################################################################
# evaluate aa_train_1 <- perf_eval_ROC_PRC(model_fit = fit,

######### partition 2 random:
aaOrigoritrain_pred <- RF_partition2_random_classWght$predictions
aatrain_pred <- RF_partition2_random_classWght$predictions
aatrain_pred[is.na(aatrain_pred)] <- 0.5
RF_partition2_random_classWght$predictions <- aatrain_pred



aa_tr_perf_random_caseweight <- perf_eval_ROC_PRC(model_fit = RF_partition2_random_classWght,
                                my_dataset = aa_train_data,
                                my_label = aa_train_data$label, 
                                train_perf = T,
                                file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_train_perf_partition2_random_casewght.png")

aa_tst_perf_random_caseweight <- perf_eval_ROC_PRC(model_fit = RF_partition2_random_classWght,
                                                  my_dataset = aa_test_data,
                                                  my_label = aa_test_data[, ncol(aa_test_data)], 
                                                  train_perf = F,
                                                  file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_test_perf_partition2_random_casewght.png")

aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_train_data), names(MM9_1kb_tiled_owner_3Mfilter))]
aa_tr_partition2_random_perclass <- perf_eval_ROC_PRC_perclass(model_fit = RF_partition2_random_classWght,
                                                my_dataset = aa_train_data,
                                                my_label = aa_train_data[, ncol(aa_train_data)], 
                                                file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_train_perf_partition2_random_perclass.png", 
                                                train_perf = T,
                                                class_name = aaclass_name)

aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_test_data), names(MM9_1kb_tiled_owner_3Mfilter))]
aa_tr_partition2_random_perclass <- perf_eval_ROC_PRC_perclass(model_fit = RF_partition2_random_classWght,
                                                               my_dataset = aa_test_data,
                                                               my_label = aa_test_data[, ncol(aa_test_data)], 
                                                               file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_test_perf_partition2_random_perclass.png", 
                                                               train_perf = F,
                                                               class_name = aaclass_name)

aaimp <- importance(RF_partition2_random_classWght)
names(aaimp) <- name_dic_p2[1:length(aaimp),1]
#names(aaimp)[sort(aaimp, decreasing = T, index.return = T)$ix[1:30]]
png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_varimp_partition2_random.png",    # create PNG for the heat map        
    width = 6*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,1), mar = c(4,9,4,4))

barplot(sort(aaimp, decreasing = T)[1:50], horiz = T, las =2)
dev.off()


aamy_pred <- predict(RF_partition2_random_classWght, data=aa_test_data)
aapreddd <-  aamy_pred$predictions[, 2]
names(aapreddd) <- rownames(aa_test_data)
aapreddd2 <-  RF_partition2_random_classWght$predictions[, 2]
names(aapreddd2) <- rownames(aa_train_data)
aapreddd <- c(aapreddd, aapreddd2)
aatst <- performance_viz_chromosome(all_tile_label=MM9_1kb_tiled_owner_labels_binary_3Mfilter,
                                    show_points = c("all"), 
                                    subset_points= c( "train", "test", "valid"),
                                    partition_dataset = Partition_2_dfs,
                                    predicted_pos_probilbilty = aapreddd,
                                    positive_thersh = 0.53,
                                    specific_chr = character(0),
                                    my_file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/performance_per_chr_random_partition2_thresh53.png",
                                    lncRNA_coor_df= lncRNA_chosen_gt1k_uniqTiles,
                                    chr_sizes_file="~/Documents/Shayan/BioInf/lncRNA/mm9_chr_sizes.txt")


######### partition 2 chunk:
aaOrigoritrain_pred <- RF_partition2_chunk_classWght$predictions
aatrain_pred <- RF_partition2_chunk_classWght$predictions
aatrain_pred[is.na(aatrain_pred)] <- 0.5
RF_partition2_chunk_classWght$predictions <- aatrain_pred

aa_tr_perf_chunk_caseweight <- perf_eval_ROC_PRC(model_fit = RF_partition2_chunk_classWght,
                                                  my_dataset = aa_train_data,
                                                  my_label = aa_train_data$label, 
                                                  train_perf = T,
                                                  file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_train_perf_partition2_chunk_casewght.png")

aa_tst_perf_chunk_caseweight <- perf_eval_ROC_PRC(model_fit = RF_partition2_chunk_classWght,
                                                   my_dataset = aa_test_data,
                                                   my_label = aa_test_data[, ncol(aa_test_data)], 
                                                   train_perf = F,
                                                   file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_test_perf_partition2_chunk_casewght.png")


aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_train_data), names(MM9_1kb_tiled_owner_3Mfilter))]
aa_tr_partition2_chunk_perclass <- perf_eval_ROC_PRC_perclass(model_fit = RF_partition2_chunk_classWght,
                                                               my_dataset = aa_train_data,
                                                               my_label = aa_train_data[, ncol(aa_train_data)], 
                                                               file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_train_perf_partition2_chunk_perclass.png", 
                                                               train_perf = T,
                                                               class_name = aaclass_name)

aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_test_data), names(MM9_1kb_tiled_owner_3Mfilter))]
aa_tst_partition2_chunk_perclass <- perf_eval_ROC_PRC_perclass(model_fit = RF_partition2_chunk_classWght,
                                                               my_dataset = aa_test_data,
                                                               my_label = aa_test_data[, ncol(aa_test_data)], 
                                                               file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_test_perf_partition2_chunk_perclass.png", 
                                                               train_perf = F,
                                                               class_name = aaclass_name)
aaimp <- importance(RF_partition2_chunk_classWght)
names(aaimp) <- name_dic_p2[1:length(aaimp),1]
#names(aaimp)[sort(aaimp, decreasing = T, index.return = T)$ix[1:30]]
png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_varimp_partition2_chunk.png",    # create PNG for the heat map        
    width = 6*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,1), mar = c(4,9,4,4))

barplot(sort(aaimp, decreasing = T)[1:50], horiz = T, las =2)
dev.off()


aamy_pred <- predict(RF_partition2_chunk_classWght, data=aa_test_data)
aapreddd <-  aamy_pred$predictions[, 2]
names(aapreddd) <- rownames(aa_test_data)
aapreddd2 <-  RF_partition2_chunk_classWght$predictions[, 2]
names(aapreddd2) <- rownames(aa_train_data)
aapreddd <- c(aapreddd, aapreddd2)
aatst <- performance_viz_chromosome(all_tile_label=MM9_1kb_tiled_owner_labels_binary_3Mfilter,
                                    show_points = c("all"), 
                                    subset_points= c( "train", "test", "valid"),
                                    partition_dataset = Partition_2_dfs,
                                    predicted_pos_probilbilty = aapreddd,
                                    positive_thersh = 0.53,
                                    specific_chr = character(0),
                                    my_file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/performance_per_chr_chunk_partition2_thresh53.png",
                                    lncRNA_coor_df= lncRNA_chosen_gt1k_uniqTiles,
                                    chr_sizes_file="~/Documents/Shayan/BioInf/lncRNA/mm9_chr_sizes.txt")

#############################################################################
# Run random forests on partition 2, only remove distance feature

# random
set.seed(1345)
aa_train_data2 <- aa_train_data
aa_train_data2 <- aa_train_data2[, -c(which(name_dic_p2[, 1] == "distance"))]
RF_partition2_random_classWght_noDist <- ranger(label ~ ., 
                                         data = aa_train_data2, 
                                         num.trees = 1000,
                                         max.depth = 8,
                                         probability = TRUE, 
                                         importance  = "impurity",
                                         case.weights = aa_case_weight)

#chunk

set.seed(1346)
aa_train_data2 <- aa_train_data
aa_train_data2 <- aa_train_data2[, -c(which(name_dic_p2[, 1] == "distance"))]
RF_partition2_chunk_classWght_noDist <- ranger(label ~ ., 
                                                data = aa_train_data2, 
                                                num.trees = 1000,
                                                max.depth = 8,
                                                probability = TRUE, 
                                                importance  = "impurity",
                                                case.weights = aa_case_weight)
##### evaluating results, partition2 -- no distance - random

aaOrigoritrain_pred2 <- RF_partition2_random_classWght_noDist$predictions
aatrain_pred <- RF_partition2_random_classWght_noDist$predictions
aatrain_pred[is.na(aatrain_pred)] <- 0.5
RF_partition2_random_classWght_noDist$predictions <- aatrain_pred



aa_tr_perf_random_caseweight <- perf_eval_ROC_PRC(model_fit = RF_partition2_random_classWght_noDist,
                                                  my_dataset = aa_train_data,
                                                  my_label = aa_train_data$label, 
                                                  train_perf = T,
                                                  file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_train_perf_partition2_random_casewght.png")

aa_tst_perf_random_caseweight <- perf_eval_ROC_PRC(model_fit = RF_partition2_random_classWght_noDist,
                                                   my_dataset = aa_test_data,
                                                   my_label = aa_test_data[, ncol(aa_test_data)], 
                                                   train_perf = F,
                                                   file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_test_perf_partition2_random_casewght.png")

aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_train_data), names(MM9_1kb_tiled_owner_3Mfilter))]
aa_tr_partition2_random_perclass <- perf_eval_ROC_PRC_perclass(model_fit = RF_partition2_random_classWght_noDist,
                                                               my_dataset = aa_train_data,
                                                               my_label = aa_train_data[, ncol(aa_train_data)], 
                                                               file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_train_perf_partition2_random_perclass.png", 
                                                               train_perf = T,
                                                               class_name = aaclass_name)

aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_test_data), names(MM9_1kb_tiled_owner_3Mfilter))]
aa_tr_partition2_random_perclass <- perf_eval_ROC_PRC_perclass(model_fit = RF_partition2_random_classWght_noDist,
                                                               my_dataset = aa_test_data,
                                                               my_label = aa_test_data[, ncol(aa_test_data)], 
                                                               file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_test_perf_partition2_random_perclass.png", 
                                                               train_perf = F,
                                                               class_name = aaclass_name)

aaimp <- importance(RF_partition2_random_classWght_noDist)
names(aaimp) <- name_dic_p2[1:length(aaimp),1]
#names(aaimp)[sort(aaimp, decreasing = T, index.return = T)$ix[1:30]]
png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_varimp_partition2_random.png",    # create PNG for the heat map        
    width = 6*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,1), mar = c(4,9,4,4))

barplot(sort(aaimp, decreasing = T)[1:50], horiz = T, las =2)
dev.off()


aamy_pred <- predict(RF_partition2_random_classWght_noDist, data=aa_test_data)
aapreddd <-  aamy_pred$predictions[, 2]
names(aapreddd) <- rownames(aa_test_data)
aapreddd2 <-  RF_partition2_random_classWght_noDist$predictions[, 2]
names(aapreddd2) <- rownames(aa_train_data)
aapreddd <- c(aapreddd, aapreddd2)
aatst <- performance_viz_chromosome(all_tile_label=MM9_1kb_tiled_owner_labels_binary_3Mfilter,
                                    show_points = c("all"), 
                                    subset_points= c( "train", "test", "valid"),
                                    partition_dataset = Partition_2_dfs,
                                    predicted_pos_probilbilty = aapreddd,
                                    positive_thersh = 0.53,
                                    specific_chr = character(0),
                                    my_file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/performance_per_chr_random_partition2_thresh53.png",
                                    lncRNA_coor_df= lncRNA_chosen_gt1k_uniqTiles,
                                    chr_sizes_file="~/Documents/Shayan/BioInf/lncRNA/mm9_chr_sizes.txt")


######### partition 2 chunk -- nodistance:
aaOrigoritrain_pred <- RF_partition2_chunk_classWght$predictions
aatrain_pred <- RF_partition2_chunk_classWght$predictions
aatrain_pred[is.na(aatrain_pred)] <- 0.5
RF_partition2_chunk_classWght$predictions <- aatrain_pred

aa_tr_perf_chunk_caseweight <- perf_eval_ROC_PRC(model_fit = RF_partition2_chunk_classWght,
                                                 my_dataset = aa_train_data,
                                                 my_label = aa_train_data$label, 
                                                 train_perf = T,
                                                 file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_train_perf_partition2_chunk_casewght.png")

aa_tst_perf_chunk_caseweight <- perf_eval_ROC_PRC(model_fit = RF_partition2_chunk_classWght,
                                                  my_dataset = aa_test_data,
                                                  my_label = aa_test_data[, ncol(aa_test_data)], 
                                                  train_perf = F,
                                                  file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_test_perf_partition2_chunk_casewght.png")


aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_train_data), names(MM9_1kb_tiled_owner_3Mfilter))]
aa_tr_partition2_chunk_perclass <- perf_eval_ROC_PRC_perclass(model_fit = RF_partition2_chunk_classWght,
                                                              my_dataset = aa_train_data,
                                                              my_label = aa_train_data[, ncol(aa_train_data)], 
                                                              file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_train_perf_partition2_chunk_perclass.png", 
                                                              train_perf = T,
                                                              class_name = aaclass_name)

aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_test_data), names(MM9_1kb_tiled_owner_3Mfilter))]
aa_tst_partition2_chunk_perclass <- perf_eval_ROC_PRC_perclass(model_fit = RF_partition2_chunk_classWght,
                                                               my_dataset = aa_test_data,
                                                               my_label = aa_test_data[, ncol(aa_test_data)], 
                                                               file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_test_perf_partition2_chunk_perclass.png", 
                                                               train_perf = F,
                                                               class_name = aaclass_name)
aaimp <- importance(RF_partition2_chunk_classWght)
names(aaimp) <- name_dic_p2[1:length(aaimp),1]
#names(aaimp)[sort(aaimp, decreasing = T, index.return = T)$ix[1:30]]
png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_varimp_partition2_chunk.png",    # create PNG for the heat map        
    width = 6*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,1), mar = c(4,9,4,4))

barplot(sort(aaimp, decreasing = T)[1:50], horiz = T, las =2)
dev.off()


aamy_pred <- predict(RF_partition2_chunk_classWght, data=aa_test_data)
aapreddd <-  aamy_pred$predictions[, 2]
names(aapreddd) <- rownames(aa_test_data)
aapreddd2 <-  RF_partition2_chunk_classWght$predictions[, 2]
names(aapreddd2) <- rownames(aa_train_data)
aapreddd <- c(aapreddd, aapreddd2)
aatst <- performance_viz_chromosome(all_tile_label=MM9_1kb_tiled_owner_labels_binary_3Mfilter,
                                    show_points = c("all"), 
                                    subset_points= c( "train", "test", "valid"),
                                    partition_dataset = Partition_2_dfs,
                                    predicted_pos_probilbilty = aapreddd,
                                    positive_thersh = 0.53,
                                    specific_chr = character(0),
                                    my_file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/performance_per_chr_chunk_partition2_thresh53.png",
                                    lncRNA_coor_df= lncRNA_chosen_gt1k_uniqTiles,
                                    chr_sizes_file="~/Documents/Shayan/BioInf/lncRNA/mm9_chr_sizes.txt")

#############################################################################
# evaluating results of models with up/down sampling or cost sensitive on P4 and P5, using multimodel viz
###################################
aa_folder_names <- list.dirs("~/Documents/Shayan/BioInf/lncRNA", full.names = T)
aa_folder_names <- aa_folder_names[grep("Learned_models", aa_folder_names)]
names(aa_folder_names) <- c("bal_cost_wght", "down_samp", "up_samp")
#Partition_4_modified_dataset_and_removedColumnNames___chunk___distance___RFmodel

aa_ch_rand <- c("rand", "chunk")
aa_dis_full <- c("full", "distance")
aa_file_names_pre <- "Partition_4_modified_dataset_and_removedColumnNames___"
# partition 4
#load("~/Documents/Shayan/BioInf/lncRNA/Partition_4_modified_dataset_and_removedColumnNames.RData")
################# rand, chunk
#test, train
aaplot_strore_address <- "~/Documents/Shayan/BioInf/lncRNA/Performance_multi/"
load("~/Documents/Shayan/BioInf/lncRNA/MM9_1kb_tiled_owner_3Mfilter.nosync.RData")
aa_pr_multimodel<- list()
aa_pr_multimodeltr <- list()
aa_pr_multimodel_multiclass <- list()
aa_pr_multimodel_multiclass_tr <- list()
for(aa_par in 1:length(aa_ch_rand)){
  print(aa_ch_rand[aa_par])
  aa_cnt <- 1
  my_model_fit_list <- list()
  if(aa_ch_rand[aa_par]  == "rand"){
    my_partition_df <- Partition_4_dfs
  }else{
    my_partition_df <- Partition_4_chunk_dfs
  }
  aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
  aa_train_data <- Partition_4_modified_dataset[aatrainind,]
  aa_test_data <- Partition_4_modified_dataset[-aatrainind,]
  
  for(aa_fl in 1:length(aa_folder_names)){
    print(aa_folder_names[aa_fl])
    for(aa_feat in 1:length(aa_dis_full)){
      print(aa_dis_full[aa_feat])
      aaaname <- paste0(aa_folder_names[aa_fl],"/",  aa_file_names_pre, aa_ch_rand[aa_par], "___", aa_dis_full[aa_feat], "___RFmodel.RData")
      load(aaaname)
      my_model_fit_list[[aa_cnt]] <- my_RF_model
      names(my_model_fit_list)[aa_cnt] <- paste0(names(aa_folder_names)[aa_fl], "__", aa_dis_full[aa_feat])
      aa_cnt <- aa_cnt + 1
      
    }
  }
  aanampas <- paste0(names(my_model_fit_list), collapse = "__vs__")
  print("plotting test all ...")
  aa_pr_multimodel[[aa_par]] <- perf_eval_ROC_PRC_multimodel(model_fit_list = my_model_fit_list,
                                           my_dataset = aa_test_data,
                                           my_label = aa_test_data[, ncol(aa_test_data)],
                                           file_name = paste0(aaplot_strore_address, aanampas,"__",aa_ch_rand[aa_par], "___test_all_P4.png"),
                                           train_perf = F,
                                           neg_value = "neg",
                                           pos_value = "pos")
  print("plotting train all ...")
  aa_pr_multimodeltr[[aa_par]] <- perf_eval_ROC_PRC_multimodel(model_fit_list = my_model_fit_list,
                                                   my_dataset = aa_train_data,
                                                   my_label = aa_train_data[, ncol(aa_train_data)],
                                                   file_name = paste0(aaplot_strore_address, aanampas,"__",aa_ch_rand[aa_par], "___train_all_P4.png"),
                                                   train_perf = T,
                                                   neg_value = "neg",
                                                   pos_value = "pos")
  aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_test_data), names(MM9_1kb_tiled_owner_3Mfilter))]
  print("plotting test perclass ...")
  aa_pr_multimodel_multiclass[[aa_par]] <- perf_eval_ROC_PRC_perclass_multimodel(model_fit_list = my_model_fit_list,
                                                            my_dataset = aa_test_data,
                                                            my_label = aa_test_data[, ncol(aa_test_data)],
                                                            file_name = paste0(aaplot_strore_address, aanampas,"__",aa_ch_rand[aa_par], "___test_all_P4_perClass.png"),
                                                            train_perf = F,
                                                            class_name = aaclass_name,
                                                            neg_value = "neg",
                                                            pos_value = "pos")
  aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_train_data), names(MM9_1kb_tiled_owner_3Mfilter))]
  print("plotting train perclass ...")
  aa_pr_multimodel_multiclass_tr[[aa_par]] <- perf_eval_ROC_PRC_perclass_multimodel(model_fit_list = my_model_fit_list,
                                                                                    my_dataset = aa_train_data,
                                                                                    my_label = aa_train_data[, ncol(aa_train_data)],
                                                                                    file_name = paste0(aaplot_strore_address, aanampas,"__",aa_ch_rand[aa_par], "___train_all_P4_perClass.png"),
                                                                                    train_perf = T,
                                                                                    class_name = aaclass_name,
                                                                                    neg_value = "neg",
                                                                                    pos_value = "pos")
  
  
}



###################################
# partition 5
#load("~/Documents/Shayan/BioInf/lncRNA/Partition_5_modified_dataset_and_removedColumnNames.RData")
################# rand
aa_folder_names <- list.dirs("~/Documents/Shayan/BioInf/lncRNA", full.names = T)
aa_folder_names <- aa_folder_names[grep("Learned_models", aa_folder_names)]
names(aa_folder_names) <- c("bal_cost_wght", "down_samp", "up_samp")
#Partition_5_modified_dataset_and_removedColumnNames___chunk___distance___RFmodel

aa_ch_rand <- c("rand", "chunk")
aa_dis_full <- c("full", "distance")
aa_file_names_pre <- "Partition_5_modified_dataset_and_removedColumnNames___"

aaplot_strore_address <- "~/Documents/Shayan/BioInf/lncRNA/Performance_multi/"
load("~/Documents/Shayan/BioInf/lncRNA/MM9_1kb_tiled_owner_3Mfilter.nosync.RData")
aa_pr_multimodel_P5<- list()
aa_pr_multimodeltr_P5 <- list()
aa_pr_multimodel_multiclass_P5 <- list()
aa_pr_multimodel_multiclass_tr_P5 <- list()

for(aa_par in 1:length(aa_ch_rand)){
  print(aa_ch_rand[aa_par])
  aa_cnt <- 1
  my_model_fit_list <- list()
  if(aa_ch_rand[aa_par]  == "rand"){
    my_partition_df <- Partition_5_dfs
  }else{
    my_partition_df <- Partition_5_chunk_dfs
  }
  aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
  aa_train_data <- Partition_5_modified_dataset[aatrainind,]
  aa_test_data <- Partition_5_modified_dataset[-aatrainind,]
  
  for(aa_fl in 1:length(aa_folder_names)){
    print(aa_folder_names[aa_fl])
    for(aa_feat in 1:length(aa_dis_full)){
      print(aa_dis_full[aa_feat])
      aaaname <- paste0(aa_folder_names[aa_fl],"/",  aa_file_names_pre, aa_ch_rand[aa_par], "___", aa_dis_full[aa_feat], "___RFmodel.RData")
      
      if(aaaname == "/Users/Shayan/Documents/Shayan/BioInf/lncRNA/Learned_models_downSample.nosync/Partition_5_modified_dataset_and_removedColumnNames___chunk___full___RFmodel.RData"){
        my_model_fit_list[[aa_cnt]] <-  my_model_fit_list[[aa_cnt - 2]]
      }else{
        load(aaaname)
        my_model_fit_list[[aa_cnt]] <- my_RF_model
      }
      
      names(my_model_fit_list)[aa_cnt] <- paste0(names(aa_folder_names)[aa_fl], "__", aa_dis_full[aa_feat])
      aa_cnt <- aa_cnt + 1
      
    }
  }
  aanampas <- paste0(names(my_model_fit_list), collapse = "__vs__")
  print("plotting test all ...")
  aa_pr_multimodel_P5[[aa_par]] <- perf_eval_ROC_PRC_multimodel(model_fit_list = my_model_fit_list,
                                                             my_dataset = aa_test_data,
                                                             my_label = aa_test_data[, ncol(aa_test_data)],
                                                             file_name = paste0(aaplot_strore_address, aanampas,"__",aa_ch_rand[aa_par], "___test_all_P5.png"),
                                                             train_perf = F,
                                                             neg_value = "neg",
                                                             pos_value = "pos")
  print("plotting train all ...")
  aa_pr_multimodeltr_P5[[aa_par]] <- perf_eval_ROC_PRC_multimodel(model_fit_list = my_model_fit_list,
                                                               my_dataset = aa_train_data,
                                                               my_label = aa_train_data[, ncol(aa_train_data)],
                                                               file_name = paste0(aaplot_strore_address, aanampas,"__",aa_ch_rand[aa_par], "___train_all_P5.png"),
                                                               train_perf = T,
                                                               neg_value = "neg",
                                                               pos_value = "pos")
  aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_test_data), names(MM9_1kb_tiled_owner_3Mfilter))]
  print("plotting test perclass ...")
  aa_pr_multimodel_multiclass_P5[[aa_par]] <- perf_eval_ROC_PRC_perclass_multimodel(model_fit_list = my_model_fit_list,
                                                                                 my_dataset = aa_test_data,
                                                                                 my_label = aa_test_data[, ncol(aa_test_data)],
                                                                                 file_name = paste0(aaplot_strore_address, aanampas,"__",aa_ch_rand[aa_par], "___test_all_P5_perClass.png"),
                                                                                 train_perf = F,
                                                                                 class_name = aaclass_name,
                                                                                 neg_value = "neg",
                                                                                 pos_value = "pos")
  aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_train_data), names(MM9_1kb_tiled_owner_3Mfilter))]
  print("plotting train perclass ...")
  aa_pr_multimodel_multiclass_tr_P5[[aa_par]] <- perf_eval_ROC_PRC_perclass_multimodel(model_fit_list = my_model_fit_list,
                                                                                    my_dataset = aa_train_data,
                                                                                    my_label = aa_train_data[, ncol(aa_train_data)],
                                                                                    file_name = paste0(aaplot_strore_address, aanampas,"__",aa_ch_rand[aa_par], "___train_all_P5_perClass.png"),
                                                                                    train_perf = T,
                                                                                    class_name = aaclass_name,
                                                                                    neg_value = "neg",
                                                                                    pos_value = "pos")
  
  
}





#############################################################################
# add new features on the number of binding sites and so on and check performance





#######
# create a small test dataset for setting up RF on hal
setupTest_dataset
aa_feature <- matrix(rnorm(30000), nrow = 1000)
aa_label <- sample(x = c("pos", "neg"), size = 1000,replace = T, prob = c(0.1, 0.9))

setupTest_dataset <- cbind(as.data.frame(aa_feature), aa_label)
colnames(setupTest_dataset)[1:(ncol(setupTest_dataset) - 1)] <- paste0("feature_", c(1:(ncol(setupTest_dataset) - 1)))
colnames(setupTest_dataset)[ncol(setupTest_dataset)] <- "label"
setupTest_dataset$label[1:5]
rownames(setupTest_dataset) <- sample(names(MM9_1kb_tiled_owner_labels_binary_3Mfilter), nrow(setupTest_dataset))
setupTest_index <- data.frame(owner = sample(Partition_2_chunk_dfs$owner, nrow(setupTest_dataset)),
                              label = setupTest_dataset$label,
                              dataset = sample(x = c("train", "valid", "test"),
                                               size =  nrow(setupTest_dataset),
                                               replace = T,
                                               prob = c(0.7,0.2,0.1)),
                              tile_name = rownames(setupTest_dataset))
# save(list = c("setupTest_index", "setupTest_dataset", "MM9_1kb_tiled_GR_filtered", "name_dic_p2"), 
#      file = "~/Documents/Shayan/BioInf/lncRNA/setupTest_dataset.RData")
# Partition_2_modified_dataset <- setupTest_dataset
# Partition_2_chunk_dfs <- setupTest_index
# Partition_2_dfs  <- setupTest_index
# MM9_1kb_tiled_GR_filtered
# name_dic_p2 <- cbind(colnames(setupTest_dataset), colnames(setupTest_dataset))
##########################
# take a look at kmers

##########################
# read the per lncRNA models one by one, write a new dataset using only the top 100 predictors
my_model_list_p6_rand <- list()
my_model_list_p6_chunk <- list()
my_model_feature_varImp_rand <- list()
my_model_feature_varImp_chunk <- list()
my_name_dic_list_p6 <- list()

aafiles <- list.files("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_tripCorr/")
aafiles2 <- list.files("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/")
aafiles2_spl <- unlist(lapply(strsplit(aafiles2, "\\."), "[[", 1))
for(i in 1:length(aafiles2_spl)){
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_tripCorr/", aafiles2_spl[i], "___rand___full___RFmodel.RData"))
  my_model_list_p6_rand[[i]] <- my_RF_model
  variable_importance <- importance(my_RF_model)
  my_model_feature_varImp_rand[[i]] <- variable_importance
  my_top_100 <- sort(variable_importance, decreasing = T, index.return = T)$ix[1:100]
  save(list = c("my_top_100"), file = paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_top_index/", aafiles2_spl[i],"___rand__topind.RData" ))
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_tripCorr/", aafiles2_spl[i], "___chunk___full___RFmodel.RData"))
  my_model_list_p6_chunk[[i]] <- my_RF_model
  variable_importance <- importance(my_RF_model)
  my_model_feature_varImp_chunk[[i]] <- variable_importance
  my_top_100 <- sort(variable_importance, decreasing = T, index.return = T)$ix[1:100]
  save(list = c("my_top_100"), file = paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_top_index/", aafiles2_spl[i],"___chunk__topind.RData" ))
  
}
names(my_model_list_p6_chunk) <- aafiles2_spl
names(my_model_feature_varImp_chunk) <- aafiles2_spl
names(my_model_list_p6_rand) <- aafiles2_spl
names(my_model_feature_varImp_rand) <- aafiles2_spl


# read the name_dics for all
for(i in 1:length(aafiles2_spl)){
  print(i)
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/", aafiles2[i]))
  my_name_dic_list_p6[[i]] <-  my_name_dic
}
names(my_name_dic_list_p6) <- aafiles2_spl
# write jobs for the 

# write job for the ones running with only the top100
aa_un_own <- unique(Partition_6_dfs$owner)
for(i in 1:length(aa_un_own)){
  cat(c("Rscript --vanilla RF_load_run_p6_subset.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_",
        aa_un_own[i],
        ".RData rand ", 
        sample(c(100:10000), 1),
        " /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_top_index/", 
        aafiles2_spl[grep(pattern = aa_un_own[i], x = aafiles2_spl)],
        "___rand__topind.RData", "\n"),
      sep = "", file = "run_rf_p6_t100.job", append = T)
  
  
  
  cat(c("Rscript --vanilla RF_load_run_p6_subset.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_",
        aa_un_own[i],
        ".RData chunk ",
        sample(c(100:10000), 1),
        " /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_top_index/",
        aafiles2_spl[grep(pattern = aa_un_own[i], x = aafiles2_spl)],
        "___chunk__topind.RData", "\n"),
      sep = "", file = "run_rf_p6_t100.job", append = T)
}

###### visualize multimodel performance for each lncRNA
aa_file_mod_t100 <- list.files("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_subset")
aa_file_mod_all <- list.files("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_tripCorr/")
aafiles2 <- list.files("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/")
aafiles2_spl <- unlist(lapply(strsplit(aafiles2, "\\."), "[[", 1))
my_model_list_p6_rand_top100 <- list()
my_model_list_p6_chunk_top100 <- list()
aafiles2_spl2 <- unlist(lapply(strsplit(aafiles2_spl, "_"), "[[", 5))
aaplot_strore_address <- "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_comp/"
aa_pr_multimodel_P6_rand <- list()
aa_pr_multimodeltr_P6_rand <- list()
aa_pr_multimodel_P6_chunk <- list()
aa_pr_multimodeltr_P6_chunk <- list()

source("~/Documents/Shayan/BioInf/lncRNA/RF_evaluation_functions.R")
for(i in 1:length(aafiles2_spl)){
  print(i)
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/", aafiles2[i]))
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_subset/", aafiles2_spl[i], "___rand___full___top100RFmodel.RData"))
  my_partition_df <- my_partition_rand
  aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
  aa_train_data <- my_Dataset[aatrainind,]
  aa_test_data <- my_Dataset[-aatrainind,]
  my_model_list_p6_rand_top100[[i]] <- my_RF_model
  my_model_fit_list <- list()
  my_model_fit_list[[1]] <- my_model_list_p6_rand[[i]]
  my_model_fit_list[[2]] <- my_RF_model
  names(my_model_fit_list) <- c("all_features", "top100")
  
  aa_pr_multimodel_P6_rand[[i]] <- perf_eval_ROC_PRC_multimodel(model_fit_list = my_model_fit_list,
                                                                     my_dataset = aa_test_data,
                                                                     my_label = aa_test_data[, ncol(aa_test_data)],
                                                                     file_name = paste0(aaplot_strore_address, aafiles2_spl2[i],"__", "rand___test_all_P6.png"),
                                                                     train_perf = F,
                                                                     neg_value = "neg",
                                                                     pos_value = "pos")
  aa_pr_multimodeltr_P6_rand[[i]] <- perf_eval_ROC_PRC_multimodel(model_fit_list = my_model_fit_list,
                                                                       my_dataset = aa_train_data,
                                                                       my_label = aa_train_data[, ncol(aa_train_data)],
                                                                       file_name = paste0(aaplot_strore_address, aafiles2_spl2[i],"__", "rand___train_all_P6.png"),
                                                                       train_perf = T,
                                                                       neg_value = "neg",
                                                                       pos_value = "pos")
  
  
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_subset/", aafiles2_spl[i], "___chunk___full___top100RFmodel.RData"))
  my_partition_df <- my_partition_chunk
  aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
  aa_train_data <- my_Dataset[aatrainind,]
  aa_test_data <- my_Dataset[-aatrainind,]
  
  my_model_list_p6_chunk_top100[[i]] <- my_RF_model
  my_model_fit_list <- list()
  my_model_fit_list[[1]] <- my_model_list_p6_chunk[[i]]
  my_model_fit_list[[2]] <- my_RF_model
  names(my_model_fit_list) <- c("all_features", "top100")
  aa_pr_multimodel_P6_chunk[[i]] <- perf_eval_ROC_PRC_multimodel(model_fit_list = my_model_fit_list,
                                                                      my_dataset = aa_test_data,
                                                                      my_label = aa_test_data[, ncol(aa_test_data)],
                                                                      file_name = paste0(aaplot_strore_address,  aafiles2_spl2[i],"__", "chunk___test_all_P6.png"),
                                                                      train_perf = F,
                                                                      neg_value = "neg",
                                                                      pos_value = "pos")
  aa_pr_multimodeltr_P6_chunk[[i]] <- perf_eval_ROC_PRC_multimodel(model_fit_list = my_model_fit_list,
                                                                        my_dataset = aa_train_data,
                                                                        my_label = aa_train_data[, ncol(aa_train_data)],
                                                                        file_name = paste0(aaplot_strore_address, aafiles2_spl2[i],"__", "chunk___train_all_P6.png"),
                                                                        train_perf = T,
                                                                        neg_value = "neg",
                                                                        pos_value = "pos")
  
  
}
# library(ranger)
# load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_tripCorr/Partition_6_modified_dataset_Gm14820___chunk___full___RFmodel.RData")
# load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Gm14820.RData")
# aaimp <- importance_pvalues(x= my_RF_model, method = "altmann", formula = label ~ ., data = my_Dataset[my_partition_chunk$dataset %in% c("train", "valid"),] )
# 
# View(aaimp)



my_model_feature_varImp_rand

my_name_dic_list_p6


aa_all_f <- sort(unique(do.call(rbind, my_name_dic_list_p6)[, 1]))
aa_feat_emp <- matrix(nrow = length(my_name_dic_list_p6), ncol = length(aa_all_f))
colnames(aa_feat_emp) <- aa_all_f
rownames(aa_feat_emp) <- unlist(lapply(strsplit(names(my_model_list_p6_rand), "_"), "[[", 5))

aa_feat_emp_rank <- matrix(nrow = length(my_name_dic_list_p6), ncol = length(aa_all_f))
colnames(aa_feat_emp_rank) <- aa_all_f
rownames(aa_feat_emp_rank) <- unlist(lapply(strsplit(names(my_model_list_p6_rand), "_"), "[[", 5))
for(i in 1:length(my_model_feature_varImp_rand)){
  print(i)
  acurimp <- my_model_feature_varImp_rand[[i]]
  names(acurimp) <- my_name_dic_list_p6[[i]][match(names(my_model_feature_varImp_rand[[i]])
                                                   ,my_name_dic_list_p6[[i]][, 2]), 1]
  acurimp_rank <- sort(acurimp, decreasing = T, index.return = T)$ix
  aa_feat_emp[i, match(names(acurimp), colnames(aa_feat_emp))] <- acurimp
  aa_feat_emp_rank[i, match(names(acurimp), colnames(aa_feat_emp))] <- acurimp_rank
}

aa_feat_emp[is.na(aa_feat_emp)] <- 0
aa_feat_emp2 <- aa_feat_emp
aa_feat_emp2[aa_feat_emp <= 0] <- 0

png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/feature_importance_heatmap_p6.png", 
    width = 20*300, 
    height = 10*300, 
    res = 300,
    pointsize = 10)
heatmap.2(aa_feat_emp2, Rowv = T, Colv = T,
          dendrogram = "both", na.rm = T, trace = "none",
          breaks = seq(from = range(aa_feat_emp2)[1], range(aa_feat_emp2)[2], length.out = 11),
          col = colorspace::sequential_hcl(n=10)
          )
dev.off()
png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/feature_importance_sum_barplot_p6.png", 
    width = 20*300, 
    height = 4*300, 
    res = 300,
    pointsize = 8)
par(mar = c(10,4,4,4))
barplot(sort(colSums(aa_feat_emp2), decreasing = T)[1:100], las = 2)
dev.off()

sum(!is.na(aa_feat_emp_rank))
aa_feat_emp_rank[is.na(aa_feat_emp_rank)] <- range(aa_feat_emp_rank, na.rm = T)[2]

heatmap.2(aa_feat_emp_rank, Rowv = T, Colv = T,
          dendrogram = "both", trace = "none",
          breaks = c(1, 10, 50, 100, 200, 500, 1000, 1500),
          col = colorspace::sequential_hcl(n=7)
)

#aavar <- apply(aa_feat_emp_rank, 2, sd)
#aa_feat_emp_rank <- aa_feat_emp_rank[, -which(aavar == 0)]

z <- zClust(x=aa_feat_emp_rank, scale="none", zlim=c(-3,3), method="average")
distCor(aa_feat_emp_rank)
# require(RColorBrewer)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(7)



aarbp <- grep(pattern = "*_mm9", x = colnames(aa_feat_emp_rank))
aachrom <- grep(pattern = '^H[0-9]{1,2}', x = colnames(aa_feat_emp_rank))
aatf <- grep(pattern = '*.bed', x = colnames(aa_feat_emp_rank))
aatf <- setdiff(aatf, aachrom)
aa_col <- character(ncol(aa_feat_emp_rank))
aa_col[] <- "grey"
aa_col[aarbp] <- "red"
aa_col[aachrom] <- "blue"
aa_col[aatf] <- "green"

png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/feature_importance_ranking_heatmap_p6.png", 
    width = 30*300, 
    height = 5*300, 
    res = 300,
    pointsize = 10)
heatmap.2(z$data,
          trace='none',
          breaks = c(0, 10, 50, 100, 200, 500, 1000, 1500),
          col=rev(cols),
          Rowv=z$Rowv,
          Colv=z$Colv, margins = c(8, 8), 
          ColSideColors =aa_col )
dev.off()

####################################################################################################################################################
# write a function to visualize top features, given the model and dataset


#Partition_6_chunk_dfs

aa_file_datast <- list.files("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected")
aa_file_model  <- list.files("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_permutation")
aaspl_datast <- unlist(lapply(strsplit(unlist(lapply(strsplit(aa_file_datast, "\\."), "[[", 1)), "_"), "[[", 5))

for(i in 1:length(aa_file_datast)){
  print(i)
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/", aa_file_datast[i]))
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_permutation/Partition_6_modified_dataset_", aaspl_datast[i], "___rand___full___RFmodel.RData"))
  AAtst <- visualize_topN_scatter(my_dataset = my_Dataset,
                                  my_partition = my_partition_rand,
                                  my_model = my_RF_model,
                                  topN=20, 
                                  my_file_name = paste0("~/Documents/Shayan/BioInf/lncRNA/plots/Feature_top20/", aaspl_datast[i], "_top20.png"), 
                                  box_only = T)
  
}

####################################################################################################################################################

# plotting Malat1 Performance after removing each feature class

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Malat1.RData")
my_partition_df <- my_partition_rand
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_MALAT1_mlist <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_tripCorr/Partition_6_modified_dataset_Malat1___rand___full___RFmodel.RData")
aa_MALAT1_mlist[[1]] <- my_RF_model
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_MALAT1/Partition_6_modified_dataset_Malat1___rand___RNA_polymerase_II__315.bed__RNAseq__CAGE___RFmodel.RData")
aa_MALAT1_mlist[[2]] <- my_RF_model
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_MALAT1/Partition_6_modified_dataset_Malat1___rand___first_229_last_1489__noexp_noChromPair___RFmodel.RData")
aa_MALAT1_mlist[[3]] <- my_RF_model
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_MALAT1/Partition_6_modified_dataset_Malat1___rand___first_229_last_1488___noexp_noChrom_noBed___RFmodel.RData")
aa_MALAT1_mlist[[4]] <- my_RF_model
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_MALAT1/Partition_6_modified_dataset_Malat1___rand___first_229_last_1487___RFmodel.RData")
aa_MALAT1_mlist[[5]] <- my_RF_model
aaplot_strore_address <- "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_MALAT1_comparison/"

names(aa_MALAT1_mlist) <- c("All", "noExp", "noExp_noChrom","noExp_noChrom_noTF", "noExp_noChrom_noTF_noRBP" )
aa_pr_multimodel_P6_Malat1_rand <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_MALAT1_mlist[1:2],
                                                               my_dataset = aa_test_data,
                                                               my_label = aa_test_data[, ncol(aa_test_data)],
                                                               file_name = paste0(aaplot_strore_address,"All_vs_noExp","__rand___test_all_P6.png"),
                                                               train_perf = F,
                                                               neg_value = "neg",
                                                               pos_value = "pos")

aa_pr_multimodel_P6_Malat1_rand <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_MALAT1_mlist[1:3],
                                                                my_dataset = aa_test_data,
                                                                my_label = aa_test_data[, ncol(aa_test_data)],
                                                                file_name = paste0(aaplot_strore_address,"All_vs_noExp_vs_noExpnoChrom","__rand___test_all_P6.png"),
                                                                train_perf = F,
                                                                neg_value = "neg",
                                                                pos_value = "pos")

aa_pr_multimodel_P6_Malat1_rand <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_MALAT1_mlist[1:4],
                                                                my_dataset = aa_test_data,
                                                                my_label = aa_test_data[, ncol(aa_test_data)],
                                                                file_name = paste0(aaplot_strore_address,"All_vs_noExp_vs_noExpnoChrom_vs_noExpnoChromnoTF","__rand___test_all_P6.png"),
                                                                train_perf = F,
                                                                neg_value = "neg",
                                                                pos_value = "pos")
aa_pr_multimodel_P6_Malat1_rand <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_MALAT1_mlist,
                                                                my_dataset = aa_test_data,
                                                                my_label = aa_test_data[, ncol(aa_test_data)],
                                                                file_name = paste0(aaplot_strore_address,"All_vs_noExp_vs_noExpnoChrom_vs_noExpnoChromnoTF_vs_noExpnoChromnoTFnoRBP","__rand___test_all_P6.png"),
                                                                train_perf = F,
                                                                neg_value = "neg",
                                                                pos_value = "pos")

##################### aligning top 100 5mers that come out of sequence only version of MALAT1 model 
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_MALAT1/Partition_6_modified_dataset_Malat1___rand___first_229_last_1487___RFmodel.RData")
aa_malat_1_SeqOnly_variable_importance <- importance(my_RF_model)
names(aa_malat_1_SeqOnly_variable_importance) <- my_name_dic[match(names(aa_malat_1_SeqOnly_variable_importance) , my_name_dic[, 2]),1]
aa_myKmer <- setdiff(names(sort(aa_malat_1_SeqOnly_variable_importance, decreasing = T))[1:104], c("triplex_TC", "triplex", "triplex_GT", "triplex_GA"))
for(i in 1:length(aa_myKmer)){
  cat(c(">seq", i, "\n"), file = "~/Documents/Shayan/BioInf/lncRNA/malat1_top100_kmer.txt", sep = "", append = T)
  cat(c(aa_myKmer[i], "NNNNN", "\n"), file = "~/Documents/Shayan/BioInf/lncRNA/malat1_top100_kmer.txt", sep = "", append = T)
}

####################################################################################################################################################

# plotting Malat1, NEAT1, Gm14820, Trerf1, Kcnq1ot1 Performance after removing each feature class

aaplot_strore_address <- "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison/"

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Malat1.RData")
my_partition_df <- my_partition_rand
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Malat1_mlist_nameDic <- my_name_dic
aa_Malat1_mlist_rand_featImp <- list()
aa_Malat1_mlist_rand <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_tripCorr/Partition_6_modified_dataset_Malat1___rand___full___RFmodel.RData")
aa_Malat1_mlist_rand[[1]] <- my_RF_model
aa_Malat1_mlist_rand_featImp[[1]] <- importance(my_RF_model)
names(aa_Malat1_mlist_rand_featImp[[1]]) <- my_name_dic[match(names(aa_Malat1_mlist_rand_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_MALAT1/Partition_6_modified_dataset_Malat1___rand___first_229_last_1488___noexp_noChrom_noBed___RFmodel.RData")
aa_Malat1_mlist_rand[[2]] <- my_RF_model
aa_Malat1_mlist_rand_featImp[[2]] <- importance(my_RF_model)
names(aa_Malat1_mlist_rand_featImp[[2]]) <- my_name_dic[match(names(aa_Malat1_mlist_rand_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_MALAT1/Partition_6_modified_dataset_Malat1___rand___first_229_last_1492___RFmodel.RData")
aa_Malat1_mlist_rand[[3]] <- my_RF_model
aa_Malat1_mlist_rand_featImp[[3]] <- importance(my_RF_model)
names(aa_Malat1_mlist_rand_featImp[[3]]) <- my_name_dic[match(names(aa_Malat1_mlist_rand_featImp[[3]]), my_name_dic[, 2]),1]

names(aa_Malat1_mlist_rand) <- c("All", "seq_only", "kmer_only")

aa_pr_multimodel_P6_Malat1_rand <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Malat1_mlist_rand,
                                                               my_dataset = aa_test_data,
                                                               my_label = aa_test_data[, ncol(aa_test_data)],
                                                               file_name = paste0(aaplot_strore_address,"Malat1_All_vs_seq_vs_kmer","__rand___test_all_P6.png"),
                                                               train_perf = F,
                                                               neg_value = "neg",
                                                               pos_value = "pos")

my_partition_df <- my_partition_chunk
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Malat1_mlist_chunk <- list()
aa_Malat1_mlist_chunk_featImp <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_tripCorr/Partition_6_modified_dataset_Malat1___chunk___full___RFmodel.RData")
aa_Malat1_mlist_chunk[[1]] <- my_RF_model
aa_Malat1_mlist_chunk_featImp[[1]] <- importance(my_RF_model)
names(aa_Malat1_mlist_chunk_featImp[[1]]) <- my_name_dic[match(names(aa_Malat1_mlist_chunk_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_MALAT1/Partition_6_modified_dataset_Malat1___chunk___first_229_last_1488___noexp_noChrom_noBed___RFmodel.RData")
aa_Malat1_mlist_chunk[[2]] <- my_RF_model
aa_Malat1_mlist_chunk_featImp[[2]] <- importance(my_RF_model)
names(aa_Malat1_mlist_chunk_featImp[[2]]) <- my_name_dic[match(names(aa_Malat1_mlist_chunk_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_MALAT1/Partition_6_modified_dataset_Malat1___chunk___first_229_last_1492___RFmodel.RData")
aa_Malat1_mlist_chunk[[3]] <- my_RF_model
aa_Malat1_mlist_chunk_featImp[[3]] <- importance(my_RF_model)
names(aa_Malat1_mlist_chunk_featImp[[3]]) <- my_name_dic[match(names(aa_Malat1_mlist_chunk_featImp[[3]]), my_name_dic[, 2]),1]

names(aa_Malat1_mlist_chunk) <- c("All", "seq_only", "kmer_only")
aa_pr_multimodel_P6_Malat1_chunk <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Malat1_mlist_chunk,
                                                                my_dataset = aa_test_data,
                                                                my_label = aa_test_data[, ncol(aa_test_data)],
                                                                file_name = paste0(aaplot_strore_address,"Malat1_All_vs_seq_vs_kmer","__chunk___test_all_P6.png"),
                                                                train_perf = F,
                                                                neg_value = "neg",
                                                                pos_value = "pos")


png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison/Malat1_feature_rand.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Malat1_mlist_rand_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Malat1_mlist_rand_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Malat1_mlist_rand_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmerOnly")
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison/Malat1_feature_chunk.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Malat1_mlist_chunk_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Malat1_mlist_chunk_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Malat1_mlist_chunk_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmerOnly")
dev.off()

######################################################
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Neat1.RData")
my_partition_df <- my_partition_rand
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Neat1_mlist_nameDic <- my_name_dic
aa_Neat1_mlist_rand_featImp <- list()
aa_Neat1_mlist_rand <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_tripCorr/Partition_6_modified_dataset_Neat1___rand___full___RFmodel.RData")
aa_Neat1_mlist_rand[[1]] <- my_RF_model
aa_Neat1_mlist_rand_featImp[[1]] <- importance(my_RF_model)
names(aa_Neat1_mlist_rand_featImp[[1]]) <- my_name_dic[match(names(aa_Neat1_mlist_rand_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation/Partition_6_modified_dataset_Neat1___rand___first_204_last_1464___RFmodel.RData")
aa_Neat1_mlist_rand[[2]] <- my_RF_model
aa_Neat1_mlist_rand_featImp[[2]] <- importance(my_RF_model)
names(aa_Neat1_mlist_rand_featImp[[2]]) <- my_name_dic[match(names(aa_Neat1_mlist_rand_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation/Partition_6_modified_dataset_Neat1___rand___first_204_last_1468___RFmodel.RData")
aa_Neat1_mlist_rand[[3]] <- my_RF_model
aa_Neat1_mlist_rand_featImp[[3]] <- importance(my_RF_model)
names(aa_Neat1_mlist_rand_featImp[[3]]) <- my_name_dic[match(names(aa_Neat1_mlist_rand_featImp[[3]]), my_name_dic[, 2]),1]

names(aa_Neat1_mlist_rand) <- c("All", "seq_only", "kmer_only")

aa_pr_multimodel_P6_Neat1_rand <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Neat1_mlist_rand,
                                                               my_dataset = aa_test_data,
                                                               my_label = aa_test_data[, ncol(aa_test_data)],
                                                               file_name = paste0(aaplot_strore_address,"Neat1_All_vs_seq_vs_kmer","__rand___test_all_P6.png"),
                                                               train_perf = F,
                                                               neg_value = "neg",
                                                               pos_value = "pos")

my_partition_df <- my_partition_chunk
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Neat1_mlist_chunk <- list()
aa_Neat1_mlist_chunk_featImp <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_tripCorr/Partition_6_modified_dataset_Neat1___chunk___full___RFmodel.RData")
aa_Neat1_mlist_chunk[[1]] <- my_RF_model
aa_Neat1_mlist_chunk_featImp[[1]] <- importance(my_RF_model)
names(aa_Neat1_mlist_chunk_featImp[[1]]) <- my_name_dic[match(names(aa_Neat1_mlist_chunk_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation/Partition_6_modified_dataset_Neat1___chunk___first_204_last_1464___RFmodel.RData")
aa_Neat1_mlist_chunk[[2]] <- my_RF_model
aa_Neat1_mlist_chunk_featImp[[2]] <- importance(my_RF_model)
names(aa_Neat1_mlist_chunk_featImp[[2]]) <- my_name_dic[match(names(aa_Neat1_mlist_chunk_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation/Partition_6_modified_dataset_Neat1___chunk___first_204_last_1468___RFmodel.RData")
aa_Neat1_mlist_chunk[[3]] <- my_RF_model
aa_Neat1_mlist_chunk_featImp[[3]] <- importance(my_RF_model)
names(aa_Neat1_mlist_chunk_featImp[[3]]) <- my_name_dic[match(names(aa_Neat1_mlist_chunk_featImp[[3]]), my_name_dic[, 2]),1]

names(aa_Neat1_mlist_chunk) <- c("All", "seq_only", "kmer_only")
aa_pr_multimodel_P6_Neat1_chunk <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Neat1_mlist_chunk,
                                                                my_dataset = aa_test_data,
                                                                my_label = aa_test_data[, ncol(aa_test_data)],
                                                                file_name = paste0(aaplot_strore_address,"Neat1_All_vs_seq_vs_kmer","__chunk___test_all_P6.png"),
                                                                train_perf = F,
                                                                neg_value = "neg",
                                                                pos_value = "pos")


png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison/NEAT1_feature_rand.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Neat1_mlist_rand_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Neat1_mlist_rand_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Neat1_mlist_rand_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmerOnly")
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison/NEAT1_feature_chunk.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Neat1_mlist_chunk_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Neat1_mlist_chunk_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Neat1_mlist_chunk_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmerOnly")
dev.off()

######################################################

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Gm14820.RData")
my_partition_df <- my_partition_rand
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Gm14820_mlist_nameDic <- my_name_dic
aa_Gm14820_mlist_rand_featImp <- list()
aa_Gm14820_mlist_rand <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_tripCorr/Partition_6_modified_dataset_Gm14820___rand___full___RFmodel.RData")
aa_Gm14820_mlist_rand[[1]] <- my_RF_model
aa_Gm14820_mlist_rand_featImp[[1]] <- importance(my_RF_model)
names(aa_Gm14820_mlist_rand_featImp[[1]]) <- my_name_dic[match(names(aa_Gm14820_mlist_rand_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation/Partition_6_modified_dataset_Gm14820___rand___first_129_last_1385___RFmodel.RData")
aa_Gm14820_mlist_rand[[2]] <- my_RF_model
aa_Gm14820_mlist_rand_featImp[[2]] <- importance(my_RF_model)
names(aa_Gm14820_mlist_rand_featImp[[2]]) <- my_name_dic[match(names(aa_Gm14820_mlist_rand_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation/Partition_6_modified_dataset_Gm14820___rand___first_129_last_1389___RFmodel.RData")
aa_Gm14820_mlist_rand[[3]] <- my_RF_model
aa_Gm14820_mlist_rand_featImp[[3]] <- importance(my_RF_model)
names(aa_Gm14820_mlist_rand_featImp[[3]]) <- my_name_dic[match(names(aa_Gm14820_mlist_rand_featImp[[3]]), my_name_dic[, 2]),1]

names(aa_Gm14820_mlist_rand) <- c("All", "seq_only", "kmer_only")

aa_pr_multimodel_P6_Gm14820_rand <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Gm14820_mlist_rand,
                                                               my_dataset = aa_test_data,
                                                               my_label = aa_test_data[, ncol(aa_test_data)],
                                                               file_name = paste0(aaplot_strore_address,"Gm14820_All_vs_seq_vs_kmer","__rand___test_all_P6.png"),
                                                               train_perf = F,
                                                               neg_value = "neg",
                                                               pos_value = "pos")

my_partition_df <- my_partition_chunk
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Gm14820_mlist_chunk <- list()
aa_Gm14820_mlist_chunk_featImp <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_tripCorr/Partition_6_modified_dataset_Gm14820___chunk___full___RFmodel.RData")
aa_Gm14820_mlist_chunk[[1]] <- my_RF_model
aa_Gm14820_mlist_chunk_featImp[[1]] <- importance(my_RF_model)
names(aa_Gm14820_mlist_chunk_featImp[[1]]) <- my_name_dic[match(names(aa_Gm14820_mlist_chunk_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation/Partition_6_modified_dataset_Gm14820___chunk___first_129_last_1385___RFmodel.RData")
aa_Gm14820_mlist_chunk[[2]] <- my_RF_model
aa_Gm14820_mlist_chunk_featImp[[2]] <- importance(my_RF_model)
names(aa_Gm14820_mlist_chunk_featImp[[2]]) <- my_name_dic[match(names(aa_Gm14820_mlist_chunk_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation/Partition_6_modified_dataset_Gm14820___chunk___first_129_last_1389___RFmodel.RData")
aa_Gm14820_mlist_chunk[[3]] <- my_RF_model
aa_Gm14820_mlist_chunk_featImp[[3]] <- importance(my_RF_model)
names(aa_Gm14820_mlist_chunk_featImp[[3]]) <- my_name_dic[match(names(aa_Gm14820_mlist_chunk_featImp[[3]]), my_name_dic[, 2]),1]

names(aa_Gm14820_mlist_chunk) <- c("All", "seq_only", "kmer_only")
aa_pr_multimodel_P6_Gm14820_chunk <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Gm14820_mlist_chunk,
                                                                my_dataset = aa_test_data,
                                                                my_label = aa_test_data[, ncol(aa_test_data)],
                                                                file_name = paste0(aaplot_strore_address,"Gm14820_All_vs_seq_vs_kmer","__chunk___test_all_P6.png"),
                                                                train_perf = F,
                                                                neg_value = "neg",
                                                                pos_value = "pos")


png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison/Gm14820_feature_rand.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Gm14820_mlist_rand_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Gm14820_mlist_rand_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Gm14820_mlist_rand_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmerOnly")
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison/Gm14820_feature_chunk.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Gm14820_mlist_chunk_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Gm14820_mlist_chunk_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Gm14820_mlist_chunk_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmerOnly")
dev.off()

######################################################

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Kcnq1ot1.RData")
my_partition_df <- my_partition_rand
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Kcnq1ot1_mlist_nameDic <- my_name_dic
aa_Kcnq1ot1_mlist_rand_featImp <- list()
aa_Kcnq1ot1_mlist_rand <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_tripCorr/Partition_6_modified_dataset_Kcnq1ot1___rand___full___RFmodel.RData")
aa_Kcnq1ot1_mlist_rand[[1]] <- my_RF_model
aa_Kcnq1ot1_mlist_rand_featImp[[1]] <- importance(my_RF_model)
names(aa_Kcnq1ot1_mlist_rand_featImp[[1]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_rand_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation/Partition_6_modified_dataset_Kcnq1ot1___rand___first_200_last_1453___RFmodel.RData")
aa_Kcnq1ot1_mlist_rand[[2]] <- my_RF_model
aa_Kcnq1ot1_mlist_rand_featImp[[2]] <- importance(my_RF_model)
names(aa_Kcnq1ot1_mlist_rand_featImp[[2]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_rand_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation/Partition_6_modified_dataset_Kcnq1ot1___rand___first_200_last_1457___RFmodel.RData")
aa_Kcnq1ot1_mlist_rand[[3]] <- my_RF_model
aa_Kcnq1ot1_mlist_rand_featImp[[3]] <- importance(my_RF_model)
names(aa_Kcnq1ot1_mlist_rand_featImp[[3]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_rand_featImp[[3]]), my_name_dic[, 2]),1]

names(aa_Kcnq1ot1_mlist_rand) <- c("All", "seq_only", "kmer_only")

aa_pr_multimodel_P6_Kcnq1ot1_rand <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Kcnq1ot1_mlist_rand,
                                                                 my_dataset = aa_test_data,
                                                                 my_label = aa_test_data[, ncol(aa_test_data)],
                                                                 file_name = paste0(aaplot_strore_address,"Kcnq1ot1_All_vs_seq_vs_kmer","__rand___test_all_P6.png"),
                                                                 train_perf = F,
                                                                 neg_value = "neg",
                                                                 pos_value = "pos")

my_partition_df <- my_partition_chunk
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Kcnq1ot1_mlist_chunk <- list()
aa_Kcnq1ot1_mlist_chunk_featImp <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_tripCorr/Partition_6_modified_dataset_Kcnq1ot1___chunk___full___RFmodel.RData")
aa_Kcnq1ot1_mlist_chunk[[1]] <- my_RF_model
aa_Kcnq1ot1_mlist_chunk_featImp[[1]] <- importance(my_RF_model)
names(aa_Kcnq1ot1_mlist_chunk_featImp[[1]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_chunk_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation/Partition_6_modified_dataset_Kcnq1ot1___chunk___first_200_last_1453___RFmodel.RData")
aa_Kcnq1ot1_mlist_chunk[[2]] <- my_RF_model
aa_Kcnq1ot1_mlist_chunk_featImp[[2]] <- importance(my_RF_model)
names(aa_Kcnq1ot1_mlist_chunk_featImp[[2]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_chunk_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation/Partition_6_modified_dataset_Kcnq1ot1___chunk___first_200_last_1457___RFmodel.RData")
aa_Kcnq1ot1_mlist_chunk[[3]] <- my_RF_model
aa_Kcnq1ot1_mlist_chunk_featImp[[3]] <- importance(my_RF_model)
names(aa_Kcnq1ot1_mlist_chunk_featImp[[3]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_chunk_featImp[[3]]), my_name_dic[, 2]),1]

names(aa_Kcnq1ot1_mlist_chunk) <- c("All", "seq_only", "kmer_only")
aa_pr_multimodel_P6_Kcnq1ot1_chunk <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Kcnq1ot1_mlist_chunk,
                                                                  my_dataset = aa_test_data,
                                                                  my_label = aa_test_data[, ncol(aa_test_data)],
                                                                  file_name = paste0(aaplot_strore_address,"Kcnq1ot1_All_vs_seq_vs_kmer","__chunk___test_all_P6.png"),
                                                                  train_perf = F,
                                                                  neg_value = "neg",
                                                                  pos_value = "pos")


png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison/Kcnq1ot1_feature_rand.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Kcnq1ot1_mlist_rand_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Kcnq1ot1_mlist_rand_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Kcnq1ot1_mlist_rand_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmerOnly")
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison/Kcnq1ot1_feature_chunk.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Kcnq1ot1_mlist_chunk_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Kcnq1ot1_mlist_chunk_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Kcnq1ot1_mlist_chunk_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmerOnly")
dev.off()


######################################################

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Trerf1.RData")
my_partition_df <- my_partition_rand
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Trerf1_mlist_nameDic <- my_name_dic
aa_Trerf1_mlist_rand_featImp <- list()
aa_Trerf1_mlist_rand <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_tripCorr/Partition_6_modified_dataset_Trerf1___rand___full___RFmodel.RData")
aa_Trerf1_mlist_rand[[1]] <- my_RF_model
aa_Trerf1_mlist_rand_featImp[[1]] <- importance(my_RF_model)
names(aa_Trerf1_mlist_rand_featImp[[1]]) <- my_name_dic[match(names(aa_Trerf1_mlist_rand_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation/Partition_6_modified_dataset_Trerf1___rand___first_173_last_1425___RFmodel.RData")
aa_Trerf1_mlist_rand[[2]] <- my_RF_model
aa_Trerf1_mlist_rand_featImp[[2]] <- importance(my_RF_model)
names(aa_Trerf1_mlist_rand_featImp[[2]]) <- my_name_dic[match(names(aa_Trerf1_mlist_rand_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation/Partition_6_modified_dataset_Trerf1___rand___first_173_last_1429___RFmodel.RData")
aa_Trerf1_mlist_rand[[3]] <- my_RF_model
aa_Trerf1_mlist_rand_featImp[[3]] <- importance(my_RF_model)
names(aa_Trerf1_mlist_rand_featImp[[3]]) <- my_name_dic[match(names(aa_Trerf1_mlist_rand_featImp[[3]]), my_name_dic[, 2]),1]

names(aa_Trerf1_mlist_rand) <- c("All", "seq_only", "kmer_only")

aa_pr_multimodel_P6_Trerf1_rand <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Trerf1_mlist_rand,
                                                                  my_dataset = aa_test_data,
                                                                  my_label = aa_test_data[, ncol(aa_test_data)],
                                                                  file_name = paste0(aaplot_strore_address,"Trerf1_All_vs_seq_vs_kmer","__rand___test_all_P6.png"),
                                                                  train_perf = F,
                                                                  neg_value = "neg",
                                                                  pos_value = "pos")

my_partition_df <- my_partition_chunk
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Trerf1_mlist_chunk <- list()
aa_Trerf1_mlist_chunk_featImp <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_tripCorr/Partition_6_modified_dataset_Trerf1___chunk___full___RFmodel.RData")
aa_Trerf1_mlist_chunk[[1]] <- my_RF_model
aa_Trerf1_mlist_chunk_featImp[[1]] <- importance(my_RF_model)
names(aa_Trerf1_mlist_chunk_featImp[[1]]) <- my_name_dic[match(names(aa_Trerf1_mlist_chunk_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation/Partition_6_modified_dataset_Trerf1___chunk___first_173_last_1425___RFmodel.RData")
aa_Trerf1_mlist_chunk[[2]] <- my_RF_model
aa_Trerf1_mlist_chunk_featImp[[2]] <- importance(my_RF_model)
names(aa_Trerf1_mlist_chunk_featImp[[2]]) <- my_name_dic[match(names(aa_Trerf1_mlist_chunk_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation/Partition_6_modified_dataset_Trerf1___chunk___first_173_last_1429___RFmodel.RData")
aa_Trerf1_mlist_chunk[[3]] <- my_RF_model
aa_Trerf1_mlist_chunk_featImp[[3]] <- importance(my_RF_model)
names(aa_Trerf1_mlist_chunk_featImp[[3]]) <- my_name_dic[match(names(aa_Trerf1_mlist_chunk_featImp[[3]]), my_name_dic[, 2]),1]

names(aa_Trerf1_mlist_chunk) <- c("All", "seq_only", "kmer_only")
aa_pr_multimodel_P6_Trerf1_chunk <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Trerf1_mlist_chunk,
                                                                   my_dataset = aa_test_data,
                                                                   my_label = aa_test_data[, ncol(aa_test_data)],
                                                                   file_name = paste0(aaplot_strore_address,"Trerf1_All_vs_seq_vs_kmer","__chunk___test_all_P6.png"),
                                                                   train_perf = F,
                                                                   neg_value = "neg",
                                                                   pos_value = "pos")


png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison/Trerf1_feature_rand.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Trerf1_mlist_rand_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Trerf1_mlist_rand_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Trerf1_mlist_rand_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmerOnly")
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison/Trerf1_feature_chunk.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Trerf1_mlist_chunk_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Trerf1_mlist_chunk_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Trerf1_mlist_chunk_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmerOnly")
dev.off()

####################################
# Repeating the above excercise after adding the U1snRNA binding site counts

####################################################################################################################################################

# plotting Malat1, NEAT1, Gm14820, Trerf1, Kcnq1ot1 Performance after removing each feature class

aaplot_strore_address <- "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_U1//"
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_U1snRNA.RData")

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Malat1.RData")
my_Dataset$U1snRNA_site <-  Partition_6_U1snRNA[match(rownames(my_Dataset), names(Partition_6_U1snRNA))]
my_Dataset <- my_Dataset[, c(c(1:(ncol(my_Dataset)-2)), (ncol(my_Dataset)), (ncol(my_Dataset) - 1))]

my_partition_df <- my_partition_rand
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

my_name_dic <- rbind(my_name_dic, c("U1snRNA_site", "U1snRNA_site"))
my_name_dic <- my_name_dic[c(c(1:(nrow(my_name_dic) - 2)), nrow(my_name_dic), (nrow(my_name_dic) - 1)), ]

aa_Malat1_mlist_nameDic <- my_name_dic

aa_Malat1_mlist_rand_featImp <- list()
aa_Malat1_mlist_rand <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_U1/Partition_6_modified_dataset_Malat1___rand___full___RFmodel.RData")
aa_Malat1_mlist_rand[[1]] <- my_RF_model
aa_Malat1_mlist_rand_featImp[[1]] <- importance(my_RF_model)
names(aa_Malat1_mlist_rand_featImp[[1]]) <- my_name_dic[match(names(aa_Malat1_mlist_rand_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation_U1/Partition_6_modified_dataset_Malat1___rand___first_229_last_1488___RFmodel.RData")
aa_Malat1_mlist_rand[[2]] <- my_RF_model
aa_Malat1_mlist_rand_featImp[[2]] <- importance(my_RF_model)
names(aa_Malat1_mlist_rand_featImp[[2]]) <- my_name_dic[match(names(aa_Malat1_mlist_rand_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation_U1/Partition_6_modified_dataset_Malat1___rand___first_229_last_1492___RFmodel.RData")
aa_Malat1_mlist_rand[[3]] <- my_RF_model
aa_Malat1_mlist_rand_featImp[[3]] <- importance(my_RF_model)
names(aa_Malat1_mlist_rand_featImp[[3]]) <- my_name_dic[match(names(aa_Malat1_mlist_rand_featImp[[3]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_OnlyU1/Partition_6_modified_dataset_Malat1___rand___onlyU1___RFmodel.RData")
aa_Malat1_mlist_rand[[4]] <- my_RF_model
aa_Malat1_mlist_rand_featImp[[4]] <- importance(my_RF_model)
names(aa_Malat1_mlist_rand_featImp[[4]]) <- my_name_dic[match(names(aa_Malat1_mlist_rand_featImp[[4]]), my_name_dic[, 2]),1]


names(aa_Malat1_mlist_rand) <- c("All", "seq_only", "kmer+U1", "U1_only")

aa_pr_multimodel_P6_Malat1_rand <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Malat1_mlist_rand,
                                                                my_dataset = aa_test_data,
                                                                my_label = aa_test_data[, ncol(aa_test_data)],
                                                                file_name = paste0(aaplot_strore_address,"Malat1_All_vs_seq_vs_kmer","__rand___test_all_P6.png"),
                                                                train_perf = F,
                                                                neg_value = "neg",
                                                                pos_value = "pos")

my_partition_df <- my_partition_chunk
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Malat1_mlist_chunk <- list()
aa_Malat1_mlist_chunk_featImp <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_U1/Partition_6_modified_dataset_Malat1___chunk___full___RFmodel.RData")
aa_Malat1_mlist_chunk[[1]] <- my_RF_model
aa_Malat1_mlist_chunk_featImp[[1]] <- importance(my_RF_model)
names(aa_Malat1_mlist_chunk_featImp[[1]]) <- my_name_dic[match(names(aa_Malat1_mlist_chunk_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation_U1/Partition_6_modified_dataset_Malat1___chunk___first_229_last_1488___RFmodel.RData")
aa_Malat1_mlist_chunk[[2]] <- my_RF_model
aa_Malat1_mlist_chunk_featImp[[2]] <- importance(my_RF_model)
names(aa_Malat1_mlist_chunk_featImp[[2]]) <- my_name_dic[match(names(aa_Malat1_mlist_chunk_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation_U1/Partition_6_modified_dataset_Malat1___chunk___first_229_last_1492___RFmodel.RData")
aa_Malat1_mlist_chunk[[3]] <- my_RF_model
aa_Malat1_mlist_chunk_featImp[[3]] <- importance(my_RF_model)
names(aa_Malat1_mlist_chunk_featImp[[3]]) <- my_name_dic[match(names(aa_Malat1_mlist_chunk_featImp[[3]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_OnlyU1/Partition_6_modified_dataset_Malat1___chunk___onlyU1___RFmodel.RData")
aa_Malat1_mlist_chunk[[4]] <- my_RF_model
aa_Malat1_mlist_chunk_featImp[[4]] <- importance(my_RF_model)
names(aa_Malat1_mlist_chunk_featImp[[4]]) <- my_name_dic[match(names(aa_Malat1_mlist_chunk_featImp[[4]]), my_name_dic[, 2]),1]

names(aa_Malat1_mlist_chunk) <- c("All", "seq_only", "kmer+U1", "U1_only")
aa_pr_multimodel_P6_Malat1_chunk <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Malat1_mlist_chunk,
                                                                 my_dataset = aa_test_data,
                                                                 my_label = aa_test_data[, ncol(aa_test_data)],
                                                                 file_name = paste0(aaplot_strore_address,"Malat1_All_vs_seq_vs_kmer","__chunk___test_all_P6.png"),
                                                                 train_perf = F,
                                                                 neg_value = "neg",
                                                                 pos_value = "pos")


png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_U1/Malat1_feature_rand.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Malat1_mlist_rand_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Malat1_mlist_rand_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Malat1_mlist_rand_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmer+U1")
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_U1/Malat1_feature_chunk.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Malat1_mlist_chunk_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Malat1_mlist_chunk_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Malat1_mlist_chunk_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmer+U1")
dev.off()

######################################################
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Neat1.RData")
my_Dataset$U1snRNA_site <-  Partition_6_U1snRNA[match(rownames(my_Dataset), names(Partition_6_U1snRNA))]
my_Dataset <- my_Dataset[, c(c(1:(ncol(my_Dataset)-2)), (ncol(my_Dataset)), (ncol(my_Dataset) - 1))]
my_partition_df <- my_partition_rand
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

my_name_dic <- rbind(my_name_dic, c("U1snRNA_site", "U1snRNA_site"))
my_name_dic <- my_name_dic[c(c(1:(nrow(my_name_dic) - 2)), nrow(my_name_dic), (nrow(my_name_dic) - 1)), ]


aa_Neat1_mlist_nameDic <- my_name_dic
aa_Neat1_mlist_rand_featImp <- list()
aa_Neat1_mlist_rand <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_U1/Partition_6_modified_dataset_Neat1___rand___full___RFmodel.RData")
aa_Neat1_mlist_rand[[1]] <- my_RF_model
aa_Neat1_mlist_rand_featImp[[1]] <- importance(my_RF_model)
names(aa_Neat1_mlist_rand_featImp[[1]]) <- my_name_dic[match(names(aa_Neat1_mlist_rand_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation_U1/Partition_6_modified_dataset_Neat1___rand___first_204_last_1464___RFmodel.RData")
aa_Neat1_mlist_rand[[2]] <- my_RF_model
aa_Neat1_mlist_rand_featImp[[2]] <- importance(my_RF_model)
names(aa_Neat1_mlist_rand_featImp[[2]]) <- my_name_dic[match(names(aa_Neat1_mlist_rand_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation_U1/Partition_6_modified_dataset_Neat1___rand___first_204_last_1468___RFmodel.RData")
aa_Neat1_mlist_rand[[3]] <- my_RF_model
aa_Neat1_mlist_rand_featImp[[3]] <- importance(my_RF_model)
names(aa_Neat1_mlist_rand_featImp[[3]]) <- my_name_dic[match(names(aa_Neat1_mlist_rand_featImp[[3]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_OnlyU1/Partition_6_modified_dataset_Neat1___rand___onlyU1___RFmodel.RData")
aa_Neat1_mlist_rand[[4]] <- my_RF_model
aa_Neat1_mlist_rand_featImp[[4]] <- importance(my_RF_model)
names(aa_Neat1_mlist_rand_featImp[[4]]) <- my_name_dic[match(names(aa_Neat1_mlist_rand_featImp[[4]]), my_name_dic[, 2]),1]

names(aa_Neat1_mlist_rand) <- c("All", "seq_only", "kmer+U1", "U1_only")

aa_pr_multimodel_P6_Neat1_rand <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Neat1_mlist_rand,
                                                               my_dataset = aa_test_data,
                                                               my_label = aa_test_data[, ncol(aa_test_data)],
                                                               file_name = paste0(aaplot_strore_address,"Neat1_All_vs_seq_vs_kmer","__rand___test_all_P6.png"),
                                                               train_perf = F,
                                                               neg_value = "neg",
                                                               pos_value = "pos")

my_partition_df <- my_partition_chunk
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Neat1_mlist_chunk <- list()
aa_Neat1_mlist_chunk_featImp <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_U1/Partition_6_modified_dataset_Neat1___chunk___full___RFmodel.RData")
aa_Neat1_mlist_chunk[[1]] <- my_RF_model
aa_Neat1_mlist_chunk_featImp[[1]] <- importance(my_RF_model)
names(aa_Neat1_mlist_chunk_featImp[[1]]) <- my_name_dic[match(names(aa_Neat1_mlist_chunk_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation_U1/Partition_6_modified_dataset_Neat1___chunk___first_204_last_1464___RFmodel.RData")
aa_Neat1_mlist_chunk[[2]] <- my_RF_model
aa_Neat1_mlist_chunk_featImp[[2]] <- importance(my_RF_model)
names(aa_Neat1_mlist_chunk_featImp[[2]]) <- my_name_dic[match(names(aa_Neat1_mlist_chunk_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation_U1/Partition_6_modified_dataset_Neat1___chunk___first_204_last_1468___RFmodel.RData")
aa_Neat1_mlist_chunk[[3]] <- my_RF_model
aa_Neat1_mlist_chunk_featImp[[3]] <- importance(my_RF_model)
names(aa_Neat1_mlist_chunk_featImp[[3]]) <- my_name_dic[match(names(aa_Neat1_mlist_chunk_featImp[[3]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_OnlyU1/Partition_6_modified_dataset_Neat1___chunk___onlyU1___RFmodel.RData")
aa_Neat1_mlist_chunk[[4]] <- my_RF_model
aa_Neat1_mlist_chunk_featImp[[4]] <- importance(my_RF_model)
names(aa_Neat1_mlist_chunk_featImp[[4]]) <- my_name_dic[match(names(aa_Neat1_mlist_chunk_featImp[[4]]), my_name_dic[, 2]),1]

names(aa_Neat1_mlist_chunk) <- c("All", "seq_only", "kmer+U1", "U1_only")
aa_pr_multimodel_P6_Neat1_chunk <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Neat1_mlist_chunk,
                                                                my_dataset = aa_test_data,
                                                                my_label = aa_test_data[, ncol(aa_test_data)],
                                                                file_name = paste0(aaplot_strore_address,"Neat1_All_vs_seq_vs_kmer","__chunk___test_all_P6.png"),
                                                                train_perf = F,
                                                                neg_value = "neg",
                                                                pos_value = "pos")


png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_U1/NEAT1_feature_rand.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Neat1_mlist_rand_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Neat1_mlist_rand_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Neat1_mlist_rand_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmer+U1")
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_U1/NEAT1_feature_chunk.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Neat1_mlist_chunk_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Neat1_mlist_chunk_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Neat1_mlist_chunk_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmer+U1")
dev.off()

######################################################

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Gm14820.RData")
my_Dataset$U1snRNA_site <-  Partition_6_U1snRNA[match(rownames(my_Dataset), names(Partition_6_U1snRNA))]
my_Dataset <- my_Dataset[, c(c(1:(ncol(my_Dataset)-2)), (ncol(my_Dataset)), (ncol(my_Dataset) - 1))]
my_partition_df <- my_partition_rand
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

my_name_dic <- rbind(my_name_dic, c("U1snRNA_site", "U1snRNA_site"))
my_name_dic <- my_name_dic[c(c(1:(nrow(my_name_dic) - 2)), nrow(my_name_dic), (nrow(my_name_dic) - 1)), ]

aa_Gm14820_mlist_nameDic <- my_name_dic
aa_Gm14820_mlist_rand_featImp <- list()
aa_Gm14820_mlist_rand <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_U1/Partition_6_modified_dataset_Gm14820___rand___full___RFmodel.RData")
aa_Gm14820_mlist_rand[[1]] <- my_RF_model
aa_Gm14820_mlist_rand_featImp[[1]] <- importance(my_RF_model)
names(aa_Gm14820_mlist_rand_featImp[[1]]) <- my_name_dic[match(names(aa_Gm14820_mlist_rand_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation_U1/Partition_6_modified_dataset_Gm14820___rand___first_129_last_1385___RFmodel.RData")
aa_Gm14820_mlist_rand[[2]] <- my_RF_model
aa_Gm14820_mlist_rand_featImp[[2]] <- importance(my_RF_model)
names(aa_Gm14820_mlist_rand_featImp[[2]]) <- my_name_dic[match(names(aa_Gm14820_mlist_rand_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation_U1/Partition_6_modified_dataset_Gm14820___rand___first_129_last_1389___RFmodel.RData")
aa_Gm14820_mlist_rand[[3]] <- my_RF_model
aa_Gm14820_mlist_rand_featImp[[3]] <- importance(my_RF_model)
names(aa_Gm14820_mlist_rand_featImp[[3]]) <- my_name_dic[match(names(aa_Gm14820_mlist_rand_featImp[[3]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_OnlyU1/Partition_6_modified_dataset_Gm14820___rand___onlyU1___RFmodel.RData")
aa_Gm14820_mlist_rand[[4]] <- my_RF_model
aa_Gm14820_mlist_rand_featImp[[4]] <- importance(my_RF_model)
names(aa_Gm14820_mlist_rand_featImp[[4]]) <- my_name_dic[match(names(aa_Gm14820_mlist_rand_featImp[[4]]), my_name_dic[, 2]),1]

names(aa_Gm14820_mlist_rand) <- c("All", "seq_only", "kmer+U1", "U1_only")

aa_pr_multimodel_P6_Gm14820_rand <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Gm14820_mlist_rand,
                                                                 my_dataset = aa_test_data,
                                                                 my_label = aa_test_data[, ncol(aa_test_data)],
                                                                 file_name = paste0(aaplot_strore_address,"Gm14820_All_vs_seq_vs_kmer","__rand___test_all_P6.png"),
                                                                 train_perf = F,
                                                                 neg_value = "neg",
                                                                 pos_value = "pos")

my_partition_df <- my_partition_chunk
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Gm14820_mlist_chunk <- list()
aa_Gm14820_mlist_chunk_featImp <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_U1/Partition_6_modified_dataset_Gm14820___chunk___full___RFmodel.RData")
aa_Gm14820_mlist_chunk[[1]] <- my_RF_model
aa_Gm14820_mlist_chunk_featImp[[1]] <- importance(my_RF_model)
names(aa_Gm14820_mlist_chunk_featImp[[1]]) <- my_name_dic[match(names(aa_Gm14820_mlist_chunk_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation_U1/Partition_6_modified_dataset_Gm14820___chunk___first_129_last_1385___RFmodel.RData")
aa_Gm14820_mlist_chunk[[2]] <- my_RF_model
aa_Gm14820_mlist_chunk_featImp[[2]] <- importance(my_RF_model)
names(aa_Gm14820_mlist_chunk_featImp[[2]]) <- my_name_dic[match(names(aa_Gm14820_mlist_chunk_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation_U1/Partition_6_modified_dataset_Gm14820___chunk___first_129_last_1389___RFmodel.RData")
aa_Gm14820_mlist_chunk[[3]] <- my_RF_model
aa_Gm14820_mlist_chunk_featImp[[3]] <- importance(my_RF_model)
names(aa_Gm14820_mlist_chunk_featImp[[3]]) <- my_name_dic[match(names(aa_Gm14820_mlist_chunk_featImp[[3]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_OnlyU1/Partition_6_modified_dataset_Gm14820___chunk___onlyU1___RFmodel.RData")
aa_Gm14820_mlist_chunk[[4]] <- my_RF_model
aa_Gm14820_mlist_chunk_featImp[[4]] <- importance(my_RF_model)
names(aa_Gm14820_mlist_chunk_featImp[[4]]) <- my_name_dic[match(names(aa_Gm14820_mlist_chunk_featImp[[4]]), my_name_dic[, 2]),1]


names(aa_Gm14820_mlist_chunk) <- c("All", "seq_only", "kmer+U1", "U1_only")
aa_pr_multimodel_P6_Gm14820_chunk <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Gm14820_mlist_chunk,
                                                                  my_dataset = aa_test_data,
                                                                  my_label = aa_test_data[, ncol(aa_test_data)],
                                                                  file_name = paste0(aaplot_strore_address,"Gm14820_All_vs_seq_vs_kmer","__chunk___test_all_P6.png"),
                                                                  train_perf = F,
                                                                  neg_value = "neg",
                                                                  pos_value = "pos")


png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_U1/Gm14820_feature_rand.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Gm14820_mlist_rand_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Gm14820_mlist_rand_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Gm14820_mlist_rand_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmer+U1")
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_U1/Gm14820_feature_chunk.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Gm14820_mlist_chunk_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Gm14820_mlist_chunk_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Gm14820_mlist_chunk_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmer+U1")
dev.off()

######################################################

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Kcnq1ot1.RData")
my_Dataset$U1snRNA_site <-  Partition_6_U1snRNA[match(rownames(my_Dataset), names(Partition_6_U1snRNA))]
my_Dataset <- my_Dataset[, c(c(1:(ncol(my_Dataset)-2)), (ncol(my_Dataset)), (ncol(my_Dataset) - 1))]

my_partition_df <- my_partition_rand
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

my_name_dic <- rbind(my_name_dic, c("U1snRNA_site", "U1snRNA_site"))
my_name_dic <- my_name_dic[c(c(1:(nrow(my_name_dic) - 2)), nrow(my_name_dic), (nrow(my_name_dic) - 1)), ]

aa_Kcnq1ot1_mlist_nameDic <- my_name_dic
aa_Kcnq1ot1_mlist_rand_featImp <- list()
aa_Kcnq1ot1_mlist_rand <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_U1/Partition_6_modified_dataset_Kcnq1ot1___rand___full___RFmodel.RData")
aa_Kcnq1ot1_mlist_rand[[1]] <- my_RF_model
aa_Kcnq1ot1_mlist_rand_featImp[[1]] <- importance(my_RF_model)
names(aa_Kcnq1ot1_mlist_rand_featImp[[1]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_rand_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation_U1/Partition_6_modified_dataset_Kcnq1ot1___rand___first_200_last_1453___RFmodel.RData")
aa_Kcnq1ot1_mlist_rand[[2]] <- my_RF_model
aa_Kcnq1ot1_mlist_rand_featImp[[2]] <- importance(my_RF_model)
names(aa_Kcnq1ot1_mlist_rand_featImp[[2]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_rand_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation_U1/Partition_6_modified_dataset_Kcnq1ot1___rand___first_200_last_1457___RFmodel.RData")
aa_Kcnq1ot1_mlist_rand[[3]] <- my_RF_model
aa_Kcnq1ot1_mlist_rand_featImp[[3]] <- importance(my_RF_model)
names(aa_Kcnq1ot1_mlist_rand_featImp[[3]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_rand_featImp[[3]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_OnlyU1/Partition_6_modified_dataset_Kcnq1ot1___rand___onlyU1___RFmodel.RData")
aa_Kcnq1ot1_mlist_rand[[4]] <- my_RF_model
aa_Kcnq1ot1_mlist_rand_featImp[[4]] <- importance(my_RF_model)
names(aa_Kcnq1ot1_mlist_rand_featImp[[4]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_rand_featImp[[4]]), my_name_dic[, 2]),1]

names(aa_Kcnq1ot1_mlist_rand) <- c("All", "seq_only", "kmer+U1", "U1_only")

aa_pr_multimodel_P6_Kcnq1ot1_rand <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Kcnq1ot1_mlist_rand,
                                                                  my_dataset = aa_test_data,
                                                                  my_label = aa_test_data[, ncol(aa_test_data)],
                                                                  file_name = paste0(aaplot_strore_address,"Kcnq1ot1_All_vs_seq_vs_kmer","__rand___test_all_P6.png"),
                                                                  train_perf = F,
                                                                  neg_value = "neg",
                                                                  pos_value = "pos")

my_partition_df <- my_partition_chunk
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Kcnq1ot1_mlist_chunk <- list()
aa_Kcnq1ot1_mlist_chunk_featImp <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_U1/Partition_6_modified_dataset_Kcnq1ot1___chunk___full___RFmodel.RData")
aa_Kcnq1ot1_mlist_chunk[[1]] <- my_RF_model
aa_Kcnq1ot1_mlist_chunk_featImp[[1]] <- importance(my_RF_model)
names(aa_Kcnq1ot1_mlist_chunk_featImp[[1]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_chunk_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation_U1/Partition_6_modified_dataset_Kcnq1ot1___chunk___first_200_last_1453___RFmodel.RData")
aa_Kcnq1ot1_mlist_chunk[[2]] <- my_RF_model
aa_Kcnq1ot1_mlist_chunk_featImp[[2]] <- importance(my_RF_model)
names(aa_Kcnq1ot1_mlist_chunk_featImp[[2]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_chunk_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation_U1/Partition_6_modified_dataset_Kcnq1ot1___chunk___first_200_last_1457___RFmodel.RData")
aa_Kcnq1ot1_mlist_chunk[[3]] <- my_RF_model
aa_Kcnq1ot1_mlist_chunk_featImp[[3]] <- importance(my_RF_model)
names(aa_Kcnq1ot1_mlist_chunk_featImp[[3]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_chunk_featImp[[3]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_OnlyU1/Partition_6_modified_dataset_Kcnq1ot1___chunk___onlyU1___RFmodel.RData")
aa_Kcnq1ot1_mlist_chunk[[4]] <- my_RF_model
aa_Kcnq1ot1_mlist_chunk_featImp[[4]] <- importance(my_RF_model)
names(aa_Kcnq1ot1_mlist_chunk_featImp[[4]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_chunk_featImp[[4]]), my_name_dic[, 2]),1]

names(aa_Kcnq1ot1_mlist_chunk) <- c("All", "seq_only", "kmer+U1", "U1_only")
aa_pr_multimodel_P6_Kcnq1ot1_chunk <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Kcnq1ot1_mlist_chunk,
                                                                   my_dataset = aa_test_data,
                                                                   my_label = aa_test_data[, ncol(aa_test_data)],
                                                                   file_name = paste0(aaplot_strore_address,"Kcnq1ot1_All_vs_seq_vs_kmer","__chunk___test_all_P6.png"),
                                                                   train_perf = F,
                                                                   neg_value = "neg",
                                                                   pos_value = "pos")


png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_U1/Kcnq1ot1_feature_rand.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Kcnq1ot1_mlist_rand_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Kcnq1ot1_mlist_rand_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Kcnq1ot1_mlist_rand_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmer+U1")
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_U1/Kcnq1ot1_feature_chunk.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Kcnq1ot1_mlist_chunk_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Kcnq1ot1_mlist_chunk_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Kcnq1ot1_mlist_chunk_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmer+U1")
dev.off()


######################################################

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Trerf1.RData")
my_Dataset$U1snRNA_site <-  Partition_6_U1snRNA[match(rownames(my_Dataset), names(Partition_6_U1snRNA))]
my_Dataset <- my_Dataset[, c(c(1:(ncol(my_Dataset)-2)), (ncol(my_Dataset)), (ncol(my_Dataset) - 1))]

my_partition_df <- my_partition_rand
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

my_name_dic <- rbind(my_name_dic, c("U1snRNA_site", "U1snRNA_site"))
my_name_dic <- my_name_dic[c(c(1:(nrow(my_name_dic) - 2)), nrow(my_name_dic), (nrow(my_name_dic) - 1)), ]


aa_Trerf1_mlist_nameDic <- my_name_dic
aa_Trerf1_mlist_rand_featImp <- list()
aa_Trerf1_mlist_rand <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_U1/Partition_6_modified_dataset_Trerf1___rand___full___RFmodel.RData")
aa_Trerf1_mlist_rand[[1]] <- my_RF_model
aa_Trerf1_mlist_rand_featImp[[1]] <- importance(my_RF_model)
names(aa_Trerf1_mlist_rand_featImp[[1]]) <- my_name_dic[match(names(aa_Trerf1_mlist_rand_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation_U1/Partition_6_modified_dataset_Trerf1___rand___first_173_last_1425___RFmodel.RData")
aa_Trerf1_mlist_rand[[2]] <- my_RF_model
aa_Trerf1_mlist_rand_featImp[[2]] <- importance(my_RF_model)
names(aa_Trerf1_mlist_rand_featImp[[2]]) <- my_name_dic[match(names(aa_Trerf1_mlist_rand_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation_U1/Partition_6_modified_dataset_Trerf1___rand___first_173_last_1429___RFmodel.RData")
aa_Trerf1_mlist_rand[[3]] <- my_RF_model
aa_Trerf1_mlist_rand_featImp[[3]] <- importance(my_RF_model)
names(aa_Trerf1_mlist_rand_featImp[[3]]) <- my_name_dic[match(names(aa_Trerf1_mlist_rand_featImp[[3]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_OnlyU1/Partition_6_modified_dataset_Trerf1___rand___onlyU1___RFmodel.RData")
aa_Trerf1_mlist_rand[[4]] <- my_RF_model
aa_Trerf1_mlist_rand_featImp[[4]] <- importance(my_RF_model)
names(aa_Trerf1_mlist_rand_featImp[[4]]) <- my_name_dic[match(names(aa_Trerf1_mlist_rand_featImp[[4]]), my_name_dic[, 2]),1]

names(aa_Trerf1_mlist_rand) <- c("All", "seq_only", "kmer+U1", "U1_only")

aa_pr_multimodel_P6_Trerf1_rand <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Trerf1_mlist_rand,
                                                                my_dataset = aa_test_data,
                                                                my_label = aa_test_data[, ncol(aa_test_data)],
                                                                file_name = paste0(aaplot_strore_address,"Trerf1_All_vs_seq_vs_kmer","__rand___test_all_P6.png"),
                                                                train_perf = F,
                                                                neg_value = "neg",
                                                                pos_value = "pos")

my_partition_df <- my_partition_chunk
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Trerf1_mlist_chunk <- list()
aa_Trerf1_mlist_chunk_featImp <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_U1/Partition_6_modified_dataset_Trerf1___chunk___full___RFmodel.RData")
aa_Trerf1_mlist_chunk[[1]] <- my_RF_model
aa_Trerf1_mlist_chunk_featImp[[1]] <- importance(my_RF_model)
names(aa_Trerf1_mlist_chunk_featImp[[1]]) <- my_name_dic[match(names(aa_Trerf1_mlist_chunk_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation_U1/Partition_6_modified_dataset_Trerf1___chunk___first_173_last_1425___RFmodel.RData")
aa_Trerf1_mlist_chunk[[2]] <- my_RF_model
aa_Trerf1_mlist_chunk_featImp[[2]] <- importance(my_RF_model)
names(aa_Trerf1_mlist_chunk_featImp[[2]]) <- my_name_dic[match(names(aa_Trerf1_mlist_chunk_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_P6_feature_ablation_U1/Partition_6_modified_dataset_Trerf1___chunk___first_173_last_1429___RFmodel.RData")
aa_Trerf1_mlist_chunk[[3]] <- my_RF_model
aa_Trerf1_mlist_chunk_featImp[[3]] <- importance(my_RF_model)
names(aa_Trerf1_mlist_chunk_featImp[[3]]) <- my_name_dic[match(names(aa_Trerf1_mlist_chunk_featImp[[3]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_OnlyU1/Partition_6_modified_dataset_Trerf1___chunk___onlyU1___RFmodel.RData")
aa_Trerf1_mlist_chunk[[4]] <- my_RF_model
aa_Trerf1_mlist_chunk_featImp[[4]] <- importance(my_RF_model)
names(aa_Trerf1_mlist_chunk_featImp[[4]]) <- my_name_dic[match(names(aa_Trerf1_mlist_chunk_featImp[[4]]), my_name_dic[, 2]),1]

names(aa_Trerf1_mlist_chunk) <- c("All", "seq_only", "kmer+U1", "U1_only")
aa_pr_multimodel_P6_Trerf1_chunk <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Trerf1_mlist_chunk,
                                                                 my_dataset = aa_test_data,
                                                                 my_label = aa_test_data[, ncol(aa_test_data)],
                                                                 file_name = paste0(aaplot_strore_address,"Trerf1_All_vs_seq_vs_kmer","__chunk___test_all_P6.png"),
                                                                 train_perf = F,
                                                                 neg_value = "neg",
                                                                 pos_value = "pos")


png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_U1/Trerf1_feature_rand.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Trerf1_mlist_rand_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Trerf1_mlist_rand_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Trerf1_mlist_rand_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmer+U1")
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_U1/Trerf1_feature_chunk.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Trerf1_mlist_chunk_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Trerf1_mlist_chunk_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Trerf1_mlist_chunk_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmer+U1")
dev.off()

####################################################################################################################################################

# plotting Malat1, NEAT1, Gm14820, Trerf1, Kcnq1ot1 Performance after removing each feature class
# after filtering the RBP scores based on expression
source("~/Documents/Shayan/BioInf/lncRNA/RBP_filter_func.R")
aaplot_strore_address <- "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_RBPfiltered/"
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_U1snRNA.RData")

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Malat1.RData")
my_Dataset$U1snRNA_site <-  Partition_6_U1snRNA[match(rownames(my_Dataset), names(Partition_6_U1snRNA))]
my_Dataset <- my_Dataset[, c(c(1:(ncol(my_Dataset)-2)), (ncol(my_Dataset)), (ncol(my_Dataset) - 1))]
my_name_dic <- rbind(my_name_dic, c("U1snRNA_site", "U1snRNA_site"))
my_name_dic <- my_name_dic[c(c(1:(nrow(my_name_dic) - 2)), nrow(my_name_dic), (nrow(my_name_dic) - 1)), ]
my_Dataset <- Filter_RBP(inp_dataset = my_Dataset,
                         inp_nameDic = my_name_dic,
                         partition_Exp_features_address = "~/Documents/Shayan/BioInf/lncRNA/partition_6_expression_features.RData")


my_partition_df <- my_partition_rand
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]



aa_Malat1_mlist_nameDic <- my_name_dic

aa_Malat1_mlist_rand_featImp <- list()
aa_Malat1_mlist_rand <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_RBP_Filtered/Partition_6_modified_dataset_Malat1___rand___full___RFmodel.RData")
aa_Malat1_mlist_rand[[1]] <- my_RF_model
aa_Malat1_mlist_rand_featImp[[1]] <- importance(my_RF_model)
names(aa_Malat1_mlist_rand_featImp[[1]]) <- my_name_dic[match(names(aa_Malat1_mlist_rand_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_RBP_Filtered_Ablation/Partition_6_modified_dataset_Malat1___rand___first_229_last_1488___RFmodel.RData")
aa_Malat1_mlist_rand[[2]] <- my_RF_model
aa_Malat1_mlist_rand_featImp[[2]] <- importance(my_RF_model)
names(aa_Malat1_mlist_rand_featImp[[2]]) <- my_name_dic[match(names(aa_Malat1_mlist_rand_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_RBP_Filtered_Ablation/Partition_6_modified_dataset_Malat1___rand___first_229_last_1492___RFmodel.RData")
aa_Malat1_mlist_rand[[3]] <- my_RF_model
aa_Malat1_mlist_rand_featImp[[3]] <- importance(my_RF_model)
names(aa_Malat1_mlist_rand_featImp[[3]]) <- my_name_dic[match(names(aa_Malat1_mlist_rand_featImp[[3]]), my_name_dic[, 2]),1]

# load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_OnlyU1/Partition_6_modified_dataset_Malat1___rand___onlyU1___RFmodel.RData")
# aa_Malat1_mlist_rand[[4]] <- my_RF_model
# aa_Malat1_mlist_rand_featImp[[4]] <- importance(my_RF_model)
# names(aa_Malat1_mlist_rand_featImp[[4]]) <- my_name_dic[match(names(aa_Malat1_mlist_rand_featImp[[4]]), my_name_dic[, 2]),1]


names(aa_Malat1_mlist_rand) <- c("All", "seq_only", "kmer+U1")

aa_pr_multimodel_P6_Malat1_rand <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Malat1_mlist_rand,
                                                                my_dataset = aa_test_data,
                                                                my_label = aa_test_data[, ncol(aa_test_data)],
                                                                file_name = paste0(aaplot_strore_address,"Malat1_All_vs_seq_vs_kmer","__rand___test_all_P6.png"),
                                                                train_perf = F,
                                                                neg_value = "neg",
                                                                pos_value = "pos")

my_partition_df <- my_partition_chunk
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Malat1_mlist_chunk <- list()
aa_Malat1_mlist_chunk_featImp <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_RBP_Filtered/Partition_6_modified_dataset_Malat1___chunk___full___RFmodel.RData")
aa_Malat1_mlist_chunk[[1]] <- my_RF_model
aa_Malat1_mlist_chunk_featImp[[1]] <- importance(my_RF_model)
names(aa_Malat1_mlist_chunk_featImp[[1]]) <- my_name_dic[match(names(aa_Malat1_mlist_chunk_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_RBP_Filtered_Ablation/Partition_6_modified_dataset_Malat1___chunk___first_229_last_1488___RFmodel.RData")
aa_Malat1_mlist_chunk[[2]] <- my_RF_model
aa_Malat1_mlist_chunk_featImp[[2]] <- importance(my_RF_model)
names(aa_Malat1_mlist_chunk_featImp[[2]]) <- my_name_dic[match(names(aa_Malat1_mlist_chunk_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_RBP_Filtered_Ablation/Partition_6_modified_dataset_Malat1___chunk___first_229_last_1492___RFmodel.RData")
aa_Malat1_mlist_chunk[[3]] <- my_RF_model
aa_Malat1_mlist_chunk_featImp[[3]] <- importance(my_RF_model)
names(aa_Malat1_mlist_chunk_featImp[[3]]) <- my_name_dic[match(names(aa_Malat1_mlist_chunk_featImp[[3]]), my_name_dic[, 2]),1]

# load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_OnlyU1/Partition_6_modified_dataset_Malat1___chunk___onlyU1___RFmodel.RData")
# aa_Malat1_mlist_chunk[[4]] <- my_RF_model
# aa_Malat1_mlist_chunk_featImp[[4]] <- importance(my_RF_model)
# names(aa_Malat1_mlist_chunk_featImp[[4]]) <- my_name_dic[match(names(aa_Malat1_mlist_chunk_featImp[[4]]), my_name_dic[, 2]),1]

names(aa_Malat1_mlist_chunk) <- c("All", "seq_only", "kmer+U1")
aa_pr_multimodel_P6_Malat1_chunk <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Malat1_mlist_chunk,
                                                                 my_dataset = aa_test_data,
                                                                 my_label = aa_test_data[, ncol(aa_test_data)],
                                                                 file_name = paste0(aaplot_strore_address,"Malat1_All_vs_seq_vs_kmer","__chunk___test_all_P6.png"),
                                                                 train_perf = F,
                                                                 neg_value = "neg",
                                                                 pos_value = "pos")


png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_RBPfiltered/Malat1_feature_rand.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Malat1_mlist_rand_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Malat1_mlist_rand_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Malat1_mlist_rand_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmer+U1")
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_RBPfiltered/Malat1_feature_chunk.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Malat1_mlist_chunk_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Malat1_mlist_chunk_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Malat1_mlist_chunk_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmer+U1")
dev.off()

######################################################
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Neat1.RData")
my_Dataset$U1snRNA_site <-  Partition_6_U1snRNA[match(rownames(my_Dataset), names(Partition_6_U1snRNA))]
my_Dataset <- my_Dataset[, c(c(1:(ncol(my_Dataset)-2)), (ncol(my_Dataset)), (ncol(my_Dataset) - 1))]
my_name_dic <- rbind(my_name_dic, c("U1snRNA_site", "U1snRNA_site"))
my_name_dic <- my_name_dic[c(c(1:(nrow(my_name_dic) - 2)), nrow(my_name_dic), (nrow(my_name_dic) - 1)), ]
my_Dataset <- Filter_RBP(inp_dataset = my_Dataset,
                         inp_nameDic = my_name_dic,
                         partition_Exp_features_address = "~/Documents/Shayan/BioInf/lncRNA/partition_6_expression_features.RData")


my_partition_df <- my_partition_rand
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Neat1_mlist_nameDic <- my_name_dic
aa_Neat1_mlist_rand_featImp <- list()
aa_Neat1_mlist_rand <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_RBP_Filtered/Partition_6_modified_dataset_Neat1___rand___full___RFmodel.RData")
aa_Neat1_mlist_rand[[1]] <- my_RF_model
aa_Neat1_mlist_rand_featImp[[1]] <- importance(my_RF_model)
names(aa_Neat1_mlist_rand_featImp[[1]]) <- my_name_dic[match(names(aa_Neat1_mlist_rand_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_RBP_Filtered_Ablation/Partition_6_modified_dataset_Neat1___rand___first_204_last_1464___RFmodel.RData")
aa_Neat1_mlist_rand[[2]] <- my_RF_model
aa_Neat1_mlist_rand_featImp[[2]] <- importance(my_RF_model)
names(aa_Neat1_mlist_rand_featImp[[2]]) <- my_name_dic[match(names(aa_Neat1_mlist_rand_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_RBP_Filtered_Ablation/Partition_6_modified_dataset_Neat1___rand___first_204_last_1468___RFmodel.RData")
aa_Neat1_mlist_rand[[3]] <- my_RF_model
aa_Neat1_mlist_rand_featImp[[3]] <- importance(my_RF_model)
names(aa_Neat1_mlist_rand_featImp[[3]]) <- my_name_dic[match(names(aa_Neat1_mlist_rand_featImp[[3]]), my_name_dic[, 2]),1]

# load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_OnlyU1/Partition_6_modified_dataset_Neat1___rand___onlyU1___RFmodel.RData")
# aa_Neat1_mlist_rand[[4]] <- my_RF_model
# aa_Neat1_mlist_rand_featImp[[4]] <- importance(my_RF_model)
# names(aa_Neat1_mlist_rand_featImp[[4]]) <- my_name_dic[match(names(aa_Neat1_mlist_rand_featImp[[4]]), my_name_dic[, 2]),1]

names(aa_Neat1_mlist_rand) <- c("All", "seq_only", "kmer+U1")

aa_pr_multimodel_P6_Neat1_rand <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Neat1_mlist_rand,
                                                               my_dataset = aa_test_data,
                                                               my_label = aa_test_data[, ncol(aa_test_data)],
                                                               file_name = paste0(aaplot_strore_address,"Neat1_All_vs_seq_vs_kmer","__rand___test_all_P6.png"),
                                                               train_perf = F,
                                                               neg_value = "neg",
                                                               pos_value = "pos")

my_partition_df <- my_partition_chunk
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Neat1_mlist_chunk <- list()
aa_Neat1_mlist_chunk_featImp <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_RBP_Filtered/Partition_6_modified_dataset_Neat1___chunk___full___RFmodel.RData")
aa_Neat1_mlist_chunk[[1]] <- my_RF_model
aa_Neat1_mlist_chunk_featImp[[1]] <- importance(my_RF_model)
names(aa_Neat1_mlist_chunk_featImp[[1]]) <- my_name_dic[match(names(aa_Neat1_mlist_chunk_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_RBP_Filtered_Ablation/Partition_6_modified_dataset_Neat1___chunk___first_204_last_1464___RFmodel.RData")
aa_Neat1_mlist_chunk[[2]] <- my_RF_model
aa_Neat1_mlist_chunk_featImp[[2]] <- importance(my_RF_model)
names(aa_Neat1_mlist_chunk_featImp[[2]]) <- my_name_dic[match(names(aa_Neat1_mlist_chunk_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_RBP_Filtered_Ablation/Partition_6_modified_dataset_Neat1___chunk___first_204_last_1468___RFmodel.RData")
aa_Neat1_mlist_chunk[[3]] <- my_RF_model
aa_Neat1_mlist_chunk_featImp[[3]] <- importance(my_RF_model)
names(aa_Neat1_mlist_chunk_featImp[[3]]) <- my_name_dic[match(names(aa_Neat1_mlist_chunk_featImp[[3]]), my_name_dic[, 2]),1]

# load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_OnlyU1/Partition_6_modified_dataset_Neat1___chunk___onlyU1___RFmodel.RData")
# aa_Neat1_mlist_chunk[[4]] <- my_RF_model
# aa_Neat1_mlist_chunk_featImp[[4]] <- importance(my_RF_model)
# names(aa_Neat1_mlist_chunk_featImp[[4]]) <- my_name_dic[match(names(aa_Neat1_mlist_chunk_featImp[[4]]), my_name_dic[, 2]),1]

names(aa_Neat1_mlist_chunk) <- c("All", "seq_only", "kmer+U1")
aa_pr_multimodel_P6_Neat1_chunk <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Neat1_mlist_chunk,
                                                                my_dataset = aa_test_data,
                                                                my_label = aa_test_data[, ncol(aa_test_data)],
                                                                file_name = paste0(aaplot_strore_address,"Neat1_All_vs_seq_vs_kmer","__chunk___test_all_P6.png"),
                                                                train_perf = F,
                                                                neg_value = "neg",
                                                                pos_value = "pos")


png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_RBPfiltered/NEAT1_feature_rand.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Neat1_mlist_rand_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Neat1_mlist_rand_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Neat1_mlist_rand_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmer+U1")
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_RBPfiltered/NEAT1_feature_chunk.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Neat1_mlist_chunk_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Neat1_mlist_chunk_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Neat1_mlist_chunk_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmer+U1")
dev.off()

######################################################

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Gm14820.RData")
my_Dataset$U1snRNA_site <-  Partition_6_U1snRNA[match(rownames(my_Dataset), names(Partition_6_U1snRNA))]
my_Dataset <- my_Dataset[, c(c(1:(ncol(my_Dataset)-2)), (ncol(my_Dataset)), (ncol(my_Dataset) - 1))]
my_name_dic <- rbind(my_name_dic, c("U1snRNA_site", "U1snRNA_site"))
my_name_dic <- my_name_dic[c(c(1:(nrow(my_name_dic) - 2)), nrow(my_name_dic), (nrow(my_name_dic) - 1)), ]
my_Dataset <- Filter_RBP(inp_dataset = my_Dataset,
                         inp_nameDic = my_name_dic,
                         partition_Exp_features_address = "~/Documents/Shayan/BioInf/lncRNA/partition_6_expression_features.RData")


my_partition_df <- my_partition_rand
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Gm14820_mlist_nameDic <- my_name_dic
aa_Gm14820_mlist_rand_featImp <- list()
aa_Gm14820_mlist_rand <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_RBP_Filtered/Partition_6_modified_dataset_Gm14820___rand___full___RFmodel.RData")
aa_Gm14820_mlist_rand[[1]] <- my_RF_model
aa_Gm14820_mlist_rand_featImp[[1]] <- importance(my_RF_model)
names(aa_Gm14820_mlist_rand_featImp[[1]]) <- my_name_dic[match(names(aa_Gm14820_mlist_rand_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_RBP_Filtered_Ablation/Partition_6_modified_dataset_Gm14820___rand___first_129_last_1385___RFmodel.RData")
aa_Gm14820_mlist_rand[[2]] <- my_RF_model
aa_Gm14820_mlist_rand_featImp[[2]] <- importance(my_RF_model)
names(aa_Gm14820_mlist_rand_featImp[[2]]) <- my_name_dic[match(names(aa_Gm14820_mlist_rand_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_RBP_Filtered_Ablation/Partition_6_modified_dataset_Gm14820___rand___first_129_last_1389___RFmodel.RData")
aa_Gm14820_mlist_rand[[3]] <- my_RF_model
aa_Gm14820_mlist_rand_featImp[[3]] <- importance(my_RF_model)
names(aa_Gm14820_mlist_rand_featImp[[3]]) <- my_name_dic[match(names(aa_Gm14820_mlist_rand_featImp[[3]]), my_name_dic[, 2]),1]

# load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_OnlyU1/Partition_6_modified_dataset_Gm14820___rand___onlyU1___RFmodel.RData")
# aa_Gm14820_mlist_rand[[4]] <- my_RF_model
# aa_Gm14820_mlist_rand_featImp[[4]] <- importance(my_RF_model)
# names(aa_Gm14820_mlist_rand_featImp[[4]]) <- my_name_dic[match(names(aa_Gm14820_mlist_rand_featImp[[4]]), my_name_dic[, 2]),1]

names(aa_Gm14820_mlist_rand) <- c("All", "seq_only", "kmer+U1")

aa_pr_multimodel_P6_Gm14820_rand <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Gm14820_mlist_rand,
                                                                 my_dataset = aa_test_data,
                                                                 my_label = aa_test_data[, ncol(aa_test_data)],
                                                                 file_name = paste0(aaplot_strore_address,"Gm14820_All_vs_seq_vs_kmer","__rand___test_all_P6.png"),
                                                                 train_perf = F,
                                                                 neg_value = "neg",
                                                                 pos_value = "pos")

my_partition_df <- my_partition_chunk
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Gm14820_mlist_chunk <- list()
aa_Gm14820_mlist_chunk_featImp <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_RBP_Filtered/Partition_6_modified_dataset_Gm14820___chunk___full___RFmodel.RData")
aa_Gm14820_mlist_chunk[[1]] <- my_RF_model
aa_Gm14820_mlist_chunk_featImp[[1]] <- importance(my_RF_model)
names(aa_Gm14820_mlist_chunk_featImp[[1]]) <- my_name_dic[match(names(aa_Gm14820_mlist_chunk_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_RBP_Filtered_Ablation/Partition_6_modified_dataset_Gm14820___chunk___first_129_last_1385___RFmodel.RData")
aa_Gm14820_mlist_chunk[[2]] <- my_RF_model
aa_Gm14820_mlist_chunk_featImp[[2]] <- importance(my_RF_model)
names(aa_Gm14820_mlist_chunk_featImp[[2]]) <- my_name_dic[match(names(aa_Gm14820_mlist_chunk_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_RBP_Filtered_Ablation/Partition_6_modified_dataset_Gm14820___chunk___first_129_last_1389___RFmodel.RData")
aa_Gm14820_mlist_chunk[[3]] <- my_RF_model
aa_Gm14820_mlist_chunk_featImp[[3]] <- importance(my_RF_model)
names(aa_Gm14820_mlist_chunk_featImp[[3]]) <- my_name_dic[match(names(aa_Gm14820_mlist_chunk_featImp[[3]]), my_name_dic[, 2]),1]

# load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_OnlyU1/Partition_6_modified_dataset_Gm14820___chunk___onlyU1___RFmodel.RData")
# aa_Gm14820_mlist_chunk[[4]] <- my_RF_model
# aa_Gm14820_mlist_chunk_featImp[[4]] <- importance(my_RF_model)
# names(aa_Gm14820_mlist_chunk_featImp[[4]]) <- my_name_dic[match(names(aa_Gm14820_mlist_chunk_featImp[[4]]), my_name_dic[, 2]),1]


names(aa_Gm14820_mlist_chunk) <- c("All", "seq_only", "kmer+U1")
aa_pr_multimodel_P6_Gm14820_chunk <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Gm14820_mlist_chunk,
                                                                  my_dataset = aa_test_data,
                                                                  my_label = aa_test_data[, ncol(aa_test_data)],
                                                                  file_name = paste0(aaplot_strore_address,"Gm14820_All_vs_seq_vs_kmer","__chunk___test_all_P6.png"),
                                                                  train_perf = F,
                                                                  neg_value = "neg",
                                                                  pos_value = "pos")


png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_RBPfiltered/Gm14820_feature_rand.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Gm14820_mlist_rand_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Gm14820_mlist_rand_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Gm14820_mlist_rand_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmer+U1")
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_RBPfiltered/Gm14820_feature_chunk.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Gm14820_mlist_chunk_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Gm14820_mlist_chunk_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Gm14820_mlist_chunk_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmer+U1")
dev.off()

######################################################

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Kcnq1ot1.RData")
my_Dataset$U1snRNA_site <-  Partition_6_U1snRNA[match(rownames(my_Dataset), names(Partition_6_U1snRNA))]
my_Dataset <- my_Dataset[, c(c(1:(ncol(my_Dataset)-2)), (ncol(my_Dataset)), (ncol(my_Dataset) - 1))]
my_name_dic <- rbind(my_name_dic, c("U1snRNA_site", "U1snRNA_site"))
my_name_dic <- my_name_dic[c(c(1:(nrow(my_name_dic) - 2)), nrow(my_name_dic), (nrow(my_name_dic) - 1)), ]
my_Dataset <- Filter_RBP(inp_dataset = my_Dataset,
                         inp_nameDic = my_name_dic,
                         partition_Exp_features_address = "~/Documents/Shayan/BioInf/lncRNA/partition_6_expression_features.RData")


my_partition_df <- my_partition_rand
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Kcnq1ot1_mlist_nameDic <- my_name_dic
aa_Kcnq1ot1_mlist_rand_featImp <- list()
aa_Kcnq1ot1_mlist_rand <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_RBP_Filtered/Partition_6_modified_dataset_Kcnq1ot1___rand___full___RFmodel.RData")
aa_Kcnq1ot1_mlist_rand[[1]] <- my_RF_model
aa_Kcnq1ot1_mlist_rand_featImp[[1]] <- importance(my_RF_model)
names(aa_Kcnq1ot1_mlist_rand_featImp[[1]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_rand_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_RBP_Filtered_Ablation/Partition_6_modified_dataset_Kcnq1ot1___rand___first_200_last_1453___RFmodel.RData")
aa_Kcnq1ot1_mlist_rand[[2]] <- my_RF_model
aa_Kcnq1ot1_mlist_rand_featImp[[2]] <- importance(my_RF_model)
names(aa_Kcnq1ot1_mlist_rand_featImp[[2]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_rand_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_RBP_Filtered_Ablation/Partition_6_modified_dataset_Kcnq1ot1___rand___first_200_last_1457___RFmodel.RData")
aa_Kcnq1ot1_mlist_rand[[3]] <- my_RF_model
aa_Kcnq1ot1_mlist_rand_featImp[[3]] <- importance(my_RF_model)
names(aa_Kcnq1ot1_mlist_rand_featImp[[3]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_rand_featImp[[3]]), my_name_dic[, 2]),1]

# load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_OnlyU1/Partition_6_modified_dataset_Kcnq1ot1___rand___onlyU1___RFmodel.RData")
# aa_Kcnq1ot1_mlist_rand[[4]] <- my_RF_model
# aa_Kcnq1ot1_mlist_rand_featImp[[4]] <- importance(my_RF_model)
# names(aa_Kcnq1ot1_mlist_rand_featImp[[4]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_rand_featImp[[4]]), my_name_dic[, 2]),1]

names(aa_Kcnq1ot1_mlist_rand) <- c("All", "seq_only", "kmer+U1")

aa_pr_multimodel_P6_Kcnq1ot1_rand <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Kcnq1ot1_mlist_rand,
                                                                  my_dataset = aa_test_data,
                                                                  my_label = aa_test_data[, ncol(aa_test_data)],
                                                                  file_name = paste0(aaplot_strore_address,"Kcnq1ot1_All_vs_seq_vs_kmer","__rand___test_all_P6.png"),
                                                                  train_perf = F,
                                                                  neg_value = "neg",
                                                                  pos_value = "pos")

my_partition_df <- my_partition_chunk
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Kcnq1ot1_mlist_chunk <- list()
aa_Kcnq1ot1_mlist_chunk_featImp <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_RBP_Filtered/Partition_6_modified_dataset_Kcnq1ot1___chunk___full___RFmodel.RData")
aa_Kcnq1ot1_mlist_chunk[[1]] <- my_RF_model
aa_Kcnq1ot1_mlist_chunk_featImp[[1]] <- importance(my_RF_model)
names(aa_Kcnq1ot1_mlist_chunk_featImp[[1]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_chunk_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_RBP_Filtered_Ablation/Partition_6_modified_dataset_Kcnq1ot1___chunk___first_200_last_1453___RFmodel.RData")
aa_Kcnq1ot1_mlist_chunk[[2]] <- my_RF_model
aa_Kcnq1ot1_mlist_chunk_featImp[[2]] <- importance(my_RF_model)
names(aa_Kcnq1ot1_mlist_chunk_featImp[[2]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_chunk_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_RBP_Filtered_Ablation/Partition_6_modified_dataset_Kcnq1ot1___chunk___first_200_last_1457___RFmodel.RData")
aa_Kcnq1ot1_mlist_chunk[[3]] <- my_RF_model
aa_Kcnq1ot1_mlist_chunk_featImp[[3]] <- importance(my_RF_model)
names(aa_Kcnq1ot1_mlist_chunk_featImp[[3]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_chunk_featImp[[3]]), my_name_dic[, 2]),1]

# load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_OnlyU1/Partition_6_modified_dataset_Kcnq1ot1___chunk___onlyU1___RFmodel.RData")
# aa_Kcnq1ot1_mlist_chunk[[4]] <- my_RF_model
# aa_Kcnq1ot1_mlist_chunk_featImp[[4]] <- importance(my_RF_model)
# names(aa_Kcnq1ot1_mlist_chunk_featImp[[4]]) <- my_name_dic[match(names(aa_Kcnq1ot1_mlist_chunk_featImp[[4]]), my_name_dic[, 2]),1]

names(aa_Kcnq1ot1_mlist_chunk) <- c("All", "seq_only", "kmer+U1")
aa_pr_multimodel_P6_Kcnq1ot1_chunk <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Kcnq1ot1_mlist_chunk,
                                                                   my_dataset = aa_test_data,
                                                                   my_label = aa_test_data[, ncol(aa_test_data)],
                                                                   file_name = paste0(aaplot_strore_address,"Kcnq1ot1_All_vs_seq_vs_kmer","__chunk___test_all_P6.png"),
                                                                   train_perf = F,
                                                                   neg_value = "neg",
                                                                   pos_value = "pos")


png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_RBPfiltered/Kcnq1ot1_feature_rand.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Kcnq1ot1_mlist_rand_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Kcnq1ot1_mlist_rand_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Kcnq1ot1_mlist_rand_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmer+U1")
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_RBPfiltered/Kcnq1ot1_feature_chunk.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Kcnq1ot1_mlist_chunk_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Kcnq1ot1_mlist_chunk_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Kcnq1ot1_mlist_chunk_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmer+U1")
dev.off()


######################################################

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Trerf1.RData")
my_Dataset$U1snRNA_site <-  Partition_6_U1snRNA[match(rownames(my_Dataset), names(Partition_6_U1snRNA))]
my_Dataset <- my_Dataset[, c(c(1:(ncol(my_Dataset)-2)), (ncol(my_Dataset)), (ncol(my_Dataset) - 1))]
my_name_dic <- rbind(my_name_dic, c("U1snRNA_site", "U1snRNA_site"))
my_name_dic <- my_name_dic[c(c(1:(nrow(my_name_dic) - 2)), nrow(my_name_dic), (nrow(my_name_dic) - 1)), ]
my_Dataset <- Filter_RBP(inp_dataset = my_Dataset,
                         inp_nameDic = my_name_dic,
                         partition_Exp_features_address = "~/Documents/Shayan/BioInf/lncRNA/partition_6_expression_features.RData")


my_partition_df <- my_partition_rand
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Trerf1_mlist_nameDic <- my_name_dic
aa_Trerf1_mlist_rand_featImp <- list()
aa_Trerf1_mlist_rand <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_RBP_Filtered/Partition_6_modified_dataset_Trerf1___rand___full___RFmodel.RData")
aa_Trerf1_mlist_rand[[1]] <- my_RF_model
aa_Trerf1_mlist_rand_featImp[[1]] <- importance(my_RF_model)
names(aa_Trerf1_mlist_rand_featImp[[1]]) <- my_name_dic[match(names(aa_Trerf1_mlist_rand_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_RBP_Filtered_Ablation/Partition_6_modified_dataset_Trerf1___rand___first_173_last_1425___RFmodel.RData")
aa_Trerf1_mlist_rand[[2]] <- my_RF_model
aa_Trerf1_mlist_rand_featImp[[2]] <- importance(my_RF_model)
names(aa_Trerf1_mlist_rand_featImp[[2]]) <- my_name_dic[match(names(aa_Trerf1_mlist_rand_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_RBP_Filtered_Ablation/Partition_6_modified_dataset_Trerf1___rand___first_173_last_1429___RFmodel.RData")
aa_Trerf1_mlist_rand[[3]] <- my_RF_model
aa_Trerf1_mlist_rand_featImp[[3]] <- importance(my_RF_model)
names(aa_Trerf1_mlist_rand_featImp[[3]]) <- my_name_dic[match(names(aa_Trerf1_mlist_rand_featImp[[3]]), my_name_dic[, 2]),1]

# load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_OnlyU1/Partition_6_modified_dataset_Trerf1___rand___onlyU1___RFmodel.RData")
# aa_Trerf1_mlist_rand[[4]] <- my_RF_model
# aa_Trerf1_mlist_rand_featImp[[4]] <- importance(my_RF_model)
# names(aa_Trerf1_mlist_rand_featImp[[4]]) <- my_name_dic[match(names(aa_Trerf1_mlist_rand_featImp[[4]]), my_name_dic[, 2]),1]

names(aa_Trerf1_mlist_rand) <- c("All", "seq_only", "kmer+U1")

aa_pr_multimodel_P6_Trerf1_rand <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Trerf1_mlist_rand,
                                                                my_dataset = aa_test_data,
                                                                my_label = aa_test_data[, ncol(aa_test_data)],
                                                                file_name = paste0(aaplot_strore_address,"Trerf1_All_vs_seq_vs_kmer","__rand___test_all_P6.png"),
                                                                train_perf = F,
                                                                neg_value = "neg",
                                                                pos_value = "pos")

my_partition_df <- my_partition_chunk
aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- my_Dataset[aatrainind,]
aa_test_data <- my_Dataset[-aatrainind,]

aa_Trerf1_mlist_chunk <- list()
aa_Trerf1_mlist_chunk_featImp <- list()
load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_RBP_Filtered/Partition_6_modified_dataset_Trerf1___chunk___full___RFmodel.RData")
aa_Trerf1_mlist_chunk[[1]] <- my_RF_model
aa_Trerf1_mlist_chunk_featImp[[1]] <- importance(my_RF_model)
names(aa_Trerf1_mlist_chunk_featImp[[1]]) <- my_name_dic[match(names(aa_Trerf1_mlist_chunk_featImp[[1]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_RBP_Filtered_Ablation/Partition_6_modified_dataset_Trerf1___chunk___first_173_last_1425___RFmodel.RData")
aa_Trerf1_mlist_chunk[[2]] <- my_RF_model
aa_Trerf1_mlist_chunk_featImp[[2]] <- importance(my_RF_model)
names(aa_Trerf1_mlist_chunk_featImp[[2]]) <- my_name_dic[match(names(aa_Trerf1_mlist_chunk_featImp[[2]]), my_name_dic[, 2]),1]

load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_RBP_Filtered_Ablation/Partition_6_modified_dataset_Trerf1___chunk___first_173_last_1429___RFmodel.RData")
aa_Trerf1_mlist_chunk[[3]] <- my_RF_model
aa_Trerf1_mlist_chunk_featImp[[3]] <- importance(my_RF_model)
names(aa_Trerf1_mlist_chunk_featImp[[3]]) <- my_name_dic[match(names(aa_Trerf1_mlist_chunk_featImp[[3]]), my_name_dic[, 2]),1]

# load("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_OnlyU1/Partition_6_modified_dataset_Trerf1___chunk___onlyU1___RFmodel.RData")
# aa_Trerf1_mlist_chunk[[4]] <- my_RF_model
# aa_Trerf1_mlist_chunk_featImp[[4]] <- importance(my_RF_model)
# names(aa_Trerf1_mlist_chunk_featImp[[4]]) <- my_name_dic[match(names(aa_Trerf1_mlist_chunk_featImp[[4]]), my_name_dic[, 2]),1]

names(aa_Trerf1_mlist_chunk) <- c("All", "seq_only", "kmer+U1")
aa_pr_multimodel_P6_Trerf1_chunk <- perf_eval_ROC_PRC_multimodel(model_fit_list = aa_Trerf1_mlist_chunk,
                                                                 my_dataset = aa_test_data,
                                                                 my_label = aa_test_data[, ncol(aa_test_data)],
                                                                 file_name = paste0(aaplot_strore_address,"Trerf1_All_vs_seq_vs_kmer","__chunk___test_all_P6.png"),
                                                                 train_perf = F,
                                                                 neg_value = "neg",
                                                                 pos_value = "pos")


png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_RBPfiltered/Trerf1_feature_rand.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Trerf1_mlist_rand_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Trerf1_mlist_rand_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Trerf1_mlist_rand_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmer+U1")
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/performance_plots_p6_FeatureAblation_comparison_RBPfiltered/Trerf1_feature_chunk.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,3), mar = c(4,9,4,4))
barplot(sort(aa_Trerf1_mlist_chunk_featImp[[1]], decreasing = T)[1:30], horiz = T, las =2, main = "Full")
barplot(sort(aa_Trerf1_mlist_chunk_featImp[[2]], decreasing = T)[1:30], horiz = T, las =2, main = "SeqOnly")
barplot(sort(aa_Trerf1_mlist_chunk_featImp[[3]], decreasing = T)[1:30], horiz = T, las =2, main = "kmer+U1")
dev.off()

#################################3

################################

###############################
#perf_eval_ROC_PRC_multimodel

#zoom
source("~/Documents/Shayan/BioInf/lncRNA/RF_evaluation_functions.R")
source("~/Documents/Shayan/BioInf/lncRNA/RBP_filter_func.R")
aadir <- list.dirs("~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Results/Learned_models")
system(paste0("cat ~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full.job ~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full2.job ~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full_TXGRO.job ~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full_TXrandomized.job run_rf_p6_dec_full_ONLYGROPRO.job ~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full_AddedGROPRO.job > ~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full_6.job"))
performance_list_lncRNA <- list()

for(i in 2:3){
  performance_list_lncRNA[[i-1]] <- perf_eval_wrapper(Learned_model_folder = aadir[i],only_imporatnce = F,save_pred = T,
                                                    lncRNA_dataset_folder = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_dataset",
                                                    model_name=character(0),
                                                    lncRNA_name = character(0),
                                                    jobfile_address = "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full_6.job",
                                                    aaplot_strore_address = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Res_comparision",
                                                    filtering_folder = "~/Documents/Shayan/BioInf/lncRNA/partition_6_Dec_TOBE_exp_FILTERED")
}

performance_pred_list_lncRNA <- list()
length(aadir)
for(i in 2:length(aadir)){
  performance_pred_list_lncRNA[[i-1]] <- perf_eval_wrapper(Learned_model_folder = aadir[i],only_imporatnce = F,save_pred = T,
                                                      lncRNA_dataset_folder = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_dataset",
                                                      model_name=character(0),
                                                      lncRNA_name = character(0),
                                                      jobfile_address = "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full_6.job",
                                                      aaplot_strore_address = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Res_comparision",
                                                      filtering_folder = "~/Documents/Shayan/BioInf/lncRNA/partition_6_Dec_TOBE_exp_FILTERED")
}
# not finished for Gm11613 and Gm53
# i <- 14 
# aatst <- perf_eval_wrapper(Learned_model_folder = aadir[i],only_imporatnce = F,
#                            lncRNA_dataset_folder = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_dataset",
#                            model_name=character(0),
#                            lncRNA_name = character(0),
#                            jobfile_address = "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full_three.job",
#                            aaplot_strore_address = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Res_comparision",
#                            filtering_folder = "~/Documents/Shayan/BioInf/lncRNA/partition_6_Dec_TOBE_exp_FILTERED")

aasp <- strsplit(aadir, "\\/")
aasp1 <- unlist(lapply(aasp, function(x) x[length(x)]))
aasp1 <- aasp1[2:length(aasp1)]
#performance_list_lncRNA <- performance_list_lncRNA[2:length(performance_list_lncRNA)]
names(performance_list_lncRNA) <- aasp1
names(performance_pred_list_lncRNA) <- aasp1

save(list = c("performance_list_lncRNA"), file = "~/Documents/Shayan/BioInf/lncRNA/performance_list_lncRNA.RData")
save(list = c("performance_pred_list_lncRNA"), file = "~/Documents/Shayan/BioInf/lncRNA/performance_pred_list_lncRNA.RData")

# aadir
# load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_dataset/Partition_6_Dec_dataset_Tsix.RData")
# load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Results/Learned_models/Tsix/Partition_6_Dec_dataset_Tsix___chunk___RBPscannedTX___RFmodel.RData")
# aamy_filter_index <- as.numeric(read.table("~/Documents/Shayan/BioInf/lncRNA/partition_6_Dec_TOBE_exp_FILTERED/my_filter_2.txt")$V1)
# 
# my_name_dic[aamy_filter_index, 1]
# aamy_Dataset <- Filter_RBP(inp_dataset = my_Dataset,
#                          inp_nameDic = my_name_dic,
#                          feature_index=aamy_filter_index,
#                          partition_Exp_features_address = "~/Documents/Shayan/BioInf/lncRNA/partition_6_expression_features.RData")
# 
# 
# rownames(my_Dataset)[which(! rownames(my_Dataset) %in% partition_6_expression_features$tile_name)]
# sum(! rownames(my_Dataset) %in% Partition_6_dfs$tile_name)
# 
# 
# 
# aaexp1 <- partition_6_expression_features$RNAseq[match(rownames(my_Dataset), partition_6_expression_features$tile_name)]
# aaexp2 <- partition_6_expression_features$CAGE[match(rownames(my_Dataset), partition_6_expression_features$tile_name)]
# aaexp3 <- partition_6_expression_features$RNA_polymerase_II__315.bed[match(rownames(my_Dataset), partition_6_expression_features$tile_name)]

# Plot Various model comparisions for Malat1


aamalat1_allnames <- list.files("~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Results/Learned_models/Malat1")
aamalat1_allnames2 <- unlist(lapply(strsplit(aamalat1_allnames, "___"), "[[", 2))
aamalat1_allnames3 <- unlist(lapply(strsplit(aamalat1_allnames, "___"), "[[", 3))
length(unique(aamalat1_allnames3))

aa_name_list <- list()
aa_name_list[[1]] <- unique(aamalat1_allnames3[grep("Kmer", aamalat1_allnames3)])
aa_name_list[[2]] <- unique(aamalat1_allnames3[grep("ChIP", aamalat1_allnames3)])
aa_name_list[[3]] <- unique(aamalat1_allnames3[grep("Chrom", aamalat1_allnames3)])
aa_name_list[[4]] <- unique(aamalat1_allnames3[grep("RBPscan", aamalat1_allnames3)])
aa_name_list[[5]] <- unique(aamalat1_allnames3[grep("TFscan", aamalat1_allnames3)])
aa_name_list[[6]] <- unique(aamalat1_allnames3[grep("Repeat", aamalat1_allnames3)])

for(i in 1:length(aa_name_list)){
  perf_eval_wrapper(Learned_model_folder = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Results/Learned_models/Malat1",
                    only_imporatnce = F,
                    lncRNA_dataset_folder = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_dataset",
                    model_name=character(0),
                    lncRNA_name = character(0),
                    jobfile_address = "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full_three.job",
                    aaplot_strore_address = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Res_comparision_malat1",
                    filtering_folder = "~/Documents/Shayan/BioInf/lncRNA/partition_6_Dec_TOBE_exp_FILTERED")
  
}



###############
## Comparing the performance of various families of classes in contrast to each other
load("~/Documents/Shayan/BioInf/lncRNA/performance_pred_list_lncRNA.RData")

#aa1 <- unlist(lapply(performance_list_lncRNA, length))
aa2_1 <- unlist(lapply(lapply(performance_pred_list_lncRNA, "[[", 1), length))
aa2_2 <- unlist(lapply(lapply(performance_pred_list_lncRNA, "[[", 2), length))

performance_auroc_rand_28_34 <- matrix(nrow = length(aa2_1), ncol = max(aa2_1))
performance_auroc_chunk_28_34 <- matrix(nrow = length(aa2_1), ncol = max(aa2_1))
performance_auprc_rand_28_34 <- matrix(nrow = length(aa2_1), ncol = max(aa2_1))
performance_auprc_chunk_28_34 <- matrix(nrow = length(aa2_1), ncol = max(aa2_1))
aawhichmax <- which.max(aa2_1)
rownames(performance_auroc_rand_28_34) <- names(aa2_1)
colnames(performance_auroc_rand_28_34) <- names(performance_pred_list_lncRNA[[aawhichmax]][[1]])
rownames(performance_auroc_chunk_28_34) <- names(aa2_1)
colnames(performance_auroc_chunk_28_34) <- names(performance_pred_list_lncRNA[[aawhichmax]][[1]])
rownames(performance_auprc_rand_28_34) <- names(aa2_1)
colnames(performance_auprc_rand_28_34) <- names(performance_pred_list_lncRNA[[aawhichmax]][[1]])
rownames(performance_auprc_chunk_28_34) <- names(aa2_1)
colnames(performance_auprc_chunk_28_34) <- names(performance_pred_list_lncRNA[[aawhichmax]][[1]])

for(i in 1:length(performance_pred_list_lncRNA)){
  #rand
  for(j in 1:aa2_1[i]){
    aatmp1 <- performance_pred_list_lncRNA[[i]][[1]][[j]]$auc
    aatmp2 <- performance_pred_list_lncRNA[[i]][[1]][[j]]$prc$auc.integral
    performance_auroc_rand_28_34[i, match(names(performance_pred_list_lncRNA[[i]][[1]])[j], colnames(performance_auroc_rand_28_34))] <- aatmp1
    performance_auprc_rand_28_34[i, match(names(performance_pred_list_lncRNA[[i]][[1]])[j], colnames(performance_auprc_rand_28_34))] <- aatmp2
  }
  #chunk
  for(j in 1:aa2_2[i]){
    aatmp1 <- performance_pred_list_lncRNA[[i]][[2]][[j]]$auc
    aatmp2 <- performance_pred_list_lncRNA[[i]][[2]][[j]]$prc$auc.integral
    
    performance_auroc_chunk_28_34[i, match(names(performance_pred_list_lncRNA[[i]][[2]])[j], colnames(performance_auroc_chunk_28_34))] <- aatmp1
    performance_auprc_chunk_28_34[i, match(names(performance_pred_list_lncRNA[[i]][[2]])[j], colnames(performance_auprc_chunk_28_34))] <- aatmp2
  }
}
# 1) plot scatter plots with full model performance on one axis, and performance of a feature family on the other axis

#  comparisons:
# RBPTFRepeatChromChIPexpKmer vs ChIP
# RBPTFRepeatChromChIPexpKmer vs Chromatin
# RBPTFRepeatChromChIPexpKmer vs Transcription
# RBPTFRepeatChromChIPexpKmer vs ChromChIPexp
# RBPTFRepeatChromChIPexpKmer vs Kmer
# RBPTFRepeatChromChIPexpKmer vs RBPscanned
# RBPTFRepeatChromChIPexpKmer vs Repeat
# RBPTFRepeatChromChIPexpKmer vs TFscanned
# RBPTFRepeatChromChIPexpKmer vs RBPTFRepeat
# RBPTFRepeatChromChIPexpKmer vs Triplex
# RBPTFRepeatChromChIPexpKmer vs RBPTFRepeatTriplex
# RBPTFRepeatChromChIPexpKmer vs Pairs
# RBPTFRepeatChromChIPexpKmer vs RBPTFRepeatChromChIPexp
# RBPTFRepeatChromChIPexpKmer vs KmerChromChipExp
# RBPTFRepeatChromChIPexpKmer vs KmerTriplexPairs

aaCompMat <- rbind(c("RBPTFRepeatChromChIPexpKmer" , "ChIP"),
                   c("RBPTFRepeatChromChIPexpKmer" , "Chromatin"),
                   c("RBPTFRepeatChromChIPexpKmer" , "Transcription"),
                   c("RBPTFRepeatChromChIPexpKmer" , "ChromChIPexp"),
                   c("RBPTFRepeatChromChIPexpKmer" , "Kmer"),
                   c("RBPTFRepeatChromChIPexpKmer" , "RBPscanned"),
                   c("RBPTFRepeatChromChIPexpKmer" , "Repeat"),
                   c("RBPTFRepeatChromChIPexpKmer" , "TFscanned"),
                   c("RBPTFRepeatChromChIPexpKmer" , "RBPTFRepeat"),
                   c("RBPTFRepeatChromChIPexpKmer" , "Triplex"),
                   c("RBPTFRepeatChromChIPexpKmer" , "RBPTFRepeatTriplex"),
                   c("RBPTFRepeatChromChIPexpKmer" , "Pairs"),
                   c("RBPTFRepeatChromChIPexpKmer" , "RBPTFRepeatChromChIPexp"),
                   c("RBPTFRepeatChromChIPexpKmer" , "KmerChromChipExp"),
                   c("RBPTFRepeatChromChIPexpKmer" , "KmerTriplexPairs"),
                   c("RBPTFRepeatChromChIPexpKmer" , "GROPRO")
                   )

library(RColorBrewer)
n <- 28
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector_Samp <- sample(col_vector, n)
pie(rep(1,n), col=col_vector_Samp)
names(col_vector_Samp) <- rownames(performance_auroc_chunk_28_34)

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Res_comparision/aggreg_family_vs_full_rand_auroc.png", 
    width = 10*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(4,4), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = performance_auroc_rand_28_34[, match(aaCompMat[i,1], colnames(performance_auroc_rand_28_34))],
       y = performance_auroc_rand_28_34[, match(aaCompMat[i,2], colnames(performance_auroc_rand_28_34))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.3,0.85), ylim = c(0.3,0.85),
       cex.lab = 0.7,
       cex= 1.2,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)

    legend(x = "top", inset=c(0,-0.7),xpd = T,ncol=3,
           legend = names(col_vector_Samp), 
           fill = col_vector_Samp,
           bty = "n",xjust = 0.5, yjust = 0,
           cex = 0.7, x.intersp = 0.5, y.intersp = 0.7) 
    
}
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Res_comparision/aggreg_family_vs_full_chunk_auroc.png", 
    width = 10*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(4,4), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = performance_auroc_chunk_28_34[, match(aaCompMat[i,1], colnames(performance_auroc_chunk_28_34))],
       y = performance_auroc_chunk_28_34[, match(aaCompMat[i,2], colnames(performance_auroc_chunk_28_34))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.3,0.85), ylim = c(0.3,0.85),
       cex.lab = 0.7,
       cex= 1.2,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.5, y.intersp = 0.7) 
  
}
dev.off()

range(performance_auprc_rand_28_34, na.rm = T)
range(performance_auprc_chunk_28_34, na.rm = T)

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Res_comparision/aggreg_family_vs_full_rand_auprc.png", 
    width = 10*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(4,4), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = performance_auprc_rand_28_34[, match(aaCompMat[i,1], colnames(performance_auprc_rand_28_34))],
       y = performance_auprc_rand_28_34[, match(aaCompMat[i,2], colnames(performance_auprc_rand_28_34))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.1,0.85), ylim = c(0.1,0.85),
       cex.lab = 0.7,
       cex= 1.2,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.5, y.intersp = 0.7) 
  
}
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Res_comparision/aggreg_family_vs_full_chunk_auprc.png", 
    width = 10*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(4,4), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = performance_auprc_chunk_28_34[, match(aaCompMat[i,1], colnames(performance_auprc_chunk_28_34))],
       y = performance_auprc_chunk_28_34[, match(aaCompMat[i,2], colnames(performance_auprc_chunk_28_34))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.1,0.85), ylim = c(0.1,0.85),
       cex.lab = 0.7,
       cex= 1.2,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.5, y.intersp = 0.7) 
  
}
dev.off()


# 2) plot performance of specific families in comparison to each other, for example sequence against cell context, or kmer against TBPTFREPEAT

#  comparisons:
# Kmer vs ChromChIPexp
# RBPTFRepeat vs ChromChIPexp
# RBPTFRepeat vs Kmer
# RBPscanned vs TFscanned
# RBPscanned vs Repeat
# TFscanned vs Repeat
# RBPscannedTX vs TFscannedTX
# RBPscannedTX vs RepeatTX
# TFscannedTX vs RepeatTX
# RBPscannedTXGRO vs TFscannedTXGRO
# RBPscannedTXGRO vs RepeatTXGRO
# TFscannedTXGRO vs RepeatTXGRO



aaCompMat <- rbind(c("Kmer" , "ChromChIPexp"),
                   c("RBPTFRepeat" , "ChromChIPexp"),
                   c("RBPTFRepeat" , "Kmer"),
                   c("RBPscanned" , "TFscanned"),
                   c("RBPscanned" , "Repeat"),
                   c("TFscanned" , "Repeat"),
                   c("RBPscanned" , "Kmer"),
                   c("TFscanned" , "Kmer"),
                   c("Repeat" , "Kmer"),
                   c("RBPscannedTX" , "TFscannedTX"),
                   c("RBPscannedTX" , "RepeatTX"),
                   c("TFscannedTX" , "RepeatTX"),
                   c("RBPscannedTXGRO" , "TFscannedTXGRO"),
                   c("RBPscannedTXGRO" , "RepeatTXGRO"),
                   c("TFscannedTXGRO" , "RepeatTXGRO")
                   )


png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Res_comparision/aggreg_family_vs_family_rand_auroc.png", 
    width = 10*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(4,4), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = performance_auroc_rand_28_34[, match(aaCompMat[i,1], colnames(performance_auroc_rand_28_34))],
       y = performance_auroc_rand_28_34[, match(aaCompMat[i,2], colnames(performance_auroc_rand_28_34))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.3,0.85), ylim = c(0.3,0.85),
       cex.lab = 0.7,
       cex= 1.2,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.5, y.intersp = 0.7) 
  
}
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Res_comparision/aggreg_family_vs_family_chunk_auroc.png", 
    width = 10*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(4,4), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = performance_auroc_chunk_28_34[, match(aaCompMat[i,1], colnames(performance_auroc_chunk_28_34))],
       y = performance_auroc_chunk_28_34[, match(aaCompMat[i,2], colnames(performance_auroc_chunk_28_34))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.3,0.85), ylim = c(0.3,0.85),
       cex.lab = 0.7,
       cex= 1.2,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.5, y.intersp = 0.7) 
  
}
dev.off()

range(performance_auprc_rand_28_34, na.rm = T)
range(performance_auprc_chunk_28_34, na.rm = T)

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Res_comparision/aggreg_family_vs_family_rand_auprc.png", 
    width = 10*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(4,4), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = performance_auprc_rand_28_34[, match(aaCompMat[i,1], colnames(performance_auprc_rand_28_34))],
       y = performance_auprc_rand_28_34[, match(aaCompMat[i,2], colnames(performance_auprc_rand_28_34))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.1,0.85), ylim = c(0.1,0.85),
       cex.lab = 0.7,
       cex= 1.2,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.5, y.intersp = 0.7) 
  
}
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Res_comparision/aggreg_family_vs_family_chunk_auprc.png", 
    width = 10*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(4,4), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = performance_auprc_chunk_28_34[, match(aaCompMat[i,1], colnames(performance_auprc_chunk_28_34))],
       y = performance_auprc_chunk_28_34[, match(aaCompMat[i,2], colnames(performance_auprc_chunk_28_34))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.1,0.85), ylim = c(0.1,0.85),
       cex.lab = 0.7,
       cex= 1.2,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.5, y.intersp = 0.7) 
  
}
dev.off()


# 3) plot expression-filtered features against their non filtered version --> for RBP-TF-REPEAT --> for both types of expression filter

#  comparisons:
# RBPscanned vs RBPscannedTX
# RBPscanned vs RBPscannedTXGRO
# RBPscannedTX vs RBPscannedTXGRO

# TFscanned vs TFscannedTX
# TFscanned vs TFscannedTXGRO
# TFscannedTX vs TFscannedTXGRO

# Repeat vs RepeatTX
# Repeat vs RepeatTXGRO
# RepeatTX vs RepeatTXGRO

# RBPTFRepeat vs RBPTFRepeatTXRBPrep
# RBPTFRepeat vs RBPTFRepeatTXGRORBPrep
# RBPTFRepeatTXRBPrep vs RBPTFRepeatTXGRORBPrep

# RBPTFRepeatTriplex vs RBPTFRepeatTriplexTXRBPrep
# RBPTFRepeatTriplex vs RBPTFRepeatTriplexTXGRORBPrep
# RBPTFRepeatTriplexTXRBPrep vs RBPTFRepeatTriplexTXGRORBPrep

# RBPTFRepeatChromChIPexp vs RBPTFRepeatTXRBPrepChromChIPexp
# RBPTFRepeatChromChIPexp vs RBPTFRepeatTXGRORBPrepChromChIPexp
# RBPTFRepeatTXRBPrepChromChIPexp vs RBPTFRepeatTXGRORBPrepChromChIPexp

# RBPTFRepeatChromChIPexpKmer vs RBPTFRepeatTXRBPrepChromChIPexpKmer
# RBPTFRepeatChromChIPexpKmer vs RBPTFRepeatTXGRORBPrepChromChIPexpKmer
# RBPTFRepeatTXRBPrepChromChIPexpKmer vs RBPTFRepeatTXGRORBPrepChromChIPexpKmer

aaCompMat <- rbind(c("RBPscanned" , "RBPscannedTX"),
                   c("RBPscanned" , "RBPscannedTXGRO"),
                   c("RBPscanned" , "RBPscannedTXRand"),
                   c("RBPscannedTX" , "RBPscannedTXGRO"),
                   c("TFscanned" , "TFscannedTX"),
                   c("TFscanned" , "TFscannedTXGRO"),
                   c("TFscanned" , "TFscannedTXRand"),
                   c("TFscannedTX" , "TFscannedTXGRO"),
                   c("Repeat" , "RepeatTX"),
                   c("Repeat" , "RepeatTXGRO"),
                   c("Repeat" , "RepeatTXRand"),
                   c("RepeatTX" , "RepeatTXGRO"),
                   c("RBPTFRepeat" , "RBPTFRepeatTXRBPrep"),
                   c("RBPTFRepeat" , "RBPTFRepeatTXGRORBPrep"),
                   c("RBPTFRepeatTXRBPrep" , "RBPTFRepeatTXGRORBPrep"),
                   c("RBPTFRepeatTriplex" , "RBPTFRepeatTriplexTXRBPrep"),
                   c("RBPTFRepeatTriplex" , "RBPTFRepeatTriplexTXGRORBPrep"),
                   c("RBPTFRepeatTriplexTXRBPrep" , "RBPTFRepeatTriplexTXGRORBPrep"),
                   c("RBPTFRepeatChromChIPexp" , "RBPTFRepeatTXRBPrepChromChIPexp"),
                   c("RBPTFRepeatChromChIPexp" , "RBPTFRepeatTXGRORBPrepChromChIPexp"),
                   c("RBPTFRepeatTXRBPrepChromChIPexp" , "RBPTFRepeatTXGRORBPrepChromChIPexp"),
                   c("RBPTFRepeatChromChIPexpKmer" , "RBPTFRepeatTXRBPrepChromChIPexpKmer"),
                   c("RBPTFRepeatChromChIPexpKmer" , "RBPTFRepeatTXGRORBPrepChromChIPexpKmer"),
                   c("RBPTFRepeatTXRBPrepChromChIPexpKmer" , "RBPTFRepeatTXGRORBPrepChromChIPexpKmer")
                   )



png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Res_comparision/aggreg_expFilter_rand_auroc.png", 
    width = 10*300,
    height = 15*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(6,4), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = performance_auroc_rand_28_34[, match(aaCompMat[i,1], colnames(performance_auroc_rand_28_34))],
       y = performance_auroc_rand_28_34[, match(aaCompMat[i,2], colnames(performance_auroc_rand_28_34))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.3,0.85), ylim = c(0.3,0.85),
       cex.lab = 0.7,
       cex= 1.2,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.5, y.intersp = 0.7) 
  
}
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Res_comparision/aggreg_expFilter_chunk_auroc.png", 
    width = 10*300,
    height = 15*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(6,4), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = performance_auroc_chunk_28_34[, match(aaCompMat[i,1], colnames(performance_auroc_chunk_28_34))],
       y = performance_auroc_chunk_28_34[, match(aaCompMat[i,2], colnames(performance_auroc_chunk_28_34))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.3,0.85), ylim = c(0.3,0.85),
       cex.lab = 0.7,
       cex= 1.2,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.5, y.intersp = 0.7) 
  
}
dev.off()

range(performance_auprc_rand_28_34, na.rm = T)
range(performance_auprc_chunk_28_34, na.rm = T)

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Res_comparision/aggreg_expFilter_rand_auprc.png", 
    width = 10*300,
    height = 15*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(6,4), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = performance_auprc_rand_28_34[, match(aaCompMat[i,1], colnames(performance_auprc_rand_28_34))],
       y = performance_auprc_rand_28_34[, match(aaCompMat[i,2], colnames(performance_auprc_rand_28_34))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.1,0.85), ylim = c(0.1,0.85),
       cex.lab = 0.7,
       cex= 1.2,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.5, y.intersp = 0.7) 
  
}
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Res_comparision/aggreg_expFilter_chunk_auprc.png", 
    width = 10*300,
    height = 15*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(6,4), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = performance_auprc_chunk_28_34[, match(aaCompMat[i,1], colnames(performance_auprc_chunk_28_34))],
       y = performance_auprc_chunk_28_34[, match(aaCompMat[i,2], colnames(performance_auprc_chunk_28_34))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.1,0.85), ylim = c(0.1,0.85),
       cex.lab = 0.7,
       cex= 1.2,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.5, y.intersp = 0.7) 
  
}
dev.off()


png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Res_comparision/aggreg_expFilter_chunk_auprc_diff.png", 
    width = 10*300,
    height = 15*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(6,4), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = performance_auprc_chunk_28_34_dif[, match(aaCompMat[i,1], colnames(performance_auprc_chunk_28_34_dif))],
       y = performance_auprc_chunk_28_34_dif[, match(aaCompMat[i,2], colnames(performance_auprc_chunk_28_34_dif))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.0,0.4), ylim = c(0.0,0.4),
       cex.lab = 0.7,
       cex= 1.2,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.5, y.intersp = 0.7) 
  
}
dev.off()


png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Res_comparision/aggreg_expFilter_rand_auprc_diff.png", 
    width = 10*300,
    height = 15*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(6,4), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = performance_auprc_rand_28_34_dif[, match(aaCompMat[i,1], colnames(performance_auprc_rand_28_34_dif))],
       y = performance_auprc_rand_28_34_dif[, match(aaCompMat[i,2], colnames(performance_auprc_rand_28_34_dif))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.0,0.45), ylim = c(0.0,0.45),
       cex.lab = 0.7,
       cex= 1.2,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.5, y.intersp = 0.7) 
  
}
dev.off()


par(mfrow = c(1,1), mar = c(18,4,4,4))
aaind <- sort(performance_auroc_chunk_28_34["Neat1",], decreasing = T, index.return = T)$ix
barplot(performance_auroc_chunk_28_34["Neat1",aaind], las = 2)


par(mfrow = c(1,1), mar = c(18,4,4,4))
aaind <- sort(performance_auroc_chunk_28_34["Malat1",], decreasing = T, index.return = T)$ix
barplot(performance_auroc_chunk_28_34["Malat1",aaind], las = 2)

aatestrand <- Partition_6_dfs[Partition_6_dfs$dataset == "test",]
aatestchun <- Partition_6_chunk_dfs[Partition_6_chunk_dfs$dataset == "test",]

aax <- as.matrix(table(aatestrand$label, aatestrand$owner))
aprcbaserand <- aax[2,]/colSums(aax)

aax2 <- as.matrix(table(aatestchun$label, aatestchun$owner))
aprcbasechun <- aax2[2,]/colSums(aax2)



performance_auprc_rand_28_34_dif <- performance_auprc_rand_28_34 - aprcbaserand[match(rownames(performance_auprc_rand_28_34), names(aprcbaserand))]
performance_auprc_chunk_28_34_dif <- performance_auprc_chunk_28_34 - aprcbasechun[match(rownames(performance_auprc_chunk_28_34), names(aprcbasechun))]
################################################################################################
# write jobs to gather performance eval




aa_un_own <- sort(unique(Partition_6_random_chunk_cv_df$owner))


for(i in 1:length(aa_un_own)){
  cat(c("Rscript --vanilla /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/gather_perf.R",
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CVfirst_Results/Learned_models/",aa_un_own[i] ),
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CVfirst_Results/Combined_results/",aa_un_own[i], "__perfMat.RData\n" )),
      sep = " ", 
      file = "~/Documents/Shayan/BioInf/lncRNA/gather_perf_all_cvfirst.job",
      append = !(i == 1))
}
################################################################################################
# write jobs to gather performance eval




aa_un_own <- sort(unique(Partition_6_random_chunk_cv_df$owner))


for(i in 1:length(aa_un_own)){
  cat(c("Rscript --vanilla /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/gather_perf.R",
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_svd1_Results/Learned_models/",aa_un_own[i] ),
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_svd1_Results/Combined_results/",aa_un_own[i], "__perfMat.RData\n" )),
      sep = " ", 
      file = "~/Documents/Shayan/BioInf/lncRNA/gather_perf_all.job",
      append = !(i == 1))
}
################################################################################################
# read combined performance for 28 lncRNAs (family svd)
aafiles <- list.files("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res")
aafn <- unlist(lapply(strsplit(aafiles, "__"), "[[", 1))
aa_all_auroc_list <- list()
aa_all_auprc_list <- list()
aa_all_impor_list <- list()

for(i in 1:length(aafiles)){
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res/", aafiles[i]))
  aa_all_auroc_list[[i]] <- my_perf[[1]][[1]]
  aa_all_auprc_list[[i]] <- my_perf[[2]][[1]]
  aa_all_impor_list[[i]] <- my_perf[[3]][[1]]
}
names(aa_all_auroc_list) <- aafn
names(aa_all_auprc_list) <- aafn
names(aa_all_impor_list) <- aafn

View(aa_all_auroc_list$Malat1)

aa_auprc_base <- matrix(nrow = length(aafn), ncol = (ncol(Partition_6_random_chunk_cv_df) - 3))
rownames(aa_auprc_base) <- aafn
colnames(aa_auprc_base) <- colnames(Partition_6_random_chunk_cv_df)[4:ncol(Partition_6_random_chunk_cv_df)]

for(i in 1:nrow(aa_auprc_base)){
  print(i)
  for(j in 1:ncol(aa_auprc_base)){
    aatst <- Partition_6_random_chunk_cv_df[Partition_6_random_chunk_cv_df[, 3+j] == 0, c(1:3)]
    
    aa_auprc_base[i, j] <- sum((aatst$label == 1) & (aatst$owner == aafn[i]))/sum(aatst$owner == aafn[i])
  }
}
View(aa_auprc_base)

aa_auprc_base <- aa_auprc_base[, c(c(6:10), c(1:5))]

aa_auprc_base_train <- matrix(nrow = length(aafn), ncol = (ncol(Partition_6_random_chunk_cv_df) - 3))
rownames(aa_auprc_base_train) <- aafn
colnames(aa_auprc_base_train) <- colnames(Partition_6_random_chunk_cv_df)[4:ncol(Partition_6_random_chunk_cv_df)]

for(i in 1:nrow(aa_auprc_base_train)){
  print(i)
  for(j in 1:ncol(aa_auprc_base_train)){
    aatst <- Partition_6_random_chunk_cv_df[Partition_6_random_chunk_cv_df[, 3+j] == 1, c(1:3)]
    
    aa_auprc_base_train[i, j] <- sum((aatst$label == 1) & (aatst$owner == aafn[i]))/sum(aatst$owner == aafn[i])
  }
}
View(aa_auprc_base_train)

aa_auprc_base_train <- aa_auprc_base_train[, c(c(6:10), c(1:5))]


family_mat_list_all <- list()

for(i in 1:length(aa_all_auroc_list)){
  print(aafn[i])
  family_mat_list <- list()
  aana <- rownames(aa_all_auroc_list[[i]])
  aana_sp <- strsplit(aana, "_")
  aana_sp_sort <- lapply(aana_sp, sort)
  aana_allfam <- sort(unique(unlist(aana_sp)))
  aana_allfam <- setdiff(aana_allfam,"sequ")
  for(j in 1:length(aana_allfam)){
    aatmp1 <- which(unlist(lapply(aana_sp, function(x) aana_allfam[j] %in% x)))
    #print(aana_allfam[j])
    #print(aana[aatmp1])
    aatmp2 <- numeric(length = length(aatmp1))
    for(k in 1:length(aatmp1)){
      if(length(aana_sp[[aatmp1[k]]]) > 1){
        aanew_name <- setdiff(aana_sp_sort[[aatmp1[k]]], aana_allfam[j])
        if((length(aanew_name) == 1) & (aanew_name[1] == "sequ")){
          aanew_name <- sort(unlist(strsplit("kmerFq_motSc_Rep_pairs_triplx", "_")))
          aastm <- which(unlist(lapply(aana_sp_sort, function(x) identical(aanew_name, x))))
          if(length(aastm) == 0){
            aanew_name <- sort(unlist(strsplit("kmerFq_motSc_Rep_triplx", "_")))
            aastm <- which(unlist(lapply(aana_sp_sort, function(x) identical(aanew_name, x))))
            if(length(aastm) == 0){
              aanew_name <- sort(unlist(strsplit("kmerFq_Rep_triplx", "_")))
            }
          }
        }
        aastm <- which(unlist(lapply(aana_sp_sort, function(x) identical(aanew_name, x))))
        if(length(aastm) == 0){
          print(paste("nomatch", names(aa_all_auroc_list)[i], aana_allfam[j],"###", paste(aana_sp_sort[[aatmp1[k]]], collapse = "_"),"###", paste(aanew_name, collapse = "||")))
          aatmp2[k] <- 0
        }else{
          aatmp2[k] <- aastm 
        }
        
      }else{
        aatmp2[k] <- 0
      }
      
    }
    family_mat_list[[j]] <- cbind(aatmp1, aatmp2)
  }
  names(family_mat_list) <- aana_allfam
  family_mat_list_all[[i]] <- family_mat_list
}

names(family_mat_list_all) <- aafn
length(family_mat_list_all[[1]])

# aa_familycomp_df
# aa_all_auroc_list
# aa_all_auprc_list

aatmp_lnc_df_list <- list()
for(i in 1:length(family_mat_list_all)){
  aatmp_lnc_fam_df_list <- list()
  print(aafn[i])
  for(j in 1:length(family_mat_list_all[[i]])){#going over families
    print(names(family_mat_list_all[[i]])[j])
    aatmptroc <- matrix(nrow = nrow(family_mat_list_all[[i]][[j]]), 
                        ncol = ncol(aa_all_auroc_list[[i]]))
    aatmptprc <- matrix(nrow = nrow(family_mat_list_all[[i]][[j]]), 
                        ncol = ncol(aa_all_auprc_list[[i]]))
    for(k in 1:nrow(family_mat_list_all[[i]][[j]])){ #going over comparisons
      if(family_mat_list_all[[i]][[j]][k,2] == 0){
        aatmpcomproc <- rep(0.5, ncol(aa_all_auroc_list[[i]]))
        aatmpcompprc <- aa_auprc_base[i, ]
      }else{
        aatmpcomproc <- aa_all_auroc_list[[i]][family_mat_list_all[[i]][[j]][k,2],]
        aatmpcompprc <- aa_all_auprc_list[[i]][family_mat_list_all[[i]][[j]][k,2],]
      }
      aatmptroc[k,] <- aa_all_auroc_list[[i]][family_mat_list_all[[i]][[j]][k,1],] - aatmpcomproc
      aatmptprc[k,] <- aa_all_auprc_list[[i]][family_mat_list_all[[i]][[j]][k,1],] - aatmpcompprc
    }
    aatmptroc_num <- colSums(!is.na(aatmptroc))
    aatmptprc_num <- colSums(!is.na(aatmptprc))
    
    aatmptroc_mean <- colMeans(aatmptroc, na.rm = T)
    aatmptprc_mean <- colMeans(aatmptprc, na.rm = T)
    aatmptroc_sd <- apply(aatmptroc, MARGIN = 2, FUN = sd, na.rm = T)
    aatmptprc_sd <- apply(aatmptprc, MARGIN = 2, FUN = sd, na.rm = T)
    aatmp_lnc_fam_df_list[[j]] <- data.frame(lncRNA = rep(aafn[i], length(aatmptroc_mean)), 
                                   feature_family = rep(names(family_mat_list_all[[i]])[j], length(aatmptroc_mean)), 
                                   partition=c(paste0("CCV", c(1:5)),paste0("RCV", c(1:5))),
                                   partition_type = c(rep("chunk", 5), rep("random", 5)),
                                   nu_comparison = aatmptroc_num, 
                                   auroc_mean_delta = aatmptroc_mean,
                                   auroc_sd_delta = aatmptroc_sd,
                                   auprc_mean_delta = aatmptprc_mean,
                                   auprc_sd_delta = aatmptprc_sd)
    
  }
  aatmp_lnc_df_list[[i]] <- do.call(rbind, aatmp_lnc_fam_df_list)
  
}

aatmp_lnc_df <- do.call(rbind, aatmp_lnc_df_list)
View(aatmp_lnc_df)



aatttgc <- summarySE(aatmp_lnc_df, measurevar=c("auroc_mean_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc_1 <-  aatttgc[aatttgc$partition_type == "random",]
aatttgc_2 <-  aatttgc[aatttgc$partition_type == "chunk",]

#aatgc <- aatgc[aatgc$Gene_number %in% c(100, 500, 1000),]
ggplot(aatttgc_1, aes(x=lncRNA, y=auroc_mean_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=auroc_mean_delta-se, ymax=auroc_mean_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auroc random")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

ggplot(aatttgc_2, aes(x=lncRNA, y=auroc_mean_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=auroc_mean_delta-se, ymax=auroc_mean_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auroc chunk")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))



aatttgc2 <- summarySE(aatmp_lnc_df, measurevar=c("auprc_mean_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc2_1 <-  aatttgc2[aatttgc2$partition_type == "random",]
aatttgc2_2 <-  aatttgc2[aatttgc2$partition_type == "chunk",]

ggplot(aatttgc2_1, aes(x=lncRNA, y=auprc_mean_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=auprc_mean_delta-se, ymax=auprc_mean_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auprc random")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

ggplot(aatttgc2_2, aes(x=lncRNA, y=auprc_mean_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=auprc_mean_delta-se, ymax=auprc_mean_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auprc chunk")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

# create for matrices delta AUROC/PRC in Chunk/Random
# each row a lncRNA
# each column a Feature family
aatttgc <- summarySE(aatmp_lnc_df, measurevar=c("auroc_mean_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc_1 <-  aatttgc[aatttgc$partition_type == "random",]
aatttgc_2 <-  aatttgc[aatttgc$partition_type == "chunk",]

aatttgc2 <- summarySE(aatmp_lnc_df, measurevar=c("auprc_mean_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc2_1 <-  aatttgc2[aatttgc2$partition_type == "random",]
aatttgc2_2 <-  aatttgc2[aatttgc2$partition_type == "chunk",]

aa_auprc_chunk <- matrix(nrow = length(unique(aatmp_lnc_df$lncRNA)), 
                         ncol = length(unique(aatmp_lnc_df$feature_family)))

rownames(aa_auprc_chunk) <- unique(aatmp_lnc_df$lncRNA)
colnames(aa_auprc_chunk) <- unique(aatmp_lnc_df$feature_family)

aa_auprc_rand <- matrix(nrow = length(unique(aatmp_lnc_df$lncRNA)), 
                         ncol = length(unique(aatmp_lnc_df$feature_family)))

rownames(aa_auprc_rand) <- unique(aatmp_lnc_df$lncRNA)
colnames(aa_auprc_rand) <- unique(aatmp_lnc_df$feature_family)

aa_auroc_chunk <- matrix(nrow = length(unique(aatmp_lnc_df$lncRNA)), 
                         ncol = length(unique(aatmp_lnc_df$feature_family)))

rownames(aa_auroc_chunk) <- unique(aatmp_lnc_df$lncRNA)
colnames(aa_auroc_chunk) <- unique(aatmp_lnc_df$feature_family)

aa_auroc_rand <- matrix(nrow = length(unique(aatmp_lnc_df$lncRNA)), 
                         ncol = length(unique(aatmp_lnc_df$feature_family)))

rownames(aa_auroc_rand) <- unique(aatmp_lnc_df$lncRNA)
colnames(aa_auroc_rand) <- unique(aatmp_lnc_df$feature_family)

for(i in 1:nrow(aatttgc_1)){
  aa_auroc_rand[match(aatttgc_1$lncRNA[i], rownames(aa_auroc_rand)), 
                match(aatttgc_1$feature_family[i], colnames(aa_auroc_rand))] <- aatttgc_1$auroc_mean_delta[i]
}
for(i in 1:nrow(aatttgc_2)){
  aa_auroc_chunk[match(aatttgc_2$lncRNA[i], rownames(aa_auroc_chunk)), 
                match(aatttgc_2$feature_family[i], colnames(aa_auroc_chunk))] <- aatttgc_2$auroc_mean_delta[i]
}
for(i in 1:nrow(aatttgc2_1)){
  aa_auprc_rand[match(aatttgc2_1$lncRNA[i], rownames(aa_auprc_rand)), 
                match(aatttgc2_1$feature_family[i], colnames(aa_auprc_rand))] <- aatttgc2_1$auprc_mean_delta[i]
}
for(i in 1:nrow(aatttgc2_2)){
  aa_auprc_chunk[match(aatttgc2_2$lncRNA[i], rownames(aa_auprc_chunk)), 
                match(aatttgc2_2$feature_family[i], colnames(aa_auprc_chunk))] <- aatttgc2_2$auprc_mean_delta[i]
}

library(gplots)
#heatmap.2(aa_auroc_rand)
require(RColorBrewer)


aa_auroc_rand[is.na(aa_auroc_rand)] <- 0
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(19)
# z <- zClust(x=aa_auroc_rand, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
heatmap.2(aa_auroc_rand, trace='none', col=(cols), breaks = seq(-0.04, 0.34, 0.02), margins = c(6,10))

aa_auroc_chunk[is.na(aa_auroc_chunk)] <- 0
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(17)
# z <- zClust(x=aa_auroc_chunk, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
heatmap.2(aa_auroc_chunk, trace='none', col=(cols), breaks = seq(-0.04, 0.30, 0.02), margins = c(6,10))

aa_auprc_rand[is.na(aa_auprc_rand)] <- 0
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(22)
# z <- zClust(x=aa_auprc_rand, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
heatmap.2(aa_auprc_rand, trace='none', col=(cols), breaks = seq(-0.06, 0.38, 0.02), margins = c(6,10))

aa_auprc_chunk[is.na(aa_auprc_chunk)] <- 0
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(19)
# z <- zClust(x=aa_auprc_chunk, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
heatmap.2(aa_auprc_chunk, trace='none', col=(cols), breaks = seq(-0.04, 0.34, 0.02), margins = c(6,10))





################################################################################################
# Redoing the above analysis excluding distance
################################################################################################


family_mat_list_all_noDist <- list()

aa_all_auroc_list_noDist <- aa_all_auroc_list
aa_all_auprc_list_noDist <- aa_all_auprc_list
for(i in 1:length(aa_all_auroc_list_noDist)){
  print(names(aa_all_auroc_list)[i])
  aagr <- grep(pattern = "dist", ignore.case = T, x = rownames(aa_all_auroc_list_noDist[[i]]))
  print(length(aagr))
  if(length(aagr) > 0){
    aa_all_auroc_list_noDist[[i]] <- aa_all_auroc_list_noDist[[i]][-aagr,]
    aa_all_auprc_list_noDist[[i]] <- aa_all_auprc_list_noDist[[i]][-aagr,]
  }
}

for(i in 1:length(aa_all_auroc_list_noDist)){
  print(aafn[i])
  family_mat_list <- list()
  aana <- rownames(aa_all_auroc_list_noDist[[i]])
  aana_sp <- strsplit(aana, "_")
  aana_sp_sort <- lapply(aana_sp, sort)
  aana_allfam <- sort(unique(unlist(aana_sp)))
  aana_allfam <- setdiff(aana_allfam,"sequ")
  for(j in 1:length(aana_allfam)){
    aatmp1 <- which(unlist(lapply(aana_sp, function(x) aana_allfam[j] %in% x)))
    #print(aana_allfam[j])
    #print(aana[aatmp1])
    aatmp2 <- numeric(length = length(aatmp1))
    for(k in 1:length(aatmp1)){
      if(length(aana_sp[[aatmp1[k]]]) > 1){
        aanew_name <- setdiff(aana_sp_sort[[aatmp1[k]]], aana_allfam[j])
        if((length(aanew_name) == 1) & (aanew_name[1] == "sequ")){
          aanew_name <- sort(unlist(strsplit("kmerFq_motSc_Rep_pairs_triplx", "_")))
          aastm <- which(unlist(lapply(aana_sp_sort, function(x) identical(aanew_name, x))))
          if(length(aastm) == 0){
            aanew_name <- sort(unlist(strsplit("kmerFq_motSc_Rep_triplx", "_")))
            aastm <- which(unlist(lapply(aana_sp_sort, function(x) identical(aanew_name, x))))
            if(length(aastm) == 0){
              aanew_name <- sort(unlist(strsplit("kmerFq_Rep_triplx", "_")))
            }
          }
        }
        aastm <- which(unlist(lapply(aana_sp_sort, function(x) identical(aanew_name, x))))
        if(length(aastm) == 0){
          print(paste("nomatch", names(aa_all_auroc_list_noDist)[i], aana_allfam[j],"###", paste(aana_sp_sort[[aatmp1[k]]], collapse = "_"),"###", paste(aanew_name, collapse = "||")))
          aatmp2[k] <- 0
        }else{
          aatmp2[k] <- aastm 
        }
        
      }else{
        aatmp2[k] <- 0
      }
      
    }
    family_mat_list[[j]] <- cbind(aatmp1, aatmp2)
  }
  names(family_mat_list) <- aana_allfam
  family_mat_list_all_noDist[[i]] <- family_mat_list
}

names(family_mat_list_all_noDist) <- aafn
length(family_mat_list_all_noDist[[1]])

# aa_familycomp_df
# aa_all_auroc_list
# aa_all_auprc_list

aatmp_lnc_df_list_noDist <- list()
for(i in 1:length(family_mat_list_all_noDist)){
  aatmp_lnc_fam_df_list <- list()
  print(aafn[i])
  for(j in 1:length(family_mat_list_all_noDist[[i]])){#going over families
    print(names(family_mat_list_all_noDist[[i]])[j])
    aatmptroc <- matrix(nrow = nrow(family_mat_list_all_noDist[[i]][[j]]), 
                        ncol = ncol(aa_all_auroc_list_noDist[[i]]))
    aatmptprc <- matrix(nrow = nrow(family_mat_list_all_noDist[[i]][[j]]), 
                        ncol = ncol(aa_all_auprc_list_noDist[[i]]))
    for(k in 1:nrow(family_mat_list_all_noDist[[i]][[j]])){ #going over comparisons
      if(family_mat_list_all_noDist[[i]][[j]][k,2] == 0){
        aatmpcomproc <- rep(0.5, ncol(aa_all_auroc_list_noDist[[i]]))
        aatmpcompprc <- aa_auprc_base[i, ]
      }else{
        aatmpcomproc <- aa_all_auroc_list_noDist[[i]][family_mat_list_all_noDist[[i]][[j]][k,2],]
        aatmpcompprc <- aa_all_auprc_list_noDist[[i]][family_mat_list_all_noDist[[i]][[j]][k,2],]
      }
      aatmptroc[k,] <- aa_all_auroc_list_noDist[[i]][family_mat_list_all_noDist[[i]][[j]][k,1],] - aatmpcomproc
      aatmptprc[k,] <- aa_all_auprc_list_noDist[[i]][family_mat_list_all_noDist[[i]][[j]][k,1],] - aatmpcompprc
    }
    aatmptroc_num <- colSums(!is.na(aatmptroc))
    aatmptprc_num <- colSums(!is.na(aatmptprc))
    
    aatmptroc_mean <- colMeans(aatmptroc, na.rm = T)
    aatmptprc_mean <- colMeans(aatmptprc, na.rm = T)
    aatmptroc_sd <- apply(aatmptroc, MARGIN = 2, FUN = sd, na.rm = T)
    aatmptprc_sd <- apply(aatmptprc, MARGIN = 2, FUN = sd, na.rm = T)
    aatmp_lnc_fam_df_list[[j]] <- data.frame(lncRNA = rep(aafn[i], length(aatmptroc_mean)), 
                                             feature_family = rep(names(family_mat_list_all_noDist[[i]])[j], length(aatmptroc_mean)), 
                                             partition=c(paste0("CCV", c(1:5)),paste0("RCV", c(1:5))),
                                             partition_type = c(rep("chunk", 5), rep("random", 5)),
                                             nu_comparison = aatmptroc_num, 
                                             auroc_mean_delta = aatmptroc_mean,
                                             auroc_sd_delta = aatmptroc_sd,
                                             auprc_mean_delta = aatmptprc_mean,
                                             auprc_sd_delta = aatmptprc_sd)
    
  }
  aatmp_lnc_df_list_noDist[[i]] <- do.call(rbind, aatmp_lnc_fam_df_list)
  
}

aatmp_lnc_df_noDist <- do.call(rbind, aatmp_lnc_df_list_noDist)
View(aatmp_lnc_df_noDist)

# create for matrices delta AUROC/PRC in Chunk/Random -- after removing distance models
# each row a lncRNA
# each column a Feature family
aatttgc_nodist <- summarySE(aatmp_lnc_df_noDist, measurevar=c("auroc_mean_delta"), 
                            groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc_1_nodist <-  aatttgc_nodist[aatttgc_nodist$partition_type == "random",]
aatttgc_2_nodist <-  aatttgc_nodist[aatttgc_nodist$partition_type == "chunk",]

aatttgc2_nodist <- summarySE(aatmp_lnc_df_noDist, measurevar=c("auprc_mean_delta"), 
                             groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc2_1_nodist <-  aatttgc2_nodist[aatttgc2_nodist$partition_type == "random",]
aatttgc2_2_nodist <-  aatttgc2_nodist[aatttgc2_nodist$partition_type == "chunk",]

aa_auprc_chunk_nodist <- matrix(nrow = length(unique(aatmp_lnc_df_noDist$lncRNA)), 
                                ncol = length(unique(aatmp_lnc_df_noDist$feature_family)))

rownames(aa_auprc_chunk_nodist) <- unique(aatmp_lnc_df_noDist$lncRNA)
colnames(aa_auprc_chunk_nodist) <- unique(aatmp_lnc_df_noDist$feature_family)

aa_auprc_rand_nodist <- matrix(nrow = length(unique(aatmp_lnc_df_noDist$lncRNA)), 
                               ncol = length(unique(aatmp_lnc_df_noDist$feature_family)))

rownames(aa_auprc_rand_nodist) <- unique(aatmp_lnc_df_noDist$lncRNA)
colnames(aa_auprc_rand_nodist) <- unique(aatmp_lnc_df_noDist$feature_family)

aa_auroc_chunk_nodist <- matrix(nrow = length(unique(aatmp_lnc_df_noDist$lncRNA)), 
                                ncol = length(unique(aatmp_lnc_df_noDist$feature_family)))

rownames(aa_auroc_chunk_nodist) <- unique(aatmp_lnc_df_noDist$lncRNA)
colnames(aa_auroc_chunk_nodist) <- unique(aatmp_lnc_df_noDist$feature_family)

aa_auroc_rand_nodist <- matrix(nrow = length(unique(aatmp_lnc_df_noDist$lncRNA)), 
                               ncol = length(unique(aatmp_lnc_df_noDist$feature_family)))

rownames(aa_auroc_rand_nodist) <- unique(aatmp_lnc_df_noDist$lncRNA)
colnames(aa_auroc_rand_nodist) <- unique(aatmp_lnc_df_noDist$feature_family)

for(i in 1:nrow(aatttgc_1_nodist)){
  aa_auroc_rand_nodist[match(aatttgc_1_nodist$lncRNA[i], rownames(aa_auroc_rand_nodist)), 
                       match(aatttgc_1_nodist$feature_family[i], colnames(aa_auroc_rand_nodist))] <- aatttgc_1_nodist$auroc_mean_delta[i]
}
for(i in 1:nrow(aatttgc_2_nodist)){
  aa_auroc_chunk_nodist[match(aatttgc_2_nodist$lncRNA[i], rownames(aa_auroc_chunk_nodist)), 
                        match(aatttgc_2_nodist$feature_family[i], colnames(aa_auroc_chunk_nodist))] <- aatttgc_2_nodist$auroc_mean_delta[i]
}
for(i in 1:nrow(aatttgc2_1)){
  aa_auprc_rand_nodist[match(aatttgc2_1_nodist$lncRNA[i], rownames(aa_auprc_rand_nodist)), 
                       match(aatttgc2_1_nodist$feature_family[i], colnames(aa_auprc_rand_nodist))] <- aatttgc2_1_nodist$auprc_mean_delta[i]
}
for(i in 1:nrow(aatttgc2_2_nodist)){
  aa_auprc_chunk_nodist[match(aatttgc2_2_nodist$lncRNA[i], rownames(aa_auprc_chunk_nodist)), 
                        match(aatttgc2_2_nodist$feature_family[i], colnames(aa_auprc_chunk_nodist))] <- aatttgc2_2_nodist$auprc_mean_delta[i]
}

library(gplots)
#heatmap.2(aa_auroc_rand)
require(RColorBrewer)


aa_auroc_rand_nodist[is.na(aa_auroc_rand_nodist)] <- 0
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(12)
# z <- zClust(x=aa_auroc_rand, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
aam1 <- heatmap.2(aa_auroc_rand_nodist, trace='none', col=(cols), breaks = seq(-0.00, 0.12, 0.01), margins = c(6,10))
aam1
aa_auroc_chunk_nodist[is.na(aa_auroc_chunk_nodist)] <- 0
range(aa_auroc_chunk_nodist)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(12)
# z <- zClust(x=aa_auroc_chunk, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
heatmap.2(aa_auroc_chunk_nodist, trace='none', col=(cols), breaks = seq(-0.00, 0.12, 0.01), margins = c(6,10)
          #, Rowv = aam1$rowInd
          , Colv = rev(aam1$colInd)
          , revC = T
          )

aa_auprc_rand_nodist[is.na(aa_auprc_rand_nodist)] <- 0
range(aa_auprc_rand_nodist)

cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(15)
# z <- zClust(x=aa_auprc_rand, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
aam <- heatmap.2(aa_auprc_rand_nodist, trace='none', col=(cols), breaks = seq(-0.00, 0.15, 0.01), margins = c(6,10))
aam
aa_auprc_chunk_nodist[is.na(aa_auprc_chunk_nodist)] <- 0
range(aa_auprc_chunk_nodist)

cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(15)
# z <- zClust(x=aa_auprc_chunk, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
heatmap.2(aa_auprc_chunk_nodist, trace='none', col=(cols), breaks = seq(-0.00, 0.15, 0.01), margins = c(6,10), Rowv = aam$rowInd)




################################################################################################
################################################################################################
# read combined performance for 28 lncRNAs (CVfirst)
aafiles <- list.files("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_Results/Combined_res/")
aafn <- unlist(lapply(strsplit(aafiles, "__"), "[[", 1))
aa_all_auroc_list_cvfirst <- list()
aa_all_auprc_list_cvfirst <- list()
aa_all_impor_list_cvfirst <- list()

for(i in 1:length(aafiles)){
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_Results/Combined_res/", aafiles[i]))
  aa_all_auroc_list_cvfirst[[i]] <- my_perf[[1]][[1]]
  aa_all_auprc_list_cvfirst[[i]] <- my_perf[[2]][[1]]
  aa_all_impor_list_cvfirst[[i]] <- my_perf[[3]][[1]]
}
names(aa_all_auroc_list_cvfirst) <- aafn
names(aa_all_auprc_list_cvfirst) <- aafn
names(aa_all_impor_list_cvfirst) <- aafn

aac <- max(unlist(lapply(aa_all_auprc_list_cvfirst, nrow)))

aa_all_auroc_rand_cvfirst <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auroc_rand_cvfirst) <- rownames(aa_all_auprc_list_cvfirst[[2]])
rownames(aa_all_auroc_rand_cvfirst) <- aafn

aa_all_auroc_chunk_cvfirst <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auroc_chunk_cvfirst) <- rownames(aa_all_auprc_list_cvfirst[[2]])
rownames(aa_all_auroc_chunk_cvfirst) <- aafn

aa_all_auprc_rand_cvfirst <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auprc_rand_cvfirst) <- rownames(aa_all_auprc_list_cvfirst[[2]])
rownames(aa_all_auprc_rand_cvfirst) <- aafn

aa_all_auprc_chunk_cvfirst <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auprc_chunk_cvfirst) <- rownames(aa_all_auprc_list_cvfirst[[2]])
rownames(aa_all_auprc_chunk_cvfirst) <- aafn


aa_all_auprc_list_cvfirst_diff <- list()
for(i in 1:length(aa_all_auprc_list_cvfirst)){
  aa_all_auprc_list_cvfirst_diff[[i]] <- aa_all_auprc_list_cvfirst[[i]]
  for(j in 1:nrow(aa_all_auprc_list_cvfirst[[i]])){
    aa_all_auprc_list_cvfirst_diff[[i]][j,] <- aa_all_auprc_list_cvfirst_diff[[i]][j,] - aa_auprc_base[i,]
  }
}
for(i in 1:length(aa_all_auroc_list_cvfirst)){
  aatmmpccv_roc <- rowMeans(aa_all_auroc_list_cvfirst[[i]][,c(1:5)], na.rm = T)
  aatmmprcv_roc <- rowMeans(aa_all_auroc_list_cvfirst[[i]][,c(6:10)], na.rm = T)
  
  aatmmpccv_prc <- rowMeans(aa_all_auprc_list_cvfirst_diff[[i]][,c(1:5)], na.rm = T)
  aatmmprcv_prc <- rowMeans(aa_all_auprc_list_cvfirst_diff[[i]][,c(6:10)], na.rm = T)
  
  aacolord <- match(rownames(aa_all_auroc_list_cvfirst[[i]]), 
                    colnames(aa_all_auroc_chunk_cvfirst))
  
  aa_all_auroc_chunk_cvfirst[i,aacolord] <- aatmmpccv_roc
  aa_all_auroc_rand_cvfirst[i,aacolord] <- aatmmprcv_roc
  aa_all_auprc_chunk_cvfirst[i,aacolord] <- aatmmpccv_prc
  aa_all_auprc_rand_cvfirst[i,aacolord] <- aatmmprcv_prc
  
  
}



#aa_all_auroc_rand_cvfirst[is.na(aa_all_auroc_rand_cvfirst)] <- 0
# cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(22)
# z <- zClust(x=t(aa_all_auroc_rand_cvfirst), scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(0.5, 1, 0.02), margins = c(10,21))
png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_Results/AUROC_rand_heatmap.png", 
    width = 10*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
range(aa_all_auroc_rand_cvfirst)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(25)
heatmap.2(t(aa_all_auroc_rand_cvfirst), trace='none', col=(cols), breaks = seq(0.5, 1, 0.02), margins = c(10,21))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_Results/AUROC_chunk_heatmap.png", 
    width = 10*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#sum(is.na(aa_all_auroc_chunk_cvfirst))
#aa_all_auroc_chunk_cvfirst[is.na(aa_all_auroc_chunk_cvfirst)] <- 0
range(aa_all_auroc_chunk_cvfirst)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(25)
heatmap.2(t(aa_all_auroc_chunk_cvfirst), trace='none', col=(cols), breaks = seq(0.5, 1, 0.02), margins = c(10,21))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_Results/AUPRC_rand_heatmap.png", 
    width = 10*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#sum(is.na(aa_all_auprc_rand_cvfirst))
#aa_all_auprc_rand_cvfirst[is.na(aa_all_auprc_rand_cvfirst)] <- 0
range(aa_all_auprc_rand_cvfirst)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(35)
heatmap.2(t(aa_all_auprc_rand_cvfirst), trace='none', col=(cols), breaks = seq(0.0, 0.7, 0.02), margins = c(10,21))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_Results/AUPRC_chunk_heatmap.png", 
    width = 10*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#sum(is.na(aa_all_auprc_chunk_cvfirst))
#aa_all_auprc_chunk_cvfirst[is.na(aa_all_auprc_chunk_cvfirst)] <- 0
range(aa_all_auprc_chunk_cvfirst)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(35)
heatmap.2(t(aa_all_auprc_chunk_cvfirst), trace='none', col=(cols), breaks = seq(0.0, 0.7, 0.02), margins = c(10,21))
dev.off()

# distance_vs_full models
aagrep <- grep(pattern = "dist", ignore.case = T, x = colnames(aa_all_auroc_rand_cvfirst))
aagr2 <- setdiff(colnames(aa_all_auroc_rand_cvfirst)[aagrep], "RBPTFRepeatTriplexChromAccessMethChIPDinucKmerExpDist")

#plot()
aaCompMat <- cbind(rep("RBPTFRepeatTriplexChromAccessMethChIPDinucKmerExpDist", length(aagr2)), aagr2)

library(RColorBrewer)
#n <- 28
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#col_vector_Samp <- sample(col_vector, n)
#pie(rep(1,n), col=col_vector_Samp)
#names(col_vector_Samp) <- rownames(performance_auroc_chunk_28_34)

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_Results/comparison_dist_auroc_rand.png", 
    width = 10*300,
    height = 5*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(2,4), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = aa_all_auroc_rand_cvfirst[, match(aaCompMat[i,1], colnames(aa_all_auroc_rand_cvfirst))],
       y = aa_all_auroc_rand_cvfirst[, match(aaCompMat[i,2], colnames(aa_all_auroc_rand_cvfirst))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.4,1), ylim = c(0.4,1),
       cex.lab = 0.7,
       cex= 1.2,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0.0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.3, y.intersp = 0.7) 
  
}
dev.off()


png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_Results/comparison_dist_auroc_chunk.png", 
    width = 10*300,
    height = 5*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(2,4), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = aa_all_auroc_chunk_cvfirst[, match(aaCompMat[i,1], colnames(aa_all_auroc_chunk_cvfirst))],
       y = aa_all_auroc_chunk_cvfirst[, match(aaCompMat[i,2], colnames(aa_all_auroc_chunk_cvfirst))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.4,1), ylim = c(0.4,1),
       cex.lab = 0.7,
       cex= 1.2,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0.0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.3, y.intersp = 0.7) 
  
}
dev.off()


png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_Results/comparison_dist_auprc_rand.png", 
    width = 10*300,
    height = 5*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(2,4), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = aa_all_auprc_rand_cvfirst[, match(aaCompMat[i,1], colnames(aa_all_auprc_rand_cvfirst))],
       y = aa_all_auprc_rand_cvfirst[, match(aaCompMat[i,2], colnames(aa_all_auprc_rand_cvfirst))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.0,0.71), ylim = c(0.0,0.71),
       cex.lab = 0.7,
       cex= 1.2,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0.0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.3, y.intersp = 0.7) 
  
}
dev.off()


png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_Results/comparison_dist_auprc_chunk.png", 
    width = 10*300,
    height = 5*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(2,4), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = aa_all_auprc_chunk_cvfirst[, match(aaCompMat[i,1], colnames(aa_all_auprc_chunk_cvfirst))],
       y = aa_all_auprc_chunk_cvfirst[, match(aaCompMat[i,2], colnames(aa_all_auprc_chunk_cvfirst))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2],  xlim = c(0.0,0.71), ylim = c(0.0,0.71),
       cex.lab = 0.7,
       cex= 1.2,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0.0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.3, y.intersp = 0.7) 
  
}
dev.off()

################################################################################################################
colnames(aa_all_auprc_chunk_cvfirst)[grep("TX", colnames(aa_all_auprc_chunk_cvfirst))]


aaCompMat <- rbind(c("RBPscanned" , "RBPscannedTX"),
                   c("TFscanned" , "TFscannedTX"),
                   c("Repeat" , "RepeatTX"),
                   c("RBPTFRepeat" , "RBPTFRepeatTXRBPrep"),
                   c("RBPTFRepeatTriplex" , "RBPTFRepeatTriplexTXRBPrep"),
                   c("RBPTFRepeatTriplexChromAccessMeth" , "RBPTFRepeatTriplexChromAccessMethTX"),
                   c("RBPTFRepeatTriplexChromAccessMethChIPDinuc" , "RBPTFRepeatTriplexChromAccessMethChIPDinucTX"),
                   c("RBPTFRepeatTriplexChromAccessMethChIPkmer" , "RBPTFRepeatTriplexChromAccessMethChIPkmerTX"),
                   c("RBPTFRepeatTriplexChromAccessMethChIPDinucKmerExp" , "RBPTFRepeatTriplexChromAccessMethChIPDinucKmerExpTX"),
                   c("RBPTFRepeatTriplexChromAccessMethChIPDinucKmerExpDist" , "RBPTFRepeatTriplexChromAccessMethChIPDinucKmerExpDistTX"))

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_Results/comparison_TX_auroc_rand.png", 
    width = 10*300,
    height = 5*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(2,5), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = aa_all_auroc_rand_cvfirst[, match(aaCompMat[i,1], colnames(aa_all_auroc_rand_cvfirst))],
       y = aa_all_auroc_rand_cvfirst[, match(aaCompMat[i,2], colnames(aa_all_auroc_rand_cvfirst))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.4,1), ylim = c(0.4,1),
       cex.lab = 0.7,
       cex= 0.8,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0.0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.3, y.intersp = 0.7) 
  
}
dev.off()


png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_Results/comparison_TX_auroc_chunk.png", 
    width = 10*300,
    height = 5*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(2,5), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = aa_all_auroc_chunk_cvfirst[, match(aaCompMat[i,1], colnames(aa_all_auroc_chunk_cvfirst))],
       y = aa_all_auroc_chunk_cvfirst[, match(aaCompMat[i,2], colnames(aa_all_auroc_chunk_cvfirst))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.4,1), ylim = c(0.4,1),
       cex.lab = 0.7,
       cex= 0.8,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0.0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.3, y.intersp = 0.7) 
  
}
dev.off()


png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_Results/comparison_TX_auprc_rand.png", 
    width = 10*300,
    height = 5*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(2,5), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = aa_all_auprc_rand_cvfirst[, match(aaCompMat[i,1], colnames(aa_all_auprc_rand_cvfirst))],
       y = aa_all_auprc_rand_cvfirst[, match(aaCompMat[i,2], colnames(aa_all_auprc_rand_cvfirst))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.0,0.5), ylim = c(0.0,0.5),
       cex.lab = 0.7,
       cex= 0.8,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0.0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.3, y.intersp = 0.7) 
  
}
dev.off()


png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_Results/comparison_TX_auprc_chunk.png", 
    width = 10*300,
    height = 5*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(2,5), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = aa_all_auprc_chunk_cvfirst[, match(aaCompMat[i,1], colnames(aa_all_auprc_chunk_cvfirst))],
       y = aa_all_auprc_chunk_cvfirst[, match(aaCompMat[i,2], colnames(aa_all_auprc_chunk_cvfirst))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2],  xlim = c(0.0,0.5), ylim = c(0.0,0.5),
       cex.lab = 0.7,
       cex= 0.8,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0.0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.3, y.intersp = 0.7) 
  
}
dev.off()

############################################################################################################################################################
colnames(aa_all_auprc_chunk_cvfirst)
aaCompMat <- rbind(c("Kmer" , "DinucFreq"),
                   c("Kmer" , "DistanceOnly"),
                   c("Kmer" , "RBPscanned"),
                   c("Kmer" , "TFscanned"),
                   c("Kmer" , "Repeat"),
                   c("Kmer" , "Triplex"),
                   c("Kmer" , "RBPTFRepeatTriplex"),
                   c("Kmer" , "AccessMethChromChIP"),
                   c("Kmer" , "Transcription"),
                   c("Kmer" , "Pairs"),
                   c("Kmer" , "RBPTFRepeatTriplexChromAccessMethChIPExp"),
                   c("RBPTFRepeat" , "AccessMethChromChIP"),
                   c("RBPscanned" , "TFscanned"),
                   c("RBPscanned" , "Repeat"),
                   c("TFscanned" , "Repeat"),
                   c("Accessibility", "Methylation"),
                   c("AccessMeth", "Chromatin"),
                   c("AccessMethChromChIP", "Transcription")
)

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_Results/comparison_FMFM_auroc_rand.png", 
    width = 10*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(4,5), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = aa_all_auroc_rand_cvfirst[, match(aaCompMat[i,1], colnames(aa_all_auroc_rand_cvfirst))],
       y = aa_all_auroc_rand_cvfirst[, match(aaCompMat[i,2], colnames(aa_all_auroc_rand_cvfirst))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.4,1), ylim = c(0.4,1),
       cex.lab = 0.7,
       cex= 0.8,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0.0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.3, y.intersp = 0.7) 
  
}
dev.off()


png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_Results/comparison_FMFM_auroc_chunk.png", 
    width = 10*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(4,5), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = aa_all_auroc_chunk_cvfirst[, match(aaCompMat[i,1], colnames(aa_all_auroc_chunk_cvfirst))],
       y = aa_all_auroc_chunk_cvfirst[, match(aaCompMat[i,2], colnames(aa_all_auroc_chunk_cvfirst))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.4,1), ylim = c(0.4,1),
       cex.lab = 0.7,
       cex= 0.8,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0.0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.3, y.intersp = 0.7) 
  
}
dev.off()


png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_Results/comparison_FMFM_auprc_rand.png", 
    width = 10*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(4,5), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = aa_all_auprc_rand_cvfirst[, match(aaCompMat[i,1], colnames(aa_all_auprc_rand_cvfirst))],
       y = aa_all_auprc_rand_cvfirst[, match(aaCompMat[i,2], colnames(aa_all_auprc_rand_cvfirst))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2], xlim = c(0.0,0.5), ylim = c(0.0,0.5),
       cex.lab = 0.7,
       cex= 0.8,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0.0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.3, y.intersp = 0.7) 
  
}
dev.off()


png(filename = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_Results/comparison_FMFM_auprc_chunk.png", 
    width = 10*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#par(mfrow = c(4,4), mar = c(4,4,4,10))
par(mfrow = c(4,5), mar = c(4,8,8,4),xpd = F)
for(i in 1:nrow(aaCompMat)){
  plot(x = aa_all_auprc_chunk_cvfirst[, match(aaCompMat[i,1], colnames(aa_all_auprc_chunk_cvfirst))],
       y = aa_all_auprc_chunk_cvfirst[, match(aaCompMat[i,2], colnames(aa_all_auprc_chunk_cvfirst))], 
       xlab = aaCompMat[i,1], ylab = aaCompMat[i,2],  xlim = c(0.0,0.5), ylim = c(0.0,0.5),
       cex.lab = 0.7,
       cex= 0.8,
       pch = 16, 
       col = col_vector_Samp)
  abline(a = 0, b =1, col = 1, lty = 2)
  
  legend(x = "top", inset=c(0.0,-0.7),xpd = T,ncol=3,
         legend = names(col_vector_Samp), 
         fill = col_vector_Samp,
         bty = "n",xjust = 0.5, yjust = 0,
         cex = 0.7, x.intersp = 0.3, y.intersp = 0.7) 
  
}
dev.off()

################################################################################################################
################################################################################################################
################################################################################################
################################################################################################
# write jobs to gather performance eval




aa_un_own <- sort(unique(Partition_6_random_chunk_cv_df$owner))


for(i in 1:length(aa_un_own)){
  cat(c("Rscript --vanilla /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/gather_perf.R",
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_svd2_lowres_Results/Learned_models/",aa_un_own[i] ),
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_svd2_lowres_Results/Combined_results/",aa_un_own[i], "__perfMat.RData\n" )),
      sep = " ", 
      file = "~/Documents/Shayan/BioInf/lncRNA/gather_perf_all_svd2.job",
      append = !(i == 1))
}
################################################################################################.
################################################################################################
# read combined performance for 28 lncRNAs (family svd2) --> distance lower resolution
aafiles <- list.files("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res2/")
aafn <- unlist(lapply(strsplit(aafiles, "__"), "[[", 1))
aa_all_auroc_list2 <- list()
aa_all_auprc_list2 <- list()
aa_all_impor_list2 <- list()

for(i in 1:length(aafiles)){
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res2/", aafiles[i]))
  aa_all_auroc_list2[[i]] <- my_perf[[1]][[1]]
  aa_all_auprc_list2[[i]] <- my_perf[[2]][[1]]
  aa_all_impor_list2[[i]] <- my_perf[[3]][[1]]
}
names(aa_all_auroc_list2) <- aafn
names(aa_all_auprc_list2) <- aafn
names(aa_all_impor_list2) <- aafn

View(aa_all_auroc_list2$Malat1)

aa_auprc_base <- matrix(nrow = length(aafn), ncol = (ncol(Partition_6_random_chunk_cv_df) - 3))
rownames(aa_auprc_base) <- aafn
colnames(aa_auprc_base) <- colnames(Partition_6_random_chunk_cv_df)[4:ncol(Partition_6_random_chunk_cv_df)]

for(i in 1:nrow(aa_auprc_base)){
  print(i)
  for(j in 1:ncol(aa_auprc_base)){
    aatst <- Partition_6_random_chunk_cv_df[Partition_6_random_chunk_cv_df[, 3+j] == 0, c(1:3)]
    
    aa_auprc_base[i, j] <- sum((aatst$label == 1) & (aatst$owner == aafn[i]))/sum(aatst$owner == aafn[i])
  }
}
View(aa_auprc_base)

aa_auprc_base <- aa_auprc_base[, c(c(6:10), c(1:5))]







aac <- max(unlist(lapply(aa_all_auprc_list2, nrow)))

aa_all_auroc_rand_2 <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auroc_rand_2) <- rownames(aa_all_auprc_list2[[2]])
rownames(aa_all_auroc_rand_2) <- aafn

aa_all_auroc_chunk_2 <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auroc_chunk_2) <- rownames(aa_all_auprc_list2[[2]])
rownames(aa_all_auroc_chunk_2) <- aafn

aa_all_auprc_rand_2 <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auprc_rand_2) <- rownames(aa_all_auprc_list2[[2]])
rownames(aa_all_auprc_rand_2) <- aafn

aa_all_auprc_chunk_2 <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auprc_chunk_2) <- rownames(aa_all_auprc_list2[[2]])
rownames(aa_all_auprc_chunk_2) <- aafn


aa_all_auprc_list2_diff <- list()
for(i in 1:length(aa_all_auprc_list2)){
  aa_all_auprc_list2_diff[[i]] <- aa_all_auprc_list2[[i]]
  for(j in 1:nrow(aa_all_auprc_list2[[i]])){
    aa_all_auprc_list2_diff[[i]][j,] <- aa_all_auprc_list2_diff[[i]][j,] - aa_auprc_base[i,]
  }
}
for(i in 1:length(aa_all_auroc_list2)){
  aatmmpccv_roc <- rowMeans(aa_all_auroc_list2[[i]][,c(1:5)], na.rm = T)
  aatmmprcv_roc <- rowMeans(aa_all_auroc_list2[[i]][,c(6:10)], na.rm = T)
  
  aatmmpccv_prc <- rowMeans(aa_all_auprc_list2_diff[[i]][,c(1:5)], na.rm = T)
  aatmmprcv_prc <- rowMeans(aa_all_auprc_list2_diff[[i]][,c(6:10)], na.rm = T)
  
  aacolord <- match(rownames(aa_all_auroc_list2[[i]]), 
                    colnames(aa_all_auroc_chunk_2))
  
  aa_all_auroc_chunk_2[i,aacolord] <- aatmmpccv_roc
  aa_all_auroc_rand_2[i,aacolord] <- aatmmprcv_roc
  aa_all_auprc_chunk_2[i,aacolord] <- aatmmpccv_prc
  aa_all_auprc_rand_2[i,aacolord] <- aatmmprcv_prc
}
  

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res2_plots/AUROC_rand_heatmap.png", 
    width = 14*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
range(aa_all_auroc_rand_2, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(25)
heatmap.2((aa_all_auroc_rand_2), trace='none', col=(cols), breaks = seq(0.5, 1, 0.02), margins = c(12,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res2_plots/AUROC_chunk_heatmap.png", 
    width = 14*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#sum(is.na(aa_all_auroc_chunk_cvfirst))
#aa_all_auroc_chunk_cvfirst[is.na(aa_all_auroc_chunk_cvfirst)] <- 0
range(aa_all_auroc_chunk_2, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(25)
heatmap.2((aa_all_auroc_chunk_2), trace='none', col=(cols), breaks = seq(0.5, 1, 0.02), margins = c(12,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res2_plots/AUPRC_rand_heatmap.png", 
    width = 14*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#sum(is.na(aa_all_auprc_rand_cvfirst))
#aa_all_auprc_rand_cvfirst[is.na(aa_all_auprc_rand_cvfirst)] <- 0
range(aa_all_auprc_rand_2, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(33)
heatmap.2((aa_all_auprc_rand_2), trace='none', col=(cols), breaks = seq(0.0, 0.66, 0.02), margins = c(12,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res2_plots/AUPRC_chunk_heatmap.png", 
    width = 14*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#sum(is.na(aa_all_auprc_chunk_cvfirst))
#aa_all_auprc_chunk_cvfirst[is.na(aa_all_auprc_chunk_cvfirst)] <- 0
range(aa_all_auprc_chunk_2, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(33)
heatmap.2((aa_all_auprc_chunk_2), trace='none', col=(cols), breaks = seq(0.0, 0.66, 0.02), margins = c(12,10))
dev.off()

family_mat_list_all2 <- list()

for(i in 1:length(aa_all_auroc_list2)){
  print(aafn[i])
  family_mat_list <- list()
  aana <- rownames(aa_all_auroc_list2[[i]])
  aana_sp <- strsplit(aana, "_")
  aana_sp_sort <- lapply(aana_sp, sort)
  aana_allfam <- sort(unique(unlist(aana_sp)))
  #if("sequ" %in% aana_allfam){
    aana_allfam <- setdiff(aana_allfam,c("sequ"))
  #}
  
  for(j in 1:length(aana_allfam)){
    aatmp1 <- which(unlist(lapply(aana_sp, function(x) aana_allfam[j] %in% x)))
    if(aana_allfam[j] == "dist"){
      aatmp111 <- which(unlist(lapply(aana_sp[aatmp1], function(x) "sequ" %in% x)))
      aatmp1 <- aatmp1[setdiff(c(1:length(aatmp1)), aatmp111)]
    }
    #print(aana_allfam[j])
    #print(aana[aatmp1])
    aatmp2 <- numeric(length = length(aatmp1))
    for(k in 1:length(aatmp1)){
      if(length(aana_sp[[aatmp1[k]]]) > 1){
        aanew_name <- setdiff(aana_sp_sort[[aatmp1[k]]], aana_allfam[j])
        if((length(aanew_name) == 2) & identical(aanew_name, c("dist", "sequ"))){
          aatmp2[k] <- 0
        }else{
          aastm <- which(unlist(lapply(aana_sp_sort, function(x) identical(aanew_name, x))))
          if(length(aastm) == 0){
            print(paste("nomatch", names(aa_all_auroc_list2)[i], aana_allfam[j],"###", paste(aana_sp_sort[[aatmp1[k]]], collapse = "_"),"###", paste(aanew_name, collapse = "||")))
            aatmp2[k] <- 0
          }else{
            aatmp2[k] <- aastm 
          }
        }
        # if((length(aanew_name) == 2) & all(aanew_name == c("sequ", "dist"))){
        #   aanew_name <- sort(unlist(strsplit("kmerFq_motSc_Rep_pairs_triplx", "_")))
        #   aastm <- which(unlist(lapply(aana_sp_sort, function(x) identical(aanew_name, x))))
        #   if(length(aastm) == 0){
        #     aanew_name <- sort(unlist(strsplit("kmerFq_motSc_Rep_triplx", "_")))
        #     aastm <- which(unlist(lapply(aana_sp_sort, function(x) identical(aanew_name, x))))
        #     if(length(aastm) == 0){
        #       aanew_name <- sort(unlist(strsplit("kmerFq_Rep_triplx", "_")))
        #     }
        #   }
        # }

        
      }else{
        aatmp2[k] <- 0
      }
      
    }
    family_mat_list[[j]] <- cbind(aatmp1, aatmp2)
  }
  names(family_mat_list) <- aana_allfam
  family_mat_list_all2[[i]] <- family_mat_list
}

names(family_mat_list_all2) <- aafn
length(family_mat_list_all2[[1]])

# aa_familycomp_df
# aa_all_auroc_list
# aa_all_auprc_list

aatmp_lnc_df_list2 <- list()
for(i in 1:length(family_mat_list_all2)){
  aatmp_lnc_fam_df_list <- list()
  print(aafn[i])
  for(j in 1:length(family_mat_list_all2[[i]])){#going over families
    print(names(family_mat_list_all2[[i]])[j])
    aatmptroc <- matrix(nrow = nrow(family_mat_list_all2[[i]][[j]]), 
                        ncol = ncol(aa_all_auroc_list2[[i]]))
    aatmptprc <- matrix(nrow = nrow(family_mat_list_all2[[i]][[j]]), 
                        ncol = ncol(aa_all_auprc_list2[[i]]))
    for(k in 1:nrow(family_mat_list_all2[[i]][[j]])){ #going over comparisons
      if(family_mat_list_all2[[i]][[j]][k,2] == 0){
        aatmpcomproc <- rep(0.5, ncol(aa_all_auroc_list2[[i]]))
        aatmpcompprc <- aa_auprc_base[i, ]
      }else{
        aatmpcomproc <- aa_all_auroc_list2[[i]][family_mat_list_all2[[i]][[j]][k,2],]
        aatmpcompprc <- aa_all_auprc_list2[[i]][family_mat_list_all2[[i]][[j]][k,2],]
      }
      aatmptroc[k,] <- aa_all_auroc_list2[[i]][family_mat_list_all2[[i]][[j]][k,1],] - aatmpcomproc
      aatmptprc[k,] <- aa_all_auprc_list2[[i]][family_mat_list_all2[[i]][[j]][k,1],] - aatmpcompprc
    }
    aatmptroc_num <- colSums(!is.na(aatmptroc))
    aatmptprc_num <- colSums(!is.na(aatmptprc))
    
    aatmptroc_mean <- colMeans(aatmptroc, na.rm = T)
    aatmptprc_mean <- colMeans(aatmptprc, na.rm = T)
    aatmptroc_sd <- apply(aatmptroc, MARGIN = 2, FUN = sd, na.rm = T)
    aatmptprc_sd <- apply(aatmptprc, MARGIN = 2, FUN = sd, na.rm = T)
    aatmptroc_max <- apply(aatmptroc, MARGIN = 2, FUN = max, na.rm = T)
    aatmptprc_max <- apply(aatmptprc, MARGIN = 2, FUN = max, na.rm = T)
    
    aatmp_lnc_fam_df_list[[j]] <- data.frame(lncRNA = rep(aafn[i], length(aatmptroc_mean)), 
                                             feature_family = rep(names(family_mat_list_all2[[i]])[j], length(aatmptroc_mean)), 
                                             partition=c(paste0("CCV", c(1:5)),paste0("RCV", c(1:5))),
                                             partition_type = c(rep("chunk", 5), rep("random", 5)),
                                             nu_comparison = aatmptroc_num, 
                                             auroc_mean_delta = aatmptroc_mean,
                                             auroc_sd_delta = aatmptroc_sd,
                                             auprc_mean_delta = aatmptprc_mean,
                                             auprc_sd_delta = aatmptprc_sd,
                                             auroc_max_delta = aatmptroc_max,
                                             auprc_max_delta = aatmptprc_max)
    
  }
  aatmp_lnc_df_list2[[i]] <- do.call(rbind, aatmp_lnc_fam_df_list)
  
}

aatmp_lnc_df2 <- do.call(rbind, aatmp_lnc_df_list2)
View(aatmp_lnc_df2)



aatttgc <- summarySE(aatmp_lnc_df2, measurevar=c("auroc_max_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc_1 <-  aatttgc[aatttgc$partition_type == "random",]
aatttgc_2 <-  aatttgc[aatttgc$partition_type == "chunk",]

#aatgc <- aatgc[aatgc$Gene_number %in% c(100, 500, 1000),]
ggplot(aatttgc_1, aes(x=lncRNA, y=auroc_max_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=auroc_max_delta-se, ymax=auroc_max_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auroc random")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

ggplot(aatttgc_2, aes(x=lncRNA, y=auroc_mean_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=auroc_mean_delta-se, ymax=auroc_mean_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auroc chunk")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))



aatttgc2 <- summarySE(aatmp_lnc_df2, measurevar=c("auprc_mean_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc2_1 <-  aatttgc2[aatttgc2$partition_type == "random",]
aatttgc2_2 <-  aatttgc2[aatttgc2$partition_type == "chunk",]

ggplot(aatttgc2_1, aes(x=lncRNA, y=auprc_mean_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=auprc_mean_delta-se, ymax=auprc_mean_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auprc random")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

ggplot(aatttgc2_2, aes(x=lncRNA, y=auprc_mean_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=auprc_mean_delta-se, ymax=auprc_mean_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auprc chunk")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

# create for matrices delta AUROC/PRC in Chunk/Random
# each row a lncRNA
# each column a Feature family
aatttgc <- summarySE(aatmp_lnc_df2, measurevar=c("auroc_mean_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc_1 <-  aatttgc[aatttgc$partition_type == "random",]
aatttgc_2 <-  aatttgc[aatttgc$partition_type == "chunk",]

aatttgc2 <- summarySE(aatmp_lnc_df2, measurevar=c("auprc_mean_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc2_1 <-  aatttgc2[aatttgc2$partition_type == "random",]
aatttgc2_2 <-  aatttgc2[aatttgc2$partition_type == "chunk",]

aa_col_seq <- c("dist", "kmerFq", "motSc","pairs",  "Rep" ,   "triplx")
aa_col_con <- setdiff(colnames(aa_auprc_chunk2), aa_col_seq)

aa_auprc_chunk2 <- matrix(nrow = length(unique(aatmp_lnc_df2$lncRNA)), 
                         ncol = length(unique(aatmp_lnc_df2$feature_family)))

rownames(aa_auprc_chunk2) <- unique(aatmp_lnc_df2$lncRNA)
colnames(aa_auprc_chunk2) <- unique(aatmp_lnc_df2$feature_family)

aa_auprc_rand2 <- matrix(nrow = length(unique(aatmp_lnc_df2$lncRNA)), 
                        ncol = length(unique(aatmp_lnc_df2$feature_family)))

rownames(aa_auprc_rand2) <- unique(aatmp_lnc_df2$lncRNA)
colnames(aa_auprc_rand2) <- unique(aatmp_lnc_df2$feature_family)

aa_auroc_chunk2 <- matrix(nrow = length(unique(aatmp_lnc_df2$lncRNA)), 
                         ncol = length(unique(aatmp_lnc_df2$feature_family)))

rownames(aa_auroc_chunk2) <- unique(aatmp_lnc_df2$lncRNA)
colnames(aa_auroc_chunk2) <- unique(aatmp_lnc_df2$feature_family)

aa_auroc_rand2 <- matrix(nrow = length(unique(aatmp_lnc_df2$lncRNA)), 
                        ncol = length(unique(aatmp_lnc_df2$feature_family)))

rownames(aa_auroc_rand2) <- unique(aatmp_lnc_df2$lncRNA)
colnames(aa_auroc_rand2) <- unique(aatmp_lnc_df2$feature_family)

for(i in 1:nrow(aatttgc_1)){
  aa_auroc_rand2[match(aatttgc_1$lncRNA[i], rownames(aa_auroc_rand2)), 
                match(aatttgc_1$feature_family[i], colnames(aa_auroc_rand2))] <- aatttgc_1$auroc_mean_delta[i]
}
for(i in 1:nrow(aatttgc_2)){
  aa_auroc_chunk2[match(aatttgc_2$lncRNA[i], rownames(aa_auroc_chunk2)), 
                 match(aatttgc_2$feature_family[i], colnames(aa_auroc_chunk2))] <- aatttgc_2$auroc_mean_delta[i]
}
for(i in 1:nrow(aatttgc2_1)){
  aa_auprc_rand2[match(aatttgc2_1$lncRNA[i], rownames(aa_auprc_rand2)), 
                match(aatttgc2_1$feature_family[i], colnames(aa_auprc_rand2))] <- aatttgc2_1$auprc_mean_delta[i]
}
for(i in 1:nrow(aatttgc2_2)){
  aa_auprc_chunk2[match(aatttgc2_2$lncRNA[i], rownames(aa_auprc_chunk2)), 
                 match(aatttgc2_2$feature_family[i], colnames(aa_auprc_chunk2))] <- aatttgc2_2$auprc_mean_delta[i]
}

aa_col_seq <- c("dist", "kmerFq", "motSc","pairs",  "Rep" ,   "triplx")
aa_col_con <- setdiff(colnames(aa_auprc_chunk2), aa_col_seq)


library(gplots)
#heatmap.2(aa_auroc_rand2)
require(RColorBrewer)


#aa_auroc_rand2[is.na(aa_auroc_rand2)] <- 0
png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res2_plots/Family_AUROC_rand_seq.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(19)
# z <- zClust(x=aa_auroc_rand2, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auroc_rand2[,aa_col_seq])
heatmap.2(aa_auroc_rand2[,aa_col_seq], trace='none', col=(cols), breaks = seq(-0.04, 0.34, 0.02), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res2_plots/Family_AUROC_rand_cont.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
range(aa_auroc_rand2[,aa_col_con])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(6)
heatmap.2(aa_auroc_rand2[,aa_col_con], trace='none', col=(cols), breaks = seq(-0.00, 0.06, 0.01), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res2_plots/Family_AUROC_chunk_seq.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(17)
# z <- zClust(x=aa_auroc_rand2, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auroc_chunk2[,aa_col_seq])
heatmap.2(aa_auroc_chunk2[,aa_col_seq], trace='none', col=(cols), breaks = seq(-0.0, 0.34, 0.02), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res2_plots/Family_AUROC_chunk_cont.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
range(aa_auroc_chunk2[,aa_col_con])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(8)
heatmap.2(aa_auroc_chunk2[,aa_col_con], trace='none', col=(cols), breaks = seq(-0.00, 0.08, 0.01), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res2_plots/Family_AUPRC_rand_seq.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(20)
# z <- zClust(x=aa_auroc_rand2, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auprc_rand2[,aa_col_seq])
heatmap.2(aa_auprc_rand2[,aa_col_seq], trace='none', col=(cols), breaks = seq(-0.00, 0.40, 0.02), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res2_plots/Family_AUPRC_rand_cont.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
range(aa_auprc_rand2[,aa_col_con])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(13)
heatmap.2(aa_auprc_rand2[,aa_col_con], trace='none', col=(cols), breaks = seq(-0.00, 0.13, 0.01), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res2_plots/Family_AUPRC_chunk_seq.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(18)
# z <- zClust(x=aa_auroc_rand2, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auprc_chunk2[,aa_col_seq])
heatmap.2(aa_auprc_chunk2[,aa_col_seq], trace='none', col=(cols), breaks = seq(-0.00, 0.36, 0.02), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res2_plots/Family_AUPRC_chunk_cont.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
range(aa_auprc_chunk2[,aa_col_con])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(13)
heatmap.2(aa_auprc_chunk2[,aa_col_con], trace='none', col=(cols), breaks = seq(-0.00, 0.13, 0.01), margins = c(6,10))
dev.off()


########################################
# replot with max instead of mean
aatttgc <- summarySE(aatmp_lnc_df2, measurevar=c("auroc_max_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc_1 <-  aatttgc[aatttgc$partition_type == "random",]
aatttgc_2 <-  aatttgc[aatttgc$partition_type == "chunk",]

aatttgc2 <- summarySE(aatmp_lnc_df2, measurevar=c("auprc_max_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc2_1 <-  aatttgc2[aatttgc2$partition_type == "random",]
aatttgc2_2 <-  aatttgc2[aatttgc2$partition_type == "chunk",]

aa_auprc_chunk2 <- matrix(nrow = length(unique(aatmp_lnc_df2$lncRNA)), 
                          ncol = length(unique(aatmp_lnc_df2$feature_family)))

rownames(aa_auprc_chunk2) <- unique(aatmp_lnc_df2$lncRNA)
colnames(aa_auprc_chunk2) <- unique(aatmp_lnc_df2$feature_family)

aa_auprc_rand2 <- matrix(nrow = length(unique(aatmp_lnc_df2$lncRNA)), 
                         ncol = length(unique(aatmp_lnc_df2$feature_family)))

rownames(aa_auprc_rand2) <- unique(aatmp_lnc_df2$lncRNA)
colnames(aa_auprc_rand2) <- unique(aatmp_lnc_df2$feature_family)

aa_auroc_chunk2 <- matrix(nrow = length(unique(aatmp_lnc_df2$lncRNA)), 
                          ncol = length(unique(aatmp_lnc_df2$feature_family)))

rownames(aa_auroc_chunk2) <- unique(aatmp_lnc_df2$lncRNA)
colnames(aa_auroc_chunk2) <- unique(aatmp_lnc_df2$feature_family)

aa_auroc_rand2 <- matrix(nrow = length(unique(aatmp_lnc_df2$lncRNA)), 
                         ncol = length(unique(aatmp_lnc_df2$feature_family)))

rownames(aa_auroc_rand2) <- unique(aatmp_lnc_df2$lncRNA)
colnames(aa_auroc_rand2) <- unique(aatmp_lnc_df2$feature_family)

for(i in 1:nrow(aatttgc_1)){
  aa_auroc_rand2[match(aatttgc_1$lncRNA[i], rownames(aa_auroc_rand2)), 
                 match(aatttgc_1$feature_family[i], colnames(aa_auroc_rand2))] <- aatttgc_1$auroc_max_delta[i]
}
for(i in 1:nrow(aatttgc_2)){
  aa_auroc_chunk2[match(aatttgc_2$lncRNA[i], rownames(aa_auroc_chunk2)), 
                  match(aatttgc_2$feature_family[i], colnames(aa_auroc_chunk2))] <- aatttgc_2$auroc_max_delta[i]
}
for(i in 1:nrow(aatttgc2_1)){
  aa_auprc_rand2[match(aatttgc2_1$lncRNA[i], rownames(aa_auprc_rand2)), 
                 match(aatttgc2_1$feature_family[i], colnames(aa_auprc_rand2))] <- aatttgc2_1$auprc_max_delta[i]
}
for(i in 1:nrow(aatttgc2_2)){
  aa_auprc_chunk2[match(aatttgc2_2$lncRNA[i], rownames(aa_auprc_chunk2)), 
                  match(aatttgc2_2$feature_family[i], colnames(aa_auprc_chunk2))] <- aatttgc2_2$auprc_max_delta[i]
}

library(gplots)
#heatmap.2(aa_auroc_rand2)
require(RColorBrewer)


png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res2_plots/Family_AUROC_rand_seq_max.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(24)
# z <- zClust(x=aa_auroc_rand2, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auroc_rand2[,aa_col_seq])
heatmap.2(aa_auroc_rand2[,aa_col_seq], trace='none', col=(cols), breaks = seq(-0.00, 0.48, 0.02), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res2_plots/Family_AUROC_rand_cont_max.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
range(aa_auroc_rand2[,aa_col_con])
aa_auroc_rand2[aa_auroc_rand2 < 0] <- 0
range(aa_auroc_rand2[,aa_col_con])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(23)
heatmap.2(aa_auroc_rand2[,aa_col_con], trace='none', col=(cols), breaks = seq(-0.00, 0.46, 0.02), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res2_plots/Family_AUROC_chunk_seq_max.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(23)
# z <- zClust(x=aa_auroc_rand2, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auroc_chunk2[,aa_col_seq])
heatmap.2(aa_auroc_chunk2[,aa_col_seq], trace='none', col=(cols), breaks = seq(-0.0, 0.46, 0.02), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res2_plots/Family_AUROC_chunk_cont_max.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
range(aa_auroc_chunk2[,aa_col_con])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(20)
heatmap.2(aa_auroc_chunk2[,aa_col_con], trace='none', col=(cols), breaks = seq(-0.00, 0.40, 0.02), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res2_plots/Family_AUPRC_rand_seq_max.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(33)
# z <- zClust(x=aa_auroc_rand2, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auprc_rand2[,aa_col_seq])
heatmap.2(aa_auprc_rand2[,aa_col_seq], trace='none', col=(cols), breaks = seq(-0.00, 0.66, 0.02), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res2_plots/Family_AUPRC_rand_cont_max.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
range(aa_auprc_rand2[,aa_col_con])
aa_auprc_rand2[aa_auprc_rand2 < 0] <- 0
range(aa_auprc_rand2[,aa_col_con])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(30)
heatmap.2(aa_auprc_rand2[,aa_col_con], trace='none', col=(cols), breaks = seq(-0.00, 0.60, 0.02), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res2_plots/Family_AUPRC_chunk_seq_max.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(27)
# z <- zClust(x=aa_auroc_rand2, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auprc_chunk2[,aa_col_seq])
heatmap.2(aa_auprc_chunk2[,aa_col_seq], trace='none', col=(cols), breaks = seq(-0.00, 0.54, 0.02), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res2_plots/Family_AUPRC_chunk_cont_max.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
range(aa_auprc_chunk2[,aa_col_con])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(28)
heatmap.2(aa_auprc_chunk2[,aa_col_con], trace='none', col=(cols), breaks = seq(-0.00, 0.56, 0.02), margins = c(6,10))
dev.off()



####################
aa_un_own <- sort(unique(Partition_6_random_chunk_cv_df$owner))


for(i in 1:length(aa_un_own)){
  cat(c("Rscript --vanilla /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/gather_perf.R",
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_svd2_lowres10MG_Results/Learned_models/",aa_un_own[i] ),
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_svd2_lowres10MG_Results/Combined_results/",aa_un_own[i], "__perfMat.RData\n" )),
      sep = " ", 
      file = "~/Documents/Shayan/BioInf/lncRNA/gather_perf_all_svd2_10MG.job",
      append = !(i == 1))
}
################################################################################################
# read combined performance for 28 lncRNAs (family svd2) --> distance lower resolution --> 10MG
aafiles <- list.files("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res3/")
aafn <- unlist(lapply(strsplit(aafiles, "__"), "[[", 1))
aa_all_auroc_list3 <- list()
aa_all_auprc_list3 <- list()
aa_all_impor_list3 <- list()

for(i in 1:length(aafiles)){
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res3/", aafiles[i]))
  aa_all_auroc_list3[[i]] <- my_perf[[1]][[1]]
  aa_all_auprc_list3[[i]] <- my_perf[[2]][[1]]
  aa_all_impor_list3[[i]] <- my_perf[[3]][[1]]
}
names(aa_all_auroc_list3) <- aafn
names(aa_all_auprc_list3) <- aafn
names(aa_all_impor_list3) <- aafn

View(aa_all_auroc_list3$Malat1)

# aa_auprc_base <- matrix(nrow = length(aafn), ncol = (ncol(Partition_6_random_chunk_cv_df) - 3))
# rownames(aa_auprc_base) <- aafn
# colnames(aa_auprc_base) <- colnames(Partition_6_random_chunk_cv_df)[4:ncol(Partition_6_random_chunk_cv_df)]
# 
# for(i in 1:nrow(aa_auprc_base)){
#   print(i)
#   for(j in 1:ncol(aa_auprc_base)){
#     aatst <- Partition_6_random_chunk_cv_df[Partition_6_random_chunk_cv_df[, 3+j] == 0, c(1:3)]
#     
#     aa_auprc_base[i, j] <- sum((aatst$label == 1) & (aatst$owner == aafn[i]))/sum(aatst$owner == aafn[i])
#   }
# }
# View(aa_auprc_base)

# aa_auprc_base <- aa_auprc_base[, c(c(6:10), c(1:5))]







aac <- max(unlist(lapply(aa_all_auprc_list3, nrow)))

aa_all_auroc_rand_3 <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auroc_rand_3) <- rownames(aa_all_auprc_list3[[2]])
rownames(aa_all_auroc_rand_3) <- aafn

aa_all_auroc_chunk_3 <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auroc_chunk_3) <- rownames(aa_all_auprc_list3[[2]])
rownames(aa_all_auroc_chunk_3) <- aafn

aa_all_auprc_rand_3 <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auprc_rand_3) <- rownames(aa_all_auprc_list3[[2]])
rownames(aa_all_auprc_rand_3) <- aafn

aa_all_auprc_chunk_3 <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auprc_chunk_3) <- rownames(aa_all_auprc_list3[[2]])
rownames(aa_all_auprc_chunk_3) <- aafn


aa_all_auprc_list3_diff <- list()
for(i in 1:length(aa_all_auprc_list3)){
  aa_all_auprc_list3_diff[[i]] <- aa_all_auprc_list3[[i]]
  for(j in 1:nrow(aa_all_auprc_list3[[i]])){
    aa_all_auprc_list3_diff[[i]][j,] <- aa_all_auprc_list3_diff[[i]][j,] - aa_auprc_base[i,]
  }
}
for(i in 1:length(aa_all_auroc_list3)){
  aatmmpccv_roc <- rowMeans(aa_all_auroc_list3[[i]][,c(1:5)], na.rm = T)
  aatmmprcv_roc <- rowMeans(aa_all_auroc_list3[[i]][,c(6:10)], na.rm = T)
  
  aatmmpccv_prc <- rowMeans(aa_all_auprc_list3_diff[[i]][,c(1:5)], na.rm = T)
  aatmmprcv_prc <- rowMeans(aa_all_auprc_list3_diff[[i]][,c(6:10)], na.rm = T)
  
  aacolord <- match(rownames(aa_all_auroc_list3[[i]]), 
                    colnames(aa_all_auroc_chunk_3))
  
  aa_all_auroc_chunk_3[i,aacolord] <- aatmmpccv_roc
  aa_all_auroc_rand_3[i,aacolord] <- aatmmprcv_roc
  aa_all_auprc_chunk_3[i,aacolord] <- aatmmpccv_prc
  aa_all_auprc_rand_3[i,aacolord] <- aatmmprcv_prc
}


png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res3_plots/AUROC_rand_heatmap.png", 
    width = 14*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
range(aa_all_auroc_rand_3, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(25)
heatmap.2((aa_all_auroc_rand_3), trace='none', col=(cols), breaks = seq(0.5, 1, 0.02), margins = c(12,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res3_plots/AUROC_chunk_heatmap.png", 
    width = 14*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#sum(is.na(aa_all_auroc_chunk_cvfirst))
#aa_all_auroc_chunk_cvfirst[is.na(aa_all_auroc_chunk_cvfirst)] <- 0
range(aa_all_auroc_chunk_3, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(25)
heatmap.2((aa_all_auroc_chunk_3), trace='none', col=(cols), breaks = seq(0.5, 1, 0.02), margins = c(12,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res3_plots/AUPRC_rand_heatmap.png", 
    width = 14*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#sum(is.na(aa_all_auprc_rand_cvfirst))
#aa_all_auprc_rand_cvfirst[is.na(aa_all_auprc_rand_cvfirst)] <- 0
range(aa_all_auprc_rand_3, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(33)
heatmap.2((aa_all_auprc_rand_3), trace='none', col=(cols), breaks = seq(0.0, 0.66, 0.02), margins = c(12,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_Res3_plots/AUPRC_chunk_heatmap.png", 
    width = 14*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#sum(is.na(aa_all_auprc_chunk_cvfirst))
#aa_all_auprc_chunk_cvfirst[is.na(aa_all_auprc_chunk_cvfirst)] <- 0
range(aa_all_auprc_chunk_3, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(33)
heatmap.2((aa_all_auprc_chunk_3), trace='none', col=(cols), breaks = seq(0.0, 0.66, 0.02), margins = c(12,10))
dev.off()
###################################################################################################################################################################################################
############################################################################################################################################################
#####################################################################################################################
# new CV partition 

for(i in 1:length(aa_un_own)){
  cat(c("Rscript --vanilla /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/gather_perf.R",
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_lowres1MG/Learned_models/",aa_un_own[i] ),
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_lowres1MG/Combined_results/",aa_un_own[i], "__perfMat.RData\n" )),
      sep = " ", 
      file = "~/Documents/Shayan/BioInf/lncRNA/gather_perf_all_CV2_svd.job",
      append = !(i == 1))
}

aafiles <- list.files("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_CV2/")
aafn <- unlist(lapply(strsplit(aafiles, "__"), "[[", 1))
aa_all_auroc_list_cv2 <- list()
aa_all_auprc_list_cv2 <- list()
aa_all_impor_list_cv2 <- list()

for(i in 1:length(aafiles)){
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_CV2/", aafiles[i]))
  aa_all_auroc_list_cv2[[i]] <- my_perf[[1]][[1]]
  aa_all_auprc_list_cv2[[i]] <- my_perf[[2]][[1]]
  aa_all_impor_list_cv2[[i]] <- my_perf[[3]][[1]]
}
names(aa_all_auroc_list_cv2) <- aafn
names(aa_all_auprc_list_cv2) <- aafn
names(aa_all_impor_list_cv2) <- aafn

View(aa_all_auroc_list_cv2$Malat1)




aac <- max(unlist(lapply(aa_all_auprc_list_cv2, nrow)))

aa_all_auroc_rand_cv2 <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auroc_rand_cv2) <- rownames(aa_all_auprc_list_cv2[[2]])
rownames(aa_all_auroc_rand_cv2) <- aafn

aa_all_auroc_chunk_cv2 <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auroc_chunk_cv2) <- rownames(aa_all_auprc_list_cv2[[2]])
rownames(aa_all_auroc_chunk_cv2) <- aafn

aa_all_auprc_rand_cv2 <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auprc_rand_cv2) <- rownames(aa_all_auprc_list_cv2[[2]])
rownames(aa_all_auprc_rand_cv2) <- aafn

aa_all_auprc_chunk_cv2 <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auprc_chunk_cv2) <- rownames(aa_all_auprc_list_cv2[[2]])
rownames(aa_all_auprc_chunk_cv2) <- aafn


aa_all_auprc_list_cv2_diff <- list()
for(i in 1:length(aa_all_auprc_list_cv2)){
  aa_all_auprc_list_cv2_diff[[i]] <- aa_all_auprc_list_cv2[[i]]
  for(j in 1:nrow(aa_all_auprc_list_cv2[[i]])){
    aa_all_auprc_list_cv2_diff[[i]][j,] <- aa_all_auprc_list_cv2_diff[[i]][j,] - aa_auprc_base[i,]
  }
}
for(i in 1:length(aa_all_auroc_list_cv2)){
  aatmmpccv_roc <- rowMeans(aa_all_auroc_list_cv2[[i]][,c(1:5)], na.rm = T)
  aatmmprcv_roc <- rowMeans(aa_all_auroc_list_cv2[[i]][,c(6:10)], na.rm = T)
  
  aatmmpccv_prc <- rowMeans(aa_all_auprc_list_cv2_diff[[i]][,c(1:5)], na.rm = T)
  aatmmprcv_prc <- rowMeans(aa_all_auprc_list_cv2_diff[[i]][,c(6:10)], na.rm = T)
  
  aacolord <- match(rownames(aa_all_auroc_list_cv2[[i]]), 
                    colnames(aa_all_auroc_chunk_cv2))
  
  aa_all_auroc_chunk_cv2[i,aacolord] <- aatmmpccv_roc
  aa_all_auroc_rand_cv2[i,aacolord] <- aatmmprcv_roc
  aa_all_auprc_chunk_cv2[i,aacolord] <- aatmmpccv_prc
  aa_all_auprc_rand_cv2[i,aacolord] <- aatmmprcv_prc
}


png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_CV2_plots/AUROC_rand_heatmap.png", 
    width = 14*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
range(aa_all_auroc_rand_cv2, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(25)
heatmap.2((aa_all_auroc_rand_cv2), trace='none', col=(cols), breaks = seq(0.5, 1, 0.02), margins = c(12,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_CV2_plots/AUROC_chunk_heatmap.png", 
    width = 14*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#sum(is.na(aa_all_auroc_chunk_cvfirst))
#aa_all_auroc_chunk_cvfirst[is.na(aa_all_auroc_chunk_cvfirst)] <- 0
range(aa_all_auroc_chunk_cv2, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(25)
heatmap.2((aa_all_auroc_chunk_cv2), trace='none', col=(cols), breaks = seq(0.5, 1, 0.02), margins = c(12,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_CV2_plots/AUPRC_rand_heatmap.png", 
    width = 14*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#sum(is.na(aa_all_auprc_rand_cvfirst))
#aa_all_auprc_rand_cvfirst[is.na(aa_all_auprc_rand_cvfirst)] <- 0
range(aa_all_auprc_rand_cv2, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(33)
heatmap.2((aa_all_auprc_rand_cv2), trace='none', col=(cols), breaks = seq(0.0, 0.66, 0.02), margins = c(12,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_CV2_plots/AUPRC_chunk_heatmap.png", 
    width = 14*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#sum(is.na(aa_all_auprc_chunk_cvfirst))
#aa_all_auprc_chunk_cvfirst[is.na(aa_all_auprc_chunk_cvfirst)] <- 0
range(aa_all_auprc_chunk_cv2, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(33)
heatmap.2((aa_all_auprc_chunk_cv2), trace='none', col=(cols), breaks = seq(0.0, 0.66, 0.02), margins = c(12,10))
dev.off()


aatmp_lnc_df_list_cv2 <- list()
for(i in 1:length(family_mat_list_all2)){
  aatmp_lnc_fam_df_list <- list()
  print(aafn[i])
  for(j in 1:length(family_mat_list_all2[[i]])){#going over families
    print(names(family_mat_list_all2[[i]])[j])
    aatmptroc <- matrix(nrow = nrow(family_mat_list_all2[[i]][[j]]), 
                        ncol = ncol(aa_all_auroc_list_cv2[[i]]))
    aatmptprc <- matrix(nrow = nrow(family_mat_list_all2[[i]][[j]]), 
                        ncol = ncol(aa_all_auprc_list_cv2[[i]]))
    for(k in 1:nrow(family_mat_list_all2[[i]][[j]])){ #going over comparisons
      if(family_mat_list_all2[[i]][[j]][k,2] == 0){
        aatmpcomproc <- rep(0.5, ncol(aa_all_auroc_list_cv2[[i]]))
        aatmpcompprc <- aa_auprc_base[i, ]
      }else{
        aatmpcomproc <- aa_all_auroc_list_cv2[[i]][family_mat_list_all2[[i]][[j]][k,2],]
        aatmpcompprc <- aa_all_auprc_list_cv2[[i]][family_mat_list_all2[[i]][[j]][k,2],]
      }
      aatmptroc[k,] <- aa_all_auroc_list_cv2[[i]][family_mat_list_all2[[i]][[j]][k,1],] - aatmpcomproc
      aatmptprc[k,] <- aa_all_auprc_list_cv2[[i]][family_mat_list_all2[[i]][[j]][k,1],] - aatmpcompprc
    }
    aatmptroc_num <- colSums(!is.na(aatmptroc))
    aatmptprc_num <- colSums(!is.na(aatmptprc))
    
    aatmptroc_mean <- colMeans(aatmptroc, na.rm = T)
    aatmptprc_mean <- colMeans(aatmptprc, na.rm = T)
    aatmptroc_sd <- apply(aatmptroc, MARGIN = 2, FUN = sd, na.rm = T)
    aatmptprc_sd <- apply(aatmptprc, MARGIN = 2, FUN = sd, na.rm = T)
    aatmptroc_max <- apply(aatmptroc, MARGIN = 2, FUN = max, na.rm = T)
    aatmptprc_max <- apply(aatmptprc, MARGIN = 2, FUN = max, na.rm = T)
    
    aatmp_lnc_fam_df_list[[j]] <- data.frame(lncRNA = rep(aafn[i], length(aatmptroc_mean)), 
                                             feature_family = rep(names(family_mat_list_all2[[i]])[j], length(aatmptroc_mean)), 
                                             partition=c(paste0("CCV", c(1:5)),paste0("RCV", c(1:5))),
                                             partition_type = c(rep("chunk", 5), rep("random", 5)),
                                             nu_comparison = aatmptroc_num, 
                                             auroc_mean_delta = aatmptroc_mean,
                                             auroc_sd_delta = aatmptroc_sd,
                                             auprc_mean_delta = aatmptprc_mean,
                                             auprc_sd_delta = aatmptprc_sd,
                                             auroc_max_delta = aatmptroc_max,
                                             auprc_max_delta = aatmptprc_max)
    
  }
  aatmp_lnc_df_list_cv2[[i]] <- do.call(rbind, aatmp_lnc_fam_df_list)
  
}

aatmp_lnc_df_cv2 <- do.call(rbind, aatmp_lnc_df_list_cv2)
View(aatmp_lnc_df_cv2)



aatttgc <- summarySE(aatmp_lnc_df_cv2, measurevar=c("auroc_mean_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc_1 <-  aatttgc[aatttgc$partition_type == "random",]
aatttgc_2 <-  aatttgc[aatttgc$partition_type == "chunk",]

#aatgc <- aatgc[aatgc$Gene_number %in% c(100, 500, 1000),]
ggplot(aatttgc_1, aes(x=lncRNA, y=auroc_mean_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=auroc_mean_delta-se, ymax=auroc_mean_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auroc random")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

ggplot(aatttgc_2, aes(x=lncRNA, y=auroc_mean_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=auroc_mean_delta-se, ymax=auroc_mean_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auroc chunk")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))



aatttgc2 <- summarySE(aatmp_lnc_df_cv2, measurevar=c("auprc_mean_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc2_1 <-  aatttgc2[aatttgc2$partition_type == "random",]
aatttgc2_2 <-  aatttgc2[aatttgc2$partition_type == "chunk",]

ggplot(aatttgc2_1, aes(x=lncRNA, y=auprc_mean_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=auprc_mean_delta-se, ymax=auprc_mean_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auprc random")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

ggplot(aatttgc2_2, aes(x=lncRNA, y=auprc_mean_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=auprc_mean_delta-se, ymax=auprc_mean_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auprc chunk")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

# create for matrices delta AUROC/PRC in Chunk/Random
# each row a lncRNA
# each column a Feature family
aatttgc <- summarySE(aatmp_lnc_df_cv2, measurevar=c("auroc_mean_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc_1 <-  aatttgc[aatttgc$partition_type == "random",]
aatttgc_2 <-  aatttgc[aatttgc$partition_type == "chunk",]

aatttgc2 <- summarySE(aatmp_lnc_df_cv2, measurevar=c("auprc_mean_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc2_1 <-  aatttgc2[aatttgc2$partition_type == "random",]
aatttgc2_2 <-  aatttgc2[aatttgc2$partition_type == "chunk",]

aa_col_seq <- c("dist", "kmerFq", "motSc","pairs",  "Rep" ,   "triplx")
aa_col_con <- setdiff(colnames(aa_auprc_chunk2), aa_col_seq)

aa_auprc_chunk2 <- matrix(nrow = length(unique(aatmp_lnc_df_cv2$lncRNA)), 
                          ncol = length(unique(aatmp_lnc_df_cv2$feature_family)))

rownames(aa_auprc_chunk2) <- unique(aatmp_lnc_df_cv2$lncRNA)
colnames(aa_auprc_chunk2) <- unique(aatmp_lnc_df_cv2$feature_family)

aa_auprc_rand2 <- matrix(nrow = length(unique(aatmp_lnc_df_cv2$lncRNA)), 
                         ncol = length(unique(aatmp_lnc_df_cv2$feature_family)))

rownames(aa_auprc_rand2) <- unique(aatmp_lnc_df_cv2$lncRNA)
colnames(aa_auprc_rand2) <- unique(aatmp_lnc_df_cv2$feature_family)

aa_auroc_chunk2 <- matrix(nrow = length(unique(aatmp_lnc_df_cv2$lncRNA)), 
                          ncol = length(unique(aatmp_lnc_df_cv2$feature_family)))

rownames(aa_auroc_chunk2) <- unique(aatmp_lnc_df_cv2$lncRNA)
colnames(aa_auroc_chunk2) <- unique(aatmp_lnc_df_cv2$feature_family)

aa_auroc_rand2 <- matrix(nrow = length(unique(aatmp_lnc_df_cv2$lncRNA)), 
                         ncol = length(unique(aatmp_lnc_df_cv2$feature_family)))

rownames(aa_auroc_rand2) <- unique(aatmp_lnc_df_cv2$lncRNA)
colnames(aa_auroc_rand2) <- unique(aatmp_lnc_df_cv2$feature_family)

for(i in 1:nrow(aatttgc_1)){
  aa_auroc_rand2[match(aatttgc_1$lncRNA[i], rownames(aa_auroc_rand2)), 
                 match(aatttgc_1$feature_family[i], colnames(aa_auroc_rand2))] <- aatttgc_1$auroc_mean_delta[i]
}
for(i in 1:nrow(aatttgc_2)){
  aa_auroc_chunk2[match(aatttgc_2$lncRNA[i], rownames(aa_auroc_chunk2)), 
                  match(aatttgc_2$feature_family[i], colnames(aa_auroc_chunk2))] <- aatttgc_2$auroc_mean_delta[i]
}
for(i in 1:nrow(aatttgc2_1)){
  aa_auprc_rand2[match(aatttgc2_1$lncRNA[i], rownames(aa_auprc_rand2)), 
                 match(aatttgc2_1$feature_family[i], colnames(aa_auprc_rand2))] <- aatttgc2_1$auprc_mean_delta[i]
}
for(i in 1:nrow(aatttgc2_2)){
  aa_auprc_chunk2[match(aatttgc2_2$lncRNA[i], rownames(aa_auprc_chunk2)), 
                  match(aatttgc2_2$feature_family[i], colnames(aa_auprc_chunk2))] <- aatttgc2_2$auprc_mean_delta[i]
}

aa_col_seq <- c("dist", "kmerFq", "motSc","pairs",  "Rep" ,   "triplx")
aa_col_con <- setdiff(colnames(aa_auprc_chunk2), aa_col_seq)


library(gplots)
#heatmap.2(aa_auroc_rand2)
require(RColorBrewer)


#aa_auroc_rand2[is.na(aa_auroc_rand2)] <- 0
png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_CV2_plots/Family_AUROC_rand_seq.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(19)
# z <- zClust(x=aa_auroc_rand2, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auroc_rand2[,aa_col_seq])
heatmap.2(aa_auroc_rand2[,aa_col_seq], trace='none', col=(cols), breaks = seq(-0.04, 0.34, 0.02), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_CV2_plots/Family_AUROC_rand_cont.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
range(aa_auroc_rand2[,aa_col_con])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(6)
heatmap.2(aa_auroc_rand2[,aa_col_con], trace='none', col=(cols), breaks = seq(-0.00, 0.06, 0.01), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_CV2_plots/Family_AUROC_chunk_seq.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(17)
# z <- zClust(x=aa_auroc_rand2, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auroc_chunk2[,aa_col_seq])
heatmap.2(aa_auroc_chunk2[,aa_col_seq], trace='none', col=(cols), breaks = seq(-0.0, 0.34, 0.02), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_CV2_plots/Family_AUROC_chunk_cont.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
range(aa_auroc_chunk2[,aa_col_con])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(8)
heatmap.2(aa_auroc_chunk2[,aa_col_con], trace='none', col=(cols), breaks = seq(-0.00, 0.08, 0.01), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_CV2_plots/Family_AUPRC_rand_seq.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(20)
# z <- zClust(x=aa_auroc_rand2, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auprc_rand2[,aa_col_seq])
heatmap.2(aa_auprc_rand2[,aa_col_seq], trace='none', col=(cols), breaks = seq(-0.00, 0.40, 0.02), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_CV2_plots/Family_AUPRC_rand_cont.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
range(aa_auprc_rand2[,aa_col_con])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(13)
heatmap.2(aa_auprc_rand2[,aa_col_con], trace='none', col=(cols), breaks = seq(-0.00, 0.13, 0.01), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_CV2_plots/Family_AUPRC_chunk_seq.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(18)
# z <- zClust(x=aa_auroc_rand2, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auprc_chunk2[,aa_col_seq])
heatmap.2(aa_auprc_chunk2[,aa_col_seq], trace='none', col=(cols), breaks = seq(-0.00, 0.36, 0.02), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_CV2_plots/Family_AUPRC_chunk_cont.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
range(aa_auprc_chunk2[,aa_col_con])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(13)
heatmap.2(aa_auprc_chunk2[,aa_col_con], trace='none', col=(cols), breaks = seq(-0.00, 0.13, 0.01), margins = c(6,10))
dev.off()


aaln <- unique(aatmp_lnc_df_cv2$lncRNA)
# write jobs to visualize lnRNA preds
# Rscript plot_gviz_lncRNA.R Meg3 C /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_lowres1MG/Gviz dist kmerFq motSc_pairs_triplx sequ_dist sequ_dist_trnsp sequ_dist_Chrom_Meth_AX_ChIP_t

for(i in 1:length(aaln)){
  cat(c("Rscript",
        "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/plot_gviz_lncRNA.R",
        aaln[i],
        "R", 
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_lowres1MG/Gviz/", aaln[i]), 
        "dist kmerFq motSc_pairs_triplx sequ_dist sequ_dist_trnsp sequ_dist_Chrom_Meth_AX_ChIP_trnsp", 
        "\n"), file = "plot_gviz_all_cv2.job", append = !(i==1))
  cat(c("Rscript",
        "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/plot_gviz_lncRNA.R",
        aaln[i],
        "C", 
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_lowres1MG/Gviz/", aaln[i]), 
        "dist kmerFq motSc_pairs_triplx sequ_dist sequ_dist_trnsp sequ_dist_Chrom_Meth_AX_ChIP_trnsp", 
        "\n"), file = "plot_gviz_all_cv2.job", append = T)
  
}
###########################
# Reading results after transcription filtering
# new CV partition 

for(i in 1:length(aa_un_own)){
  cat(c("Rscript --vanilla /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/gather_perf.R",
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_lowres1MG/Learned_models/",aa_un_own[i] ),
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_lowres1MG/Combined_results/",aa_un_own[i], "__perfMat.RData\n" )),
      sep = " ", 
      file = "~/Documents/Shayan/BioInf/lncRNA/gather_perf_all_CV2_svd.job",
      append = !(i == 1))
}

aafiles <- list.files("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX/")
aafn <- unlist(lapply(strsplit(aafiles, "__"), "[[", 1))
aa_all_auroc_list_cv2TX <- list()
aa_all_auprc_list_cv2TX <- list()
aa_all_impor_list_cv2TX<- list()

for(i in 1:length(aafiles)){
  print(i)
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX/", aafiles[i]))
  aa_all_auroc_list_cv2TX[[i]] <- my_perf[[1]][[1]]
  aa_all_auprc_list_cv2TX[[i]] <- my_perf[[2]][[1]]
  aa_all_impor_list_cv2TX[[i]] <- my_perf[[3]][[1]]
}
names(aa_all_auroc_list_cv2TX) <- aafn
names(aa_all_auprc_list_cv2TX) <- aafn
names(aa_all_impor_list_cv2TX) <- aafn

View(aa_all_auroc_list_cv2TX$Malat1)




aac <- max(unlist(lapply(aa_all_auprc_list_cv2TX, nrow)))

aa_all_auroc_rand_cv2TX <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auroc_rand_cv2TX) <- rownames(aa_all_auprc_list_cv2TX[[2]])
rownames(aa_all_auroc_rand_cv2TX) <- aafn

aa_all_auroc_chunk_cv2TX <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auroc_chunk_cv2TX) <- rownames(aa_all_auprc_list_cv2TX[[2]])
rownames(aa_all_auroc_chunk_cv2TX) <- aafn

aa_all_auprc_rand_cv2TX <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auprc_rand_cv2TX) <- rownames(aa_all_auprc_list_cv2TX[[2]])
rownames(aa_all_auprc_rand_cv2TX) <- aafn

aa_all_auprc_chunk_cv2TX <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auprc_chunk_cv2TX) <- rownames(aa_all_auprc_list_cv2TX[[2]])
rownames(aa_all_auprc_chunk_cv2TX) <- aafn


aa_all_auprc_list_cv2_diffTX <- list()
for(i in 1:length(aa_all_auprc_list_cv2TX)){
  aa_all_auprc_list_cv2_diffTX[[i]] <- aa_all_auprc_list_cv2TX[[i]]
  for(j in 1:nrow(aa_all_auprc_list_cv2TX[[i]])){
    aa_all_auprc_list_cv2_diffTX[[i]][j,] <- aa_all_auprc_list_cv2_diffTX[[i]][j,] - aa_auprc_base[i,]
  }
}
for(i in 1:length(aa_all_auroc_list_cv2TX)){
  aatmmpccv_roc <- rowMeans(aa_all_auroc_list_cv2TX[[i]][,c(1:5)], na.rm = T)
  aatmmprcv_roc <- rowMeans(aa_all_auroc_list_cv2TX[[i]][,c(6:10)], na.rm = T)
  
  aatmmpccv_prc <- rowMeans(aa_all_auprc_list_cv2_diffTX[[i]][,c(1:5)], na.rm = T)
  aatmmprcv_prc <- rowMeans(aa_all_auprc_list_cv2_diffTX[[i]][,c(6:10)], na.rm = T)
  
  aacolord <- match(rownames(aa_all_auroc_list_cv2TX[[i]]), 
                    colnames(aa_all_auroc_chunk_cv2TX))
  
  aa_all_auroc_chunk_cv2TX[i,aacolord] <- aatmmpccv_roc
  aa_all_auroc_rand_cv2TX[i,aacolord] <- aatmmprcv_roc
  aa_all_auprc_chunk_cv2TX[i,aacolord] <- aatmmpccv_prc
  aa_all_auprc_rand_cv2TX[i,aacolord] <- aatmmprcv_prc
}


png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/AUROC_rand_heatmap.png", 
    width = 17*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
range(aa_all_auroc_rand_cv2TX, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(25)
heatmap.2((aa_all_auroc_rand_cv2TX), trace='none', col=(cols), breaks = seq(0.5, 1, 0.02), margins = c(12,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/AUROC_chunk_heatmap.png", 
    width = 17*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#sum(is.na(aa_all_auroc_chunk_cvfirst))
#aa_all_auroc_chunk_cvfirst[is.na(aa_all_auroc_chunk_cvfirst)] <- 0
range(aa_all_auroc_chunk_cv2TX, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(25)
heatmap.2((aa_all_auroc_chunk_cv2TX), trace='none', col=(cols), breaks = seq(0.5, 1, 0.02), margins = c(12,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/AUPRC_rand_heatmap.png", 
    width = 17*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#sum(is.na(aa_all_auprc_rand_cvfirst))
#aa_all_auprc_rand_cvfirst[is.na(aa_all_auprc_rand_cvfirst)] <- 0
range(aa_all_auprc_rand_cv2TX, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(33)
heatmap.2((aa_all_auprc_rand_cv2TX), trace='none', col=(cols), breaks = seq(0.0, 0.66, 0.02), margins = c(12,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/AUPRC_chunk_heatmap.png", 
    width = 17*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#sum(is.na(aa_all_auprc_chunk_cvfirst))
#aa_all_auprc_chunk_cvfirst[is.na(aa_all_auprc_chunk_cvfirst)] <- 0
range(aa_all_auprc_chunk_cv2TX, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(33)
heatmap.2((aa_all_auprc_chunk_cv2TX), trace='none', col=(cols), breaks = seq(0.0, 0.66, 0.02), margins = c(12,10))
dev.off()

#######
# make barplots of best models with distance and best ones without distance
# make a dataframe with the following columns: lncRNA, partition_type, AUROC_dist, AUROC_nodist, AUPRC_dist, AUPRC_nodist
aa_perf_df_o <- data.frame(lncRNA = character(0), 
                           model_partition= character(0), # can be "Wdist_chunk", "Wdist_rand", "WOdist_chunk", "WOdist_rand"
                           AUROC = numeric(0), 
                           AUPRC = numeric(0))

aad <- grep(pattern = "dist", x = colnames(aa_all_auroc_rand_cv2TX))
#colnames(aa_all_auroc_rand_cv2TX)
aa_all_c <- which(colnames(aa_all_auroc_rand_cv2TX) %in% "sequ_dist_Chrom_Meth_AX_ChIP_trnsp")
for(i in 1:nrow(aa_all_auroc_rand_cv2TX)){
  
  aamd1 <- "Full+dist_rand"
  # aamaxrd_roc <- max(aa_all_auroc_rand_cv2TX[i,aad], na.rm = T)
  # aamaxrd_prc <- max(aa_all_auprc_rand_cv2TX[i,aad], na.rm = T)
  aamaxrd_roc <- aa_all_auroc_rand_cv2TX[i,aa_all_c]
  aamaxrd_prc <- aa_all_auprc_rand_cv2TX[i,aa_all_c]
  
  aaaxx1 <- c(rownames(aa_all_auroc_rand_cv2TX)[i], aamd1, aamaxrd_roc, aamaxrd_prc)
  aamd2 <- "Full+dist_chunk"
  # aamaxcd_roc <- max(aa_all_auroc_chunk_cv2TX[i,aad], na.rm = T)
  # aamaxcd_prc <- max(aa_all_auprc_chunk_cv2TX[i,aad], na.rm = T)
  aamaxcd_roc <- aa_all_auroc_chunk_cv2TX[i,aa_all_c]
  aamaxcd_prc <- aa_all_auprc_chunk_cv2TX[i,aa_all_c]
  
  aaaxx2 <- c(rownames(aa_all_auroc_rand_cv2TX)[i], aamd2, aamaxcd_roc, aamaxcd_prc)
  aa_perf_df_o <- rbind(aa_perf_df_o, rbind(aaaxx1,aaaxx2))
}
#######

family_mat_list_all2TX <- list()

for(i in 1:length(aa_all_auroc_list_cv2TX)){
  print(aafn[i])
  family_mat_list <- list()
  aana <- rownames(aa_all_auroc_list_cv2TX[[i]])
  aana_sp <- strsplit(aana, "_")
  aana_sp_sort <- lapply(aana_sp, sort)
  aana_allfam <- sort(unique(unlist(aana_sp)))
  if("sequ" %in% aana_allfam){
   aana_allfam <- setdiff(aana_allfam,c("sequ"))
  }
  
  for(j in 1:length(aana_allfam)){
    aatmp1 <- which(unlist(lapply(aana_sp, function(x) aana_allfam[j] %in% x)))
    if(!(aana_allfam[j] %in% c("TXRB", "TXRE"))){
      aatmp111 <- which(unlist(lapply(aana_sp[aatmp1], function(x) "TXRB" %in% x)))
      aatmp1112 <- which(unlist(lapply(aana_sp[aatmp1], function(x) "TXRE" %in% x)))
      aatmp1 <- aatmp1[setdiff(c(1:length(aatmp1)), union(aatmp111, aatmp1112))]
    }
    #print(aana_allfam[j])
    #print(aana[aatmp1])
    aatmp2 <- numeric(length = length(aatmp1))
    for(k in 1:length(aatmp1)){
      if(length(aana_sp[[aatmp1[k]]]) > 1){
        aanew_name <- setdiff(aana_sp_sort[[aatmp1[k]]], aana_allfam[j])
        #if((length(aanew_name) == 2) & identical(aanew_name, c("dist", "sequ"))){
          #aatmp2[k] <- 0
        #}else{
          aastm <- which(unlist(lapply(aana_sp_sort, function(x) identical(aanew_name, x))))
          if(length(aastm) == 0){
            print(paste("nomatch", names(aa_all_auroc_list_cv2TX)[i], aana_allfam[j],"###", paste(aana_sp_sort[[aatmp1[k]]], collapse = "_"),"###", paste(aanew_name, collapse = "||")))
            aatmp2[k] <- 0
          }else{
            aatmp2[k] <- aastm 
          }
        #}
        # if((length(aanew_name) == 2) & all(aanew_name == c("sequ", "dist"))){
        #   aanew_name <- sort(unlist(strsplit("kmerFq_motSc_Rep_pairs_triplx", "_")))
        #   aastm <- which(unlist(lapply(aana_sp_sort, function(x) identical(aanew_name, x))))
        #   if(length(aastm) == 0){
        #     aanew_name <- sort(unlist(strsplit("kmerFq_motSc_Rep_triplx", "_")))
        #     aastm <- which(unlist(lapply(aana_sp_sort, function(x) identical(aanew_name, x))))
        #     if(length(aastm) == 0){
        #       aanew_name <- sort(unlist(strsplit("kmerFq_Rep_triplx", "_")))
        #     }
        #   }
        # }
        
        
      }else{
        aatmp2[k] <- 0
      }
      
    }
    family_mat_list[[j]] <- cbind(aatmp1, aatmp2)
  }
  names(family_mat_list) <- aana_allfam
  family_mat_list_all2TX[[i]] <- family_mat_list
}

names(family_mat_list_all2TX) <- aafn
length(family_mat_list_all2TX[[1]])





aatmp_lnc_df_list_cv2TX <- list()
for(i in 1:length(family_mat_list_all2TX)){
  aatmp_lnc_fam_df_list <- list()
  print(aafn[i])
  for(j in 1:length(family_mat_list_all2TX[[i]])){#going over families
    print(names(family_mat_list_all2TX[[i]])[j])
    aatmptroc <- matrix(nrow = nrow(family_mat_list_all2TX[[i]][[j]]), 
                        ncol = ncol(aa_all_auroc_list_cv2TX[[i]]))
    aatmptprc <- matrix(nrow = nrow(family_mat_list_all2TX[[i]][[j]]), 
                        ncol = ncol(aa_all_auprc_list_cv2TX[[i]]))
    for(k in 1:nrow(family_mat_list_all2TX[[i]][[j]])){ #going over comparisons
      if(family_mat_list_all2TX[[i]][[j]][k,2] == 0){
        aatmpcomproc <- rep(0.5, ncol(aa_all_auroc_list_cv2TX[[i]]))
        aatmpcompprc <- aa_auprc_base[i, ]
      }else{
        aatmpcomproc <- aa_all_auroc_list_cv2TX[[i]][family_mat_list_all2TX[[i]][[j]][k,2],]
        aatmpcompprc <- aa_all_auprc_list_cv2TX[[i]][family_mat_list_all2TX[[i]][[j]][k,2],]
      }
      aatmptroc[k,] <- aa_all_auroc_list_cv2TX[[i]][family_mat_list_all2TX[[i]][[j]][k,1],] - aatmpcomproc
      aatmptprc[k,] <- aa_all_auprc_list_cv2TX[[i]][family_mat_list_all2TX[[i]][[j]][k,1],] - aatmpcompprc
    }
    aatmptroc_num <- colSums(!is.na(aatmptroc))
    aatmptprc_num <- colSums(!is.na(aatmptprc))
    
    aatmptroc_mean <- colMeans(aatmptroc, na.rm = T)
    aatmptprc_mean <- colMeans(aatmptprc, na.rm = T)
    aatmptroc_sd <- apply(aatmptroc, MARGIN = 2, FUN = sd, na.rm = T)
    aatmptprc_sd <- apply(aatmptprc, MARGIN = 2, FUN = sd, na.rm = T)
    aatmptroc_max <- apply(aatmptroc, MARGIN = 2, FUN = max, na.rm = T)
    aatmptprc_max <- apply(aatmptprc, MARGIN = 2, FUN = max, na.rm = T)
    
    aatmp_lnc_fam_df_list[[j]] <- data.frame(lncRNA = rep(aafn[i], length(aatmptroc_mean)), 
                                             feature_family = rep(names(family_mat_list_all2TX[[i]])[j], length(aatmptroc_mean)), 
                                             partition=c(paste0("CCV", c(1:5)),paste0("RCV", c(1:5))),
                                             partition_type = c(rep("chunk", 5), rep("random", 5)),
                                             nu_comparison = aatmptroc_num, 
                                             auroc_mean_delta = aatmptroc_mean,
                                             auroc_sd_delta = aatmptroc_sd,
                                             auprc_mean_delta = aatmptprc_mean,
                                             auprc_sd_delta = aatmptprc_sd,
                                             auroc_max_delta = aatmptroc_max,
                                             auprc_max_delta = aatmptprc_max)
    
  }
  aatmp_lnc_df_list_cv2TX[[i]] <- do.call(rbind, aatmp_lnc_fam_df_list)
  
}

aatmp_lnc_df_cv2TX <- do.call(rbind, aatmp_lnc_df_list_cv2TX)
View(aatmp_lnc_df_cv2TX)



aatttgc <- summarySE(data = aatmp_lnc_df_cv2TX, 
                     measurevar=c("auroc_mean_delta"), 
                     groupvars=c("lncRNA", "feature_family", "partition_type"), 
                     na.rm=T)

aatttgc_1 <-  aatttgc[aatttgc$partition_type == "random",]
aatttgc_2 <-  aatttgc[aatttgc$partition_type == "chunk",]

#aatgc <- aatgc[aatgc$Gene_number %in% c(100, 500, 1000),]
ggplot(aatttgc_1, aes(x=lncRNA, y=auroc_mean_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=auroc_mean_delta-se, ymax=auroc_mean_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auroc random")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

ggplot(aatttgc_2, aes(x=lncRNA, y=auroc_mean_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=auroc_mean_delta-se, ymax=auroc_mean_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auroc chunk")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))



aatttgc2 <- summarySE(aatmp_lnc_df_cv2TX, measurevar=c("auprc_mean_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc2_1 <-  aatttgc2[aatttgc2$partition_type == "random",]
aatttgc2_2 <-  aatttgc2[aatttgc2$partition_type == "chunk",]

ggplot(aatttgc2_1, aes(x=lncRNA, y=auprc_mean_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=auprc_mean_delta-se, ymax=auprc_mean_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auprc random")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

ggplot(aatttgc2_2, aes(x=lncRNA, y=auprc_mean_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=auprc_mean_delta-se, ymax=auprc_mean_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auprc chunk")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

# create for matrices delta AUROC/PRC in Chunk/Random
# each row a lncRNA
# each column a Feature family
aatttgc <- summarySE(aatmp_lnc_df_cv2TX, measurevar=c("auroc_mean_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc_1 <-  aatttgc[aatttgc$partition_type == "random",]
aatttgc_2 <-  aatttgc[aatttgc$partition_type == "chunk",]

aatttgc2 <- summarySE(aatmp_lnc_df_cv2TX, measurevar=c("auprc_mean_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc2_1 <-  aatttgc2[aatttgc2$partition_type == "random",]
aatttgc2_2 <-  aatttgc2[aatttgc2$partition_type == "chunk",]

aa_col_seq <- c("dist", "kmerFq", "motSc","pairs",  "Rep" ,   "triplx")
#aa_col_con <- setdiff(colnames(aa_auprc_chunk2), aa_col_seq)

aa_auprc_chunk2 <- matrix(nrow = length(unique(aatmp_lnc_df_cv2TX$lncRNA)), 
                          ncol = length(unique(aatmp_lnc_df_cv2TX$feature_family)))

rownames(aa_auprc_chunk2) <- unique(aatmp_lnc_df_cv2TX$lncRNA)
colnames(aa_auprc_chunk2) <- unique(aatmp_lnc_df_cv2TX$feature_family)

aa_auprc_rand2 <- matrix(nrow = length(unique(aatmp_lnc_df_cv2TX$lncRNA)), 
                         ncol = length(unique(aatmp_lnc_df_cv2TX$feature_family)))

rownames(aa_auprc_rand2) <- unique(aatmp_lnc_df_cv2TX$lncRNA)
colnames(aa_auprc_rand2) <- unique(aatmp_lnc_df_cv2TX$feature_family)

aa_auroc_chunk2 <- matrix(nrow = length(unique(aatmp_lnc_df_cv2TX$lncRNA)), 
                          ncol = length(unique(aatmp_lnc_df_cv2TX$feature_family)))

rownames(aa_auroc_chunk2) <- unique(aatmp_lnc_df_cv2TX$lncRNA)
colnames(aa_auroc_chunk2) <- unique(aatmp_lnc_df_cv2TX$feature_family)

aa_auroc_rand2 <- matrix(nrow = length(unique(aatmp_lnc_df_cv2TX$lncRNA)), 
                         ncol = length(unique(aatmp_lnc_df_cv2TX$feature_family)))

rownames(aa_auroc_rand2) <- unique(aatmp_lnc_df_cv2TX$lncRNA)
colnames(aa_auroc_rand2) <- unique(aatmp_lnc_df_cv2TX$feature_family)

for(i in 1:nrow(aatttgc_1)){
  aa_auroc_rand2[match(aatttgc_1$lncRNA[i], rownames(aa_auroc_rand2)), 
                 match(aatttgc_1$feature_family[i], colnames(aa_auroc_rand2))] <- aatttgc_1$auroc_mean_delta[i]
}
for(i in 1:nrow(aatttgc_2)){
  aa_auroc_chunk2[match(aatttgc_2$lncRNA[i], rownames(aa_auroc_chunk2)), 
                  match(aatttgc_2$feature_family[i], colnames(aa_auroc_chunk2))] <- aatttgc_2$auroc_mean_delta[i]
}
for(i in 1:nrow(aatttgc2_1)){
  aa_auprc_rand2[match(aatttgc2_1$lncRNA[i], rownames(aa_auprc_rand2)), 
                 match(aatttgc2_1$feature_family[i], colnames(aa_auprc_rand2))] <- aatttgc2_1$auprc_mean_delta[i]
}
for(i in 1:nrow(aatttgc2_2)){
  aa_auprc_chunk2[match(aatttgc2_2$lncRNA[i], rownames(aa_auprc_chunk2)), 
                  match(aatttgc2_2$feature_family[i], colnames(aa_auprc_chunk2))] <- aatttgc2_2$auprc_mean_delta[i]
}

aa_col_seq <- c("dist", "kmerFq", "motSc","pairs",  "Rep" ,   "triplx")
aa_col_con <- setdiff(colnames(aa_auprc_chunk2), aa_col_seq)


library(gplots)
#heatmap.2(aa_auroc_rand2)
require(RColorBrewer)


#aa_auroc_rand2[is.na(aa_auroc_rand2)] <- 0
png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Family_AUROC_rand_seq.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(19)
# z <- zClust(x=aa_auroc_rand2, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auroc_rand2[,aa_col_seq])
heatmap.2(aa_auroc_rand2[,aa_col_seq], trace='none', col=(cols), breaks = seq(-0.04, 0.34, 0.02), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Family_AUROC_rand_cont.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
range(aa_auroc_rand2[,aa_col_con])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(6)
heatmap.2(aa_auroc_rand2[,aa_col_con], trace='none', col=(cols), breaks = seq(-0.00, 0.06, 0.01), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Family_AUROC_chunk_seq.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(19)
# z <- zClust(x=aa_auroc_rand2, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auroc_chunk2[,aa_col_seq])
heatmap.2(aa_auroc_chunk2[,aa_col_seq], trace='none', col=(cols), breaks = seq(-0.0, 0.38, 0.02), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Family_AUROC_chunk_cont.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
range(aa_auroc_chunk2[,aa_col_con])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(8)
heatmap.2(aa_auroc_chunk2[,aa_col_con], trace='none', col=(cols), breaks = seq(-0.00, 0.08, 0.01), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Family_AUPRC_rand_seq.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(20)
# z <- zClust(x=aa_auroc_rand2, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auprc_rand2[,aa_col_seq])
heatmap.2(aa_auprc_rand2[,aa_col_seq], trace='none', col=(cols), breaks = seq(-0.00, 0.40, 0.02), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Family_AUPRC_rand_cont.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
range(aa_auprc_rand2[,aa_col_con])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(13)
heatmap.2(aa_auprc_rand2[,aa_col_con], trace='none', col=(cols), breaks = seq(-0.00, 0.13, 0.01), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Family_AUPRC_chunk_seq.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(19)
# z <- zClust(x=aa_auroc_rand2, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auprc_chunk2[,aa_col_seq])
heatmap.2(aa_auprc_chunk2[,aa_col_seq], trace='none', col=(cols), breaks = seq(-0.00, 0.38, 0.02), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Family_AUPRC_chunk_cont.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
range(aa_auprc_chunk2[,aa_col_con])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(13)
heatmap.2(aa_auprc_chunk2[,aa_col_con], trace='none', col=(cols), breaks = seq(-0.00, 0.13, 0.01), margins = c(6,10))
dev.off()


aaln <- unique(aatmp_lnc_df_cv2TX$lncRNA)
# write jobs to visualize lnRNA preds
# Rscript plot_gviz_lncRNA.R Meg3 C /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_lowres1MG/Gviz dist kmerFq motSc_pairs_triplx sequ_dist sequ_dist_trnsp sequ_dist_Chrom_Meth_AX_ChIP_t

for(i in 1:length(aaln)){
  cat(c("Rscript",
        "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/plot_gviz_lncRNA.R",
        aaln[i],
        "R", 
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_lowres1MG/Gviz/", aaln[i]), 
        "dist kmerFq motSc_pairs_triplx sequ_dist sequ_dist_trnsp sequ_dist_Chrom_Meth_AX_ChIP_trnsp sequ_dist_Chrom_Meth_AX_ChIP_trnsp_TXRE", 
        "\n"), file = "plot_gviz_all_cv2.job", append = !(i==1))
  cat(c("Rscript",
        "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/plot_gviz_lncRNA.R",
        aaln[i],
        "C", 
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_lowres1MG/Gviz/", aaln[i]), 
        "dist kmerFq motSc_pairs_triplx sequ_dist sequ_dist_trnsp sequ_dist_Chrom_Meth_AX_ChIP_trnsp sequ_dist_Chrom_Meth_AX_ChIP_trnsp_TXRE", 
        "\n"), file = "plot_gviz_all_cv2.job", append = T)
  
}

# plot the best model and 20 most important fearutes

for(i in 1:length(aaln)){
  cat(c("Rscript",
        "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/plot_gviz_OneModel_topFeat.R",
        aaln[i],
        "R", 
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_lowres1MG/Gviz/", aaln[i]), 
        "20", 
        colnames(aa_all_auprc_rand_cv2TX)[which.max(aa_all_auprc_rand_cv2TX[match(aaln[i], rownames(aa_all_auprc_rand_cv2TX)),])], 
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_CVfirst_dataset_", aaln[i], ".RData"),
        "\n"), file = "plot_gviz_top20_cv2TX.job", append = !(i==1))
  cat(c("Rscript",
        "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/plot_gviz_OneModel_topFeat.R",
        aaln[i],
        "C", 
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_lowres1MG/Gviz/", aaln[i]), 
        "20", 
        colnames(aa_all_auprc_chunk_cv2TX)[which.max(aa_all_auprc_chunk_cv2TX[match(aaln[i], rownames(aa_all_auprc_chunk_cv2TX)),])], 
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_CVfirst_dataset_", aaln[i], ".RData"),
        "\n"), file = "plot_gviz_top20_cv2TX.job", append = T)
  
}



colnames(aa_all_auprc_chunk_cv2TX)[which.max(aa_all_auprc_chunk_cv2TX[match(aaln[i], rownames(aa_all_auprc_chunk_cv2TX)),])]


################################################################################################
# Redoing the above analysis excluding distance
################################################################################################


family_mat_list_all2TX_noDist <- list()

aa_all_auroc_list_noDist <- aa_all_auroc_list_cv2TX
aa_all_auprc_list_noDist <- aa_all_auprc_list_cv2TX
for(i in 1:length(aa_all_auroc_list_noDist)){
  print(names(aa_all_auroc_list_cv2TX)[i])
  aagr <- grep(pattern = "dist", ignore.case = T, x = rownames(aa_all_auroc_list_noDist[[i]]))
  print(length(aagr))
  if(length(aagr) > 0){
    aa_all_auroc_list_noDist[[i]] <- aa_all_auroc_list_noDist[[i]][-aagr,]
    aa_all_auprc_list_noDist[[i]] <- aa_all_auprc_list_noDist[[i]][-aagr,]
  }
}

for(i in 1:length(aa_all_auroc_list_noDist)){
  print(aafn[i])
  family_mat_list <- list()
  aana <- rownames(aa_all_auroc_list_noDist[[i]])
  aana_sp <- strsplit(aana, "_")
  aana_sp_sort <- lapply(aana_sp, sort)
  aana_allfam <- sort(unique(unlist(aana_sp)))
  aana_allfam <- setdiff(aana_allfam,"sequ")
  for(j in 1:length(aana_allfam)){
    aatmp1 <- which(unlist(lapply(aana_sp, function(x) aana_allfam[j] %in% x)))
    if(!(aana_allfam[j] %in% c("TXRB", "TXRE"))){
      aatmp111 <- which(unlist(lapply(aana_sp[aatmp1], function(x) "TXRB" %in% x)))
      aatmp1112 <- which(unlist(lapply(aana_sp[aatmp1], function(x) "TXRE" %in% x)))
      aatmp1 <- aatmp1[setdiff(c(1:length(aatmp1)), union(aatmp111, aatmp1112))]
    }
    #print(aana_allfam[j])
    #print(aana[aatmp1])
    aatmp2 <- numeric(length = length(aatmp1))
    for(k in 1:length(aatmp1)){
      if(length(aana_sp[[aatmp1[k]]]) > 1){
        aanew_name <- setdiff(aana_sp_sort[[aatmp1[k]]], aana_allfam[j])
        # if((length(aanew_name) == 1) & (aanew_name[1] == "sequ")){
        #   aanew_name <- sort(unlist(strsplit("kmerFq_motSc_Rep_pairs_triplx", "_")))
        #   aastm <- which(unlist(lapply(aana_sp_sort, function(x) identical(aanew_name, x))))
        #   if(length(aastm) == 0){
        #     aanew_name <- sort(unlist(strsplit("kmerFq_motSc_Rep_triplx", "_")))
        #     aastm <- which(unlist(lapply(aana_sp_sort, function(x) identical(aanew_name, x))))
        #     if(length(aastm) == 0){
        #       aanew_name <- sort(unlist(strsplit("kmerFq_Rep_triplx", "_")))
        #     }
        #   }
        # }
        aastm <- which(unlist(lapply(aana_sp_sort, function(x) identical(aanew_name, x))))
        if(length(aastm) == 0){
          print(paste("nomatch", names(aa_all_auroc_list_noDist)[i], aana_allfam[j],"###", paste(aana_sp_sort[[aatmp1[k]]], 
                                                                                                 collapse = "_"),"###", 
                      paste(aanew_name, collapse = "||")))
          aatmp2[k] <- 0
        }else{
          aatmp2[k] <- aastm 
        }
        
      }else{
        aatmp2[k] <- 0
      }
      
    }
    family_mat_list[[j]] <- cbind(aatmp1, aatmp2)
  }
  names(family_mat_list) <- aana_allfam
  family_mat_list_all2TX_noDist[[i]] <- family_mat_list
}

names(family_mat_list_all2TX_noDist) <- aafn
length(family_mat_list_all2TX_noDist[[1]])

# aa_familycomp_df
# aa_all_auroc_list
# aa_all_auprc_list

aatmp_lnc_df_list_cv2TX_noDist <- list()
for(i in 1:length(family_mat_list_all2TX_noDist)){
  aatmp_lnc_fam_df_list <- list()
  print(aafn[i])
  for(j in 1:length(family_mat_list_all2TX_noDist[[i]])){#going over families
    print(names(family_mat_list_all2TX_noDist[[i]])[j])
    aatmptroc <- matrix(nrow = nrow(family_mat_list_all2TX_noDist[[i]][[j]]), 
                        ncol = ncol(aa_all_auroc_list_noDist[[i]]))
    aatmptprc <- matrix(nrow = nrow(family_mat_list_all2TX_noDist[[i]][[j]]), 
                        ncol = ncol(aa_all_auprc_list_noDist[[i]]))
    for(k in 1:nrow(family_mat_list_all2TX_noDist[[i]][[j]])){ #going over comparisons
      if(family_mat_list_all2TX_noDist[[i]][[j]][k,2] == 0){
        aatmpcomproc <- rep(0.5, ncol(aa_all_auroc_list_noDist[[i]]))
        aatmpcompprc <- aa_auprc_base[i, ]
      }else{
        aatmpcomproc <- aa_all_auroc_list_noDist[[i]][family_mat_list_all2TX_noDist[[i]][[j]][k,2],]
        aatmpcompprc <- aa_all_auprc_list_noDist[[i]][family_mat_list_all2TX_noDist[[i]][[j]][k,2],]
      }
      aatmptroc[k,] <- aa_all_auroc_list_noDist[[i]][family_mat_list_all2TX_noDist[[i]][[j]][k,1],] - aatmpcomproc
      aatmptprc[k,] <- aa_all_auprc_list_noDist[[i]][family_mat_list_all2TX_noDist[[i]][[j]][k,1],] - aatmpcompprc
    }
    aatmptroc_num <- colSums(!is.na(aatmptroc))
    aatmptprc_num <- colSums(!is.na(aatmptprc))
    
    aatmptroc_mean <- colMeans(aatmptroc, na.rm = T)
    aatmptprc_mean <- colMeans(aatmptprc, na.rm = T)
    aatmptroc_sd <- apply(aatmptroc, MARGIN = 2, FUN = sd, na.rm = T)
    aatmptprc_sd <- apply(aatmptprc, MARGIN = 2, FUN = sd, na.rm = T)
    aatmp_lnc_fam_df_list[[j]] <- data.frame(lncRNA = rep(aafn[i], length(aatmptroc_mean)), 
                                             feature_family = rep(names(family_mat_list_all2TX_noDist[[i]])[j], length(aatmptroc_mean)), 
                                             partition=c(paste0("CCV", c(1:5)),paste0("RCV", c(1:5))),
                                             partition_type = c(rep("chunk", 5), rep("random", 5)),
                                             nu_comparison = aatmptroc_num, 
                                             auroc_mean_delta = aatmptroc_mean,
                                             auroc_sd_delta = aatmptroc_sd,
                                             auprc_mean_delta = aatmptprc_mean,
                                             auprc_sd_delta = aatmptprc_sd)
    
  }
  aatmp_lnc_df_list_cv2TX_noDist[[i]] <- do.call(rbind, aatmp_lnc_fam_df_list)
  
}

aatmp_lnc_df_cv2TX_noDist <- do.call(rbind, aatmp_lnc_df_list_cv2TX_noDist)
View(aatmp_lnc_df_cv2TX_noDist)

# create for matrices delta AUROC/PRC in Chunk/Random -- after removing distance models
# each row a lncRNA
# each column a Feature family
aatttgc_nodist <- summarySE(aatmp_lnc_df_cv2TX_noDist, measurevar=c("auroc_mean_delta"), 
                            groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc_1_nodist <-  aatttgc_nodist[aatttgc_nodist$partition_type == "random",]
aatttgc_2_nodist <-  aatttgc_nodist[aatttgc_nodist$partition_type == "chunk",]

aatttgc2_nodist <- summarySE(aatmp_lnc_df_cv2TX_noDist, measurevar=c("auprc_mean_delta"), 
                             groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc2_1_nodist <-  aatttgc2_nodist[aatttgc2_nodist$partition_type == "random",]
aatttgc2_2_nodist <-  aatttgc2_nodist[aatttgc2_nodist$partition_type == "chunk",]

aa_auprc_chunk_nodist <- matrix(nrow = length(unique(aatmp_lnc_df_cv2TX_noDist$lncRNA)), 
                                ncol = length(unique(aatmp_lnc_df_cv2TX_noDist$feature_family)))

rownames(aa_auprc_chunk_nodist) <- unique(aatmp_lnc_df_cv2TX_noDist$lncRNA)
colnames(aa_auprc_chunk_nodist) <- unique(aatmp_lnc_df_cv2TX_noDist$feature_family)

aa_auprc_rand_nodist <- matrix(nrow = length(unique(aatmp_lnc_df_cv2TX_noDist$lncRNA)), 
                               ncol = length(unique(aatmp_lnc_df_cv2TX_noDist$feature_family)))

rownames(aa_auprc_rand_nodist) <- unique(aatmp_lnc_df_cv2TX_noDist$lncRNA)
colnames(aa_auprc_rand_nodist) <- unique(aatmp_lnc_df_cv2TX_noDist$feature_family)

aa_auroc_chunk_nodist <- matrix(nrow = length(unique(aatmp_lnc_df_cv2TX_noDist$lncRNA)), 
                                ncol = length(unique(aatmp_lnc_df_cv2TX_noDist$feature_family)))

rownames(aa_auroc_chunk_nodist) <- unique(aatmp_lnc_df_cv2TX_noDist$lncRNA)
colnames(aa_auroc_chunk_nodist) <- unique(aatmp_lnc_df_cv2TX_noDist$feature_family)

aa_auroc_rand_nodist <- matrix(nrow = length(unique(aatmp_lnc_df_cv2TX_noDist$lncRNA)), 
                               ncol = length(unique(aatmp_lnc_df_cv2TX_noDist$feature_family)))

rownames(aa_auroc_rand_nodist) <- unique(aatmp_lnc_df_cv2TX_noDist$lncRNA)
colnames(aa_auroc_rand_nodist) <- unique(aatmp_lnc_df_cv2TX_noDist$feature_family)

for(i in 1:nrow(aatttgc_1_nodist)){
  aa_auroc_rand_nodist[match(aatttgc_1_nodist$lncRNA[i], rownames(aa_auroc_rand_nodist)), 
                       match(aatttgc_1_nodist$feature_family[i], colnames(aa_auroc_rand_nodist))] <- aatttgc_1_nodist$auroc_mean_delta[i]
}
for(i in 1:nrow(aatttgc_2_nodist)){
  aa_auroc_chunk_nodist[match(aatttgc_2_nodist$lncRNA[i], rownames(aa_auroc_chunk_nodist)), 
                        match(aatttgc_2_nodist$feature_family[i], colnames(aa_auroc_chunk_nodist))] <- aatttgc_2_nodist$auroc_mean_delta[i]
}
for(i in 1:nrow(aatttgc2_1)){
  aa_auprc_rand_nodist[match(aatttgc2_1_nodist$lncRNA[i], rownames(aa_auprc_rand_nodist)), 
                       match(aatttgc2_1_nodist$feature_family[i], colnames(aa_auprc_rand_nodist))] <- aatttgc2_1_nodist$auprc_mean_delta[i]
}
for(i in 1:nrow(aatttgc2_2_nodist)){
  aa_auprc_chunk_nodist[match(aatttgc2_2_nodist$lncRNA[i], rownames(aa_auprc_chunk_nodist)), 
                        match(aatttgc2_2_nodist$feature_family[i], colnames(aa_auprc_chunk_nodist))] <- aatttgc2_2_nodist$auprc_mean_delta[i]
}

library(gplots)
#heatmap.2(aa_auroc_rand)
require(RColorBrewer)


#aa_auroc_rand_nodist[is.na(aa_auroc_rand_nodist)] <- 0
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(12)
# z <- zClust(x=aa_auroc_rand, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auroc_rand_nodist)
aam1 <- heatmap.2(aa_auroc_rand_nodist, trace='none', col=(cols), breaks = seq(-0.00, 0.12, 0.01), margins = c(6,10))
aam1
#aa_auroc_chunk_nodist[is.na(aa_auroc_chunk_nodist)] <- 0
range(aa_auroc_chunk_nodist)
cols <- colorRampPalette(brewer.pal(11, "RdBu"))(26)
# z <- zClust(x=aa_auroc_chunk, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
heatmap.2(aa_auroc_chunk_nodist, trace='none', col=(cols),
          breaks = seq(-0.13, 0.13, 0.01), margins = c(6,10)
          #, Rowv = aam1$rowInd
          , Colv = rev(aam1$colInd)
          , revC = T
)

#aa_auprc_rand_nodist[is.na(aa_auprc_rand_nodist)] <- 0
range(aa_auprc_rand_nodist)

cols <- colorRampPalette(brewer.pal(13, "YlOrRd"))(13)
# z <- zClust(x=aa_auprc_rand, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
aam <- heatmap.2(aa_auprc_rand_nodist, trace='none', col=(cols), breaks = seq(-0.00, 0.13, 0.01), margins = c(6,10))
aam
#aa_auprc_chunk_nodist[is.na(aa_auprc_chunk_nodist)] <- 0
range(aa_auprc_chunk_nodist)

cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(13)
# z <- zClust(x=aa_auprc_chunk, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
heatmap.2(aa_auprc_chunk_nodist, trace='none', col=(cols), breaks = seq(-0.00, 0.13, 0.01), margins = c(6,10), Rowv = aam$rowInd)




################################################################################################
################################################################################################
# comprehensive evaluation of feature importances --> only using models in the top 10 percentile for each lncRNA

aa_all_impor_list_cv2TX_nodist <- list()
for(i in 1:length(aa_all_impor_list_cv2TX)){
  aa_all_impor_list_cv2TX_nodist[[i]] <- list()
  for(j in 1:length(aa_all_impor_list_cv2TX[[i]])){
    aa_all_impor_list_cv2TX_nodist[[i]][[j]] <- aa_all_impor_list_cv2TX[[i]][[j]][-grep(pattern = "dist",
                                                                                        x = names(aa_all_impor_list_cv2TX[[i]][[j]]))]
  }
  names(aa_all_impor_list_cv2TX_nodist[[i]]) <- names(aa_all_impor_list_cv2TX[[i]])
}
names(aa_all_impor_list_cv2TX_nodist) <- names(aa_all_impor_list_cv2TX)




#aa_all_auroc_list_noDist


aac <- max(unlist(lapply(aa_all_auprc_list_noDist, nrow)))

aa_all_auroc_rand_cv2TX_noDist <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auroc_rand_cv2TX_noDist) <- rownames(aa_all_auprc_list_noDist[[2]])
rownames(aa_all_auroc_rand_cv2TX_noDist) <- aafn

aa_all_auroc_chunk_cv2TX_noDist <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auroc_chunk_cv2TX_noDist) <- rownames(aa_all_auprc_list_noDist[[2]])
rownames(aa_all_auroc_chunk_cv2TX_noDist) <- aafn

aa_all_auprc_rand_cv2TX_noDist <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auprc_rand_cv2TX_noDist) <- rownames(aa_all_auprc_list_noDist[[2]])
rownames(aa_all_auprc_rand_cv2TX_noDist) <- aafn

aa_all_auprc_chunk_cv2TX_noDist <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auprc_chunk_cv2TX_noDist) <- rownames(aa_all_auprc_list_noDist[[2]])
rownames(aa_all_auprc_chunk_cv2TX_noDist) <- aafn


aa_all_auprc_list_cv2_diffTX_NoDist <- list()
for(i in 1:length(aa_all_auprc_list_noDist)){
  aa_all_auprc_list_cv2_diffTX_NoDist[[i]] <- aa_all_auprc_list_noDist[[i]]
  for(j in 1:nrow(aa_all_auprc_list_noDist[[i]])){
    aa_all_auprc_list_cv2_diffTX_NoDist[[i]][j,] <- aa_all_auprc_list_cv2_diffTX_NoDist[[i]][j,] - aa_auprc_base[i,]
  }
}
for(i in 1:length(aa_all_auroc_list_noDist)){
  aatmmpccv_roc <- rowMeans(aa_all_auroc_list_noDist[[i]][,c(1:5)], na.rm = T)
  aatmmprcv_roc <- rowMeans(aa_all_auroc_list_noDist[[i]][,c(6:10)], na.rm = T)
  
  aatmmpccv_prc <- rowMeans(aa_all_auprc_list_cv2_diffTX_NoDist[[i]][,c(1:5)], na.rm = T)
  aatmmprcv_prc <- rowMeans(aa_all_auprc_list_cv2_diffTX_NoDist[[i]][,c(6:10)], na.rm = T)
  
  aacolord <- match(rownames(aa_all_auroc_list_noDist[[i]]), 
                    colnames(aa_all_auroc_chunk_cv2TX_noDist))
  
  aa_all_auroc_chunk_cv2TX_noDist[i,aacolord] <- aatmmpccv_roc
  aa_all_auroc_rand_cv2TX_noDist[i,aacolord] <- aatmmprcv_roc
  aa_all_auprc_chunk_cv2TX_noDist[i,aacolord] <- aatmmpccv_prc
  aa_all_auprc_rand_cv2TX_noDist[i,aacolord] <- aatmmprcv_prc
}


png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/AUROC_rand_heatmap_noDist.png", 
    width = 17*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
range(aa_all_auroc_rand_cv2TX_noDist, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(25)
heatmap.2((aa_all_auroc_rand_cv2TX_noDist), trace='none', col=(cols), breaks = seq(0.5, 1, 0.02), margins = c(12,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/AUROC_chunk_heatmap_noDist.png", 
    width = 17*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
range(aa_all_auroc_chunk_cv2TX_noDist, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(25)
heatmap.2((aa_all_auroc_chunk_cv2TX_noDist), trace='none', col=(cols), breaks = seq(0.5, 1, 0.02), margins = c(12,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/AUPRC_rand_heatmap_noDist.png", 
    width = 17*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
range(aa_all_auprc_rand_cv2TX_noDist, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(20)
heatmap.2((aa_all_auprc_rand_cv2TX_noDist), trace='none', col=(cols), breaks = seq(0.0, 0.4, 0.02), margins = c(12,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/AUPRC_chunk_heatmap_noDist.png", 
    width = 17*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
range(aa_all_auprc_chunk_cv2TX_noDist, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(20)
heatmap.2((aa_all_auprc_chunk_cv2TX_noDist), trace='none', col=(cols), breaks = seq(0.0, 0.4, 0.02), margins = c(12,10))
dev.off()

aatop10perlnc_rand_sequ <- list()
aatop10perlnc_rand_cont <- list()

aatop10perlnc_chun_sequ <- list()
aatop10perlnc_chun_cont <- list()

aagghh <- grep("sequ", colnames(aa_all_auroc_rand_cv2TX_noDist))

aa_all_auroc_rand_cv2TX_noDist_sequ <- aa_all_auroc_rand_cv2TX_noDist[,-aagghh]
aa_all_auroc_rand_cv2TX_noDist_cont <- aa_all_auroc_rand_cv2TX_noDist[,aagghh]

aaxxt <- aa_all_auprc_rand_cv2TX_noDist_sequ
rownames(aaxxt)[sort(apply(aaxxt,1,max), decreasing = T, index.return=T)$ix]

aa_all_auprc_rand_cv2TX_noDist_sequ <- aa_all_auprc_rand_cv2TX_noDist[,-aagghh]
aa_all_auprc_rand_cv2TX_noDist_cont <- aa_all_auprc_rand_cv2TX_noDist[,aagghh]

aa_all_auroc_chun_cv2TX_noDist_sequ <- aa_all_auroc_chunk_cv2TX_noDist[,-aagghh]
aa_all_auroc_chun_cv2TX_noDist_cont <- aa_all_auroc_chunk_cv2TX_noDist[,aagghh]

aa_all_auprc_chun_cv2TX_noDist_sequ <- aa_all_auprc_chunk_cv2TX_noDist[,-aagghh]
aa_all_auprc_chun_cv2TX_noDist_cont <- aa_all_auprc_chunk_cv2TX_noDist[,aagghh]


aagghh1 <- grep("_TXRB", colnames(aa_all_auprc_rand_cv2TX_noDist_sequ))
aagghh2 <- grep("_TXRE", colnames(aa_all_auprc_rand_cv2TX_noDist_sequ))
aa_all_auroc_rand_cv2TX_noDist_sequ_noTX <- aa_all_auroc_rand_cv2TX_noDist_sequ[,-union(aagghh1,aagghh2)]
aa_all_auprc_rand_cv2TX_noDist_sequ_noTX <- aa_all_auprc_rand_cv2TX_noDist_sequ[,-union(aagghh1,aagghh2)]
aa_all_auroc_chun_cv2TX_noDist_sequ_noTX <- aa_all_auroc_chun_cv2TX_noDist_sequ[,-union(aagghh1,aagghh2)]
aa_all_auprc_chun_cv2TX_noDist_sequ_noTX <- aa_all_auprc_chun_cv2TX_noDist_sequ[,-union(aagghh1,aagghh2)]

for(i in 1:nrow(aa_all_auprc_chunk_cv2TX_noDist)){
  
  aamyind1RS <- quantile(aa_all_auroc_rand_cv2TX_noDist_sequ[i, ],  seq(0,1,0.1))
  aamyind2RS <- quantile(aa_all_auprc_rand_cv2TX_noDist_sequ[i, ],  seq(0,1,0.1))
  aamyind1CS <- quantile(aa_all_auroc_chun_cv2TX_noDist_sequ[i, ],  seq(0,1,0.1))
  aamyind2CS <- quantile(aa_all_auprc_chun_cv2TX_noDist_sequ[i, ],  seq(0,1,0.1))
  
  aamyind1RC <- quantile(aa_all_auroc_rand_cv2TX_noDist_cont[i, ],  seq(0,1,0.1))
  aamyind2RC <- quantile(aa_all_auprc_rand_cv2TX_noDist_cont[i, ],  seq(0,1,0.1))
  aamyind1CC <- quantile(aa_all_auroc_chun_cv2TX_noDist_cont[i, ],  seq(0,1,0.1))
  aamyind2CC <- quantile(aa_all_auprc_chun_cv2TX_noDist_cont[i, ],  seq(0,1,0.1))
  

  aatop10perlnc_rand_sequ[[i]] <- union(colnames(aa_all_auroc_rand_cv2TX_noDist_sequ)[aa_all_auroc_rand_cv2TX_noDist_sequ[i, ] >= aamyind1RS[10]],
                                   colnames(aa_all_auprc_rand_cv2TX_noDist_sequ)[aa_all_auprc_rand_cv2TX_noDist_sequ[i, ] >= aamyind2RS[10]])
  aatop10perlnc_chun_sequ[[i]] <- union(colnames(aa_all_auroc_chun_cv2TX_noDist_sequ)[aa_all_auroc_chun_cv2TX_noDist_sequ[i, ] >= aamyind1CS[10]],
                                   colnames(aa_all_auprc_chun_cv2TX_noDist_sequ)[aa_all_auprc_chun_cv2TX_noDist_sequ[i, ] >=  aamyind2CS[10]])
  
  aatop10perlnc_rand_cont[[i]] <- union(colnames(aa_all_auroc_rand_cv2TX_noDist_cont)[aa_all_auroc_rand_cv2TX_noDist_cont[i, ] >= aamyind1RC[10]],
                                   colnames(aa_all_auprc_rand_cv2TX_noDist_cont)[aa_all_auprc_rand_cv2TX_noDist_cont[i, ] >= aamyind2RC[10]])
  aatop10perlnc_chun_cont[[i]] <- union(colnames(aa_all_auroc_chun_cv2TX_noDist_cont)[aa_all_auroc_chun_cv2TX_noDist_cont[i, ] >= aamyind1CC[10]],
                                   colnames(aa_all_auprc_chun_cv2TX_noDist_cont)[aa_all_auprc_chun_cv2TX_noDist_cont[i, ] >=  aamyind2CC[10]])
  
}
names(aatop10perlnc_rand_sequ) <- rownames(aa_all_auprc_chunk_cv2TX_noDist)
names(aatop10perlnc_chun_sequ) <- rownames(aa_all_auprc_chunk_cv2TX_noDist)
names(aatop10perlnc_rand_cont) <- rownames(aa_all_auprc_chunk_cv2TX_noDist)
names(aatop10perlnc_chun_cont) <- rownames(aa_all_auprc_chunk_cv2TX_noDist)

#make a dataframe with the following columns: lncRNA, feature, model, feature_ranking
aa_all_feat <- c(my_name_dic[,1], paste0("sequ_TR_", c(1:73)))


aaimpdf_rand_s <- list()
aaimpdf_rand_c <- list()
aaimpdf_chun_s <- list()
aaimpdf_chun_c <- list()

for(i in 1:length(aatop10perlnc_rand_sequ)){
  aaimpdf_rand_s[[i]] <- list()
  aaimpdf_rand_c[[i]] <- list()
  aaimpdf_chun_s[[i]] <- list()
  aaimpdf_chun_c[[i]] <- list()
  aacur_feat_c <- aa_all_impor_list_cv2TX_nodist[[i]][1:5]
  aacur_feat_r <- aa_all_impor_list_cv2TX_nodist[[i]][6:10]
  
  for(j in 1:length(aatop10perlnc_rand_sequ[[i]])){
    aatmplist <- list()
    for(k in 1:length(aacur_feat_r)){
      aayytmp <- aacur_feat_r[[k]][[match(aatop10perlnc_rand_sequ[[i]][j], 
                                          names( aacur_feat_r[[k]]))]]

      aatmplist[[k]] <- array(dim = length(aa_all_feat))
      names(aatmplist[[k]]) <- aa_all_feat
      aatmplist[[k]][match(names(aayytmp), names(aatmplist[[k]]))] <- aayytmp
      
      # if(k > 1){
      #   stopifnot(all(names(aatmplist[[k]]) == names(aatmplist[[k-1]])))
      # }
    }
    #stopifnot(length(unique(unlist(lapply(aatmplist, length)))) == 1)
    aatmplist_df <- do.call(rbind, aatmplist)
    aatmplist_df_mean <- colMeans(aatmplist_df, na.rm = F)
    aasortind <- sort(aatmplist_df_mean, decreasing = T, index.return = T, na.last =T)$ix
    aasortind2 <- c(1:length(aasortind))[sort(aasortind, index.return = T)$ix]
    aasortind2[is.na(aatmplist_df_mean)] <- NA
    aaimpdf_rand_s[[i]][[j]] <- data.frame(lncRNA = rep(names(aatop10perlnc_rand_sequ)[i], length(aatmplist_df_mean)),
                                           model =  rep(aatop10perlnc_rand_sequ[[i]][j], length(aatmplist_df_mean)),
                                           feature = names(aatmplist_df_mean),
                                           rank = aasortind2)
  }
  
  for(j in 1:length(aatop10perlnc_rand_cont[[i]])){
    aatmplist <- list()
    for(k in 1:length(aacur_feat_r)){
      aayytmp <- aacur_feat_r[[k]][[match(aatop10perlnc_rand_cont[[i]][j], 
                                          names( aacur_feat_r[[k]]))]]
      
      aatmplist[[k]] <- array(dim = length(aa_all_feat))
      names(aatmplist[[k]]) <- aa_all_feat
      aatmplist[[k]][match(names(aayytmp), names(aatmplist[[k]]))] <- aayytmp
      
      # if(k > 1){
      #   stopifnot(all(names(aatmplist[[k]]) == names(aatmplist[[k-1]])))
      # }
    }
    #stopifnot(length(unique(unlist(lapply(aatmplist, length)))) == 1)
    aatmplist_df <- do.call(rbind, aatmplist)
    aatmplist_df_mean <- colMeans(aatmplist_df, na.rm = F)
    aasortind <- sort(aatmplist_df_mean, decreasing = T, index.return = T, na.last =T)$ix
    aasortind2 <- c(1:length(aasortind))[sort(aasortind, index.return = T)$ix]
    aasortind2[is.na(aatmplist_df_mean)] <- NA
    aaimpdf_rand_c[[i]][[j]] <- data.frame(lncRNA = rep(names(aatop10perlnc_rand_cont)[i], length(aatmplist_df_mean)),
                                           model =  rep(aatop10perlnc_rand_cont[[i]][j], length(aatmplist_df_mean)),
                                           feature = names(aatmplist_df_mean),
                                           rank = aasortind2)
  }
  
  
  for(j in 1:length(aatop10perlnc_chun_sequ[[i]])){
    aatmplist <- list()
    for(k in 1:length(aacur_feat_c)){
      aayytmp <- aacur_feat_c[[k]][[match(aatop10perlnc_chun_sequ[[i]][j], 
                                          names( aacur_feat_c[[k]]))]]
      
      aatmplist[[k]] <- array(dim = length(aa_all_feat))
      names(aatmplist[[k]]) <- aa_all_feat
      aatmplist[[k]][match(names(aayytmp), names(aatmplist[[k]]))] <- aayytmp
      
      # if(k > 1){
      #   stopifnot(all(names(aatmplist[[k]]) == names(aatmplist[[k-1]])))
      # }
    }
    #stopifnot(length(unique(unlist(lapply(aatmplist, length)))) == 1)
    aatmplist_df <- do.call(rbind, aatmplist)
    aatmplist_df_mean <- colMeans(aatmplist_df, na.rm = F)
    aasortind <- sort(aatmplist_df_mean, decreasing = T, index.return = T, na.last =T)$ix
    aasortind2 <- c(1:length(aasortind))[sort(aasortind, index.return = T)$ix]
    aasortind2[is.na(aatmplist_df_mean)] <- NA
    aaimpdf_chun_s[[i]][[j]] <- data.frame(lncRNA = rep(names(aatop10perlnc_chun_sequ)[i], length(aatmplist_df_mean)),
                                           model =  rep(aatop10perlnc_chun_sequ[[i]][j], length(aatmplist_df_mean)),
                                           feature = names(aatmplist_df_mean),
                                           rank = aasortind2)
  }
  
  
  for(j in 1:length(aatop10perlnc_chun_cont[[i]])){
    aatmplist <- list()
    for(k in 1:length(aacur_feat_c)){
      aayytmp <- aacur_feat_c[[k]][[match(aatop10perlnc_chun_cont[[i]][j], 
                                          names( aacur_feat_c[[k]]))]]
      
      aatmplist[[k]] <- array(dim = length(aa_all_feat))
      names(aatmplist[[k]]) <- aa_all_feat
      aatmplist[[k]][match(names(aayytmp), names(aatmplist[[k]]))] <- aayytmp
      
      # if(k > 1){
      #   stopifnot(all(names(aatmplist[[k]]) == names(aatmplist[[k-1]])))
      # }
    }
    #stopifnot(length(unique(unlist(lapply(aatmplist, length)))) == 1)
    aatmplist_df <- do.call(rbind, aatmplist)
    aatmplist_df_mean <- colMeans(aatmplist_df, na.rm = F)
    aasortind <- sort(aatmplist_df_mean, decreasing = T, index.return = T, na.last =T)$ix
    aasortind2 <- c(1:length(aasortind))[sort(aasortind, index.return = T)$ix]
    aasortind2[is.na(aatmplist_df_mean)] <- NA
    aaimpdf_chun_c[[i]][[j]] <- data.frame(lncRNA = rep(names(aatop10perlnc_chun_cont)[i], length(aatmplist_df_mean)),
                                           model =  rep(aatop10perlnc_chun_cont[[i]][j], length(aatmplist_df_mean)),
                                           feature = names(aatmplist_df_mean),
                                           rank = aasortind2)
  }
  # aatmp1
  # aatmp2
  # aatmp3
  # aatmp4
  
}

aalncFeat_seq_R_list <- list()
aalncFeat_seq_C_list <- list()

aalncFeat_con_R_list <- list()
aalncFeat_con_C_list <- list()

for(i in 1:length(aaimpdf_chun_c)){
  aalncFeat_seq_R_list[[i]] <- do.call(rbind, aaimpdf_rand_s[[i]])
  aalncFeat_seq_C_list[[i]] <- do.call(rbind, aaimpdf_chun_s[[i]])
  aalncFeat_con_R_list[[i]] <- do.call(rbind, aaimpdf_rand_c[[i]])
  aalncFeat_con_C_list[[i]] <- do.call(rbind, aaimpdf_chun_c[[i]])
}

aalncFeat_seq_R_df <- do.call(rbind, aalncFeat_seq_R_list)
aalncFeat_seq_C_df <- do.call(rbind, aalncFeat_seq_C_list)
aalncFeat_con_R_df <- do.call(rbind, aalncFeat_con_R_list)
aalncFeat_con_C_df <- do.call(rbind, aalncFeat_con_C_list)



lncRNA_feature_seq_rand_summary_list <- list()
for(i in 1:length(aalncFeat_seq_R_list)){
  print(i)
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_violin/Random/Seq/rand_seq_",names(aatop10perlnc_rand_sequ)[i], ".png"),    # create PNG for the heat map        
      width = 12*300,        # 5 x 300 pixels
      height = 8*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  
  aatgc <- summarySE(aalncFeat_seq_R_list[[i]], measurevar=c("rank"), groupvars=c("feature"), na.rm=T)
  lncRNA_feature_seq_rand_summary_list[[i]] <- aatgc
  aatgc <- aatgc[aatgc$N >0 ,]
  aatgc2 <- aatgc[sort(aatgc$median, index.return=T)$ix,]
  aaplx <- aalncFeat_seq_R_list[[i]]
  aaplx$feature <- factor(aalncFeat_seq_R_list[[i]]$feature, levels = aatgc2$feature)
  aaxplt <- ggplot(aaplx[aalncFeat_seq_R_list[[i]]$feature %in% aatgc2$feature[1:100],], aes(x=feature, y=rank)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),scale = "count") +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(paste0("ranking of sequence features for ",names(aatop10perlnc_rand_sequ)[i]))+
    ylab("rank")+
    xlab("seq feature") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}

lncRNA_feature_seq_chunk_summary_list <- list()
for(i in 1:length(aalncFeat_seq_C_list)){
  print(i)
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_violin/Chunk/Seq/chunk_seq_",
                        names(aatop10perlnc_chun_sequ)[i], ".png"),    # create PNG for the heat map        
      width = 12*300,        # 5 x 300 pixels
      height = 8*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  
  aatgc <- summarySE(aalncFeat_seq_C_list[[i]], measurevar=c("rank"), groupvars=c("feature"), na.rm=T)
  lncRNA_feature_seq_chunk_summary_list[[i]] <- aatgc
  aatgc <- aatgc[aatgc$N >0 ,]
  aatgc2 <- aatgc[sort(aatgc$median, index.return=T)$ix,]
  aaplx <- aalncFeat_seq_C_list[[i]]
  aaplx$feature <- factor(aalncFeat_seq_C_list[[i]]$feature, levels = aatgc2$feature)
  aaxplt <- ggplot(aaplx[aalncFeat_seq_C_list[[i]]$feature %in% aatgc2$feature[1:100],], aes(x=feature, y=rank)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),scale = "count") +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(paste0("ranking of sequence features for ",names(aatop10perlnc_chun_sequ)[i]))+
    ylab("rank")+
    xlab("seq feature") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}


lncRNA_feature_context_rand_summary_list <- list()
for(i in 1:length(aalncFeat_con_R_list)){
  print(i)
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_violin/Random/Context/rand_context_",names(aatop10perlnc_rand_sequ)[i], ".png"),    # create PNG for the heat map        
      width = 18*300,        # 5 x 300 pixels
      height = 8*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  
  aatgc <- summarySE(aalncFeat_con_R_list[[i]], measurevar=c("rank"), groupvars=c("feature"), na.rm=T)
  lncRNA_feature_context_rand_summary_list[[i]] <- aatgc
  aatgc <- aatgc[aatgc$N >0 ,]
  aatgc2 <- aatgc[sort(aatgc$median, index.return=T)$ix,]
  aaplx <- aalncFeat_con_R_list[[i]]
  aaplx$feature <- factor(aalncFeat_con_R_list[[i]]$feature, levels = aatgc2$feature)
  aaxplt <- ggplot(aaplx[aalncFeat_con_R_list[[i]]$feature %in% aatgc2$feature[1:173],], aes(x=feature, y=rank)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),scale = "count") +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(paste0("ranking of context features for ",names(aatop10perlnc_rand_sequ)[i]))+
    ylab("rank")+
    xlab("context feature") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}

lncRNA_feature_context_chunk_summary_list <- list()
for(i in 1:length(aalncFeat_con_C_list)){
  print(i)
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_violin/Chunk/Context/chunk_context_",
                        names(aatop10perlnc_chun_sequ)[i], ".png"),    # create PNG for the heat map        
      width = 18*300,        # 5 x 300 pixels
      height = 8*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  
  aatgc <- summarySE(aalncFeat_con_C_list[[i]], measurevar=c("rank"), groupvars=c("feature"), na.rm=T)
  lncRNA_feature_context_chunk_summary_list[[i]] <- aatgc
  aatgc <- aatgc[aatgc$N >0 ,]
  aatgc2 <- aatgc[sort(aatgc$median, index.return=T)$ix,]
  aaplx <- aalncFeat_con_C_list[[i]]
  aaplx$feature <- factor(aalncFeat_con_C_list[[i]]$feature, levels = aatgc2$feature)
  aaxplt <- ggplot(aaplx[aalncFeat_con_C_list[[i]]$feature %in% aatgc2$feature[1:173],], aes(x=feature, y=rank)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),scale = "count") +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(paste0("ranking of context features for ",names(aatop10perlnc_chun_sequ)[i]))+
    ylab("rank")+
    xlab("context feature") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}

#########
# make boxplot of the same plots
#lncRNA_feature_seq_rand_summary_list <- list()
for(i in 1:length(aalncFeat_seq_R_list)){
  print(i)
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_violin/Random/Seq/rand_seq_",
                        names(aatop10perlnc_rand_sequ)[i], "_boxplot.png"),    # create PNG for the heat map        
      width = 12*300,        # 5 x 300 pixels
      height = 8*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  
  aatgc <- summarySE(aalncFeat_seq_R_list[[i]], measurevar=c("rank"), groupvars=c("feature"), na.rm=T)
  #lncRNA_feature_seq_rand_summary_list[[i]] <- aatgc
  aatgc <- aatgc[aatgc$N >0 ,]
  aatgc2 <- aatgc[sort(aatgc$median, index.return=T)$ix,]
  aaplx <- aalncFeat_seq_R_list[[i]]
  aaplx$feature <- factor(aalncFeat_seq_R_list[[i]]$feature, levels = aatgc2$feature)
  aaxplt <- ggplot(aaplx[aalncFeat_seq_R_list[[i]]$feature %in% aatgc2$feature[1:100],], aes(x=feature, y=rank)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(paste0("ranking of sequence features for ",names(aatop10perlnc_rand_sequ)[i]))+
    ylab("rank")+
    xlab("seq feature") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}

#lncRNA_feature_seq_chunk_summary_list <- list()
for(i in 1:length(aalncFeat_seq_C_list)){
  print(i)
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_violin/Chunk/Seq/chunk_seq_",
                        names(aatop10perlnc_chun_sequ)[i], "_boxplot.png"),    # create PNG for the heat map        
      width = 12*300,        # 5 x 300 pixels
      height = 8*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  
  aatgc <- summarySE(aalncFeat_seq_C_list[[i]], measurevar=c("rank"), groupvars=c("feature"), na.rm=T)
  #lncRNA_feature_seq_chunk_summary_list[[i]] <- aatgc
  aatgc <- aatgc[aatgc$N >0 ,]
  aatgc2 <- aatgc[sort(aatgc$median, index.return=T)$ix,]
  aaplx <- aalncFeat_seq_C_list[[i]]
  aaplx$feature <- factor(aalncFeat_seq_C_list[[i]]$feature, levels = aatgc2$feature)
  aaxplt <- ggplot(aaplx[aalncFeat_seq_C_list[[i]]$feature %in% aatgc2$feature[1:100],], aes(x=feature, y=rank)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(paste0("ranking of sequence features for ",names(aatop10perlnc_chun_sequ)[i]))+
    ylab("rank")+
    xlab("seq feature") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}


#lncRNA_feature_context_rand_summary_list <- list()
for(i in 1:length(aalncFeat_con_R_list)){
  print(i)
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_violin/Random/Context/rand_context_",
                        names(aatop10perlnc_rand_sequ)[i], "_boxplot.png"),    # create PNG for the heat map        
      width = 18*300,        # 5 x 300 pixels
      height = 8*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  
  aatgc <- summarySE(aalncFeat_con_R_list[[i]], measurevar=c("rank"), groupvars=c("feature"), na.rm=T)
 # lncRNA_feature_context_rand_summary_list[[i]] <- aatgc
  aatgc <- aatgc[aatgc$N >0 ,]
  aatgc2 <- aatgc[sort(aatgc$median, index.return=T)$ix,]
  aaplx <- aalncFeat_con_R_list[[i]]
  aaplx$feature <- factor(aalncFeat_con_R_list[[i]]$feature, levels = aatgc2$feature)
  aaxplt <- ggplot(aaplx[aalncFeat_con_R_list[[i]]$feature %in% aatgc2$feature[1:173],], aes(x=feature, y=rank)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(paste0("ranking of context features for ",names(aatop10perlnc_rand_sequ)[i]))+
    ylab("rank")+
    xlab("context feature") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}

#lncRNA_feature_context_chunk_summary_list <- list()
for(i in 1:length(aalncFeat_con_C_list)){
  print(i)
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_violin/Chunk/Context/chunk_context_",
                        names(aatop10perlnc_chun_sequ)[i], "_boxplot.png"),    # create PNG for the heat map        
      width = 18*300,        # 5 x 300 pixels
      height = 8*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  
  aatgc <- summarySE(aalncFeat_con_C_list[[i]], measurevar=c("rank"), groupvars=c("feature"), na.rm=T)
#  lncRNA_feature_context_chunk_summary_list[[i]] <- aatgc
  aatgc <- aatgc[aatgc$N >0 ,]
  aatgc2 <- aatgc[sort(aatgc$median, index.return=T)$ix,]
  aaplx <- aalncFeat_con_C_list[[i]]
  aaplx$feature <- factor(aalncFeat_con_C_list[[i]]$feature, levels = aatgc2$feature)
  aaxplt <- ggplot(aaplx[aalncFeat_con_C_list[[i]]$feature %in% aatgc2$feature[1:173],], aes(x=feature, y=rank)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(paste0("ranking of context features for ",names(aatop10perlnc_chun_sequ)[i]))+
    ylab("rank")+
    xlab("context feature") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}

##########################################################################################
##############################################################################
# draw heatmaps using the median value



lncRNA_feature_seq_rand_summary_list
lncRNA_feature_seq_chunk_summary_list
lncRNA_feature_context_rand_summary_list
lncRNA_feature_context_chunk_summary_list


#unlist(lapply(lncRNA_feature_context_rand_summary_list, nrow))
all(lncRNA_feature_context_rand_summary_list[[1]]$feature == lncRNA_feature_context_rand_summary_list[[20]]$feature)



aaxlsi <- lapply(lncRNA_feature_context_rand_summary_list, "[[", "median")
aaxlsi <- do.call(rbind, aaxlsi)
colnames(aaxlsi) <- lncRNA_feature_context_rand_summary_list[[1]]$feature
rownames(aaxlsi) <- rownames(aa_all_auprc_chunk_cv2TX_noDist)

aaxlsi <- aaxlsi[,!(colSums(is.na(aaxlsi)) == nrow(aaxlsi))]
aamin <- apply(aaxlsi, 2, min, na.rm=T)
aaxlsi <- aaxlsi[, aamin <= 100]
aaarr <- range(aaxlsi, na.rm = T)
aabreak <- c(0,5,10,20,50,75,100,aaarr[2])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(length(aabreak)-1)
png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_heatmap_Context_rand.png",    # create PNG for the heat map        
    width = 30*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11)
heatmap.2((aaxlsi), trace='none', col=rev(cols), breaks = aabreak, margins = c(12,8))
dev.off()

aaxlsi <- aaxlsi[, -grep("sequ_TR_", colnames(aaxlsi))]
png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_heatmap_Context_rand_noSEQU.png",    # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11)
heatmap.2((aaxlsi), trace='none', col=rev(cols), breaks = aabreak, margins = c(12,8))
dev.off()



#################

aaxlsi <- lapply(lncRNA_feature_context_chunk_summary_list  , "[[", "median")
aaxlsi <- do.call(rbind, aaxlsi)
colnames(aaxlsi) <- lncRNA_feature_context_chunk_summary_list[[1]]$feature
rownames(aaxlsi) <- rownames(aa_all_auprc_chunk_cv2TX_noDist)

aaxlsi <- aaxlsi[,!(colSums(is.na(aaxlsi)) == nrow(aaxlsi))]
aamin <- apply(aaxlsi, 2, min, na.rm=T)
aaxlsi <- aaxlsi[, aamin <= 100]
aaarr <- range(aaxlsi, na.rm = T)
aabreak <- c(0,5,10,20,50,75,100,aaarr[2])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(length(aabreak)-1)
png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_heatmap_Context_chunk.png",    # create PNG for the heat map        
    width = 30*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11)
heatmap.2((aaxlsi), trace='none', col=rev(cols), breaks = aabreak, margins = c(12,8))
dev.off()

aaxlsi <- aaxlsi[, -grep("sequ_TR_", colnames(aaxlsi))]
png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_heatmap_Context_chunk_noSEQU.png",    # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11)
heatmap.2((aaxlsi), trace='none', col=rev(cols), breaks = aabreak, margins = c(12,8))
dev.off()



#################3


aaxlsi <- lapply(lncRNA_feature_seq_rand_summary_list  , "[[", "median")
aaxlsi <- do.call(rbind, aaxlsi)
colnames(aaxlsi) <- lncRNA_feature_seq_rand_summary_list[[1]]$feature
rownames(aaxlsi) <- rownames(aa_all_auprc_chunk_cv2TX_noDist)

aaxlsi <- aaxlsi[,!(colSums(is.na(aaxlsi)) == nrow(aaxlsi))]
aamin <- apply(aaxlsi, 2, min, na.rm=T)
aaxlsi <- aaxlsi[, aamin <= 50]
aaarr <- range(aaxlsi, na.rm = T)
aabreak <- c(0,5,10,20,50,75,100,aaarr[2])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(length(aabreak)-1)
png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_heatmap_Seq_rand.png",    # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11)
heatmap.2((aaxlsi), trace='none', col=rev(cols), breaks = aabreak, margins = c(12,8))
dev.off()

#################


aaxlsi <- lapply(lncRNA_feature_seq_chunk_summary_list  , "[[", "median")
aaxlsi <- do.call(rbind, aaxlsi)
colnames(aaxlsi) <- lncRNA_feature_seq_chunk_summary_list[[1]]$feature
rownames(aaxlsi) <- rownames(aa_all_auprc_chunk_cv2TX_noDist)

aaxlsi <- aaxlsi[,!(colSums(is.na(aaxlsi)) == nrow(aaxlsi))]
aamin <- apply(aaxlsi, 2, min, na.rm=T)
aaxlsi <- aaxlsi[, aamin <= 100]
aaarr <- range(aaxlsi, na.rm = T)
aabreak <- c(0,5,10,20,50,75,100,aaarr[2])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(length(aabreak)-1)
png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_heatmap_Seq_chunk.png",    # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11)
heatmap.2((aaxlsi), trace='none', col=rev(cols), breaks = aabreak, margins = c(12,8))
dev.off()

#######################################################################################################################
#######################################################################################################################
# for each lncRNA 
# do two 100 boxplots of the top 100 sequence or context feaures: full dataset, positive against negative
# also note down the positive tiles with values in top 10 percentile for the feature

aalnc_feature_top_pos_tile_list <- list()
aaln <- unique(aatmp_lnc_df_cv2TX$lncRNA)

lncRNA_feature_seq_chunk_summary_list
lncRNA_feature_seq_rand_summary_list
lncRNA_feature_context_chunk_summary_list
lncRNA_feature_context_rand_summary_list


aafile <- list.files("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/")
for(i in 1:length(aafile)){
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/", aafile[i]))
  aa_myf1 <- lncRNA_feature_seq_chunk_summary_list[[i]][lncRNA_feature_seq_chunk_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("sequ_TR_", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_boxplot/fb_",aaln[i],"_Seq_chunk.png"),    # create PNG for the heat map        
      width = 16*300,        # 5 x 300 pixels
      height = 16*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(2,2,2,2))
  for(j in 1:length(aa_myf11)){
    boxplot(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])]~my_Dataset$label, main = aa_myf11[j])
  }
  
  aa_myf1 <- lncRNA_feature_seq_rand_summary_list[[i]][lncRNA_feature_seq_rand_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("sequ_TR_", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  dev.off()
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_boxplot/fb_",aaln[i],"_Seq_rand.png"),    # create PNG for the heat map        
      width = 16*300,        # 5 x 300 pixels
      height = 16*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(2,2,2,2))
  for(j in 1:length(aa_myf11)){
    boxplot(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])]~my_Dataset$label, main = aa_myf11[j])
  }
  dev.off()
  aa_myf1 <- lncRNA_feature_context_rand_summary_list[[i]][lncRNA_feature_context_rand_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("sequ_TR_", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_boxplot/fb_",aaln[i],"_Context_rand.png"),    # create PNG for the heat map        
      width = 16*300,        # 5 x 300 pixels
      height = 16*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(2,2,2,2))
  for(j in 1:length(aa_myf11)){
    boxplot(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])]~my_Dataset$label, main = aa_myf11[j])
  }
  dev.off()
  aa_myf1 <- lncRNA_feature_context_chunk_summary_list[[i]][lncRNA_feature_context_chunk_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("sequ_TR_", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_boxplot/fb_",aaln[i],"_Context_chunk.png"),    # create PNG for the heat map        
      width = 16*300,        # 5 x 300 pixels
      height = 16*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(2,2,2,2))
  for(j in 1:length(aa_myf11)){
    boxplot(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])]~my_Dataset$label, main = aa_myf11[j])
  }
  dev.off()
  
  
}

######## plot overlapping histograms instead of the above boxplots

aafile <- list.files("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/")

c1 <- rgb(0,0,255,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,0,0, max = 255, alpha = 80, names = "lt.red")
length(aafile)
for(i in 1:1){
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/", aafile[i]))
  aa_myf1 <- lncRNA_feature_seq_chunk_summary_list[[i]][lncRNA_feature_seq_chunk_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("sequ_TR_", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_Hist_ovl/fh_",aaln[i],"_Seq_chunk.png"),    # create PNG for the heat map        
      width = 16*300,        # 5 x 300 pixels
      height = 16*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(2,2,2,2))
  for(j in 1:length(aa_myf11)){
    aa_myneg <- my_Dataset[my_Dataset$label == "neg", match(aa_myf11[j], my_name_dic[,1])]
    aa_mypos <- my_Dataset[my_Dataset$label == "pos", match(aa_myf11[j], my_name_dic[,1])]
    aa_f_range <- range(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])], na.rm = T)
    ax <- pretty(floor(aa_f_range[1]):ceiling(aa_f_range[2]), n = 20)
    histneg <- hist(aa_myneg,  breaks = ax , plot = F)
    histpos <- hist(aa_mypos,  breaks = ax , plot = F)
    plot(histneg, col = c1, freq = F, main = aa_myf11[j]) # Plot 1st histogram using a transparent color
    plot(histpos, col = c2, add = TRUE, freq = F) # Add 2nd histogram using different color
    #boxplot(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])]~my_Dataset$label, main = aa_myf11[j])
  }
  
  aa_myf1 <- lncRNA_feature_seq_rand_summary_list[[i]][lncRNA_feature_seq_rand_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("sequ_TR_", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  dev.off()
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_Hist_ovl/fh_",aaln[i],"_Seq_rand.png"),    # create PNG for the heat map        
      width = 16*300,        # 5 x 300 pixels
      height = 16*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(2,2,2,2))
  for(j in 1:length(aa_myf11)){
    aa_myneg <- my_Dataset[my_Dataset$label == "neg", match(aa_myf11[j], my_name_dic[,1])]
    aa_mypos <- my_Dataset[my_Dataset$label == "pos", match(aa_myf11[j], my_name_dic[,1])]
    aa_f_range <- range(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])], na.rm = T)
    ax <- pretty(floor(aa_f_range[1]):ceiling(aa_f_range[2]), n = 20)
    histneg <- hist(aa_myneg,  breaks = ax , plot = F)
    histpos <- hist(aa_mypos,  breaks = ax , plot = F)
    plot(histneg, col = c1, freq = F, main = aa_myf11[j]) # Plot 1st histogram using a transparent color
    plot(histpos, col = c2, add = TRUE, freq = F) # Add 2nd histogram using different color
    
  }
  dev.off()
  aa_myf1 <- lncRNA_feature_context_rand_summary_list[[i]][lncRNA_feature_context_rand_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("sequ_TR_", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_Hist_ovl/fh_",aaln[i],"_Context_rand.png"),    # create PNG for the heat map        
      width = 16*300,        # 5 x 300 pixels
      height = 16*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(2,2,2,2))
  for(j in 1:length(aa_myf11)){
    aa_myneg <- my_Dataset[my_Dataset$label == "neg", match(aa_myf11[j], my_name_dic[,1])]
    aa_mypos <- my_Dataset[my_Dataset$label == "pos", match(aa_myf11[j], my_name_dic[,1])]
    aa_f_range <- range(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])], na.rm = T)
    ax <- pretty(floor(aa_f_range[1]):ceiling(aa_f_range[2]), n = 20)
    histneg <- hist(aa_myneg,  breaks = ax , plot = F)
    histpos <- hist(aa_mypos,  breaks = ax , plot = F)
    plot(histneg, col = c1, freq = F, main = aa_myf11[j]) # Plot 1st histogram using a transparent color
    plot(histpos, col = c2, add = TRUE, freq = F) # Add 2nd histogram using different color
    
  }
  dev.off()
  aa_myf1 <- lncRNA_feature_context_chunk_summary_list[[i]][lncRNA_feature_context_chunk_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("sequ_TR_", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_Hist_ovl/fh_",aaln[i],"_Context_chunk.png"),    # create PNG for the heat map        
      width = 16*300,        # 5 x 300 pixels
      height = 16*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(2,2,2,2))
  for(j in 1:length(aa_myf11)){
    aa_myneg <- my_Dataset[my_Dataset$label == "neg", match(aa_myf11[j], my_name_dic[,1])]
    aa_mypos <- my_Dataset[my_Dataset$label == "pos", match(aa_myf11[j], my_name_dic[,1])]
    aa_f_range <- range(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])], na.rm = T)
    ax <- pretty(floor(aa_f_range[1]):ceiling(aa_f_range[2]), n = 20)
    histneg <- hist(aa_myneg,  breaks = ax , plot = F)
    histpos <- hist(aa_mypos,  breaks = ax , plot = F)
    plot(histneg, col = c1, freq = F, main = aa_myf11[j]) # Plot 1st histogram using a transparent color
    plot(histpos, col = c2, add = TRUE, freq = F) # Add 2nd histogram using different color
    
  }
  dev.off()
  
  
}
######## plot qqplots instead of the above overlapping histograms

aafile <- list.files("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/")

c1 <- rgb(0,0,255,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,0,0, max = 255, alpha = 80, names = "lt.red")

for(i in 1:length(aafile)){
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/", aafile[i]))
  aa_myf1 <- lncRNA_feature_seq_chunk_summary_list[[i]][lncRNA_feature_seq_chunk_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("sequ_TR_", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_qq/fq_",aaln[i],"_Seq_chunk.png"),    # create PNG for the heat map        
      width = 18*300,        # 5 x 300 pixels
      height = 18*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(4,4,3,2))
  for(j in 1:length(aa_myf11)){
    aa_myneg <- my_Dataset[my_Dataset$label == "neg", match(aa_myf11[j], my_name_dic[,1])]
    aa_mypos <- my_Dataset[my_Dataset$label == "pos", match(aa_myf11[j], my_name_dic[,1])]
    # aa_f_range <- range(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])], na.rm = T)
    # ax <- pretty(floor(aa_f_range[1]):ceiling(aa_f_range[2]), n = 20)
    aar <- range(c(aa_myneg,aa_mypos))
    
    qqplot(x = aa_myneg, y = aa_mypos,plot.it = T,xlab = "-",ylab = "+", main = aa_myf11[j],pch=19,cex = 0.7, xlim = aar, ylim = aar)
    abline(a = 0,b = 1,col=2)
    # histneg <- hist(aa_myneg,  breaks = ax , plot = F)
    # histpos <- hist(aa_mypos,  breaks = ax , plot = F)
    # plot(histneg, col = c1, freq = F, main = aa_myf11[j]) # Plot 1st histogram using a transparent color
    # plot(histpos, col = c2, add = TRUE, freq = F) # Add 2nd histogram using different color
    #boxplot(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])]~my_Dataset$label, main = aa_myf11[j])
  }
  
  aa_myf1 <- lncRNA_feature_seq_rand_summary_list[[i]][lncRNA_feature_seq_rand_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("sequ_TR_", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  dev.off()
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_qq/fq_",aaln[i],"_Seq_rand.png"),    # create PNG for the heat map        
      width = 18*300,        # 5 x 300 pixels
      height = 18*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(4,4,3,2))
  for(j in 1:length(aa_myf11)){
    aa_myneg <- my_Dataset[my_Dataset$label == "neg", match(aa_myf11[j], my_name_dic[,1])]
    aa_mypos <- my_Dataset[my_Dataset$label == "pos", match(aa_myf11[j], my_name_dic[,1])]
    aar <- range(c(aa_myneg,aa_mypos))
    
    qqplot(x = aa_myneg, y = aa_mypos,plot.it = T,xlab = "-",ylab = "+", main = aa_myf11[j],pch=19,cex = 0.7, xlim = aar, ylim = aar)
    abline(a = 0,b = 1,col=2)
    # aa_f_range <- range(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])], na.rm = T)
    # ax <- pretty(floor(aa_f_range[1]):ceiling(aa_f_range[2]), n = 20)
    # histneg <- hist(aa_myneg,  breaks = ax , plot = F)
    # histpos <- hist(aa_mypos,  breaks = ax , plot = F)
    # plot(histneg, col = c1, freq = F, main = aa_myf11[j]) # Plot 1st histogram using a transparent color
    # plot(histpos, col = c2, add = TRUE, freq = F) # Add 2nd histogram using different color
    
  }
  dev.off()
  aa_myf1 <- lncRNA_feature_context_rand_summary_list[[i]][lncRNA_feature_context_rand_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("sequ_TR_", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_qq/fq_",aaln[i],"_Context_rand.png"),    # create PNG for the heat map        
      width = 18*300,        # 5 x 300 pixels
      height = 18*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(4,4,3,2))
  for(j in 1:length(aa_myf11)){
    aa_myneg <- my_Dataset[my_Dataset$label == "neg", match(aa_myf11[j], my_name_dic[,1])]
    aa_mypos <- my_Dataset[my_Dataset$label == "pos", match(aa_myf11[j], my_name_dic[,1])]
    aar <- range(c(aa_myneg,aa_mypos))
    qqplot(x = aa_myneg, y = aa_mypos,plot.it = T,xlab = "-",ylab = "+", main = aa_myf11[j],pch=19,cex = 0.7, xlim = aar, ylim = aar)
    abline(a = 0,b = 1,col=2)
    # aa_f_range <- range(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])], na.rm = T)
    # ax <- pretty(floor(aa_f_range[1]):ceiling(aa_f_range[2]), n = 20)
    # histneg <- hist(aa_myneg,  breaks = ax , plot = F)
    # histpos <- hist(aa_mypos,  breaks = ax , plot = F)
    # plot(histneg, col = c1, freq = F, main = aa_myf11[j]) # Plot 1st histogram using a transparent color
    # plot(histpos, col = c2, add = TRUE, freq = F) # Add 2nd histogram using different color
    
  }
  dev.off()
  aa_myf1 <- lncRNA_feature_context_chunk_summary_list[[i]][lncRNA_feature_context_chunk_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("sequ_TR_", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res4TX_plots/Feature_qq/fq_",aaln[i],"_Context_chunk.png"),    # create PNG for the heat map        
      width = 18*300,        # 5 x 300 pixels
      height = 18*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(4,4,3,2))
  for(j in 1:length(aa_myf11)){
    aa_myneg <- my_Dataset[my_Dataset$label == "neg", match(aa_myf11[j], my_name_dic[,1])]
    aa_mypos <- my_Dataset[my_Dataset$label == "pos", match(aa_myf11[j], my_name_dic[,1])]
    aar <- range(c(aa_myneg,aa_mypos))
    qqplot(x = aa_myneg, y = aa_mypos,plot.it = T,xlab = "-",ylab = "+", main = aa_myf11[j],pch=19,cex = 0.7, xlim = aar, ylim = aar)
    abline(a = 0,b = 1,col=2)
    # aa_f_range <- range(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])], na.rm = T)
    # ax <- pretty(floor(aa_f_range[1]):ceiling(aa_f_range[2]), n = 20)
    # histneg <- hist(aa_myneg,  breaks = ax , plot = F)
    # histpos <- hist(aa_mypos,  breaks = ax , plot = F)
    # plot(histneg, col = c1, freq = F, main = aa_myf11[j]) # Plot 1st histogram using a transparent color
    # plot(histpos, col = c2, add = TRUE, freq = F) # Add 2nd histogram using different color
    
  }
  dev.off()
  
  
}
######## get feature index for running a model for each lncRNA

aafile <- list.files("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/")

aa_model_feat_list_rand <- list()
aa_model_feat_list_chunk <- list()
aaln <- unique(aatmp_lnc_df_cv2TX$lncRNA)

load(paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/", aafile[1]))
for(i in 1:length(aafile)){
  
  aa_myf1 <- lncRNA_feature_seq_chunk_summary_list[[i]][lncRNA_feature_seq_chunk_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("sequ_TR_", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  aa_model_feat_list_chunk[[i]] <- aa_myf11 # which(my_name_dic[,1] %in% aa_myf11)
  

  
  aa_myf1 <- lncRNA_feature_seq_rand_summary_list[[i]][lncRNA_feature_seq_rand_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("sequ_TR_", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  aa_model_feat_list_rand[[i]] <- aa_myf11 #which(my_name_dic[,1] %in% aa_myf11)
  
  
  aa_myf1 <- lncRNA_feature_context_rand_summary_list[[i]][lncRNA_feature_context_rand_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("sequ_TR_", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  aa_model_feat_list_rand[[i]] <- c(aa_model_feat_list_rand[[i]], aa_myf11)# which(my_name_dic[,1] %in% ))
  

  aa_myf1 <- lncRNA_feature_context_chunk_summary_list[[i]][lncRNA_feature_context_chunk_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("sequ_TR_", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  aa_model_feat_list_chunk[[i]] <- c(aa_model_feat_list_chunk[[i]], aa_myf11) # which(my_name_dic[,1] %in% aa_myf11))
  my_200_feat <- aa_model_feat_list_rand[[i]]
  save(list = c("my_200_feat"), file = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Feature200/perlncIndex/index_", aaln[i], "_R.RData"))
  my_200_feat <- aa_model_feat_list_chunk[[i]]
  save(list = c("my_200_feat"), file = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Feature200/perlncIndex/index_", aaln[i], "_C.RData"))
  
  
}
names(aa_model_feat_list_chunk) <- aaln
names(aa_model_feat_list_rand) <- aaln


# write jobs visualization
# plot the best model and 20 most important fearutes

aa_filt_dic <- c("no_filter",
                 "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_CV2_filt/my_filter_RBP.txt",
                 "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_CV2_filt/my_filter_repeat.txt",
                 "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_CV2_filt/my_filter_both.txt")

for(i in 1:length(aaln)){
  cat(c("Rscript",
        "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/plot_gviz_200Feat.R",
        aaln[i],
        "R", 
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_200feat/Gviz/", aaln[i]), 
        "20", 
        c("top200Features_RBPFilter"),#, "top200Features_BothFilter", "top200Features_RepFilter","top200Features_noFilter"), 
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_CVfirst_dataset_", aaln[i], ".RData"),
        aa_filt_dic[2],
        "\n"), file = "plot_gviz_top20_200Feat.job", append = !(i==1))
  cat(c("Rscript",
        "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/plot_gviz_200Feat.R",
        aaln[i],
        "C", 
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_200feat/Gviz/", aaln[i]), 
        "20", 
        c("top200Features_RBPFilter"),#, "top200Features_BothFilter", "top200Features_RepFilter","top200Features_noFilter"),
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_CVfirst_dataset_", aaln[i], ".RData"),
        aa_filt_dic[2],
        "\n"), file = "plot_gviz_top20_200Feat.job", append = T)
  
}




###########################
# Reading results of top200_feature model

for(i in 1:length(aa_un_own)){
  cat(c("Rscript --vanilla /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/gather_perf.R",
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_200feat/Learned_models/",aa_un_own[i] ),
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_200feat/Combined_results/",aa_un_own[i], "__perfMat.RData\n" )),
      sep = " ", 
      file = "~/Documents/Shayan/BioInf/lncRNA/gather_perf_all_CV2_200feat.job",
      append = !(i == 1))
}

aafiles <- list.files("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Feature200/Combined_res")
aafn <- unlist(lapply(strsplit(aafiles, "__"), "[[", 1))
aa_all_auroc_list_cv2TX200feat <- list()
aa_all_auprc_list_cv2TX200feat <- list()
aa_all_impor_list_cv2TX200feat<- list()

for(i in 1:length(aafiles)){
  print(i)
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Feature200/Combined_res/", aafiles[i]))
  aa_all_auroc_list_cv2TX200feat[[i]] <- my_perf[[1]][[1]]
  aa_all_auprc_list_cv2TX200feat[[i]] <- my_perf[[2]][[1]]
  aa_all_impor_list_cv2TX200feat[[i]] <- my_perf[[3]][[1]]
}
names(aa_all_auroc_list_cv2TX200feat) <- aafn
names(aa_all_auprc_list_cv2TX200feat) <- aafn
names(aa_all_impor_list_cv2TX200feat) <- aafn

View(aa_all_auroc_list_cv2TX200feat$Malat1)




aac <- max(unlist(lapply(aa_all_auprc_list_cv2TX200feat, nrow)))

aa_all_auroc_rand_cv2TX200feat <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auroc_rand_cv2TX200feat) <- rownames(aa_all_auprc_list_cv2TX200feat[[2]])
rownames(aa_all_auroc_rand_cv2TX200feat) <- aafn

aa_all_auroc_chunk_cv2TX200feat <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auroc_chunk_cv2TX200feat) <- rownames(aa_all_auprc_list_cv2TX200feat[[2]])
rownames(aa_all_auroc_chunk_cv2TX200feat) <- aafn

aa_all_auprc_rand_cv2TX200feat <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auprc_rand_cv2TX200feat) <- rownames(aa_all_auprc_list_cv2TX200feat[[2]])
rownames(aa_all_auprc_rand_cv2TX200feat) <- aafn

aa_all_auprc_chunk_cv2TX200feat <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auprc_chunk_cv2TX200feat) <- rownames(aa_all_auprc_list_cv2TX200feat[[2]])
rownames(aa_all_auprc_chunk_cv2TX) <- aafn


aa_all_auprc_list_cv2_diffTX200feat <- list()
for(i in 1:length(aa_all_auprc_list_cv2TX200feat)){
  aa_all_auprc_list_cv2_diffTX200feat[[i]] <- aa_all_auprc_list_cv2TX200feat[[i]]
  for(j in 1:nrow(aa_all_auprc_list_cv2TX200feat[[i]])){
    aa_all_auprc_list_cv2_diffTX200feat[[i]][j,] <- aa_all_auprc_list_cv2_diffTX200feat[[i]][j,] - aa_auprc_base[i,]
  }
}
for(i in 1:length(aa_all_auroc_list_cv2TX200feat)){
  aatmmpccv_roc <- rowMeans(aa_all_auroc_list_cv2TX200feat[[i]][,c(1:5)], na.rm = T)
  aatmmprcv_roc <- rowMeans(aa_all_auroc_list_cv2TX200feat[[i]][,c(6:10)], na.rm = T)
  
  aatmmpccv_prc <- rowMeans(aa_all_auprc_list_cv2_diffTX200feat[[i]][,c(1:5)], na.rm = T)
  aatmmprcv_prc <- rowMeans(aa_all_auprc_list_cv2_diffTX200feat[[i]][,c(6:10)], na.rm = T)
  
  aacolord <- match(rownames(aa_all_auroc_list_cv2TX200feat[[i]]), 
                    colnames(aa_all_auroc_chunk_cv2TX200feat))
  
  aa_all_auroc_chunk_cv2TX200feat[i,aacolord] <- aatmmpccv_roc
  aa_all_auroc_rand_cv2TX200feat[i,aacolord] <- aatmmprcv_roc
  aa_all_auprc_chunk_cv2TX200feat[i,aacolord] <- aatmmpccv_prc
  aa_all_auprc_rand_cv2TX200feat[i,aacolord] <- aatmmprcv_prc
}

View(aa_all_auroc_rand_cv2TX200feat)
View(aa_all_auroc_chunk_cv2TX200feat)

aa_all_auprc_rand_cv2TX200feat
aa_all_auprc_chunk_cv2TX200feat




# biomTrack <- BiomartGeneRegionTrack(genome = "mm9",
#                                     chromosome = chr, start = 20000000, end = 21000000,
#                                     name = "ENSEMBL")
# 
# 
# plotTracks(biomTrack, col.line = NULL, col = NULL)

##################################################################################################################################################################
##################################################################################################################################################################
# why are some regions zero in all motif scans? --> there was a bug because of duplicated names--> reran the analysis
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/Partition_6_CVfirst_dataset_Neat1.RData")
my_Dataset
aax1 <- which(my_name_dic[, 1] %in% "RBP_pair") 
aax2 <- which(my_name_dic[, 1] %in% "PPI_pair") 
aax3 <- which(my_name_dic[, 1] %in% "TFscan_pair") 
aax4 <- which(my_name_dic[, 1] %in% "Rbmx__RBPscan") 

aaxc <- which(rowSums(my_Dataset[,c(aax1, aax2, aax3, aax4)]) == 0)
table(my_Dataset$label[aaxc])


aax1 <- which(my_name_dic[, 1] %in% "GG__DiFreq") 
aax2 <- which(my_name_dic[, 1] %in% "GC__DiFreq") 
aax3 <- which(my_name_dic[, 1] %in% "CC__DiFreq") 
aax4 <- which(my_name_dic[, 1] %in% "AT__DiFreq") 

aaxc <- which(rowSums(my_Dataset[,c(aax1, aax2, aax3, aax4)]) == 0)

table(my_Dataset$label[aaxc])


Partition_6_RBP_scanned <- read.table("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/All_RBP_profile_P6.txt", header = T)
#load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_RBP_scanned.RData")

dim(Partition_6_RBP_scanned)
aalz1 <- rownames(Partition_6_RBP_scanned)[(rowSums(Partition_6_RBP_scanned) == 0)]

table(Partition_6_dfs$owner[Partition_6_dfs$tile_name %in% aalz])/table(Partition_6_dfs$owner)
table(Partition_6_dfs$owner[Partition_6_dfs$tile_name %in% aalz], Partition_6_dfs$label[Partition_6_dfs$tile_name %in% aalz]) / table(Partition_6_dfs$owner, Partition_6_dfs$label)


Partition_6_TF_scanned <- read.table("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/All_TF_profile_P6.txt", header = T)
dim(Partition_6_TF_scanned)
aalz2 <- rownames(Partition_6_RBP_scanned)[(rowSums(Partition_6_TF_scanned) == 0)]

table(Partition_6_dfs$owner[Partition_6_dfs$tile_name %in% aalz])/table(Partition_6_dfs$owner)
table(Partition_6_dfs$owner[Partition_6_dfs$tile_name %in% aalz], Partition_6_dfs$label[Partition_6_dfs$tile_name %in% aalz]) / table(Partition_6_dfs$owner, Partition_6_dfs$label)


aamx <- Partition_6_dfs$tile_name[as.numeric(intersect(aalz1, aalz2))]
aamxGR <- Partition_6_dfs_GR[match(aamx, Partition_6_dfs_GR$tile)]

aamxseq <- getSeq(x = BSgenome.Mmusculus.UCSC.mm9, names = aamxGR)
toberemoved_tiles_partition6 <- unique(aamxGR$tile)



##################################################################################################################################################################
##################################################################################################################################################################
# gathering performance of new runs --> after fixing issue with 0 motif scan scores for about 50000 regions

aa_un_own <- sort(unique(Partition_6_random_chunk_cv_df$owner))


for(i in 1:length(aa_un_own)){
  cat(c("Rscript --vanilla /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/gather_perf_wt.R",
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_NOD/Learned_models/",aa_un_own[i] ),
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_NOD/Combined_results/",aa_un_own[i], "__perfMat.RData\n" )),
      sep = " ", 
      file = "~/Documents/Shayan/BioInf/lncRNA/gather_perf_all_CV2_svd_noD.job",
      append = !(i == 1))
}

##################################################################################################################################################################
##################################################################################################################################################################
# reading results
aafiles <- list.files("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX/")
aafn <- unlist(lapply(strsplit(aafiles, "__"), "[[", 1))
aa_all_auroc_list_cv2TX_NoD <- list()
aa_all_auprc_list_cv2TX_NoD <- list()
aa_all_impor_list_cv2TX_NoD <- list()

aa_all_auroc_list_cv2TX_NoD_train <- list()
aa_all_auprc_list_cv2TX_NoD_train <- list()

aallnamm <- character(0)
for(i in 1:length(aafiles)){
  print(i)
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX/", aafiles[i]))
  aa_all_auroc_list_cv2TX_NoD[[i]] <- my_perf[[1]][[1]]
  aa_all_auprc_list_cv2TX_NoD[[i]] <- my_perf[[2]][[1]]
  aa_all_impor_list_cv2TX_NoD[[i]] <- my_perf[[3]][[1]]
  aa_all_auroc_list_cv2TX_NoD_train[[i]] <- my_perf[[4]][[1]]
  aa_all_auprc_list_cv2TX_NoD_train[[i]] <- my_perf[[5]][[1]]
  
  
  aallnamm <- union(aallnamm, rownames(aa_all_auprc_list_cv2TX_NoD[[i]]))
}
names(aa_all_auroc_list_cv2TX_NoD) <- aafn
names(aa_all_auprc_list_cv2TX_NoD) <- aafn
names(aa_all_impor_list_cv2TX_NoD) <- aafn
names(aa_all_auroc_list_cv2TX_NoD_train) <- aafn
names(aa_all_auprc_list_cv2TX_NoD_train) <- aafn

View(aa_all_auroc_list_cv2TX_NoD$Neat1)




aac <- length(aallnamm)

aa_all_auroc_rand_cv2TX_NoD <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auroc_rand_cv2TX_NoD) <- aallnamm
rownames(aa_all_auroc_rand_cv2TX_NoD) <- aafn

aa_all_auroc_chunk_cv2TX_NoD <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auroc_chunk_cv2TX_NoD) <- aallnamm
rownames(aa_all_auroc_chunk_cv2TX_NoD) <- aafn

aa_all_auprc_rand_cv2TX_NoD <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auprc_rand_cv2TX_NoD) <- aallnamm
rownames(aa_all_auprc_rand_cv2TX_NoD) <- aafn

aa_all_auprc_chunk_cv2TX_NoD <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auprc_chunk_cv2TX_NoD) <- aallnamm
rownames(aa_all_auprc_chunk_cv2TX_NoD) <- aafn


aa_all_auprc_list_cv2_diffTX_NoD <- list()
for(i in 1:length(aa_all_auprc_list_cv2TX_NoD)){
  aa_all_auprc_list_cv2_diffTX_NoD[[i]] <- aa_all_auprc_list_cv2TX_NoD[[i]]
  for(j in 1:nrow(aa_all_auprc_list_cv2TX_NoD[[i]])){
    aa_all_auprc_list_cv2_diffTX_NoD[[i]][j,] <- aa_all_auprc_list_cv2_diffTX_NoD[[i]][j,] - aa_auprc_base[i,]
  }
}
for(i in 1:length(aa_all_auroc_list_cv2TX_NoD)){
  aatmmpccv_roc <- rowMeans(aa_all_auroc_list_cv2TX_NoD[[i]][,c(1:5)], na.rm = T)
  aatmmprcv_roc <- rowMeans(aa_all_auroc_list_cv2TX_NoD[[i]][,c(6:10)], na.rm = T)
  
  aatmmpccv_prc <- rowMeans(aa_all_auprc_list_cv2_diffTX_NoD[[i]][,c(1:5)], na.rm = T)
  aatmmprcv_prc <- rowMeans(aa_all_auprc_list_cv2_diffTX_NoD[[i]][,c(6:10)], na.rm = T)
  
  aacolord <- match(rownames(aa_all_auroc_list_cv2TX_NoD[[i]]), 
                    colnames(aa_all_auroc_chunk_cv2TX_NoD))
  
  aa_all_auroc_chunk_cv2TX_NoD[i,aacolord] <- aatmmpccv_roc
  aa_all_auroc_rand_cv2TX_NoD[i,aacolord] <- aatmmprcv_roc
  aa_all_auprc_chunk_cv2TX_NoD[i,aacolord] <- aatmmpccv_prc
  aa_all_auprc_rand_cv2TX_NoD[i,aacolord] <- aatmmprcv_prc
}

library(RColorBrewer)
library(gplots)
png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/AUROC_rand_heatmap.png", 
    width = 17*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
range(aa_all_auroc_rand_cv2TX_NoD, na.rm = T)
#aa_all_auroc_rand_cv2TX_NoD[is.na(aa_all_auroc_rand_cv2TX_NoD)] <- 0.5
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(25)
heatmap.2((aa_all_auroc_rand_cv2TX_NoD), trace='none', col=(cols), breaks = seq(0.5, 1, 0.02), margins = c(12,10), na.rm = F)
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/AUROC_chunk_heatmap.png", 
    width = 17*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#sum(is.na(aa_all_auroc_chunk_cvfirst))
#aa_all_auroc_chunk_cvfirst[is.na(aa_all_auroc_chunk_cvfirst)] <- 0
range(aa_all_auroc_chunk_cv2TX_NoD, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(25)
#aa_all_auroc_chunk_cv2TX_NoD[is.na(aa_all_auroc_chunk_cv2TX_NoD)] <- 0.5

heatmap.2((aa_all_auroc_chunk_cv2TX_NoD), trace='none', col=(cols), breaks = seq(0.5, 1, 0.02), margins = c(12,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/AUPRC_rand_heatmap.png", 
    width = 17*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#sum(is.na(aa_all_auprc_rand_cvfirst))
#aa_all_auprc_rand_cvfirst[is.na(aa_all_auprc_rand_cvfirst)] <- 0
range(aa_all_auprc_rand_cv2TX_NoD, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(20)
heatmap.2((aa_all_auprc_rand_cv2TX_NoD), trace='none', col=(cols), breaks = seq(0.0, 0.40, 0.02), margins = c(12,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/AUPRC_chunk_heatmap.png", 
    width = 17*300,
    height = 10*300,
    res = 300,
    pointsize = 11)
#sum(is.na(aa_all_auprc_chunk_cvfirst))
#aa_all_auprc_chunk_cvfirst[is.na(aa_all_auprc_chunk_cvfirst)] <- 0
range(aa_all_auprc_chunk_cv2TX_NoD, na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(20)
heatmap.2((aa_all_auprc_chunk_cv2TX_NoD), trace='none', col=(cols), breaks = seq(0.0, 0.40, 0.02), margins = c(12,10))
dev.off()

#######
# make barplots of best models with distance and best ones without distance
# make a dataframe with the following columns: lncRNA, partition_type, AUROC_dist, AUROC_nodist, AUPRC_dist, AUPRC_nodist
aa_perf_df_o2 <- data.frame(lncRNA = character(0), 
                           model_partition= character(0), # can be "Wdist_chunk", "Wdist_rand", "WOdist_chunk", "WOdist_rand"
                           AUROC = numeric(0), 
                           AUPRC = numeric(0))

aa_all_col <- which(colnames(aa_all_auroc_rand_cv2TX_NoD) %in% "sequ_Chrom_Meth_AX_ChIP_trnsp")
for(i in 1:nrow(aa_all_auroc_rand_cv2TX_NoD)){
  
  aamd1 <- "Full_rand"
  #aamaxrd_roc <- max(aa_all_auroc_rand_cv2TX_NoD[i,], na.rm = T)
  #aamaxrd_prc <- max(aa_all_auprc_rand_cv2TX_NoD[i,], na.rm = T)
  aamaxrd_roc <- aa_all_auroc_rand_cv2TX_NoD[i,aa_all_col]
  aamaxrd_prc <- aa_all_auprc_rand_cv2TX_NoD[i,aa_all_col]
  aaaxx1 <- c(rownames(aa_all_auroc_rand_cv2TX_NoD)[i], aamd1, aamaxrd_roc, aamaxrd_prc)
  aamd2 <- "Full_chunk"
  #aamaxcd_roc <- max(aa_all_auroc_chunk_cv2TX_NoD[i,], na.rm = T)
  #aamaxcd_prc <- max(aa_all_auprc_chunk_cv2TX_NoD[i,], na.rm = T)
  aamaxcd_roc <- aa_all_auroc_chunk_cv2TX_NoD[i,aa_all_col]
  aamaxcd_prc <- aa_all_auprc_chunk_cv2TX_NoD[i,aa_all_col]
  aaaxx2 <- c(rownames(aa_all_auroc_rand_cv2TX_NoD)[i], aamd2, aamaxcd_roc, aamaxcd_prc)
  aa_perf_df_o2 <- rbind(aa_perf_df_o2, rbind(aaaxx1,aaaxx2))
}

aa_perf_df_all <- rbind(aa_perf_df_o, aa_perf_df_o2)
colnames(aa_perf_df_all) <- c("lncRNA", "Type", "AUROC", "AUPRC")
aa_perf_df_all$Type <- factor(aa_perf_df_all$Type, levels = unique(aa_perf_df_all$Type))
aa_perf_df_all$AUROC <- as.numeric(aa_perf_df_all$AUROC)
aa_perf_df_all$AUPRC <- as.numeric(aa_perf_df_all$AUPRC)

aaakkls <- summarySE(data = aa_perf_df_all, 
                     measurevar=c("AUROC"), 
                     groupvars=c( "Type"), 
                     na.rm=T)
aaakkls2 <- summarySE(data = aa_perf_df_all, 
                     measurevar=c("AUPRC"), 
                     groupvars=c( "Type"), 
                     na.rm=T)

aa_perf_df_o3 <- data.frame(lncRNA = rep("Average", 4), 
                            Type= aaakkls$Type, 
                            AUROC = aaakkls$AUROC, 
                            AUPRC = aaakkls2$AUPRC)
aa_perf_df_all <- rbind(aa_perf_df_all, aa_perf_df_o3)
aa_perf_df_all$lncRNA <- factor(aa_perf_df_all$lncRNA, levels = c(setdiff(sort(unique(aa_perf_df_all$lncRNA)), "Average"),"Average" ))

ggplot(aa_perf_df_all, aes(x=lncRNA, y=AUROC, fill=Type)) + 
  # geom_errorbar(aes(ymin=auroc_mean_delta-se, ymax=auroc_mean_delta+se),
  #               width=.2,                    # Width of the error bars
  #               position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("AUROC")+
  xlab("lncRNA") + 
  #ylim(0.5, 1)+
  coord_cartesian(ylim = c(0.5, 1))+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 11), 
        axis.title.y=element_text(size = 11),
        axis.text.x = element_text(size = 10, angle = 90),
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"))

ggplot(aa_perf_df_all, aes(x=lncRNA, y=AUPRC, fill=Type)) + 
  # geom_errorbar(aes(ymin=auroc_mean_delta-se, ymax=auroc_mean_delta+se),
  #               width=.2,                    # Width of the error bars
  #               position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("AUPRC")+
  xlab("lncRNA") + 
  #ylim(0.5, 1)+
  #coord_cartesian(ylim = c(0.5, 1))+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 11), 
        axis.title.y=element_text(size = 11),
        axis.text.x = element_text(size = 10, angle = 90),
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"))
#######


family_mat_list_all2TX_NoD <- list()

for(i in 1:length(aa_all_auroc_list_cv2TX_NoD)){
  print(aafn[i])
  family_mat_list <- list()
  aana <- rownames(aa_all_auroc_list_cv2TX_NoD[[i]])
  aana_sp <- strsplit(aana, "_")
  aana_sp_sort <- lapply(aana_sp, sort)
  aana_allfam <- sort(unique(unlist(aana_sp)))
  if("sequ" %in% aana_allfam){
    aana_allfam <- setdiff(aana_allfam,c("sequ"))
  }
  
  for(j in 1:length(aana_allfam)){
    aatmp1 <- which(unlist(lapply(aana_sp, function(x) aana_allfam[j] %in% x)))
    if(!(aana_allfam[j] %in% c("TXRB", "TXRE"))){
      aatmp111 <- which(unlist(lapply(aana_sp[aatmp1], function(x) "TXRB" %in% x)))
      aatmp1112 <- which(unlist(lapply(aana_sp[aatmp1], function(x) "TXRE" %in% x)))
      aatmp1 <- aatmp1[setdiff(c(1:length(aatmp1)), union(aatmp111, aatmp1112))]
    }else{
      aatmp111 <- which(unlist(lapply(aana_sp[aatmp1], function(x) "sequ" %in% x)))
      aatmp1 <- aatmp1[aatmp111] # evaluating the importance of transcriptional filtering only in models that already use the transcription information separately --> so the improvement is not due to introduction of new information
      
    }
    #print(aana_allfam[j])
    #print(aana[aatmp1])
    aatmp2 <- numeric(length = length(aatmp1))
    for(k in 1:length(aatmp1)){
      if(length(aana_sp[[aatmp1[k]]]) > 1){
        aanew_name <- setdiff(aana_sp_sort[[aatmp1[k]]], aana_allfam[j])
        #if((length(aanew_name) == 2) & identical(aanew_name, c("dist", "sequ"))){
        #aatmp2[k] <- 0
        #}else{
        aastm <- which(unlist(lapply(aana_sp_sort, function(x) identical(aanew_name, x))))
        if(length(aastm) == 0){
          print(paste("nomatch", names(aa_all_auroc_list_cv2TX_NoD)[i], aana_allfam[j],"###", paste(aana_sp_sort[[aatmp1[k]]], collapse = "_"),"###",
                      paste(aanew_name, collapse = "||")))
          aatmp2[k] <- 0
        }else{
          aatmp2[k] <- aastm 
        }
        #}
        # if((length(aanew_name) == 2) & all(aanew_name == c("sequ", "dist"))){
        #   aanew_name <- sort(unlist(strsplit("kmerFq_motSc_Rep_pairs_triplx", "_")))
        #   aastm <- which(unlist(lapply(aana_sp_sort, function(x) identical(aanew_name, x))))
        #   if(length(aastm) == 0){
        #     aanew_name <- sort(unlist(strsplit("kmerFq_motSc_Rep_triplx", "_")))
        #     aastm <- which(unlist(lapply(aana_sp_sort, function(x) identical(aanew_name, x))))
        #     if(length(aastm) == 0){
        #       aanew_name <- sort(unlist(strsplit("kmerFq_Rep_triplx", "_")))
        #     }
        #   }
        # }
        
        
      }else{
        aatmp2[k] <- 0
      }
      
    }
    family_mat_list[[j]] <- cbind(aatmp1, aatmp2)
  }
  names(family_mat_list) <- aana_allfam
  family_mat_list_all2TX_NoD[[i]] <- family_mat_list
}

names(family_mat_list_all2TX_NoD) <- aafn
length(family_mat_list_all2TX_NoD[[1]])





aatmp_lnc_df_list_cv2TX_NoD <- list()
for(i in 1:length(family_mat_list_all2TX_NoD)){
  aatmp_lnc_fam_df_list <- list()
  print(aafn[i])
  for(j in 1:length(family_mat_list_all2TX_NoD[[i]])){#going over families
    print(names(family_mat_list_all2TX_NoD[[i]])[j])
    aatmptroc <- matrix(nrow = nrow(family_mat_list_all2TX_NoD[[i]][[j]]), 
                        ncol = ncol(aa_all_auroc_list_cv2TX_NoD[[i]]))
    aatmptprc <- matrix(nrow = nrow(family_mat_list_all2TX_NoD[[i]][[j]]), 
                        ncol = ncol(aa_all_auprc_list_cv2TX_NoD[[i]]))
    for(k in 1:nrow(family_mat_list_all2TX_NoD[[i]][[j]])){ #going over comparisons
      if(family_mat_list_all2TX_NoD[[i]][[j]][k,2] == 0){
        aatmpcomproc <- rep(0.5, ncol(aa_all_auroc_list_cv2TX_NoD[[i]]))
        aatmpcompprc <- aa_auprc_base[i, ]
      }else{
        aatmpcomproc <- aa_all_auroc_list_cv2TX_NoD[[i]][family_mat_list_all2TX_NoD[[i]][[j]][k,2],]
        aatmpcompprc <- aa_all_auprc_list_cv2TX_NoD[[i]][family_mat_list_all2TX_NoD[[i]][[j]][k,2],]
      }
      aatmptroc[k,] <- aa_all_auroc_list_cv2TX_NoD[[i]][family_mat_list_all2TX_NoD[[i]][[j]][k,1],] - aatmpcomproc
      aatmptprc[k,] <- aa_all_auprc_list_cv2TX_NoD[[i]][family_mat_list_all2TX_NoD[[i]][[j]][k,1],] - aatmpcompprc
    }
    aatmptroc_num <- colSums(!is.na(aatmptroc))
    aatmptprc_num <- colSums(!is.na(aatmptprc))
    
    aatmptroc_mean <- colMeans(aatmptroc, na.rm = T)
    aatmptprc_mean <- colMeans(aatmptprc, na.rm = T)
    aatmptroc_sd <- apply(aatmptroc, MARGIN = 2, FUN = sd, na.rm = T)
    aatmptprc_sd <- apply(aatmptprc, MARGIN = 2, FUN = sd, na.rm = T)
    aatmptroc_max <- apply(aatmptroc, MARGIN = 2, FUN = max, na.rm = T)
    aatmptprc_max <- apply(aatmptprc, MARGIN = 2, FUN = max, na.rm = T)
    
    aatmp_lnc_fam_df_list[[j]] <- data.frame(lncRNA = rep(aafn[i], length(aatmptroc_mean)), 
                                             feature_family = rep(names(family_mat_list_all2TX_NoD[[i]])[j], length(aatmptroc_mean)), 
                                             partition=c(paste0("CCV", c(1:5)),paste0("RCV", c(1:5))),
                                             partition_type = c(rep("chunk", 5), rep("random", 5)),
                                             nu_comparison = aatmptroc_num, 
                                             auroc_mean_delta = aatmptroc_mean,
                                             auroc_sd_delta = aatmptroc_sd,
                                             auprc_mean_delta = aatmptprc_mean,
                                             auprc_sd_delta = aatmptprc_sd,
                                             auroc_max_delta = aatmptroc_max,
                                             auprc_max_delta = aatmptprc_max)
    
  }
  aatmp_lnc_df_list_cv2TX_NoD[[i]] <- do.call(rbind, aatmp_lnc_fam_df_list)
  
}

aatmp_lnc_df_cv2TX_NoD <- do.call(rbind, aatmp_lnc_df_list_cv2TX_NoD)
View(aatmp_lnc_df_cv2TX_NoD)

###########


aatttgc <- summarySE(data = aatmp_lnc_df_cv2TX_NoD, 
                     measurevar=c("auroc_mean_delta"), 
                     groupvars=c("lncRNA", "feature_family", "partition_type"), 
                     na.rm=T)

aatttgc_1 <-  aatttgc[aatttgc$partition_type == "random",]
aatttgc_2 <-  aatttgc[aatttgc$partition_type == "chunk",]

#aatgc <- aatgc[aatgc$Gene_number %in% c(100, 500, 1000),]
ggplot(aatttgc_1, aes(x=lncRNA, y=auroc_mean_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=auroc_mean_delta-se, ymax=auroc_mean_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auroc random")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

ggplot(aatttgc_2, aes(x=lncRNA, y=auroc_mean_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=auroc_mean_delta-se, ymax=auroc_mean_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auroc chunk")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))



aatttgc2 <- summarySE(aatmp_lnc_df_cv2TX_NoD, measurevar=c("auprc_mean_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc2_1 <-  aatttgc2[aatttgc2$partition_type == "random",]
aatttgc2_2 <-  aatttgc2[aatttgc2$partition_type == "chunk",]

ggplot(aatttgc2_1, aes(x=lncRNA, y=auprc_mean_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=auprc_mean_delta-se, ymax=auprc_mean_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auprc random")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

ggplot(aatttgc2_2, aes(x=lncRNA, y=auprc_mean_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=auprc_mean_delta-se, ymax=auprc_mean_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auprc chunk")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

# create for matrices delta AUROC/PRC in Chunk/Random
# each row a lncRNA
# each column a Feature family
aatttgc <- summarySE(aatmp_lnc_df_cv2TX_NoD, measurevar=c("auroc_mean_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc_1 <-  aatttgc[aatttgc$partition_type == "random",]
aatttgc_2 <-  aatttgc[aatttgc$partition_type == "chunk",]

aatttgc2 <- summarySE(aatmp_lnc_df_cv2TX_NoD, measurevar=c("auprc_mean_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc2_1 <-  aatttgc2[aatttgc2$partition_type == "random",]
aatttgc2_2 <-  aatttgc2[aatttgc2$partition_type == "chunk",]

aa_col_seq <- c( "kmerFq", "motSc","pairs",  "Rep" ,   "triplx")
#aa_col_con <- setdiff(colnames(aa_auprc_chunk2), aa_col_seq)

aa_auprc_chunk2 <- matrix(nrow = length(unique(aatmp_lnc_df_cv2TX_NoD$lncRNA)), 
                          ncol = length(unique(aatmp_lnc_df_cv2TX_NoD$feature_family)))

rownames(aa_auprc_chunk2) <- unique(aatmp_lnc_df_cv2TX_NoD$lncRNA)
colnames(aa_auprc_chunk2) <- unique(aatmp_lnc_df_cv2TX_NoD$feature_family)

aa_auprc_rand2 <- matrix(nrow = length(unique(aatmp_lnc_df_cv2TX_NoD$lncRNA)), 
                         ncol = length(unique(aatmp_lnc_df_cv2TX_NoD$feature_family)))

rownames(aa_auprc_rand2) <- unique(aatmp_lnc_df_cv2TX_NoD$lncRNA)
colnames(aa_auprc_rand2) <- unique(aatmp_lnc_df_cv2TX_NoD$feature_family)

aa_auroc_chunk2 <- matrix(nrow = length(unique(aatmp_lnc_df_cv2TX_NoD$lncRNA)), 
                          ncol = length(unique(aatmp_lnc_df_cv2TX_NoD$feature_family)))

rownames(aa_auroc_chunk2) <- unique(aatmp_lnc_df_cv2TX_NoD$lncRNA)
colnames(aa_auroc_chunk2) <- unique(aatmp_lnc_df_cv2TX_NoD$feature_family)

aa_auroc_rand2 <- matrix(nrow = length(unique(aatmp_lnc_df_cv2TX_NoD$lncRNA)), 
                         ncol = length(unique(aatmp_lnc_df_cv2TX_NoD$feature_family)))

rownames(aa_auroc_rand2) <- unique(aatmp_lnc_df_cv2TX_NoD$lncRNA)
colnames(aa_auroc_rand2) <- unique(aatmp_lnc_df_cv2TX_NoD$feature_family)

for(i in 1:nrow(aatttgc_1)){
  aa_auroc_rand2[match(aatttgc_1$lncRNA[i], rownames(aa_auroc_rand2)), 
                 match(aatttgc_1$feature_family[i], colnames(aa_auroc_rand2))] <- aatttgc_1$auroc_mean_delta[i]
}
for(i in 1:nrow(aatttgc_2)){
  aa_auroc_chunk2[match(aatttgc_2$lncRNA[i], rownames(aa_auroc_chunk2)), 
                  match(aatttgc_2$feature_family[i], colnames(aa_auroc_chunk2))] <- aatttgc_2$auroc_mean_delta[i]
}
for(i in 1:nrow(aatttgc2_1)){
  aa_auprc_rand2[match(aatttgc2_1$lncRNA[i], rownames(aa_auprc_rand2)), 
                 match(aatttgc2_1$feature_family[i], colnames(aa_auprc_rand2))] <- aatttgc2_1$auprc_mean_delta[i]
}
for(i in 1:nrow(aatttgc2_2)){
  aa_auprc_chunk2[match(aatttgc2_2$lncRNA[i], rownames(aa_auprc_chunk2)), 
                  match(aatttgc2_2$feature_family[i], colnames(aa_auprc_chunk2))] <- aatttgc2_2$auprc_mean_delta[i]
}

aa_col_seq <- c("kmerFq", "motSc","pairs",  "Rep" ,   "triplx")
aa_col_con <- setdiff(colnames(aa_auprc_chunk2), aa_col_seq)


library(gplots)
#heatmap.2(aa_auroc_rand2)
require(RColorBrewer)


#aa_auroc_rand2[is.na(aa_auroc_rand2)] <- 0
png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Family_AUROC_rand_seq.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(7)
# z <- zClust(x=aa_auroc_rand2, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auroc_rand2[,aa_col_seq])
heatmap.2(aa_auroc_rand2[,aa_col_seq], 
          trace='none', col=(cols),
          breaks = seq(-0.00, 0.07, 0.01),
          margins = c(10,10)
         )
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Family_AUROC_rand_cont.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
range(aa_auroc_rand2[,aa_col_con])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(9)
heatmap.2(aa_auroc_rand2[,aa_col_con], trace='none', col=(cols), breaks = seq(-0.00, 0.09, 0.01), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Family_AUROC_chunk_seq.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
cols <- colorRampPalette(brewer.pal(9, "RdYlBu"))(14)
# z <- zClust(x=aa_auroc_rand2, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auroc_chunk2[,aa_col_seq])
heatmap.2(aa_auroc_chunk2[,aa_col_seq], trace='none', col=rev(cols), breaks = seq(-0.07, 0.07, 0.01), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Family_AUROC_chunk_cont.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
range(aa_auroc_chunk2[,aa_col_con])
cols <- colorRampPalette(brewer.pal(9, "RdYlBu"))(22)
heatmap.2(aa_auroc_chunk2[,aa_col_con], trace='none', col=rev(cols), breaks = seq(-0.11, 0.11, 0.01), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Family_AUPRC_rand_seq.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(9)
# z <- zClust(x=aa_auroc_rand2, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auprc_rand2[,aa_col_seq])
heatmap.2(aa_auprc_rand2[,aa_col_seq], trace='none', col=(cols), breaks = seq(-0.00, 0.09, 0.01), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Family_AUPRC_rand_cont.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
range(aa_auprc_rand2[,aa_col_con])
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(10)
heatmap.2(aa_auprc_rand2[,aa_col_con], trace='none', col=(cols), breaks = seq(-0.00, 0.10, 0.01), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Family_AUPRC_chunk_seq.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
cols <- colorRampPalette(brewer.pal(9, "RdYlBu"))(18)
# z <- zClust(x=aa_auroc_rand2, scale="none", zlim=c(-3,3), method="average")
# heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv, breaks = seq(-0.06, 0.36, 0.02), margins = c(6,10))
range(aa_auprc_chunk2[,aa_col_seq])
heatmap.2(aa_auprc_chunk2[,aa_col_seq], trace='none', col=rev(cols), breaks = seq(-0.09, 0.09, 0.01), margins = c(6,10))
dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Family_AUPRC_chunk_cont.png", 
    width = 8*300,
    height = 8*300,
    res = 300,
    pointsize = 11)
range(aa_auprc_chunk2[,aa_col_con])
cols <- colorRampPalette(brewer.pal(9, "RdYlBu"))(26)
heatmap.2(aa_auprc_chunk2[,aa_col_con], trace='none', col=rev(cols), breaks = seq(-0.13, 0.13, 0.01), margins = c(6,10))
dev.off()
##################################################################################################################################################################
aatstss2 <- marginal_perf_impro(perf_list = aa_all_auprc_list_cv2TX_NoD,
                                default_mat = aadefmat,return_comp = T)

aaxna <- unique(aatstss2$marginal_perf_df$lncRNA)

for(i in 1:length(aaxna)){
  aamyd <- aatstss2$marginal_perf_df[aatstss2$marginal_perf_df$lncRNA %in% aaxna[i],]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/per_comp_boxplot_auprc/", aaxna[i], ".png"), 
      width = 8*300,
      height = 8*300,
      res = 300,
      pointsize = 11)
  aaxplt <- ggplot(aamyd, aes(x=feature_family, y=perf_delta, fill = partition_type)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(aaxna[i])+
    ylab("delta perf")+
    xlab("feature family") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}
######################################################
# for training prc
aatstss2_tr <- marginal_perf_impro(perf_list = aa_all_auprc_list_cv2TX_NoD_train,
                                default_mat = aa_auprc_base_train,
                                return_comp = T)

aaxna <- unique(aatstss2_tr$marginal_perf_df$lncRNA)

for(i in 1:length(aaxna)){
  aamyd <- aatstss2_tr$marginal_perf_df[aatstss2_tr$marginal_perf_df$lncRNA %in% aaxna[i],]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/per_comp_boxplot_auprc_train/", aaxna[i], ".png"), 
      width = 8*300,
      height = 8*300,
      res = 300,
      pointsize = 11)
  aaxplt <- ggplot(aamyd, aes(x=feature_family, y=perf_delta, fill = partition_type)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(aaxna[i])+
    ylab("delta perf")+
    xlab("feature family") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}

######################################################

aadefmat2 <- matrix(nrow = nrow(aa_auprc_base),
                    ncol = ncol(aa_auprc_base))
aadefmat2[] <- 0.5
aatstss3 <- marginal_perf_impro(perf_list = aa_all_auroc_list_cv2TX_NoD,
                                default_mat = aadefmat2,return_comp = T)

aaxna <- unique(aatstss3$marginal_perf_df$lncRNA)

for(i in 1:length(aaxna)){
  aamyd <- aatstss3$marginal_perf_df[aatstss3$marginal_perf_df$lncRNA %in% aaxna[i],]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/per_comp_boxplot_auroc/", aaxna[i], ".png"), 
      width = 8*300,
      height = 8*300,
      res = 300,
      pointsize = 11)
  aaxplt <- ggplot(aamyd, aes(x=feature_family, y=perf_delta, fill = partition_type)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(aaxna[i])+
    ylab("delta perf")+
    xlab("feature family") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}
######################################################
# for training roc
aatstss3_tr <- marginal_perf_impro(perf_list = aa_all_auroc_list_cv2TX_NoD_train,
                                   default_mat = aadefmat2,
                                   return_comp = T)

aaxna <- unique(aatstss3_tr$marginal_perf_df$lncRNA)

for(i in 1:length(aaxna)){
  aamyd <- aatstss3_tr$marginal_perf_df[aatstss3_tr$marginal_perf_df$lncRNA %in% aaxna[i],]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/per_comp_boxplot_auroc_train/", aaxna[i], ".png"), 
      width = 8*300,
      height = 8*300,
      res = 300,
      pointsize = 11)
  aaxplt <- ggplot(aamyd, aes(x=feature_family, y=perf_delta, fill = partition_type)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(aaxna[i])+
    ylab("delta perf")+
    xlab("feature family") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}

######################################################
####################################
aatyty3 <- marginal_perf_impro_pair(perf_list = aa_all_auprc_list_cv2TX_NoD,
                                    default_mat = aadefmat)

aatyty4 <- marginal_perf_impro_pair(perf_list = aa_all_auprc_list_cv2TX_NoD,
                                    default_mat = aadefmat,
                                    return_comp = T)

aatttgc2 <- summarySE(aatyty3$marginal_perf_df, measurevar=c("perf_mean_delta"), groupvars=c("lncRNA", "feature_family", "partition_type"), na.rm=T)

aatttgc2_1 <-  aatttgc2[aatttgc2$partition_type == "random",]
aatttgc2_2 <-  aatttgc2[aatttgc2$partition_type == "chunk",]

ggplot(aatttgc2_2, aes(x=lncRNA, y=perf_mean_delta, fill=feature_family)) + 
  geom_errorbar(aes(ymin=perf_mean_delta-se, ymax=perf_mean_delta+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  ylab("delta auprc random")+
  xlab("lncRNA") + 
  #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12), 
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

aaxna <- unique(aatyty4$marginal_perf_df$lncRNA)

for(i in 1:length(aaxna)){
  aamyd <- aatyty4$marginal_perf_df[aatyty4$marginal_perf_df$lncRNA %in% aaxna[i],]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/pair_per_comp_boxplot_auprc/", aaxna[i], ".png"), 
      width = 8*300,
      height = 8*300,
      res = 300,
      pointsize = 11)
  aaxplt <- ggplot(aamyd, aes(x=feature_family, y=perf_delta, fill = partition_type)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(aaxna[i])+
    ylab("delta perf")+
    xlab("feature family pair") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}

aamytsspr <- pair_indiv_compare(pair_df = aatyty4$marginal_perf_df, indiv_df = aatstss2$marginal_perf_df)
aamytsspr$diff <- aamytsspr$perf_pair - (aamytsspr$perf_1 + aamytsspr$perf_2)

aaxna <- unique(aamytsspr$lncRNA)

for(i in 1:length(aaxna)){
  aamyd <- aamytsspr[aamytsspr$lncRNA %in% aaxna[i],]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/pair_indDiff_per_comp_boxplot_auprc/", aaxna[i], ".png"), 
      width = 8*300,
      height = 8*300,
      res = 300,
      pointsize = 11)
  aaxplt <- ggplot(aamyd, aes(x=feature_family, y=diff, fill = partition_type)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(aaxna[i])+
    ylab("delta perf")+
    xlab("feature family") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}
# draw heatmap
aamytss_agg <- aggregate(x = aamytsspr[c("diff")],
                         by = aamytsspr[c("lncRNA","partition_type", "partition", "feature_family")], 
                         FUN = mean)
aamytsspr_agg2 <- aggregate(x = aamytss_agg[c("diff")],
                          by = aamytss_agg[c("lncRNA","partition_type", "feature_family")], 
                          FUN = mean)
aalncna <- unique(aamytsspr_agg2$lncRNA)
aapairr <- unique(aamytsspr_agg2$feature_family)
aa_redun_syn_mat_prc_C <- matrix(nrow = length(aalncna),
                                 ncol = length(aapairr))
rownames(aa_redun_syn_mat_prc_C) <- aalncna
colnames(aa_redun_syn_mat_prc_C) <- aapairr
aa_redun_syn_mat_prc_R <- matrix(nrow = length(aalncna),
                                 ncol = length(aapairr))
rownames(aa_redun_syn_mat_prc_R) <- aalncna
colnames(aa_redun_syn_mat_prc_R) <- aapairr
for(i in 1:nrow(aa_redun_syn_mat_prc_C)){
  aa1 <- aamytsspr_agg2[((aamytsspr_agg2$partition_type == "chunk") & (aamytsspr_agg2$lncRNA == rownames(aa_redun_syn_mat_prc_C)[i])),]
  aa2 <- aamytsspr_agg2[((aamytsspr_agg2$partition_type == "random") & (aamytsspr_agg2$lncRNA == rownames(aa_redun_syn_mat_prc_C)[i])),]
  
  aa_redun_syn_mat_prc_C[i, ] <- aa1$diff[match(aa1$feature_family, colnames(aa_redun_syn_mat_prc_C))]
  aa_redun_syn_mat_prc_R[i, ] <- aa2$diff[match(aa2$feature_family, colnames(aa_redun_syn_mat_prc_R))]
}

aabreak <- quantile(aa_redun_syn_mat_prc_C, prob = seq(0,1,0.01),na.rm = T)
aa_redun_syn_mat_prc_C[is.na(aa_redun_syn_mat_prc_C)] <- 0
cols <- colorRampPalette(brewer.pal(9, "RdBu"))(length(aabreak)-1)
png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/FeatureFamily_synergy_chunk_prc.png",    # create PNG for the heat map        
    width = 10*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11)
heatmap.2((aa_redun_syn_mat_prc_C), 
          trace='none', Colv = T,
          dendrogram = "both",
          #ColSideColors = aa_col_color,
          col=rev(cols), 
          breaks = aabreak,
          margins = c(10,10))
dev.off()

aabreak <- quantile(aa_redun_syn_mat_prc_R, prob = seq(0,1,0.01),na.rm = T)
aa_redun_syn_mat_prc_R[is.na(aa_redun_syn_mat_prc_R)] <- 0
cols <- colorRampPalette(brewer.pal(9, "RdBu"))(length(aabreak)-1)
png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/FeatureFamily_synergy_random_prc.png",    # create PNG for the heat map        
    width = 10*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11)
heatmap.2((aa_redun_syn_mat_prc_R), 
          trace='none', Colv = T,
          dendrogram = "both",
          #ColSideColors = aa_col_color,
          col=rev(cols), 
          breaks = aabreak,
          margins = c(10,10))
dev.off()




# for auroc
aadefmat2 <- matrix(nrow = nrow(aa_auprc_base),
                    ncol = ncol(aa_auprc_base))
aadefmat2[] <- 0.5
aatstss31 <- marginal_perf_impro_pair(perf_list = aa_all_auroc_list_cv2TX_NoD,
                                     default_mat = aadefmat2,
                                     return_comp = T)

aaxna <- unique(aatstss31$marginal_perf_df$lncRNA)

for(i in 1:length(aaxna)){
  aamyd <- aatstss31$marginal_perf_df[aatstss31$marginal_perf_df$lncRNA %in% aaxna[i],]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/pair_per_comp_boxplot_auroc/", aaxna[i], ".png"), 
      width = 8*300,
      height = 8*300,
      res = 300,
      pointsize = 11)
  aaxplt <- ggplot(aamyd, aes(x=feature_family, y=perf_delta, fill = partition_type)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(aaxna[i])+
    ylab("delta perf")+
    xlab("feature family") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}
####################################
aamytss <- pair_indiv_compare(pair_df = aatstss31$marginal_perf_df, indiv_df = aatstss3$marginal_perf_df)
aamytss$diff <- aamytss$perf_pair - (aamytss$perf_1 + aamytss$perf_2)

aaxna <- unique(aamytss$lncRNA)

for(i in 1:length(aaxna)){
  aamyd <- aamytss[aamytss$lncRNA %in% aaxna[i],]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/pair_indDiff_per_comp_boxplot_auroc/", aaxna[i], ".png"), 
      width = 8*300,
      height = 8*300,
      res = 300,
      pointsize = 11)
  aaxplt <- ggplot(aamyd, aes(x=feature_family, y=diff, fill = partition_type)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(aaxna[i])+
    ylab("delta perf")+
    xlab("feature family") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}
# draw heatmap
aamytss_agg <- aggregate(x = aamytss[c("diff")],
                         by = aamytss[c("lncRNA","partition_type", "partition", "feature_family")], 
                         FUN = mean)
aamytss_agg2 <- aggregate(x = aamytss_agg[c("diff")],
                         by = aamytss_agg[c("lncRNA","partition_type", "feature_family")], 
                         FUN = mean)
aalncna <- unique(aamytss_agg2$lncRNA)
aapairr <- unique(aamytss_agg2$feature_family)
aa_redun_syn_mat_roc_C <- matrix(nrow = length(aalncna),
                               ncol = length(aapairr))
rownames(aa_redun_syn_mat_roc_C) <- aalncna
colnames(aa_redun_syn_mat_roc_C) <- aapairr
aa_redun_syn_mat_roc_R <- matrix(nrow = length(aalncna),
                                 ncol = length(aapairr))
rownames(aa_redun_syn_mat_roc_R) <- aalncna
colnames(aa_redun_syn_mat_roc_R) <- aapairr
for(i in 1:nrow(aa_redun_syn_mat_roc_C)){
  aa1 <- aamytss_agg2[((aamytss_agg2$partition_type == "chunk") & (aamytss_agg2$lncRNA == rownames(aa_redun_syn_mat_roc_C)[i])),]
  aa2 <- aamytss_agg2[((aamytss_agg2$partition_type == "random") & (aamytss_agg2$lncRNA == rownames(aa_redun_syn_mat_roc_C)[i])),]
  
  aa_redun_syn_mat_roc_C[i, ] <- aa1$diff[match(aa1$feature_family, colnames(aa_redun_syn_mat_roc_C))]
  aa_redun_syn_mat_roc_R[i, ] <- aa2$diff[match(aa2$feature_family, colnames(aa_redun_syn_mat_roc_R))]
}

aabreak <- quantile(aa_redun_syn_mat_roc_C, prob = seq(0,1,0.01),na.rm = T)
aa_redun_syn_mat_roc_C[is.na(aa_redun_syn_mat_roc_C)] <- 0
cols <- colorRampPalette(brewer.pal(9, "RdBu"))(length(aabreak)-1)
png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/FeatureFamily_synergy_chunk_roc.png",    # create PNG for the heat map        
    width = 10*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11)
heatmap.2((aa_redun_syn_mat_roc_C), 
          trace='none', Colv = T,
          dendrogram = "both",
          #ColSideColors = aa_col_color,
          col=rev(cols), 
          breaks = aabreak,
          margins = c(10,10))
dev.off()

aabreak <- quantile(aa_redun_syn_mat_roc_R, prob = seq(0,1,0.01),na.rm = T)
aa_redun_syn_mat_roc_R[is.na(aa_redun_syn_mat_roc_R)] <- 0
cols <- colorRampPalette(brewer.pal(9, "RdBu"))(length(aabreak)-1)
png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/FeatureFamily_synergy_random_roc.png",    # create PNG for the heat map        
    width = 10*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11)
heatmap.2((aa_redun_syn_mat_roc_R), 
          trace='none', Colv = T,
          dendrogram = "both",
          #ColSideColors = aa_col_color,
          col=rev(cols), 
          breaks = aabreak,
          margins = c(10,10))
dev.off()

####################################



aa2r <- aa_auroc_chunk2[,aa_col_con[1:5]]
aa2p <- aa_auprc_chunk2[,aa_col_con[1:5]]
aa1r <- aa_auroc_chunk2[,aa_col_seq]
aa1p <- aa_auprc_chunk2[,aa_col_seq]

aa2rrm <- apply(aa2r, 1, function(x) sum(x[x>0]))
aa2prm <- apply(aa2p, 1, function(x) sum(x[x>0]))

aa1rrm <- apply(aa1r, 1, function(x) sum(x[x>0]))
aa1prm <- apply(aa1p, 1, function(x) sum(x[x>0]))

cbind(lncRNA_mESC_vs_mOPC_sameBS_ratio, aa2rrm[match(names(lncRNA_mESC_vs_mOPC_sameBS_ratio), names(aa2rrm))])
plot(lncRNA_mESC_vs_mOPC_sameBS_ratio, aa2rrm[match(names(lncRNA_mESC_vs_mOPC_sameBS_ratio), names(aa2rrm))], ylab = "marginal increase context AUROC")
cbind(lncRNA_mESC_vs_mOPC_sameBS_ratio, aa2prm[match(names(lncRNA_mESC_vs_mOPC_sameBS_ratio), names(aa2prm))])
plot(lncRNA_mESC_vs_mOPC_sameBS_ratio, aa2prm[match(names(lncRNA_mESC_vs_mOPC_sameBS_ratio), names(aa2prm))], ylab = "marginal increase context AUPRC")

cor.test(lncRNA_mESC_vs_mOPC_sameBS_ratio,aa2prm[match(names(lncRNA_mESC_vs_mOPC_sameBS_ratio), names(aa2prm))], method = "spearman")
cor.test(lncRNA_mESC_vs_mOPC_sameBS_ratio,aa2rrm[match(names(lncRNA_mESC_vs_mOPC_sameBS_ratio), names(aa2rrm))], method = "spearman")

cbind(lncRNA_mESC_vs_mOPC_sameBS_ratio, aa1rrm[match(names(lncRNA_mESC_vs_mOPC_sameBS_ratio), names(aa1rrm))])
plot(lncRNA_mESC_vs_mOPC_sameBS_ratio, aa1rrm[match(names(lncRNA_mESC_vs_mOPC_sameBS_ratio), names(aa1rrm))])
cbind(lncRNA_mESC_vs_mOPC_sameBS_ratio, aa1prm[match(names(lncRNA_mESC_vs_mOPC_sameBS_ratio), names(aa1prm))])
plot(lncRNA_mESC_vs_mOPC_sameBS_ratio, aa1prm[match(names(lncRNA_mESC_vs_mOPC_sameBS_ratio), names(aa1prm))])

plot(aaperc2, aa2rrm[match(names(aaperc2), names(aa1rrm))])
plot(aaperc2, aa2prm[match(names(aaperc2), names(aa1prm))])

lncRNA_mESC_vs_mOPC_sameBS_ratio


plot(lncRNA_mESC_vs_mOPC_sameBS_ratio, (aa2prm / aa1prm)[match(names(lncRNA_mESC_vs_mOPC_sameBS_ratio), names(aa2prm))])
##############################
#get the performance difference between the best sequence and best context models for each lncRNA
aagghh <- grep("sequ", colnames(aa_all_auroc_rand_cv2TX_NoD))
aagghh2 <- grep("TXR", colnames(aa_all_auroc_rand_cv2TX_NoD))
aagghh3 <- union(aagghh, aagghh2)
  
aa_all_auroc_rand_cv2TX_noDist_sequ <- aa_all_auroc_rand_cv2TX_NoD[,-aagghh3]
aa_all_auroc_rand_cv2TX_noDist_cont <- aa_all_auroc_rand_cv2TX_NoD[,aagghh]

#aaxxt <- aa_all_auprc_rand_cv2TX_noDist_sequ
#rownames(aaxxt)[sort(apply(aaxxt,1,max), decreasing = T, index.return=T)$ix]

aa_all_auprc_rand_cv2TX_noDist_sequ <- aa_all_auprc_rand_cv2TX_NoD[,-aagghh3]
aa_all_auprc_rand_cv2TX_noDist_cont <- aa_all_auprc_rand_cv2TX_NoD[,aagghh]

aa_all_auroc_chun_cv2TX_noDist_sequ <- aa_all_auroc_chunk_cv2TX_NoD[,-aagghh3]
aa_all_auroc_chun_cv2TX_noDist_cont <- aa_all_auroc_chunk_cv2TX_NoD[,aagghh]

aa_all_auprc_chun_cv2TX_noDist_sequ <- aa_all_auprc_chunk_cv2TX_NoD[,-aagghh3]
aa_all_auprc_chun_cv2TX_noDist_cont <- aa_all_auprc_chunk_cv2TX_NoD[,aagghh]


aamax_cont_pr_chun <- apply(aa_all_auprc_chun_cv2TX_noDist_cont, 1, max, na.rm =T)
aamax_sequ_pr_chun <- apply(aa_all_auprc_chun_cv2TX_noDist_sequ, 1, max, na.rm =T)
aa1 <-aamax_cont_pr_chun- aamax_sequ_pr_chun
plot(lncRNA_mESC_vs_mOPC_sameBS_ratio, aa1[match(names(lncRNA_mESC_vs_mOPC_sameBS_ratio), names(aa1))], 
     ylab = " AUPRC improvement due to context", xlab = "Binding pattern conservation (mESC vs mOPC)", pch = 19, cex = 1)
cbind(lncRNA_mESC_vs_mOPC_sameBS_ratio, aa1[match(names(lncRNA_mESC_vs_mOPC_sameBS_ratio), names(aa1))])

cor.test(lncRNA_mESC_vs_mOPC_sameBS_ratio,aa1[match(names(lncRNA_mESC_vs_mOPC_sameBS_ratio), names(aa1))], method = "spearman")


aamax_cont_ro_chun <- apply(aa_all_auroc_chun_cv2TX_noDist_cont, 1, max, na.rm =T)
aamax_sequ_ro_chun <- apply(aa_all_auroc_chun_cv2TX_noDist_sequ, 1, max, na.rm =T)
aa2 <- aamax_cont_ro_chun - aamax_sequ_ro_chun
plot(lncRNA_mESC_vs_mOPC_sameBS_ratio, aa2[match(names(lncRNA_mESC_vs_mOPC_sameBS_ratio), names(aa2))])



aamax_cont_pr_rand <- apply(aa_all_auprc_rand_cv2TX_noDist_cont, 1, max, na.rm =T)
aamax_sequ_pr_rand <- apply(aa_all_auprc_rand_cv2TX_noDist_sequ, 1, max, na.rm =T)
aa3 <- aamax_cont_pr_rand - aamax_sequ_pr_rand
plot(lncRNA_mESC_vs_mOPC_sameBS_ratio, aa3[match(names(lncRNA_mESC_vs_mOPC_sameBS_ratio), names(aa3))])

aamax_cont_ro_rand <- apply(aa_all_auroc_rand_cv2TX_noDist_cont, 1, max, na.rm =T)
aamax_sequ_ro_rand <- apply(aa_all_auroc_rand_cv2TX_noDist_sequ, 1, max, na.rm =T)
aa4 <- aamax_cont_ro_rand  - aamax_sequ_ro_rand
plot(lncRNA_mESC_vs_mOPC_sameBS_ratio, aa4[match(names(lncRNA_mESC_vs_mOPC_sameBS_ratio), names(aa4))])



##############################
# feature importance analysis

aatop10perlnc_rand_sequ <- list()
aatop10perlnc_rand_cont <- list()

aatop10perlnc_chun_sequ <- list()
aatop10perlnc_chun_cont <- list()

aagghh <- grep("sequ", colnames(aa_all_auroc_rand_cv2TX_NoD))

aa_all_auroc_rand_cv2TX_noDist_sequ <- aa_all_auroc_rand_cv2TX_NoD[,-aagghh]
aa_all_auroc_rand_cv2TX_noDist_cont <- aa_all_auroc_rand_cv2TX_NoD[,aagghh]

#aaxxt <- aa_all_auprc_rand_cv2TX_noDist_sequ
#rownames(aaxxt)[sort(apply(aaxxt,1,max), decreasing = T, index.return=T)$ix]

aa_all_auprc_rand_cv2TX_noDist_sequ <- aa_all_auprc_rand_cv2TX_NoD[,-aagghh]
aa_all_auprc_rand_cv2TX_noDist_cont <- aa_all_auprc_rand_cv2TX_NoD[,aagghh]

aa_all_auroc_chun_cv2TX_noDist_sequ <- aa_all_auroc_chunk_cv2TX_NoD[,-aagghh]
aa_all_auroc_chun_cv2TX_noDist_cont <- aa_all_auroc_chunk_cv2TX_NoD[,aagghh]

aa_all_auprc_chun_cv2TX_noDist_sequ <- aa_all_auprc_chunk_cv2TX_NoD[,-aagghh]
aa_all_auprc_chun_cv2TX_noDist_cont <- aa_all_auprc_chunk_cv2TX_NoD[,aagghh]

aagghh1 <- grep("_TXRB", colnames(aa_all_auprc_rand_cv2TX_noDist_sequ))
aagghh2 <- grep("_TXRE", colnames(aa_all_auprc_rand_cv2TX_noDist_sequ))
aa_all_auroc_rand_cv2TX_noDist_sequ_noTX <- aa_all_auroc_rand_cv2TX_noDist_sequ[,-union(aagghh1,aagghh2)]
aa_all_auprc_rand_cv2TX_noDist_sequ_noTX <- aa_all_auprc_rand_cv2TX_noDist_sequ[,-union(aagghh1,aagghh2)]
aa_all_auroc_chun_cv2TX_noDist_sequ_noTX <- aa_all_auroc_chun_cv2TX_noDist_sequ[,-union(aagghh1,aagghh2)]
aa_all_auprc_chun_cv2TX_noDist_sequ_noTX <- aa_all_auprc_chun_cv2TX_noDist_sequ[,-union(aagghh1,aagghh2)]

for(i in 1:nrow(aa_all_auprc_chunk_cv2TX_NoD)){
  aamyind1RS <- quantile(aa_all_auroc_rand_cv2TX_noDist_sequ_noTX[i, ],  seq(0,1,0.1))
  aamyind2RS <- quantile(aa_all_auprc_rand_cv2TX_noDist_sequ_noTX[i, ],  seq(0,1,0.1))
  aamyind1CS <- quantile(aa_all_auroc_chun_cv2TX_noDist_sequ_noTX[i, ],  seq(0,1,0.1))
  aamyind2CS <- quantile(aa_all_auprc_chun_cv2TX_noDist_sequ_noTX[i, ],  seq(0,1,0.1))
  
  aamyind1RC <- quantile(aa_all_auroc_rand_cv2TX_noDist_cont[i, ],  seq(0,1,0.1), na.rm =T)
  aamyind2RC <- quantile(aa_all_auprc_rand_cv2TX_noDist_cont[i, ],  seq(0,1,0.1), na.rm =T)
  aamyind1CC <- quantile(aa_all_auroc_chun_cv2TX_noDist_cont[i, ],  seq(0,1,0.1), na.rm =T)
  aamyind2CC <- quantile(aa_all_auprc_chun_cv2TX_noDist_cont[i, ],  seq(0,1,0.1), na.rm =T)
  
  
  aatop10perlnc_rand_sequ[[i]] <- union(colnames(aa_all_auroc_rand_cv2TX_noDist_sequ_noTX)[aa_all_auroc_rand_cv2TX_noDist_sequ_noTX[i, ] >= aamyind1RS[9]],
                                        colnames(aa_all_auprc_rand_cv2TX_noDist_sequ_noTX)[aa_all_auprc_rand_cv2TX_noDist_sequ_noTX[i, ] >= aamyind2RS[9]])
  aatop10perlnc_chun_sequ[[i]] <- union(colnames(aa_all_auroc_chun_cv2TX_noDist_sequ_noTX)[aa_all_auroc_chun_cv2TX_noDist_sequ_noTX[i, ] >= aamyind1CS[9]],
                                        colnames(aa_all_auprc_chun_cv2TX_noDist_sequ_noTX)[aa_all_auprc_chun_cv2TX_noDist_sequ_noTX[i, ] >=  aamyind2CS[9]])
  
  aatop10perlnc_rand_cont[[i]] <- union(colnames(aa_all_auroc_rand_cv2TX_noDist_cont)[aa_all_auroc_rand_cv2TX_noDist_cont[i, ] >= aamyind1RC[9]],
                                        colnames(aa_all_auprc_rand_cv2TX_noDist_cont)[aa_all_auprc_rand_cv2TX_noDist_cont[i, ] >= aamyind2RC[9]])
  aatop10perlnc_chun_cont[[i]] <- union(colnames(aa_all_auroc_chun_cv2TX_noDist_cont)[aa_all_auroc_chun_cv2TX_noDist_cont[i, ] >= aamyind1CC[9]],
                                        colnames(aa_all_auprc_chun_cv2TX_noDist_cont)[aa_all_auprc_chun_cv2TX_noDist_cont[i, ] >=  aamyind2CC[9]])
  
}
names(aatop10perlnc_rand_sequ) <- rownames(aa_all_auprc_chunk_cv2TX_NoD)
names(aatop10perlnc_chun_sequ) <- rownames(aa_all_auprc_chunk_cv2TX_NoD)
names(aatop10perlnc_rand_cont) <- rownames(aa_all_auprc_chunk_cv2TX_NoD)
names(aatop10perlnc_chun_cont) <- rownames(aa_all_auprc_chunk_cv2TX_NoD)

#make a dataframe with the following columns: lncRNA, feature, model, feature_ranking
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/Partition_6_CVfirst_dataset_2410003L11Rik.RData")
aa_all_feat <- c(my_name_dic[,1], paste0("TR_sequ_", c(1:150)))


aaimpdf_rand_s <- list()
aaimpdf_rand_c <- list()
aaimpdf_chun_s <- list()
aaimpdf_chun_c <- list()

for(i in 1:length(aatop10perlnc_rand_sequ)){
  aaimpdf_rand_s[[i]] <- list()
  aaimpdf_rand_c[[i]] <- list()
  aaimpdf_chun_s[[i]] <- list()
  aaimpdf_chun_c[[i]] <- list()
  aacur_feat_c <- aa_all_impor_list_cv2TX_NoD[[i]][1:5]
  aacur_feat_r <- aa_all_impor_list_cv2TX_NoD[[i]][6:10]
  
  for(j in 1:length(aatop10perlnc_rand_sequ[[i]])){
    aatmplist <- list()
    aatmplist_rank <- list()
    for(k in 1:length(aacur_feat_r)){
      aayytmp <- aacur_feat_r[[k]][[match(aatop10perlnc_rand_sequ[[i]][j], 
                                            names( aacur_feat_r[[k]]))]]
      
      aatmplist[[k]] <- array(dim = length(aa_all_feat))
      names(aatmplist[[k]]) <- aa_all_feat
      aatmplist[[k]][match(names(aayytmp), names(aatmplist[[k]]))] <- aayytmp
      aatmplist_rank[[k]] <- array(dim = length(aa_all_feat))
      
      aasortind <- sort(aatmplist[[k]], decreasing = T, index.return = T, na.last =T)$ix
      asortind2 <- c(1:length(aasortind))[sort(aasortind, index.return = T)$ix]
      aatmplist_rank[[k]] <- asortind2
      aatmplist_rank[[k]][is.na(aatmplist[[k]])] <- NA
      names(aatmplist_rank[[k]]) <- aa_all_feat
      # if(k > 1){
      #   stopifnot(all(names(aatmplist[[k]]) == names(aatmplist[[k-1]])))
      # }
    }
    #stopifnot(length(unique(unlist(lapply(aatmplist, length)))) == 1)
    aatmplist_df <- do.call(rbind, aatmplist_rank)
    aatmplist_df_mean <- apply(aatmplist_df, 2, median, na.rm = F) #colMeans(aatmplist_df, na.rm = F)
    #aasortind <- sort(aatmplist_df_mean, decreasing = T, index.return = T, na.last =T)$ix
    #aasortind2 <- c(1:length(aasortind))[sort(aasortind, index.return = T)$ix]
    #aasortind2[is.na(aatmplist_df_mean)] <- NA
    aaimpdf_rand_s[[i]][[j]] <- data.frame(lncRNA = rep(names(aatop10perlnc_rand_sequ)[i], length(aatmplist_df_mean)),
                                           model =  rep(aatop10perlnc_rand_sequ[[i]][j], length(aatmplist_df_mean)),
                                           feature = names(aatmplist_df_mean),
                                           rank = aatmplist_df_mean)
  }
  
  for(j in 1:length(aatop10perlnc_rand_cont[[i]])){
    aatmplist <- list()
    aatmplist_rank <- list()
    aatmplist_noTR <- list()
    aatmplist_rank_noTR <- list()
    
    for(k in 1:length(aacur_feat_r)){
      aayytmp <- aacur_feat_r[[k]][[match(aatop10perlnc_rand_cont[[i]][j], 
                                          names( aacur_feat_r[[k]]))]]
      aayytmp_noTR <- aayytmp[-grep("TR_sequ", names(aayytmp))]
      aatmplist[[k]] <- array(dim = length(aa_all_feat))
      names(aatmplist[[k]]) <- aa_all_feat
      
      aatmplist_noTR[[k]] <- array(dim = length(aa_all_feat))
      names(aatmplist_noTR[[k]]) <- aa_all_feat
      
      
      aatmplist[[k]][match(names(aayytmp), names(aatmplist[[k]]))] <- aayytmp
      aatmplist_noTR[[k]][match(names(aayytmp_noTR), names(aatmplist_noTR[[k]]))] <- aayytmp_noTR
      
      
      aatmplist_rank[[k]] <- array(dim = length(aa_all_feat))
      aatmplist_rank_noTR[[k]] <- array(dim = length(aa_all_feat))
      
      aasortind <- sort(aatmplist[[k]], decreasing = T, index.return = T, na.last =T)$ix
      asortind2 <- c(1:length(aasortind))[sort(aasortind, index.return = T)$ix]
      aatmplist_rank[[k]] <- asortind2
      aatmplist_rank[[k]][is.na(aatmplist[[k]])] <- NA
      names(aatmplist_rank[[k]]) <- aa_all_feat
     
      aasortind <- sort(aatmplist_noTR[[k]], decreasing = T, index.return = T, na.last =T)$ix
      asortind2 <- c(1:length(aasortind))[sort(aasortind, index.return = T)$ix]
      aatmplist_rank_noTR[[k]] <- asortind2
      aatmplist_rank_noTR[[k]][is.na(aatmplist_noTR[[k]])] <- NA
      names(aatmplist_rank_noTR[[k]]) <- aa_all_feat
      
      
      # if(k > 1){
      #   stopifnot(all(names(aatmplist[[k]]) == names(aatmplist[[k-1]])))
      # }
    }
    #stopifnot(length(unique(unlist(lapply(aatmplist, length)))) == 1)
    aatmplist_df <- do.call(rbind, aatmplist_rank)
    aatmplist_df_mean <- apply(aatmplist_df, 2, median, na.rm = F) #colMeans(aatmplist_df, na.rm = F)
    
    aatmplist_df_noTR <- do.call(rbind, aatmplist_rank_noTR)
    aatmplist_df_mean_noTR <- apply(aatmplist_df_noTR, 2, median, na.rm = F) #colMeans(aatmplist_df, na.rm = F)
    stopifnot(all(names(aatmplist_df_mean) == names(aatmplist_df_mean_noTR)) )
    #aasortind <- sort(aatmplist_df_mean, decreasing = T, index.return = T, na.last =T)$ix
    #aasortind2 <- c(1:length(aasortind))[sort(aasortind, index.return = T)$ix]
    #aasortind2[is.na(aatmplist_df_mean)] <- NA
    aaimpdf_rand_c[[i]][[j]] <- data.frame(lncRNA = rep(names(aatop10perlnc_rand_cont)[i], length(aatmplist_df_mean)),
                                           model =  rep(aatop10perlnc_rand_cont[[i]][j], length(aatmplist_df_mean)),
                                           feature = names(aatmplist_df_mean),
                                           rank = aatmplist_df_mean, 
                                           rank_WOTR = aatmplist_df_mean_noTR)
  }
  
  
  for(j in 1:length(aatop10perlnc_chun_sequ[[i]])){
    aatmplist <- list()
    aatmplist_rank <- list()
    for(k in 1:length(aacur_feat_c)){
      aayytmp <- aacur_feat_c[[k]][[match(aatop10perlnc_chun_sequ[[i]][j], 
                                          names( aacur_feat_c[[k]]))]]
      
      aatmplist[[k]] <- array(dim = length(aa_all_feat))
      names(aatmplist[[k]]) <- aa_all_feat
      aatmplist[[k]][match(names(aayytmp), names(aatmplist[[k]]))] <- aayytmp
      aatmplist_rank[[k]] <- array(dim = length(aa_all_feat))
      
      aasortind <- sort(aatmplist[[k]], decreasing = T, index.return = T, na.last =T)$ix
      asortind2 <- c(1:length(aasortind))[sort(aasortind, index.return = T)$ix]
      aatmplist_rank[[k]] <- asortind2
      aatmplist_rank[[k]][is.na(aatmplist[[k]])] <- NA
      names(aatmplist_rank[[k]]) <- aa_all_feat
      # if(k > 1){
      #   stopifnot(all(names(aatmplist[[k]]) == names(aatmplist[[k-1]])))
      # }
    }
    #stopifnot(length(unique(unlist(lapply(aatmplist, length)))) == 1)
    aatmplist_df <- do.call(rbind, aatmplist_rank)
    aatmplist_df_mean <- apply(aatmplist_df, 2, median, na.rm = F) #colMeans(aatmplist_df, na.rm = F)
    #aasortind <- sort(aatmplist_df_mean, decreasing = T, index.return = T, na.last =T)$ix
    #aasortind2 <- c(1:length(aasortind))[sort(aasortind, index.return = T)$ix]
    #aasortind2[is.na(aatmplist_df_mean)] <- NA
    aaimpdf_chun_s[[i]][[j]] <- data.frame(lncRNA = rep(names(aatop10perlnc_chun_sequ)[i], length(aatmplist_df_mean)),
                                           model =  rep(aatop10perlnc_chun_sequ[[i]][j], length(aatmplist_df_mean)),
                                           feature = names(aatmplist_df_mean),
                                           rank = aatmplist_df_mean)
  }
  
  
  for(j in 1:length(aatop10perlnc_chun_cont[[i]])){
    aatmplist <- list()
    aatmplist_rank <- list()
    aatmplist_noTR <- list()
    aatmplist_rank_noTR <- list()
    
    for(k in 1:length(aacur_feat_c)){
      aayytmp <- aacur_feat_c[[k]][[match(aatop10perlnc_chun_cont[[i]][j], 
                                          names( aacur_feat_c[[k]]))]]
      aayytmp_noTR <- aayytmp[-grep("TR_sequ", names(aayytmp))]
      
      aatmplist[[k]] <- array(dim = length(aa_all_feat))
      names(aatmplist[[k]]) <- aa_all_feat
      
      aatmplist_noTR[[k]] <- array(dim = length(aa_all_feat))
      names(aatmplist_noTR[[k]]) <- aa_all_feat
      
      aatmplist[[k]][match(names(aayytmp), names(aatmplist[[k]]))] <- aayytmp
      aatmplist_noTR[[k]][match(names(aayytmp_noTR), names(aatmplist_noTR[[k]]))] <- aayytmp_noTR
      
      aatmplist_rank[[k]] <- array(dim = length(aa_all_feat))
      aatmplist_rank_noTR[[k]] <- array(dim = length(aa_all_feat))
      
      
      aasortind <- sort(aatmplist[[k]], decreasing = T, index.return = T, na.last =T)$ix
      asortind2 <- c(1:length(aasortind))[sort(aasortind, index.return = T)$ix]
      aatmplist_rank[[k]] <- asortind2
      aatmplist_rank[[k]][is.na(aatmplist[[k]])] <- NA
      names(aatmplist_rank[[k]]) <- aa_all_feat
      
      aasortind <- sort(aatmplist_noTR[[k]], decreasing = T, index.return = T, na.last =T)$ix
      asortind2 <- c(1:length(aasortind))[sort(aasortind, index.return = T)$ix]
      aatmplist_rank_noTR[[k]] <- asortind2
      aatmplist_rank_noTR[[k]][is.na(aatmplist_noTR[[k]])] <- NA
      names(aatmplist_rank_noTR[[k]]) <- aa_all_feat
      
      # if(k > 1){
      #   stopifnot(all(names(aatmplist[[k]]) == names(aatmplist[[k-1]])))
      # }
    }
    #stopifnot(length(unique(unlist(lapply(aatmplist, length)))) == 1)
    aatmplist_df <- do.call(rbind, aatmplist_rank)
    aatmplist_df_mean <- apply(aatmplist_df, 2, median, na.rm = F) #colMeans(aatmplist_df, na.rm = F)
    
    aatmplist_df_noTR <- do.call(rbind, aatmplist_rank_noTR)
    aatmplist_df_mean_noTR <- apply(aatmplist_df_noTR, 2, median, na.rm = F) #colMeans(aatmplist_df, na.rm = F)
    stopifnot(all(names(aatmplist_df_mean) == names(aatmplist_df_mean_noTR)) )
    
    # aasortind <- sort(aatmplist_df_mean, decreasing = T, index.return = T, na.last =T)$ix
    # aasortind2 <- c(1:length(aasortind))[sort(aasortind, index.return = T)$ix]
    # aasortind2[is.na(aatmplist_df_mean)] <- NA
    aaimpdf_chun_c[[i]][[j]] <- data.frame(lncRNA = rep(names(aatop10perlnc_chun_cont)[i], length(aatmplist_df_mean)),
                                           model =  rep(aatop10perlnc_chun_cont[[i]][j], length(aatmplist_df_mean)),
                                           feature = names(aatmplist_df_mean),
                                           rank = aatmplist_df_mean, 
                                           rank_WOTR = aatmplist_df_mean_noTR)
  }
  # aatmp1
  # aatmp2
  # aatmp3
  # aatmp4
  
}

aalncFeat_seq_R_list <- list()
aalncFeat_seq_C_list <- list()

aalncFeat_con_R_list <- list()
aalncFeat_con_C_list <- list()

for(i in 1:length(aaimpdf_chun_c)){
  aalncFeat_seq_R_list[[i]] <- do.call(rbind, aaimpdf_rand_s[[i]])
  aalncFeat_seq_C_list[[i]] <- do.call(rbind, aaimpdf_chun_s[[i]])
  aalncFeat_con_R_list[[i]] <- do.call(rbind, aaimpdf_rand_c[[i]])
  aalncFeat_con_C_list[[i]] <- do.call(rbind, aaimpdf_chun_c[[i]])
}

aalncFeat_seq_R_df <- do.call(rbind, aalncFeat_seq_R_list)
aalncFeat_seq_C_df <- do.call(rbind, aalncFeat_seq_C_list)
aalncFeat_con_R_df <- do.call(rbind, aalncFeat_con_R_list)
aalncFeat_con_C_df <- do.call(rbind, aalncFeat_con_C_list)



lncRNA_feature_seq_rand_summary_list <- list()
for(i in 1:length(aalncFeat_seq_R_list)){
  print(i)
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_violin/Random/Seq/rand_seq_",names(aatop10perlnc_rand_sequ)[i], ".png"),    # create PNG for the heat map        
      width = 12*300,        # 5 x 300 pixels
      height = 8*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  
  aatgc <- summarySE(aalncFeat_seq_R_list[[i]], measurevar=c("rank"), groupvars=c("feature"), na.rm=T)
  lncRNA_feature_seq_rand_summary_list[[i]] <- aatgc
  aatgc <- aatgc[aatgc$N >0 ,]
  aatgc2 <- aatgc[sort(aatgc$median, index.return=T)$ix,]
  aaplx <- aalncFeat_seq_R_list[[i]]
  aaplx$feature <- factor(aalncFeat_seq_R_list[[i]]$feature, levels = aatgc2$feature)
  aaxplt <- ggplot(aaplx[aalncFeat_seq_R_list[[i]]$feature %in% aatgc2$feature[1:100],], aes(x=feature, y=rank)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),scale = "count") +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(paste0("ranking of sequence features for ",names(aatop10perlnc_rand_sequ)[i]))+
    ylab("rank")+
    xlab("seq feature") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}

lncRNA_feature_seq_chunk_summary_list <- list()
for(i in 1:length(aalncFeat_seq_C_list)){
  print(i)
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_violin/Chunk/Seq/chunk_seq_",
                        names(aatop10perlnc_chun_sequ)[i], ".png"),    # create PNG for the heat map        
      width = 12*300,        # 5 x 300 pixels
      height = 8*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  
  aatgc <- summarySE(aalncFeat_seq_C_list[[i]], measurevar=c("rank"), groupvars=c("feature"), na.rm=T)
  lncRNA_feature_seq_chunk_summary_list[[i]] <- aatgc
  aatgc <- aatgc[aatgc$N >0 ,]
  aatgc2 <- aatgc[sort(aatgc$median, index.return=T)$ix,]
  aaplx <- aalncFeat_seq_C_list[[i]]
  aaplx$feature <- factor(aalncFeat_seq_C_list[[i]]$feature, levels = aatgc2$feature)
  aaxplt <- ggplot(aaplx[aalncFeat_seq_C_list[[i]]$feature %in% aatgc2$feature[1:100],], aes(x=feature, y=rank)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),scale = "count") +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(paste0("ranking of sequence features for ",names(aatop10perlnc_chun_sequ)[i]))+
    ylab("rank")+
    xlab("seq feature") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}


lncRNA_feature_context_rand_summary_list <- list()
lncRNA_feature_context_rand_summary_list_NOTR <- list()

for(i in 1:length(aalncFeat_con_R_list)){
  print(i)
  # png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_violin/Random/Context/rand_context_",names(aatop10perlnc_rand_sequ)[i], ".png"),    # create PNG for the heat map
  #     width = 18*300,        # 5 x 300 pixels
  #     height = 8*300,
  #     res = 300,            # 300 pixels per inch
  #     pointsize = 11)

  aatgc <- summarySE(aalncFeat_con_R_list[[i]], measurevar=c("rank"), groupvars=c("feature"), na.rm=T)
  lncRNA_feature_context_rand_summary_list[[i]] <- aatgc
  # aatgc <- aatgc[aatgc$N >0 ,]
  # aatgc2 <- aatgc[sort(aatgc$median, index.return=T)$ix,]
  # aaplx <- aalncFeat_con_R_list[[i]]
  # aaplx$feature <- factor(aalncFeat_con_R_list[[i]]$feature, levels = aatgc2$feature)
  # aaxplt <- ggplot(aaplx[aalncFeat_con_R_list[[i]]$feature %in% aatgc2$feature[1:200],], aes(x=feature, y=rank)) +
  #   # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
  #   #               labels = trans_format("log10", math_format(10^.x))) +
  #   geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),scale = "count") +
  #   #ggtitle("all points (2116029) distance to owner origin")+
  #   ggtitle(paste0("ranking of context features for ",names(aatop10perlnc_rand_sequ)[i]))+
  #   ylab("rank")+
  #   xlab("context feature") +
  #   theme(plot.title = element_text(hjust = 0.5),
  #         panel.grid.major.x = element_blank(),
  #         axis.text = element_text(angle = 90))
  # print(aaxplt)
  # dev.off()
  
  ### NOTR
  # png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_violin/Random/Context/rand_context_noTR_",names(aatop10perlnc_rand_sequ)[i], ".png"),    # create PNG for the heat map        
  #     width = 18*300,        # 5 x 300 pixels
  #     height = 8*300,
  #     res = 300,            # 300 pixels per inch
  #     pointsize = 11)
  
  aatgc <- summarySE(aalncFeat_con_R_list[[i]], measurevar=c("rank_WOTR"), groupvars=c("feature"), na.rm=T)
  lncRNA_feature_context_rand_summary_list_NOTR[[i]] <- aatgc
  # aatgc <- aatgc[aatgc$N >0 ,]
  # aatgc2 <- aatgc[sort(aatgc$median, index.return=T)$ix,]
  # aaplx <- aalncFeat_con_R_list[[i]]
  # aaplx$feature <- factor(aalncFeat_con_R_list[[i]]$feature, levels = aatgc2$feature)
  # aaxplt <- ggplot(aaplx[aalncFeat_con_R_list[[i]]$feature %in% aatgc2$feature[1:200],], aes(x=feature, y=rank_WOTR)) + 
  #   # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
  #   #               labels = trans_format("log10", math_format(10^.x))) +
  #   geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),scale = "count") +
  #   #ggtitle("all points (2116029) distance to owner origin")+
  #   ggtitle(paste0("ranking of context features for ",names(aatop10perlnc_rand_sequ)[i]))+
  #   ylab("rank")+
  #   xlab("context feature") + 
  #   theme(plot.title = element_text(hjust = 0.5),
  #         panel.grid.major.x = element_blank(),
  #         axis.text = element_text(angle = 90))
  # print(aaxplt)
  # dev.off()
}

lncRNA_feature_context_chunk_summary_list <- list()
lncRNA_feature_context_chunk_summary_list_NoTR <- list()

for(i in 1:length(aalncFeat_con_C_list)){
  print(i)
  # png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_violin/Chunk/Context/chunk_context_",
  #                       names(aatop10perlnc_chun_sequ)[i], ".png"),    # create PNG for the heat map        
  #     width = 18*300,        # 5 x 300 pixels
  #     height = 8*300,
  #     res = 300,            # 300 pixels per inch
  #     pointsize = 11)
  # 
  # aatgc <- summarySE(aalncFeat_con_C_list[[i]], measurevar=c("rank"), groupvars=c("feature"), na.rm=T)
  # lncRNA_feature_context_chunk_summary_list[[i]] <- aatgc
  # aatgc <- aatgc[aatgc$N >0 ,]
  # aatgc2 <- aatgc[sort(aatgc$median, index.return=T)$ix,]
  # aaplx <- aalncFeat_con_C_list[[i]]
  # aaplx$feature <- factor(aalncFeat_con_C_list[[i]]$feature, levels = aatgc2$feature)
  # aaxplt <- ggplot(aaplx[aalncFeat_con_C_list[[i]]$feature %in% aatgc2$feature[1:200],], aes(x=feature, y=rank)) + 
  #   # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
  #   #               labels = trans_format("log10", math_format(10^.x))) +
  #   geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),scale = "count") +
  #   #ggtitle("all points (2116029) distance to owner origin")+
  #   ggtitle(paste0("ranking of context features for ",names(aatop10perlnc_chun_sequ)[i]))+
  #   ylab("rank")+
  #   xlab("context feature") + 
  #   theme(plot.title = element_text(hjust = 0.5),
  #         panel.grid.major.x = element_blank(),
  #         axis.text = element_text(angle = 90))
  # print(aaxplt)
  # dev.off()
  
  #####NOTR
  print(i)
  # png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_violin/Chunk/Context/chunk_context_",
  #                       names(aatop10perlnc_chun_sequ)[i], ".png"),    # create PNG for the heat map        
  #     width = 18*300,        # 5 x 300 pixels
  #     height = 8*300,
  #     res = 300,            # 300 pixels per inch
  #     pointsize = 11)
  
  aatgc <- summarySE(aalncFeat_con_C_list[[i]], measurevar=c("rank_WOTR"), groupvars=c("feature"), na.rm=T)
  lncRNA_feature_context_chunk_summary_list_NoTR[[i]] <- aatgc
  # aatgc <- aatgc[aatgc$N >0 ,]
  # aatgc2 <- aatgc[sort(aatgc$median, index.return=T)$ix,]
  # aaplx <- aalncFeat_con_C_list[[i]]
  # aaplx$feature <- factor(aalncFeat_con_C_list[[i]]$feature, levels = aatgc2$feature)
  # aaxplt <- ggplot(aaplx[aalncFeat_con_C_list[[i]]$feature %in% aatgc2$feature[1:200],], aes(x=feature, y=rank_WOTR)) + 
  #   # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
  #   #               labels = trans_format("log10", math_format(10^.x))) +
  #   geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),scale = "count") +
  #   #ggtitle("all points (2116029) distance to owner origin")+
  #   ggtitle(paste0("ranking of context features for ",names(aatop10perlnc_chun_sequ)[i]))+
  #   ylab("rank")+
  #   xlab("context feature") + 
  #   theme(plot.title = element_text(hjust = 0.5),
  #         panel.grid.major.x = element_blank(),
  #         axis.text = element_text(angle = 90))
  # print(aaxplt)
  # dev.off()
}

#########
# make boxplot of the same plots
#lncRNA_feature_seq_rand_summary_list <- list()
for(i in 1:length(aalncFeat_seq_R_list)){
  print(i)
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_violin/Random/Seq/rand_seq_",
                        names(aatop10perlnc_rand_sequ)[i], "_boxplot.png"),    # create PNG for the heat map        
      width = 12*300,        # 5 x 300 pixels
      height = 8*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  
  aatgc <- summarySE(aalncFeat_seq_R_list[[i]], measurevar=c("rank"), groupvars=c("feature"), na.rm=T)
  #lncRNA_feature_seq_rand_summary_list[[i]] <- aatgc
  aatgc <- aatgc[aatgc$N >0 ,]
  aatgc2 <- aatgc[sort(aatgc$median, index.return=T)$ix,]
  aaplx <- aalncFeat_seq_R_list[[i]]
  aaplx$feature <- factor(aalncFeat_seq_R_list[[i]]$feature, levels = aatgc2$feature)
  aaxplt <- ggplot(aaplx[aalncFeat_seq_R_list[[i]]$feature %in% aatgc2$feature[1:100],], aes(x=feature, y=rank)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(paste0("ranking of sequence features for ",names(aatop10perlnc_rand_sequ)[i]))+
    ylab("rank")+
    xlab("seq feature") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}

#lncRNA_feature_seq_chunk_summary_list <- list()
for(i in 1:length(aalncFeat_seq_C_list)){
  print(i)
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_violin/Chunk/Seq/chunk_seq_",
                        names(aatop10perlnc_chun_sequ)[i], "_boxplot.png"),    # create PNG for the heat map        
      width = 12*300,        # 5 x 300 pixels
      height = 8*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  
  aatgc <- summarySE(aalncFeat_seq_C_list[[i]], measurevar=c("rank"), groupvars=c("feature"), na.rm=T)
  #lncRNA_feature_seq_chunk_summary_list[[i]] <- aatgc
  aatgc <- aatgc[aatgc$N >0 ,]
  aatgc2 <- aatgc[sort(aatgc$median, index.return=T)$ix,]
  aaplx <- aalncFeat_seq_C_list[[i]]
  aaplx$feature <- factor(aalncFeat_seq_C_list[[i]]$feature, levels = aatgc2$feature)
  aaxplt <- ggplot(aaplx[aalncFeat_seq_C_list[[i]]$feature %in% aatgc2$feature[1:100],], aes(x=feature, y=rank)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(paste0("ranking of sequence features for ",names(aatop10perlnc_chun_sequ)[i]))+
    ylab("rank")+
    xlab("seq feature") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}


#lncRNA_feature_context_rand_summary_list <- list()
for(i in 1:length(aalncFeat_con_R_list)){
  print(i)
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_violin/Random/Context/rand_context_",
                        names(aatop10perlnc_rand_sequ)[i], "_boxplot.png"),    # create PNG for the heat map        
      width = 18*300,        # 5 x 300 pixels
      height = 8*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  
  aatgc <- summarySE(aalncFeat_con_R_list[[i]], measurevar=c("rank"), groupvars=c("feature"), na.rm=T)
  # lncRNA_feature_context_rand_summary_list[[i]] <- aatgc
  aatgc <- aatgc[aatgc$N >0 ,]
  aatgc2 <- aatgc[sort(aatgc$median, index.return=T)$ix,]
  aaplx <- aalncFeat_con_R_list[[i]]
  aaplx$feature <- factor(aalncFeat_con_R_list[[i]]$feature, levels = aatgc2$feature)
  aaxplt <- ggplot(aaplx[aalncFeat_con_R_list[[i]]$feature %in% aatgc2$feature[1:200],], aes(x=feature, y=rank)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(paste0("ranking of context features for ",names(aatop10perlnc_rand_sequ)[i]))+
    ylab("rank")+
    xlab("context feature") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}

#lncRNA_feature_context_chunk_summary_list <- list()
for(i in 1:length(aalncFeat_con_C_list)){
  print(i)
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_violin/Chunk/Context/chunk_context_",
                        names(aatop10perlnc_chun_sequ)[i], "_boxplot.png"),    # create PNG for the heat map        
      width = 18*300,        # 5 x 300 pixels
      height = 8*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  
  aatgc <- summarySE(aalncFeat_con_C_list[[i]], measurevar=c("rank"), groupvars=c("feature"), na.rm=T)
  #  lncRNA_feature_context_chunk_summary_list[[i]] <- aatgc
  aatgc <- aatgc[aatgc$N >0 ,]
  aatgc2 <- aatgc[sort(aatgc$median, index.return=T)$ix,]
  aaplx <- aalncFeat_con_C_list[[i]]
  aaplx$feature <- factor(aalncFeat_con_C_list[[i]]$feature, levels = aatgc2$feature)
  aaxplt <- ggplot(aaplx[aalncFeat_con_C_list[[i]]$feature %in% aatgc2$feature[1:200],], aes(x=feature, y=rank)) + 
    # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
    #               labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() +
    #ggtitle("all points (2116029) distance to owner origin")+
    ggtitle(paste0("ranking of context features for ",names(aatop10perlnc_chun_sequ)[i]))+
    ylab("rank")+
    xlab("context feature") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(angle = 90))
  print(aaxplt)
  dev.off()
}

##########################################################################################
##############################################################################
# draw heatmaps using the median value



lncRNA_feature_seq_rand_summary_list
lncRNA_feature_seq_chunk_summary_list
lncRNA_feature_context_rand_summary_list
lncRNA_feature_context_chunk_summary_list


#unlist(lapply(lncRNA_feature_context_rand_summary_list, nrow))
all(lncRNA_feature_context_rand_summary_list[[1]]$feature == lncRNA_feature_context_rand_summary_list[[20]]$feature)



aaxlsi <- lapply(lncRNA_feature_context_rand_summary_list, "[[", "median")
aaxlsi <- do.call(rbind, aaxlsi)
colnames(aaxlsi) <- lncRNA_feature_context_rand_summary_list[[1]]$feature
rownames(aaxlsi) <- rownames(aa_all_auprc_chunk_cv2TX_noDist)

aaxlsi <- aaxlsi[,!(colSums(is.na(aaxlsi)) == nrow(aaxlsi))]
aamin <- apply(aaxlsi, 2, min, na.rm=T)
aaquant <- quantile(aamin, prob = seq(0,1,0.1), na.rm=T)

aaxlsi <- aaxlsi[, aamin <= aaquant[6]]
aaarr <- range(aaxlsi, na.rm = T)
#aabreak <- c(0,5,10,20,50,75,100,155,200,250,aaarr[2])
aapostname <- unlist(lapply(strsplit(colnames(aaxlsi), "__"), function(x) x[[length(x)]]))
aapostname2 <- strsplit(aapostname, "_")
aapostname2_2 <- unlist(lapply(aapostname2, function(x) x[[min(2, length(x))]]))
table(aapostname2_2)
aa_co_sort <- sort(aapostname2_2, index.return=T)$ix
aa_col_color <- brewer.pal(9, "Pastel1")[factor(aapostname2_2[aa_co_sort])]
#aapostname2_2[aa_co_sort]

aabreak <- quantile(aaxlsi, prob = seq(0,1,0.01),na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(length(aabreak)-1)
png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_heatmap_Context_rand.png",    # create PNG for the heat map        
    width = 30*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11)
heatmap.2((aaxlsi[,aa_co_sort]), 
          trace='none', Colv = F,
          dendrogram = "none",
          Rowv = T,
          ColSideColors = aa_col_color,
          col=rev(cols), 
          breaks = aabreak,
          margins = c(12,8),
          density.info="none",
          key.title = "",
          key.xlab = "",
          key.ylab = "",
          keysize = 0.75)
dev.off()




aaxlsi <- aaxlsi[, -grep("TR_sequ", colnames(aaxlsi))]
aapostname <- unlist(lapply(strsplit(colnames(aaxlsi), "__"), function(x) x[[length(x)]]))
aapostname2 <- strsplit(aapostname, "_")
aapostname2_2 <- unlist(lapply(aapostname2, function(x) x[[min(2, length(x))]]))
table(aapostname2_2)
aa_co_sort <- sort(aapostname2_2, index.return=T)$ix
aa_col_color <- brewer.pal(9, "Pastel1")[factor(aapostname2_2[aa_co_sort])]

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_heatmap_Context_rand_noSEQU.png",    # create PNG for the heat map        
    width = 12*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 12)
heatmap.2((aaxlsi[,aa_co_sort]), 
          trace='none', Colv = F,
          dendrogram = "none",Rowv = T,
          ColSideColors = aa_col_color,
          col=rev(cols), 
          breaks = aabreak,
          margins = c(13,8),
          density.info="none",
          key.title = "",
          key.xlab = "",
          key.ylab = "",
          keysize = 0.75)
dev.off()



#################

aaxlsi <- lapply(lncRNA_feature_context_chunk_summary_list  , "[[", "median")
aaxlsi <- do.call(rbind, aaxlsi)
colnames(aaxlsi) <- lncRNA_feature_context_chunk_summary_list[[1]]$feature
rownames(aaxlsi) <- rownames(aa_all_auprc_chunk_cv2TX_noDist)

aaxlsi <- aaxlsi[,!(colSums(is.na(aaxlsi)) == nrow(aaxlsi))]
aamin <- apply(aaxlsi, 2, min, na.rm=T)
aaquant <- quantile(aamin, prob = seq(0,1,0.1), na.rm=T)
aaxlsi <- aaxlsi[, aamin <= aaquant[6]]
aaarr <- range(aaxlsi, na.rm = T)
#aabreak <- c(0,5,10,20,50,75,100,155,200,250,aaarr[2])
aapostname <- unlist(lapply(strsplit(colnames(aaxlsi), "__"), function(x) x[[length(x)]]))
aapostname2 <- strsplit(aapostname, "_")
aapostname2_2 <- unlist(lapply(aapostname2, function(x) x[[min(2, length(x))]]))
table(aapostname2_2)
aa_co_sort <- sort(aapostname2_2, index.return=T)$ix
aa_col_color <- brewer.pal(9, "Pastel1")[factor(aapostname2_2[aa_co_sort])]
#aapostname2_2[aa_co_sort]

aabreak <- quantile(aaxlsi, prob = seq(0,1,0.01),na.rm = T)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(length(aabreak)-1)
png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_heatmap_Context_chunk.png",    # create PNG for the heat map        
    width = 30*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11)
heatmap.2((aaxlsi[,aa_co_sort]), 
          trace='none', Colv = F,
          dendrogram = "none",Rowv = T,
          ColSideColors = aa_col_color,
          col=rev(cols), 
          breaks = aabreak,
          margins = c(12,8),
          density.info="none",
          key.title = "",
          key.xlab = "",
          key.ylab = "",
          keysize = 0.75)
dev.off()

aaxlsi <- aaxlsi[, -grep("TR_sequ", colnames(aaxlsi))]
aapostname <- unlist(lapply(strsplit(colnames(aaxlsi), "__"), function(x) x[[length(x)]]))
aapostname2 <- strsplit(aapostname, "_")
aapostname2_2 <- unlist(lapply(aapostname2, function(x) x[[min(2, length(x))]]))
table(aapostname2_2)
aa_co_sort <- sort(aapostname2_2, index.return=T)$ix
aa_col_color <- brewer.pal(9, "Pastel1")[factor(aapostname2_2[aa_co_sort])]
#aapostname2_2[aa_co_sort]


png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_heatmap_Context_chunk_noSEQU.png",    # create PNG for the heat map        
    width = 12*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 12)
heatmap.2((aaxlsi[,aa_co_sort]), 
          trace='none', Colv = F,Rowv = T,
          dendrogram = "none",
          ColSideColors = aa_col_color,
          col=rev(cols), 
          breaks = aabreak,
          margins = c(12,8),   
          density.info="none",
          key.title = "",
          key.xlab = "",
          key.ylab = "",
          keysize = 0.75)
dev.off()
#################3 NOTR

#################

aaxlsi <- lapply(lncRNA_feature_context_chunk_summary_list_NoTR  , "[[", "median")
aaxlsi <- do.call(rbind, aaxlsi)
colnames(aaxlsi) <- lncRNA_feature_context_chunk_summary_list_NoTR[[1]]$feature
rownames(aaxlsi) <- rownames(aa_all_auprc_chunk_cv2TX_noDist)

aaxlsi <- aaxlsi[,!(colSums(is.na(aaxlsi)) == nrow(aaxlsi))]
aamin <- apply(aaxlsi, 2, min, na.rm=T)
aaquant <- quantile(aamin, prob = seq(0,1,0.1), na.rm=T)
aaxlsi <- aaxlsi[, aamin <= aaquant[6]]
aaarr <- range(aaxlsi, na.rm = T)
#aabreak <- c(0,5,10,20,50,75,100,155,200,250,aaarr[2])
aapostname <- unlist(lapply(strsplit(colnames(aaxlsi), "__"), function(x) x[[length(x)]]))
aapostname2 <- strsplit(aapostname, "_")
aapostname2_2 <- unlist(lapply(aapostname2, function(x) x[[min(2, length(x))]]))
table(aapostname2_2)
aa_co_sort <- sort(aapostname2_2, index.return=T)$ix
aa_col_color <- brewer.pal(9, "Pastel1")[factor(aapostname2_2[aa_co_sort])]
#aapostname2_2[aa_co_sort]

aabreak <- quantile(aaxlsi, prob = seq(0,1,0.01),na.rm = T)
aarng <- range(aaxlsi, na.rm = T)
aabreak <- c(seq(1,100,5), c(150,200,300))

cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(length(aabreak)-1)
png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_heatmap_Context_chunk_NOTR.png",    # create PNG for the heat map        
    width = 30*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11)
heatmap.2((aaxlsi[,aa_co_sort]), 
          trace='none', Colv = F,
          dendrogram = "none",Rowv = F,
          ColSideColors = aa_col_color,
          col=rev(cols), 
          breaks = aabreak,
          margins = c(12,8),
          density.info="none",
          key.title = "",
          key.xlab = "",
          key.ylab = "",
          keysize = 0.75)
dev.off()



#################3


aaxlsi <- lapply(lncRNA_feature_seq_rand_summary_list  , "[[", "median")
aaxlsi <- do.call(rbind, aaxlsi)
colnames(aaxlsi) <- lncRNA_feature_seq_rand_summary_list[[1]]$feature
rownames(aaxlsi) <- rownames(aa_all_auprc_chunk_cv2TX_noDist)

aaxlsi <- aaxlsi[,!(colSums(is.na(aaxlsi)) == nrow(aaxlsi))]
aamin <- apply(aaxlsi, 2, min, na.rm=T)
aaquant <- quantile(aamin, prob = seq(0,1,0.01), na.rm=T)
aaxlsi <- aaxlsi[, aamin <= aaquant[10]]
#aaarr <- range(aaxlsi, na.rm = T)
#aabreak <- c(0,5,10,20,50,75,100,155,200,250,aaarr[2])
aapostname <- unlist(lapply(strsplit(colnames(aaxlsi), "__"), function(x) x[[length(x)]]))
aapostname2 <- strsplit(aapostname, "_")
aapostname2_2 <- unlist(lapply(aapostname2, function(x) ifelse((x[[1]] != "triplex"), yes = x[[min(2, length(x))]], no = "triplex") ))
table(aapostname2_2)
aa_co_sort <- sort(aapostname2_2, index.return=T)$ix
aa_col_color <- brewer.pal(9, "Pastel1")[factor(aapostname2_2[aa_co_sort])]
#aapostname2_2[aa_co_sort]

aabreak <- unique(quantile(aaxlsi, prob = seq(0,1,0.02),na.rm = T))
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(length(aabreak)-1)

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_heatmap_Seq_rand.png",    # create PNG for the heat map        
    width = 12*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11)
heatmap.2((aaxlsi[,aa_co_sort]), 
          trace='none', Colv = T,Rowv = T,
          dendrogram = "none",
          ColSideColors = aa_col_color,
          col=rev(cols), 
          breaks = aabreak,
          margins = c(12,8),      
          density.info="none",
          key.title = "",
          key.xlab = "",
          key.ylab = "",
          keysize = 0.75)
dev.off()

#################


aaxlsi <- lapply(lncRNA_feature_seq_chunk_summary_list  , "[[", "median")
aaxlsi <- do.call(rbind, aaxlsi)
colnames(aaxlsi) <- lncRNA_feature_seq_chunk_summary_list[[1]]$feature
rownames(aaxlsi) <- rownames(aa_all_auprc_chunk_cv2TX_noDist)


aaxlsi <- aaxlsi[,!(colSums(is.na(aaxlsi)) == nrow(aaxlsi))]
aamin <- apply(aaxlsi, 2, min, na.rm=T)
aaquant <- quantile(aamin, prob = seq(0,1,0.01), na.rm=T)
aaxlsi <- aaxlsi[, aamin <= aaquant[10]]
#aaarr <- range(aaxlsi, na.rm = T)
#aabreak <- c(0,5,10,20,50,75,100,155,200,250,aaarr[2])
aapostname <- unlist(lapply(strsplit(colnames(aaxlsi), "__"), function(x) x[[length(x)]]))
aapostname2 <- strsplit(aapostname, "_")
aapostname2_2 <- unlist(lapply(aapostname2, function(x) ifelse((x[[1]] != "triplex"), yes = x[[min(2, length(x))]], no = "triplex") ))
table(aapostname2_2)
aa_co_sort <- sort(aapostname2_2, index.return=T)$ix
aa_col_color <- brewer.pal(9, "Pastel1")[factor(aapostname2_2[aa_co_sort])]
#aapostname2_2[aa_co_sort]

aabreak <- unique(quantile(aaxlsi, prob = seq(0,1,0.02),na.rm = T))
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(length(aabreak)-1)

png(filename = "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_heatmap_Seq_chunk.png",    # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 12)
heatmap.2((aaxlsi[,aa_co_sort]), 
          trace='none', Colv = F,
          dendrogram = "none",
          ColSideColors = aa_col_color,
          col=rev(cols), 
          breaks = aabreak,
          margins = c(12,8),
          density.info="none",
          key.title = "",
          key.xlab = "",
          key.ylab = "",
          keysize = 0.75)
dev.off()


#######################################################################################################################
#######################################################################################################################
# for each lncRNA 
# do two 100 boxplots of the top 100 sequence or context feaures: full dataset, positive against negative
# also note down the positive tiles with values in top 10 percentile for the feature

aalnc_feature_top_pos_tile_list <- list()
aaln <- unique(aatmp_lnc_df_cv2TX_NoD$lncRNA)

lncRNA_feature_seq_chunk_summary_list
lncRNA_feature_seq_rand_summary_list
lncRNA_feature_context_chunk_summary_list
lncRNA_feature_context_rand_summary_list


aafile <- list.files("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/")
for(i in 1:length(aafile)){
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/", aafile[i]))
  aa_myf1 <- lncRNA_feature_seq_chunk_summary_list[[i]][lncRNA_feature_seq_chunk_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("TR_sequ", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_boxplot/fb_",aaln[i],"_Seq_chunk.png"),    # create PNG for the heat map        
      width = 16*300,        # 5 x 300 pixels
      height = 16*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(2,2,2,2))
  for(j in 1:length(aa_myf11)){
    boxplot(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])]~my_Dataset$label, main = aa_myf11[j])
  }
  
  aa_myf1 <- lncRNA_feature_seq_rand_summary_list[[i]][lncRNA_feature_seq_rand_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("TR_sequ", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  dev.off()
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_boxplot/fb_",aaln[i],"_Seq_rand.png"),    # create PNG for the heat map        
      width = 16*300,        # 5 x 300 pixels
      height = 16*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(2,2,2,2))
  for(j in 1:length(aa_myf11)){
    boxplot(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])]~my_Dataset$label, main = aa_myf11[j])
  }
  dev.off()
  aa_myf1 <- lncRNA_feature_context_rand_summary_list[[i]][lncRNA_feature_context_rand_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("TR_sequ", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_boxplot/fb_",aaln[i],"_Context_rand.png"),    # create PNG for the heat map        
      width = 16*300,        # 5 x 300 pixels
      height = 16*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(2,2,2,2))
  for(j in 1:length(aa_myf11)){
    boxplot(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])]~my_Dataset$label, main = aa_myf11[j])
  }
  dev.off()
  aa_myf1 <- lncRNA_feature_context_chunk_summary_list[[i]][lncRNA_feature_context_chunk_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("TR_sequ", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_boxplot/fb_",aaln[i],"_Context_chunk.png"),    # create PNG for the heat map        
      width = 16*300,        # 5 x 300 pixels
      height = 16*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(2,2,2,2))
  for(j in 1:length(aa_myf11)){
    boxplot(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])]~my_Dataset$label, main = aa_myf11[j])
  }
  dev.off()
  
  
}

######## plot overlapping histograms instead of the above boxplots

aafile <- list.files("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/")

c1 <- rgb(0,0,255,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,0,0, max = 255, alpha = 80, names = "lt.red")
length(aafile)
for(i in 1:1){
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/", aafile[i]))
  aa_myf1 <- lncRNA_feature_seq_chunk_summary_list[[i]][lncRNA_feature_seq_chunk_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("TR_sequ", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_Hist_ovl/fh_",aaln[i],"_Seq_chunk.png"),    # create PNG for the heat map        
      width = 16*300,        # 5 x 300 pixels
      height = 16*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(2,2,2,2))
  for(j in 1:length(aa_myf11)){
    aa_myneg <- my_Dataset[my_Dataset$label == "neg", match(aa_myf11[j], my_name_dic[,1])]
    aa_mypos <- my_Dataset[my_Dataset$label == "pos", match(aa_myf11[j], my_name_dic[,1])]
    aa_f_range <- range(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])], na.rm = T)
    ax <- pretty(floor(aa_f_range[1]):ceiling(aa_f_range[2]), n = 20)
    histneg <- hist(aa_myneg,  breaks = ax , plot = F)
    histpos <- hist(aa_mypos,  breaks = ax , plot = F)
    plot(histneg, col = c1, freq = F, main = aa_myf11[j]) # Plot 1st histogram using a transparent color
    plot(histpos, col = c2, add = TRUE, freq = F) # Add 2nd histogram using different color
    #boxplot(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])]~my_Dataset$label, main = aa_myf11[j])
  }
  
  aa_myf1 <- lncRNA_feature_seq_rand_summary_list[[i]][lncRNA_feature_seq_rand_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("TR_sequ", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  dev.off()
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_Hist_ovl/fh_",aaln[i],"_Seq_rand.png"),    # create PNG for the heat map        
      width = 16*300,        # 5 x 300 pixels
      height = 16*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(2,2,2,2))
  for(j in 1:length(aa_myf11)){
    aa_myneg <- my_Dataset[my_Dataset$label == "neg", match(aa_myf11[j], my_name_dic[,1])]
    aa_mypos <- my_Dataset[my_Dataset$label == "pos", match(aa_myf11[j], my_name_dic[,1])]
    aa_f_range <- range(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])], na.rm = T)
    ax <- pretty(floor(aa_f_range[1]):ceiling(aa_f_range[2]), n = 20)
    histneg <- hist(aa_myneg,  breaks = ax , plot = F)
    histpos <- hist(aa_mypos,  breaks = ax , plot = F)
    plot(histneg, col = c1, freq = F, main = aa_myf11[j]) # Plot 1st histogram using a transparent color
    plot(histpos, col = c2, add = TRUE, freq = F) # Add 2nd histogram using different color
    
  }
  dev.off()
  aa_myf1 <- lncRNA_feature_context_rand_summary_list[[i]][lncRNA_feature_context_rand_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("TR_sequ", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_Hist_ovl/fh_",aaln[i],"_Context_rand.png"),    # create PNG for the heat map        
      width = 16*300,        # 5 x 300 pixels
      height = 16*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(2,2,2,2))
  for(j in 1:length(aa_myf11)){
    aa_myneg <- my_Dataset[my_Dataset$label == "neg", match(aa_myf11[j], my_name_dic[,1])]
    aa_mypos <- my_Dataset[my_Dataset$label == "pos", match(aa_myf11[j], my_name_dic[,1])]
    aa_f_range <- range(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])], na.rm = T)
    ax <- pretty(floor(aa_f_range[1]):ceiling(aa_f_range[2]), n = 20)
    histneg <- hist(aa_myneg,  breaks = ax , plot = F)
    histpos <- hist(aa_mypos,  breaks = ax , plot = F)
    plot(histneg, col = c1, freq = F, main = aa_myf11[j]) # Plot 1st histogram using a transparent color
    plot(histpos, col = c2, add = TRUE, freq = F) # Add 2nd histogram using different color
    
  }
  dev.off()
  aa_myf1 <- lncRNA_feature_context_chunk_summary_list[[i]][lncRNA_feature_context_chunk_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("TR_sequ", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_Hist_ovl/fh_",aaln[i],"_Context_chunk.png"),    # create PNG for the heat map        
      width = 16*300,        # 5 x 300 pixels
      height = 16*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(2,2,2,2))
  for(j in 1:length(aa_myf11)){
    aa_myneg <- my_Dataset[my_Dataset$label == "neg", match(aa_myf11[j], my_name_dic[,1])]
    aa_mypos <- my_Dataset[my_Dataset$label == "pos", match(aa_myf11[j], my_name_dic[,1])]
    aa_f_range <- range(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])], na.rm = T)
    ax <- pretty(floor(aa_f_range[1]):ceiling(aa_f_range[2]), n = 20)
    histneg <- hist(aa_myneg,  breaks = ax , plot = F)
    histpos <- hist(aa_mypos,  breaks = ax , plot = F)
    plot(histneg, col = c1, freq = F, main = aa_myf11[j]) # Plot 1st histogram using a transparent color
    plot(histpos, col = c2, add = TRUE, freq = F) # Add 2nd histogram using different color
    
  }
  dev.off()
  
  
}
######## plot qqplots instead of the above overlapping histograms

aafile <- list.files("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/")
aaln <- unique(aatmp_lnc_df_cv2TX_NoD$lncRNA)

c1 <- rgb(0,0,255,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,0,0, max = 255, alpha = 80, names = "lt.red")

load("Partition_6_CVfirst_dataset/Partition_6_CVfirst_dataset_2410003L11Rik.RData")
aa_colnn <- my_name_dic[1:(nrow(my_name_dic) - 1),1]  #c(colnames(Partition_6_feature_mat_tile_2), colnames(Partition_6_feature_mat_pair_2))

#aa_filter_list3 <- list()
aass1 <- grep(pattern = "__RBPscan", x = aa_colnn)
aass2 <- grep(pattern = "RBP_pair", x = aa_colnn)
#aa_filter_list3[[1]] <- union(aass1, aass2)
aass1r <- grep(pattern = "__rep", x = aa_colnn)
aass2r <- grep(pattern = "repeat_pair", x = aa_colnn)
#aa_filter_list3[[2]] <- union(aass1r, aass2r)

aaaalfil <- aa_colnn[unique(c(aass1, aass2, aass1r, aass2r))]
for(i in 1:length(aafile)){
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/", aafile[i]))
  aa_myf1 <- lncRNA_feature_seq_chunk_summary_list[[i]][lncRNA_feature_seq_chunk_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("TR_sequ", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }

  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]

  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_qq/fq_",aaln[i],"_Seq_chunk.png"),    # create PNG for the heat map
      width = 18*300,        # 5 x 300 pixels
      height = 18*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(4,4,3,2))
  for(j in 1:length(aa_myf11)){
    aa_myneg <- my_Dataset[my_Dataset$label == "neg", match(aa_myf11[j], my_name_dic[,1])]
    aa_mypos <- my_Dataset[my_Dataset$label == "pos", match(aa_myf11[j], my_name_dic[,1])]
    # aa_f_range <- range(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])], na.rm = T)
    # ax <- pretty(floor(aa_f_range[1]):ceiling(aa_f_range[2]), n = 20)
    aar <- range(c(aa_myneg,aa_mypos))

    qqplot(x = aa_myneg, y = aa_mypos,plot.it = T,xlab = "-",ylab = "+", main = aa_myf11[j],pch=19,cex = 0.7, xlim = aar, ylim = aar)
    abline(a = 0,b = 1,col=2)
    # histneg <- hist(aa_myneg,  breaks = ax , plot = F)
    # histpos <- hist(aa_mypos,  breaks = ax , plot = F)
    # plot(histneg, col = c1, freq = F, main = aa_myf11[j]) # Plot 1st histogram using a transparent color
    # plot(histpos, col = c2, add = TRUE, freq = F) # Add 2nd histogram using different color
    #boxplot(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])]~my_Dataset$label, main = aa_myf11[j])
  }

  aa_myf1 <- lncRNA_feature_seq_rand_summary_list[[i]][lncRNA_feature_seq_rand_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("TR_sequ", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  dev.off()

  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_qq/fq_",aaln[i],"_Seq_rand.png"),    # create PNG for the heat map
      width = 18*300,        # 5 x 300 pixels
      height = 18*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(4,4,3,2))
  for(j in 1:length(aa_myf11)){
    aa_myneg <- my_Dataset[my_Dataset$label == "neg", match(aa_myf11[j], my_name_dic[,1])]
    aa_mypos <- my_Dataset[my_Dataset$label == "pos", match(aa_myf11[j], my_name_dic[,1])]
    aar <- range(c(aa_myneg,aa_mypos))

    qqplot(x = aa_myneg, y = aa_mypos,plot.it = T,xlab = "-",ylab = "+", main = aa_myf11[j],pch=19,cex = 0.7, xlim = aar, ylim = aar)
    abline(a = 0,b = 1,col=2)
    # aa_f_range <- range(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])], na.rm = T)
    # ax <- pretty(floor(aa_f_range[1]):ceiling(aa_f_range[2]), n = 20)
    # histneg <- hist(aa_myneg,  breaks = ax , plot = F)
    # histpos <- hist(aa_mypos,  breaks = ax , plot = F)
    # plot(histneg, col = c1, freq = F, main = aa_myf11[j]) # Plot 1st histogram using a transparent color
    # plot(histpos, col = c2, add = TRUE, freq = F) # Add 2nd histogram using different color

  }
  dev.off()
  
  #############################################
  ############################################# adding  filtered features
  
  
  aa_myf1 <- lncRNA_feature_seq_chunk_summary_list[[i]][lncRNA_feature_seq_chunk_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("TR_sequ", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_qq/fq_",aaln[i],"_Seq_filtered_chunk.png"),    # create PNG for the heat map        
      width = 18*300,        # 5 x 300 pixels
      height = 18*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(4,4,3,2))
  for(j in 1:length(aa_myf11)){
    aamyma1 <- match(rownames(my_Dataset), partition_6_expression_features$tile_name)
    aamyma2 <- match(c("GROseq", "PROseq"), colnames(partition_6_expression_features))
    aa_cur_exp <- partition_6_expression_features[aamyma1, aamyma2]
    scale_vec <- log2(apply(aa_cur_exp[, c("GROseq", "PROseq")], 1, max) + 1)
    if(aa_myf11[j] %in% aaaalfil){
      aacursss <- my_Dataset[, match(aa_myf11[j], my_name_dic[,1])]
      aacursss <-  aacursss * scale_vec
    }else{
      aacursss <- my_Dataset[, match(aa_myf11[j], my_name_dic[,1])]
    }
    aa_myneg <- aacursss[my_Dataset$label == "neg"]
    aa_mypos <- aacursss[my_Dataset$label == "pos"]
    

    # aa_f_range <- range(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])], na.rm = T)
    # ax <- pretty(floor(aa_f_range[1]):ceiling(aa_f_range[2]), n = 20)
    aar <- range(c(aa_myneg,aa_mypos))
    
    qqplot(x = aa_myneg, y = aa_mypos,plot.it = T,xlab = "-",ylab = "+", main = aa_myf11[j],pch=19,cex = 0.7, xlim = aar, ylim = aar)
    abline(a = 0,b = 1,col=2)
    # histneg <- hist(aa_myneg,  breaks = ax , plot = F)
    # histpos <- hist(aa_mypos,  breaks = ax , plot = F)
    # plot(histneg, col = c1, freq = F, main = aa_myf11[j]) # Plot 1st histogram using a transparent color
    # plot(histpos, col = c2, add = TRUE, freq = F) # Add 2nd histogram using different color
    #boxplot(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])]~my_Dataset$label, main = aa_myf11[j])
  }
  
  aa_myf1 <- lncRNA_feature_seq_rand_summary_list[[i]][lncRNA_feature_seq_rand_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("TR_sequ", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  dev.off()
  
  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_qq/fq_",aaln[i],"_Seq_filtered_rand.png"),    # create PNG for the heat map        
      width = 18*300,        # 5 x 300 pixels
      height = 18*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(4,4,3,2))
  for(j in 1:length(aa_myf11)){
    aamyma1 <- match(rownames(my_Dataset), partition_6_expression_features$tile_name)
    aamyma2 <- match(c("GROseq", "PROseq"), colnames(partition_6_expression_features))
    aa_cur_exp <- partition_6_expression_features[aamyma1, aamyma2]
    scale_vec <- log2(apply(aa_cur_exp[, c("GROseq", "PROseq")], 1, max) + 1)
    if(aa_myf11[j] %in% aaaalfil){
      aacursss <- my_Dataset[, match(aa_myf11[j], my_name_dic[,1])]
      aacursss <-  aacursss * scale_vec
    }else{
      aacursss <- my_Dataset[, match(aa_myf11[j], my_name_dic[,1])]
    }
    aa_myneg <- aacursss[my_Dataset$label == "neg"]
    aa_mypos <- aacursss[my_Dataset$label == "pos"]
    aar <- range(c(aa_myneg,aa_mypos))
    
    qqplot(x = aa_myneg, y = aa_mypos,plot.it = T,xlab = "-",ylab = "+", main = aa_myf11[j],pch=19,cex = 0.7, xlim = aar, ylim = aar)
    abline(a = 0,b = 1,col=2)
    # aa_f_range <- range(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])], na.rm = T)
    # ax <- pretty(floor(aa_f_range[1]):ceiling(aa_f_range[2]), n = 20)
    # histneg <- hist(aa_myneg,  breaks = ax , plot = F)
    # histpos <- hist(aa_mypos,  breaks = ax , plot = F)
    # plot(histneg, col = c1, freq = F, main = aa_myf11[j]) # Plot 1st histogram using a transparent color
    # plot(histpos, col = c2, add = TRUE, freq = F) # Add 2nd histogram using different color
    
  }
  dev.off()
  ######################################################################
  
  aa_myf1 <- lncRNA_feature_context_rand_summary_list[[i]][lncRNA_feature_context_rand_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("TR_sequ", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]

  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_qq/fq_",aaln[i],"_Context_rand.png"),    # create PNG for the heat map
      width = 18*300,        # 5 x 300 pixels
      height = 18*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(4,4,3,2))
  for(j in 1:length(aa_myf11)){
    aa_myneg <- my_Dataset[my_Dataset$label == "neg", match(aa_myf11[j], my_name_dic[,1])]
    aa_mypos <- my_Dataset[my_Dataset$label == "pos", match(aa_myf11[j], my_name_dic[,1])]
    if(length(grep("TRNSCRP", aa_myf11[j])) > 0){
      aa_myneg <- log2(aa_myneg + 1)
      aa_mypos <- log2(aa_mypos + 1)
    }
    aar <- range(c(aa_myneg,aa_mypos))
    qqplot(x = aa_myneg, y = aa_mypos,plot.it = T,xlab = "-",ylab = "+", main = aa_myf11[j],pch=19,cex = 0.7, xlim = aar, ylim = aar)
    abline(a = 0,b = 1,col=2)
    # aa_f_range <- range(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])], na.rm = T)
    # ax <- pretty(floor(aa_f_range[1]):ceiling(aa_f_range[2]), n = 20)
    # histneg <- hist(aa_myneg,  breaks = ax , plot = F)
    # histpos <- hist(aa_mypos,  breaks = ax , plot = F)
    # plot(histneg, col = c1, freq = F, main = aa_myf11[j]) # Plot 1st histogram using a transparent color
    # plot(histpos, col = c2, add = TRUE, freq = F) # Add 2nd histogram using different color

  }
  dev.off()
  aa_myf1 <- lncRNA_feature_context_chunk_summary_list[[i]][lncRNA_feature_context_chunk_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("TR_sequ", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]

  png(filename = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Combined_res_noD_TX_plots/Feature_qq/fq_",aaln[i],"_Context_chunk.png"),    # create PNG for the heat map
      width = 18*300,        # 5 x 300 pixels
      height = 18*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(10,10), mar = c(4,4,3,2))
  for(j in 1:length(aa_myf11)){
    aa_myneg <- my_Dataset[my_Dataset$label == "neg", match(aa_myf11[j], my_name_dic[,1])]
    aa_mypos <- my_Dataset[my_Dataset$label == "pos", match(aa_myf11[j], my_name_dic[,1])]
    if(length(grep("TRNSCRP", aa_myf11[j])) > 0){
      aa_myneg <- log2(aa_myneg + 1)
      aa_mypos <- log2(aa_mypos + 1)
    }
    aar <- range(c(aa_myneg,aa_mypos))
    qqplot(x = aa_myneg, y = aa_mypos,plot.it = T,xlab = "-",ylab = "+", main = aa_myf11[j],pch=19,cex = 0.7, xlim = aar, ylim = aar)
    abline(a = 0,b = 1,col=2)
    # aa_f_range <- range(my_Dataset[, match(aa_myf11[j], my_name_dic[,1])], na.rm = T)
    # ax <- pretty(floor(aa_f_range[1]):ceiling(aa_f_range[2]), n = 20)
    # histneg <- hist(aa_myneg,  breaks = ax , plot = F)
    # histpos <- hist(aa_mypos,  breaks = ax , plot = F)
    # plot(histneg, col = c1, freq = F, main = aa_myf11[j]) # Plot 1st histogram using a transparent color
    # plot(histpos, col = c2, add = TRUE, freq = F) # Add 2nd histogram using different color

  }
  dev.off()
  
  
}
######## get feature index for running a model for each lncRNA

aafile <- list.files("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/")

aa_model_feat_list_rand <- list()
aa_model_feat_list_chunk <- list()
aaln <- unique(aatmp_lnc_df_cv2TX_NoD$lncRNA)
load(paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/", aafile[1]))

for(i in 1:length(aafile)){
  
  aa_myf1 <- lncRNA_feature_seq_chunk_summary_list[[i]][lncRNA_feature_seq_chunk_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("TR_sequ", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  aa_model_feat_list_chunk[[i]] <- aa_myf11 # which(my_name_dic[,1] %in% aa_myf11)
  
  
  
  aa_myf1 <- lncRNA_feature_seq_rand_summary_list[[i]][lncRNA_feature_seq_rand_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("TR_sequ", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  aa_model_feat_list_rand[[i]] <- aa_myf11 #which(my_name_dic[,1] %in% aa_myf11)
  
  
  aa_myf1 <- lncRNA_feature_context_rand_summary_list[[i]][lncRNA_feature_context_rand_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("TR_sequ", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  aa_model_feat_list_rand[[i]] <- c(aa_model_feat_list_rand[[i]], aa_myf11)# which(my_name_dic[,1] %in% ))
  
  
  aa_myf1 <- lncRNA_feature_context_chunk_summary_list[[i]][lncRNA_feature_context_chunk_summary_list[[i]]$N > 1,]
  aa_myf1 <- aa_myf1[!(is.na(aa_myf1$median)),]
  aaggr <- grep("TR_sequ", aa_myf1$feature)
  if(length(aaggr) > 0){
    aa_myf1 <- aa_myf1[-aaggr,]
  }
  aa_myf11 <- aa_myf1$feature[sort(aa_myf1$median, index.return=T)$ix[1:min(100, nrow(aa_myf1))]]
  aa_model_feat_list_chunk[[i]] <- c(aa_model_feat_list_chunk[[i]], aa_myf11) # which(my_name_dic[,1] %in% aa_myf11))
  my_200_feat <- aa_model_feat_list_rand[[i]]
  save(list = c("my_200_feat"), file = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Feature200_rerun1/perlncIndex/index_", aaln[i], "_R.RData"))
  my_200_feat <- aa_model_feat_list_chunk[[i]]
  save(list = c("my_200_feat"), file = paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Feature200_rerun1/perlncIndex/index_", aaln[i], "_C.RData"))
  
  
}
names(aa_model_feat_list_chunk) <- aaln
names(aa_model_feat_list_rand) <- aaln

unlist(lapply(aa_model_feat_list_chunk, length))




#############################3
# Reading results of top200_feature model

for(i in 1:length(aa_un_own)){
  cat(c("Rscript --vanilla /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/gather_perf.R",
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_200feat_rerun1/Learned_models/",aa_un_own[i] ),
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_200feat_rerun1/Combined_results/",aa_un_own[i], "__perfMat.RData\n" )),
      sep = " ", 
      file = "~/Documents/Shayan/BioInf/lncRNA/gather_perf_all_CV2_200feat_rerun1.job",
      append = !(i == 1))
}

aafiles <- list.files("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Feature200_rerun1/Combined_res")
aafn <- unlist(lapply(strsplit(aafiles, "__"), "[[", 1))
aa_all_auroc_list_cv2TX200feat <- list()
aa_all_auprc_list_cv2TX200feat <- list()
aa_all_impor_list_cv2TX200feat<- list()

for(i in 1:length(aafiles)){
  print(i)
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Feature200_rerun1/Combined_res/", aafiles[i]))
  aa_all_auroc_list_cv2TX200feat[[i]] <- my_perf[[1]][[1]]
  aa_all_auprc_list_cv2TX200feat[[i]] <- my_perf[[2]][[1]]
  aa_all_impor_list_cv2TX200feat[[i]] <- my_perf[[3]][[1]]
}
names(aa_all_auroc_list_cv2TX200feat) <- aafn
names(aa_all_auprc_list_cv2TX200feat) <- aafn
names(aa_all_impor_list_cv2TX200feat) <- aafn

View(aa_all_auroc_list_cv2TX200feat$Malat1)




aac <- max(unlist(lapply(aa_all_auprc_list_cv2TX200feat, nrow)))

aa_all_auroc_rand_cv2TX200feat <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auroc_rand_cv2TX200feat) <- rownames(aa_all_auprc_list_cv2TX200feat[[2]])
rownames(aa_all_auroc_rand_cv2TX200feat) <- aafn

aa_all_auroc_chunk_cv2TX200feat <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auroc_chunk_cv2TX200feat) <- rownames(aa_all_auprc_list_cv2TX200feat[[2]])
rownames(aa_all_auroc_chunk_cv2TX200feat) <- aafn

aa_all_auprc_rand_cv2TX200feat <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auprc_rand_cv2TX200feat) <- rownames(aa_all_auprc_list_cv2TX200feat[[2]])
rownames(aa_all_auprc_rand_cv2TX200feat) <- aafn

aa_all_auprc_chunk_cv2TX200feat <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auprc_chunk_cv2TX200feat) <- rownames(aa_all_auprc_list_cv2TX200feat[[2]])
rownames(aa_all_auprc_chunk_cv2TX) <- aafn


aa_all_auprc_list_cv2_diffTX200feat <- list()
for(i in 1:length(aa_all_auprc_list_cv2TX200feat)){
  aa_all_auprc_list_cv2_diffTX200feat[[i]] <- aa_all_auprc_list_cv2TX200feat[[i]]
  for(j in 1:nrow(aa_all_auprc_list_cv2TX200feat[[i]])){
    aa_all_auprc_list_cv2_diffTX200feat[[i]][j,] <- aa_all_auprc_list_cv2_diffTX200feat[[i]][j,] - aa_auprc_base[i,]
  }
}
for(i in 1:length(aa_all_auroc_list_cv2TX200feat)){
  aatmmpccv_roc <- rowMeans(aa_all_auroc_list_cv2TX200feat[[i]][,c(1:5)], na.rm = T)
  aatmmprcv_roc <- rowMeans(aa_all_auroc_list_cv2TX200feat[[i]][,c(6:10)], na.rm = T)
  
  aatmmpccv_prc <- rowMeans(aa_all_auprc_list_cv2_diffTX200feat[[i]][,c(1:5)], na.rm = T)
  aatmmprcv_prc <- rowMeans(aa_all_auprc_list_cv2_diffTX200feat[[i]][,c(6:10)], na.rm = T)
  
  aacolord <- match(rownames(aa_all_auroc_list_cv2TX200feat[[i]]), 
                    colnames(aa_all_auroc_chunk_cv2TX200feat))
  
  aa_all_auroc_chunk_cv2TX200feat[i,aacolord] <- aatmmpccv_roc
  aa_all_auroc_rand_cv2TX200feat[i,aacolord] <- aatmmprcv_roc
  aa_all_auprc_chunk_cv2TX200feat[i,aacolord] <- aatmmpccv_prc
  aa_all_auprc_rand_cv2TX200feat[i,aacolord] <- aatmmprcv_prc
}

View(aa_all_auroc_rand_cv2TX200feat)
View(aa_all_auroc_chunk_cv2TX200feat)

aa_all_auprc_rand_cv2TX200feat
aa_all_auprc_chunk_cv2TX200feat




# biomTrack <- BiomartGeneRegionTrack(genome = "mm9",
#                                     chromosome = chr, start = 20000000, end = 21000000,
#                                     name = "ENSEMBL")
# 

##############################
# write visulization jobs --> Gviz and feature importance

aa_filt_dic <- c("no_filter",
                 "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_CV2_filt/my_filter_RBP.txt",
                 "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_CV2_filt/my_filter_repeat.txt",
                 "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_CV2_filt/my_filter_both.txt")

for(i in 1:length(aaln)){
  cat(c("Rscript",
        "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/plot_gviz_200Feat.R",
        aaln[i],
        "R", 
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_200feat_rerun1/Gviz/", aaln[i]), 
        "20", 
        c("top200Features_RBPFilter"),#, "top200Features_BothFilter", "top200Features_RepFilter","top200Features_noFilter"), 
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_CVfirst_dataset_", aaln[i], ".RData"),
        aa_filt_dic[2],
        "\n"), file = "plot_gviz_top20_200Feat_rerun1.job", append = !(i==1))
  cat(c("Rscript",
        "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/plot_gviz_200Feat.R",
        aaln[i],
        "C", 
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_200feat_rerun1/Gviz/", aaln[i]), 
        "20", 
        c("top200Features_RBPFilter"),#, "top200Features_BothFilter", "top200Features_RepFilter","top200Features_noFilter"),
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_CVfirst_dataset_", aaln[i], ".RData"),
        aa_filt_dic[2],
        "\n"), file = "plot_gviz_top20_200Feat_rerun1.job", append = T)
  
}


#######################################################################################################################################
#######################################################################################################################################

# write jobs to visualize lnRNA preds
# Rscript plot_gviz_lncRNA.R Meg3 C /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_lowres1MG/Gviz dist kmerFq motSc_pairs_triplx sequ_dist sequ_dist_trnsp sequ_dist_Chrom_Meth_AX_ChIP_t

for(i in 1:length(aaln)){
  cat(c("Rscript",
        "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/plot_gviz_RandVsChunk.R",
        aaln[i],
        
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Random_vs_chunk_plots/", aaln[i]), 
 
        "\n"), file = "plot_gviz_all_randomVschunk.job", append = !(i==1))

  
}
#######################################################################################################################################
#######################################################################################################################################
# check the quantiles of Rbmx
aafile <- list.files("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/")
aaqua <- list()
for(i in 1:length(aafile)){
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/",aafile[i]))
  aaw <- grep("Rbmx", my_name_dic[, 1])
  aaqua[[i]] <- quantile(my_Dataset[, aaw], prob = seq(0,1,0.1))
}

#######################################################################################################################################
#######################################################################################################################################
# find genes overlapping with Gm14820 and Neat1 binding sites

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/Partition_6_CVfirst_dataset_Gm14820.RData")
aa_my_reg <- Partition_6_dfs_GR[match(rownames(my_Dataset)[my_Dataset$label == "pos"],  Partition_6_dfs_GR$tile )]


#######################################################################################################################################
# gather performance for pair runs

aa_un_own <- sort(unique(Partition_6_random_chunk_cv_df$owner))


for(i in 1:length(aa_un_own)){
  cat(c("Rscript --vanilla /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/gather_perf_wt.R",
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_pair/Learned_models/",aa_un_own[i] ),
        paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_pair/Combined_results/",aa_un_own[i], "__perfMat.RData\n" )),
      sep = " ", 
      file = "~/Documents/Shayan/BioInf/lncRNA/gather_perf_all_CV_pair.job",
      append = !(i == 1))
}



aafiles <- list.files("~/Documents/Shayan/BioInf/lncRNA/Family_svd/pair_run/")
aafn <- unlist(lapply(strsplit(aafiles, "__"), "[[", 1))
aa_all_auroc_list_cvpair <- list()
aa_all_auprc_list_cvpair <- list()
aa_all_impor_list_cvpair <- list()

for(i in 1:length(aafiles)){
  print(i)
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Family_svd/pair_run/", aafiles[i]))
  aa_all_auroc_list_cvpair[[i]] <- my_perf[[1]][[1]]
  aa_all_auprc_list_cvpair[[i]] <- my_perf[[2]][[1]]
  aa_all_impor_list_cvpair[[i]] <- my_perf[[3]][[1]]
}
names(aa_all_auroc_list_cvpair) <- aafn
names(aa_all_auprc_list_cvpair) <- aafn
names(aa_all_impor_list_cvpair) <- aafn

View(aa_all_auroc_list_cvpair$Malat1)




aac <- max(unlist(lapply(aa_all_auprc_list_cvpair, nrow)))

aa_all_auroc_rand_cv_pair <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auroc_rand_cv_pair) <- rownames(aa_all_auprc_list_cvpair[[2]])
rownames(aa_all_auroc_rand_cv_pair) <- aafn

aa_all_auroc_chunk_cv_pair <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auroc_chunk_cv_pair) <- rownames(aa_all_auprc_list_cvpair[[2]])
rownames(aa_all_auroc_chunk_cv_pair) <- aafn

aa_all_auprc_rand_cv_pair <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auprc_rand_cv_pair) <- rownames(aa_all_auprc_list_cvpair[[2]])
rownames(aa_all_auprc_rand_cv_pair) <- aafn

aa_all_auprc_chunk_cv_pair <- matrix(nrow = length(aafn), ncol = aac)
colnames(aa_all_auprc_chunk_cv_pair) <- rownames(aa_all_auprc_list_cvpair[[2]])
rownames(aa_all_auprc_chunk_cv_pair) <- aafn


aa_all_auprc_list_cv2_diff_pair <- list()
for(i in 1:length(aa_all_auprc_list_cvpair)){
  aa_all_auprc_list_cv2_diff_pair[[i]] <- aa_all_auprc_list_cvpair[[i]]
  for(j in 1:nrow(aa_all_auprc_list_cvpair[[i]])){
    aa_all_auprc_list_cv2_diff_pair[[i]][j,] <- aa_all_auprc_list_cv2_diff_pair[[i]][j,] - aa_auprc_base[i,]
  }
}
for(i in 1:length(aa_all_auroc_list_cvpair)){
  aatmmpccv_roc <- rowMeans(aa_all_auroc_list_cvpair[[i]][,c(1:5)], na.rm = T)
  aatmmprcv_roc <- rowMeans(aa_all_auroc_list_cvpair[[i]][,c(6:10)], na.rm = T)
  
  aatmmpccv_prc <- rowMeans(aa_all_auprc_list_cv2_diff_pair[[i]][,c(1:5)], na.rm = T)
  aatmmprcv_prc <- rowMeans(aa_all_auprc_list_cv2_diff_pair[[i]][,c(6:10)], na.rm = T)
  
  aacolord <- match(rownames(aa_all_auroc_list_cvpair[[i]]), 
                    colnames(aa_all_auroc_chunk_cv_pair))
  
  aa_all_auroc_chunk_cv_pair[i,aacolord] <- aatmmpccv_roc
  aa_all_auroc_rand_cv_pair[i,aacolord] <- aatmmprcv_roc
  aa_all_auprc_chunk_cv_pair[i,aacolord] <- aatmmpccv_prc
  aa_all_auprc_rand_cv_pair[i,aacolord] <- aatmmprcv_prc
}

View(aa_all_auroc_rand_cv_pair)
View(aa_all_auroc_chunk_cv_pair)

aa_all_auprc_rand_cv_pair
aa_all_auprc_chunk_cv_pair


# create a matrix to store performance of pairs: nrow= number of feature fams, ncol =  number of feature fams
aa_indiv <- c(colnames(aa_all_auroc_rand_cv_pair)[-grep("_a_", colnames(aa_all_auroc_rand_cv_pair))],
              "RBP_pair_a_TFscan_pair_a_PPI_pair_a_repeat_pair", "HISTONE_a_Chromatin_pair"   ,"__RBPscan_a___TFscan","ChIP_a_ChIP_pair")


aa_pair_auroc_list_r <- list()
aa_pair_auroc_list_c <- list()
aa_pair_auprc_list_r <- list()
aa_pair_auprc_list_c <- list()
aa_ind_auroc_list_r <- list()
aa_ind_auprc_list_r <- list()
aa_ind_auroc_list_c <- list()
aa_ind_auprc_list_c <- list()

for(i in 1:5){
  aa_pair_auroc_list_r[[i]] <- matrix(nrow = length(all_families),
                                    ncol = length(all_families))
  rownames(aa_pair_auroc_list_r[[i]]) <- aa_indiv
  colnames(aa_pair_auroc_list_r[[i]]) <- aa_indiv
  
  aa_pair_auroc_list_c[[i]] <- matrix(nrow = length(all_families),
                                      ncol = length(all_families))
  rownames(aa_pair_auroc_list_c[[i]]) <- aa_indiv
  colnames(aa_pair_auroc_list_c[[i]]) <- aa_indiv
  
  aa_pair_auprc_list_r[[i]] <- matrix(nrow = length(all_families),
                                      ncol = length(all_families))
  rownames(aa_pair_auprc_list_r[[i]]) <- aa_indiv
  colnames(aa_pair_auprc_list_r[[i]]) <- aa_indiv
  
  aa_pair_auprc_list_c[[i]] <- matrix(nrow = length(all_families),
                                      ncol = length(all_families))
  rownames(aa_pair_auprc_list_c[[i]]) <- aa_indiv
  colnames(aa_pair_auprc_list_c[[i]]) <- aa_indiv
  
  aa_ind_auroc_list_r[[i]] <- matrix(nrow = length(all_families),
                                     ncol = 1)
  rownames(aa_ind_auroc_list_r[[i]]) <- aa_indiv
  
  aa_ind_auroc_list_c[[i]] <- matrix(nrow = length(all_families),
                                     ncol = 1)
  rownames(aa_ind_auroc_list_c[[i]]) <- aa_indiv
  aa_ind_auprc_list_r[[i]] <- matrix(nrow = length(all_families),
                                     ncol = 1)
  rownames(aa_ind_auprc_list_r[[i]]) <- aa_indiv
  aa_ind_auprc_list_c[[i]] <- matrix(nrow = length(all_families),
                                     ncol = 1)
  rownames(aa_ind_auprc_list_c[[i]]) <- aa_indiv
}


names(aa_pair_auroc_list_r) <- paste0("RCV", c(1:5))
names(aa_pair_auprc_list_r) <- paste0("RCV", c(1:5))
names(aa_pair_auroc_list_c) <- paste0("CCV", c(1:5))
names(aa_pair_auprc_list_c) <- paste0("CCV", c(1:5))

names(aa_ind_auroc_list_r) <- paste0("RCV", c(1:5))
names(aa_ind_auprc_list_r) <- paste0("RCV", c(1:5))
names(aa_ind_auroc_list_c) <- paste0("CCV", c(1:5))
names(aa_ind_auprc_list_c) <- paste0("CCV", c(1:5))

for(i in 1:length(aa_indiv)){
  aa_1col <- match(aa_indiv[i], colnames(aa_all_auroc_chunk_cv_pair))
  for(aacurcv in 1:5){
    
    
    for(j in 1:length(aa_indiv)){
      
    }
  }
}



#colnames(aa_all_auroc_rand_cv_pair)[-grep("_a_", colnames(aa_all_auroc_rand_cv_pair))]

#aa_all_auroc_chunk_cv_pair
aa_comp_red_list <- list()
aa_indiv <- c(colnames(aa_all_auroc_rand_cv_pair)[-grep("_a_", colnames(aa_all_auroc_rand_cv_pair))],
              "RBP_pair_a_TFscan_pair_a_PPI_pair_a_repeat_pair", "HISTONE_a_Chromatin_pair"   ,"__RBPscan_a___TFscan","ChIP_a_ChIP_pair")
for(k in 1:nrow(aa_all_auroc_chunk_cv_pair)){
  aa_comp_red_list[[k]] <- matrix(nrow=length(aa_indiv), ncol = length(aa_indiv))
  rownames(aa_comp_red_list[[k]]) <- aa_indiv
  colnames(aa_comp_red_list[[k]]) <- aa_indiv
  
  for( i in 1:(length(aa_indiv)-1)){
    aa_1col <- match(aa_indiv[i], colnames(aa_all_auroc_chunk_cv_pair))
    stopifnot(length(aa_1col) == 1)
    for(j in (i+1):length(aa_indiv)){
      aa_2col <- match(aa_indiv[j], colnames(aa_all_auroc_chunk_cv_pair))
      stopifnot(length(aa_2col) == 1)
      aa_both_col <- intersect(grep(aa_indiv[i],colnames(aa_all_auroc_chunk_cv_pair)),
                               grep(aa_indiv[j],colnames(aa_all_auroc_chunk_cv_pair)))
      stopifnot(length(aa_both_col) == 1)
      aa_comp_red_list[[k]][i, j] <- aa_all_auroc_chunk_cv_pair[k,aa_both_col] - max(aa_all_auroc_chunk_cv_pair[k,aa_2col],
                                                                                     aa_all_auroc_chunk_cv_pair[k,aa_1col])
    }
  }
  
}
names(aa_comp_red_list) <- rownames(aa_all_auroc_chunk_cv_pair)

View(aa_comp_red_list$Neat1)

aa_name_mat <- aa_comp_red_list$Neat1
for(i in 1:ncol(aa_name_mat)){
  for(j in 1:nrow(aa_name_mat)){
    aa_name_mat[j, i] <- paste(rownames(aa_name_mat)[j], "__Vs__", colnames(aa_name_mat)[i])
  }
}

aa_full_com <- do.call(rbind, lapply(aa_comp_red_list, function(x) x[upper.tri(x)]))
rownames(aa_full_com) <- names(aa_comp_red_list)
colnames(aa_full_com) <- aa_name_mat[upper.tri(aa_name_mat)]
library(gplots)

aabreak <- c(seq(-0.1,0.1,0.001))
# "Blues"   "BuGn"    "BuPu"    "GnBu"    "Greens"  "Greys"   "Oranges" "OrRd"    "PuBu"    "PuBuGn"  "PuRd"    "Purples" "RdPu"    "Reds"    "YlGn"   
#  "YlGnBu"  "YlOrBr"  "YlOrRd" 
cols <- colorRampPalette(brewer.pal(9, "RdBu"))(length(aabreak)-1)


aaxxvc <- colnames(aa_full_com)

aa_indiv <- c(colnames(aa_all_auroc_rand_cv_pair)[-grep("_a_", colnames(aa_all_auroc_rand_cv_pair))],
              "RBP_pair_a_TFscan_pair_a_PPI_pair_a_repeat_pair", "HISTONE_a_Chromatin_pair"   ,"__RBPscan_a___TFscan","ChIP_a_ChIP_pair")
aa_indiv2 <- c("kmer", "Repeat", "AXSB", "Methylation","Triplex","TRNSCRP","PairSeq","Histone","MotifScan","ChIP")
for(i in 1:length(aa_indiv)){
  aaxxvc <- gsub(pattern = aa_indiv[i],replacement = aa_indiv2[i],x = aaxxvc)
}
aaxxvc <- gsub(pattern = "__",replacement = "",x = aaxxvc)
colnames(aa_full_com) <- aaxxvc
heatmap.2(aa_full_com[-match(aa_TOBEREMOVED, rownames(aa_full_com)),], 
          trace='none', Colv = T,Rowv = T,
          dendrogram = "both",
          #ColSideColors = aa_col_color,
          col=rev(cols), 
          breaks = aabreak,
          margins = c(12,8),   
          density.info="none",
          key.title = "",
          key.xlab = "",
          key.ylab = "",
          keysize = 1)


sort(aa_all_impor_list_cvpair$Trerf1$CCV1$HISTONE_a_Chromatin_pair_a_Methylation__WGBS, decreasing = T)[1:10]
sort(aa_all_impor_list_cvpair$Trerf1$CCV2$HISTONE_a_Chromatin_pair_a_Methylation__WGBS, decreasing = T)[1:10]
sort(aa_all_impor_list_cvpair$Trerf1$CCV3$HISTONE_a_Chromatin_pair_a_Methylation__WGBS, decreasing = T)[1:10]
sort(aa_all_impor_list_cvpair$Trerf1$CCV4$HISTONE_a_Chromatin_pair_a_Methylation__WGBS, decreasing = T)[1:10]

barplot(sort(aa_all_impor_list_cvpair$Trerf1$CCV3$Methylation__WGBS_a_ChIP_a_ChIP_pair, decreasing = T)[1:20],las = 2)
barplot(sort(aa_all_impor_list_cvpair$Neat1$CCV3$Methylation__WGBS_a_ChIP_a_ChIP_pair, decreasing = T)[1:20],las = 2)
barplot(sort(aa_all_impor_list_cvpair$Gm14820$CCV3$Methylation__WGBS_a_ChIP_a_ChIP_pair, decreasing = T)[1:20],las = 2)



