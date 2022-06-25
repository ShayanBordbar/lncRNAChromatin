# RF evaluation functions
require(Gviz)
require(IRanges)
############################################################################################################################################
# write a function that gets model fit, test data,  creates 
# two plots, auprc and auroc for training and test samples

perf_eval_ROC_PRC <- function(model_fit,
                              my_dataset,
                              my_label, 
                              file_name,
                              train_perf = F,
                              neg_value = "neg",
                              pos_value = "pos"){
  stopifnot(length(my_label) == nrow(my_dataset))
  require(ROCR)
  require(PRROC)
  if(train_perf){
    my_pred <- list()
    my_pred$predictions <- model_fit$predictions
  }else{
    my_pred <- predict(model_fit, data=my_dataset)
  }
  #my_pred <- predict(model_fit, data=my_dataset)
  my_prediction <- prediction(my_pred$predictions[, 2],my_label)
  prc_rec <- performance(my_prediction, measure="prec", x.measure="rec")
  tpr_fpr <- performance(my_prediction, measure = "tpr", x.measure = "fpr")
  auc <- performance(my_prediction, measure = "auc")
  auc <- auc@y.values[[1]]
  print("after auc")
  index_Pos <- my_label == pos_value
  index_Neg <- my_label == neg_value
  aaprc <- pr.curve(my_pred$predictions[index_Pos, 2], my_pred$predictions[index_Neg, 2], curve = T)
  #plot(aaprc)
  print("after prc")
  png(filename = file_name,    # create PNG for the heat map        
      width = 8*300,        # 5 x 300 pixels
      height = 4*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11) 
  par(mfrow = c(1,2), mar = c(4,4,4,4))
  plot(tpr_fpr, ylab="TPR", xlab = "FPR",
       main = paste0("auroc: ", (format(auc, digits = 2))),
       xlim = c(0,1), ylim = c(0,1))
  plot(aaprc$curve[,1], aaprc$curve[,2], type = "l", ylab="Precision", xlab = "Recall",
       main = paste0("auprc: ", (format(aaprc$auc.integral, digits = 2))),
       xlim = c(0,1), 
       ylim = c(0,1))
  #plot(aaprc)
  dev.off()
  return(list(Prec_Rec = prc_rec,
              TPR_FPR = tpr_fpr,
              AUROC = auc, 
              AUPRC = aaprc$auc.integral))
}
############################################################################################################################################
############################################################################################################################################

perf_eval_ROC_PRC_multimodel <- function(model_fit_list,
                              my_dataset_list,
                              my_label, 
                              file_name,
                              train_perf = F,
                              neg_value = "neg",
                              pos_value = "pos"){
  stopifnot(#length(my_label) == nrow(my_dataset), 
            length(names(model_fit_list)) == length(model_fit_list))
  require(ROCR)
  require(PRROC)
  my_pred_list <- list()
  index_Pos <- my_label == pos_value
  index_Neg <- my_label == neg_value
  
  perf_list <- list()
  for(i in 1:length(model_fit_list)){
    if(train_perf){
      my_pred_list[[i]] <- list()
      my_pred_list[[i]]$predictions <- model_fit_list[[i]]$predictions
      aana <- sum(is.na(my_pred_list[[i]]$predictions))
      if(aana > 0){
        print(paste0("There are ", aana, " NA predicted values in the trained model ", names(model_fit_list)[i], " . Will be set to zero..."))
        my_pred_list[[i]]$predictions[is.na(my_pred_list[[i]]$predictions)]<- 0
      }
    }else{
      my_pred_list[[i]] <- predict(model_fit_list[[i]], data=my_dataset_list[[i]])
    }
    perf_list[[i]] <- list()
    my_prediction <- prediction(my_pred_list[[i]]$predictions[, 2],my_label)
    perf_list[[i]]$prc_rec <- performance(my_prediction, measure="prec", x.measure="rec")
    perf_list[[i]]$tpr_fpr <- performance(my_prediction, measure = "tpr", x.measure = "fpr")
    auc <- performance(my_prediction, measure = "auc")
    perf_list[[i]]$auc <- auc@y.values[[1]]
    perf_list[[i]]$prc <- pr.curve(my_pred_list[[i]]$predictions[index_Pos, 2], my_pred_list[[i]]$predictions[index_Neg, 2], curve = T)
    
  }
  names(perf_list) <- names(model_fit_list)
  names(my_pred_list) <- names(model_fit_list)
  
  png(filename = file_name,    # create PNG for the heat map
      width = 16*300,        # 5 x 300 pixels
      height = 4*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11)
  par(mfrow = c(1,2), mar = c(4,5,4,15), xpd = T)
  for(i in 1:length(perf_list)){
    if(i == 1){
      plot(perf_list[[i]]$tpr_fpr, ylab="TPR", xlab = "FPR",
          # main = paste0("auroc: ", (format(perf_list[[i]]$auc, digits = 2))),
           xlim = c(0,1), ylim = c(0,1))
      }else{
        #print("i")
        #print(i)
        #print(perf_list[[i]]$tpr_fpr)
        lines(x = perf_list[[i]]$tpr_fpr@x.values[[1]],
              y = perf_list[[i]]$tpr_fpr@y.values[[1]], col = i)
      }
  }
  aa_all_auc<- unlist(lapply(perf_list, "[[", 3))
  legend(x = "topright", inset=c(-0.65,0),xpd = T,
         legend = paste(names(perf_list),
                                           format(aa_all_auc, digits = 2), 
                                           sep = "___"), 
         fill = c(1:i),
         bty = "n",
         cex = 0.7, x.intersp = 0.5, y.intersp = 0.6)
  for(i in 1:length(perf_list)){
    if(i == 1){
      plot(perf_list[[i]]$prc_rec, ylab="Precision", xlab = "Recall",
           #main = paste0("auprc: ", (format(aaprc$auc.integral, digits = 2))),
           xlim = c(0,1),
           ylim = c(0,1))
    }else{
      lines(x = perf_list[[i]]$prc_rec@x.values[[1]],
            y=perf_list[[i]]$prc_rec@y.values[[1]] , col = i)
    }
  }
  aa_all_prc<- lapply(perf_list, "[[", 4)
  aa_all_prc <- unlist(lapply(aa_all_prc, "[[", "auc.integral") )
  legend(x = "topright", inset=c(-0.65,0), xpd = T,
         legend = paste(names(perf_list), 
                        format(aa_all_prc, digits = 2), sep = "___"),
         
         bty = "n",
         fill = c(1:i),
         cex = 0.7, x.intersp = 0.5, y.intersp = 0.6)
  dev.off()
  return(perf_list)
}
############################################################################################################################################
# my_model_fit_list <- list()
# load("RF_partition2_random_caseWght.RData")
# load("Learned_models/Partition_2_modified_dataset_and_removedColumnNames___chunk___distance___RFmodel.RData")
# my_model_fit_list[[1]] <- RF_partition2_chunk_classWght
# my_model_fit_list[[2]] <- my_RF_model
# names(my_model_fit_list) <- c("random","random_noDis")
# aanampas <- paste0(names(my_model_fit_list), collapse = "__vs__")
# aa_pr_multimodel <- perf_eval_ROC_PRC_multimodel(model_fit_list = my_model_fit_list,
#                                          my_dataset = aa_test_data,
#                                          my_label = aa_test_data[, ncol(aa_test_data)],
#                                          file_name = paste0(plot_strore_address, aanampas, "___test_all.png"),
#                                          train_perf = F,
#                                          neg_value = "neg",
#                                          pos_value = "pos")
# aa_pr_multimodeltr <- perf_eval_ROC_PRC_multimodel(model_fit_list = my_model_fit_list,
#                                                  my_dataset = aa_train_data,
#                                                  my_label = aa_train_data[, ncol(aa_train_data)],
#                                                  file_name = paste0(plot_strore_address, aanampas, "___train_all.png"),
#                                                  train_perf = T,
#                                                  neg_value = "neg",
#                                                  pos_value = "pos")
# 

############################################################################################################################################
############################################################################################################################################
# write a function that gets model fit, test data, and the owner for each label and creates 
# two plots, auprc and auroc for training and test samples for each of the lncRNAs
perf_eval_ROC_PRC_perclass <- function(model_fit,
                                       my_dataset,
                                       my_label, 
                                       file_name, 
                                       neg_value = "neg",
                                       pos_value = "pos",
                                       class_name , 
                                       train_perf = F){
  # class_name is a character vector indicating the class (Tile owner) of each example 
  # train_perf : if true it will use the saved model predictions to get training performance (won't run the model on the full training data)
  stopifnot(length(my_label) == nrow(my_dataset),
            length(my_label) == length(class_name))
  require(ROCR)
  require(PRROC)
  if(train_perf){
    my_pred <- list()
    my_pred$predictions <- model_fit$predictions
  }else{
    my_pred <- predict(model_fit, data=my_dataset)
  }
  my_pred_list <- list()
  class_uniq <- sort(unique(class_name))
  aaskip <- numeric(0)
  for(i in 1:length(class_uniq)){
    if(sum(!is.na(my_pred$predictions[class_name %in% class_uniq[i] , 2])) > 2){
      my_prediction <- prediction(my_pred$predictions[class_name %in% class_uniq[i] , 2],
                                  my_label[class_name %in% class_uniq[i]])
      
      prc_rec <- performance(my_prediction, measure="prec", x.measure="rec")
      tpr_fpr <- performance(my_prediction, measure = "tpr", x.measure = "fpr")
      auc <- performance(my_prediction, measure = "auc")
      auc <- auc@y.values[[1]]
      index_Pos <- my_label[class_name %in% class_uniq[i]] == pos_value
      index_Neg <- my_label[class_name %in% class_uniq[i]] == neg_value
      aatestpred <- my_pred$predictions[class_name %in% class_uniq[i], 2]
      aaprc <- pr.curve(aatestpred[index_Pos], aatestpred[index_Neg], curve = T)
      my_pred_list[[i]] <- list()
      my_pred_list[[i]]$PRC_REC <- prc_rec
      my_pred_list[[i]]$TPR_FPR <- tpr_fpr
      my_pred_list[[i]]$AUROC <- auc
      my_pred_list[[i]]$AUPRC <- aaprc$auc.integral
    }else{
      my_pred_list[[i]] <- list()
      aaskip <- c(aaskip, i)
    }


  }
  names(my_pred_list) <- class_uniq
  #plot(aaprc)
  png(filename = file_name,    # create PNG for the heat map        
      width = 20*300,        # 5 x 300 pixels
      height = 20*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11) 
  par(mfrow = c(7,8), mar = c(3,4,4,2))
  for(i in 1:length(my_pred_list)){
    if(! (i %in% aaskip)){
      plot(my_pred_list[[i]]$TPR_FPR, ylab="TPR", xlab = "FPR", xlim = c(0,1), ylim = c(0,1),
           main = paste0(names(my_pred_list)[i], "\nauroc: ", (format(my_pred_list[[i]]$AUROC, digits = 2))))
      plot(my_pred_list[[i]]$PRC_REC, ylab="Precision", xlab = "Recall", xlim = c(0,1), ylim = c(0,1),
           main = paste0(names(my_pred_list)[i], "\nauprc: ", (format(my_pred_list[[i]]$AUPRC, digits = 2))))
    }else{
      plot(cbind(c(0,1), c(1,0)), ylab="", xlab = "", xlim = c(0,1), ylim = c(0,1),
           main = paste0(names(my_pred_list)[i], "\nauroc: ", "NOT enough observations"))
      plot(cbind(c(0,1), c(1,0)), ylab="", xlab = "", xlim = c(0,1), ylim = c(0,1),
           main = paste0(names(my_pred_list)[i], "\nauprc: ", "NOT enough observations"))
    }

  }
  #plot(aaprc)
  dev.off()
  return(my_pred_list)
}




############################################################################################################################################
# example
# for random partition eval
# aa_tst_1 <- perf_eval_ROC_PRC(model_fit = fit,
#                               my_dataset = aa_test_data,
#                               my_label = aa_test_data[, ncol(aa_test_data)], 
#                               file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/RF_test_perf_partition1_random.png")
############################################################################################################################################
################################################################################################################################################################
perf_eval_ROC_PRC_perclass_multimodel <- function(model_fit_list,
                                                  my_dataset,
                                                  my_label, 
                                                  file_name, 
                                                  neg_value = "neg",
                                                  pos_value = "pos",
                                                  class_name , 
                                                  train_perf = F){
  # class_name is a character vector indicating the class (Tile owner) of each example 
  # train_perf : if true it will use the saved model predictions to get training performance (won't run the model on the full training data)
  stopifnot(length(my_label) == nrow(my_dataset),
            length(my_label) == length(class_name))
  require(ROCR)
  require(PRROC)
  
  # index_Pos <- my_label == pos_value
  # index_Neg <- my_label == neg_value
  my_predict_list <- list()
  my_perf_list <- list()
  aaskip_list <- list()
  class_uniq <- sort(unique(class_name))
  for(cur_model in 1:length(model_fit_list)){
    my_pred_list <- list()
    if(train_perf){
      my_predict_list[[cur_model]] <- list()
      my_predict_list[[cur_model]]$predictions <- model_fit_list[[cur_model]]$predictions
      aana <- sum(is.na(my_predict_list[[cur_model]]$predictions))
      if(aana > 0){
        print(paste0("There are ", aana, " NA predicted values in the trained model ", names(model_fit_list)[cur_model], " . Will be set to zero..."))
        my_predict_list[[cur_model]]$predictions[is.na(my_predict_list[[cur_model]]$predictions)]<- 0
      }
      
    }else{
      my_predict_list[[cur_model]] <- predict(model_fit_list[[cur_model]], data=my_dataset)
    }
    aaskip_list[[cur_model]] <- numeric(0)
    for(i in 1:length(class_uniq)){
      if(sum(!is.na(my_predict_list[[cur_model]]$predictions[class_name %in% class_uniq[i] , 2])) > 2){
        my_prediction <- prediction(my_predict_list[[cur_model]]$predictions[class_name %in% class_uniq[i] , 2],
                                    my_label[class_name %in% class_uniq[i]])
        
        prc_rec <- performance(my_prediction, measure="prec", x.measure="rec")
        tpr_fpr <- performance(my_prediction, measure = "tpr", x.measure = "fpr")
        auc <- performance(my_prediction, measure = "auc")
        auc <- auc@y.values[[1]]
        index_Pos <- my_label[class_name %in% class_uniq[i]] == pos_value
        index_Neg <- my_label[class_name %in% class_uniq[i]] == neg_value
        aatestpred <- my_predict_list[[cur_model]]$predictions[class_name %in% class_uniq[i], 2]
        aaprc <- pr.curve(aatestpred[index_Pos], aatestpred[index_Neg], curve = T)
        my_pred_list[[i]] <- list()
        my_pred_list[[i]]$PRC_REC <- prc_rec
        my_pred_list[[i]]$TPR_FPR <- tpr_fpr
        my_pred_list[[i]]$AUROC <- auc
        my_pred_list[[i]]$AUPRC <- aaprc$auc.integral
      }else{
        my_pred_list[[i]] <- list()
        aaskip_list[[cur_model]] <- c(aaskip_list[[cur_model]], i)
      }
      
      
    }
    names(my_pred_list) <- class_uniq
    my_perf_list[[cur_model]] <- my_pred_list
  }
  names(my_perf_list) <- names(model_fit_list)
  names(aaskip_list) <- names(model_fit_list)
  

  #plot(aaprc)
  png(filename = file_name,    # create PNG for the heat map        
      width = 20*300,        # 5 x 300 pixels
      height = 20*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11) 
  par(mfrow = c(7,8), mar = c(3,4,4,2))
  for(i in 1:length(class_uniq)){
    for(cur_model in 1:length(model_fit_list)){ # going through models and plotting auroc
      if(! (i %in% aaskip_list[[cur_model]])){ # if there are enough predictions for this class
        if(cur_model == 1){ # plot if the first
          plot(my_perf_list[[cur_model]][[i]]$TPR_FPR,
               ylab="TPR", 
               xlab = "FPR",
               xlim = c(0,1), ylim = c(0,1),
               main = names(my_perf_list[[cur_model]])[i])
        }else{# add to plot if it's not the first model 
          lines(x = my_perf_list[[cur_model]][[i]]$TPR_FPR@x.values[[1]], 
                y =  my_perf_list[[cur_model]][[i]]$TPR_FPR@y.values[[1]], 
                col = cur_model)
        }
      }else{ # if there are not enough predictions for this class
        if(cur_model == 1){
          plot(cbind(c(0,1), c(1,0)), ylab="TPR", xlab = "FPR", xlim = c(0,1), ylim = c(0,1),
               main = names(my_perf_list[[cur_model]])[i])
        }
      }
    }
    if(i == 1){ # add legend only for the first panel
      legend(x = "bottomright",
             legend = names(model_fit_list),
             bty = "n",
             fill = c(1:length(model_fit_list)),
             cex = 0.7, x.intersp = 0.7, y.intersp = 0.7)
    }
    
    for(cur_model in 1:length(model_fit_list)){ # going through models and plotting auprc
      if(! (i %in% aaskip_list[[cur_model]])){ # if there are enough predictions for this class
        if(cur_model == 1){ # plot if the first
          plot(my_perf_list[[cur_model]][[i]]$PRC_REC,
               ylab="Precision", 
               xlab = "Recall",
               xlim = c(0,1), ylim = c(0,1),
               main = names(my_perf_list[[cur_model]])[i])
        }else{# add to plot if it's not the first model 
          lines(x = my_perf_list[[cur_model]][[i]]$PRC_REC@x.values[[1]], 
                y =  my_perf_list[[cur_model]][[i]]$PRC_REC@y.values[[1]], 
                col = cur_model)
        }
      }else{ # if there are not enough predictions for this class
        if(cur_model == 1){
          plot(cbind(c(0,1), c(1,0)), ylab="Precision", xlab = "Recall", xlim = c(0,1), ylim = c(0,1),
               main = names(my_perf_list[[cur_model]])[i])
        }
      }
    }
    if(i == 1){ # add legend only for the first panel
      legend(x = "topright",
             legend = names(model_fit_list),
             bty = "n",
             fill = c(1:length(model_fit_list)),
             cex = 0.7, x.intersp = 0.7, y.intersp = 0.7)
    }
    

  }
  #plot(aaprc)
  dev.off()
  return(my_perf_list)
}
############################################################################################################################################
############################################################################################################################################
# example
# my_model_fit_list <- list()
# load("RF_partition2_random_caseWght.RData")
# load("Learned_models/Partition_2_modified_dataset_and_removedColumnNames___chunk___distance___RFmodel.RData")
# my_model_fit_list[[1]] <- RF_partition2_chunk_classWght
# my_model_fit_list[[2]] <- my_RF_model
# names(my_model_fit_list) <- c("random","random_noDis")
# aanampas <- paste0(names(my_model_fit_list), collapse = "__vs__")
# load("~/Documents/Shayan/BioInf/lncRNA/MM9_1kb_tiled_owner_3Mfilter.nosync.RData")
# aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_test_data), names(MM9_1kb_tiled_owner_3Mfilter))]
# aa_pr_multimodel <- perf_eval_ROC_PRC_perclass_multimodel(model_fit_list = my_model_fit_list,
#                                          my_dataset = aa_test_data,
#                                          my_label = aa_test_data[, ncol(aa_test_data)],
#                                          file_name = paste0(plot_strore_address, aanampas, "___test_all_perclass.png"),
#                                          train_perf = F,
#                                          class_name = aaclass_name,
#                                          neg_value = "neg",
#                                          pos_value = "pos")
# aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_train_data), names(MM9_1kb_tiled_owner_3Mfilter))]
# aa_pr_multimodeltr <- perf_eval_ROC_PRC_perclass_multimodel(model_fit_list = my_model_fit_list,
#                                                  my_dataset = aa_train_data,
#                                                  my_label = aa_train_data[, ncol(aa_train_data)],
#                                                  file_name = paste0(plot_strore_address, aanampas, "___train_all_perclass.png"),
#                                                  train_perf = T,
#                                                  class_name = aaclass_name,
#                                                  neg_value = "neg",
#                                                  pos_value = "pos")


############################################################################################################################################
############################################################################################################################################
performance_viz_chromosome <- function(all_tile_label=MM9_1kb_tiled_owner_labels_binary_3Mfilter,
                                       show_points = c("all"), 
                                       subset_points= c( "train", "test", "valid"),
                                       partition_dataset,
                                       predicted_pos_probilbilty,
                                       positive_thersh = 0.5,
                                       specific_chr = character(0),
                                       my_file_name = "performance_per_chr.png",
                                       lncRNA_coor_df= lncRNA_chosen_gt1k_uniqTiles,
                                       chr_sizes_file="~/Documents/Shayan/BioInf/lncRNA/mm9_chr_sizes.txt"){
  # all_tile_label : is a binary vector named with the names of the universe of all tiles available
  # show_points : is character: either "all", or "used" , if all chosen, the points will be shown regardless of if they have been used by the model,
  #  but if 'used' is chosen, points are shown only if they appear in one of the training, test, or validation sets.
  # subset_points: is a character vector if length greater than or equal to 1, can be any subset of these c( "train", "test", "valid")
  # partition_dataset: is a dataframe that has at least the following columns: tile_name,owner, label, dataset
  # is a nemed numeric vector , containing the pobability of positive label predicition. (all numbers between 0, 1)
  # positive_thersh, the threshold on positive probability below which is considered to be prediced negative
  # specific_chr : a character vector containing the names of chromosomes for which you want the plots, if left empty uses all chrs
  stopifnot(length(names(predicted_pos_probilbilty)) == length(predicted_pos_probilbilty),
            all(partition_dataset$tile_name  %in% names(all_tile_label)),
            all(names(predicted_pos_probilbilty) %in% partition_dataset$tile_name),
            all(predicted_pos_probilbilty <= 1), 
            all(predicted_pos_probilbilty >= 0),
            is.character(partition_dataset$tile_name), 
            all(subset_points %in% c( "train", "test", "valid")), 
            length(show_points) == 1,
            show_points %in% c("all", "used")  )
  label_pos <- "pos"
  label_neg <- "neg"
  
  partition_dataset$label[partition_dataset$label == 0] <- label_neg
  partition_dataset$label[partition_dataset$label == 1] <- label_pos
  
  mm9_chr_size <-read.table(chr_sizes_file, 
                            header = F,
                            stringsAsFactors = F)
  train_test_valid_heightdf <- data.frame(name = c("train", "valid", "test"),
                                          neg = c(0, 0.4, 0.8), 
                                          pos = c(1.2, 1.6, 2), 
                                          col = c(4,3,2),
                                          pch = c(6,6,6),
                                          cex = c(0.5, 0.5, 0.5))
  lncRNA_height <- 2.4
  mu_shift <- 0.1 # shift due to misclassification
  
  ########  
  # adding predicions o the partition dataset
  partition_dataset$predicted_label <- array(dim = nrow(partition_dataset))
  
  aapredpos <- names(predicted_pos_probilbilty)[predicted_pos_probilbilty >= positive_thersh]
  aapredneg <- names(predicted_pos_probilbilty)[predicted_pos_probilbilty < positive_thersh]
  
  partition_dataset$predicted_label[match(aapredpos, partition_dataset$tile_name)] <- label_pos
  partition_dataset$predicted_label[match(aapredneg, partition_dataset$tile_name)] <- label_neg
  aag <- partition_dataset$predicted_label == partition_dataset$label
  
  partition_dataset$label_agreement[!aag] <- -1
  partition_dataset$label_agreement[aag] <- 1
  partition_dataset$label_agreement[is.na(aag)] <- 0
  
  
  #############
  # constructing the base for all case
  full_pos_neg_layer <- list()
  
  all_tilename_neg <- names(all_tile_label)[all_tile_label == 0]
  aa_tilenamesp_neg <- strsplit(all_tilename_neg, "_")
  aa_tilenamesp1_neg<- unlist(lapply(aa_tilenamesp_neg, "[[", 1))
  aa_tilenamesp2_neg <- as.numeric(unlist(lapply(aa_tilenamesp_neg, "[[", 2)))
  print(sum(is.na(aa_tilenamesp2_neg)))
  
  all_tilename_pos <- names(all_tile_label)[all_tile_label == 1]
  aa_tilenamesp <- strsplit(all_tilename_pos, "_")
  aa_tilenamesp1<- unlist(lapply(aa_tilenamesp, "[[", 1))
  aa_tilenamesp2<- as.numeric(unlist(lapply(aa_tilenamesp, "[[", 2)))
  print(sum(is.na(aa_tilenamesp2)))
  
  if(length(specific_chr) >0){
    stopifnot(all(specific_chr %in% aa_tilenamesp1))
    aachr_size2 <- mm9_chr_size[mm9_chr_size$V1 %in% specific_chr,]
  }else{
    aachr_size2 <- mm9_chr_size[mm9_chr_size$V1 %in% aa_tilenamesp1,]
  }
  
  for(i in 1:nrow(aachr_size2)){
    full_pos_neg_layer[[i]] <- array(dim = ceiling(aachr_size2$V2[i]/1000))
    names(full_pos_neg_layer[[i]]) <- c(1:length(full_pos_neg_layer[[i]]))
    full_pos_neg_layer[[i]][aa_tilenamesp2[which(aa_tilenamesp1 %in% aachr_size2$V1[i])]] <- 1
    full_pos_neg_layer[[i]][aa_tilenamesp2_neg[which(aa_tilenamesp1_neg %in% aachr_size2$V1[i])]] <- 0
  }
  names(full_pos_neg_layer) <- aachr_size2$V1
  ###########
  # constructing the second layer for all cases: the one specified by subset_points
  partition_sp <- strsplit(partition_dataset$tile_name, "_")
  partition_sp1 <- unlist(lapply(partition_sp, "[[", 1))
  partition_sp2 <- as.numeric(unlist(lapply(partition_sp, "[[", 2)))
  print(sum(is.na(partition_sp2)))
  layer_2 <- list() # for the height of points in datasets
  for(i in 1:nrow(aachr_size2)){
    layer_2[[i]] <-   array(dim =  ceiling(aachr_size2$V2[i]/1000))
    names(layer_2[[i]]) <- c(1:length(layer_2[[i]]))
  }
  names(layer_2) <- aachr_size2$V1
  
  layer_3 <- list() # for the color of points in datasets
  for(i in 1:nrow(aachr_size2)){
    layer_3[[i]] <-   array(dim =  ceiling(aachr_size2$V2[i]/1000))
    names(layer_3[[i]]) <- c(1:length(layer_3[[i]]))
  }
  names(layer_3) <- aachr_size2$V1
  
  
  performance_matrix <- matrix(nrow = 4, ncol = length(layer_2))
  rownames(performance_matrix) <- c("TP", "TN", "FP", "FN")
  colnames(performance_matrix) <- names(layer_2)
  for(i in 1:length(layer_2)){
    #print(sum((partition_sp1 %in%  names(layer_2)[i]) ))
    #print(sum(partition_dataset$label == label_pos))
    #print(sum((partition_dataset$predicted_label== label_pos), na.rm = T ))
    #print("\n###############\n")
    aatp <- sum(((partition_sp1 %in%  names(layer_2)[i]) & (partition_dataset$label == label_pos) & (partition_dataset$predicted_label== label_pos)), na.rm = T)
    aatn <- sum(((partition_sp1 %in%  names(layer_2)[i]) & (partition_dataset$label == label_neg) & (partition_dataset$predicted_label== label_neg)), na.rm = T)
    aafp <- sum(((partition_sp1 %in%  names(layer_2)[i]) & (partition_dataset$label == label_neg) & (partition_dataset$predicted_label== label_pos)), na.rm = T)
    aafn <- sum(((partition_sp1 %in%  names(layer_2)[i]) & (partition_dataset$label == label_pos) & (partition_dataset$predicted_label== label_neg)), na.rm = T)
    
    performance_matrix[1, i] <- aatp
    performance_matrix[2, i] <- aatn
    performance_matrix[3, i] <- aafp
    performance_matrix[4, i] <- aafn
  }
  ########
  # construct the second layer
  for(cur_subset in 1:length(subset_points)){ # between the chosen subsets: train, test, valid
    for(cur_chr in 1:nrow(aachr_size2)){
      aa_whch_pos <-  which((partition_sp1 %in% aachr_size2$V1[cur_chr]) &
                              (partition_dataset$label %in% label_pos) &
                              (partition_dataset$dataset %in% subset_points[cur_subset])) 
      #print(length(aa_whch_pos))
      layer_2[[cur_chr]][partition_sp2[aa_whch_pos]] <- rnorm(n = length(aa_whch_pos), mean = train_test_valid_heightdf$pos[train_test_valid_heightdf$name == subset_points[cur_subset]], sd = 0.02)
      layer_2[[cur_chr]][partition_sp2[aa_whch_pos]] <- layer_2[[cur_chr]][partition_sp2[aa_whch_pos]] + mu_shift * partition_dataset$label_agreement[aa_whch_pos]
      layer_3[[cur_chr]][partition_sp2[aa_whch_pos]] <- train_test_valid_heightdf$col[train_test_valid_heightdf$name == subset_points[cur_subset]]
      
      aa_whch_neg <-  which((partition_sp1 %in% aachr_size2$V1[cur_chr]) & 
                              (partition_dataset$label == label_neg) & 
                              (partition_dataset$dataset %in% subset_points[cur_subset])) 
      layer_2[[cur_chr]][partition_sp2[aa_whch_neg]] <-   + rnorm(n = length(aa_whch_neg), mean =train_test_valid_heightdf$neg[train_test_valid_heightdf$name == subset_points[cur_subset]] , sd = 0.02)
      layer_2[[cur_chr]][partition_sp2[aa_whch_neg]] <- layer_2[[cur_chr]][partition_sp2[aa_whch_neg]] + mu_shift * partition_dataset$label_agreement[aa_whch_neg]
      layer_3[[cur_chr]][partition_sp2[aa_whch_neg]] <- train_test_valid_heightdf$col[train_test_valid_heightdf$name == subset_points[cur_subset]]
      
    }
  }
  
  ##### add lncRNA position
  for(i in 1:nrow(aachr_size2)){
    aaw <- which(lncRNA_coor_df$chromosome_name %in% aachr_size2$V1[i])
    for(j in 1:length(aaw)){
      aawch <- c(floor((lncRNA_coor_df$start_position[aaw[j]])/1000):ceiling((lncRNA_coor_df$end_position[aaw[j]])/1000))
      layer_2[[i]][aawch] <- lncRNA_height
      layer_3[[i]][aawch] <- "pink"
    }
  }
  #########
  
  png(filename = my_file_name,       
      width = 30*300,        # 5 x 300 pixels
      height = (3 * length(layer_2)*300),
      res = 300,            # 300 pixels per inch
      pointsize = 10) 
  par(mfrow = c(length(layer_2) ,1), mar = c(3,4,3,4))
  
  if(show_points == "all"){
    for(i in 1:length(full_pos_neg_layer)){
      plot(full_pos_neg_layer[[i]],
           main = paste0(names(full_pos_neg_layer)[i],"_" ,aachr_size2$V2[i],
                         "\nTP: ",performance_matrix[1,i],  
                         " ,TN: ",performance_matrix[2,i], 
                         " ,FP: ",performance_matrix[3,i],
                         " ,FN: ",performance_matrix[4,i]),
           las = 2,
           xlab = "genomic position", 
           xaxt = "n", ylim = c(-0.2,2.5),
           pch = 19, 
           cex = 0.3)
      
      points(layer_2[[i]], pch = 8, col = layer_3[[i]], cex = 0.4)
      abline(h = c(train_test_valid_heightdf$neg, train_test_valid_heightdf$pos, lncRNA_height), col ="lightgrey", lwd = 0.6, lty = 3)
      # abline(h = c(train_test_valid_heightdf$neg + mu_shift, train_test_valid_heightdf$pos, lncRNA_height), col ="lightgrey", lwd = 0.3, lty = 3)
    }
    
  }else{
    for(i in 1:length(layer_2)){
      plot(layer_2[[i]],
           main = paste0(names(layer_2)[i],"_" ,aachr_size2$V2[i],
                         "\nTP: ",performance_matrix[1,i],  
                         " ,TN: ",performance_matrix[2,i], 
                         " ,FP: ",performance_matrix[3,i],
                         " ,FN: ",performance_matrix[4,i]),
           las = 2, xlab = "genomic position", xaxt = "n", ylim = c(-0.2,2.5),
           pch = 8, cex = 0.4, col = layer_3[[i]])
      abline(h = c(train_test_valid_heightdf$neg, train_test_valid_heightdf$pos, lncRNA_height), col ="lightgrey", lwd = 0.6, lty = c(2,3,4,2,3,4,5))
    }
    
  }
  dev.off()
  return(list(Layer_1 = full_pos_neg_layer,
              Layer_2 = layer_2,
              Layer_3 = layer_3,
              Perf_mat = performance_matrix))
  
}
################################################################################################################################################################################################
################################################################################################################################################################################################
# example
# aamy_pred <- predict(rf_model_partition1_chunk_varimp_purity, data=aa_test_data)
# aapreddd <-  aamy_pred$predictions[, 2]
# names(aapreddd) <- rownames(aa_test_data)
# aatst <- performance_viz_chromosome(all_tile_label=MM9_1kb_tiled_owner_labels_binary_3Mfilter,
#                                     show_points = c("all"), 
#                                     subset_points= c( "train", "test", "valid"),
#                                     partition_dataset = Partition_1_chunk_dfs,
#                                     predicted_pos_probilbilty = aapreddd,
#                                     positive_thersh = 0.2,
#                                     specific_chr = character(0),
#                                     my_file_name = "~/Documents/Shayan/BioInf/lncRNA/plots/performance_per_chr_chunk_partition1_thresh2.png",
#                                     lncRNA_coor_df= lncRNA_chosen_gt1k_uniqTiles,
#                                     chr_sizes_file="~/Documents/Shayan/BioInf/lncRNA/mm9_chr_sizes.txt")
################################################################################################################################################################################################
visualize_topN_scatter <- function(my_dataset, my_partition, my_model, topN=10, my_file_name = "top_feature_inspection.png", box_only = F){
  require(ranger)
  # form pos-neg map for scatter plot
  aa_all_name <- my_partition$tile_name
  aa_tile_split <- strsplit(aa_all_name, "_")
  stopifnot(all(unlist(lapply(aa_tile_split, length)) == 2))
  aachr <- unlist(lapply(aa_tile_split, "[[", 1))
  aanu <- as.numeric(unlist(lapply(aa_tile_split, "[[", 2)))
  aachr_uniq <- unique(aachr)
  pos_neg_map <- matrix(nrow = 0, ncol = 2)
  colnames(pos_neg_map) <- c("pos", "neg")
  for(cur_chr in 1:length(aachr_uniq)){
    print("cur_chr")
    print(cur_chr)
    aa_w_pos <- which((aachr %in% aachr_uniq[cur_chr]) & (my_partition$label == 1))
    aa_w_neg <- which((aachr %in% aachr_uniq[cur_chr]) & (my_partition$label == 0))
    if((length(aa_w_pos)) == 0 | (length(aa_w_neg) == 0)){
      next()
    }
    print("length(aa_w_neg)")
    print(length(aa_w_neg))
    print("length(aa_w_pos)")
    print(length(aa_w_pos))
    
    aa_neg_clos <- numeric(length = length(aa_w_neg))
    for(cur_neg in 1:length(aa_w_neg)){
      aa_neg_clos[cur_neg] <- aa_w_pos[which.min(abs(aanu[cur_neg] - aanu[aa_w_pos]))]
    }
    
    aarem_pos <- setdiff(aa_w_pos, aa_neg_clos)
    if(length(aarem_pos) > 1){
      print("there are unmatched positives")
      print("length(aarem_pos)")
      print(length(aarem_pos))
      print("length(rep(aarem_pos[length(aarem_pos)], length(aarem_pos)))")
      print(length(rep(aarem_pos[length(aarem_pos)], length(aarem_pos))))
      aa1 <- cbind(aa_neg_clos, aa_w_neg)
      aa2 <- cbind(aarem_pos ,rep(aarem_pos[length(aarem_pos)], length(aarem_pos)))
      print("ncol(aa1)")
      print(ncol(aa1))
      print("ncol(aa2)")
      print(ncol(aa2))
      aa_curmap <- rbind(aa1,aa2)
      print("after_map")
      
    }else{
      aa_curmap <- cbind(aa_neg_clos, 
                         aa_w_neg)
    }
    pos_neg_map <- rbind(pos_neg_map, aa_curmap)
  }
  
  # extract important features
  
  variable_importance <- importance(my_model)
  my_top_N <- sort(variable_importance, decreasing = T, index.return = T)$ix[1:topN]
  topN_name <- my_name_dic[my_top_N, 1]
  
  # plot
  if((floor(sqrt(topN)) * ceiling(sqrt(topN))) >= topN){
    aaprow <- c(ceiling(sqrt(topN)), floor(sqrt(topN)))
  }else{
    aaprow <- c(ceiling(sqrt(topN)), ceiling(sqrt(topN)))
  }
  if(!box_only){
    png(filename = my_file_name,    # create PNG for the heat map        
        width = 10*300,        # 5 x 300 pixels
        height = 10*300,
        res = 300,            # 300 pixels per inch
        pointsize = 11) 
    par(mfrow = aaprow, mar = c(4,4,3,3))
    on.exit(dev.off())
    for(cur_f in 1:length(my_top_N)){
      aaplot <- cbind(my_dataset[pos_neg_map[, 1],my_top_N[cur_f]], 
                      my_dataset[pos_neg_map[, 2],my_top_N[cur_f]])
      if(length(aarem_pos) > 1){
        aaplot[pos_neg_map[, 2] %in% aarem_pos[length(aarem_pos)], 2] <- -1
      }
      
      plot(aaplot,
           xlab = "+", ylab = "-",
           main = paste0(topN_name[cur_f], "___",my_partition$owner[1]),
           pch = 19, 
           cex = 0.4, 
           xlim = range(aaplot),
           ylim =range(aaplot) )
      abline(a = 0, b = 1, col = 2)
      abline(h = 0, col = 3)
    }
  }
  
  print("########################")
  
  png(filename = paste0(my_file_name, "__boxplot"),    # create PNG for the heat map        
      width = 10*300,        # 5 x 300 pixels
      height = 10*300,
      res = 300,            # 300 pixels per inch
      pointsize = 11) 
  par(mfrow = aaprow, mar = c(4,4,3,3))
  on.exit(dev.off())
  
  for(cur_f in 1:length(my_top_N)){
    aaplot <- cbind(my_dataset[pos_neg_map[, 1],my_top_N[cur_f]], 
                    my_dataset[pos_neg_map[, 2],my_top_N[cur_f]])
    if(length(aarem_pos) > 1){
      aaplot[pos_neg_map[, 2] %in% aarem_pos[length(aarem_pos)], 2] <- -1
    }
    
    boxplot(my_dataset[, my_top_N[cur_f]]~my_partition$label,
            main = paste0(topN_name[cur_f], "___",my_partition$owner[1]), outline = F)
  }
  return(pos_neg_map)
}

#####
# example

# 
# aa_file_datast <- list.files("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected")
# aa_file_model  <- list.files("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_permutation")
# aaspl_datast <- unlist(lapply(strsplit(unlist(lapply(strsplit(aa_file_datast, "\\."), "[[", 1)), "_"), "[[", 5))
# 
# for(i in 1:length(aa_file_datast)){
#   print(i)
#   load(paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/", aa_file_datast[i]))
#   load(paste0("~/Documents/Shayan/BioInf/lncRNA/Learned_models_p6_permutation/Partition_6_modified_dataset_", aaspl_datast[i], "___rand___full___RFmodel.RData"))
#   AAtst <- visualize_topN_scatter(my_dataset = my_Dataset,
#                                   my_partition = my_partition_rand,
#                                   my_model = my_RF_model,
#                                   topN=20, 
#                                   my_file_name = paste0("~/Documents/Shayan/BioInf/lncRNA/plots/Feature_top20/", aaspl_datast[i], "_top20.png"), 
#                                   box_only = T)
#   
# }


perf_eval_wrapper <- function(Learned_model_folder,
                              lncRNA_dataset_folder,
                              model_name=character(0),
                              lncRNA_name = character(0),
                              jobfile_address,
                              aaplot_strore_address,
                              prename = "all_comp",
                              only_imporatnce = F,
                              save_pred = F,
                              filtering_folder = "~/Documents/Shayan/BioInf/lncRNA/partition_6_Dec_TOBE_exp_FILTERED",
                              plot_importance = T,
                              plot_roc_prc=T,
                              my_partition_Exp_features_address = "~/Documents/Shayan/BioInf/lncRNA/partition_6_expression_features.RData"){
  require(ranger)
  require(RColorBrewer)
  require(ROCR)
  require(PRROC)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  aa_dataset <- list.files(lncRNA_dataset_folder)
  aa_files <- list.files(Learned_model_folder)
  aa_lnc_list <- list()
  aa_dataset_list <- list()
  if(length(lncRNA_name) > 0){
    
    for(i in 1:length(lncRNA_name)){
      aa_lnc_list[[i]] <- aa_files[grep(lncRNA_name[i], aa_files)]
      aa_dataset_list[[i]] <- aa_dataset[grep(lncRNA_name, aa_dataset)]
    }
  }else{ 
    aa_lnc_list[[1]] <- aa_files
    aattt <- unlist(lapply(strsplit(unlist(lapply(strsplit(aa_lnc_list[[1]], "___"), "[[", 1)), "_"), "[[", 5))
    lncRNA_name <- unique(aattt)
    print(lncRNA_name)
    stopifnot(length(lncRNA_name) == 1)
    #print(lncRNA_name)
    aa_dataset_list[[1]] <- aa_dataset[grep(lncRNA_name, aa_dataset)]
    print("aa_dataset")
    print(aa_dataset_list[[1]])
  }
  
  my_model_names <- list()
  my_model_type <- list() # chunk or rand
  for(i in 1:length(aa_lnc_list)){
    my_model_names[[i]] <- unlist(lapply(strsplit(aa_lnc_list[[i]], "___"), "[[", 3))
    my_model_type[[i]] <- unlist(lapply(strsplit(aa_lnc_list[[i]], "___"), "[[", 2))
    
    
  }
  
  if(length(model_name) > 0){
    for(i in 1:length(aa_lnc_list)){
      aa_lnc_list[[i]] <- aa_lnc_list[[i]][my_model_names[[i]] %in% model_name]
      my_model_type[[i]] <- my_model_type[[i]][my_model_names[[i]] %in% model_name]
      my_model_names[[i]] <- my_model_names[[i]][my_model_names[[i]] %in% model_name]
      
    }
  }
  
  
  aa_all_job <- readLines(jobfile_address)
  # for each lncRNA divide the models into rand and chunk groups and 
  for(cur_lnc in 1:length(aa_lnc_list)){
    #if(!only_imporatnce){
    load(paste0(lncRNA_dataset_folder,"/",aa_dataset_list[[i]]))
    aanammm <- rownames(my_Dataset)
    aanammmsp <- unlist(lapply(strsplit(aanammm, "\\."), "[[", 1))
    rownames(my_Dataset) <- aanammmsp
    #}
    
    all_partitions <- unique(my_model_type[[cur_lnc]])
    my_pt_model_list <- list()
    my_pt_model_name_list <- list()
    for(cur_pt in 1:length(all_partitions)){
      my_pt_model_list[[cur_pt]] <- aa_lnc_list[[cur_lnc]][my_model_type[[cur_lnc]] %in% all_partitions[cur_pt]]
      my_pt_model_name_list[[cur_pt]] <- my_model_names[[cur_lnc]][my_model_type[[cur_lnc]] %in% all_partitions[cur_pt]]
      stopifnot(length(my_pt_model_name_list[[cur_pt]]) == length(my_pt_model_list[[cur_pt]]))
    }

    
    
    aa_curjob <- aa_all_job[grep(pattern = lncRNA_name, x = aa_all_job)]
    my_pred_list_perpt <- list()

    perf_list_perpt <- list()
    
    if(save_pred){
      pred_list_perpt <- list()
    }
    featImp_perpt_list <- list()
    
    for(cur_pt in 1:length(all_partitions)){
      print(all_partitions[cur_pt])
      perf_list_perpt[[cur_pt]] <- list()
      pred_list_perpt[[cur_pt]] <- list()
      featImp_perpt_list[[cur_pt]] <- list()
      if(length(my_pt_model_list[[cur_pt]]) > 0){
        aa_curjob_pt <- aa_curjob[grep(all_partitions[cur_pt], aa_curjob)]
        #print("length(aa_curjob_pt)")
        #print(length(aa_curjob_pt))
        #print(aa_curjob_pt)
        aawh <- which(colnames(my_partition) == all_partitions[cur_pt])
        my_partition_df <- my_partition[, c(1,3,2, aawh)]
        colnames(my_partition_df)[4] <- "dataset"
        aatrainind <- which(my_partition_df$dataset %in% 1)
        
        for(cur_mod in 1:length(my_pt_model_list[[cur_pt]])){
          print(my_pt_model_list[[cur_pt]][cur_mod])
          load(paste0(Learned_model_folder, "/", my_pt_model_list[[cur_pt]][cur_mod]))
          aagrep <- paste0("\\b", my_pt_model_name_list[[cur_pt]][cur_mod], "\\b")
         # print(aagrep)
          aa_curjob_model <- aa_curjob_pt[grep(pattern = aagrep, x = aa_curjob_pt )]
         # print(aa_curjob_model)
          aaspl <- unlist(strsplit(aa_curjob_model, " "))
          aa_RBP_filtering <- aaspl[grep("filter", aaspl)]
          #print(aa_RBP_filtering)
          if(!only_imporatnce){
            if(aa_RBP_filtering != "no_filter"){
              aa_RBP_filtering2 <- unlist(strsplit(aa_RBP_filtering, "\\/"))
              aa_RBP_filtering2 <- aa_RBP_filtering2[length(aa_RBP_filtering2)]
              aa_RBP_filtering2 <- paste0(filtering_folder, "/", aa_RBP_filtering2)
              my_filter_index <- as.numeric(read.table(aa_RBP_filtering2)$V1)
              #aagn <- grep(pattern = "TXGRO", x = my_chunk_model_names[cur_mod])
              feat_to_filterby <- c("GROseq","PROseq")
              my_Dataset_tmp <- Filter_RBP(inp_dataset = my_Dataset,
                                           exp_feature_name = feat_to_filterby,
                                           inp_nameDic = my_name_dic,
                                           feature_index=my_filter_index,
                                           partition_Exp_features_address = my_partition_Exp_features_address)
              #print("after Filter_RBP")
              #print(my_chunk_model_names[cur_mod])
            }else{
              my_Dataset_tmp <- my_Dataset
            } # end of going over filters
            tmp_test_dataset <- my_Dataset_tmp[-aatrainind,]
            rm(my_Dataset_tmp)
            my_label <- tmp_test_dataset[, ncol(tmp_test_dataset)]
            index_Pos <- my_label == "pos"
            index_Neg <- my_label == "neg"
            #my_pred_list_chunk[[i]]
            #print("before predict")
            tmp_pred <- predict(my_RF_model, data=tmp_test_dataset)
            if(save_pred){
              pred_list_perpt[[cur_pt]][[cur_mod]] <- tmp_pred
            }
            perf_list_perpt[[cur_pt]][[cur_mod]] <- list()
            #print("before prediction")
            my_prediction <- prediction(tmp_pred$predictions[, 2],my_label)
            #print("after prediction")
            
            perf_list_perpt[[cur_pt]][[cur_mod]]$prc_rec <- performance(my_prediction, measure="prec", x.measure="rec")
            #print("after performance prc_rec")
            perf_list_perpt[[cur_pt]][[cur_mod]]$tpr_fpr <- performance(my_prediction, measure = "tpr", x.measure = "fpr")
            #print("after performance tpr_fpr")
            auc <- performance(my_prediction, measure = "auc")
            #print("after performance auroc")
            perf_list_perpt[[cur_pt]][[cur_mod]]$auc <- auc@y.values[[1]]

            perf_list_perpt[[cur_pt]][[cur_mod]]$prc <- pr.curve(tmp_pred$predictions[index_Pos, 2], tmp_pred$predictions[index_Neg, 2], curve = F)
          } # end of only importance
          featImp_perpt_list[[cur_pt]][[cur_mod]] <- importance(my_RF_model)
          names(featImp_perpt_list[[cur_pt]][[cur_mod]]) <- my_name_dic[match(names(featImp_perpt_list[[cur_pt]][[cur_mod]]), my_name_dic[, 2]),1]

        } # end of going models in the current partition
        if(!only_imporatnce){
          names(perf_list_perpt[[cur_pt]]) <- my_pt_model_name_list[[cur_pt]]
          if(save_pred){
            names(pred_list_perpt[[cur_pt]]) <- my_pt_model_name_list[[cur_pt]]
          }
        }
      } # end of if partition has any members
    } # end of going over partitions
    
    if(plot_importance){
      for(cur_pt in 1:length(all_partitions)){
        print("length(featImp_perpt_list[[cur_pt]])")
        print(length(featImp_perpt_list[[cur_pt]]))
        png(filename = paste0(aaplot_strore_address,"/Importance_",lncRNA_name[cur_lnc],"_", all_partitions[cur_pt], "_",prename,".png"),         
            width = length(featImp_perpt_list[[cur_pt]])* 1 *300,        # 5 x 300 pixels
            height = 40*300,
            res = 300,            # 300 pixels per inch
            pointsize = 10)
        par(mfrow = c(5,(ceiling(length(featImp_perpt_list[[cur_pt]])/5))), mar = c(4,9,4,4))
        for(i in 1:length(featImp_perpt_list[[cur_pt]]) ){
          barplot(sort(featImp_perpt_list[[cur_pt]][[i]], decreasing = T)[1:50], horiz = T, las =2, main = my_pt_model_name_list[[cur_pt]][i])
        }
        dev.off()
      }

    }
    ###############
    if(plot_roc_prc){
      for(cur_pt in 1:length(all_partitions)){
        png(filename = paste0(aaplot_strore_address,"/perf_",lncRNA_name[cur_lnc],"_", all_partitions[cur_pt], "_",prename,"__test.png"),    
            width = 16*300,        
            height = 4*300,
            res = 300,           
            pointsize = 11)
        par(mfrow = c(1,2), mar = c(4,5,4,15), xpd = T)
        for(ii in 1:length(perf_list_perpt[[cur_pt]])){
          if(ii == 1){
            plot(perf_list_perpt[[cur_pt]][[ii]]$tpr_fpr, ylab="TPR", xlab = "FPR",
                 # main = paste0("auroc: ", (format(perf_list[[i]]$auc, digits = 2))),
                 xlim = c(0,1), ylim = c(0,1), col = col_vector[ii])
          }else{

            lines(x = perf_list_perpt[[cur_pt]][[ii]]$tpr_fpr@x.values[[1]],
                  y = perf_list_perpt[[cur_pt]][[ii]]$tpr_fpr@y.values[[1]], col = col_vector[ii])
          }
        }
        aa_all_auc<- unlist(lapply(perf_list_perpt[[cur_pt]], "[[", 3))
        legend(x = "topright", inset=c(-0.6,0),xpd = T,
               legend = paste(names(perf_list_perpt[[cur_pt]]),
                              format(aa_all_auc, digits = 2), 
                              sep = "___"), 
               fill = col_vector[c(1:ii)],
               bty = "n",
               cex = 0.5, x.intersp = 0.6, y.intersp = 0.6)
        
        
        
        for(ii in 1:length(perf_list_perpt[[cur_pt]])){
          if(ii == 1){
            plot(perf_list_perpt[[cur_pt]][[ii]]$prc_rec, ylab="Precision", xlab = "Recall",
                 #main = paste0("auprc: ", (format(aaprc$auc.integral, digits = 2))),
                 xlim = c(0,1),
                 ylim = c(0,1), col = col_vector[ii])
          }else{
            lines(x = perf_list_perpt[[cur_pt]][[ii]]$prc_rec@x.values[[1]],
                  y=perf_list_perpt[[cur_pt]][[ii]]$prc_rec@y.values[[1]] , col = col_vector[ii])
          }
        }
        aa_all_prc<- lapply(perf_list_perpt[[cur_pt]], "[[", 4)
        aa_all_prc <- unlist(lapply(aa_all_prc, "[[", "auc.integral") )
        legend(x = "topright", inset=c(-0.6,0), xpd = T,
               legend = paste(names(perf_list_perpt[[cur_pt]]), 
                              format(aa_all_prc, digits = 2), sep = "___"),
               
               bty = "n",
               fill = col_vector[c(1:ii)],
               cex = 0.5, x.intersp = 0.6, y.intersp = 0.6)
        dev.off()
      } # end of loop over paritions

    }

    
  } # end of loop over lncRNAs
  ret_lst <- list(performance=perf_list_perpt ,
                  importance=featImp_perpt_list)
  if(save_pred){
    ret_lst$prediction <- pred_list_perpt

  }


  return(ret_lst)
}




######### example
# source("~/Documents/Shayan/BioInf/lncRNA/RF_evaluation_functions.R")
# perf_eval_wrapper(Learned_model_folder = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Results/Learned_models/Gm14820",
#                   lncRNA_dataset_folder = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_dataset",
#                   model_name=character(0),
#                   lncRNA_name = character(0),
#                   jobfile_address = "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full.job",
#                   aaplot_strore_address = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_Res_comparision",
#                   filtering_folder = "~/Documents/Shayan/BioInf/lncRNA/partition_6_Dec_TOBE_exp_FILTERED")




# function to gather performnace data for various models

perf_gather <- function(Learned_model_folder, 
                        model_name = character(0), 
                        lncRNA_name = character(0),
                        return_imp = T,
                        return_pred=F,
                        return_per_full=F,
                        train_perf_included=F,
                        file_pattern = "*___perf_pred_varIMP.RData"){
#  require(ranger)
  #require(RColorBrewer)
 # require(ROCR)
 # require(PRROC)
  # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  # Learned_model_folder : address to folder containing performance results for a particular lncRNA
  aa_files <- list.files(Learned_model_folder, pattern = file_pattern)
  aa_lnc_list <- list()
  #aa_dataset_list <- list()
  if(length(lncRNA_name) > 0){
    
    for(i in 1:length(lncRNA_name)){
      aa_lnc_list[[i]] <- aa_files[grep(lncRNA_name[i], aa_files)]
      #aa_dataset_list[[i]] <- aa_dataset[grep(lncRNA_name, aa_dataset)]
    }
  }else{ 
    aa_lnc_list[[1]] <- aa_files
    aattt <- unlist(lapply(strsplit(unlist(lapply(strsplit(aa_lnc_list[[1]], "___"), "[[", 1)), "_"), "[[", 5))
    lncRNA_name <- unique(aattt)
    print(lncRNA_name)
    stopifnot(length(lncRNA_name) == 1)
    # aa_dataset_list[[1]] <- aa_dataset[grep(lncRNA_name, aa_dataset)]
    # print("aa_dataset")
    # print(aa_dataset_list[[1]])
  }
  
  my_model_names <- list()
  my_model_type <- list() # RCV1:5, CCV1:5
  for(i in 1:length(aa_lnc_list)){
    my_model_names[[i]] <- unlist(lapply(strsplit(aa_lnc_list[[i]], "___"), "[[", 3))
    my_model_type[[i]] <- unlist(lapply(strsplit(aa_lnc_list[[i]], "___"), "[[", 2))
    print(length(my_model_names[[i]]))
    print(length(my_model_type[[i]]))
    
  }
  
  if(length(model_name) > 0){
    for(i in 1:length(aa_lnc_list)){
      aa_lnc_list[[i]] <- aa_lnc_list[[i]][my_model_names[[i]] %in% model_name]
      my_model_type[[i]] <- my_model_type[[i]][my_model_names[[i]] %in% model_name]
      my_model_names[[i]] <- my_model_names[[i]][my_model_names[[i]] %in% model_name]
      
    }
  }
  my_model_type_uniq  <- lapply(my_model_type, unique)
  my_model_name_uniq  <- lapply(my_model_names, unique)
  
  aa_CV_len_lis <- unlist(lapply(my_model_type_uniq, length))
  aa_name_len_lis <- unlist(lapply(my_model_name_uniq, length))
  print("aa_CV_len_lis")
  print(aa_CV_len_lis)
  print("aa_name_len_lis")
  print(aa_name_len_lis)
  
  aa_my_col_num <- max(aa_CV_len_lis)
  aa_my_row_num <- max(aa_name_len_lis)
  AUROC_mat_list <- list()
  AUPRC_mat_list <- list()
  Importance_list <- list()
  Prediction_list <- list()
  performance_list <- list()
  if(train_perf_included){
    performance_list_train <- list()
    AUROC_mat_list_train <- list()
    AUPRC_mat_list_train <- list()
  }

  for(cur_lnc in 1:length(aa_lnc_list)){
    AUROC_mat_list[[cur_lnc]] <- matrix(nrow = aa_my_row_num, 
                                  ncol = aa_my_col_num)
    AUPRC_mat_list[[cur_lnc]] <- matrix(nrow = aa_my_row_num, 
                                  ncol = aa_my_col_num)
    rownames(AUROC_mat_list[[cur_lnc]]) <- sort(my_model_name_uniq[[which.max(aa_name_len_lis)]])
    colnames(AUROC_mat_list[[cur_lnc]]) <- sort(my_model_type_uniq[[which.max(aa_CV_len_lis)]])
    rownames(AUPRC_mat_list[[cur_lnc]]) <-  rownames(AUROC_mat_list[[cur_lnc]]) 
    colnames(AUPRC_mat_list[[cur_lnc]]) <-  colnames(AUROC_mat_list[[cur_lnc]]) 
    if(train_perf_included){
      AUROC_mat_list_train[[cur_lnc]] <- matrix(nrow = aa_my_row_num, 
                                          ncol = aa_my_col_num)
      AUPRC_mat_list_train[[cur_lnc]] <- matrix(nrow = aa_my_row_num, 
                                          ncol = aa_my_col_num)
      rownames(AUROC_mat_list_train[[cur_lnc]]) <- sort(my_model_name_uniq[[which.max(aa_name_len_lis)]])
      colnames(AUROC_mat_list_train[[cur_lnc]]) <- sort(my_model_type_uniq[[which.max(aa_CV_len_lis)]])
      rownames(AUPRC_mat_list_train[[cur_lnc]]) <-  rownames(AUROC_mat_list_train[[cur_lnc]]) 
      colnames(AUPRC_mat_list_train[[cur_lnc]]) <-  colnames(AUROC_mat_list_train[[cur_lnc]]) 
      performance_list_train[[cur_lnc]] <- list()
    }
    
    Importance_list[[cur_lnc]] <- list()
    Prediction_list[[cur_lnc]] <- list()
    performance_list[[cur_lnc]] <- list()
    my_pt_model_list <- list()
    my_pt_model_name_list <- list()
    aa_allpart <- sort(unique(my_model_type[[cur_lnc]]))
    for(cur_pt in 1:length(aa_allpart)){
      my_pt_model_list[[cur_pt]] <- aa_lnc_list[[cur_lnc]][my_model_type[[cur_lnc]] %in% aa_allpart[cur_pt]]
      my_pt_model_name_list[[cur_pt]] <- my_model_names[[cur_lnc]][my_model_type[[cur_lnc]] %in% aa_allpart[cur_pt]]
      Importance_list[[cur_lnc]][[cur_pt]] <- list()
      Prediction_list[[cur_lnc]][[cur_pt]] <- list()
      performance_list[[cur_lnc]][[cur_pt]] <- list()
      if(train_perf_included){
        performance_list_train[[cur_lnc]][[cur_pt]] <- list()
      }
      stopifnot(length(my_pt_model_name_list[[cur_pt]]) == length(my_pt_model_list[[cur_pt]]))
      for(cur_mod in 1:length(my_pt_model_list[[cur_pt]])){
        load(paste0(Learned_model_folder, "/", my_pt_model_list[[cur_pt]][cur_mod]))
        Importance_list[[cur_lnc]][[cur_pt]][[cur_mod]] <- variable_importance
        Prediction_list[[cur_lnc]][[cur_pt]][[cur_mod]] <- test_predictions
        performance_list[[cur_lnc]][[cur_pt]][[cur_mod]] <- test_perf

        AUROC_mat_list[[cur_lnc]][match(my_pt_model_name_list[[cur_pt]][cur_mod], rownames(AUROC_mat_list[[cur_lnc]])),
                                  match(aa_allpart[cur_pt], colnames(AUROC_mat_list[[cur_lnc]]))] <- performance_list[[cur_lnc]][[cur_pt]][[cur_mod]][[3]]
        AUPRC_mat_list[[cur_lnc]][match(my_pt_model_name_list[[cur_pt]][cur_mod], rownames(AUPRC_mat_list[[cur_lnc]])),
                                  match(aa_allpart[cur_pt], colnames(AUPRC_mat_list[[cur_lnc]]))] <- performance_list[[cur_lnc]][[cur_pt]][[cur_mod]][[4]]
        if(train_perf_included){
          performance_list_train[[cur_lnc]][[cur_pt]][[cur_mod]] <- training_perf
          AUROC_mat_list_train[[cur_lnc]][match(my_pt_model_name_list[[cur_pt]][cur_mod], rownames(AUROC_mat_list_train[[cur_lnc]])),
                                    match(aa_allpart[cur_pt], colnames(AUROC_mat_list_train[[cur_lnc]]))] <- performance_list_train[[cur_lnc]][[cur_pt]][[cur_mod]][[3]]
          AUPRC_mat_list_train[[cur_lnc]][match(my_pt_model_name_list[[cur_pt]][cur_mod], rownames(AUPRC_mat_list_train[[cur_lnc]])),
                                    match(aa_allpart[cur_pt], colnames(AUPRC_mat_list_train[[cur_lnc]]))] <- performance_list_train[[cur_lnc]][[cur_pt]][[cur_mod]][[4]]
          
        }
        
      }
      names(Importance_list[[cur_lnc]][[cur_pt]]) <- my_pt_model_name_list[[cur_pt]]
      names(Prediction_list[[cur_lnc]][[cur_pt]]) <- my_pt_model_name_list[[cur_pt]]
      names(performance_list[[cur_lnc]][[cur_pt]]) <- my_pt_model_name_list[[cur_pt]]
      if(train_perf_included){
        names(performance_list_train[[cur_lnc]][[cur_pt]]) <- my_pt_model_name_list[[cur_pt]]
      }
    }
    names(Importance_list[[cur_lnc]]) <- aa_allpart
    names(Prediction_list[[cur_lnc]]) <- aa_allpart
    names(performance_list[[cur_lnc]]) <- aa_allpart
    if(train_perf_included){
      names(performance_list_train[[cur_lnc]]) <- aa_allpart
    }
  }
  my_results <- list(auroc_mat=AUROC_mat_list, auprc_mat=AUPRC_mat_list)
  if(return_imp){
    my_results$importance <- Importance_list
  }
  if(return_pred){
    my_results$predictions <- Prediction_list
  }
  if(return_per_full){
    my_results$performance_full <- performance_list
  }
  if(train_perf_included){
    my_results$auroc_mat_train <- AUROC_mat_list_train
    my_results$auprc_mat_train <- AUPRC_mat_list_train
  }
  return(my_results)
}
############################################################################################################################################

perf_gather_pair <- function(Learned_model_folder, 
                        model_name = character(0), 
                        lncRNA_name = character(0),
                        return_imp = T,
                        return_pred=F,
                        return_per_full=F,
                        train_perf_included=F,
                        file_pattern = "*___perf_pred_varIMP.RData"){
  #  require(ranger)
  #require(RColorBrewer)
  # require(ROCR)
  # require(PRROC)
  # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  # Learned_model_folder : address to folder containing performance results for a particular lncRNA
  aa_files <- list.files(Learned_model_folder, pattern = file_pattern)
  aa_lnc_list <- list()
  #aa_dataset_list <- list()
  if(length(lncRNA_name) > 0){
    
    for(i in 1:length(lncRNA_name)){
      aa_lnc_list[[i]] <- aa_files[grep(lncRNA_name[i], aa_files)]
      #aa_dataset_list[[i]] <- aa_dataset[grep(lncRNA_name, aa_dataset)]
    }
  }else{ 
    aa_lnc_list[[1]] <- aa_files
    aattt <- unlist(lapply(strsplit(unlist(lapply(strsplit(aa_lnc_list[[1]], "___"), "[[", 1)), "_"), "[[", 5))
    lncRNA_name <- unique(aattt)
    print(lncRNA_name)
    stopifnot(length(lncRNA_name) == 1)
    # aa_dataset_list[[1]] <- aa_dataset[grep(lncRNA_name, aa_dataset)]
    # print("aa_dataset")
    # print(aa_dataset_list[[1]])
  }
  
  my_model_names <- list()
  my_model_type <- list() # RCV1:5, CCV1:5
  for(i in 1:length(aa_lnc_list)){
    #my_model_names[[i]] <- unlist(lapply(strsplit(aa_lnc_list[[i]], "___"), "[[", 3))
    aaxxxc <- unlist(lapply(strsplit(aa_lnc_list[[i]], "_feat_"), "[[", 2))
    my_model_names[[i]] <- unlist(lapply(strsplit(aaxxxc, "___perf_pred"), "[[", 1))
    my_model_type[[i]] <- unlist(lapply(strsplit(aa_lnc_list[[i]], "___"), "[[", 2))
    print(length(my_model_names[[i]]))
    print(length(my_model_type[[i]]))
    
  }
  
  if(length(model_name) > 0){
    for(i in 1:length(aa_lnc_list)){
      aa_lnc_list[[i]] <- aa_lnc_list[[i]][my_model_names[[i]] %in% model_name]
      my_model_type[[i]] <- my_model_type[[i]][my_model_names[[i]] %in% model_name]
      my_model_names[[i]] <- my_model_names[[i]][my_model_names[[i]] %in% model_name]
      
    }
  }
  my_model_type_uniq  <- lapply(my_model_type, unique)
  my_model_name_uniq  <- lapply(my_model_names, unique)
  
  aa_CV_len_lis <- unlist(lapply(my_model_type_uniq, length))
  aa_name_len_lis <- unlist(lapply(my_model_name_uniq, length))
  print("aa_CV_len_lis")
  print(aa_CV_len_lis)
  print("aa_name_len_lis")
  print(aa_name_len_lis)
  
  aa_my_col_num <- max(aa_CV_len_lis)
  aa_my_row_num <- max(aa_name_len_lis)
  AUROC_mat_list <- list()
  AUPRC_mat_list <- list()
  Importance_list <- list()
  Prediction_list <- list()
  performance_list <- list()
  if(train_perf_included){
    performance_list_train <- list()
    AUROC_mat_list_train <- list()
    AUPRC_mat_list_train <- list()
  }
  
  for(cur_lnc in 1:length(aa_lnc_list)){
    AUROC_mat_list[[cur_lnc]] <- matrix(nrow = aa_my_row_num, 
                                        ncol = aa_my_col_num)
    AUPRC_mat_list[[cur_lnc]] <- matrix(nrow = aa_my_row_num, 
                                        ncol = aa_my_col_num)
    rownames(AUROC_mat_list[[cur_lnc]]) <- sort(my_model_name_uniq[[which.max(aa_name_len_lis)]])
    colnames(AUROC_mat_list[[cur_lnc]]) <- sort(my_model_type_uniq[[which.max(aa_CV_len_lis)]])
    rownames(AUPRC_mat_list[[cur_lnc]]) <-  rownames(AUROC_mat_list[[cur_lnc]]) 
    colnames(AUPRC_mat_list[[cur_lnc]]) <-  colnames(AUROC_mat_list[[cur_lnc]]) 
    if(train_perf_included){
      AUROC_mat_list_train[[cur_lnc]] <- matrix(nrow = aa_my_row_num, 
                                                ncol = aa_my_col_num)
      AUPRC_mat_list_train[[cur_lnc]] <- matrix(nrow = aa_my_row_num, 
                                                ncol = aa_my_col_num)
      rownames(AUROC_mat_list_train[[cur_lnc]]) <- sort(my_model_name_uniq[[which.max(aa_name_len_lis)]])
      colnames(AUROC_mat_list_train[[cur_lnc]]) <- sort(my_model_type_uniq[[which.max(aa_CV_len_lis)]])
      rownames(AUPRC_mat_list_train[[cur_lnc]]) <-  rownames(AUROC_mat_list_train[[cur_lnc]]) 
      colnames(AUPRC_mat_list_train[[cur_lnc]]) <-  colnames(AUROC_mat_list_train[[cur_lnc]]) 
      performance_list_train[[cur_lnc]] <- list()
    }
    
    Importance_list[[cur_lnc]] <- list()
    Prediction_list[[cur_lnc]] <- list()
    performance_list[[cur_lnc]] <- list()
    my_pt_model_list <- list()
    my_pt_model_name_list <- list()
    aa_allpart <- sort(unique(my_model_type[[cur_lnc]]))
    for(cur_pt in 1:length(aa_allpart)){
      my_pt_model_list[[cur_pt]] <- aa_lnc_list[[cur_lnc]][my_model_type[[cur_lnc]] %in% aa_allpart[cur_pt]]
      my_pt_model_name_list[[cur_pt]] <- my_model_names[[cur_lnc]][my_model_type[[cur_lnc]] %in% aa_allpart[cur_pt]]
      Importance_list[[cur_lnc]][[cur_pt]] <- list()
      Prediction_list[[cur_lnc]][[cur_pt]] <- list()
      performance_list[[cur_lnc]][[cur_pt]] <- list()
      if(train_perf_included){
        performance_list_train[[cur_lnc]][[cur_pt]] <- list()
      }
      stopifnot(length(my_pt_model_name_list[[cur_pt]]) == length(my_pt_model_list[[cur_pt]]))
      for(cur_mod in 1:length(my_pt_model_list[[cur_pt]])){
        load(paste0(Learned_model_folder, "/", my_pt_model_list[[cur_pt]][cur_mod]))
        Importance_list[[cur_lnc]][[cur_pt]][[cur_mod]] <- variable_importance
        Prediction_list[[cur_lnc]][[cur_pt]][[cur_mod]] <- test_predictions
        performance_list[[cur_lnc]][[cur_pt]][[cur_mod]] <- test_perf
        
        AUROC_mat_list[[cur_lnc]][match(my_pt_model_name_list[[cur_pt]][cur_mod], rownames(AUROC_mat_list[[cur_lnc]])),
                                  match(aa_allpart[cur_pt], colnames(AUROC_mat_list[[cur_lnc]]))] <- performance_list[[cur_lnc]][[cur_pt]][[cur_mod]][[3]]
        AUPRC_mat_list[[cur_lnc]][match(my_pt_model_name_list[[cur_pt]][cur_mod], rownames(AUPRC_mat_list[[cur_lnc]])),
                                  match(aa_allpart[cur_pt], colnames(AUPRC_mat_list[[cur_lnc]]))] <- performance_list[[cur_lnc]][[cur_pt]][[cur_mod]][[4]]
        if(train_perf_included){
          performance_list_train[[cur_lnc]][[cur_pt]][[cur_mod]] <- training_perf
          AUROC_mat_list_train[[cur_lnc]][match(my_pt_model_name_list[[cur_pt]][cur_mod], rownames(AUROC_mat_list_train[[cur_lnc]])),
                                          match(aa_allpart[cur_pt], colnames(AUROC_mat_list_train[[cur_lnc]]))] <- performance_list_train[[cur_lnc]][[cur_pt]][[cur_mod]][[3]]
          AUPRC_mat_list_train[[cur_lnc]][match(my_pt_model_name_list[[cur_pt]][cur_mod], rownames(AUPRC_mat_list_train[[cur_lnc]])),
                                          match(aa_allpart[cur_pt], colnames(AUPRC_mat_list_train[[cur_lnc]]))] <- performance_list_train[[cur_lnc]][[cur_pt]][[cur_mod]][[4]]
          
        }
        
      }
      names(Importance_list[[cur_lnc]][[cur_pt]]) <- my_pt_model_name_list[[cur_pt]]
      names(Prediction_list[[cur_lnc]][[cur_pt]]) <- my_pt_model_name_list[[cur_pt]]
      names(performance_list[[cur_lnc]][[cur_pt]]) <- my_pt_model_name_list[[cur_pt]]
      if(train_perf_included){
        names(performance_list_train[[cur_lnc]][[cur_pt]]) <- my_pt_model_name_list[[cur_pt]]
      }
    }
    names(Importance_list[[cur_lnc]]) <- aa_allpart
    names(Prediction_list[[cur_lnc]]) <- aa_allpart
    names(performance_list[[cur_lnc]]) <- aa_allpart
    if(train_perf_included){
      names(performance_list_train[[cur_lnc]]) <- aa_allpart
    }
  }
  my_results <- list(auroc_mat=AUROC_mat_list, auprc_mat=AUPRC_mat_list)
  if(return_imp){
    my_results$importance <- Importance_list
  }
  if(return_pred){
    my_results$predictions <- Prediction_list
  }
  if(return_per_full){
    my_results$performance_full <- performance_list
  }
  if(train_perf_included){
    my_results$auroc_mat_train <- AUROC_mat_list_train
    my_results$auprc_mat_train <- AUPRC_mat_list_train
  }
  return(my_results)
}
############################################################################################################################################
performance_to_Granges <- function(tile_GR=Partition_6_dfs_GR,
                                   partition_type= c("R"),
                                   partition_dataset,
                                   model_folder,
                                   model_name,
                                   nu_partitions=5){
  # tile_GR is a granges object that contains coordinates for all tiles considered here, with a metadata column called "tile" containing the name of the tile
  # partition_type: is a character vector with length 1, either "C" for chunk or "R" for random partitioning
  # model_folder: address of the folder where all models for this lncRNA reside
  # model_name: a character vector of length 1 or higher. is the name of the models to be visualized 
  # partition_dataset: is a dataframe that has at least the following columns: tile_name,owner, label, RCV1 RCV2 RCV3 RCV4 RCV5 CCV1 CCV2 CCV3 CCV4 CCV5
  stopifnot(all(partition_dataset$tile_name %in% tile_GR$tile))
  nu_models <- length(model_name)
  aa_files <- list.files(model_folder, pattern = "*_perf_pred_varIMP.RData")
  model_file_list <- list()
  model_part_list <- list()
  for(i in 1:nu_models){ # get the name of models to work with
    model_file_list[[i]] <- aa_files[grep(pattern = paste0("___",model_name[i],"___"), x = aa_files)]
    model_file_list[[i]] <- model_file_list[[i]][grep(pattern = paste0("___",partition_type,"CV"), x = model_file_list[[i]])]
    model_part_list[[i]] <- unlist(lapply(strsplit(model_file_list[[i]], "___"), "[[", 2))
  }
  names(model_file_list) <- model_name
  
  aa_all_part <- paste0(partition_type,"CV", c(1:nu_partitions))
  # load files one by one and form the full prediction vector
  my_prediction_list_byModel_by_part <- list()
  my_prediction_list_byModel <- list()
  my_prediction_list_byModelGR <- list()
  for(cur_model in 1:nu_models){
    my_prediction_list_byModel_by_part[[cur_model]] <- list()
    for(cur_part in 1:length(model_file_list[[cur_model]])){
      load(paste0(model_folder, "/", model_file_list[[cur_model]][cur_part]))
      cur_prediction <-  test_predictions$predictions[, 2]
      aawtest <- which(partition_dataset[,match(model_part_list[[cur_model]][cur_part], colnames(partition_dataset))] == 0)
      cur_name <- partition_dataset$tile_name[aawtest]
      names(cur_prediction) <- cur_name
      my_prediction_list_byModel_by_part[[cur_model]][[cur_part]] <- cur_prediction
    }
    my_prediction_list_byModel[[cur_model]] <- do.call(c, my_prediction_list_byModel_by_part[[cur_model]])
    my_prediction_list_byModel[[cur_model]] <- my_prediction_list_byModel[[cur_model]][match(partition_dataset$tile_name, names(my_prediction_list_byModel[[cur_model]]))]
    stopifnot(((length(model_part_list[[cur_model]]) < nu_partitions) | (sum(is.na(my_prediction_list_byModel[[cur_model]])) == 0)))
    if(length(model_part_list[[cur_model]]) < nu_partitions){
      print(paste0("Predictions missing for model ", model_name, " partitions ", setdiff(aa_all_part, model_part_list[[cur_model]])))
    }
   # my_prediction_list_byModelGR[[cur_model]] <- tile_GR[match(,tile_GR$tile)]
  }
  # Put each model prediction in one Granges
  my_init_GR <- tile_GR[match(partition_dataset$tile_name,tile_GR$tile)]
  my_out_Granges_list <- list()
  for(i in 1:length(my_prediction_list_byModel)){
    my_out_Granges_list[[i]] <- my_init_GR
    my_out_Granges_list[[i]]$model_score <-  my_prediction_list_byModel[[i]]
    my_out_Granges_list[[i]]$label <- partition_dataset$label
  }
  names(my_out_Granges_list) <- model_name
  return(my_out_Granges_list)
}
################################################################################################################################################################################################
# example
# save(list = c('Partition_6_dfs_GR'), file = "Partition_6_dfs_GR.RData")
# 
# load('input_RData/Partition_6_CVfirst_dataset_Malat1.RData')
# aatst <- performance_to_Granges(tile_GR=Partition_6_dfs_GR,
#                                 partition_type= c("R"),
#                                 partition_dataset = my_partition,
#                                 model_folder = "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_svd2_lowres_Results/Learned_models/Malat1",
#                                 model_name = c( "dist", "kmerFq", "motSc_pairs_triplx","sequ_dist_trnsp", "sequ_dist_Chrom_Meth_AX_ChIP_trnsp"))
################################################################################################################################################################################################
# get Granges and make Gviz figure

# can visualize features as well --> like transcription for Malat1
# Add each model predictions as a datatrack
# add a track for labels
# add track for feature vectors if provided
Gviz_modelPred <- function(prediction_Granges_list,
                           feature_Granges_list = list(),
                           partition_df,
                           nu_partitions=5,
                           partition_type= c("R"),
                           my_genome = "mm9"){
  # prediction_Granges_list is the output of performance_to_Granges function
  #feature_Granges_list is a list similar to prediction_Granges_list but with features as input
  #partition_df is a dataframe that has at least the following columns: tile_name,owner, label, RCV1 RCV2 RCV3 RCV4 RCV5 CCV1 CCV2 CCV3 CCV4 CCV5
  # partition_type: is a character vector with length 1, either "C" for chunk or "R" for random partitioning
  
  require(Gviz)
  require(RColorBrewer)
  nu_models <- length(prediction_Granges_list)
  model_names <- names(prediction_Granges_list)
  # create a factor for partiton of each datapoint
  
  
  aa_all_part <- paste0(partition_type,"CV", c(1:nu_partitions))
  aa_part_as <- character(length = nrow(partition_df))
  for(i in 1:length(aa_all_part)){
    aa_part_as[partition_df[,match(aa_all_part[i], colnames(partition_df))] == 0] <- aa_all_part[i]
  }
  aa_part_as2 <- factor(aa_part_as, levels = aa_all_part)
  print(table(aa_part_as2))
  print(length(aa_part_as2))
  #aa_part_as3 <- factor(aa_all_part)
  #prediction_Granges_list[[i]]$my_grouping <- aa_part_as2
  
  qual_col_pals = brewer.pal.info[rownames(brewer.pal.info) == 'Set3',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  #prediction_Granges_list[[i]]$my_color <- col_vector[aa_part_as2]
  
  agtrack <- GenomeAxisTrack(cex = 0.55)
  aitrack <- IdeogramTrack(genome = my_genome)
  
  prediction_track_list <- list()
  for(i in 1:nu_models){
    #print(length(prediction_Granges_list[[i]]))
    # add each partition as a sample
    #tmp_plt_list <- list()
    #print("before")
    tmp_dataset <- matrix(0L, nrow = (length(aa_all_part) + 4), 
                          ncol = length(prediction_Granges_list[[i]]))
    #print("after")
    rownames(tmp_dataset) <- c(aa_all_part, "TP", "FP", "TN", "FN")
    
    for(j in 1:length(aa_all_part)){
      #print(j)
      aa_cur_ind <- which(partition_df[,match(aa_all_part[j], colnames(partition_df))] == 0)
      tmp_dataset[j, aa_cur_ind] <- prediction_Granges_list[[i]]$model_score[aa_cur_ind]
      aamycurp <- prediction_Granges_list[[i]]$model_score[aa_cur_ind]
      aamycurl <- prediction_Granges_list[[i]]$label[aa_cur_ind]
      aamycurp_quan <- quantile(aamycurp, seq(0,1,0.1))
      aawwww <- which((aamycurp > aamycurp_quan[6]) & (aamycurl >= 1))
      #print("length(aawwww) TP")
      #print(length(aawwww))
      
      if(length(aawwww) > 0){
        tmp_dataset[length(aa_all_part) + 1, aa_cur_ind[aawwww]] <- aamycurp[aawwww]
        tmp_dataset[j                      , aa_cur_ind[aawwww]] <- rep(0, length(aawwww))
      }

      aawwww <- which((aamycurp > aamycurp_quan[6]) & (aamycurl < 1))
      #print("length(aawwww) FP")
      #print(length(aawwww))
      
      if(length(aawwww) > 0){
        tmp_dataset[length(aa_all_part) + 2, aa_cur_ind[aawwww]] <- aamycurp[aawwww]
        tmp_dataset[j                      , aa_cur_ind[aawwww]] <- rep(0, length(aawwww))
      }
      
      aawwww <- which((aamycurp <= aamycurp_quan[6]) & (aamycurl < 1))
      #print("length(aawwww) TN")
      #print(length(aawwww))
      
      if(length(aawwww) > 0){
        tmp_dataset[length(aa_all_part) + 3, aa_cur_ind[aawwww]] <- aamycurp[aawwww]
        tmp_dataset[j                      , aa_cur_ind[aawwww]] <- rep(0, length(aawwww))
      }
      
      aawwww <- which((aamycurp <= aamycurp_quan[6]) & (aamycurl >= 1))
      #print("length(aawwww) FN")
      #print(length(aawwww))
      
      if(length(aawwww) > 0){
        tmp_dataset[length(aa_all_part) + 4, aa_cur_ind[aawwww]] <- aamycurp[aawwww]
        tmp_dataset[j                      , aa_cur_ind[aawwww]] <- rep(0, length(aawwww))
      }
      
    }
    
    tmp_dataset <- tmp_dataset[((nrow(tmp_dataset)-3):(nrow(tmp_dataset))),]
    prediction_track_list[[i]] <-DataTrack(range = prediction_Granges_list[[i]],
                                           # start = start(prediction_Granges_list[[i]]),
                                           #                              end = end(prediction_Granges_list[[i]]),
                                           #                              chromosome = seqnames(prediction_Granges_list[[i]]),
                                           data = tmp_dataset,
                                           genome = my_genome, 
                                           name = model_names[i],  
                                           window=100,aggregation="mean",
                                           groups = factor(rownames(tmp_dataset)),
                                           col = c(col_vector[4], col_vector[3], col_vector[2], col_vector[1]),#col_vector[c(1:(nu_partitions+1))+3],
                                           legend = TRUE, 
                                           ylim = c(-0.02, 1.02),
                                           cex.legend = 1, 
                                           type= "hist") #, degree = 1, span = 0.66, family = "gaussian"
    
    displayPars(prediction_track_list[[i]]) <- list(groups = factor(rownames(tmp_dataset)), legend = TRUE)
    
  }
  prediction_Granges_list[[1]]$label <- prediction_Granges_list[[1]]$label + 0.2
  tmp_dataset <- matrix(0L, nrow = length(aa_all_part), ncol = length(prediction_Granges_list[[i]]))
  rownames(tmp_dataset) <- aa_all_part
  for(j in 1:length(aa_all_part)){
    aa_cur_ind <- which(partition_df[,match(aa_all_part[j], colnames(partition_df))] == 0)
    tmp_dataset[j, aa_cur_ind] <- prediction_Granges_list[[1]]$label[aa_cur_ind]
  }
  label_track <- DataTrack(range = prediction_Granges_list[[1]],
                           data = tmp_dataset,
                           genome = my_genome, 
                           name = "Label",  
                           groups = factor(aa_all_part),
                           col = col_vector[c(1:(nu_partitions))+5],
                           window = 100,aggregation="mean",
                           legend = TRUE, 
                           ylim = c(-0.00, 1.22),
                          cex.legend =1, 
                           type= "hist")
  if(length(feature_Granges_list) > 0){
    #stopifnot(length(feature_Granges_list[[1]]) == nrow(partition_df))
    feature_track_list <- list()
    for(cur_ff in 1:length(feature_Granges_list)){
      #print(i)
      feature_track_list[[cur_ff]] <- DataTrack(range = feature_Granges_list[[cur_ff]],
                                              data = "score",
                                              genome = my_genome, 
                                              name = names(feature_Granges_list)[cur_ff], 
                                              window = 100,aggregation="mean",
                                              #groups = aa_part_as2,
                                              #col = col_vector[1:nu_partitions],
                                              # legend = TRUE, 
                                              #ylim = c(-0.02, 1.02),
                                              #cex.legend = 0.7, 
                                              type= "hist")
    }
    my_res <- list(aitrack, agtrack, label_track, unlist(prediction_track_list),  unlist(feature_track_list))
  }else{
    my_res <- list(aitrack, agtrack, label_track, unlist(prediction_track_list))
  }
  
  return(my_res)
}
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # #v# # # # #
plot_gviz_labelOnly <- function(tile_GR,
                                my_partition_dataset,
                                lnc_name,
                                nu_partitions=5, 
                                my_genome="mm9", 
                                output_address){
  require(Gviz)
  require(RColorBrewer)
  require(colorspace)
  # to compare chunk and random CV for a particular lncRNA
  my_tiles <- my_partition_dataset$tile_name[my_partition_dataset$owner %in% lnc_name]
  my_part <- my_partition_dataset[my_partition_dataset$owner %in% lnc_name,]
  my_gr <- tile_GR[tile_GR$tile %in% my_tiles]
  my_gr <- my_gr[match(my_part$tile_name, my_gr$tile)]
  # create datatracks for labels of each partitioning type:Chunk, Random
  my_rand_parts <- paste0("RCV", c(1:nu_partitions))
  my_chun_parts <- paste0("CCV", c(1:nu_partitions))
  my_part$label <- my_part$label + 0.2
  print("before making colors")
  qual_col_pals = brewer.pal.info[rownames(brewer.pal.info) == 'Set3',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- colorspace::darken(col = col_vector,amount = 0.3)
  print("before making tmp dataset")
  tmp_dataset <- matrix(0L, nrow = length(my_rand_parts), ncol = length(my_gr))
  rownames(tmp_dataset) <- my_rand_parts
  for(j in 1:length(my_rand_parts)){
    aa_cur_ind <- which(my_part[,match(my_rand_parts[j], colnames(my_part))] == 0)
    tmp_dataset[j, aa_cur_ind] <- my_part$label[aa_cur_ind]
  }
  print("before making rand track")
  label_track_rand <- DataTrack(range = my_gr,
                           data = tmp_dataset,
                           genome = my_genome, 
                           name = "Random CV",  
                           groups = factor(my_rand_parts),
                           col = col_vector[c(1:(nu_partitions))+5],
                           #window = 1000,aggregation="mean",
                           legend = TRUE, 
                           ylim = c(-0.00, 1.22),
                           cex.legend = 1, 
                           type= "histogram",
                           col.histogram="#FFFFFF")
  
  tmp_dataset <- matrix(0L, nrow = length(my_chun_parts), ncol = length(my_gr))
  rownames(tmp_dataset) <- my_chun_parts
  for(j in 1:length(my_chun_parts)){
    aa_cur_ind <- which(my_part[,match(my_chun_parts[j], colnames(my_part))] == 0)
    tmp_dataset[j, aa_cur_ind] <- my_part$label[aa_cur_ind]
  }
  print("before making chunk track")
  label_track_chun <- DataTrack(range = my_gr,
                                data = tmp_dataset,
                                genome = my_genome, 
                                name = "Chunk CV",  
                                groups = factor(my_chun_parts),
                                col = col_vector[c(1:(nu_partitions))+5],
                                #window = 1000,aggregation="mean",
                                legend = TRUE, 
                                ylim = c(-0.00, 1.22),
                                cex.legend = 1, 
                                type= "histogram",
                                col.histogram="#FFFFFF")
  my_label_tracks <- list(label_track_rand, label_track_chun)
  
  
  if(is.factor(seqnames(my_gr))){
    my_chr <- unique(levels(seqnames(my_gr))[as.numeric(seqnames(my_gr))])
  }else{
    my_chr <- unique(seqnames(my_gr))
  }
  if(is.factor(my_chr)){
    my_chr <- as.character(levels(my_chr)[as.numeric(my_chr)])
  }
  print("my_chr")
  print(my_chr)
  aa_min_max <- data.frame(Chr = my_chr, Start = numeric(length(my_chr)), end = numeric(length(my_chr)))
  for(cur_chr in 1:length(my_chr)){
    aa_min_max[cur_chr, 2] <- min(start(my_gr[seqnames(my_gr) %in% my_chr[cur_chr]]))
    aa_min_max[cur_chr, 3] <- max(end(my_gr[seqnames(my_gr) %in% my_chr[cur_chr]]))
  }
  print("aa_min_max")
  print(aa_min_max)
  my_coord_grange <- makeGRangesFromDataFrame(aa_min_max)
  
  
  # my_coord_grange_wind <- my_coord_grange
  # my_chr_new <- my_chr
  
  my_coord_grange_wind <- unlist(slidingWindows(my_coord_grange, width=5000000, step=5000000))
  print("my_coord_grange_wind")
  print(my_coord_grange_wind)
  my_chr_new <- seqnames(my_coord_grange_wind)
  #if(is.factor(my_chr_new)){
  my_chr_new <- as.character(levels(my_chr_new)[as.numeric(my_chr_new)])
  
  #output_address
  my_file_names_coord <- paste(seqnames(my_coord_grange_wind), start(my_coord_grange_wind), end(my_coord_grange_wind), sep = "_")
  #my_file_names_model <- paste0(my_model_name, collapse = ".")
  my_file_names_lnc <- lnc_name
  file_pre_name <- "Rand_vs_Chunk"
  my_file_names <- paste0(file_pre_name,"__",my_file_names_lnc, "__", my_file_names_coord, ".png")
  
  print("my_file_names")
  print(my_file_names)
  
  for(cur_coor in 1:length(my_coord_grange_wind)){
    #aa_cleen <- ceiling((aa_min_max[cur_coor, 3] - aa_min_max[cur_coor, 2])/1000)
    #aa_my_wd <- max(20, ceiling(aa_cleen/600))
    # print(aa_my_wd)

      aamych <- unlist(strsplit(my_chr_new[cur_coor], split = "chr"))
      #print(aamych)
      aa_trrlist <- unlist(my_label_tracks)
      png(filename = paste0(output_address, "/", my_file_names[cur_coor]),
          width = 10*300,
          height = 2 * 300,# max(10, length(aa_trrlist)*1.5)*300,
          res = 300,
          #units = "in",
          pointsize = 11)
      #print(aa_min_max[cur_coor, 1])
      
      aitrack <- IdeogramTrack(genome = my_genome, chromosome = aamych[length(aamych)])
      aa_trrlist2 <- list()
      aa_trrlist2[[1]] <- aitrack
      aa_trrlist2[2:(length(aa_trrlist)+1)] <- aa_trrlist
      plotTracks(trackList=aa_trrlist2,
                 from = start(my_coord_grange_wind)[cur_coor] - 1000, # aa_min_max[cur_coor, 2] - 1000,#start(my_coord_grange)[cur_coor],
                 to = end(my_coord_grange_wind)[cur_coor] + 1000, #aa_min_max[cur_coor, 3] + 1000,#end(my_coord_grange)[cur_coor],
                 chromosome = aamych[length(aamych)]#aa_min_max[cur_coor, 1]#seqnames(my_coord_grange)[cur_coor]
                 # ylim = aaylims
                 # type = "hist", legend = T
                 #, window = 30,
                 #,background.panel = "#FFFEDB",
                 ,background.title = "darkblue"
      )
      dev.off()
    
    
  }
  print("completed plotting ...")
}
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # #v# # # # #

# # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # #
# # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # ## # # # #v# # # # #

plot_gviz_lncRNA <- function(my_tile_GR=Partition_6_dfs_GR,
                             my_partition_type= c("R"),
                             my_partition_dataset,
                             my_model_folder,
                             my_model_name,
                             my_nu_partitions=5,
                             my_feature_Granges_list= list(),
                             my_Gnome = "mm9",
                             output_address,
                             file_pre_name = "Gvisual",
                             my_coord_grange = numeric(0),
                             lnc_name,
                             sliding_wind_length=5000000){
  require(Gviz)
  # my_coord_grange is a granges object, each rntey of which will be plotted separately
  # sliding_wind_length : is the lenth of sliding window being used, if set as 0, it will not use sliding windows and plots the originl thing
  prediction_gr <- performance_to_Granges(tile_GR=my_tile_GR,
                                  partition_type= my_partition_type,
                                  partition_dataset = my_partition_dataset,
                                  model_folder = my_model_folder,
                                  model_name = my_model_name,
                                  nu_partitions = my_nu_partitions)
  #print(prediction_gr)
  gviz_plot <- Gviz_modelPred(prediction_Granges_list = prediction_gr,
                         feature_Granges_list = my_feature_Granges_list,
                         partition_df = my_partition_dataset,
                         nu_partitions=my_nu_partitions,
                         partition_type= my_partition_type,
                         my_genome = my_Gnome)
  # # getting the coordinates to plot id not provided
  if(length(my_coord_grange) == 0){
    if(is.factor(seqnames(prediction_gr[[1]]))){
      my_chr <- unique(levels(seqnames(prediction_gr[[1]]))[as.numeric(seqnames(prediction_gr[[1]]))])
    }else{
      my_chr <- unique(seqnames(prediction_gr[[1]]))
    }
    if(is.factor(my_chr)){
      my_chr <- as.character(levels(my_chr)[as.numeric(my_chr)])
    }
    print("my_chr")
    print(my_chr)
    aa_min_max <- data.frame(Chr = my_chr, Start = numeric(length(my_chr)), end = numeric(length(my_chr)))
    for(cur_chr in 1:length(my_chr)){
      aa_min_max[cur_chr, 2] <- min(start(prediction_gr[[1]][seqnames(prediction_gr[[1]]) %in% my_chr[cur_chr]]))
      aa_min_max[cur_chr, 3] <- max(end(prediction_gr[[1]][seqnames(prediction_gr[[1]]) %in% my_chr[cur_chr]]))
    }
    print("aa_min_max")
    print(aa_min_max)
    my_coord_grange <- makeGRangesFromDataFrame(aa_min_max)
  }
 print("my_coord_grange")
 print(my_coord_grange)
 if(sliding_wind_length != 0){
   my_coord_grange_wind <- unlist(slidingWindows(my_coord_grange, width=sliding_wind_length, step=sliding_wind_length/2))
   print("my_coord_grange_wind")
   print(my_coord_grange_wind)
   my_chr_new <- seqnames(my_coord_grange_wind)
   #if(is.factor(my_chr_new)){
   my_chr_new <- as.character(levels(my_chr_new)[as.numeric(my_chr_new)])
   #}
   print(my_chr_new)
 }else{
   my_coord_grange_wind <- my_coord_grange
   my_chr_new <- seqnames(my_coord_grange_wind)
   my_chr_new <- as.character(levels(my_chr_new)[as.numeric(my_chr_new)])
   #my_chr_new <- my_chr
 }

  #output_address
  my_file_names_coord <- paste(seqnames(my_coord_grange_wind), start(my_coord_grange_wind), end(my_coord_grange_wind), sep = "_")
  my_file_names_model <- paste0(my_model_name, collapse = ".")
  my_file_names_lnc <- lnc_name
  if(length(my_coord_grange_wind$auprc) > 0){
    aup <- paste(format(my_coord_grange_wind$auroc, digits = 3),format(my_coord_grange_wind$auprc, digits = 3), sep="_" )
    my_file_names_coord <- paste(my_file_names_coord, aup, sep = "_")
  }
  my_file_names <- paste0(file_pre_name,"__",my_file_names_lnc,"__", my_partition_type, "__", my_file_names_model, "__", my_file_names_coord, ".png")
  
  print("my_file_names")
  print(my_file_names)

  for(cur_coor in 1:length(my_coord_grange_wind)){
    #aa_cleen <- ceiling((aa_min_max[cur_coor, 3] - aa_min_max[cur_coor, 2])/1000)
    #aa_my_wd <- max(20, ceiling(aa_cleen/600))
   # print(aa_my_wd)
    aaov <- findOverlaps(query = my_coord_grange_wind[cur_coor], subject = prediction_gr[[1]])
    if(length(aaov@from) > 0){
      aamych <- unlist(strsplit(my_chr_new[cur_coor], split = "chr"))
      #print(aamych)
      aa_trrlist <- unlist(gviz_plot)
      png(filename = paste0(output_address, "/", my_file_names[cur_coor]),
          width = 20*300,
          height = max(10, length(aa_trrlist)*1.5)*300,
          res = 300,
          #units = "in",
          pointsize = 11)
      #print(aa_min_max[cur_coor, 1])

      aitrack <- IdeogramTrack(genome = my_Gnome, chromosome = aamych[length(aamych)])
      aa_trrlist[[1]] <- aitrack
      plotTracks(trackList=aa_trrlist,
                 from = start(my_coord_grange_wind)[cur_coor] - 1000, # aa_min_max[cur_coor, 2] - 1000,#start(my_coord_grange)[cur_coor],
                 to = end(my_coord_grange_wind)[cur_coor] + 1000, #aa_min_max[cur_coor, 3] + 1000,#end(my_coord_grange)[cur_coor],
                 chromosome = aamych[length(aamych)]#aa_min_max[cur_coor, 1]#seqnames(my_coord_grange)[cur_coor]
                 # ylim = aaylims
                 # type = "hist", legend = T
                 #, window = 30,
                 #,background.panel = "#FFFEDB",
                 ,background.title = "darkblue"
      )
      dev.off()
    }

  }
  print("completed plotting ...")
  
}



# example
# my_lnc_name <- "Gm14820"
# source("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/RF_evaluation_functions.R")
# load("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/distance_feature_partition_6_lowRes_dfGRs.RData")
# aanmm <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_CVfirst_dataset_", my_lnc_name, ".RData")
# load(aanmm)
# load("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_random_chunk_cv_df_updated.RData")
# load("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_dfs_GR.RData")
# load("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/distance_feature_partition_6_lowRes_df.RData")
# 
# aaw <- which(Partition_6_random_chunk_cv_df$owner %in% my_lnc_name)
# my_partition <- Partition_6_random_chunk_cv_df[aaw,]
# distance_feature_partition_6_lowRes_dfGR$score <- distance_feature_partition_6_lowRes_dfGR$lowerres1Mg
# 
# aammmfl <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_lowres1MG/Learned_models/", my_lnc_name)
# aamy_feature_Granges_list <- list()
# my_distance_feature_partition_6_lowRes_df <- distance_feature_partition_6_lowRes_dfGR[distance_feature_partition_6_lowRes_df$owner == my_lnc_name]
# #my_new_dist <- my_distance_feature_partition_6_lowRes_df$lowerres1Mg[match(aanammmsp, my_distance_feature_partition_6_lowRes_df$tile_name)]
# aamy_feature_Granges_list[[1]] <- my_distance_feature_partition_6_lowRes_df
# names(aamy_feature_Granges_list) <- "Genomic_distance"
# plot_gviz_lncRNA(my_tile_GR=Partition_6_dfs_GR,
#                  my_partition_type= c("R"),
#                  my_partition_dataset = my_partition,
#                  my_model_folder = aammmfl,
#                  my_model_name = c("dist", "kmerFq", "motSc_pairs_triplx", "sequ_dist", "sequ_dist_trnsp","sequ_dist_Chrom_Meth_AX_ChIP_trnsp"),
#                  my_nu_partitions=5,
#                  my_feature_Granges_list= aamy_feature_Granges_list,
#                  my_Gnome = "mm9",
#                  output_address = "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_lowres1MG/Gviz",
#                  file_pre_name = "Gvisual",
#                  my_coord_grange = numeric(0),
#                  lnc_name = my_lnc_name)



# source("RF_evaluation_functions.R")
# 
# aatst <- performance_to_Granges(tile_GR=Partition_6_dfs_GR,
#                                 partition_type= c("C"),
#                                 partition_dataset = my_partition,
#                                 model_folder = "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_svd2_lowres_Results/Learned_models/Gm14820",
#                                 model_name = c( "dist", "kmerFq", "motSc_pairs_triplx","sequ_dist_trnsp", "sequ_dist_Chrom_Meth_AX_ChIP_trnsp"))
# # 
# # 
# aapl <- Gviz_modelPred(prediction_Granges_list = aatst,
#                        feature_Granges_list = list(),
#                        partition_df = my_partition,
#                        nu_partitions=5,
#                        partition_type= c("C"),
#                        my_genome = "mm9")
# # #
# # # # #
# # # # #
# # # # #
# # # # #
# png(filename = "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/tst2gm14.png",
#     width = 20*300,
#     height = 10*300,
#     res = 300,
#     pointsize = 11)
# plotTracks(unlist(aapl),
#            from = 5880000,
#            to = 9015000,
#            chromosome = "chrX"
#            # ylim = aaylims
#            # type = "hist", legend = T
#            #, window = 30,
#            #,background.panel = "#FFFEDB",
#            ,background.title = "darkblue"
# )
# dev.off()
# 
# 
# 
# 
# 
# 
# aagr <- Partition_6_dfs_GR[match(Partition_6_chunk_cv[[1]]$tile_name, Partition_6_dfs_GR$tile)]
# aagr$model_score <- rep(0.5, length(aagr))
# aagr$label <- Partition_6_chunk_cv[[1]]$label
# aapp <- Partition_6_chunk_cv[[1]]
# colnames(aapp)[c(4:8)] <- paste0("C",colnames(aapp)[c(4:8)])
# aa33x <- list(rand_mod = aagr, randmod2 = aagr)
# 
# aapl <- Gviz_modelPred(prediction_Granges_list = aa33x,
#                        feature_Granges_list = list(),
#                        partition_df = aapp,
#                        nu_partitions=5,
#                        partition_type= c("C"),
#                        my_genome = "mm9")
# 
# plotTracks(unlist(aapl),
#            from = 35880000,
#            to = 75015000,
#            chromosome = "chr19"
#            # ylim = aaylims
#            # type = "hist", legend = T
#            #, window = 30,
#            #,background.panel = "#FFFEDB",
#            ,background.title = "darkblue"
# )

evaluate_byWindow <- function(my_tile_GR,
                              my_partition_dataset,
                              my_model_folder,
                              my_model_name){
  library(ranger, lib.loc = "~/R/library/")
  require(PRROC)
  require(ROCR)
  
  # aaw <- which(Partition_6_random_chunk_cv_df$owner %in% my_lnc_name)
  # my_partition <- Partition_6_random_chunk_cv_df[aaw,]
  
  
  prediction_gr <- performance_to_Granges(tile_GR=my_tile_GR,
                                          partition_type= "C",
                                          partition_dataset = my_partition_dataset,
                                          model_folder = my_model_folder,
                                          model_name = my_model_name,
                                          nu_partitions = 5)
  if(is.factor(seqnames(prediction_gr[[1]]))){
    my_chr <- unique(levels(seqnames(prediction_gr[[1]]))[as.numeric(seqnames(prediction_gr[[1]]))])
  }else{
    my_chr <- unique(seqnames(prediction_gr[[1]]))
  }
  if(is.factor(my_chr)){
    my_chr <- as.character(levels(my_chr)[as.numeric(my_chr)])
  }
  print("my_chr")
  print(my_chr)
  aa_min_max <- data.frame(Chr = my_chr, Start = numeric(length(my_chr)), end = numeric(length(my_chr)))
  for(cur_chr in 1:length(my_chr)){
    aa_min_max[cur_chr, 2] <- min(start(prediction_gr[[1]][seqnames(prediction_gr[[1]]) %in% my_chr[cur_chr]]))
    aa_min_max[cur_chr, 3] <- max(end(prediction_gr[[1]][seqnames(prediction_gr[[1]]) %in% my_chr[cur_chr]]))
  }
  print("aa_min_max")
  print(aa_min_max)
  my_coord_grange <- makeGRangesFromDataFrame(aa_min_max)
  my_coord_grange_wind <- unlist(slidingWindows(my_coord_grange, width=1000000, step=250000))
  print("my_coord_grange_wind")
  print(my_coord_grange_wind)
  my_chr_new <- seqnames(my_coord_grange_wind)
  #if(is.factor(my_chr_new)){
  my_chr_new <- as.character(levels(my_chr_new)[as.numeric(my_chr_new)])
  
  aaov <- findOverlaps(query = my_coord_grange_wind, subject = prediction_gr[[1]])
  aaovdf <- as.data.frame(cbind(aaov@from, aaov@to))
  aaovdf_agg1 <- aggregate(aaovdf["V2"], aaovdf["V1"], c, simplify=F)
  
  aadupd <- duplicated(aaovdf_agg1$V2)
  if(sum(aadupd) > 0){
    aa_to_rem <- aaovdf_agg1$V1[which(aadupd)]
  }else{
    aa_to_rem <- numeric(0)
  }
  aamun <- setdiff(unique(aaov@from), aa_to_rem)
  
  aaovdf_agg2 <- aggregate(aaovdf["V2"], aaovdf["V1"], length, simplify=T)
  aa_to_rem2 <- aaovdf_agg2$V1[aaovdf_agg2$V2 < 100]
  aamun <- setdiff(aamun, aa_to_rem2)
  my_coord_grange_wind_filtered <- my_coord_grange_wind[aamun]
  
  my_auroc_list <- list() 
  my_auprc_list <- list()  
  my_scored_grange <- list()
  # evaluate predictions per window and output a granges object with two metadata columns indicating the AUROC and AUPRC_diff
  for(cur_mod in 1:length(prediction_gr)){
    my_auroc_list[[cur_mod]] <- numeric(length = length(my_coord_grange_wind_filtered))
    my_auprc_list[[cur_mod]] <- numeric(length = length(my_coord_grange_wind_filtered))
    my_scored_grange[[cur_mod]] <- my_coord_grange_wind_filtered
    for(cur_reg in 1:length(my_coord_grange_wind_filtered)){
      aaov <- findOverlaps(query = my_coord_grange_wind_filtered[cur_reg], subject = prediction_gr[[cur_mod]])
      cur_pred <- prediction_gr[[cur_mod]]$model_score[unique(aaov@to)]
      cur_label <- prediction_gr[[cur_mod]]$label[unique(aaov@to)]
      print(table(cur_label))
      if(length(table(cur_label)) != 2){
        my_auroc_list[[cur_mod]][cur_reg] <- NA
        my_auprc_list[[cur_mod]][cur_reg] <- NA
      }else{
        my_prediction <- prediction(cur_pred,cur_label)
        auc <- performance(my_prediction, measure = "auc")
        auc <- auc@y.values[[1]]
        index_Pos <- cur_label == 1
        index_Neg <- cur_label == 0
        base_prc <- sum(cur_label == 1)/length(cur_label)
        aaprc <- pr.curve(cur_pred[index_Pos], cur_pred[index_Neg], curve = F)
        aaprc <- aaprc$auc.integral
        my_auroc_list[[cur_mod]][cur_reg] <- auc
        my_auprc_list[[cur_mod]][cur_reg] <- aaprc-base_prc
      }

    }
    my_scored_grange[[cur_mod]]$auroc <- my_auroc_list[[cur_mod]]
    my_scored_grange[[cur_mod]]$auprc <- my_auprc_list[[cur_mod]]
  }
  names(my_scored_grange) <- my_model_name
  return(my_scored_grange)
  
}


#example
# my_lnc_name <- "Gm14820"
# source("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/RF_evaluation_functions.R")
# load("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_random_chunk_cv_df_updated.RData")
# load("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_dfs_GR.RData")
# aaw <- which(Partition_6_random_chunk_cv_df$owner %in% my_lnc_name)
# my_partition <- Partition_6_random_chunk_cv_df[aaw,]
# aammmfl <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_200feat/Learned_models/", my_lnc_name)
# 
# aattsst <- evaluate_byWindow(my_tile_GR = Partition_6_dfs_GR,
#                              my_partition_dataset = my_partition,
#                              my_model_folder = aammmfl,
#                              my_model_name = c("top200Features_RBPFilter", "top200Features_BothFilter", "top200Features_RepFilter", 
#                                                "top200Features_noFilter"))
