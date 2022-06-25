args <- commandArgs(trailingOnly = TRUE)
# first argument is the address of RData file containing the features and partitionings
# second feature is either "rand", or "chunk"
# third feature is the number to be used for random seed
# fourth and after that are features that want to be removed before training dataset
library(ranger)
source("~/Documents/Shayan/BioInf/lncRNA/RF_evaluation_functions.R")

model_store_address <- "~/Documents/Shayan/BioInf/lncRNA/Learned_models/"
plot_strore_address <- "~/Documents/Shayan/BioInf/lncRNA/performance_plots/"
# 
input_rdata <- args[1]
partition_mode <- args[2] # can be either "rand" or "chunk"
my_random_seed <- as.numeric(args[3])
if(length(args) > 3){
	avoid_columns <- as.numeric(args[4:length(args)])
	}else{
		avoid_columns <- numeric(0)
	}

stopifnot(sum(is.na(avoid_columns)) == 0,
          args[2] %in% c("rand", "chunk"))
print("loading input data ...")
load(input_rdata)
print("data loaded")

if(args[2]  == "rand"){
	my_partition_df <- Partition_2_dfs
	}else{
		my_partition_df <- Partition_2_chunk_dfs
	}

aatrainind <- which(my_partition_df$dataset %in% c("train", "valid"))
aa_train_data <- Partition_2_modified_dataset[aatrainind, -avoid_columns]
aa_test_data <-  Partition_2_modified_dataset[-aatrainind, -avoid_columns]
rm(Partition_2_modified_dataset)

aa_table <- table(aa_train_data$label)
aa_class_weight <- c(aa_table[2]/nrow(aa_train_data), aa_table[1]/nrow(aa_train_data))
names(aa_class_weight) <- c("neg", "pos")
aa_case_weight <- numeric(length = nrow(aa_train_data))
aa_case_weight[aa_train_data$label == "pos"] <- aa_class_weight[2]
aa_case_weight[aa_train_data$label == "neg"] <- aa_class_weight[1]
set.seed(my_random_seed)
print("training ...")
my_RF_model <- ranger(label ~ ., 
                      data = aa_train_data, 
                      num.trees = 1000,
                      max.depth = 8,
                      probability = TRUE, 
                      importance  = "impurity",
                      num.threads = 3,
                      save.memory= F,
                      case.weights = aa_case_weight)

print("training done")
if(length(avoid_columns) > 0){
	if(length(avoid_columns) <= 5){
		avoided <- paste(name_dic_p2[avoid_columns, 1], collapse = "__")
		}else{
			avoided <- paste(avoid_columns, collapse = "__")
		}
	
}else{
	avoided <- "full"
}
print(avoided)
pre_name <- unlist(strsplit(input_rdata, "\\."))[1]
pre_name_last <- strsplit(pre_name[1], "\\/")
pre_name_last <- unlist(pre_name_last)
pre_name_last <- pre_name_last[length(pre_name_last)]

aa_outputname <- paste0(model_store_address, pre_name_last, "___", partition_mode, "___", avoided, "___RFmodel.RData")
print(paste0("saving model to ", aa_outputname))
save(list = c("my_RF_model"), file = aa_outputname)





aaOrigoritrain_pred <- my_RF_model$predictions
aatrain_pred <- my_RF_model$predictions
aatrain_pred[is.na(aatrain_pred)] <- 0.5
my_RF_model$predictions <- aatrain_pred




training_perf <- perf_eval_ROC_PRC(model_fit = my_RF_model,
                                my_dataset = aa_train_data,
                                my_label = aa_train_data$label, 
                                train_perf = T,
                                file_name = paste0(plot_strore_address, pre_name_last, "___", partition_mode, "___", avoided, "___train_all.png"))

test_perf <- perf_eval_ROC_PRC(model_fit = my_RF_model,
                                                  my_dataset = aa_test_data,
                                                  my_label = aa_test_data[, ncol(aa_test_data)], 
                                                  train_perf = F,
                                                  file_name = paste0(plot_strore_address, pre_name_last, "___", partition_mode, "___", avoided, "___test_all.png"))



load("~/Documents/Shayan/BioInf/lncRNA/MM9_1kb_tiled_owner_3Mfilter.nosync.RData")

aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_train_data), names(MM9_1kb_tiled_owner_3Mfilter))]
training_perf_perclass <- perf_eval_ROC_PRC_perclass(model_fit = my_RF_model,
                                                my_dataset = aa_train_data,
                                                my_label = aa_train_data[, ncol(aa_train_data)], 
                                                file_name = paste0(plot_strore_address, pre_name_last, "___", partition_mode, "___", avoided, "___train_perClass.png"), 
                                                train_perf = T,
                                                class_name = aaclass_name)

aaclass_name <- MM9_1kb_tiled_owner_3Mfilter[match(rownames(aa_test_data), names(MM9_1kb_tiled_owner_3Mfilter))]
test_perf_perclass <- perf_eval_ROC_PRC_perclass(model_fit = my_RF_model,
                                                               my_dataset = aa_test_data,
                                                               my_label = aa_test_data[, ncol(aa_test_data)], 
                                                               file_name = paste0(plot_strore_address, pre_name_last, "___", partition_mode, "___", avoided, "___test_perClass.png"), 
                                                               train_perf = F,
                                                               class_name = aaclass_name)

variable_importance <- importance(my_RF_model)
names(variable_importance) <- name_dic_p2[match(names(variable_importance), name_dic_p2[, 2]),1]

png(filename = paste0(plot_strore_address, pre_name_last, "___", partition_mode, "___", avoided, "___variableImportanceTop50.png"),    # create PNG for the heat map        
    width = 6*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
par(mfrow = c(1,1), mar = c(4,9,4,4))
barplot(sort(variable_importance, decreasing = T)[1:50], horiz = T, las =2)
dev.off()

load("~/Documents/Shayan/BioInf/lncRNA/lncRNA_chosen_gt1k_uniqTiles.nosync.RData")
load("~/Documents/Shayan/BioInf/lncRNA/MM9_1kb_tiled_owner_labels_binary_3Mfilter.nosync.RData")


aamy_pred <- predict(my_RF_model, data=aa_test_data)
aapreddd <-  aamy_pred$predictions[, 2]
names(aapreddd) <- rownames(aa_test_data)
aapreddd2 <-  my_RF_model$predictions[, 2]
names(aapreddd2) <- rownames(aa_train_data)
aapreddd <- c(aapreddd, aapreddd2)
if(is.factor(my_partition_df$tile_name)){
  my_partition_df$tile_name <- levels(my_partition_df$tile_name)[as.numeric(my_partition_df$tile_name)]
}

aatst <- performance_viz_chromosome(all_tile_label=MM9_1kb_tiled_owner_labels_binary_3Mfilter,
                                    show_points = c("all"), 
                                    subset_points= c( "train", "test", "valid"),
                                    partition_dataset = my_partition_df,
                                    predicted_pos_probilbilty = aapreddd,
                                    positive_thersh = 0.5,
                                    specific_chr = character(0),
                                    my_file_name = paste0(plot_strore_address, pre_name_last, "___", partition_mode, "___", avoided, "___perChr_thr50.png"),
                                    lncRNA_coor_df= lncRNA_chosen_gt1k_uniqTiles,
                                    chr_sizes_file="~/Documents/Shayan/BioInf/lncRNA/mm9_chr_sizes.nosync.txt")

