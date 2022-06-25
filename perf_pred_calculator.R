args <- commandArgs(trailingOnly = TRUE)
# get the predictions and performance
my_lncRNA_dataset_folder <- args[1]
my_Learned_model_folder <- args[2]
my_jobfile_address <- args[3]
my_aaplot_strore_address <-  args[4]
my_filtering_folder <- args[5]
my_output_file <- args[6]


performance_pred_list_lncRNA <- list()
aadir <- list.dirs(my_Learned_model_folder)
source("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/RF_evaluation_functions.R")
source("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/RBP_filter_func.R")

library(ranger, lib.loc = "~/R/library/")
length(aadir)
for(i in 2:length(aadir)){
  performance_pred_list_lncRNA[[i-1]] <- perf_eval_wrapper(Learned_model_folder = aadir[i],
                                                           only_imporatnce = F,
                                                           save_pred = T,
                                                           lncRNA_dataset_folder = my_lncRNA_dataset_folder,
                                                           model_name=character(0),
                                                           lncRNA_name = character(0),
                                                           jobfile_address = my_jobfile_address,
                                                           aaplot_strore_address = my_aaplot_strore_address,
                                                           filtering_folder = my_filtering_folder)
}

save(list = performance_pred_list_lncRNA, file = my_output_file)
#Rscript --vanilla perf_pred_calculator.R 
#/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CVfirst_dataset 
#/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CVfirst_Results/Learned_models
#/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/run_rf_p6_CVfirst_full_rem.job
#/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CVfirst_Results_comparison
#/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_CVfirst_TOBE_exp_FILTERED
#/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CVfirst_Results_comparison/Partition_6_CVfirst_all_Res.RData