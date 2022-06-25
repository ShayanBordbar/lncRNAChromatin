args <- commandArgs(trailingOnly = TRUE)
my_Learned_model_folder <- args[1]
my_output_file  <- args[2]
source("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/RF_evaluation_functions.R")
my_perf <- perf_gather(Learned_model_folder=my_Learned_model_folder, 
                        model_name = character(0), 
                        return_imp = T,
                        return_pred=F,
                        return_per_full=F)

save(list = c("my_perf"), file= my_output_file)