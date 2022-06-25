# starting to work on the lncRNA project
setwd("~/Documents/Shayan/BioInf/lncRNA")
library(GenomicRanges)
library(biomaRt)
library(biomartr)
library(rtracklayer)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(Vennerable)
library(ChIPseeker)
########################################################################################################################################################################
########################################################################################################################################################################
##############################################        FUNCTIONS        #############################################################################################
########################################################################################################################################################################
########################################################################################################################################################################
DisChaLDis <- function(DisMat, k_closest = 1, dist_thresh = numeric(0), tie_break=F){
  # 'Dis'tance-based 'Cha-L'ncRNA:'D'NA C'is' interaction predictor
  #   uses the distance matrix
  #   to perform this prediction task (to predict a binary matrix with the same size as
  #   lncRNA_interaction_label_matrix) under various parameter settings and evaluate the
  #   results generally and for each lncRNA separately in terms of FP and FN 
  # it can use the number of lncRNAs on each chromosome to provide a random baseline to compare with.
  # it will also get the number of trans ints for each lncRNA and reports them, 
  #  during evaluation these interactions will not be considered
  # tie_break if True it will check if there are ties and deals with them.
  
  # DisMat is the distance matrix, where each row represents a tile, and each column represents a lncRNA
  #  and the values represent bp distance between the tile and lncRNA origin of biosynthesis
  # k_closest: integer, defining the number of k closest lncRNA to be assigned to each tile
  # dist_thresh: is a numeric value representing the max distance between a tile and its assosciated lncRNA eg. 1e6
  # default behavior is using 1 closest
  # output is a matrix with the same dimensions as DisMat
  
  pred_int <- matrix(0L, nrow = nrow(DisMat), ncol = ncol(DisMat))
  rownames(pred_int) <- rownames(DisMat)
  colnames(pred_int) <- colnames(DisMat)
  
  max_k <- rowSums(!is.na(DisMat))
  
  if(length(dist_thresh) == 0){
    # k closest mode
    for(i in 1:nrow(pred_int)){
      if(sum(!is.na(DisMat[i, ])) > 0){ # if there are any distances provided for the tile
        aa_c_dis <- DisMat[i, ]
        aa_c_dis[is.na(aa_c_dis)] <- max(aa_c_dis, na.rm = T) + 1
        aa_which_min <- sort(aa_c_dis, index.return=T, decreasing = F)$ix[1:(min(max_k[i],k_closest))]
        if(tie_break){ # if tie breaks need to be explored
          aatst <- which(aa_c_dis <= aa_c_dis[aa_which_min[min(max_k[i],k_closest)]])
          aatst2 <- (min(max_k[i],k_closest) >= length(aatst))
          if(!aatst2){
            aa_which_min <- sample(x = aatst, replace = F, size = min(max_k[i],k_closest))
          }
        }
        pred_int[i, aa_which_min] <- 1
      }
    }
  }else if(length(dist_thresh) == 1){
    # distance mode
    for(i in 1:nrow(pred_int)){
      if(sum(!is.na(DisMat[i, ])) > 0){ # if there are any distances provided for the tile
        aa_which_min <- which(DisMat[i, ] <= dist_thresh)
        pred_int[i, aa_which_min] <- 1
      }
    }
  }
  return(pred_int)
}
########################################################################################################################################################################
IntMatEval <- function(real_intMat, 
                       pred_intMat,
                       DisMat,
                       metaData = lncRNA_chosen_gt1k_uniqTiles){
  # function evaluates FP and FN rate for each lncRNA. 
  # metaData :is  a dataframe where each row corresponds to a column in DisMat and it
  #  contains information  chromosome position and  namings
  
  # it also uses DisMat to only look at interactions that a distance matric was
  #  provieded for (i.e. filter for cis interactions)
  stopifnot(all(dim(real_intMat) == dim(pred_intMat)),
            all(dim(real_intMat) == dim(DisMat)),
            nrow(metaData) == ncol(DisMat), 
            all(colnames(DisMat) == colnames(pred_intMat)))
  nu_data_pts_total <- sum(!is.na(DisMat))
  nu_data_pts_bycol <- colSums(!is.na(DisMat))
  
  # create a cis interaction matrix
  real_intMat_Cis <- real_intMat
  real_intMat_Cis[is.na(DisMat)] <- 20
  nu_cis_ints <- sum(real_intMat_Cis == 1)
  #TP
  TP_per_lnc <- numeric(ncol(real_intMat))
  
  names(TP_per_lnc) <- metaData$gene_name[match(colnames(DisMat),
                                                metaData$geneID)]
  TP_total <- sum((pred_intMat + real_intMat_Cis) == 2)
  #FP
  FP_per_lnc <- numeric(ncol(real_intMat))
  names(FP_per_lnc) <- names(TP_per_lnc)
  FP_total <- sum((pred_intMat == 1) & (real_intMat_Cis == 0))
  #FN
  FN_per_lnc <- numeric(ncol(real_intMat))
  names(FN_per_lnc) <-  names(TP_per_lnc)
  FN_total <- sum((pred_intMat == 0) & (real_intMat_Cis == 1))
  #TN
  TN_per_lnc <- numeric(ncol(real_intMat))
  names(TN_per_lnc) <- names(TP_per_lnc)
  TN_total <- sum((pred_intMat + real_intMat_Cis) == 0)
  
  # estimate random baseline (choosing randomly from the same chromosome) for TP, TN, FP and FN
  TP_random <- numeric(ncol(real_intMat))
  TN_random <- numeric(ncol(real_intMat))
  FP_random <- numeric(ncol(real_intMat))
  FN_random <- numeric(ncol(real_intMat))
  
  names(TP_random) <- names(TP_per_lnc)
  names(TN_random) <- names(TP_random)
  names(FP_random) <- names(TP_random)
  names(FN_random) <- names(TP_random)
  
  
  # get the number of lncRNAs on same chromosome for all lncRNAs
  nu_all_lnc_per_chr <- numeric(ncol(real_intMat))
  names(nu_all_lnc_per_chr) <- names(TP_per_lnc)
  for(i in 1:length(nu_all_lnc_per_chr)){
    nu_all_lnc_per_chr[i] <- sum(metaData$chromosome_name %in% metaData$chromosome_name[metaData$gene_name == names(nu_all_lnc_per_chr)[i]])
  }
  
  
  for(i in 1:ncol(real_intMat_Cis)){
    #random
    TP_random[i] <- (1/nu_all_lnc_per_chr[i])*(sum(real_intMat_Cis[, i] == 1))
    TN_random[i] <- ((nu_all_lnc_per_chr[i] - 1)/nu_all_lnc_per_chr[i])*sum(real_intMat_Cis[, i] == 0)
    FP_random[i] <- (1/nu_all_lnc_per_chr[i])*(sum(!is.na(DisMat[, i]) & real_intMat_Cis[, i] == 0))
    FN_random[i] <- ((nu_all_lnc_per_chr[i] - 1)/nu_all_lnc_per_chr[i])*sum(real_intMat_Cis[, i] == 1)
    #distance-based
    
    TP_per_lnc[i] <- sum((pred_intMat[, i] + real_intMat_Cis[, i]) == 2)
    TN_per_lnc[i] <- sum((pred_intMat[, i] + real_intMat_Cis[, i]) == 0)
    FP_per_lnc[i] <- sum((pred_intMat[, i] == 1) & (real_intMat_Cis[, i] == 0))
    FN_per_lnc[i] <- sum((pred_intMat[, i] == 0) & (real_intMat_Cis[, i] == 1))
  }
  myresults <- list(General = list(TP = TP_total,
                                   TN = TN_total, 
                                   FP = FP_total, 
                                   FN = FN_total,
                                   sens = (TP_total)/(TP_total + FN_total), 
                                   spec = (TN_total)/(TN_total + FP_total + 1), 
                                   prec = TP_total/(TP_total + FP_total), 
                                   accu = (TP_total + TN_total)/(TP_total + FN_total + TN_total + FP_total), 
                                   bal_accu= ((TP_total)/(TP_total + FN_total) + (TN_total)/(TN_total + FP_total + 1))/2,
                                   F1_score = (2*TP_total)/(2*TP_total + FN_total + FP_total)
  ),
  PerlncRNA = list(TP = TP_per_lnc,
                   TN = TN_per_lnc,
                   FP = FP_per_lnc, 
                   FN = FN_per_lnc, 
                   sens = (TP_per_lnc)/(TP_per_lnc + FN_per_lnc), 
                   spec = (TN_per_lnc)/(TN_per_lnc + FP_per_lnc + 1),
                   prec = TP_per_lnc/(TP_per_lnc + FP_per_lnc),
                   accu = (TP_per_lnc + TN_per_lnc)/(TP_per_lnc + FN_per_lnc + TN_per_lnc + FP_per_lnc),
                   bal_accu= ((TP_per_lnc)/(TP_per_lnc + FN_per_lnc) + (TN_per_lnc)/(TN_per_lnc + FP_per_lnc + 1))/2,
                   F1_score = (2*TP_per_lnc)/(2*TP_per_lnc + FN_per_lnc + FP_per_lnc),
                   TP_rand = TP_random,
                   TN_rand = TN_random,
                   FP_rand = FP_random, 
                   FN_rand = FN_random, 
                   sens_rand = (TP_random)/(TP_random + FN_random), 
                   spec_rand = (TN_random)/(TN_random + FP_random + 1),
                   prec_rand = TP_random/(TP_random + FP_random),
                   accu_rand = (TP_random + TN_random)/(TP_random + FN_random + TN_random + FP_random),
                   bal_accu_rand= ( ((TP_random)/(TP_random + FN_random)) + ((TN_random)/(TN_random + FP_random + 1)) )/2,
                   F1_score_rand = (2*TP_random)/(2*TP_random + FN_random + FP_random)
  ),
  Nu_per_chr_lncRNA = nu_all_lnc_per_chr)
  return(myresults)
}
########################################################################################################################################################################
# example
aatsttt <- DisChaLDis(DisMat = lncRNA_interaction_distance_matrix_0Removed_cis, 
                      k_closest = 1)
aacmp2 <- IntMatEval(real_intMat = lncRNA_interaction_label_matrix_cis,
                     pred_intMat = aatsttt,
                     DisMat = lncRNA_interaction_distance_matrix_0Removed_cis)



########################################################################################################################################################################
IntMatEval_plot <- function(eval_out, my_filename){
  # eval_out : is the output of IntMatEval funtion
  # filename : is the path to the plot file
  aa_bp_sens <- rbind(eval_out$PerlncRNA$sens, 
                      eval_out$PerlncRNA$sens_rand)[,sort(eval_out$Nu_per_chr_lncRNA,
                                                          decreasing = F, 
                                                          index.return = T)$ix]
  aa_bp_spec <- rbind(eval_out$PerlncRNA$spec, 
                      eval_out$PerlncRNA$spec_rand)[,sort(eval_out$Nu_per_chr_lncRNA,
                                                          decreasing = F, 
                                                          index.return = T)$ix]
  aa_bp_prec <- rbind(eval_out$PerlncRNA$prec, 
                      eval_out$PerlncRNA$prec_rand)[,sort(eval_out$Nu_per_chr_lncRNA,
                                                          decreasing = F, 
                                                          index.return = T)$ix]
  aa_bp_accu <- rbind(eval_out$PerlncRNA$accu, 
                      eval_out$PerlncRNA$accu_rand)[,sort(eval_out$Nu_per_chr_lncRNA,
                                                          decreasing = F, 
                                                          index.return = T)$ix]
  aa_bp_bal_accu <- rbind(eval_out$PerlncRNA$bal_accu, 
                          eval_out$PerlncRNA$bal_accu_rand)[,sort(eval_out$Nu_per_chr_lncRNA,
                                                                  decreasing = F, 
                                                                  index.return = T)$ix]
  aa_bp_F1_score <- rbind(eval_out$PerlncRNA$F1_score, 
                          eval_out$PerlncRNA$F1_score_rand)[,sort(eval_out$Nu_per_chr_lncRNA,
                                                                  decreasing = F, 
                                                                  index.return = T)$ix]
#  par(mfrow = c(6,1), mar = c(2,4,2,2))
  png(filename = my_filename,    # create PNG for the heat map        
      width = 8*300,        # 5 x 300 pixels
      height = 12*300,
      res = 300,            # 300 pixels per inch
      pointsize = 10)        # smaller font size
  
  par(mfrow = c(5,2), mar = c(3,6,5,3))
  barplot(rbind(eval_out$PerlncRNA$TP, eval_out$PerlncRNA$TP_rand)[,sort(eval_out$Nu_per_chr_lncRNA,
                                                                     decreasing = F,
                                                                     index.return = T)$ix], las = 2, main= "TP", beside = T
                 ,legend.text = c("distance", "random")
  )
  barplot(rbind(eval_out$PerlncRNA$FP, eval_out$PerlncRNA$FP_rand)[,sort(eval_out$Nu_per_chr_lncRNA,
                                                                     decreasing = F, 
                                                                     index.return = T)$ix], las = 2, main= "FP", beside = T)
  barplot(rbind(eval_out$PerlncRNA$TN, eval_out$PerlncRNA$TN_rand)[,sort(eval_out$Nu_per_chr_lncRNA,
                                                                     decreasing = F, 
                                                                     index.return = T)$ix], las = 2, main= "TN", beside = T)
  barplot(rbind(eval_out$PerlncRNA$FN, eval_out$PerlncRNA$FN_rand)[,sort(eval_out$Nu_per_chr_lncRNA,
                                                                     decreasing = F, 
                                                                     index.return = T)$ix], las = 2, main= "FN", beside = T)
  barplot(aa_bp_sens,
          beside = T,
         # legend.text = c("dist-based", "random"), 
          las= 2, 
          ylim = c(0,1),
          main = "Sensitivity")
  abline(h=seq(0,1,0.1), lwd = 0.6, lty = 2, col= "grey")
  barplot(aa_bp_spec,
          beside = T,
          #legend.text = c("dist-based", "random"), 
          las= 2, 
          ylim = c(0,1),
          main = "Specificity")
  abline(h=seq(0,1,0.1), lwd = 0.6, lty = 2, col = "grey")
  barplot(aa_bp_prec,
          beside = T,
          #legend.text = c("dist-based", "random"), 
          las= 2, 
          ylim = c(0,1),
          main = "Precision")
  abline(h=seq(0,1,0.1), lwd = 0.6, lty = 2, col = "grey")
  barplot(aa_bp_accu,
          beside = T,
          #legend.text = c("dist-based", "random"), 
          las= 2, 
          ylim = c(0,1),
          main = "Accuracy")
  abline(h=seq(0,1,0.1), lwd = 0.6, lty = 2, col = "grey")
  barplot(aa_bp_bal_accu,
          beside = T,
          #legend.text = c("dist-based", "random"), 
          las= 2, 
          ylim = c(0,1),
          main = "Balanced Accuracy")
  abline(h=seq(0,1,0.1), lwd = 0.6, lty = 2, col = "grey")
  barplot(aa_bp_F1_score,
          beside = T,
          #legend.text = c("dist-based", "random"), 
          las= 2, 
          ylim = c(0,1),
          main = "F1 score")
  abline(h=seq(0,1,0.1), lwd = 0.6, lty = 2, col = "grey")
  dev.off()
}
########################################################################################################################################################################
########################################################################################################################################################################
distance_lnc_tile_TAD  <- function(lncRNA_TAD_Assgined = lncRNA_chosen_gt1k_TADnu, 
                                   Tile_TAD_Assigned = MM9_1kb_tiled_TADnu, 
                                   lncRNA_metaData = lncRNA_chosen_gt1k_uniqTiles){
  # lncRNA_TAD_Assgined is matrix or dataframe, with one row per lncRNA and one column per TAD assignment scheme: entires indicate the TAD number under the column assigment map
  # Tile_TAD_Assigned is matrix or dataframe, with one row per tile and one column per TAD assignment scheme: entires indicate the TAD number under the column assigment map
  # lncRNA_metaData : is a dataframe that notes the chromosome on which each lncRNA resides
  # outputs a distance matrix that can then be furhter evaluated
  out_dis_list <- list()
  
  
  aa_tile_char <- rownames(MM9_1kb_tiled_TADnu)
  aa_tilenames_split <- strsplit(aa_tile_char, "_")
  aa_tilenames_split_2 <- lapply(aa_tilenames_split, function(x) x[-length(x)])
  aa_tilenames_split_3 <- mapply(paste, aa_tilenames_split_2, collapse = "_")
  
  for(i in 1:ncol(Tile_TAD_Assigned)){ # going through different TAD assignments
    print(i)
    out_dis_list[[i]] <- matrix(nrow = nrow(Tile_TAD_Assigned),
                                ncol = nrow(lncRNA_TAD_Assgined))
    rownames(out_dis_list[[i]]) <- rownames(Tile_TAD_Assigned)
    colnames(out_dis_list[[i]]) <- rownames(lncRNA_TAD_Assgined)
    for(j in 1:ncol(out_dis_list[[i]])){
      
      aa_cur_chr <- lncRNA_chosen_gt1k_uniqTiles$chromosome_name[lncRNA_chosen_gt1k_uniqTiles$gene_name %in% colnames(out_dis_list[[i]])[j]]
      out_dis_list[[i]][aa_tilenames_split_3 %in% aa_cur_chr, j] <- abs(Tile_TAD_Assigned[aa_tilenames_split_3 %in% aa_cur_chr, i] - lncRNA_TAD_Assgined[j,i])
      
    }
  }
  return(out_dis_list)
}


aatst <- distance_lnc_tile_TAD(lncRNA_TAD_Assgined = lncRNA_chosen_gt1k_TADnu, 
                               Tile_TAD_Assigned = MM9_1kb_tiled_TADnu, 
                               lncRNA_metaData = lncRNA_chosen_gt1k_uniqTiles)

########################################################################################################################################################################
########################################################################################################################################################################
add_pair_features <- function(lncRNA_feat, tile_feat, my_network_list = numeric(0), presence_thr=0){
  # lncRNA_feat is the feature set for lncRNA of interest --> numeric
  # tile_feat is the feature set for the tiles (matrix with ncol(tile_feat) == length(lncRNA_feat))
  # my_network_list is a set of edges over which the aggregation takes place, if nothing is porvided (and if name name_pattern is not provided) homo pairing of all feature will be done
  stopifnot(ncol(tile_feat) == length(lncRNA_feat))
  # if(length(my_network_list) == 0){
  #   my_network_list = list()
  #   my_network_list[[1]] <- cbind(c(1: length(lncRNA_feat)), 
  #                                 c(1: length(lncRNA_feat)))
  #   names(my_network_list) <- "All_Homo_pairs"
  # }
  # pair_feat <- matrix(nrow = nrow(tile_feat), ncol = length(my_network_list))
  # colnames(pair_feat) <- names(my_network_list)
  pair_feat <- numeric(length = nrow(tile_feat))
  for(cur_tile in 1:nrow(tile_feat)){
    print(cur_tile)
    #for(cur_net in 1:length(my_network_list)){
      # pair_feat[cur_tile, cur_net] <- sum((lncRNA_feat[my_network_list[[cur_net]][, 1]] > presence_thr) & 
      #                                       (tile_feat[cur_tile, my_network_list[[cur_net]][, 2]] > presence_thr))
    pair_feat[cur_tile] <- sum((lncRNA_feat > presence_thr) & 
                                 (tile_feat[cur_tile,] > presence_thr))
    #}
  }
  return(pair_feat)
}

########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################



# Getting genomic coordinates and names for all the lnc RNAs from UCSC: downleaded separately for mm9(mouse_lncRNA_ucsc_mm9) and mm10 (mouse_lncRNA_ucsc_mm10) using table browser with the following filters:
#  In the "name" field, set "does" to match "NR_*".
# In the "Free-form query:" box, enter "txEnd - txStart >=200" 
aa_lnc_mouse_mm10 <- read.table("~/Documents/Shayan/BioInf/lncRNA/mouse_lncRNA_ucsc_mm10",stringsAsFactors = F)
#colanmes(aa_lnc_mouse_mm10)<- 
aana <-   unlist(strsplit("bin name chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds score name2 cdsStartStat cdsEndStat exonFrames", split=" "))
colnames(aa_lnc_mouse_mm10) <- aana
#write.table(aa_lnc_mouse_mm10$name, file = "~/Documents/Shayan/BioInf/lncRNA/lncRNA_names_refseq.txt", quote = F, row.names = F, col.names = F)
sum( aa_lnc_mouse_mm10$name2 %in% "Malat1")
aa_lnc_mouse_mm10_mod <- aa_lnc_mouse_mm10[,c("chrom", "txStart", "txEnd", "strand", "name2", "name")]
names(aa_lnc_mouse_mm10_mod)[c(2,3)] <- c("start", "end")
aa_lnc_mouse_mm10_mod$name3 <- unlist(lapply(strsplit(aa_lnc_mouse_mm10_mod$name, "\\."), "[[", 1))
#write.table(aa_lnc_mouse_mm10_mod, file = "~/Documents/Shayan/BioInf/lncRNA/lncRNA_mousemm10.bed", quote = F, row.names = F, col.names = F)

aa_lnc_mouse_mm9 <- read.table("~/Documents/Shayan/BioInf/lncRNA/mouse_lncRNA_ucsc_mm9",stringsAsFactors = F)
colnames(aa_lnc_mouse_mm9) <- aana
sum( aa_lnc_mouse_mm9$name2 %in% "Malat1")
aa_lnc_mouse_mm9_mod <- aa_lnc_mouse_mm9[,c("chrom", "txStart", "txEnd", "strand", "name2", "name")]
names(aa_lnc_mouse_mm9_mod)[c(2,3)] <- c("start", "end")
#write.table(aa_lnc_mouse_mm9_mod, file = "~/Documents/Shayan/BioInf/lncRNA/lncRNA_mousemm9.bed", quote = F, row.names = F, col.names = F)

sum(! aa_lnc_mouse_mm9_mod$name %in%  aa_lnc_mouse_mm10_mod$name3)
aa_lnc_mouse_mm9_mod$name2[! aa_lnc_mouse_mm9_mod$name %in%  aa_lnc_mouse_mm10_mod$name3]
aa_lnc_mouse_mm10_mod$name2[(! aa_lnc_mouse_mm10_mod$name3 %in%  aa_lnc_mouse_mm9_mod$name)]

sum(! aa_lnc_mouse_mm10_mod$name2 %in%  aa_lnc_mouse_mm9_mod$name2)

lncRNA_mouse_mm9 <- aa_lnc_mouse_mm9_mod
lncRNA_mouse_mm10 <- aa_lnc_mouse_mm10_mod
lncRNA_mouse_mm9_GR <- makeGRangesFromDataFrame(lncRNA_mouse_mm9, keep.extra.columns = T)
lncRNA_mouse_mm10_GR <- makeGRangesFromDataFrame(lncRNA_mouse_mm10, keep.extra.columns = T)

rm(list = c("aa_lnc_mouse_mm10", "aana", "aa_lnc_mouse_mm9", "aa_lnc_mouse_mm9_mod", "aa_lnc_mouse_mm10_mod"))
########################################################################################################################################################################
########################################################################################################################################################################
# read the GRID-seq results for mESC, downloaded from: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2396700&format=file&file=GSM2396700%5FmESC%5Fmerged%2Eghits%2Epkbin%2Enet%2Etxt%2Egz
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2396700

GRIDseq_mm9 <- read.table( "~/Documents/Shayan/BioInf/lncRNA/GRID-seq/GSM2396700_mESC_merged.ghits.pkbin.net.txt",stringsAsFactors = F)
colnames(GRIDseq_mm9) <- unlist(strsplit("RNA_Chrom, RNA_Start, RNA_End, GeneID, RNA_Total_Reads, RNA_Strand, DNA_Chrom, DNA_Bin, Value", ", "))

GRIDseq_mm9_GR <- makeGRangesFromDataFrame(GRIDseq_mm9, keep.extra.columns = T)
GRIDseq_mm9_GR_uniq <- unique(GRIDseq_mm9_GR)
aaov1 <- findOverlaps(query = (GRIDseq_mm9_GR), subject = lncRNA_mouse_mm9_GR)
length(aaov1@from)
aaov2 <- findOverlaps(query = (GRIDseq_mm9_GR_uniq), subject = lncRNA_mouse_mm9_GR)
length(unique(aaov2@from))
length((aaov2@from))
length((aaov2@to))


# read chromatin enriched RNA
GRIDseq_chromatin_enriched_RNA_mm9 <- read.csv(file = "~/Documents/Shayan/BioInf/lncRNA/GRID-seq/mm9_chromatin_enriched_RNA.csv", header = T, stringsAsFactors = F)
sum(GRIDseq_chromatin_enriched_RNA_mm9$Gene.ID %in% GRIDseq_mm9_GR$GeneID)

sum(GRIDseq_mm9_GR_uniq$GeneID %in% GRIDseq_chromatin_enriched_RNA_mm9$Gene.ID)
head(GRIDseq_mm9_GR_uniq$GeneID)

GRIDseq_mm9_lnc <- GRIDseq_mm9[GRIDseq_mm9$GeneID %in% aa_only_long_nc,]
GRIDseq_mm9_lnc_filtered <- GRIDseq_mm9_lnc[GRIDseq_mm9_lnc$Value >= 2,]
#setwd("~/Documents/Shayan/BioInf/lncRNA/GRID-seq/")
# write.table(x = GRIDseq_mm9_lnc_filtered,file = "~/Documents/Shayan/BioInf/lncRNA/GRID-seq/GRID_lncRNA_value_gt2.txt", 
#             col.names = F, row.names = F, quote = F)
#system("sort -k7,8 GRID_lncRNA_value_gt2.txt > GRID_lncRNA_value_gt2_sorted.txt")
#system("join -1 7,8 GRID_lncRNA_value_gt2.txt > GRID_lncRNA_value_gt2_sorted.txt")
MM9_1kb_tiled <- read.table("~/Documents/Shayan/BioInf/lncRNA/mm9_1kb_tiles.bed")
names()
GRIDseq_mm9_lnc_filtered_tiled <- left_join( x = GRIDseq_mm9_lnc_filtered,
                                             y = MM9_1kb_tiled, 
                                            by = c("DNA_Chrom" = "V1", "DNA_Bin" = "V2"))
GRIDseq_mm9_lnc_filtered_tiled <- GRIDseq_mm9_lnc_filtered_tiled[,-c(10)]
sum(is.na(GRIDseq_mm9_lnc_filtered_tiled$V4))

########################################################################################################################################################################
# using biomart for annotating
library(biomartr)

listEnsembl()
ensembl <- useEnsembl(biomart = "ensembl")
#(mart = ensembl, pattern = "hsapiens")
aamart <- biomartr::getDatasets("ENSEMBL_MART_ENSEMBL")
"mmusculus_gene_ensembl"

mouse_mart <- biomaRt::useDataset(  dataset = "mmusculus_gene_ensembl",         
                                    mart    = useMart("ENSEMBL_MART_ENSEMBL",      
                                                      host    = "www.ensembl.org"))
aa_all_attrib <- biomaRt::listAttributes(mouse_mart)    

aa_all_filters <- biomaRt::listFilters(mouse_mart)  

organismAttributes(organism = "Mus musculus", update = T)
                                                
aa_my_att <- c("ensembl_gene_id", "ensembl_transcript_id", "chromosome_name", 
               "start_position", "end_position", "strand","external_gene_name",
               "external_gene_source", "gene_biotype", "transcript_biotype")

aats <- getBM(attributes = aa_my_att, 
      filters = "ensembl_gene_id", 
      values = GRIDseq_mm9_GR_uniq$GeneID, 
      mart = mouse_mart)

# get only lncRNAs
aaallow <- intersect(unique(aats$transcript_biotype), unique(RADICLseq_mm10$V9))

aats2 <- aats[((aats$transcript_biotype %in% aaallow | aats$gene_biotype %in% aaallow) & aats$gene_biotype != "protein_coding") ,]
nrow(aats2)
aa_only_long_nc <- unique(aats2$ensembl_gene_id)


aaradic <- unlist(lapply(strsplit(unique(RADICLseq_mm10$V7), "\\."), "[[", 1))

aats3 <- getBM(attributes = aa_my_att, 
              filters = "ensembl_gene_id", 
              values = aaradic, 
              mart = mouse_mart)

sum(aaradic %in% aats3$ensembl_gene_id)
table(aats3$gene_biotype)
table(aats3$transcript_biotype)

########################################################################################################################################################################
########################################################################################################################################################################
# read the results of RADICL-Seq

aafiles_full <- list.files("~/Documents/Shayan/BioInf/lncRNA/RADICL-seq/", pattern = "*.gz", full.names = T)
aafiles <- list.files("~/Documents/Shayan/BioInf/lncRNA/RADICL-seq/", pattern = "*.gz", full.names = F)
aafiles_2 <- unlist(strsplit(aafiles, split = ".gz"))

# write bash file to get only lncRNA related interactions in all of them
setwd("~/Documents/Shayan/BioInf/lncRNA/RADICL-seq/")
cat(c("#!/bin/bash\n"), file = "getlncRNA_ints.sh",append = F)
for(i in 1:length(aafiles_full)){
  cat(c("gunzip ", aafiles_full[i], "\n"), file = "getlncRNA_ints.sh",sep = " ", append = T)
  cat(c("awk '$8 == %long_ncRNA% {print}' ", aafiles_2[i], " > ", paste0(aafiles_2[i], "_lnc"), "\n"), file = "getlncRNA_ints.sh",sep = " ", append = T)
  cat(c("gzip ", aafiles_2[i], "\n"), file = "getlncRNA_ints.sh",sep = " ", append = T)
}

setwd("~/Documents/Shayan/BioInf/lncRNA")

RADICLseq_mm10 <-  read.table("GSE132190_all_mESC_lncRNA_interactions.txt", header = F, stringsAsFactors = F)
#colnames(RADICLseq_mm10) <- ""
table(RADICLseq_mm10$V9)

head(RADICLseq_mm10)



# convert the DNA part of RADICL from mm10 to mm9 

RADICLseq_mm10_DNA_GR <- makeGRangesFromDataFrame(RADICLseq_mm10[, c("V11", "V12", "V13")],
                                                  start.field = "V12", 
                                                  seqnames.field = "V11", 
                                                  end.field = "V13")
names(RADICLseq_mm10_DNA_GR)<- c(1:length(RADICLseq_mm10_DNA_GR))
aa_chain <- import.chain("~/Documents/Shayan/BioInf/Liftover_chain_files/mm10ToMm9.over.chain")
aaRADICLseq_mm9_DNA_GR_liftover  <- liftOver(x = RADICLseq_mm10_DNA_GR, chain = aa_chain)
#aal <- unlist(lapply(RADICLseq_mm9_DNA_GR, length))
RADICLseq_mm9_DNA_GR <- RADICLseq_mm10_DNA_GR
RADICLseq_mm10_to_mm9_stat <- numeric(length(RADICLseq_mm10_DNA_GR))
save(list = c("RADICLseq_mm9_DNA_GR", "aaRADICLseq_mm9_DNA_GR_liftover", "RADICLseq_mm10_to_mm9_stat"),
     file = "~/Documents/Shayan/BioInf/lncRNA/RADICL-seq/conversion_to_mm9_pre.RData")
# for(i in 1:length(RADICLseq_mm9_DNA_GR)){
#   print(i)
#   if(length(aaRADICLseq_mm9_DNA_GR_liftover[[i]]) > 0){
#     RADICLseq_mm9_DNA_GR[i] <- aaRADICLseq_mm9_DNA_GR_liftover[[i]][1]
#   }else{
#     print("no")
#     RADICLseq_mm10_to_mm9_stat[i] <- 1
#   }
# }



names(aaRADICLseq_mm9_DNA_GR_liftover) <- c(1:length(aaRADICLseq_mm9_DNA_GR_liftover))
names(RADICLseq_mm10_DNA_GR) <- c(1:length(aaRADICLseq_mm9_DNA_GR_liftover))

aaRADICLseq_mm9_DNA_GR_liftover_un <- unlist(aaRADICLseq_mm9_DNA_GR_liftover, use.names = T)


RADICLseq_mm10_DNA_df <- as.data.frame(RADICLseq_mm10_DNA_GR)
RADICLseq_mm10_DNA_df$index <- rownames(RADICLseq_mm10_DNA_df)
aaRADICLseq_mm9_DNA_GR_liftover_df <- as.data.frame(aaRADICLseq_mm9_DNA_GR_liftover_un)
aaRADICLseq_mm9_DNA_GR_liftover_df$index <- names(aaRADICLseq_mm9_DNA_GR_liftover_un)
head(aaRADICLseq_mm9_DNA_GR_liftover_df)

RADICLseq_mm10_mm9_DNA <- left_join( x = RADICLseq_mm10_DNA_df,
                                             y = aaRADICLseq_mm9_DNA_GR_liftover_df, 
                                             by = c("index"))

sum(is.na(RADICLseq_mm10_mm9_DNA$seqnames.y))

RADICLseq_mm10$index <- as.character(c(1:nrow(RADICLseq_mm10)))


RADICLseq_mm10_mm9 <-  left_join( x = RADICLseq_mm10,
                                  y = aaRADICLseq_mm9_DNA_GR_liftover_df, 
                                  by = c("index"))
ncol(RADICLseq_mm10_mm9)
head(RADICLseq_mm10_mm9)


RADICLseq_mm10_mm9 <- RADICLseq_mm10_mm9[!is.na(RADICLseq_mm10_mm9$seqnames), ] 

aaRADICLseq_mm10_mm9_GR <- makeGRangesFromDataFrame(df = RADICLseq_mm10_mm9[, c("seqnames", "start", "end", "index")], 
                         keep.extra.columns = T)

names(MM9_1kb_tiled) <- c("chr", "start", "end", "tile")
MM9_1kb_tiled_GR <- makeGRangesFromDataFrame(df = MM9_1kb_tiled, keep.extra.columns = T)
end(MM9_1kb_tiled_GR) <- end(MM9_1kb_tiled_GR) - 1
aaov4 <- findOverlaps(query = aaRADICLseq_mm10_mm9_GR, subject = MM9_1kb_tiled_GR)

all(aaov4@from == c(1:length(aaov4@from)))

aanewtile <- MM9_1kb_tiled_GR$tile[aaov4@to]

RADICLseq_mm10_mm9$kbTile <- aanewtile
########################################################################################################################
# RADICLseq_mm10_mm9 contains all information for RADICL-seq
head(RADICLseq_mm10_mm9)
max((RADICLseq_mm10_mm9$V18))
########################################################################################################################
# GRIDseq_mm9_lnc_filtered_tiled contains all information for GRID-seq
head(GRIDseq_mm9_lnc_filtered_tiled)
#####################################################################################################################
######################################################################################################################
#####################################################################################################################

# combine the two experiments and aggregate based on geneID, or tile

aarad <- RADICLseq_mm10_mm9[, c("seqnames", "start", "end", "V6", "V7", "V16", "V18", "kbTile")]
aagrid <- GRIDseq_mm9_lnc_filtered_tiled[, c("RNA_Chrom", "RNA_Start", "RNA_End",
                                             "RNA_Strand", "GeneID", "Value", "V4")]

aaradic <- unlist(lapply(strsplit((RADICLseq_mm10_mm9$V7), "\\."), "[[", 1))
aarad$V7 <- aaradic
names(aarad) <- c("chr","start", "end", "strand", "geneID", "experiment", "score", "kbTile")

aagrid$experiment <- rep("GRIDseq", nrow(aagrid))
aagrid <- aagrid[, c(1,2,3,4,5,8,6,7)]
names(aagrid) <- names(aarad)

RADICL_GRID_Tiled_df <- rbind(aarad, aagrid)
RADICL_GRID_Tiled_df_short <- RADICL_GRID_Tiled_df[,c(5,6,7,8)]
RADICL_GRID_Tiled_df_short$kbTile <- as.character(levels(RADICL_GRID_Tiled_df_short$kbTile)[as.numeric(RADICL_GRID_Tiled_df_short$kbTile)])

########################################################################################################################################################################
########################################################################################################################################################################
# aggregate RADICL_GRID_Tiled_df_short based on GeneID and then sort based on the number of interactions
head(RADICL_GRID_Tiled_df_short)
RADICL_GRID_Tiled_df_short_noGRID <- RADICL_GRID_Tiled_df_short[RADICL_GRID_Tiled_df_short$experiment != "GRIDseq",]
RADICL_GRID_Tiled_df_short_onlygrid
RADICL_GRID_Tiled_df_short_byGeneID <- aggregate(x = RADICL_GRID_Tiled_df_short,
                                                      by = RADICL_GRID_Tiled_df_short[c("geneID")], 
                                                        FUN = c)

RADICL_GRID_Tiled_df_short_noGRID_byGeneID <- aggregate(x = RADICL_GRID_Tiled_df_short_noGRID,
                                                 by = RADICL_GRID_Tiled_df_short_noGRID[c("geneID")], 
                                                 FUN = c)

aa_tab <- lapply(RADICL_GRID_Tiled_df_short_byGeneID$experiment, table)
RADICL_GRID_Tiled_df_short_byGeneID$Exp_table <- aa_tab
aa_len <- unlist(lapply(RADICL_GRID_Tiled_df_short_byGeneID$experiment, length))
aa_len_uniq <- unlist(lapply(lapply(RADICL_GRID_Tiled_df_short_byGeneID$experiment, unique), length))

RADICL_GRID_Tiled_df_short_byGeneID_modified <- RADICL_GRID_Tiled_df_short_byGeneID[, c("geneID", "kbTile", "Exp_table")]
aa_tab <- lapply(RADICL_GRID_Tiled_df_short_byGeneID$kbTile, table)
RADICL_GRID_Tiled_df_short_byGeneID_modified$kbTile_table <- aa_tab
RADICL_GRID_Tiled_df_short_byGeneID_modified <- RADICL_GRID_Tiled_df_short_byGeneID_modified[, c("geneID", "Exp_table", "kbTile_table")]
RADICL_GRID_Tiled_df_short_byGeneID_modified$nu_ints_total <- aa_len
RADICL_GRID_Tiled_df_short_byGeneID_modified$nu_exp_types <- aa_len_uniq
aa_tile_uniq <- unlist(lapply(lapply(RADICL_GRID_Tiled_df_short_byGeneID$kbTile, unique), length))
RADICL_GRID_Tiled_df_short_byGeneID_modified$nu_uniq_tiles <- aa_tile_uniq
aasind <- sort(aa_len, index.return=T, decreasing = T)$ix
RADICL_GRID_Tiled_df_short_byGeneID_modified_sorted <- RADICL_GRID_Tiled_df_short_byGeneID_modified[aasind,]

aa_tab <- lapply(RADICL_GRID_Tiled_df_short_noGRID_byGeneID$experiment, table)
RADICL_GRID_Tiled_df_short_noGRID_byGeneID$Exp_table <- aa_tab
aa_len <- unlist(lapply(RADICL_GRID_Tiled_df_short_noGRID_byGeneID$experiment, length))
aa_len_uniq <- unlist(lapply(lapply(RADICL_GRID_Tiled_df_short_noGRID_byGeneID$experiment, unique), length))

RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified <- RADICL_GRID_Tiled_df_short_noGRID_byGeneID[, c("geneID", "kbTile", "Exp_table")]
aa_tab <- lapply(RADICL_GRID_Tiled_df_short_noGRID_byGeneID$kbTile, table)
RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified$kbTile_table <- aa_tab
RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified <- RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified[, c("geneID", "Exp_table", "kbTile_table")]
RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified$nu_ints_total <- aa_len
RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified$nu_exp_types <- aa_len_uniq
aa_tile_uniq <- unlist(lapply(lapply(RADICL_GRID_Tiled_df_short_noGRID_byGeneID$kbTile, unique), length))
RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified$nu_uniq_tiles <- aa_tile_uniq
aasind <- sort(aa_len, index.return=T, decreasing = T)$ix
RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified_sorted <- RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified[aasind,]


par(mfrow = c(2,1))
hist(log10(aa_len), main = "histogram:log10 of total number of interactions")
hist(log10(aa_tile_uniq), main = "histogram:log10 of number of unique interacting 1kb tiles")


aa_my_att <- c("chromosome_name", "start_position", "end_position",
               "ensembl_gene_id","external_gene_name", "gene_biotype")

aats4 <- getBM(attributes = aa_my_att, 
              filters = "ensembl_gene_id", 
              values = RADICL_GRID_Tiled_df_short_byGeneID_modified_sorted$geneID, 
              mart = mouse_mart)
lncRNA_all_info_biomart <- aats4
aaa_name <-aats4$external_gene_name[match(RADICL_GRID_Tiled_df_short_byGeneID_modified_sorted$geneID, aats4$ensembl_gene_id)]
aaa_btyp <-aats4$gene_biotype[match(RADICL_GRID_Tiled_df_short_byGeneID_modified_sorted$geneID, aats4$ensembl_gene_id)]

RADICL_GRID_Tiled_df_short_byGeneID_modified_sorted$gene_name <- aaa_name
RADICL_GRID_Tiled_df_short_byGeneID_modified_sorted$gene_biotype <- aaa_btyp

aats4 <- getBM(attributes = aa_my_att, 
               filters = "ensembl_gene_id", 
               values = RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified_sorted$geneID, 
               mart = mouse_mart)
#lncRNA_all_info_biomart <- aats4
aaa_name <-aats4$external_gene_name[match(RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified_sorted$geneID, aats4$ensembl_gene_id)]
aaa_btyp <-aats4$gene_biotype[match(RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified_sorted$geneID, aats4$ensembl_gene_id)]

RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified_sorted$gene_name <- aaa_name
RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified_sorted$gene_biotype <- aaa_btyp



View(RADICL_GRID_Tiled_df_short_byGeneID_modified_sorted[1:100,c("nu_ints_total", "nu_exp_types", 
                                                                 "nu_uniq_tiles",
                                                                 "gene_name", "gene_biotype")])
View(RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified_sorted[1:100,c("nu_ints_total", "nu_exp_types", 
                                                                 "nu_uniq_tiles",
                                                                 "gene_name", "gene_biotype")])


aa_potentially_interesting <- c("Xist", "Neat1", "Malat1", "Linc-YY1", "Rmst", 
                                "Lunar1", "Trp53cor1", "Sammson", "Firre", "Meg3", 
                                "Khps1", "Coolair", "Airn", "Hotair", "Kcnq1ot1", 
                                "Hottip", "Bvht", "Fendrr", "Anril",
                                "Lncpint", "Chaer1", "pRNA")
aa_present <- aa_potentially_interesting[which(aa_potentially_interesting %in% RADICL_GRID_Tiled_df_short_byGeneID_modified_sorted$gene_name) ]

View(cbind(aa_present, match(aa_present, RADICL_GRID_Tiled_df_short_byGeneID_modified_sorted$gene_name)))



"ENSMUSG00000098393" %in% RADICL_GRID_Tiled_df_short_byGeneID_modified_sorted$geneID
"Chaer1" %in% RADICL_GRID_Tiled_df_short_byGeneID_modified_sorted$gene_name

View(RADICL_GRID_Tiled_df_short_byGeneID_modified_sorted[match(aa_present, RADICL_GRID_Tiled_df_short_byGeneID_modified_sorted$gene_name),
                                                         c("nu_ints_total", "nu_exp_types","nu_uniq_tiles",
                                                            "gene_name", "gene_biotype")])
# choosing the lncRNAs:  Ones with more than 1000 unique interacting tiles
sum(RADICL_GRID_Tiled_df_short_byGeneID_modified_sorted$nu_uniq_tiles > 1000)

lncRNA_chosen_gt1k_uniqTiles <- RADICL_GRID_Tiled_df_short_byGeneID_modified_sorted[(RADICL_GRID_Tiled_df_short_byGeneID_modified_sorted$nu_uniq_tiles > 1000),
                                                                                    c("geneID", "gene_name", "nu_uniq_tiles", "nu_ints_total", "nu_exp_types","gene_biotype"  )]
View(lncRNA_chosen_gt1k_uniqTiles)
write.table(x = lncRNA_chosen_gt1k_uniqTiles[,c("gene_name","geneID","nu_uniq_tiles","nu_ints_total")], file = "~/Desktop/lncRNA_paper/iScience/suppTable1.csv",
            append = F, quote = F, sep = ",", row.names = F, col.names = T)
lncRNA_chosen_gt1k_uniqTiles <- lncRNA_chosen_gt1k_uniqTiles[-c(which(lncRNA_chosen_gt1k_uniqTiles$gene_name == "4930519F16Rik")),]
View(lncRNA_chosen_gt1k_uniqTiles)

aaachr <- lncRNA_all_info_biomart[match(lncRNA_chosen_gt1k_uniqTiles$geneID, lncRNA_all_info_biomart$ensembl_gene_id), c("chromosome_name", "start_position", "end_position", "external_gene_name")]

lncRNA_chosen_gt1k_uniqTiles <- cbind(aaachr[, 1:3], lncRNA_chosen_gt1k_uniqTiles)

barplot(RADICL_GRID_Tiled_df_short_byGeneID_modified_sorted$nu_ints_total)
aas <-cumsum(RADICL_GRID_Tiled_df_short_byGeneID_modified_sorted$nu_ints_total)
png(filename = "~/Desktop/lncRNA_paper/Figs/cumulative_dist.png",    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 3*300,
    res = 300,            # 300 pixels per inch
    pointsize = 12)       
par(mfrow = c(1,1), mar = c(4,4,4,4),family = "Arial")
plot(aas/aas[length(aas)] * 100, type = "l", main = "Cumulative distribution of \nlncRNA-chromatin interactions", 
     xlab = "No. lncRNAs", ylab = "% interactions", lwd = 1.5)
abline(v = 28, col =2, lty = 2, lwd = 1 )
dev.off()

aas2 <-cumsum(RADICL_GRID_Tiled_df_short_byGeneID_modified_sorted$nu_uniq_tiles)
plot(aas2/aas2[length(aas2)] * 100, type = "l", main = "Cumulative distribution of \n unique lncRNA-chromatin interactions", 
     xlab = "No. lncRNAs", ylab = "% interactions", lwd = 1.5)

plot(aas2, type = "l", main = "Cumulative distribution of \n unique lncRNA-chromatin interactions", 
     xlab = "No. lncRNAs", ylab = "% interactions", lwd = 1.5)

plot(aas, type = "l", main = "Cumulative distribution of \n unique lncRNA-chromatin interactions", 
     xlab = "No. lncRNAs", ylab = "% interactions", lwd = 1.5)

####################################################################################################
# Continue with the 28 chosen lncRNA, aggreagate lncRNA,tile pairs
# then look at the total number of unique tiles
# create label profile for each tile
# basically create a matrix with one row for each tile, and one column for each lncRNA, then do heatmap, ...
# look at correlation between lncRNAs
length(unique(RADICL_GRID_Tiled_df_short_byGeneID__tile_pairs$kbTile))


RADICL_GRID_Tiled_df_short_ChosenLNC <- RADICL_GRID_Tiled_df_short[RADICL_GRID_Tiled_df_short$geneID %in% lncRNA_chosen_gt1k_uniqTiles$geneID, ]
RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs <- aggregate(x = RADICL_GRID_Tiled_df_short_ChosenLNC[c("geneID", "kbTile","experiment")],
                                                                       by = RADICL_GRID_Tiled_df_short_ChosenLNC[c("geneID", "kbTile")], 
                                                                       FUN = c)

nrow(RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs)
aatab <- lapply(RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs$experiment, table)
aalnp <- unlist(lapply(RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs$experiment, length))
aalnp_unq <- unlist(lapply(lapply(RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs$experiment, unique), length))
RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs$exp_table <- aatab
RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs$exp_number <- aalnp
RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs$exp_type_number <- aalnp_unq




hist(log10(RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs$exp_number),
     main = "histogram: log10 of number of experimental\n evidence supporting each interaction")
hist(RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs$exp_type_number ,
     main = "histogram: number of uniq experiment types\n supporting each interaction")
table(RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs$exp_type_number)

RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat <- RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs[, c("geneID", "kbTile", "exp_number", "exp_type_number")]

length(unique(RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat$kbTile))

aaunq_tile <- unique(RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat$kbTile)
aaunq_geneid <- unique(RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat$geneID)
lncRNA_interaction_label_matrix <- matrix(0L,nrow = length(aaunq_tile),
                                          ncol = length(aaunq_geneid))
rownames(lncRNA_interaction_label_matrix) <- aaunq_tile
colnames(lncRNA_interaction_label_matrix) <- aaunq_geneid
for(i in 1:ncol(lncRNA_interaction_label_matrix)){
  print(i)
  aa_cur_tiles <- RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat$kbTile[RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat$geneID == colnames(lncRNA_interaction_label_matrix)[i]]
  lncRNA_interaction_label_matrix[match(aa_cur_tiles, rownames(lncRNA_interaction_label_matrix)),i] <- 1
}
sum(lncRNA_interaction_label_matrix)

# compute correlation between lncRNAs binding profiles
aa_lnc_correlation <- matrix(nrow = ncol(lncRNA_interaction_label_matrix),
                             ncol = ncol(lncRNA_interaction_label_matrix))
# colnames(aa_lnc_correlation) <- colnames(lncRNA_interaction_label_matrix)
# rownames(aa_lnc_correlation) <- colnames(lncRNA_interaction_label_matrix)

par(mfrow = c(1,1))
table(rowSums(lncRNA_interaction_label_matrix))

for(i in 1:(ncol(aa_lnc_correlation) - 1)){
  for(j in (i+1):ncol(aa_lnc_correlation)){
    aaac <- cor.test(lncRNA_interaction_label_matrix[, i], lncRNA_interaction_label_matrix[, j])
    aa_lnc_correlation[i, j] <- aaac$estimate
    aa_lnc_correlation[j, i] <- aaac$estimate
  }
}


colnames(aa_lnc_correlation) <- lncRNA_chosen_gt1k_uniqTiles$gene_name[match(colnames(lncRNA_interaction_label_matrix), 
                                                                             lncRNA_chosen_gt1k_uniqTiles$geneID)]
rownames(aa_lnc_correlation) <- colnames(aa_lnc_correlation)
aa_lnc_correlation[is.na(aa_lnc_correlation)] <- 1
View(aa_lnc_correlation)
heatmap.2(aa_lnc_correlation, 
          trace = "none",
          Rowv = T,
          Colv = T, 
          dendrogram = "col",
          margins = c(8,8), 
          scale = "none",
          main = "lncRNA binding correlation")
# compute the same plot but only for cis interactions between lncRNAs on the same chr
aa_Ag <- aggregate(lncRNA_chosen_gt1k_uniqTiles[c("geneID", "gene_name")],
                   by = lncRNA_chosen_gt1k_uniqTiles[c("chromosome_name")], 
                   FUN = c)
aa_cormat <- matrix(nrow = nrow(lncRNA_chosen_gt1k_uniqTiles),
                    ncol = nrow(lncRNA_chosen_gt1k_uniqTiles))
rownames(aa_cormat) <- lncRNA_chosen_gt1k_uniqTiles$gene_name
colnames(aa_cormat) <- lncRNA_chosen_gt1k_uniqTiles$gene_name
for(i in 1:nrow(aa_Ag)){
  aa_cur_chr <- paste0("chr", aa_Ag$chromosome_name[i])
  aa_ind <- which(aa_tilenames_split_3 %in% aa_cur_chr)
  if(length(aa_Ag$geneID[[i]]) > 1){
    for(j in 1:(length(aa_Ag$geneID[[i]]))){
      for(k in j:(length(aa_Ag$geneID[[i]]))){
        aaac <- cor.test(lncRNA_interaction_label_matrix[aa_ind, aa_Ag$geneID[[i]][j]],
                         lncRNA_interaction_label_matrix[aa_ind, aa_Ag$geneID[[i]][k]])
        aa_cormat[aa_Ag$gene_name[[i]][j], aa_Ag$gene_name[[i]][k]] <- aaac$estimate
        aa_cormat[aa_Ag$gene_name[[i]][k], aa_Ag$gene_name[[i]][j]] <- aaac$estimate
      }
    }
  }
}
View(aa_cormat)
aa_cormat[is.na(aa_cormat)] <- 0

# choosing color
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
aaaaaacolor = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
aamycoltest <- sample(x = aaaaaacolor, size = n, replace = F)
barplot(rep(5, 17), names.arg = c(1:17), col = aamycoltest, las = 2)
####

aaccl <- lncRNA_chosen_gt1k_uniqTiles$chromosome_name
aaccl[aaccl == "X"] <- 20
aaccl <- aamycoltest[as.numeric(aaccl)]

heatmap.2(aa_cormat, 
          trace = "none",
          Rowv = T,
          Colv = T, 
          dendrogram = "both",
          margins = c(8,8), 
          scale = "none",
          ColSideColors = aaccl,
          RowSideColors = aaccl,
          main = "lncRNA cis binding correlation_per_chr")

##########################################################################################
# make a barplot colored for different chromosomes

aa_barmat <- matrix(nrow = length(levels(MM9_1kb_tiled$chr)), 
                    ncol= ncol(lncRNA_interaction_label_matrix))
colnames(aa_barmat) <- lncRNA_chosen_gt1k_uniqTiles$gene_name[match(colnames(lncRNA_interaction_label_matrix), 
                                                                    lncRNA_chosen_gt1k_uniqTiles$geneID)]
rownames(aa_barmat) <- levels(MM9_1kb_tiled$chr)

aa_tile_char <- levels(MM9_1kb_tiled$tile)[as.numeric(MM9_1kb_tiled$tile)]
aa_tilenames_split <- strsplit(aa_tile_char, "_")
aa_tilenames_split_2 <- lapply(aa_tilenames_split, function(x) x[-length(x)])
aa_tilenames_split_3 <- mapply(paste, aa_tilenames_split_2, collapse = "_")
par(mar = c(2,2,2,2))
pie(table(aa_tilenames_split_3), 
    main = "Distribution of chromosomes tiles mm9")

aa_tile_char <- rownames(lncRNA_interaction_label_matrix)
aa_tilenames_split <- strsplit(aa_tile_char, "_")
aa_tilenames_split_2 <- lapply(aa_tilenames_split, function(x) x[-length(x)])
aa_tilenames_split_3 <- mapply(paste, aa_tilenames_split_2, collapse = "_")
par(mar = c(2,2,2,2))
pie(table(aa_tilenames_split_3), 
    main = "Distribution of chromosomes\nbound at least once")


for(i in 1:ncol(lncRNA_interaction_label_matrix)){
  aa_cur_tiles <- aa_tilenames_split_3[which(lncRNA_interaction_label_matrix[,i] > 0)]
  for(j in 1:nrow(aa_barmat)){
    aa_barmat[j, i] <- sum( aa_cur_tiles %in% rownames(aa_barmat)[j])
  }
}

barplot(aa_barmat, las=2, main = "chromosome distribution of lncRNA binding sites")

aa_barmat_norm <- aa_barmat
for(i in 1:ncol(aa_barmat_norm)){
  aa_barmat_norm[,i] <- aa_barmat_norm[,i] * 100 / (sum(aa_barmat_norm[,i]))
}
par(mar = c(7,4,4,4))
barplot(aa_barmat_norm, las=2, main = "chromosome distribution of lncRNA binding sites\nNormalized")


aa_barmat <- aa_barmat[-which(rowSums(aa_barmat) == 0), ]

par(mfrow = c(7,4),mar = c(3,2,2,2))
for(i in 1:ncol(aa_barmat_norm)){
#  aa_cur_geneID <- lncRNA_chosen_gt1k_uniqTiles$geneID[lncRNA_chosen_gt1k_uniqTiles$gene_name == colnames(aa_barmat_norm)[i]]
  
  aa_cur_chr <- unique(lncRNA_all_info_biomart$chromosome_name[lncRNA_all_info_biomart$external_gene_name == colnames(aa_barmat_norm)[i]])
  
  barplot(aa_barmat[,i], las = 2, main = paste0(colnames(aa_barmat_norm)[i], "\nchr", aa_cur_chr))
}



par(mfrow = c(1,1),mar = c(4,4,4,4))
barplot(rowSums(aa_barmat), las = 2, 
        main = "Total lncRNA binding site distribution\nfor the 28 chosen ones")

par(mfrow = c(1,1),mar = c(4,4,4,4))
barplot(table(lncRNA_chosen_gt1k_uniqTiles$chromosome_name), las = 2,
        main = "chromosome distribution of origin of biosynthesis\n of the 28 considered lncRNAs")

##########################################################################################

# find the distance between each tile and the lncRNAs (out of 28) on the same chromosome, 
# store in a structure similar to feature matrix



aamtch <- match(lncRNA_chosen_gt1k_uniqTiles$gene_name, lncRNA_mouse_mm9$name2)
lncRNA_chosen_gt1k_uniqTiles$gene_name[which(is.na(aamtch))]

aa_coor_mm9_28 <- cbind(lncRNA_chosen_gt1k_uniqTiles$gene_name, lncRNA_mouse_mm9[match(lncRNA_chosen_gt1k_uniqTiles$gene_name, lncRNA_mouse_mm9$name2),1:3])


aalncRNA_chosen_gt1k_uniqTiles <- lncRNA_chosen_gt1k_uniqTiles
aalncRNA_chosen_gt1k_uniqTiles$chromosome_name <- paste0("chr", aalncRNA_chosen_gt1k_uniqTiles$chromosome_name)
aalncRNA_chosen_gt1k_uniqTiles <- makeGRangesFromDataFrame(aalncRNA_chosen_gt1k_uniqTiles, keep.extra.columns = T, 
                                                           seqnames.field = "chromosome_name",
                                                           start.field = "start_position",
                                                           end.field = "end_position")
aa_chain <- import.chain("~/Documents/Shayan/BioInf/Liftover_chain_files/mm10ToMm9.over.chain")
aalncRNA_chosen_gt1k_uniqTiles_mm9 <- liftOver(aalncRNA_chosen_gt1k_uniqTiles, chain = aa_chain)
aalncRNA_chosen_gt1k_uniqTiles_mm9 <- unlist(aalncRNA_chosen_gt1k_uniqTiles_mm9)

# now calculate the distance between each lncRNA and all tiles on the same chromosome

lncRNA_interaction_distance_matrix <- matrix(nrow = nrow(lncRNA_interaction_label_matrix),
                                             ncol = ncol(lncRNA_interaction_label_matrix))
rownames(lncRNA_interaction_distance_matrix) <- rownames(lncRNA_interaction_label_matrix)
colnames(lncRNA_interaction_distance_matrix) <- colnames(lncRNA_interaction_label_matrix)

aa_tile_char <- rownames(lncRNA_interaction_label_matrix)
aa_tilenames_split <- strsplit(aa_tile_char, "_")
aa_tilenames_split_2 <- lapply(aa_tilenames_split, function(x) x[-length(x)])
aa_tilenames_split_3 <- mapply(paste, aa_tilenames_split_2, collapse = "_")
MM9_1kb_tiled_GR$tile <- levels(MM9_1kb_tiled_GR$tile)[as.numeric(MM9_1kb_tiled_GR$tile)]
MM9_1kb_tiled$tile <- levels(MM9_1kb_tiled$tile)[as.numeric(MM9_1kb_tiled$tile)]
for(i in 1:ncol(lncRNA_interaction_distance_matrix)){
  aa_cur_qu <- aalncRNA_chosen_gt1k_uniqTiles_mm9[aalncRNA_chosen_gt1k_uniqTiles_mm9$geneID == colnames(lncRNA_interaction_distance_matrix)[i]]
  aa_cur_tn <- aa_tile_char[aa_tilenames_split_3 %in% as.character(seqnames(aa_cur_qu))]
  aa_cur_su <- MM9_1kb_tiled_GR[MM9_1kb_tiled_GR$tile %in% aa_cur_tn]
  aa_cur_dist <- GenomicRanges::distance(x = aa_cur_qu, y = aa_cur_su, select = "all")
  lncRNA_interaction_distance_matrix[match(aa_cur_su$tile, rownames(lncRNA_interaction_distance_matrix)), i] <- aa_cur_dist
}
########################################################################################################################################################
# make lncRNA_interaction_distance_matrix for all the tiles not just the interacting ones
lncRNA_interaction_distance_matrix_full <- matrix(nrow = nrow(MM9_1kb_tiled),
                                                  ncol = ncol(lncRNA_interaction_label_matrix))
rownames(lncRNA_interaction_distance_matrix_full) <- MM9_1kb_tiled$tile
colnames(lncRNA_interaction_distance_matrix_full) <-lncRNA_chosen_gt1k_uniqTiles$gene_name[match(colnames(lncRNA_interaction_label_matrix), lncRNA_chosen_gt1k_uniqTiles$geneID)]

aa_tilenames_split <- strsplit(MM9_1kb_tiled$tile, "_")
aa_tilenames_split_2 <- lapply(aa_tilenames_split, function(x) x[-length(x)])
aa_tilenames_split_3 <- mapply(paste, aa_tilenames_split_2, collapse = "_")
MM9_1kb_tiled_chrname <- aa_tilenames_split_3
for(i in 1:ncol(lncRNA_interaction_distance_matrix_full)){
  print(i)
  aa_cur_qu <- aalncRNA_chosen_gt1k_uniqTiles_mm9[aalncRNA_chosen_gt1k_uniqTiles_mm9$gene_name ==
                                                    colnames(lncRNA_interaction_distance_matrix_full)[i]]
  #aa_cur_tn <- aa_tile_char[aa_tilenames_split_3 %in% as.character(seqnames(aa_cur_qu))]
  aa_cur_su <- MM9_1kb_tiled_GR[aa_tilenames_split_3 %in% as.character(seqnames(aa_cur_qu))]
  aa_cur_dist <- GenomicRanges::distance(x = aa_cur_qu, y = aa_cur_su, select = "all")
  lncRNA_interaction_distance_matrix_full[match(aa_cur_su$tile, rownames(lncRNA_interaction_distance_matrix_full)), i] <- aa_cur_dist
}
########################################################################################################################################################
################################################################################################################################################
###############################################       Prediction       ####################################################################
########################################################################################################################################################
########################################################################################################################################################
# Distance-based prediction of cis (same-chromosome) interactions
# assign tiles to lncRNAs only based on distance and see how many you get right for each lncRNA

sum(!is.na(lncRNA_interaction_distance_matrix))
#314773 total number of data-points
sum(lncRNA_interaction_label_matrix == 1)
# 189989 total number of true interactions --> actually a lot of these are trans ints that we can't predict using genomic-distance-based methods
# excluding trans interactions we have:
sum(!is.na(lncRNA_interaction_distance_matrix) & lncRNA_interaction_label_matrix == 1)
# 84787 true cis interactions (total number of TPs) that we have a chance to predict correctly and
sum(!is.na(lncRNA_interaction_distance_matrix) & lncRNA_interaction_label_matrix == 0)
# 229986 are potential FPs--> these are same chromsosome interactions that are not reported
# Hence for a random method the general 84787/314773 = 0.2693592 is the expected TPR rate

# number of trans interactions:
sum(is.na(lncRNA_interaction_distance_matrix) & lncRNA_interaction_label_matrix == 1)
# 105202
# number of trans interaction for Malat1 only:
sum(is.na(lncRNA_interaction_distance_matrix[,1]) & lncRNA_interaction_label_matrix[, 1] == 1)
# 104553 --> Hence 99% of all trans interactions among these lncRNAs are mediated by Malat1 only
sum(is.na(lncRNA_interaction_distance_matrix[,1]) & lncRNA_interaction_label_matrix[, 1] == 1)/sum(lncRNA_interaction_label_matrix[, 1] == 1)
# 0.9461382 of all Malat1 interactions are trans -- > this is very different from other TFs (close to 1% for most)
aa_nu_zero_dist <- numeric(length = ncol(lncRNA_interaction_distance_matrix))
aa_nu_total_cis_ints <- numeric(length = ncol(lncRNA_interaction_distance_matrix))
aa_nu_total_all_ints <- numeric(length = ncol(lncRNA_interaction_distance_matrix))

for(i in 1:ncol(lncRNA_interaction_distance_matrix)){
  #print(lncRNA_chosen_gt1k_uniqTiles$gene_name[lncRNA_chosen_gt1k_uniqTiles$geneID == colnames(lncRNA_interaction_distance_matrix)[i]])
  #print("Total number of ints:")
  aa_nu_total_all_ints[i] <- sum(lncRNA_interaction_label_matrix[,i])
  #print(sum(lncRNA_interaction_label_matrix[,i]))
  #print("Percent Cis ints:")
  aa_nu_total_cis_ints[i] <- sum((lncRNA_interaction_label_matrix[,i] == 1) & !is.na(lncRNA_interaction_distance_matrix[,i]))
  #print(aa_nu_total_cis_ints[i] / sum(lncRNA_interaction_label_matrix[,i]))
  aa_nu_zero_dist[i] <- sum(lncRNA_interaction_distance_matrix[,i] == 0, na.rm = T)
  
  #print("##########")
  
}

names(aa_nu_total_all_ints) <- lncRNA_chosen_gt1k_uniqTiles$gene_name[match(colnames(lncRNA_interaction_distance_matrix), lncRNA_chosen_gt1k_uniqTiles$geneID)]
names(aa_nu_total_cis_ints) <- names(aa_nu_total_all_ints)
names(aa_nu_zero_dist) <- names(aa_nu_total_all_ints)

par(mfrow = c(1,1), mar = c(8,4,4,4))
barplot(aa_nu_zero_dist, las = 2, 
        main = "total number of 0 distance interactions")
barplot(aa_nu_zero_dist/aa_nu_total_cis_ints, 
        las = 2, main = "% of cis interactions that are 0 distance")


# make a version of the distance matrix that removes the 0 distance ones as
# the probably show nascent transcription

lncRNA_interaction_distance_matrix_0Removed <- lncRNA_interaction_distance_matrix
lncRNA_interaction_distance_matrix_0Removed[lncRNA_interaction_distance_matrix == 0] <- NA

lncRNA_interaction_distance_matrix_full_0Removed <- lncRNA_interaction_distance_matrix_full
lncRNA_interaction_distance_matrix_full_0Removed[lncRNA_interaction_distance_matrix_full == 0] <- NA



#  potential FPs  can be looked at, with varying thresholds
#  for calling interactions for each tile (eg 1, 2, 3 closest lncRNA
#  being associated with each tile, or distance threshold by increaments 
#  of 1MB or some percentage of chromosome's length, Or using thresholds on some combination 
#  of distance and relative expression mesure (like distance (D) * relative expression(E) * beta (some constant factor that can be changed in a range))
# Idea is to use some measure of the expression of that lncRNA to reflect how far it can be found (like amplifier strength) as well).

aa_Expected_TPR_per_lncRNA <-  numeric(length = ncol(lncRNA_interaction_distance_matrix))
for(i in 1:ncol(lncRNA_interaction_distance_matrix)){
  aa_Expected_TPR_per_lncRNA[i] <- sum((lncRNA_interaction_label_matrix[,i] == 1) & !is.na(lncRNA_interaction_distance_matrix_0Removed[,i]))/sum(!is.na(lncRNA_interaction_distance_matrix_0Removed[,i]))
}
names(aa_Expected_TPR_per_lncRNA) <- lncRNA_chosen_gt1k_uniqTiles$gene_name[match(colnames(lncRNA_interaction_distance_matrix), lncRNA_chosen_gt1k_uniqTiles$geneID)]
par(mfrow = c(1,1), mar = c(8,4,4,4))
barplot(aa_Expected_TPR_per_lncRNA, las = 2, 
        main = "expected random TPR\nbased on chromosome assignment\nremoved 0 distances")
abline(h = seq(0.1,1, 0.1), lty = 2, lwd = 0.5, col = "grey")


sum((rowSums(lncRNA_interaction_label_matrix)) > 1)
# 8581 are the number of tiles that interact with more than one lncRNA (among 180220 total tiles that interact with at least one lncRNA)
#sum(((lncRNA_interaction_label_matrix) == 1) & !is.na(lncRNA_interaction_distance_matrix)) - sum(rowSums(!is.na(lncRNA_interaction_distance_matrix))> 0) 
# 9769 this the total number of potential false negatives, due to one-to-one assignment of tiles to lncRNAs -- > wrong
# of course the total number of potential FN is  sum(lncRNA_interaction_label_matrix) which is equal to 189989
# hence only about 5% imptovement in FN avoidance can be potentially achieved through optimizing the k, or distance threshold
# This is in contrast to the 39% of FPs that might be detected by assigning looser thresholds --> 
#  hence it is likely that stricter thresholds lead to better results --> this depends on the similarity of
#  the bound tiles between various lncRNAs residing on the same chromososme.
nrow(lncRNA_interaction_label_matrix)
ncol(lncRNA_interaction_label_matrix)


sum(is.na(lncRNA_interaction_distance_matrix) & (lncRNA_interaction_label_matrix == 1))
# 105202 this is the total number of trans interactions -> which makes 0.5537268 of all interactions


# looking at distance based separation per lncRNA:
par(mfrow = c(7,4), mar = c(2,2,3,2))
for(i in 1:ncol(lncRNA_interaction_distance_matrix)){
  aa_cname <- lncRNA_chosen_gt1k_uniqTiles$gene_name[lncRNA_chosen_gt1k_uniqTiles$geneID == colnames(lncRNA_interaction_distance_matrix)[i]]
  aa_one_nu <- sum(!is.na(lncRNA_interaction_distance_matrix_0Removed[, i]) & lncRNA_interaction_label_matrix[,i] == 1)
  aa_zer_nu <- sum(!is.na(lncRNA_interaction_distance_matrix_0Removed[, i]) & lncRNA_interaction_label_matrix[,i] == 0)
  boxplot(lncRNA_interaction_distance_matrix_0Removed[,i]~lncRNA_interaction_label_matrix[,i],
          main = paste0(aa_cname, "\n#1: ",aa_one_nu,  " #0: ", aa_zer_nu))
}
par(mfrow = c(7,4), mar = c(2,2,3,2))
for(i in 1:ncol(lncRNA_interaction_distance_matrix)){
  aa_cname <- lncRNA_chosen_gt1k_uniqTiles$gene_name[lncRNA_chosen_gt1k_uniqTiles$geneID == colnames(lncRNA_interaction_distance_matrix)[i]]
  aa_one_nu <- sum(!is.na(lncRNA_interaction_distance_matrix_0Removed[, i]) & lncRNA_interaction_label_matrix[,i] == 1)
  aa_zer_nu <- sum(!is.na(lncRNA_interaction_distance_matrix_0Removed[, i]) & lncRNA_interaction_label_matrix[,i] == 0)
  boxplot(log10(lncRNA_interaction_distance_matrix_0Removed[,i])~lncRNA_interaction_label_matrix[,i],
          main = paste0("log10 ", aa_cname, "\n#1: ",aa_one_nu,  " #0: ", aa_zer_nu), ylim= c(2,8.5))
  abline (h=seq(0,10,1), lty = 2, lwd = 0.7, col = "grey")
}

########################################################################################################################################################################
# create a cis label matrix and a cis distance matrix --> remove extra tiles that are only trans
aachn <- paste0("chr" ,lncRNA_chosen_gt1k_uniqTiles$chromosome_name)
#aa_tilenames_split_3
lncRNA_interaction_label_matrix_cis <- lncRNA_interaction_label_matrix
lncRNA_interaction_distance_matrix_0Removed_cis <- lncRNA_interaction_distance_matrix_0Removed
for(i in 1:ncol(lncRNA_interaction_label_matrix_cis)){
  aa_cur_chr <- aachn[lncRNA_chosen_gt1k_uniqTiles$geneID == colnames(lncRNA_interaction_label_matrix_cis)[i]]
  print("nu_trans of")
  print(lncRNA_chosen_gt1k_uniqTiles$gene_name[lncRNA_chosen_gt1k_uniqTiles$geneID == colnames(lncRNA_interaction_label_matrix_cis)[i]])
  print(sum(lncRNA_interaction_label_matrix_cis[!(aa_tilenames_split_3 %in% aa_cur_chr),i]))
  lncRNA_interaction_label_matrix_cis[!(aa_tilenames_split_3 %in% aa_cur_chr),i] <- 0
}
aaz <- which(rowSums(lncRNA_interaction_label_matrix_cis == 1) == 0)
lncRNA_interaction_label_matrix_cis <- lncRNA_interaction_label_matrix_cis[-aaz,]
lncRNA_interaction_distance_matrix_0Removed_cis <- lncRNA_interaction_distance_matrix_0Removed_cis[-aaz,]



#######################################################################################################################
#######################################################################################################################

aacmp$General

#Sensitivity (TPR) : TP/P = TP/(TP + FN)
#Specificity (TNR) : TN/N = TN/(TN + FP)
#Precision (PPV)   : TP/(TP + FP)
#False Disc Rate (FDR): FP/(TP + FP)
#accuracy: (TP + TN)/(P + N) = (TP + TN)/(TP + FN + TN + FP)
# balanced accuracy: (TPR + TNR)/2
# F1 score : is the harmonic mean of precision and sensitivity = 2/((1/Prec) + (1/Sens)) = (2*TP)/(2*TP + FN + FP)

aatsttt <- DisChaLDis(DisMat = lncRNA_interaction_distance_matrix_0Removed_cis, 
                      k_closest = 1)
aacmp2 <- IntMatEval(real_intMat = lncRNA_interaction_label_matrix,
                     pred_intMat = aatsttt,
                     DisMat = lncRNA_interaction_distance_matrix_0Removed)
aacmp3 <- IntMatEval(real_intMat = lncRNA_interaction_label_matrix_cis,
                     pred_intMat = aatsttt,
                     DisMat = lncRNA_interaction_distance_matrix_0Removed_cis)


aacmp2
aatsttt


aaareal_intMat_Cis <- lncRNA_interaction_label_matrix
aaareal_intMat_Cis[aaareal_intMat_Cis == 1 & is.na(lncRNA_interaction_distance_matrix_0Removed)] <- 0

IntMatEval_plot(eval_out = aacmp2, my_filename = "~/Documents/Shayan/BioInf/lncRNA/plots/eval_k1.png")
IntMatEval_plot(eval_out = aacmp3, my_filename = "~/Documents/Shayan/BioInf/lncRNA/plots/eval_k1_cisOnly3.png")


aa_res <- list()
aa_cmp <- list()
aa_dis <- c(1e3, 1e4, 1e5, 1e6, 1e7, 1e8)
for(i in 1:length(aa_dis)){
  print(i)
  aa_res[[i]] <-  DisChaLDis(DisMat = lncRNA_interaction_distance_matrix_0Removed_cis,
                             dist_thresh = aa_dis[i])
  aa_cmp[[i]] <- IntMatEval(real_intMat = lncRNA_interaction_label_matrix_cis,
                            pred_intMat = aa_res[[i]],
                            DisMat = lncRNA_interaction_distance_matrix_0Removed_cis)
  IntMatEval_plot(eval_out = aa_cmp[[i]], my_filename = paste0("~/Documents/Shayan/BioInf/lncRNA/plots/eval_dist_cisOnly_",aa_dis[i], ".png"))
}

aa_res2 <- list()
aa_cmp2 <- list()
aa_k <- c(1, 2, 3, 4)
for(i in 1:length(aa_k)){
  print(i)
  aa_res2[[i]] <-  DisChaLDis(DisMat = lncRNA_interaction_distance_matrix_0Removed_cis,
                             k_closest = aa_k[i])
  aa_cmp2[[i]] <- IntMatEval(real_intMat = lncRNA_interaction_label_matrix_cis,
                            pred_intMat = aa_res2[[i]],
                            DisMat = lncRNA_interaction_distance_matrix_0Removed_cis)
  IntMatEval_plot(eval_out = aa_cmp2[[i]], my_filename = paste0("~/Documents/Shayan/BioInf/lncRNA/plots/eval_cisOnly_k_",aa_k[i], ".png"))
}
#######################################################################################################################

# create venn diagram per chromosome of bound vs unbound tiles
# create venn diagram per chromosome of bound by different lncRNAs in cis, trans and both settings



######################################################################################################################
######################################################################################################################
######################################################################################################################
# Creating positive and negative sets for cis lncRNA interacting vs non-interacting regions.
# I will first look at the overlap between the positive class and available accessibility data
# Create a GRanges object for Cis-interacting regions:

rownames(lncRNA_interaction_label_matrix_cis)[1:5]
Cis_RNA28_interacting_positive <- MM9_1kb_tiled_GR[match(rownames(lncRNA_interaction_label_matrix_cis), MM9_1kb_tiled_GR$tile)]

options(scipen=999)
write.table(Cis_RNA28_interacting_positive, 
            file="~/Documents/Shayan/BioInf/lncRNA/Data_sets/Cis_RNA_interaction/Cis_RNA28_interacting_positive.bed",
            quote=F, sep="\t", row.names=F, col.names=F)
options(scipen=0)





############################
# read Cage nuclear experession # this is mm9
aa_cage_nuclear <- read.table("~/Documents/Shayan/BioInf/lncRNA/NG-A36034R2-Carninci_Supplementary_DataSets2/Carninci_Mm1_CAGE_nuc.txt",
                              stringsAsFactors = F)
aanaame <- strsplit("Chr	Start	End	ID	.	Strand	miPSF_Nucl	miPSB_Nucl	miPST_Nucl	mESR08_Nucl	mESFVB1_Nucl	mESB6G2_Nucl	MEF_Nucl	mLymB_Nucl	mLymT_Nucl	logConc	logFC	P.Value	FDR", "\t")
colnames(aa_cage_nuclear) <- unlist(aanaame)

View(head(aa_cage_nuclear))
aa_cage_nuclear2 <- aa_cage_nuclear[, c(c(1:4), 6, c(10:12), c(16:19))]

sum(duplicated(aa_cage_nuclear2$ID))

par(mfrow = c(1,1), mar = c(4,4,4,4))
boxplot.matrix(log10(as.matrix(aa_cage_nuclear2[, c(6:8)]) + 1))

aa_cage_nuclear2_num <- aa_cage_nuclear2[c("mESR08_Nucl", "mESFVB1_Nucl", "mESB6G2_Nucl")]
rownames(aa_cage_nuclear2_num) <- aa_cage_nuclear2$ID
aa_max <- apply(aa_cage_nuclear2_num, 1, max)
aa_cage_nuclear3 <- cbind(aa_cage_nuclear2[, c(1,2,3,4,5)], aa_max)
aa_cage_nuclear3 <- aa_cage_nuclear3[aa_cage_nuclear3$aa_max >= 1,]
mESC_CAGE_mm9 <- aa_cage_nuclear3
mESC_CAGE_mm9_GR <- makeGRangesFromDataFrame(mESC_CAGE_mm9, keep.extra.columns = T)


aaov <- findOverlaps(MM9_1kb_tiled_GR, mESC_CAGE_mm9_GR)
aadf <- data.frame(tile_index = aaov@from,
                   atac_index = aaov@to, 
                   CAGE_score = mESC_CAGE_mm9_GR$aa_max[aaov@to])
aadf_agg <- aggregate(aadf[c("CAGE_score")],
                      by = aadf[c("tile_index")],
                      FUN = max)

mESC_CAGE_mm9_tiled <- cbind(MM9_1kb_tiled[aadf_agg$tile_index,], aadf_agg$CAGE_score)
mESC_CAGE_mm9_tiledGR <- makeGRangesFromDataFrame(mESC_CAGE_mm9_tiled, keep.extra.columns = T)


vennplot(Sets = list(RNA_interacting = mESC_cis_interacting_mm9_tiled,
                     ATAC_seq = mESC_ATAC_Seq_tiledGR,
                     CAGE = mESC_CAGE_mm9_tiledGR
                     #,eRNA_intersect = li_erna_intersect_gr
                     #,eRNA_intersect_minus_etoh = li_erna_intersect_minus_etoh_intersect_gr
),
by = "Vennerable")

####################################################################################
# read mESC ATAC-seq # this is converted mm9
mESC_ATAC_Seq <- read.table(file = "~/Documents/Shayan/BioInf/lncRNA/mESC_ATACseq/GSM3109355_HH_229_peaks_5cols_mm9.bed",
                            stringsAsFactors = F)
colnames(mESC_ATAC_Seq) <- c("chr", "start", "end", "name", "score")
hist(log10(mESC_ATAC_Seq$score))
mESC_ATAC_Seq_GR <- makeGRangesFromDataFrame(df = mESC_ATAC_Seq, keep.extra.columns = T)

# get which tiles overlap with ATAC-seq data
aaov <- findOverlaps(MM9_1kb_tiled_GR, mESC_ATAC_Seq_GR)
aadf <- data.frame(tile_index = aaov@from,
                   atac_index = aaov@to, 
                   atac_score = mESC_ATAC_Seq_GR$score[aaov@to])
aadf_agg <- aggregate(aadf[c("atac_score")],
                      by = aadf[c("tile_index")],
                      FUN = max)
#MM9_1kb_tiled$end <- MM9_1kb_tiled$end - 1
mESC_ATAC_Seq_tiled <- cbind(MM9_1kb_tiled[aadf_agg$tile_index,], aadf_agg$atac_score)
mESC_ATAC_Seq_tiledGR <- makeGRangesFromDataFrame(mESC_ATAC_Seq_tiled, keep.extra.columns = T)

mESC_cis_interacting_mm9_tiled <- MM9_1kb_tiled_GR[MM9_1kb_tiled_GR$tile %in% rownames(lncRNA_interaction_label_matrix_cis)]

vennplot(Sets = list(RNA_interacting = mESC_cis_interacting_mm9_tiled,
                     ATAC_seq = mESC_ATAC_Seq_tiledGR
                     #eRNA_union = li_erna_union_gr
                     #,eRNA_intersect = li_erna_intersect_gr
                     #,eRNA_intersect_minus_etoh = li_erna_intersect_minus_etoh_intersect_gr
),
by = "Vennerable")
####################################################################################
# read data from encode chromatin marks 
####################################################################################
# read  chromatin mark data from encode 
# write a script to:
# 1. unzip the files that are from both replicates
# 2. convert them to mm9 using liftover and name appropriately
# 3. read them in and convert to mm9 tile format
# 4. look at intersection with the other features

aa_meta <- read.delim("~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/Downloaded_from_Encode/metadata.tsv",
                      header = T,
                      stringsAsFactors = F,
                      strip.white = T)

aa_fname <- aa_meta$File.accession[aa_meta$Biological.replicate.s. == "1, 2"]

aafiles_full <- paste0("~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/Downloaded_from_Encode/",aa_fname, ".bed.gz" )
aafiles_full_2 <- paste0("~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/Downloaded_from_Encode/",aa_fname, ".bed" )
system("mkdir ~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/Downloaded_from_Encode/processed_mm9/")
cat(c("#!/bin/bash\n"), file = "~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/Downloaded_from_Encode/getChromatinMarks.sh",
    append = F)
for(i in 1:length(aafiles_full)){
  cat(c("gunzip ", aafiles_full[i], "\n"), 
      file = "~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/Downloaded_from_Encode/getChromatinMarks.sh",
      sep = " ",
      append = T)
  cat(c("~/Documents/Shayan/BioInf/liftOver -bedPlus=6", 
        aafiles_full_2[i],
        " ~/Documents/Shayan/BioInf/Liftover_chain_files/mm10ToMm9.over.chain ", 
        paste0("~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/Downloaded_from_Encode/processed_mm9/",
               aa_fname[i],
               "_",
               aa_meta$Experiment.target[aa_meta$File.accession == aa_fname[i]], 
               "_",
               aa_meta$Biosample.term.name[aa_meta$File.accession == aa_fname[i]],
               ".bed"),
        paste0("aa_unconverted_", i),
        "\n"),
      file = "~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/Downloaded_from_Encode/getChromatinMarks.sh",
      sep = " ", append = T)
  cat(c("gzip ", aafiles_full_2[i], "\n"),
      file = "~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/Downloaded_from_Encode/getChromatinMarks.sh",
      sep = " ", append = T)
}

# read them in, and tile them
aafiles_full <- list.files("~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/Downloaded_from_Encode/processed_mm9/", 
                      full.names = T)
aafiles <- list.files("~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/Downloaded_from_Encode/processed_mm9/")
aa_chrom <- list()
for(i in 1:length(aafiles)){
  aa_chrom[[i]] <- read.table(aafiles_full[i], stringsAsFactors = F)
}
names(aa_chrom) <- aafiles
aasp <- strsplit(aafiles, "_")
aasp_mark <- unlist(lapply(aasp, "[[", 2))
aasp_mark_unq <- unique(aasp_mark)
aa_chrom_combind <- list()
for(i in 1:length(aasp_mark_unq)){
  aa_chrom_combind[[i]] <- do.call(rbind, aa_chrom[aasp_mark %in% aasp_mark_unq[i]])
}
names(aa_chrom_combind) <- aasp_mark_unq
head(aa_chrom_combind$`H3K27me3-mouse`)
mESC_Chromatin_mark_list <- aa_chrom_combind
mESC_Chromatin_mark_GR_list <- lapply(mESC_Chromatin_mark_list,
                                      makeGRangesFromDataFrame,
                                      keep.extra.columns = T,
                                      seqnames.field = "V1",
                                      start.field = "V2",
                                      end.field = "V3")

mESC_Chromatin_mark_tiled_list <- list()
mESC_Chromatin_mark_tiled_list_GR <- list()

for(i in 1:length(mESC_Chromatin_mark_GR_list)){
  aaov <- findOverlaps(MM9_1kb_tiled_GR, mESC_Chromatin_mark_GR_list[[i]])
  aadf <- data.frame(tile_index = aaov@from,
                     chrom_index = aaov@to, 
                     chrom_score = mESC_Chromatin_mark_GR_list[[i]]$V5[aaov@to])
  aadf_agg <- aggregate(aadf[c("chrom_score")],
                        by = aadf[c("tile_index")],
                        FUN = max)
  mESC_Chromatin_mark_tiled_list[[i]] <- cbind(MM9_1kb_tiled[aadf_agg$tile_index,], aadf_agg$chrom_score)
  mESC_Chromatin_mark_tiled_list_GR[[i]] <- makeGRangesFromDataFrame(mESC_Chromatin_mark_tiled_list[[i]], keep.extra.columns = T)
  
}
names(mESC_Chromatin_mark_tiled_list) <- names(mESC_Chromatin_mark_GR_list)
names(mESC_Chromatin_mark_tiled_list_GR) <- names(mESC_Chromatin_mark_GR_list)







vennplot(Sets = list(RNA_interacting = mESC_cis_interacting_mm9_tiled,
                     ATAC_seq = mESC_ATAC_Seq_tiledGR,
                     H3K27ac=mESC_Chromatin_mark_tiled_list_GR$`H3K27ac-mouse`
                     #eRNA_union = li_erna_union_gr
                     #,eRNA_intersect = li_erna_intersect_gr
                     #,eRNA_intersect_minus_etoh = li_erna_intersect_minus_etoh_intersect_gr
),
by = "Vennerable")

vennplot(Sets = list(RNA_interacting = mESC_cis_interacting_mm9_tiled,
                     H3K27ac=mESC_Chromatin_mark_tiled_list_GR$`H3K27ac-mouse`,
                     H3K27me3=mESC_Chromatin_mark_tiled_list_GR$`H3K27me3-mouse`
                     #eRNA_union = li_erna_union_gr
                     #,eRNA_intersect = li_erna_intersect_gr
                     #,eRNA_intersect_minus_etoh = li_erna_intersect_minus_etoh_intersect_gr
),
by = "Vennerable")

vennplot(Sets = list(RNA_interacting = mESC_cis_interacting_mm9_tiled,
                     H3K4me1=mESC_Chromatin_mark_tiled_list_GR$`H3K9ac-mouse`#,
                     #H3K27me3=mESC_Chromatin_mark_tiled_list_GR$`H3K27me3-mouse`
                     #eRNA_union = li_erna_union_gr
                     #,eRNA_intersect = li_erna_intersect_gr
                     #,eRNA_intersect_minus_etoh = li_erna_intersect_minus_etoh_intersect_gr
),
by = "Vennerable")

##################################################################################################
# read data from Ren Lab

aafiles_full <- list.files("~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/RenLab", pattern = "*.bigwig", full.names = T )
aafiles_full_wig <- paste0(unlist(lapply(strsplit(aafiles_full, "\\."), "[[", 1)), ".wig")
aafiles_full_bed <- paste0(unlist(lapply(strsplit(aafiles_full, "\\."), "[[", 1)), ".bed")

aafiles <- list.files("~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/RenLab/", pattern = "*.bigwig", full.names = F )

cat(c("#!/bin/bash\n"),
    file = "~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/RenLab/convert_bigwig_to_bed.sh",
    append = F)

for(i in 1:length(aafiles_full)){
  cat(c("~/Documents/Shayan/BioInf/bigWigToWig ", aafiles_full[i], aafiles_full_wig[i], "\n"), 
      file = "~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/RenLab/convert_bigwig_to_bed.sh",
      sep = " ",
      append = T)
  cat(c("wig2bed < ", 
        aafiles_full_wig[i], " > ",
        aafiles_full_bed[i],
        "\n"),
      file = "~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/RenLab/convert_bigwig_to_bed.sh",
      sep = " ", append = T)
  cat(c("rm ", aafiles_full_wig[i], "\n"),
      file = "~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/RenLab/convert_bigwig_to_bed.sh",
      sep = " ", append = T)
}

aafiles_full <- list.files("~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/RenLab",
                           pattern = "*.bed",
                           full.names = T)[2:7]
aafiles <- list.files("~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/RenLab",
                      pattern = "*.bed")[2:7]
aa_Renlab <- list()
for(i in 1:length(aafiles)){
  print(aafiles[i])
  aa_Renlab[[i]] <- read.table(aafiles_full[i], stringsAsFactors = F)
}
names(aa_Renlab) <- aafiles

mESC_Renlab_list <- aa_Renlab
rm(aa_Renlab)
mESC_Renlab_list_FILTERED <- mESC_Renlab_list
# filter the bed files for significant ones
mESC_Renlab_list_FILTERED$`GSM723015_RenLab-CTCF-mESC.bed` <- mESC_Renlab_list$`GSM723015_RenLab-CTCF-mESC.bed`[mESC_Renlab_list$`GSM723015_RenLab-CTCF-mESC.bed`$V5 > -log10(0.05),]
mESC_Renlab_list_FILTERED$`GSM723016_RenLab-H3K4me1-mESC.bed` <- mESC_Renlab_list$`GSM723016_RenLab-H3K4me1-mESC.bed`[mESC_Renlab_list$`GSM723016_RenLab-H3K4me1-mESC.bed`$V5 > -log10(0.05),]
mESC_Renlab_list_FILTERED$`GSM723017_RenLab-H3K4me3-mESC.bed` <- mESC_Renlab_list$`GSM723017_RenLab-H3K4me3-mESC.bed`[mESC_Renlab_list$`GSM723017_RenLab-H3K4me3-mESC.bed`$V5 > -log10(0.05),]
mESC_Renlab_list_FILTERED$`GSM723018_RenLab-P300-mESC.bed` <- mESC_Renlab_list$`GSM723018_RenLab-P300-mESC.bed`[mESC_Renlab_list$`GSM723018_RenLab-P300-mESC.bed`$V5 > -log10(0.05),]
mESC_Renlab_list_FILTERED$`GSM723019_RenLab-Pol2-mESC.bed` <- mESC_Renlab_list$`GSM723019_RenLab-Pol2-mESC.bed`[mESC_Renlab_list$`GSM723019_RenLab-Pol2-mESC.bed`$V5 > -log10(0.05),]
mESC_Renlab_list_FILTERED$`GSM723776_RenLab-RNA-Seq-mESC-ZY28.bed` <- mESC_Renlab_list$`GSM723776_RenLab-RNA-Seq-mESC-ZY28.bed`[mESC_Renlab_list$`GSM723776_RenLab-RNA-Seq-mESC-ZY28.bed`$V5 >= 5,]
summary(mESC_Renlab_list$`GSM723015_RenLab-CTCF-mESC.bed`$V5)
lapply(mESC_Renlab_list_FILTERED, nrow)


mESC_Renlab_GR_list <- lapply(mESC_Renlab_list_FILTERED,
                              makeGRangesFromDataFrame,
                              keep.extra.columns = T,
                              seqnames.field = "V1",
                              start.field = "V2",
                              end.field = "V3")

mESC_Renlab_tiled_list <- list()
mESC_Renlab_tiled_list_GR <- list()

for(i in 1:length(mESC_Renlab_GR_list)){
  print(i)
  aaov <- findOverlaps(MM9_1kb_tiled_GR, mESC_Renlab_GR_list[[i]])
  aadf <- data.frame(tile_index = aaov@from,
                     RL_index = aaov@to, 
                     score = mESC_Renlab_GR_list[[i]]$V5[aaov@to])
  if(i == 6){ # use sum for RNA seq and max for others
    aadf_agg <- aggregate(aadf[c("score")],
                          by = aadf[c("tile_index")],
                          FUN = sum)
  }else{
    aadf_agg <- aggregate(aadf[c("score")],
                          by = aadf[c("tile_index")],
                          FUN = max)
  }

  mESC_Renlab_tiled_list[[i]] <- cbind(MM9_1kb_tiled[aadf_agg$tile_index,], aadf_agg$score)
  mESC_Renlab_tiled_list_GR[[i]] <- makeGRangesFromDataFrame(mESC_Renlab_tiled_list[[i]],
                                                             keep.extra.columns = T)
  
}
names(mESC_Renlab_tiled_list) <- names(mESC_Renlab_GR_list)
names(mESC_Renlab_tiled_list_GR) <- names(mESC_Renlab_GR_list)



vennplot(Sets = list(RNA_interacting = mESC_cis_interacting_mm9_tiled,
                     ATAC_seq = mESC_ATAC_Seq_tiledGR,
                     RNA_Seq_gt10=mESC_Renlab_tiled_list_GR$`GSM723776_RenLab-RNA-Seq-mESC-ZY28.bed`[mESC_Renlab_tiled_list_GR$`GSM723776_RenLab-RNA-Seq-mESC-ZY28.bed`$`aadf_agg$score` > 10]
                     #eRNA_union = li_erna_union_gr
                     #,eRNA_intersect = li_erna_intersect_gr
                     #,eRNA_intersect_minus_etoh = li_erna_intersect_minus_etoh_intersect_gr
),
by = "Vennerable")
####################################################################################
# read DNAase data

aa_dnase1 <- read.table("~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/DNAse/GSM1014154_mm9_wgEncodeUwDnaseEse14129olaME0HotspotsRep1.broadPeak", stringsAsFactors = F)
aa_dnase2 <- read.table("~/Documents/Shayan/BioInf/lncRNA/mESC_encode_chip/DNAse/GSM1014154_mm9_wgEncodeUwDnaseEse14129olaME0HotspotsRep2.broadPeak", stringsAsFactors = F)
mESC_DNAse <- rbind(aa_dnase1, aa_dnase2)
mESC_DNAse <- mESC_DNAse[, c(1,2,3,7)]
mESC_DNAse_GR <- makeGRangesFromDataFrame(mESC_DNAse,
                                          keep.extra.columns = T,
                                          seqnames.field = "V1",
                                          start.field = "V2",
                                          end.field = "V3" )

aaov <- findOverlaps(MM9_1kb_tiled_GR, mESC_DNAse_GR)
aadf <- data.frame(tile_index = aaov@from,
                   dnase_index = aaov@to, 
                   dnase_score = mESC_DNAse_GR$V7[aaov@to])
aadf_agg <- aggregate(aadf[c("dnase_score")],
                      by = aadf[c("tile_index")],
                      FUN = max)

mESC_DNAse_mm9_tiled <- cbind(MM9_1kb_tiled[aadf_agg$tile_index,], aadf_agg$dnase_score)
mESC_DNAse_mm9_tiledGR <- makeGRangesFromDataFrame(mESC_DNAse_mm9_tiled, keep.extra.columns = T)


vennplot(Sets = list(RNA_interacting = mESC_cis_interacting_mm9_tiled,
                     ATAC_seq = mESC_ATAC_Seq_tiledGR,
                     DNAse = mESC_DNAse_mm9_tiledGR
                     #,eRNA_intersect = li_erna_intersect_gr
                     #,eRNA_intersect_minus_etoh = li_erna_intersect_minus_etoh_intersect_gr
),
by = "Vennerable")
################################################################################################
# read repeats
MM9_repeats <- read.table("~/Documents/Shayan/BioInf/lncRNA/Repeat_elements/mm9.fa.out", skip = 2)
aa_repeat_types <- unique(MM9_repeats$V11)
aa_repeat_types_modified <- gsub(pattern = "\\?", replacement = "", x = aa_repeat_types)
aa_repeat_types_modified <- gsub(pattern = "\\/", replacement = "__", x = aa_repeat_types_modified)
aa_repeat_types_modified <- gsub(pattern = "-", replacement = "_", x = aa_repeat_types_modified)

MM9_repeats_list <- list()
for(i in 1:length(aa_repeat_types)){
  print(i)
  MM9_repeats_list[[i]] <- MM9_repeats[MM9_repeats$V11 %in% aa_repeat_types[i],]
  MM9_repeats_list[[i]]$V11 <- gsub(pattern = "\\?", replacement = "", x = MM9_repeats_list[[i]]$V11)
  MM9_repeats_list[[i]]$V11 <- gsub(pattern = "\\/", replacement = "__", x = MM9_repeats_list[[i]]$V11)
  MM9_repeats_list[[i]]$V11 <- gsub(pattern = "-", replacement = "_", x = MM9_repeats_list[[i]]$V11)
}
names(MM9_repeats_list) <- aa_repeat_types_modified

MM9_repeats_list_GR <- lapply(MM9_repeats_list, makeGRangesFromDataFrame, seqnames.field = "V5",
                              start.field = "V6", end.field = "V7", keep.extra.columns = T)
MM9_repeats_list_tiled <- list()
load("~/Documents/Shayan/BioInf/lncRNA/MM9_1kb_tiled_GR_filtered.RData")
for(i in 1:length(MM9_repeats_list_GR)){
  print(i)
  aaov <- findOverlaps(MM9_1kb_tiled_GR_filtered, MM9_repeats_list_GR[[i]])
  aadf <- data.frame(tile_index = aaov@from,
                     repeat_index = aaov@to, 
                     repeat_class = MM9_repeats_list_GR[[i]]$V11[aaov@to])
  aadf_agg <- aggregate(aadf[c("repeat_class")],
                        by = aadf[c("tile_index")],
                        FUN = length)
  
  MM9_repeats_list_tiled[[i]] <- data.frame(tile_name = MM9_1kb_tiled_GR_filtered$tile[aadf_agg$tile_index],
                                            nu_hits = aadf_agg$repeat_class)
}
names(MM9_repeats_list_tiled) <- aa_repeat_types_modified
save(list = c("MM9_repeats_list_tiled"), file = "~/Documents/Shayan/BioInf/lncRNA/MM9_repeats_list_tiled.RData")
# mESC_repeat_mm9_tiled <- cbind(MM9_1kb_tiled[aadf_agg$tile_index,], aal)
# mESC_repeat_mm9_tiledGR <- makeGRangesFromDataFrame(mESC_repeat_mm9_tiled, keep.extra.columns = T)
# 
# 
# vennplot(Sets = list(RNA_interacting = mESC_cis_interacting_mm9_tiled,
#                      repeats = mESC_repeat_mm9_tiledGR,
#                      DNAse = mESC_DNAse_mm9_tiledGR
#                      #,eRNA_intersect = li_erna_intersect_gr
#                      #,eRNA_intersect_minus_etoh = li_erna_intersect_minus_etoh_intersect_gr
# ),
# by = "Vennerable")

####################################################################################################
####################################################################################################
# write a function that assigns TADs to tiles and lncRNAs
MM9_1kb_tiled_GR
mm9_TAD_list <- list()
aafiles <- list.files("~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/mES", full.names = T, pattern = ".txt")

# sort bed files
cat(c("#!/bin/bash\n"),
    file = "~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/mES/sortBedfiles.sh", 
    append = F)
for(i in 1:length(aafiles)){
  cat(c("cat ", aafiles[i], " | tr ' ' '\\t' > ", paste0(aafiles[i], "_tabdel\n")),
      sep = " ", 
      file = "~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/mES/sortBedfiles.sh", 
      append = T)
  cat(c("sortBed -i ", paste0(aafiles[i], "_tabdel"), " > ", paste0(aafiles[i], "_sorted\n")), 
      sep = " ", 
      file = "~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/mES/sortBedfiles.sh", 
      append = T)
}

aafiles2 <- list.files("~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/mES", full.names = F, pattern = ".txt$")
aaf <- unlist(lapply(strsplit(aafiles2, "\\."), "[[", 1))

for(i in 1:length(aafiles)){
  mm9_TAD_list[[i]] <- read.table(paste0(aafiles[i],
                                         "_sorted"),
                                  stringsAsFactors = F)
  mm9_TAD_list[[i]]$V3 <- mm9_TAD_list[[i]]$V3 - 1
  aa_chr_un <- unique(mm9_TAD_list[[i]]$V1)
  aa_number_all <- numeric(length = nrow(mm9_TAD_list[[i]]))
  aa_name_all <- character(length = nrow(mm9_TAD_list[[i]]))
  for(j in 1:length(aa_chr_un)){
    aaw <- which(mm9_TAD_list[[i]]$V1 %in% aa_chr_un[j])
    aa_name_all[aaw] <- paste(aa_chr_un[j], c(1:length(aaw)),sep = "__")
    aa_number_all[aaw] <- c(1:length(aaw))
  }
  mm9_TAD_list[[i]]$TAD_number <- aa_number_all
  mm9_TAD_list[[i]]$TAD_name <- aa_name_all
}
names(mm9_TAD_list) <- aaf
mm9_TAD_list_GR <- lapply(mm9_TAD_list,makeGRangesFromDataFrame,keep.extra.columns = T, 
                          seqnames.field = "V1", start.field = "V2", end.field = "V3")
for(i in 1:length(mm9_TAD_list_GR)){
  write.table(x = mm9_TAD_list_GR[[i]],
              file = paste0("~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/mES/TAD_", names(mm9_TAD_list_GR)[i], ".bed"),
              quote = F, row.names = F, col.names = F, sep = "\t")
  
}


lncRNA_chosen_gt1k_uniqTiles 
write.table(x = MM9_1kb_tiled_GR,
            file = "~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/mES/MM9_1kb_tile.bed",
            quote = F, row.names = F, col.names = F, sep = "\t")
lncRNA_chosen_gt1k_uniqTiles$chromosome_name <- paste0("chr", lncRNA_chosen_gt1k_uniqTiles$chromosome_name)
lncRNA_chosen_gt1k_uniqTiles_GR <- makeGRangesFromDataFrame(lncRNA_chosen_gt1k_uniqTiles[, c(1:5)],
                                                            keep.extra.columns = T, 
                                                            seqnames.field = "chromosome_name", 
                                                            start.field = "start_position", 
                                                            end.field = "end_position")
write.table(x = lncRNA_chosen_gt1k_uniqTiles_GR, file = "~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/mES/MM9_lncRNA_28.bed",
            quote = F, row.names = F, col.names = F, sep = "\t")

aa_tad_bed <- list.files("~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/mES", pattern = "bed$")
aa_tad_bed <- aa_tad_bed[-c(1,2)]

cat("#!/bin/bash\n",
    file = "~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/mES/get_closest_TAD.sh", 
    append = F)
for(i in 1:length(aa_tad_bed)){
  cat(c("bedtools closest -t all -a MM9_1kb_tile_sorted.bed -b ", aa_tad_bed[i], " | cut -f 1,6,12,13 > ", paste0("MM9_1kb_tile_", aa_tad_bed[i], ".closest\n")),
      file = "~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/mES/get_closest_TAD.sh",
      sep = " ",
      append = T)
  cat(c("bedtools closest -t all -a MM9_lncRNA_28_sorted.bed -b ", aa_tad_bed[i], " > ", paste0("MM9_lncRNA_28_", aa_tad_bed[i], ".closest\n")),
      file = "~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/mES/get_closest_TAD.sh",
      sep = " ",
      append = T)
}

# gather the TAD number from each dataset in one dataframe, for tiles and lncRNAs
aa_tile_closest_files <- list.files(path = "~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/mES",
                                    pattern = "MM9_1kb_tile_TAD_*")
aa_tile_tad_number <- matrix(nrow = length(MM9_1kb_tiled_GR),
                             ncol = length(aa_tile_closest_files))
rownames(aa_tile_tad_number) <- MM9_1kb_tiled_GR$tile
colnames(aa_tile_tad_number) <- aa_tile_closest_files
for(i in 1:length(aa_tile_closest_files)){
  print(i)
  aa_tab <- read.table(file = paste0("~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/mES/", aa_tile_closest_files[i]), 
                       stringsAsFactors = F)
  aa_tile_tad_number[,i] <- aa_tab$V3[match(rownames(aa_tile_tad_number), aa_tab$V2)]
}

MM9_1kb_tiled_TADnu <- data.frame(aa_tile_tad_number, stringsAsFactors = F)
for(i in 1:ncol(MM9_1kb_tiled_TADnu)){
  MM9_1kb_tiled_TADnu[, i] <- as.numeric(MM9_1kb_tiled_TADnu[, i] )
}
#MM9_1kb_tiled_TADnu <- as.matrix(MM9_1kb_tiled_TADnu)

aa_lncr_closest_files <- list.files(path = "~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/mES",
                                    pattern = "MM9_lncRNA_28_TAD_*")

aa_lnc28_tad_number <- matrix(nrow = length(lncRNA_chosen_gt1k_uniqTiles_GR),
                             ncol = length(aa_lncr_closest_files))
rownames(aa_lnc28_tad_number) <- lncRNA_chosen_gt1k_uniqTiles_GR$gene_name
colnames(aa_lnc28_tad_number) <- aa_lncr_closest_files
for(i in 1:length(aa_lncr_closest_files)){
  print(i)
  aa_tab <- read.table(file = paste0("~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/mES/", aa_lncr_closest_files[i]), 
                       stringsAsFactors = F)
  aa_lnc28_tad_number[,i] <- as.numeric(aa_tab$V13[match(rownames(aa_lnc28_tad_number), aa_tab$V7)])
}
lncRNA_chosen_gt1k_TADnu <- aa_lnc28_tad_number



aatst <- distance_lnc_tile_TAD(lncRNA_TAD_Assgined = lncRNA_chosen_gt1k_TADnu, 
                               Tile_TAD_Assigned = MM9_1kb_tiled_TADnu, 
                               lncRNA_metaData = lncRNA_chosen_gt1k_uniqTiles)

TAD_based_distance_list <- aatst
names(TAD_based_distance_list) <- colnames(MM9_1kb_tiled_TADnu)
TAD_based_territory_assign_list <- list()

save(list = c("DisChaLDis", "TAD_based_distance_list"), file = "Territory_assign_input.RData")


# for(i in 1:length(TAD_based_distance_list)){
#   print(i)
#   TAD_based_territory_assign_list[[i]]  <- DisChaLDis(DisMat = TAD_based_distance_list[[i]], 
#                                                                 k_closest = 1)
# }
# names(TAD_based_territory_assign_list) <- names(TAD_based_distance_list)
# save(list = c("TAD_based_territory_assign_list"), file = "Territory_assign_output.RData")
load("~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/TAD_based_territory_assign_list.RData")




#create barplots for all FN counts of all 18 TAD based together with k=1 

# first evaluate the 18 TAD based

aa_nearest_1_Pred <- DisChaLDis(DisMat = lncRNA_interaction_distance_matrix_0Removed_cis, 
                                k_closest = 1)
aa_nearest_1_eval <- IntMatEval(real_intMat = lncRNA_interaction_label_matrix_cis,
                                pred_intMat = aa_nearest_1_Pred,
                                DisMat = lncRNA_interaction_distance_matrix_0Removed_cis)

aa_tad_based_eval <- list()
for(i in 1:length(TAD_based_territory_assign_list)){
  print(i)
  aatst <- TAD_based_territory_assign_list[[i]][match(rownames(lncRNA_interaction_distance_matrix_0Removed_cis),
                                                      rownames(TAD_based_territory_assign_list[[i]])),]
  colnames(aatst) <- lncRNA_chosen_gt1k_uniqTiles$geneID[match(colnames(aatst),
                                                               lncRNA_chosen_gt1k_uniqTiles$gene_name)]
  aatst2 <- aatst[, colnames(lncRNA_interaction_label_matrix_cis)]
  
  aa_tad_based_eval[[i]] <- IntMatEval(real_intMat = lncRNA_interaction_label_matrix_cis,
                       pred_intMat = aatst2,
                       DisMat = lncRNA_interaction_distance_matrix_0Removed_cis)
  IntMatEval_plot(eval_out = aa_tad_based_eval[[i]], my_filename = paste0("plots/TAD_based_assignment_eval/eval_", 
                                                          names(TAD_based_territory_assign_list)[i],
                                                          ".png"))
  
}
names(aa_tad_based_eval) <- names(TAD_based_territory_assign_list)

aa_tad_based_eval$MM9_1kb_tile_TAD_HiC_ES_DI_10kb.bed.closest$General$FN

aafn <- lapply(aa_tad_based_eval, "[[", 1)

aafn2 <- unlist(lapply(aafn, "[[", 4))
aanam <- names(aafn2)
aanam2 <- unlist(mapply(paste, lapply(strsplit(unlist(lapply(strsplit(aanam, "\\."), "[[", 1)), "_"),
                               "[", c(5,6,7,8)),
                           collapse = "_"))
names(aafn2) <- aanam2

aatp2 <- unlist(lapply(aafn, "[[", 1))
names(aatp2) <- aanam2

par(mar = c(12,3,4,4))
barplot(c(aa_nearest_1_eval$General$FN, aafn2), las =2 , main = "# FN per method")

par(mar = c(12,6,4,4))
barplot(c(aa_nearest_1_eval$General$FN, aafn2), las =2 , main = "# FN per method")

aa_FN_general <- c(aa_nearest_1_eval$General$FN, aafn2)
names(aa_FN_general)[1] <- "nearest_1"
aa_TP_general <- c(aa_nearest_1_eval$General$TP, aatp2)
names(aa_TP_general)[1] <- "nearest_1"
barplot(aa_TP_general, las =2 , main = "# TP per method")
abline(h  = max(aa_TP_general), col = 2, lty = 2)

aa_tad_fn <- do.call(rbind, lapply(lapply(aa_tad_based_eval, "[[", 2), "[[", 4))
aa_tad_fn <- rbind(aa_nearest_1_eval$PerlncRNA$FN, aa_tad_fn)
rownames(aa_tad_fn) <- c("nearest_1", aanam2)

par(mar = c(10,3,4,4))
barplot(log10(aa_tad_fn + 1), beside = T, las =2 )

png(filename = "plots/TAD_based_assignment_eval/FN_eval_perlncRNA.png",    # create PNG for the heat map        
    width = 12*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10)        # smaller font size


par(mfrow = c(4,7),mar = c(9,3,6,2))
for(i in 1:ncol(aa_tad_fn)){
  barplot(aa_tad_fn[, i], las = 2, main = colnames(aa_tad_fn)[i])
}
dev.off()

aa_tad_fn_filt <- aa_tad_fn[, colSums(aa_tad_fn) > 0]

barplot(aa_tad_fn_filt, beside = T)

which.min(aa_FN_general)
aa_tad_Chosen <- TAD_based_distance_list$MM9_1kb_tile_TAD_HiChIP_ES_GMAP_50kb.bed.closest


aa_tad_Chosen2 <- aa_tad_Chosen
colnames(aa_tad_Chosen2) <- lncRNA_chosen_gt1k_uniqTiles$geneID[match(colnames(aa_tad_Chosen2),
                                                                      lncRNA_chosen_gt1k_uniqTiles$gene_name )]


aa_tad_Chosen2 <- aa_tad_Chosen2[, colnames(lncRNA_interaction_distance_matrix_0Removed_cis)]
aa_tad_Chosen2 <- aa_tad_Chosen2[match(rownames(lncRNA_interaction_distance_matrix_0Removed_cis), 
                                       rownames(aa_tad_Chosen2)), ]
png(filename = "plots/TAD_HiChIP_ES_GMAP_50kb_distance_boxplot.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 12*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10)        # smaller font size

par(mfrow = c(7,4), mar = c(4,4,4,4))
for(i in 1:ncol(lncRNA_interaction_distance_matrix_0Removed_cis)){
  aa_cname <- lncRNA_chosen_gt1k_uniqTiles$gene_name[lncRNA_chosen_gt1k_uniqTiles$geneID == colnames(aa_tad_Chosen2)[i]]
  aa_one_nu <- sum(!is.na(aa_tad_Chosen2[, i]) & lncRNA_interaction_label_matrix_cis[,i] == 1)
  aa_zer_nu <- sum(!is.na(aa_tad_Chosen2[, i]) & lncRNA_interaction_label_matrix_cis[,i] == 0)
  boxplot(aa_tad_Chosen2[,i]~lncRNA_interaction_label_matrix_cis[,i],
          main = paste0(aa_cname, "\n#1: ",aa_one_nu,  " #0: ", aa_zer_nu))
}
dev.off()

aa_near1_Chosen <- aa_nearest_1_Pred


# aa_tad_Chosen2 <- aa_near1_Chosen
# colnames(aa_tad_Chosen2) <- lncRNA_chosen_gt1k_uniqTiles$geneID[match(colnames(aa_tad_Chosen2),
#                                                                       lncRNA_chosen_gt1k_uniqTiles$gene_name )]
# 
# 
# aa_tad_Chosen2 <- aa_tad_Chosen2[, colnames(lncRNA_interaction_distance_matrix_0Removed_cis)]
# aa_tad_Chosen2 <- aa_tad_Chosen2[match(rownames(lncRNA_interaction_distance_matrix_0Removed_cis), 
#                                        rownames(aa_tad_Chosen2)), ]
png(filename = "plots/Nearest_1_assigned_boxplot.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 12*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10)        # smaller font size

par(mfrow = c(7,4), mar = c(4,4,4,4))
for(i in 1:ncol(lncRNA_interaction_distance_matrix_0Removed_cis)){
  aa_cname <- lncRNA_chosen_gt1k_uniqTiles$gene_name[lncRNA_chosen_gt1k_uniqTiles$geneID == colnames(aa_near1_Chosen)[i]]
  aa_one_nu <- sum(!is.na(aa_near1_Chosen[, i]) & lncRNA_interaction_label_matrix_cis[,i] == 1)
  aa_zer_nu <- sum(!is.na(aa_near1_Chosen[, i]) & lncRNA_interaction_label_matrix_cis[,i] == 0)
  boxplot(aa_near1_Chosen[,i]~lncRNA_interaction_label_matrix_cis[,i],
          main = paste0(aa_cname, "\n#1: ",aa_one_nu,  " #0: ", aa_zer_nu))
}
dev.off()

####################################################################################################
# repeat this with correcting for the ties
save(list = c("DisChaLDis"), file = "DisChaLDis_tie_fixed.RData")
TAD_based_territory_assign_list_tie_cor <- list()
# for(i in 1:length(TAD_based_distance_list)){
#   print(i)
#   TAD_based_territory_assign_list_tie_cor[[i]]  <- DisChaLDis(DisMat = TAD_based_distance_list[[i]],
#                                                       k_closest = 1, 
#                                                       tie_break = T)
# }
# names(TAD_based_territory_assign_list_tie_cor) <- names(TAD_based_distance_list)
# save(list = c("TAD_based_territory_assign_list_tie_cor"), file = "Territory_assign_output_tie_cor.RData")
#load("/Users/Shayan/Territory_assign_output_tie_cor.RData")
aa_tad_based_eval_tie <- list()
for(i in 1:length(TAD_based_territory_assign_list_tie_cor)){
  print(i)
  aatst <- TAD_based_territory_assign_list_tie_cor[[i]][match(rownames(lncRNA_interaction_distance_matrix_0Removed_cis),
                                                      rownames(TAD_based_territory_assign_list_tie_cor[[i]])),]
  colnames(aatst) <- lncRNA_chosen_gt1k_uniqTiles$geneID[match(colnames(aatst),
                                                               lncRNA_chosen_gt1k_uniqTiles$gene_name)]
  aatst2 <- aatst[, colnames(lncRNA_interaction_label_matrix_cis)]
  
  aa_tad_based_eval_tie[[i]] <- IntMatEval(real_intMat = lncRNA_interaction_label_matrix_cis,
                                       pred_intMat = aatst2,
                                       DisMat = lncRNA_interaction_distance_matrix_0Removed_cis)
  IntMatEval_plot(eval_out = aa_tad_based_eval_tie[[i]], my_filename = paste0("plots/TAD_based_assignment_eval/eval_tieCor_", 
                                                          names(TAD_based_territory_assign_list_tie_cor)[i],
                                                          ".png"))
  
}
names(aa_tad_based_eval_tie) <- names(TAD_based_territory_assign_list)

aa_tad_based_eval_tie$MM9_1kb_tile_TAD_HiC_ES_DI_10kb.bed.closest$General$FN

aafn <- lapply(aa_tad_based_eval_tie, "[[", 1)

aafn2 <- unlist(lapply(aafn, "[[", 4))
aanam <- names(aafn2)
aanam2 <- unlist(mapply(paste, lapply(strsplit(unlist(lapply(strsplit(aanam, "\\."), "[[", 1)), "_"),
                                      "[", c(5,6,7,8)),
                        collapse = "_"))
names(aafn2) <- aanam2

aatp2 <- unlist(lapply(aafn, "[[", 1))
names(aatp2) <- aanam2

par(mar = c(12,3,4,4))
barplot(c(aa_nearest_1_eval$General$FN, aafn2), las =2 , main = "# FN per method")

par(mar = c(12,6,4,4), xpd = F)

aa_FN_general <- c(aa_nearest_1_eval$General$FN, aafn2)
names(aa_FN_general)[1] <- "nearest_1"
aa_TP_general <- c(aa_nearest_1_eval$General$TP, aatp2)
names(aa_TP_general)[1] <- "nearest_1"
barplot(aa_FN_general, las =2 , main = "# FN per method")
abline(h  = min(aa_FN_general), col = 2, lty = 2)

barplot(aa_TP_general, las =2 , main = "# TP per method")
abline(h  = max(aa_TP_general), col = 2, lty = 2)

aa_tad_fn <- do.call(rbind, lapply(lapply(aa_tad_based_eval_tie, "[[", 2), "[[", 4))
aa_tad_fn <- rbind(aa_nearest_1_eval$PerlncRNA$FN, aa_tad_fn)
rownames(aa_tad_fn) <- c("nearest_1", aanam2)

par(mar = c(10,3,4,4))
barplot(log10(aa_tad_fn + 1), beside = T, las =2 )

png(filename = "plots/TAD_based_assignment_eval/FN_eval_perlncRNA_tieFixed.png",    # create PNG for the heat map        
    width = 12*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10)        # smaller font size


par(mfrow = c(4,7),mar = c(9,3,6,2))
for(i in 1:ncol(aa_tad_fn)){
  barplot(aa_tad_fn[, i], las = 2, main = colnames(aa_tad_fn)[i])
}
dev.off()

aa_tad_fn_filt <- aa_tad_fn[, colSums(aa_tad_fn) > 0]

barplot(aa_tad_fn_filt, beside = T, las = 2)

which.min(aa_FN_general)
aa_tad_Chosen <- TAD_based_distance_list$MM9_1kb_tile_TAD_HiChIP_ES_GMAP_50kb.bed.closest


aa_tad_Chosen2 <- aa_tad_Chosen
colnames(aa_tad_Chosen2) <- lncRNA_chosen_gt1k_uniqTiles$geneID[match(colnames(aa_tad_Chosen2),
                                                                      lncRNA_chosen_gt1k_uniqTiles$gene_name )]


aa_tad_Chosen2 <- aa_tad_Chosen2[, colnames(lncRNA_interaction_distance_matrix_0Removed_cis)]
aa_tad_Chosen2 <- aa_tad_Chosen2[match(rownames(lncRNA_interaction_distance_matrix_0Removed_cis), 
                                       rownames(aa_tad_Chosen2)), ]
png(filename = "plots/TAD_HiChIP_ES_GMAP_50kb_distance_boxplot_tie_corrected.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 12*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10)        # smaller font size

par(mfrow = c(7,4), mar = c(4,4,4,4))
for(i in 1:ncol(lncRNA_interaction_distance_matrix_0Removed_cis)){
  aa_cname <- lncRNA_chosen_gt1k_uniqTiles$gene_name[lncRNA_chosen_gt1k_uniqTiles$geneID == colnames(aa_tad_Chosen2)[i]]
  aa_one_nu <- sum(!is.na(aa_tad_Chosen2[, i]) & lncRNA_interaction_label_matrix_cis[,i] == 1)
  aa_zer_nu <- sum(!is.na(aa_tad_Chosen2[, i]) & lncRNA_interaction_label_matrix_cis[,i] == 0)
  boxplot(aa_tad_Chosen2[,i]~lncRNA_interaction_label_matrix_cis[,i],
          main = paste0(aa_cname, "\n#1: ",aa_one_nu,  " #0: ", aa_zer_nu))
}
dev.off()

aa_tad_Chosen <- TAD_based_territory_assign_list_tie_cor$MM9_1kb_tile_TAD_HiChIP_ES_GMAP_50kb.bed.closest


aa_tad_Chosen2 <- aa_tad_Chosen
colnames(aa_tad_Chosen2) <- lncRNA_chosen_gt1k_uniqTiles$geneID[match(colnames(aa_tad_Chosen2),
                                                                      lncRNA_chosen_gt1k_uniqTiles$gene_name )]


aa_tad_Chosen2 <- aa_tad_Chosen2[, colnames(lncRNA_interaction_distance_matrix_0Removed_cis)]
aa_tad_Chosen2 <- aa_tad_Chosen2[match(rownames(lncRNA_interaction_distance_matrix_0Removed_cis), 
                                       rownames(aa_tad_Chosen2)), ]
png(filename = "plots/TAD_HiChIP_ES_GMAP_50kb_assigned_boxplot_tie_corrected.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 12*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10)        # smaller font size

par(mfrow = c(7,4), mar = c(4,4,4,4))
for(i in 1:ncol(lncRNA_interaction_distance_matrix_0Removed_cis)){
  aa_cname <- lncRNA_chosen_gt1k_uniqTiles$gene_name[lncRNA_chosen_gt1k_uniqTiles$geneID == colnames(aa_tad_Chosen2)[i]]
  aa_one_nu <- sum(!is.na(aa_tad_Chosen2[, i]) & lncRNA_interaction_label_matrix_cis[,i] == 1)
  aa_zer_nu <- sum(!is.na(aa_tad_Chosen2[, i]) & lncRNA_interaction_label_matrix_cis[,i] == 0)
  boxplot(aa_tad_Chosen2[,i]~lncRNA_interaction_label_matrix_cis[,i],
          main = paste0(aa_cname, "\n#1: ",aa_one_nu,  " #0: ", aa_zer_nu))
}
dev.off()




# how many tiles have equal distance to MALAT1 and NEAT1
nrow(TAD_based_distance_list$MM9_1kb_tile_TAD_HiChIP_ES_GMAP_50kb.bed.closest)
colnames(TAD_based_distance_list$MM9_1kb_tile_TAD_HiChIP_ES_GMAP_50kb.bed.closest)

sum(TAD_based_distance_list$MM9_1kb_tile_TAD_HiChIP_ES_GMAP_50kb.bed.closest[, "Platr28"] == 
      TAD_based_distance_list$MM9_1kb_tile_TAD_HiChIP_ES_GMAP_50kb.bed.closest[, "Kcnq1ot1"], 
    na.rm = T )
sum(lncRNA_chosen_gt1k_TADnu["Firre",] == 
      lncRNA_chosen_gt1k_TADnu["Gm14820",], 
    na.rm = F )
View(lncRNA_chosen_gt1k_TADnu)
####################################################################################################
# evaluate a new assignment method

aa_dist_new <- TAD_based_distance_list$MM9_1kb_tile_TAD_HiChIP_ES_GMAP_50kb.bed.closest
#aa_dist_new[, c("Malat1", "Neat1", "2410003L11Rik", "2610035D17Rik", "Gm53", "Gm11613", "Meg3", "Rian")] <- 
aa_dist_new[is.na(lncRNA_interaction_distance_matrix_full_0Removed)] <- NA

aa_old_dist <- lncRNA_interaction_distance_matrix_full_0Removed
all(is.na(aa_dist_new[,1]) == is.na(aa_old_dist[, 1]))


# colnames(aa_old_dist) <- lncRNA_chosen_gt1k_uniqTiles$gene_name[match(colnames(aa_old_dist),
#                                                                       lncRNA_chosen_gt1k_uniqTiles$geneID)]
# aa_old_dist2 <- TAD_based_distance_list$MM9_1kb_tile_TAD_HiChIP_ES_GMAP_50kb.bed.closest

aa_dist_new[,
             c("Malat1", "Neat1", 
               "2410003L11Rik", "2610035D17Rik",
               "Gm53", "Gm11613",
               "Meg3", "Rian")] <- aa_old_dist[, c("Malat1", "Neat1", 
                                                   "2410003L11Rik", "2610035D17Rik",
                                                   "Gm53", "Gm11613", "Meg3", "Rian")]
sum((aa_dist_new[,1] == aa_old_dist[, 1]), na.rm = T) == sum(!is.na(aa_old_dist[, 1]))

Hybrid_distance_matrix <-  aa_dist_new

# assign territories using the new 

Hybrid_based_territory_assignment  <- DisChaLDis(DisMat = Hybrid_distance_matrix,
                                                            k_closest = 1,
                                                            tie_break = T)
aatst <- Hybrid_based_territory_assignment[match(rownames(lncRNA_interaction_distance_matrix_0Removed_cis),
                                                            rownames(Hybrid_based_territory_assignment)),]
colnames(aatst) <- lncRNA_chosen_gt1k_uniqTiles$geneID[match(colnames(aatst),
                                                             lncRNA_chosen_gt1k_uniqTiles$gene_name)]
aatst2 <- aatst[, colnames(lncRNA_interaction_label_matrix_cis)]

aa_hybrid_eval <- IntMatEval(real_intMat = lncRNA_interaction_label_matrix_cis,
                                         pred_intMat = aatst2,
                                         DisMat = lncRNA_interaction_distance_matrix_0Removed_cis)
IntMatEval_plot(eval_out = aa_hybrid_eval,
                my_filename = "~/Documents/Shayan/BioInf/lncRNA/plots/Hybrid_based_eval.png")
barplot(rbind(aa_hybrid_eval$PerlncRNA$FN, 
              aa_tad_based_eval_tie$MM9_1kb_tile_TAD_HiChIP_ES_GMAP_50kb.bed.closest$PerlncRNA$FN),
        las = 2, beside = T)
##### it didn't work, first compute assignments based on their own matrices and then do the change at the assignment level
aa_dist_new <- TAD_based_distance_list$MM9_1kb_tile_TAD_HiChIP_ES_GMAP_50kb.bed.closest
aa_dist_new <- aa_dist_new[, match(colnames(lncRNA_interaction_distance_matrix_full), colnames(aa_dist_new))]


#aa_dist_new[, c("Malat1", "Neat1", "2410003L11Rik", "2610035D17Rik", "Gm53", "Gm11613", "Meg3", "Rian")] <- 
for(i in 1:ncol(aa_dist_new)){
  aa_dist_new[lncRNA_interaction_distance_matrix_full[,i]==0, i] <- NA
  
}

aa_old_dist <- lncRNA_interaction_distance_matrix_full_0Removed

aaHybrid_based_territory_assignment1  <- DisChaLDis(DisMat = aa_old_dist,
                                                 k_closest = 1,
                                                 tie_break = T)
aaHybrid_based_territory_assignment2  <- DisChaLDis(DisMat = aa_dist_new,
                                                 k_closest = 1,
                                                 tie_break = T)

Hybrid_based_territory_assignment <- aaHybrid_based_territory_assignment2
Hybrid_based_territory_assignment[, c("Malat1", "Neat1", "2410003L11Rik", "2610035D17Rik", "Gm53", "Gm11613", "Meg3", "Rian", "Gm28373")] <- aaHybrid_based_territory_assignment1[, c("Malat1", "Neat1", "2410003L11Rik", "2610035D17Rik", "Gm53", "Gm11613", "Meg3", "Rian", "Gm28373")]

Distance_based_territory_assignment <- aaHybrid_based_territory_assignment1
aatst <- Hybrid_based_territory_assignment[match(rownames(lncRNA_interaction_distance_matrix_0Removed_cis),
                                                 rownames(Hybrid_based_territory_assignment)),]
colnames(aatst) <- lncRNA_chosen_gt1k_uniqTiles$geneID[match(colnames(aatst),
                                                             lncRNA_chosen_gt1k_uniqTiles$gene_name)]
aatst2 <- aatst[, colnames(lncRNA_interaction_label_matrix_cis)]

aa_hybrid_eval2 <- IntMatEval(real_intMat = lncRNA_interaction_label_matrix_cis,
                             pred_intMat = aatst2,
                             DisMat = lncRNA_interaction_distance_matrix_0Removed_cis)
aatst <- aaHybrid_based_territory_assignment1[match(rownames(lncRNA_interaction_distance_matrix_0Removed_cis),
                                                 rownames(aaHybrid_based_territory_assignment1)),]
colnames(aatst) <- lncRNA_chosen_gt1k_uniqTiles$geneID[match(colnames(aatst),
                                                             lncRNA_chosen_gt1k_uniqTiles$gene_name)]
aatst2 <- aatst[, colnames(lncRNA_interaction_label_matrix_cis)]


aa_dist_eval2 <- IntMatEval(real_intMat = lncRNA_interaction_label_matrix_cis,
                              pred_intMat = aatst2,
                              DisMat = lncRNA_interaction_distance_matrix_0Removed_cis)

aatst <- aaHybrid_based_territory_assignment2[match(rownames(lncRNA_interaction_distance_matrix_0Removed_cis),
                                                    rownames(aaHybrid_based_territory_assignment2)),]
colnames(aatst) <- lncRNA_chosen_gt1k_uniqTiles$geneID[match(colnames(aatst),
                                                             lncRNA_chosen_gt1k_uniqTiles$gene_name)]
aatst2 <- aatst[, colnames(lncRNA_interaction_label_matrix_cis)]


aa_tad_eval2 <- IntMatEval(real_intMat = lncRNA_interaction_label_matrix_cis,
                            pred_intMat = aatst2,
                            DisMat = lncRNA_interaction_distance_matrix_0Removed_cis)


IntMatEval_plot(eval_out = aa_hybrid_eval2,
                my_filename = "~/Documents/Shayan/BioInf/lncRNA/plots/Hybrid_based_eval_2.png")

aabpl <- rbind(aa_hybrid_eval2$PerlncRNA$FN, 
               aa_dist_eval2$PerlncRNA$FN, 
               aa_tad_eval2$PerlncRNA$FN)
rownames(aabpl) <- c("Hybrid", "Distance-based", "TAD-based")
colnames(aabpl) <- names(aa_hybrid_eval2$PerlncRNA$FN)

aabpl <- aabpl[, sort(colSums(aabpl), decreasing = T, index.return = T)$ix]
barplot(aabpl,
        las = 2, beside = T, ylab = "# FN"
        , legend.text = rownames(aabpl),    
        args.legend=list(x = "topright", bty = "n", 
                         inset=c(0.2, 0),
                         cex = 0.7, 
                         y.intersp = 0.6,
                         x.intersp = 0.5,
                         text.width = 1)
)
# based on this plot I choose to work with the bp distance-based method, while using same tad membership as a pair feature


####################################################################################################
# choose the territory assignment method
# based on the number of FN: I choose HiChIP_ES_GMAP_50kb # NO THAT WAS A BUG (tie breaking)
# here are number of False negatives per method 
# barplot(c(aa_nearest_1_eval$General$FN, aafn2), las =2 , main = "# FN per method")

# TAD_based_territory_assign_list$MM9_1kb_tile_TAD_HiChIP_ES_GMAP_50kb.bed.closest


# this is the chosen territory assignment:
Distance_based_territory_assignment
# add territorial owner to MM9_1kb_tiled
MM9_1kb_tiled


aatst <- apply(Distance_based_territory_assignment, 
               MARGIN = 1, 
               function(x) colnames(Distance_based_territory_assignment)[x == 1])
aaln <- unlist(lapply(aatst, length))
aatst2 <- aatst
aatst2[aaln == 0] <- NA
MM9_1kb_tiled_owner <- unlist(aatst2)
#MM9_1kb_tiled_owner[is.na(MM9_1kb_tiled_owner)] <- "unassigned"
par(mfrow = c(1,1), mar = c(10,5,4,4))
barplot( sort(table(MM9_1kb_tiled_owner), decreasing = T),
         las = 2,
         main = "Territory members", 
         ylab = "")
names(MM9_1kb_tiled_owner) <- rownames(Distance_based_territory_assignment)

MM9_1kb_tiled_owner_labels <- numeric(length = length(MM9_1kb_tiled_owner))
names(MM9_1kb_tiled_owner_labels) <- names(MM9_1kb_tiled_owner)

# which ones are the real interactions among the defined ones
territory_assigned_interaction <- Distance_based_territory_assignment
territory_assigned_interaction_label <- territory_assigned_interaction
territory_assigned_interaction_label[] <- 0
aalncRNA_interaction_label_matrix_cis <- lncRNA_interaction_label_matrix_cis
colnames(aalncRNA_interaction_label_matrix_cis) <- lncRNA_chosen_gt1k_uniqTiles$gene_name[match(colnames(lncRNA_interaction_label_matrix_cis), lncRNA_chosen_gt1k_uniqTiles$geneID)]
aalncRNA_interaction_label_matrix_cis <- aalncRNA_interaction_label_matrix_cis[, colnames(territory_assigned_interaction_label)]



territory_assigned_interaction_label[match(rownames(lncRNA_interaction_label_matrix_cis), 
                                           rownames(territory_assigned_interaction_label)), ] <- aalncRNA_interaction_label_matrix_cis
territory_assigned_interaction_label_filtered <- territory_assigned_interaction_label
territory_assigned_interaction_label_filtered[territory_assigned_interaction == 0] <- 0


sum(territory_assigned_interaction)
sum(territory_assigned_interaction_label)
sum(territory_assigned_interaction_label_filtered)


colSums(territory_assigned_interaction)
colSums(territory_assigned_interaction_label)
colSums(territory_assigned_interaction_label_filtered)

aa_colsums <- rbind(rbind(colSums(territory_assigned_interaction),
            colSums(territory_assigned_interaction_label)), 
      colSums(territory_assigned_interaction_label_filtered))
aa_colsums <- aa_colsums[,sort(aa_colsums[3,], decreasing = T, index.return = T)$ix]

par(mfrow = c(1,1), mar = c(8,5,6,6), xpd = T)
barplot(log10(aa_colsums + 1),
        beside = T, las = 2, 
        legend.text = c("# tile_owner_pairs", "TP + FN", "TP"),
        main = "",
        ylab = "log10(# of interactions)",    
       #legend.text=TRUE,
       args.legend=list(x = "topright", bty = "n", 
                        inset=c(-0.05, 0),
                        cex = 0.7, 
                        y.intersp = 0.6,
                        x.intersp = 0.5,
                        text.width = 1)
       )
#legend("top right", legend = c("# tile_owner_pairs", "TP + FN", "TP"), fill = )

aa_colsums <- aa_colsums[,sort(aa_colsums[3,], decreasing = T, index.return = T)$ix]
par(xpd = F)
barplot(aa_colsums,
        beside = T, las = 2, 
        legend.text = c("# tile_owner_pairs", "TP + FN", "TP"),
        main = "",
        ylab = "# of interactions",    
        #legend.text=TRUE,
        args.legend=list(x = "topright", bty = "n", 
                         inset=c(-0.05, 0),
                         cex = 0.7, 
                         y.intersp = 0.6,
                         x.intersp = 0.5,
                         text.width = 1)
)
abline(h= seq(0,180000, 10000), col = "red", lty = 2, lwd = 0.6)

barplot(aa_colsums[3,], las = 2, main = "# TP per lncRNA")
abline(h= seq(0,15000, 1000), col = "red", lty = 2, lwd = 0.6)

barplot(aa_colsums[2,] - aa_colsums[3,], las = 2, main = "# FN per lncRNA")
abline(h= seq(0,5000, 500), col = "red", lty = 2, lwd = 0.6)

aatst <- apply(territory_assigned_interaction_label_filtered, 
               MARGIN = 1, 
               function(x) ifelse(test = (sum(x) == 0), yes = NA, no = colnames(territory_assigned_interaction_label_filtered)[x == 1]) )
#aaln <- unlist(lapply(aatst, length))
aatst2 <- aatst
#aatst2[aaln == 0] <- NA
MM9_1kb_tiled_owner_labels <- unlist(aatst2)
names(MM9_1kb_tiled_owner_labels) <- names(MM9_1kb_tiled_owner) 
MM9_1kb_tiled_owner_labels_binary <- !is.na(MM9_1kb_tiled_owner_labels)
names(MM9_1kb_tiled_owner_labels_binary) <- names(MM9_1kb_tiled_owner) 

MM9_1kb_tiled_owner_labels_binary_3Mfilter_noNA

aatst2 <- apply(territory_assigned_interaction_3Mfilered, 
               MARGIN = 1, 
               function(x) ifelse(test = (sum(x) == 0), yes = NA, no = colnames(territory_assigned_interaction_3Mfilered)[x == 1]) )
sum(!is.na(MM9_1kb_tiled_owner_3Mfilter_noNA))


######################################################################################################
######################################################################################################
# make tiles or hg38
system("bedtools makewindows -g ~/Documents/Shayan/BioInf/lncRNA/hg38.genome -w 1000 > ~/Documents/Shayan/BioInf/lncRNA/hg38_1kb_tiles.bed")
HG38_1kb_tiled <- read.table(file = "~/Documents/Shayan/BioInf/lncRNA/hg38_1kb_tiles.bed", stringsAsFactors = F)
aa_unq <- unique(HG38_1kb_tiled$V1)
aa_st_ps <- numeric(length(aa_unq))
aa_st_nu <- numeric(length(aa_unq))

for(i in 1:length(aa_st_ps)){
  aax <- which(HG38_1kb_tiled$V1 %in% aa_unq[i])
  aa_st_ps[i] <- aax[1]
  aa_st_nu[i] <- length(aax)
}
aa_nms <- character(0)
for(i in 1:length(aa_st_ps)){
  aa_nms <- c(aa_nms, paste0(aa_unq[i], "_", c(1:aa_st_nu[i])))
}
HG38_1kb_tiled$tile <- aa_nms
names(HG38_1kb_tiled) <- c("chr", "start", "end", "tile")
HG38_1kb_tiled$start <- HG38_1kb_tiled$start + 1
HG38_1kb_tiled_GR <- makeGRangesFromDataFrame(df = HG38_1kb_tiled, keep.extra.columns = T)

#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)

which(HG38_1kb_tiled_GR$tile == "chr2_11")
aa_tile_Seq <- getSeq(x = BSgenome.Hsapiens.UCSC.hg38, names = HG38_1kb_tiled_GR[248968])
aa_tile_Seq

# IT SEEMS LIKE THE first 10 tiles if each chromosome are NNNs
aaHG38_1kb_tiled_GR_filtered <- HG38_1kb_tiled_GR
aarem <- numeric(0)
aa_unq_chr <- unique(HG38_1kb_tiled$chr)
for(i in 1:length(aa_unq_chr)){
  aarem <-  c(aarem, which(HG38_1kb_tiled$chr == aa_unq_chr[i])[1:min(10, sum(HG38_1kb_tiled$chr == aa_unq_chr[i]))])
}
aaHG38_1kb_tiled_GR_filtered <- aaHG38_1kb_tiled_GR_filtered[-c(aarem)]
index_first_10_perchr <- aarem
HG38_1kb_tiled_GR_filtered <- aaHG38_1kb_tiled_GR_filtered
save(list = c("HG38_1kb_tiled_GR_filtered"), file = "~/Documents/Shayan/BioInf/lncRNA/HG38_1kb_tiled_GR_filtered.RData")

aa_tile_Seq <- getSeq(x = BSgenome.Hsapiens.UCSC.hg38, names = HG38_1kb_tiled_GR_filtered[1:1000])
names(aa_tile_Seq) <- HG38_1kb_tiled_GR_filtered$tile[1:1000]

# getting kmers for hg38

# write all the tiles in one chromsome to a file

aa_tile_Seq <- getSeq(x = BSgenome.Hsapiens.UCSC.hg38, names = HG38_1kb_tiled_GR_filtered[3001])
aa_tile_Seq_kmer <- get.kmers(as.character(aa_tile_Seq), .k = 5)
aa_tile_Seq_kmer

aa_myseq <- getSeq(x = BSgenome.Hsapiens.UCSC.hg38, names = HG38_1kb_tiled_GR_filtered )
names(aa_myseq) <- HG38_1kb_tiled_GR_filtered$tile
aasppp <- unlist(lapply(strsplit(names(aa_myseq), "_"), "[[", 1))
aaunq1 <- unique(aasppp)
#aa_unq_chr <-levels(seqnames(HG38_1kb_tiled_GR_filtered))
for(i in 1:length(aaunq1)){
  print(aaunq1[i])
  aaw <- which(aasppp %in% aaunq1[i])
  ShortRead::writeFasta(object = aa_myseq[aaw], 
                        file = paste0("~/Documents/Shayan/BioInf/lncRNA/HG38_tiles_perChr/hg38_", aaunq1[i], ".fasta"))
  
}

# write bash to move these to hal
cat("#!/bin/bash\n", file = "~/Documents/Shayan/BioInf/lncRNA/move_chr_hg38.sh", append = F)
for(i in 1:length(aaunq1)){
  cat(c("rsync",
        paste0("~/Documents/Shayan/BioInf/lncRNA/HG38_tiles_perChr/hg38_", aaunq1[i], ".fasta"),
        paste0("/Users/Shayan/remote/hg38_tiles/", aaunq1[i], "\n")),
      file = "~/Documents/Shayan/BioInf/lncRNA/move_chr_hg38.sh", sep = " ", append = T)
}



# using Jellyfish

cat("#!/bin/bash\n", file = "~/Documents/Shayan/BioInf/lncRNA/break_to_tiles_hg38.sh", append = F)
for(i in 1:length(aaunq1)){
  cat(c("cd",paste0("/shared-mounts/sinhas/tabebor2/hg38_tiles/",aaunq1[i]), "\n"),
      file = "~/Documents/Shayan/BioInf/lncRNA/break_to_tiles_hg38.sh", append = T, sep = " ")
  cat(c("/shared-mounts/sinhas/tabebor2/lncRNA/break_fasta.sh <", paste0("hg38_", aaunq1[i], ".fasta"), "\n"),
      file = "~/Documents/Shayan/BioInf/lncRNA/break_to_tiles_hg38.sh", append = T, sep = " " )
 # cat(c("rm", paste0("hg38_", aaunq1[i], ".fasta"), "\n"),
  #    file = "~/Documents/Shayan/BioInf/lncRNA/break_to_tiles_hg38.sh", append = T, sep = " " )
}

# write jellyfish jobs
aasppp <- unlist(lapply(strsplit(HG38_1kb_tiled_GR_filtered$tile, "_"), "[[", 1))

for(i in 1:length(HG38_1kb_tiled_GR_filtered)){
  if((i %% 1000) == 0){
    print(i)
  }
  cat(c("jellyfish count -m 5 -t 1 --text  -s 1024 -o", 
        paste0("/shared-mounts/sinhas/tabebor2/hg38_tiles/hg_38_kmers/",aasppp[i], "/", HG38_1kb_tiled_GR_filtered$tile[i], ".5mer"), 
        paste0("/shared-mounts/sinhas/tabebor2/hg38_tiles/",aasppp[i], "/", HG38_1kb_tiled_GR_filtered$tile[i], ".fa\n")),
      file = "~/Documents/Shayan/BioInf/lncRNA/Jellyfish_hg39.job", append = !(i == 1), sep = " " )
}

#################################################################################################################################################################
# read RADICL-seq data on mOPC --> to see variation in binding for different lncRNAs



aafiles_full <- list.files("~/Documents/Shayan/BioInf/lncRNA/RADICL-seq/mOPC", pattern = "*.gz", full.names = T)
aafiles <- list.files("~/Documents/Shayan/BioInf/lncRNA/RADICL-seq/mOPC", pattern = "*.gz", full.names = F)
aafiles_2 <- unlist(strsplit(aafiles, split = ".gz"))

# write bash file to get only lncRNA related interactions in all of them
setwd("~/Documents/Shayan/BioInf/lncRNA/RADICL-seq/mOPC")
cat(c("#!/bin/bash\n"), file = "getlncRNA_ints_mOPC.sh",append = F)
for(i in 1:length(aafiles_full)){
  cat(c("gunzip ", aafiles_full[i], "\n"), file = "getlncRNA_ints_mOPC.sh",sep = " ", append = T)
  cat(c("awk '$8 == %long_ncRNA% {print}' ", aafiles_2[i], " > ", paste0(aafiles_2[i], "_lnc"), "\n"), file = "getlncRNA_ints_mOPC.sh",sep = " ", append = T)
  cat(c("gzip ", aafiles_2[i], "\n"), file = "getlncRNA_ints_mOPC.sh",sep = " ", append = T)
}

setwd("~/Documents/Shayan/BioInf/lncRNA")

RADICLseq_mm10_mOPC <-  read.table("~/Documents/Shayan/BioInf/lncRNA/RADICL-seq/mOPC/GSE132190_all_mOPC_lncRNA_interactions.txt", header = F, stringsAsFactors = F)
#colnames(RADICLseq_mm10) <- ""
table(RADICLseq_mm10_mOPC$V9)

head(RADICLseq_mm10_mOPC)



# convert the DNA part of RADICL from mm10 to mm9 

RADICLseq_mm10_DNA_GR_mOPC <- makeGRangesFromDataFrame(RADICLseq_mm10_mOPC[, c("V11", "V12", "V13")],
                                                  start.field = "V12", 
                                                  seqnames.field = "V11", 
                                                  end.field = "V13")
names(RADICLseq_mm10_DNA_GR_mOPC)<- c(1:length(RADICLseq_mm10_DNA_GR_mOPC))
aa_chain <- import.chain("~/Documents/Shayan/BioInf/Liftover_chain_files/mm10ToMm9.over.chain")
aaRADICLseq_mm9_DNA_GR_liftover_mOPC  <- liftOver(x = RADICLseq_mm10_DNA_GR_mOPC, chain = aa_chain)
#aal <- unlist(lapply(RADICLseq_mm9_DNA_GR, length))
RADICLseq_mm9_DNA_GR_mOPC <- RADICLseq_mm10_DNA_GR_mOPC
RADICLseq_mm10_to_mm9_stat_mOPC <- numeric(length(RADICLseq_mm10_DNA_GR_mOPC))
save(list = c("RADICLseq_mm9_DNA_GR_mOPC", "aaRADICLseq_mm9_DNA_GR_liftover_mOPC", "RADICLseq_mm10_to_mm9_stat_mOPC"),
     file = "~/Documents/Shayan/BioInf/lncRNA/RADICL-seq/mOPC/conversion_to_mm9_pre.RData")
# 
# for(i in 1:length(RADICLseq_mm9_DNA_GR_mOPC)){
#   print(i)
#   if(length(aaRADICLseq_mm9_DNA_GR_liftover_mOPC[[i]]) > 0){
#     RADICLseq_mm9_DNA_GR_mOPC[i] <- aaRADICLseq_mm9_DNA_GR_liftover_mOPC[[i]][1]
#   }else{
#     print("no")
#     RADICLseq_mm10_to_mm9_stat_mOPC[i] <- 1
#   }
# }



names(aaRADICLseq_mm9_DNA_GR_liftover_mOPC) <- c(1:length(aaRADICLseq_mm9_DNA_GR_liftover_mOPC))
names(RADICLseq_mm10_DNA_GR_mOPC) <- c(1:length(aaRADICLseq_mm9_DNA_GR_liftover_mOPC))

aaRADICLseq_mm9_DNA_GR_liftover_un_mOPC <- unlist(aaRADICLseq_mm9_DNA_GR_liftover_mOPC, use.names = T)


RADICLseq_mm10_DNA_df_mOPC <- as.data.frame(RADICLseq_mm10_DNA_GR_mOPC)
RADICLseq_mm10_DNA_df_mOPC$index <- rownames(RADICLseq_mm10_DNA_df_mOPC)
aaRADICLseq_mm9_DNA_GR_liftover_df_mOPC <- as.data.frame(aaRADICLseq_mm9_DNA_GR_liftover_un_mOPC)
aaRADICLseq_mm9_DNA_GR_liftover_df_mOPC$index <- names(aaRADICLseq_mm9_DNA_GR_liftover_un_mOPC)
head(aaRADICLseq_mm9_DNA_GR_liftover_df_mOPC)

RADICLseq_mm10_mm9_DNA_mOPC <- left_join( x = RADICLseq_mm10_DNA_df_mOPC,
                                     y = aaRADICLseq_mm9_DNA_GR_liftover_df_mOPC, 
                                     by = c("index"))

sum(is.na(RADICLseq_mm10_mm9_DNA_mOPC$seqnames.y))

RADICLseq_mm10_mOPC$index <- as.character(c(1:nrow(RADICLseq_mm10_mOPC)))


RADICLseq_mm10_mm9_mOPC <-  left_join( x = RADICLseq_mm10_mOPC,
                                  y = aaRADICLseq_mm9_DNA_GR_liftover_df_mOPC, 
                                  by = c("index"))
ncol(RADICLseq_mm10_mm9_mOPC)
head(RADICLseq_mm10_mm9_mOPC)


RADICLseq_mm10_mm9_mOPC <- RADICLseq_mm10_mm9_mOPC[!is.na(RADICLseq_mm10_mm9_mOPC$seqnames), ] 

aaRADICLseq_mm10_mm9_GR_mOPC <- makeGRangesFromDataFrame(df = RADICLseq_mm10_mm9_mOPC[, c("seqnames", "start", "end", "index")], 
                                                    keep.extra.columns = T)

# names(MM9_1kb_tiled) <- c("chr", "start", "end", "tile")
# MM9_1kb_tiled_GR <- makeGRangesFromDataFrame(df = MM9_1kb_tiled, keep.extra.columns = T)
# end(MM9_1kb_tiled_GR) <- end(MM9_1kb_tiled_GR) - 1
aaov4 <- findOverlaps(query = aaRADICLseq_mm10_mm9_GR_mOPC, subject = MM9_1kb_tiled_GR)

all(aaov4@from == c(1:length(aaov4@from)))

aanewtile <- MM9_1kb_tiled_GR$tile[aaov4@to]

RADICLseq_mm10_mm9_mOPC$kbTile <- aanewtile
####################################
aarad <- RADICLseq_mm10_mm9_mOPC[, c("seqnames", "start", "end", "V6", "V7", "V16", "V18", "kbTile")]
aaradic <- unlist(lapply(strsplit((RADICLseq_mm10_mm9_mOPC$V7), "\\."), "[[", 1))
aarad$V7 <- aaradic
names(aarad) <- c("chr","start", "end", "strand", "geneID", "experiment", "score", "kbTile")
RADICL_Tiled_df_mOPC <- aarad
RADICL_Tiled_df_mOPC_short <- RADICL_Tiled_df_mOPC[,c(5,6,7,8)]
#RADICL_Tiled_df_mOPC_short$kbTile <- as.character(levels(RADICL_Tiled_df_mOPC_short$kbTile)[as.numeric(RADICL_Tiled_df_mOPC_short$kbTile)])

# RADICL_GRID_Tiled_df <- rbind(aarad, aagrid)
# RADICL_GRID_Tiled_df_short <- RADICL_GRID_Tiled_df[,c(5,6,7,8)]
# RADICL_GRID_Tiled_df_short$kbTile <- as.character(levels(RADICL_GRID_Tiled_df_short$kbTile)[as.numeric(RADICL_GRID_Tiled_df_short$kbTile)])

########################################################################################################################################################################
########################################################################################################################################################################
# aggregate RADICL_GRID_Tiled_df_short based on GeneID and then sort based on the number of interactions
head(RADICL_Tiled_df_mOPC_short)

RADICL_Tiled_df_mOPC_short_byGeneID <- aggregate(x = RADICL_Tiled_df_mOPC_short,
                                                 by = RADICL_Tiled_df_mOPC_short[c("geneID")], 
                                                 FUN = c)
aa_tab <- lapply(RADICL_Tiled_df_mOPC_short_byGeneID$experiment, table)
RADICL_Tiled_df_mOPC_short_byGeneID$Exp_table <- aa_tab
aa_len <- unlist(lapply(RADICL_Tiled_df_mOPC_short_byGeneID$experiment, length))
aa_len_uniq <- unlist(lapply(lapply(RADICL_Tiled_df_mOPC_short_byGeneID$experiment, unique), length))

RADICL_Tiled_df_mOPC_short_byGeneID_modified <- RADICL_Tiled_df_mOPC_short_byGeneID[, c("geneID", "kbTile", "Exp_table")]
aa_tab <- lapply(RADICL_Tiled_df_mOPC_short_byGeneID$kbTile, table)
RADICL_Tiled_df_mOPC_short_byGeneID_modified$kbTile_table <- aa_tab
RADICL_Tiled_df_mOPC_short_byGeneID_modified <- RADICL_Tiled_df_mOPC_short_byGeneID_modified[, c("geneID", "Exp_table", "kbTile_table")]
RADICL_Tiled_df_mOPC_short_byGeneID_modified$nu_ints_total <- aa_len
RADICL_Tiled_df_mOPC_short_byGeneID_modified$nu_exp_types <- aa_len_uniq
aa_tile_uniq <- unlist(lapply(lapply(RADICL_Tiled_df_mOPC_short_byGeneID$kbTile, unique), length))
RADICL_Tiled_df_mOPC_short_byGeneID_modified$nu_uniq_tiles <- aa_tile_uniq
aasind <- sort(aa_len, index.return=T, decreasing = T)$ix
RADICL_Tiled_df_mOPC_short_byGeneID_modified_sorted <- RADICL_Tiled_df_mOPC_short_byGeneID_modified[aasind,]



par(mfrow = c(2,1), mar = c(4,4,4,4))
hist(log10(aa_len), main = "histogram:log10 of total number of interactions")
hist(log10(aa_tile_uniq), main = "histogram:log10 of number of unique interacting 1kb tiles")


aa_my_att <- c("chromosome_name", "start_position", "end_position",
               "ensembl_gene_id","external_gene_name", "gene_biotype")

aats4 <- getBM(attributes = aa_my_att, 
               filters = "ensembl_gene_id", 
               values = RADICL_Tiled_df_mOPC_short_byGeneID_modified_sorted$geneID, 
               mart = mouse_mart)

lncRNA_all_info_biomart_mOPC <- aats4
aaa_name <-aats4$external_gene_name[match(RADICL_Tiled_df_mOPC_short_byGeneID_modified_sorted$geneID, aats4$ensembl_gene_id)]
aaa_btyp <-aats4$gene_biotype[match(RADICL_Tiled_df_mOPC_short_byGeneID_modified_sorted$geneID, aats4$ensembl_gene_id)]

RADICL_Tiled_df_mOPC_short_byGeneID_modified_sorted$gene_name <- aaa_name
RADICL_Tiled_df_mOPC_short_byGeneID_modified_sorted$gene_biotype <- aaa_btyp

View(RADICL_Tiled_df_mOPC_short_byGeneID_modified_sorted[1:100,c("nu_ints_total", "nu_exp_types", 
                                                                 "nu_uniq_tiles",
                                                                 "gene_name", "gene_biotype")])


# choosing the lncRNAs:  Ones with more than 1000 unique interacting tiles
sum(RADICL_Tiled_df_mOPC_short_byGeneID_modified_sorted$nu_uniq_tiles > 1000)
save(list = c("RADICL_Tiled_df_mOPC_short_byGeneID_modified_sorted"), 
     file = "RADICL_Tiled_df_mOPC_short_byGeneID_modified_sorted.RData")
# lncRNA_chosen_gt1k_uniqTiles <- RADICL_Tiled_df_mOPC_short_byGeneID_modified_sorted[(RADICL_Tiled_df_mOPC_short_byGeneID_modified_sorted$nu_uniq_tiles > 1000),
#                                                                                     c("geneID", "gene_name", "nu_uniq_tiles", "nu_ints_total", "nu_exp_types","gene_biotype"  )]
# View(lncRNA_chosen_gt1k_uniqTiles)
# lncRNA_chosen_gt1k_uniqTiles <- lncRNA_chosen_gt1k_uniqTiles[-c(which(lncRNA_chosen_gt1k_uniqTiles$gene_name == "4930519F16Rik")),]
# View(lncRNA_chosen_gt1k_uniqTiles)
# 
# aaachr <- lncRNA_all_info_biomart_mOPC[match(lncRNA_chosen_gt1k_uniqTiles$geneID, lncRNA_all_info_biomart_mOPC$ensembl_gene_id), c("chromosome_name", "start_position", "end_position", "external_gene_name")]
# 
# lncRNA_chosen_gt1k_uniqTiles <- cbind(aaachr[, 1:3], lncRNA_chosen_gt1k_uniqTiles)


aa_mutual_link<- lncRNA_chosen_gt1k_uniqTiles$gene_name[which(lncRNA_chosen_gt1k_uniqTiles$geneID %in% RADICL_Tiled_df_mOPC_short_byGeneID_modified_sorted$geneID)] 
#aa_mutual_link <- setdiff(aa_mutual_link, "Trerf1") # removinf Trerf1 since it seems to be  protein coding 
RADICL_Tiled_df_mOPC_short_byGeneID_modified_mutual <-  RADICL_Tiled_df_mOPC_short_byGeneID_modified_sorted[RADICL_Tiled_df_mOPC_short_byGeneID_modified_sorted$gene_name %in% aa_mutual_link,]
RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified_mutual <- RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified_sorted[match(RADICL_Tiled_df_mOPC_short_byGeneID_modified_mutual$gene_name, RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified_sorted$gene_name),]

aaperc <- numeric(nrow(RADICL_Tiled_df_mOPC_short_byGeneID_modified_mutual))

names(aaperc) <- RADICL_Tiled_df_mOPC_short_byGeneID_modified_mutual$gene_name
aaperc2 <- numeric(nrow(RADICL_Tiled_df_mOPC_short_byGeneID_modified_mutual))

names(aaperc2) <- RADICL_Tiled_df_mOPC_short_byGeneID_modified_mutual$gene_name

for(i in 1:nrow(RADICL_Tiled_df_mOPC_short_byGeneID_modified_mutual)){
  print(RADICL_Tiled_df_mOPC_short_byGeneID_modified_mutual$gene_name[i])
  aa1 <- names(RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified_mutual$kbTile_table[[i]])
  aa2 <- names(RADICL_Tiled_df_mOPC_short_byGeneID_modified_mutual$kbTile_table[[i]])
  print("mESC")
  print(length(aa1))
  print("mOPC")
  print(length(aa2))
  print("intersect")
  print(length(intersect(aa1,aa2)))
  print("union")
  print(length(union(aa1,aa2)))
  aaperc[i] <- length(intersect(aa1,aa2))/length(union(aa1,aa2))
  aaperc2[i] <- length(intersect(aa1,aa2))/ min(length(aa1), length(aa2))
  if((length(aa1) < 100) | (length(aa2) < 100)){
    aaperc[i] <- NA
    aaperc2[i] <- NA
  }
  # if((length(aa1) < 100) | (length(aa2) < 100)){
  #   aaperc[i] <- NA
  #   aaperc2[i] <- NA
  # }
  # if(((length(aa1)/length(aa2)) > 4) | ((length(aa2)/length(aa1)) > 4)){
  #   aaperc[i] <- NA
  #   aaperc2[i] <- NA
  # }
}
par(mfrow = c(1,1), mar = c(8,4,4,4))
barplot(sort(aaperc, decreasing = T), las = 2, ylab = "length(intersect)/length(union)")

lncRNA_mESC_vs_mOPC_sameBS_ratio <- aaperc
save(list = c("lncRNA_mESC_vs_mOPC_sameBS_ratio"), file = "lncRNA_mESC_vs_mOPC_sameBS_ratio.RData")



# create a df wit labels for the Partition6 tiles
aa_df_opc_mesc <- list()
for(i in 1:length(RADICL_Tiled_df_mOPC_short_byGeneID_modified_mutual$geneID)){
  aa_df_opc_mesc[[i]] <- Partition_6_random_chunk_cv_df[Partition_6_random_chunk_cv_df$owner %in% RADICL_Tiled_df_mOPC_short_byGeneID_modified_mutual$gene_name[i], c(1:3)]
  colnames(aa_df_opc_mesc[[i]]) <- c("tile_name", "label_mESC",     "owner" )
  aa_df_opc_mesc[[i]]$label_mOPC <- numeric(nrow(aa_df_opc_mesc[[i]]))
  aa_df_opc_mesc[[i]]$label_mOPC[aa_df_opc_mesc[[i]]$tile_name %in% names(RADICL_Tiled_df_mOPC_short_byGeneID_modified_mutual$kbTile_table[[i]])] <- 1
}
names(aa_df_opc_mesc) <- RADICL_Tiled_df_mOPC_short_byGeneID_modified_mutual$gene_name

aa_whlnc <- RADICL_Tiled_df_mOPC_short_byGeneID_modified_mutual$gene_name[unlist(lapply(RADICL_Tiled_df_mOPC_short_byGeneID_modified_mutual$kbTile_table, length)) > 1000]
aa_whlnc2 <- RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified_mutual$gene_name[unlist(lapply(RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified_mutual$kbTile_table, length)) > 1000]
intersect(aa_whlnc, aa_whlnc2)
aa_df_opc_mesc_filt <- aa_df_opc_mesc[which(names(aa_df_opc_mesc) %in% aa_whlnc)]
rbind(RADICL_Tiled_df_mOPC_short_byGeneID_modified_mutual$gene_name, 
      unlist(lapply(RADICL_Tiled_df_mOPC_short_byGeneID_modified_mutual$kbTile_table, length)),
      unlist(lapply(RADICL_GRID_Tiled_df_short_noGRID_byGeneID_modified_mutual$kbTile_table, length)))
##########################################################################################
##########################################################
# CREATE  a heatmap where rows are lncRNAs, columns are 1Mb regions of DNA --> to show they mostly bind near their origin of biosynthesis

RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat
MM9_1kb_tiled
MM9_1kb_tiled_chrname
lncRNA_chosen_gt1k_uniqTiles
aaun <- unique(MM9_1kb_tiled_chrname)
aaun <- aaun[-grep("Random", aaun, ignore.case = T)]
aaun_list <- list()
for(i in 1:length(aaun)){
  print(i)
  print(aaun[i])
  aax <- which(MM9_1kb_tiled_chrname %in% aaun[i])
  aanut <- ceiling(length(aax)/1000)
  aaun_list[[i]] <- list()
  for(j in 1:aanut){
    aamynam <- paste0(aaun[i], "_", c(((j-1)*1000 + 1):(min((j * 1000), (length(aax))))))
    aaun_list[[i]][[j]] <- MM9_1kb_tiled_GR[MM9_1kb_tiled_GR$tile %in% aamynam]
  }
  
}
names(aaun_list[[i]]) <- aaun


aaun_mat_list <- list()
for(i in 1:length(aaun_list)){
  print(i)
  aaun_mat_list[[i]] <- matrix(0L,
                               nrow = nrow(lncRNA_chosen_gt1k_uniqTiles),
                               ncol = length(aaun_list[[i]]))
  rownames(aaun_mat_list[[i]]) <- lncRNA_chosen_gt1k_uniqTiles$gene_name
  for(j in 1:length(aaun_list[[i]])){
    for(k in 1:nrow(lncRNA_chosen_gt1k_uniqTiles)){
      aamykb <- RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat$kbTile[RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat$geneID %in% lncRNA_chosen_gt1k_uniqTiles$geneID[k]]
      
      aaun_mat_list[[i]][k,j] <- length(intersect(aamykb, aaun_list[[i]][[j]]$tile))
    }
  }
}

names(aaun_mat_list) <- aaun
lncBinding_1MBres <- aaun_mat_list
save(list = c("lncBinding_1MBres"),
     file = "~/Documents/Shayan/BioInf/lncRNA/lncBinding_1MBres.RData")

# get the postition of each lncRNA on the same coords
aawco <- numeric(nrow(lncRNA_chosen_gt1k_uniqTiles))
for(i in 1:nrow(lncRNA_chosen_gt1k_uniqTiles)){
  print(i)
  j <- which(names(aaun_list)==lncRNA_chosen_gt1k_uniqTiles$chromosome_name[i])
  if(j == 1){
    aaccnt <- 0
  }else{
    aaccnt <- cumsum(unlist(lapply(aaun_mat_list, ncol)))[j-1]
  }
  
    aafound <- 0
    for(k in 1:length(aaun_list[[j]])){
      aatss <- findOverlaps(lncRNA_chosen_gt1k_uniqTiles_GR[i], aaun_list[[j]][[k]])
      if(length(aatss@from) > 0){
        aawco[i] <- aaccnt + k
        aafound <- 1
        break
      }
    }
}
aawco
aaun_mat_list_all <- do.call(cbind, aaun_mat_list)

dim(aaun_mat_list_all)
# get the coordinates of all miRNAs on the same coords
miRNA_mmu_mm10 <- read.delim("~/Documents/Shayan/BioInf/lncRNA/miRNA/miRNA_mmu_coords_mm10.txt", skip =  13, header = F)[,c(1,4,5,7,3,9)]
colnames(miRNA_mmu_mm10) <- c("chr", "start", "end", "strand", "type", "Desc")
aax <- unlist(lapply(strsplit(miRNA_mmu_mm10$Desc, ";"), "[[", 3))
miRNA_mmu_mm10$name <- gsub(pattern = "Name=mmu-", replacement = "", aax)
miRNA_mmu_mm10_Gr <- makeGRangesFromDataFrame(df = miRNA_mmu_mm10, keep.extra.columns = T)
aa_chain <- import.chain("~/Documents/Shayan/BioInf/Liftover_chain_files/mm10ToMm9.over.chain")
miRNA_mmu_mm9_Gr  <- unlist(liftOver(x = miRNA_mmu_mm10_Gr, chain = aa_chain),use.names = T)
miRNA_mmu_mm9_Gr_primary <- miRNA_mmu_mm9_Gr[miRNA_mmu_mm9_Gr$type == "miRNA_primary_transcript"]
miRNA_mmu_mm9_Gr_mature <- miRNA_mmu_mm9_Gr[miRNA_mmu_mm9_Gr$type == "miRNA"]

amigg <- miRNA_mmu_mm9_Gr_mature[grep(("miR-361-|miR-296|miR-423|miR-92a|miR-331|miR-330|miR-326|miR-24|miR-149|miR-484|miR-339|miR-328|miR-874|miR-22-|miR-214|miR-370"),miRNA_mmu_mm9_Gr_mature$name, ignore.case = T)]
miRNA_mmu_mm9_Gr_mature[grep("miR-296", miRNA_mmu_mm9_Gr_mature$name)]



aawco_mi <- numeric(length(miRNA_mmu_mm9_Gr_mature))
for(i in 1:length(miRNA_mmu_mm9_Gr_mature)){
  print(i)
  j <- which(names(aaun_list)==seqnames(miRNA_mmu_mm9_Gr_mature)[i]@values)
  if(j == 1){
    aaccnt <- 0
  }else{
    aaccnt <- cumsum(unlist(lapply(aaun_mat_list, ncol)))[j-1]
  }
  
  aafound <- 0
  for(k in 1:length(aaun_list[[j]])){
    aatss <- findOverlaps(miRNA_mmu_mm9_Gr_mature[i], aaun_list[[j]][[k]])
    if(length(aatss@from) > 0){
      aawco_mi[i] <- aaccnt + k
      aafound <- 1
      break
    }
  }
}
aawco_mi
names(aawco_mi) <- miRNA_mmu_mm9_Gr_mature$name

save(list = c("aawco_mi"),
     file = "~/Documents/Shayan/BioInf/lncRNA/aawco_mi.RData")

# for each lncRNA get the miRNAs whose genes are located within that lncRNA's interacting tiles
RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat$gene_name <- lncRNA_chosen_gt1k_uniqTiles$gene_name[match(RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat$geneID, 
                                                                                                                        lncRNA_chosen_gt1k_uniqTiles$geneID)]
RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat_conf2 <- RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat[RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat$exp_type_number > 1,]
par(mfrow = c(1,1), mar = c(8,4,4,4))
barplot(rbind(table(RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat$gene_name) - table(RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat_conf2$gene_name),
              table(RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat_conf2$gene_name)), las = 2)

lncRNA_mRNA_int_list <- list()
lncRNA_mRNA_int_list_conf2 <- list()
lncRNA_target_Grange_list <- list()
lncRNA_target_Grange_list_conf2 <- list() 
for(i in 1:nrow(lncRNA_chosen_gt1k_uniqTiles)){
  lncRNA_target_Grange_list[[i]] <- MM9_1kb_tiled_GR[match(RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat$kbTile[RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat$gene_name == lncRNA_chosen_gt1k_uniqTiles$gene_name[i]],
                                                           MM9_1kb_tiled_GR$tile)]
  lncRNA_target_Grange_list_conf2[[i]] <- MM9_1kb_tiled_GR[match(RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat_conf2$kbTile[RADICL_GRID_Tiled_df_short_ChosenLNC_byGeneID__tile_pairs_mat_conf2$gene_name == lncRNA_chosen_gt1k_uniqTiles$gene_name[i]],
                                                                 MM9_1kb_tiled_GR$tile)]
  aaov1 <- findOverlaps(query = miRNA_mmu_mm9_Gr_mature, subject = lncRNA_target_Grange_list[[i]])
  aaov2 <- findOverlaps(query = miRNA_mmu_mm9_Gr_mature, subject = lncRNA_target_Grange_list_conf2[[i]])
  lncRNA_mRNA_int_list[[i]] <- miRNA_mmu_mm9_Gr_mature$name[unique(aaov1@from)]
  lncRNA_mRNA_int_list_conf2[[i]] <- miRNA_mmu_mm9_Gr_mature$name[unique(aaov2@from)]
}

names(lncRNA_mRNA_int_list) <- lncRNA_chosen_gt1k_uniqTiles$gene_name
names(lncRNA_mRNA_int_list_conf2) <- lncRNA_chosen_gt1k_uniqTiles$gene_name

names(lncRNA_target_Grange_list) <- lncRNA_chosen_gt1k_uniqTiles$gene_name
names(lncRNA_target_Grange_list_conf2) <- lncRNA_chosen_gt1k_uniqTiles$gene_name

save(list = c("lncRNA_target_Grange_list", "lncRNA_target_Grange_list_conf2"),
     file = "~/Documents/Shayan/BioInf/lncRNA/lncRNA_target_Grange_list.RData")

save(list = c("lncRNA_mRNA_int_list", "lncRNA_mRNA_int_list_conf2"),
     file = "~/Documents/Shayan/BioInf/lncRNA/lncRNA_mRNA_int_list.RData")

library(BSgenome.Mmusculus.UCSC.mm9)
aa_miRNA_Seq <- getSeq(x = BSgenome.Mmusculus.UCSC.mm9, names = miRNA_mmu_mm9_Gr_mature)
miRNA_mmu_mm9_Gr_mature$Seq <- as.character(aa_miRNA_Seq)
miRNA_mmu_mm9_Gr_mature$Seed <- substr(as.character(aa_miRNA_Seq),2,8)
save(list = c("miRNA_mmu_mm9_Gr_mature"), file ="~/Documents/Shayan/BioInf/lncRNA/miRNA_mmu_mm9_Gr_mature.RData" )

load("~/Documents/Shayan/BioInf/lncRNA/kmer_qq_posneg.RData")
aatstsdf <- data.frame(miR.family = miRNA_mmu_mm9_Gr_mature$name,Seed.m8=miRNA_mmu_mm9_Gr_mature$Seed)
miRNA_seed_lncRNA_qq <- lncRNA_qq_miRNA(miRNA_df = aatstsdf, QQ_df_pos = kmer_qq_posneg$positive, QQ_df_neg = kmer_qq_posneg$negative)

aatstsdf <- data.frame(miR.family = miRNA_mmu_mm9_Gr_mature$name,Seed.m8=miRNA_mmu_mm9_Gr_mature$Seq)
miRNA_full_lncRNA_qq <- lncRNA_qq_miRNA(miRNA_df = aatstsdf, QQ_df_pos = kmer_qq_posneg$positive, QQ_df_neg = kmer_qq_posneg$negative)

save(list = c("miRNA_seed_lncRNA_qq","miRNA_full_lncRNA_qq"), file = "~/Documents/Shayan/BioInf/lncRNA/miRNA_seed_lncRNA_qq.RData")

#aacsep <- cumsum(unlist(lapply(aaun_mat_list, ncol)))

aabreak <- c(0, quantile(aaun_mat_list_all[aaun_mat_list_all > 0], prob = seq(0,1,0.1)))
aacol <- colorRampPalette(brewer.pal(9, "YlOrRd"))(length(aabreak)-1)
aachraa <- character(0)
for(i in 1:length(aaun_mat_list)){
  aachraa <- c(aachraa, rep(names(aaun_mat_list)[i], ncol(aaun_mat_list[[i]])))
}

aaccooll <- colorRampPalette(brewer.pal(12,"Set3"))(length(table(aachraa)))
aaccooll2 <- sample(aaccooll, length(aaccooll), replace = F)
barplot(rep(10,22), col = aaccooll2)
aaColSideColors <- aaccooll2[factor(aachraa)]

aadro <- numeric(ncol(aaun_mat_list_all))
aadro[aawco] <- 200
aaun_mat_list_all_2 <- rbind(aaun_mat_list_all, aadro)
rownames(aaun_mat_list_all_2)[nrow(aaun_mat_list_all_2)] <- "Origin"
aaun_mat_list_all_2_sorted <- aaun_mat_list_all_2[c(sort(aawco, index.return=T)$ix, 29),]
png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/binding_heatmap.png",    # create PNG for the heat map        
    width = 12*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11)
heatmap.2(aaun_mat_list_all_2_sorted, 
          Rowv = F, 
          Colv = F, 
          dendrogram = "none",
          trace = "none", 
          colsep = aacsep,
          sepcolor = "grey",
          sepwidth = c(0.4,0.4),
          breaks = aabreak,
          col = (aacol), 
          labCol = "",
          margins = c(10,10),
          key.title=NA,
          key.xlab=NA,
          keysize=0.5,
          density.info="none"
          #,ColSideColors =  aaColSideColors
          )
dev.off()
####################################################################################################################################################

