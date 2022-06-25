#lncRNA data partitioning
setwd("~/Documents/shayan/BioInf/lncRNA/")
library(zoo)
#####################

# tile features:
# ChIPATLAS_features
# RBP_features
# Partition_1_5mer_mat

# owner features:
# kmer_features_lncRNA28_normalized
# mESC_CHIPATLAS_lncRNA_features
# mESC_RBP_lncRNA_features

library(BSgenome.Mmusculus.UCSC.mm9)
# 
# Distance_based_territory_assignment
# 
# MM9_1kb_tiled
# MM9_1kb_tiled_owner
# MM9_1kb_tiled_owner_labels
# MM9_1kb_tiled_owner_labels_binary
# MM9_1kb_tiled_TADnu
# 
# 
# # lncRNA df
# lncRNA_chosen_gt1k_uniqTiles

MM9_1kb_tiled_owner_3Mfilter <- MM9_1kb_tiled_owner[-c(index_first_3000_perchr)]
MM9_1kb_tiled_owner_labels_3Mfilter <- MM9_1kb_tiled_owner_labels[-c(index_first_3000_perchr)]
MM9_1kb_tiled_owner_labels_binary_3Mfilter <- MM9_1kb_tiled_owner_labels_binary[-c(index_first_3000_perchr)]
MM9_1kb_tiled_TADnu_3Mfilter <- MM9_1kb_tiled_TADnu[-c(index_first_3000_perchr),]


MM9_1kb_tiled_owner_3Mfilter_noNA <- MM9_1kb_tiled_owner_3Mfilter[!is.na(MM9_1kb_tiled_owner_3Mfilter)]
MM9_1kb_tiled_owner_labels_binary_3Mfilter_noNA <- MM9_1kb_tiled_owner_labels_binary_3Mfilter[!is.na(MM9_1kb_tiled_owner_3Mfilter)]
##############
# find out the comulative distribution of - examples in various distances from positive examples for each lncRNA


colSums(territory_assigned_interaction)
colSums(territory_assigned_interaction_label)
colSums(territory_assigned_interaction_label_filtered)


table(MM9_1kb_tiled_chrname[which(is.na(MM9_1kb_tiled_owner))])

territory_assigned_interaction_3Mfilered <- territory_assigned_interaction[-c(index_first_3000_perchr),]
territory_assigned_interaction_label_3Mfiltered <- territory_assigned_interaction_label_filtered[-c(index_first_3000_perchr),]
MM9_1kb_tiled_chrname_3Mfiltered <- MM9_1kb_tiled_chrname[-c(index_first_3000_perchr)]


# for each lnc RNA compute the number and percentage of negative examples that are in (1kb, 10kb, 50kb, 100kb, 200kb, 300kb, 500kb, 1Mb) distance a
#  of a positive example for that lncRNA
aa_dist_Set <- c(1, 10, 50, 100, 200, 300, 500, 1000, 5000)
comulative_distance_neg_cnt <- matrix(nrow = ncol(territory_assigned_interaction_3Mfilered), ncol = length(aa_dist_Set))
rownames(comulative_distance_neg_cnt) <- colnames(territory_assigned_interaction_3Mfilered)
colnames(comulative_distance_neg_cnt) <- aa_dist_Set

comulative_distance_neg_pct <- matrix(nrow = ncol(territory_assigned_interaction_3Mfilered), ncol = length(aa_dist_Set))
rownames(comulative_distance_neg_pct) <- colnames(territory_assigned_interaction_3Mfilered)
colnames(comulative_distance_neg_pct) <- aa_dist_Set

Total_positive_perlncRNA <- numeric(length = ncol(territory_assigned_interaction_3Mfilered))
Total_negative_perlncRNA <- numeric(length = ncol(territory_assigned_interaction_3Mfilered))
names(Total_positive_perlncRNA) <- colnames(territory_assigned_interaction_3Mfilered)
names(Total_negative_perlncRNA) <- colnames(territory_assigned_interaction_3Mfilered)
lncRNA_chromosome <- character(length = ncol(territory_assigned_interaction_3Mfilered))
names(lncRNA_chromosome) <- colnames(territory_assigned_interaction_3Mfilered)
for(aa_lnc in 1:ncol(territory_assigned_interaction_3Mfilered)){
  print(colnames(territory_assigned_interaction_3Mfilered)[aa_lnc])
  aa_ter_members <- which(territory_assigned_interaction_3Mfilered[, aa_lnc] == 1)
  aa_ter_pos <- which(territory_assigned_interaction_label_3Mfiltered[, aa_lnc] == 1)
  Total_positive_perlncRNA[aa_lnc] <- length(aa_ter_pos)
  aa_ter_neg <- setdiff(aa_ter_members, aa_ter_pos)
  Total_negative_perlncRNA[aa_lnc] <- length(aa_ter_neg)
  lncRNA_chromosome[aa_lnc] <- MM9_1kb_tiled_chrname_3Mfiltered[aa_ter_members[1]]
  for(aa_dis in 1:length(aa_dist_Set)){
    print(aa_dist_Set[aa_dis])
    aa_ptn_negb <-  unique(unlist(lapply(aa_ter_pos, function(x, d) c((x-d):(x+d)), aa_dist_Set[aa_dis])))
    aa_ptn_negb <- aa_ptn_negb[aa_ptn_negb > 0]
    comulative_distance_neg_cnt[aa_lnc,aa_dis] <- sum(aa_ptn_negb %in% aa_ter_neg)
    comulative_distance_neg_pct[aa_lnc,aa_dis] <- comulative_distance_neg_cnt[aa_lnc,aa_dis]/length(aa_ter_neg)
  }
}


par(mfrow = c(7,4), mar = c(3,4,4,3))
for(i in 1:nrow(comulative_distance_neg_cnt)){
  barplot(comulative_distance_neg_cnt[i,], las = 2, 
          main =paste0(rownames(comulative_distance_neg_cnt)[i], " ",lncRNA_chromosome[i] ,"\n",
                       "# pos: ", Total_positive_perlncRNA[i], "\n",
                       "# neg: ", Total_negative_perlncRNA[i]))
  
}
par(mfrow = c(7,4), mar = c(3,3,3,3))
for(i in 1:nrow(comulative_distance_neg_cnt)){
  barplot(comulative_distance_neg_pct[i,], las = 2, main = rownames(comulative_distance_neg_cnt)[i], ylim = c(0,1))
  
}
par(mfrow = c(1,1), mar = c(4,4,4,4))

barplot(log10(colSums(comulative_distance_neg_cnt)),
        main = "log10 #total - examples by kb distance from +", ylim = c(0,6), las = 2)
abline(h = seq(4,6,1), col = 2, lwd = 0.6, lty = 2)


####### get the ratio between available pos and neg per lncRNA
par(mar = c(8,4,4,4))
barplot(sort(Total_negative_perlncRNA/Total_positive_perlncRNA, decreasing = T), las = 2, main= '#-/#+')
abline(h = c(1,10,100), col = 2)


# chose at most a 1:10 ratio --> check how many data points does this make
aa_chosen_neg_nu <- numeric(length = length(Total_negative_perlncRNA))

for(i in 1:length(aa_chosen_neg_nu)){
  aa_chosen_neg_nu[i] <- min(Total_negative_perlncRNA[i], Total_positive_perlncRNA[i]*10)
}
names(aa_chosen_neg_nu) <- names(Total_negative_perlncRNA)
boxplot(aa_chosen_neg_nu, las = 2, main= '#-')


# get the distance from each negative to its closest positive
distance_neg_closest_pos <- matrix(nrow = nrow(territory_assigned_interaction_3Mfilered), 
                                   ncol = ncol(territory_assigned_interaction_3Mfilered))
rownames(distance_neg_closest_pos) <- rownames(territory_assigned_interaction_3Mfilered)
colnames(distance_neg_closest_pos) <- colnames(territory_assigned_interaction_3Mfilered)

distance_neg_closest_pos_vec <- numeric(nrow(territory_assigned_interaction_3Mfilered))
names(distance_neg_closest_pos_vec) <-  rownames(territory_assigned_interaction_3Mfilered)
for(aa_lnc in 1:ncol(territory_assigned_interaction_3Mfilered)){
  print(colnames(territory_assigned_interaction_3Mfilered)[aa_lnc])
  aa_ter_members <- which(territory_assigned_interaction_3Mfilered[, aa_lnc] == 1)
  aa_ter_pos <- which(territory_assigned_interaction_label_3Mfiltered[, aa_lnc] == 1)
#  Total_positive_perlncRNA[aa_lnc] <- length(aa_ter_pos)
  aa_ter_neg <- setdiff(aa_ter_members, aa_ter_pos)
#  Total_negative_perlncRNA[aa_lnc] <- length(aa_ter_neg)
  lncRNA_chromosome[aa_lnc] <- MM9_1kb_tiled_chrname_3Mfiltered[aa_ter_members[1]]
  for(aa_neg in 1:length(aa_ter_neg)){
#    print(aa_dist_Set[aa_dis])
    distance_neg_closest_pos[aa_ter_neg[aa_neg],aa_lnc] <- min(abs(aa_ter_neg[aa_neg] - aa_ter_pos))
    distance_neg_closest_pos_vec[aa_ter_neg[aa_neg]] <- distance_neg_closest_pos[aa_ter_neg[aa_neg],aa_lnc]
  }
}






# setup pipeline such that negative examples are sampled with prob inversely proportional to their distance to the closest positive

# divide into train/test/validation --> chuncks --> how?

#in order to decide on training test modeling, i want to see how positive examples are located on each chromosome

aa_tilename <- names(MM9_1kb_tiled_owner_labels_binary_3Mfilter)[MM9_1kb_tiled_owner_labels_binary_3Mfilter == 1]

aa_tilenamesp <- strsplit(aa_tilename, "_")
aa_tilenamesp1<- unlist(lapply(aa_tilenamesp, "[[", 1))
aa_tilenamesp2<- as.numeric(unlist(lapply(aa_tilenamesp, "[[", 2)))
aachr_size <-read.table("~/Documents/Shayan/BioInf/lncRNA/mm9_chr_sizes.txt", header = F, stringsAsFactors = F)
aachr_size2 <- aachr_size[aachr_size$V1 %in% aa_tilenamesp1,]
aachr_vecs <- list()
for(i in 1:nrow(aachr_size2)){
  aachr_vecs[[i]] <- numeric(length = ceiling(aachr_size2$V2[i]/1000))
  names(aachr_vecs[[i]]) <- c(1:length(aachr_vecs[[i]]))
  aachr_vecs[[i]][aa_tilenamesp2[which(aa_tilenamesp1 %in% aachr_size2$V1[i])]] <- 1
}
names(aachr_vecs) <- aachr_size2$V1



for(i in 1:nrow(aachr_size2)){
  aaw <- which(lncRNA_chosen_gt1k_uniqTiles$chromosome_name %in% aachr_size2$V1[i])

  for(j in 1:length(aaw)){
    print(lncRNA_chosen_gt1k_uniqTiles$gene_name[aaw[j]])
    aawch <- c(floor((lncRNA_chosen_gt1k_uniqTiles$start_position[aaw[j]])/1000):ceiling((lncRNA_chosen_gt1k_uniqTiles$end_position[aaw[j]])/1000))
    print(range(aawch))
    aachr_vecs[[i]][aawch] <- 3
    
  }
}


png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/chr_binding_map.png",    # create PNG for the heat map        
    width = 12*300,        # 5 x 300 pixels
    height = 12*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10) 
par(mfrow = c(4,4), mar = c(4,4,4,4))
for(i in 1:length(aachr_vecs)){
  print(i)
  barplot(aachr_vecs[[i]], main = paste(names(aachr_vecs)[i], aachr_size2$V2[i]),
          las = 2, xlab = "genomic position", xaxt = "n", yaxt  = "n", space = 2)
}
dev.off()

################################################################################
# for main training/validation/test choose a number of positives per lncRNA min(3000, Total_positive_perlncRNA[i])
#  For negatives min(10 * min(3000, Total_positive_perlncRNA[i]), Total_negative_perlncRNA)
# Then mark the chosen negative and the chosen positives and plot them --> showing territories
distance_neg_closest_pos_vec
territory_assigned_interaction_3Mfilered
territory_assigned_interaction_label_3Mfiltered
MM9_1kb_tiled_owner_labels_binary_3Mfilter
lncRNA_interaction_label_matrix
#find lncRNA overlapping tiles and remove them from pos and neg sets



set.seed(seed = 1234)
Partition_1 <- list()
aa_ovl <- findOverlaps(query = MM9_1kb_tiled_GR_filtered, subject = lncRNA_chosen_gt1k_uniqTiles_GR)

aa_thresh <- 3000
neg_to_pos_ratio <- 10
for(i in 1:ncol(territory_assigned_interaction_label_3Mfiltered)){
  aa_w_pos <- which(territory_assigned_interaction_label_3Mfiltered[, i] == 1)
  aa_w_pos <- setdiff(aa_w_pos,unique(aa_ovl@from))
  aa_w_Terr <- which(territory_assigned_interaction_3Mfilered[, i] == 1)
  aa_w_neg <- setdiff(aa_w_Terr,aa_w_pos)
  aa_w_neg <- setdiff(aa_w_neg,unique(aa_ovl@from))
  if(length(aa_w_pos) > aa_thresh){
    aa_w_pos_sampled <- sort(sample(x = aa_w_pos, size = aa_thresh, replace = F))
  }else{
    aa_w_pos_sampled <- sort(aa_w_pos)
  }
  if(length(aa_w_neg) > (neg_to_pos_ratio * length(aa_w_pos_sampled))){
    aa_w_neg_sampled <- sort(sample(x = aa_w_neg,
                                    size = (neg_to_pos_ratio * length(aa_w_pos_sampled)),
                                    replace = F,
                                    prob = 1/(distance_neg_closest_pos_vec[aa_w_neg]^(2))))
  }else{
    aa_w_neg_sampled <- aa_w_neg
  }
  Partition_1[[i]] <- list(positive=aa_w_pos_sampled, 
                           negative=aa_w_neg_sampled,
                           extra_positive = setdiff(aa_w_pos,aa_w_pos_sampled), 
                           extra_negative = setdiff(aa_w_neg,aa_w_neg_sampled))
  
}
names(Partition_1) <- colnames(territory_assigned_interaction_label_3Mfiltered)

aa_ps <- lapply(Partition_1, "[[", 1)
aa_ng <- lapply(Partition_1, "[[", 2)
aa_data_cnt <- rbind(unlist(lapply(aa_ps, length)),unlist(lapply(aa_ng, length)))
barplot(aa_data_cnt, las = 2, legend.text = c("+", "-")   ,
        main = "data points per lncRNA",
        args.legend=list(x = "topright", bty = "n", 
                         inset=c(0, 0),
                         cex = 1, 
                         y.intersp = 0.6,
                         x.intersp = 0.5,
                         text.width = 1))

aa_ps2 <- lapply(Partition_1, "[[", 3)
aa_ng2 <- lapply(Partition_1, "[[", 4)
aa_data_cnt2 <- rbind(unlist(lapply(aa_ps2, length)),unlist(lapply(aa_ng2, length)))
barplot(aa_data_cnt2, las = 2, legend.text = c("+", "-")   ,
        main = "extra data points per lncRNA",
        args.legend=list(x = "topright", bty = "n", 
                         inset=c(0, 0),
                         cex = 1, 
                         y.intersp = 0.6,
                         x.intersp = 0.5,
                         text.width = 1))


aachr_vecs2 <- list()
for(i in 1:nrow(aachr_size2)){
  aachr_vecs2[[i]] <-   array(dim =  ceiling(aachr_size2$V2[i]/1000))
  names(aachr_vecs2[[i]]) <- c(1:length(aachr_vecs2[[i]]))
  #aachr_vecs[[i]][aa_tilenamesp2[which(aa_tilenamesp1 %in% aachr_size2$V1[i])]] <- 1
}
names(aachr_vecs2) <- aachr_size2$V1

for(i in 1:length(aachr_vecs2)){
  aaw_lnc <- lncRNA_chosen_gt1k_uniqTiles$gene_name[lncRNA_chosen_gt1k_uniqTiles$chromosome_name  == names(aachr_vecs2)[i]]
  aa_wp <- which(names(Partition_1) %in% aaw_lnc)
  for(j in 1:length(aa_wp)){
    aa_tilename_pos <- rownames(territory_assigned_interaction_3Mfilered)[Partition_1[[aa_wp[j]]]$positive]
    aa_tilename_pos2<- as.numeric(unlist(lapply(strsplit(aa_tilename_pos, "_"), "[[", 2)))
    aa_tilename_neg <- rownames(territory_assigned_interaction_3Mfilered)[Partition_1[[aa_wp[j]]]$negative]
    aa_tilename_neg2<- as.numeric(unlist(lapply(strsplit(aa_tilename_neg, "_"), "[[", 2)))
    
    aachr_vecs2[[i]][aa_tilename_pos2] <- 1.25
    aachr_vecs2[[i]][aa_tilename_neg2] <- 0.25
    
    
  }
}

png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/chr_pos_neg_chosen2.png",       
    width = 40*300,        # 5 x 300 pixels
    height = 20*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10) 
par(mfrow = c(16,1), mar = c(4,4,4,4))
for(i in 1:length(aachr_vecs)){
  print(i)
  plot(aachr_vecs[[i]], main = paste(names(aachr_vecs)[i], aachr_size2$V2[i]),
          las = 2, xlab = "genomic position", xaxt = "n", ylim = c(0,4), pch = 19, cex = 0.2)
  points(aachr_vecs2[[i]], pch = 16, col = 2, cex = 0.2)
}
dev.off()


distance_neg_closest_pos_vec
aa_all_pos <- unlist( lapply(Partition_1, "[[", 1))
aa_all_posextra <- unlist( lapply(Partition_1, "[[", 3))


aa_all_neg <- unlist( lapply(Partition_1, "[[", 2))
aa_all_negextra <- unlist( lapply(Partition_1, "[[", 4))

summary(log10(distance_neg_closest_pos_vec[aa_all_neg]))
summary(log10(distance_neg_closest_pos_vec[aa_all_negextra]))
aa_all <- cbind(c(distance_neg_closest_pos_vec[aa_all_neg],
                  distance_neg_closest_pos_vec[aa_all_negextra]),
                c(rep(1, length(aa_all_neg)), 
                  rep(0, length(aa_all_negextra))))
boxplot(log10(aa_all[, 1]) ~ aa_all[, 2],
        ylab ="distance in log10(kb)",
        main ="Distance between chosen (1) and extra (0) negative points\n and closest positive point")

#########################################################################################################
#########################################################################################################
# making a second partition by 1/exp(Distance) sampling

Partition_2 <- list()
aa_ovl <- findOverlaps(query = MM9_1kb_tiled_GR_filtered, subject = lncRNA_chosen_gt1k_uniqTiles_GR)
set.seed(seed = 12345)

aa_thresh <- 3000
neg_to_pos_ratio <- 10
for(i in 1:ncol(territory_assigned_interaction_label_3Mfiltered)){
  aa_w_pos <- which(territory_assigned_interaction_label_3Mfiltered[, i] == 1)
  aa_w_pos <- setdiff(aa_w_pos,unique(aa_ovl@from))
  aa_w_Terr <- which(territory_assigned_interaction_3Mfilered[, i] == 1)
  aa_w_neg <- setdiff(aa_w_Terr,aa_w_pos)
  aa_w_neg <- setdiff(aa_w_neg,unique(aa_ovl@from))
  if(length(aa_w_pos) > aa_thresh){
    aa_w_pos_sampled <- sort(sample(x = aa_w_pos, size = aa_thresh, replace = F))
  }else{
    aa_w_pos_sampled <- sort(aa_w_pos)
  }
  if(length(aa_w_neg) > (neg_to_pos_ratio * length(aa_w_pos_sampled))){
    aa_exp <- exp(1) + 0.1
    aa_prob_sum <- 0
    while(aa_prob_sum < (neg_to_pos_ratio * length(aa_w_pos_sampled))){
      aa_exp <- aa_exp - 0.1
      aa_prob <- (1/(aa_exp^distance_neg_closest_pos_vec[aa_w_neg]))
      aa_prob_sum <- sum(aa_prob > 0)
      stopifnot(aa_exp > 1)
    }
    aa_w_neg_sampled <- sort(sample(x = aa_w_neg,
                                    size = (neg_to_pos_ratio * length(aa_w_pos_sampled)),
                                    replace = F,
                                    prob = aa_prob))
  }else{
    aa_w_neg_sampled <- aa_w_neg
  }
  Partition_2[[i]] <- list(positive=aa_w_pos_sampled, 
                           negative=aa_w_neg_sampled,
                           extra_positive = setdiff(aa_w_pos,aa_w_pos_sampled), 
                           extra_negative = setdiff(aa_w_neg,aa_w_neg_sampled))
  
}
names(Partition_2) <- colnames(territory_assigned_interaction_label_3Mfiltered)

aa_ps <- lapply(Partition_2, "[[", 1)
aa_ng <- lapply(Partition_2, "[[", 2)
aa_data_cnt <- rbind(unlist(lapply(aa_ps, length)),unlist(lapply(aa_ng, length)))
barplot(aa_data_cnt, las = 2, legend.text = c("+", "-")   ,
        main = "data points per lncRNA",
        args.legend=list(x = "topright", bty = "n", 
                         inset=c(0, 0),
                         cex = 1, 
                         y.intersp = 0.6,
                         x.intersp = 0.5,
                         text.width = 1))

aa_ps2 <- lapply(Partition_2, "[[", 3)
aa_ng2 <- lapply(Partition_2, "[[", 4)
aa_data_cnt2 <- rbind(unlist(lapply(aa_ps2, length)),unlist(lapply(aa_ng2, length)))
barplot(aa_data_cnt2, las = 2, legend.text = c("+", "-")   ,
        main = "extra data points per lncRNA",
        args.legend=list(x = "topright", bty = "n", 
                         inset=c(0, 0),
                         cex = 1, 
                         y.intersp = 0.6,
                         x.intersp = 0.5,
                         text.width = 1))


aachr_vecs2 <- list()
for(i in 1:nrow(aachr_size2)){
  aachr_vecs2[[i]] <-   array(dim =  ceiling(aachr_size2$V2[i]/1000))
  names(aachr_vecs2[[i]]) <- c(1:length(aachr_vecs2[[i]]))
  #aachr_vecs[[i]][aa_tilenamesp2[which(aa_tilenamesp1 %in% aachr_size2$V1[i])]] <- 1
}
names(aachr_vecs2) <- aachr_size2$V1

for(i in 1:length(aachr_vecs2)){
  aaw_lnc <- lncRNA_chosen_gt1k_uniqTiles$gene_name[lncRNA_chosen_gt1k_uniqTiles$chromosome_name  == names(aachr_vecs2)[i]]
  aa_wp <- which(names(Partition_2) %in% aaw_lnc)
  for(j in 1:length(aa_wp)){
    aa_tilename_pos <- rownames(territory_assigned_interaction_3Mfilered)[Partition_2[[aa_wp[j]]]$positive]
    aa_tilename_pos2<- as.numeric(unlist(lapply(strsplit(aa_tilename_pos, "_"), "[[", 2)))
    aa_tilename_neg <- rownames(territory_assigned_interaction_3Mfilered)[Partition_2[[aa_wp[j]]]$negative]
    aa_tilename_neg2<- as.numeric(unlist(lapply(strsplit(aa_tilename_neg, "_"), "[[", 2)))
    
    aachr_vecs2[[i]][aa_tilename_pos2] <- 1.25 + rnorm(n = length(aa_tilename_pos2), mean = 0, sd = 0.01)
    aachr_vecs2[[i]][aa_tilename_neg2] <- 0.25 + rnorm(n = length(aa_tilename_neg2), mean = 0, sd = 0.01)
    
    
  }
}

png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/chr_pos_neg_chosen_partition2.png",       
    width = 40*300,        # 5 x 300 pixels
    height = 30*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10) 
par(mfrow = c(16,1), mar = c(3,4,3,4))
for(i in 1:length(aachr_vecs)){
  print(i)
  plot(aachr_vecs[[i]], main = paste(names(aachr_vecs)[i], aachr_size2$V2[i]),
       las = 2, xlab = "genomic position", xaxt = "n", ylim = c(0,3.1), pch = 19, cex = 0.2)
  points(aachr_vecs2[[i]], pch = 8, col = 2, cex = 0.5)
}
dev.off()

#########################################################################################################





###############
# make pairs of datapoints (tile, owner, label), for positive and negatives 


aa_df_list <- list()
for(i in 1:length(Partition_1)){
  aa_df_pos <- data.frame(tile_name = rownames(territory_assigned_interaction_3Mfilered)[Partition_1[[i]]$positive],
                          owner = rep(names(Partition_1)[i], length(Partition_1[[i]]$positive)), 
                          label = rep(1, length(Partition_1[[i]]$positive)), 
                          dataset = sample(c("train", "valid", "test"), 
                                           size = length(Partition_1[[i]]$positive),
                                           replace = T, 
                                           prob = c(0.7, 0.2, 0.1)))
  aa_df_neg <- data.frame(tile_name = rownames(territory_assigned_interaction_3Mfilered)[Partition_1[[i]]$negative],
                          owner = rep(names(Partition_1)[i], length(Partition_1[[i]]$negative)), 
                          label = rep(0, length(Partition_1[[i]]$negative)),
                          dataset = sample(c("train", "valid", "test"), 
                                           size = length(Partition_1[[i]]$negative),
                                           replace = T, 
                                           prob = c(0.7, 0.2, 0.1))
                          )
  
  aa_df_list[[i]] <- rbind(aa_df_pos, aa_df_neg)
  
}

names(aa_df_list) <- names(Partition_1)

Partition_1_dfs <- do.call(rbind, aa_df_list)
Partition_1_df_list <- aa_df_list
for(i in 1:length(Partition_1_df_list)){
  Partition_1_df_list[[i]]$tile_name <- as.character(levels(Partition_1_df_list[[i]]$tile_name)[as.numeric(Partition_1_df_list[[i]]$tile_name)])
}



aachr_vecs3 <- list()
for(i in 1:nrow(aachr_size2)){
  aachr_vecs3[[i]] <-   array(dim =  ceiling(aachr_size2$V2[i]/1000))
  names(aachr_vecs3[[i]]) <- c(1:length(aachr_vecs3[[i]]))
  #aachr_vecs[[i]][aa_tilenamesp2[which(aa_tilenamesp1 %in% aachr_size2$V1[i])]] <- 1
}
names(aachr_vecs3) <- aachr_size2$V1

aachr_vecs4 <- list()
for(i in 1:nrow(aachr_size2)){
  aachr_vecs4[[i]] <-   array(dim =  ceiling(aachr_size2$V2[i]/1000))
  names(aachr_vecs4[[i]]) <- c(1:length(aachr_vecs4[[i]]))
  #aachr_vecs[[i]][aa_tilenamesp2[which(aa_tilenamesp1 %in% aachr_size2$V1[i])]] <- 1
}
names(aachr_vecs4) <- aachr_size2$V1


for(i in 1:length(aachr_vecs3)){
  print(i)
  aaw_lnc <- lncRNA_chosen_gt1k_uniqTiles$gene_name[lncRNA_chosen_gt1k_uniqTiles$chromosome_name  == names(aachr_vecs3)[i]]
  aa_wp <- which(names(Partition_1_df_list) %in% aaw_lnc)
  for(j in 1:length(aa_wp)){
    aa_tilename_pos_train <- Partition_1_df_list[[aa_wp[j]]]$tile_name[(Partition_1_df_list[[aa_wp[j]]]$dataset == "train") & (Partition_1_df_list[[aa_wp[j]]]$label == 1) ]
    aa_tilename_pos_test <- Partition_1_df_list[[aa_wp[j]]]$tile_name[(Partition_1_df_list[[aa_wp[j]]]$dataset == "test") & (Partition_1_df_list[[aa_wp[j]]]$label == 1) ]
    aa_tilename_pos_valid <- Partition_1_df_list[[aa_wp[j]]]$tile_name[(Partition_1_df_list[[aa_wp[j]]]$dataset == "valid") & (Partition_1_df_list[[aa_wp[j]]]$label == 1) ]
    
    
    aa_tilename_pos_train2<- as.numeric(unlist(lapply(strsplit(aa_tilename_pos_train, "_"), "[[", 2)))
    aa_tilename_pos_test2<- as.numeric(unlist(lapply(strsplit(aa_tilename_pos_test, "_"), "[[", 2)))
    aa_tilename_pos_valid2<- as.numeric(unlist(lapply(strsplit(aa_tilename_pos_valid, "_"), "[[", 2)))
    
    aa_tilename_neg_train <- Partition_1_df_list[[aa_wp[j]]]$tile_name[(Partition_1_df_list[[aa_wp[j]]]$dataset == "train") & (Partition_1_df_list[[aa_wp[j]]]$label == 0) ]
    aa_tilename_neg_test <- Partition_1_df_list[[aa_wp[j]]]$tile_name[(Partition_1_df_list[[aa_wp[j]]]$dataset == "test") & (Partition_1_df_list[[aa_wp[j]]]$label == 0) ]
    aa_tilename_neg_valid <- Partition_1_df_list[[aa_wp[j]]]$tile_name[(Partition_1_df_list[[aa_wp[j]]]$dataset == "valid") & (Partition_1_df_list[[aa_wp[j]]]$label == 0) ]
    
    
    aa_tilename_neg_train2<- as.numeric(unlist(lapply(strsplit(aa_tilename_neg_train, "_"), "[[", 2)))
    aa_tilename_neg_test2<- as.numeric(unlist(lapply(strsplit(aa_tilename_neg_test, "_"), "[[", 2)))
    aa_tilename_neg_valid2<- as.numeric(unlist(lapply(strsplit(aa_tilename_neg_valid, "_"), "[[", 2)))
    

    aachr_vecs3[[i]][aa_tilename_pos_test2] <- 1.8
    aachr_vecs3[[i]][aa_tilename_pos_valid2] <- 1.6
    aachr_vecs3[[i]][aa_tilename_pos_train2] <- 1.4
    
    aachr_vecs3[[i]][aa_tilename_neg_test2] <- 0.8
    aachr_vecs3[[i]][aa_tilename_neg_valid2] <- 0.6
    aachr_vecs3[[i]][aa_tilename_neg_train2] <- 0.4
    
    aachr_vecs4[[i]][aa_tilename_pos_train2] <- 5
    aachr_vecs4[[i]][aa_tilename_pos_test2] <- 2
    aachr_vecs4[[i]][aa_tilename_pos_valid2] <- 3
    
    aachr_vecs4[[i]][aa_tilename_neg_train2] <- 5
    aachr_vecs4[[i]][aa_tilename_neg_test2] <- 2
    aachr_vecs4[[i]][aa_tilename_neg_valid2] <- 3
    
  }
}

png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/chr_pos_neg_chosen_random_partitioning.png",       
    width = 40*300,        # 5 x 300 pixels
    height = 35*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10) 
par(mfrow = c(16,1), mar = c(3,4,3,4))
for(i in 1:length(aachr_vecs)){
  print(i)
  plot(aachr_vecs[[i]], main = paste(names(aachr_vecs)[i], aachr_size2$V2[i]),
       las = 2, xlab = "genomic position", xaxt = "n", ylim = c(0,3.1), pch = 19, cex = 0.3)
  points(aachr_vecs2[[i]], pch = 19, col = 4, cex = 0.3)
  points(aachr_vecs3[[i]], pch = 8, col = aachr_vecs4[[i]], cex = 0.3)
}
dev.off()
##########################################################################################################################################
write.table(Partition_1_dfs$tile_name, 
            row.names = F, 
            col.names = F,
            quote = F, 
            file = "~/Documents/Shayan/BioInf/lncRNA/partition_1_tile_list.txt")
##########################################################################################################################################



aatst <- read.table("chr10_100000.5mer", stringsAsFactors = F)
all_5mers <- sort(generate.kmers(.k = 5))


aax <- c("chr10_100000", "chr10_100001", "chr10_100002")
aa_kmer_list <- list()
for(i in 1:length(aax)){
  aa_kmer_list[[i]] <- read.table(file = paste0(aax[i], ".5mer") ,skip = 1, stringsAsFactors = F)
}

aa_kmer_list2 <- lapply(aa_kmer_list, map_to_all5mer, all_5mers)
aa_kmer_all <- do.call(rbind, aa_kmer_list2)
colnames(aa_kmer_all) <- all_5mers

save(list = c("all_5mers", "map_to_all5mer"),
     file = "~/Documents/Shayan/BioInf/lncRNA/read_merge.RData" )
write.table(aax, row.names = F, col.names = F, quote = F, file = "test_tile_name.txt")
load("test_tile_name.txt_5mer.RData")
nrow(my_tile_kmers)
ncol(my_tile_kmers)

##########################################################################################################################################
##########################################################################################################################################
# Chunk partitioning
Partition_1_chunk <- list()
aa_df_list2 <- list()
aa_nu_start <- 5
aa_test_frac <- 0.1
aa_Valid_frac <- 0.2

for(i in 1:length(Partition_1)){
  # choosing test pos
  print(i)
  aa_nu_test_pos <- floor(length(Partition_1[[i]]$positive)*aa_test_frac)
  aa_nu_vali_pos <- floor(length(Partition_1[[i]]$positive)*aa_Valid_frac)
  
  aa_nu_test_neg <- floor(length(Partition_1[[i]]$negative)*aa_test_frac)
  aa_nu_vali_neg <- floor(length(Partition_1[[i]]$negative)*aa_Valid_frac)
  
  aa_distcheck_test <- 0
  aa_pos_ind_univ <- c(1:(length(Partition_1[[i]]$positive) - ceiling(aa_nu_test_pos/aa_nu_start)))
  while(aa_distcheck_test <= ceiling(aa_nu_test_pos/aa_nu_start)){ # making sure the chosen test positives do not overlap
    aa_test_pos_st <- sample(x = aa_pos_ind_univ,
                             size = aa_nu_start,replace = F)
    aa_distcheck_test <- min(dist(aa_test_pos_st))
  }
  aa_test_pos_index <- unlist(lapply(aa_test_pos_st, function(x) return(c(x:(x + floor(aa_nu_test_pos/aa_nu_start))))))
  stopifnot(sum(duplicated(aa_test_pos_index)) == 0)
  
  # choosing positives for validation set
  
  aa_pos_ind_univ_valid <- setdiff(c(1:(length(Partition_1[[i]]$positive) - ceiling(aa_nu_vali_pos/aa_nu_start))), aa_test_pos_index)
  aa_distcheck_valid <- 0
  aa_check_ovl <- 1
  while((aa_distcheck_valid <= ceiling(aa_nu_vali_pos/aa_nu_start)) | aa_check_ovl > 0 ){ # making sure the chosen valid positives do not overlap
    aa_vali_pos_st <- sample(x = aa_pos_ind_univ_valid,
                             size = aa_nu_start,replace = F)
    aa_distcheck_valid <- min(dist(aa_vali_pos_st))
    aa_valid_pos_index <- unlist(lapply(aa_vali_pos_st, function(x) return(c(x:(x + floor(aa_nu_vali_pos/aa_nu_start))))))
   
    aa_check_ovl <- length(intersect(aa_test_pos_index,aa_valid_pos_index))
    
  }
  stopifnot(length(intersect(aa_test_pos_index,aa_valid_pos_index)) == 0)
  stopifnot(sum(duplicated(aa_valid_pos_index)) == 0)
  # positives for training set:
  aa_train_pos_index <- setdiff(c(1:length(Partition_1[[i]]$positive)), c(aa_test_pos_index,aa_valid_pos_index))
  
  #negatives for test set
  aa_neg_univ_ind <- c(1:length(Partition_1[[i]]$negative))
  aa_neg_univ_pos_test_dis <- numeric(length(aa_neg_univ_ind))
  for(aa_neg in 1:length(aa_neg_univ_ind)){
    # calc distance for  negatives to test positives
    aa_neg_univ_pos_test_dis[aa_neg] <- min(abs(Partition_1[[i]]$negative[aa_neg_univ_ind[aa_neg]] - Partition_1[[i]]$positive[aa_test_pos_index]))
  }
  aa_exp <- exp(1) + 0.1
  aa_prob_sum <- 0
  while(aa_prob_sum < aa_nu_test_neg){
    aa_exp <- aa_exp - 0.1
    aa_prob <- (1/(aa_exp^aa_neg_univ_pos_test_dis))
    aa_prob_sum <- sum(aa_prob > 0)
    stopifnot(aa_exp > 1)
  }
  aa_test_neg_ind <- sample(x = aa_neg_univ_ind,
                           size = aa_nu_test_neg,
                           prob = aa_prob,
                           replace = F)
  aa_neg_univ_ind_valid <- setdiff(aa_neg_univ_ind, aa_test_neg_ind)
  aa_neg_univ_pos_valid_dis <- numeric(length(aa_neg_univ_ind_valid))
  for(aa_neg in 1:length(aa_neg_univ_ind_valid)){
    # calc distance for remaining negatives to valid positives
    aa_neg_univ_pos_valid_dis[aa_neg] <- min(abs(Partition_1[[i]]$negative[aa_neg_univ_ind_valid[aa_neg]] - Partition_1[[i]]$positive[aa_valid_pos_index]))
  }
  aa_exp <- exp(1) + 0.1
  aa_prob_sum <- 0
  while(aa_prob_sum < aa_nu_vali_neg){
    aa_exp <- aa_exp - 0.1
    aa_prob <- (1/(aa_exp^aa_neg_univ_pos_valid_dis))
    aa_prob_sum <- sum(aa_prob > 0)
    stopifnot(aa_exp > 1)
  }
 

  aa_valid_neg_ind <- sample(x = aa_neg_univ_ind_valid,
                            size = aa_nu_vali_neg,
                            prob = aa_prob,
                            replace = F)
  aa_train_neg_ind <- setdiff(aa_neg_univ_ind, 
                              c(aa_test_neg_ind, aa_valid_neg_ind))
  # forming the dataframes
  aa_Sets <- c("train", "valid", "test")
  aa_pos_frame_set <- character(length = length(Partition_1[[i]]$positive))
  aa_neg_frame_set <- character(length = length(Partition_1[[i]]$negative))
  
  aa_pos_frame_set[aa_test_pos_index] <- aa_Sets[3]
  aa_pos_frame_set[aa_valid_pos_index] <- aa_Sets[2]
  aa_pos_frame_set[aa_train_pos_index] <- aa_Sets[1]
  aa_neg_frame_set[aa_test_neg_ind] <- aa_Sets[3]
  aa_neg_frame_set[aa_valid_neg_ind] <- aa_Sets[2]
  aa_neg_frame_set[aa_train_neg_ind] <- aa_Sets[1]
  
  
  
  
  aa_df_pos <- data.frame(tile_name = rownames(territory_assigned_interaction_3Mfilered)[Partition_1[[i]]$positive],
                          owner = rep(names(Partition_1)[i], length(Partition_1[[i]]$positive)), 
                          label = rep(1, length(Partition_1[[i]]$positive)), 
                          dataset = aa_pos_frame_set)
  aa_df_neg <- data.frame(tile_name = rownames(territory_assigned_interaction_3Mfilered)[Partition_1[[i]]$negative],
                          owner = rep(names(Partition_1)[i], length(Partition_1[[i]]$negative)), 
                          label = rep(0, length(Partition_1[[i]]$negative)),
                          dataset = aa_neg_frame_set)
  
  aa_df_list2[[i]] <- rbind(aa_df_pos, aa_df_neg)
  
}
Partition_1_chunk <- aa_df_list2
names(Partition_1_chunk) <- names(Partition_1)
Partition_1_chunk_dfs <- do.call(rbind, Partition_1_chunk)
table(Partition_1_chunk$Malat1$dataset)

for(i in 1:length(Partition_1_chunk)){
  Partition_1_chunk[[i]]$tile_name <- as.character(levels(Partition_1_chunk[[i]]$tile_name)[as.numeric(Partition_1_chunk[[i]]$tile_name)])
}

# check if tiles are sorted correctly by their numeric part
aatts <- rownames(territory_assigned_interaction_3Mfilered)[Partition_1[[i]]$negative]

aanu <- as.numeric(unlist(lapply(strsplit(aatts, "_"), "[[", 2)))
identical(aanu, sort(aanu)) # True -> hence I can use the indecis for distance




aachr_vecs3 <- list()
for(i in 1:nrow(aachr_size2)){
  aachr_vecs3[[i]] <-   array(dim =  ceiling(aachr_size2$V2[i]/1000))
  names(aachr_vecs3[[i]]) <- c(1:length(aachr_vecs3[[i]]))
  #aachr_vecs[[i]][aa_tilenamesp2[which(aa_tilenamesp1 %in% aachr_size2$V1[i])]] <- 1
}
names(aachr_vecs3) <- aachr_size2$V1

aachr_vecs4 <- list()
for(i in 1:nrow(aachr_size2)){
  aachr_vecs4[[i]] <-   array(dim =  ceiling(aachr_size2$V2[i]/1000))
  names(aachr_vecs4[[i]]) <- c(1:length(aachr_vecs4[[i]]))
  #aachr_vecs[[i]][aa_tilenamesp2[which(aa_tilenamesp1 %in% aachr_size2$V1[i])]] <- 1
}
names(aachr_vecs4) <- aachr_size2$V1


for(i in 1:length(aachr_vecs3)){
  print(i)
  aaw_lnc <- lncRNA_chosen_gt1k_uniqTiles$gene_name[lncRNA_chosen_gt1k_uniqTiles$chromosome_name  == names(aachr_vecs3)[i]]
  aa_wp <- which(names(Partition_1_chunk) %in% aaw_lnc)
  for(j in 1:length(aa_wp)){
    aa_tilename_pos_train <- Partition_1_chunk[[aa_wp[j]]]$tile_name[(Partition_1_chunk[[aa_wp[j]]]$dataset == "train") & (Partition_1_chunk[[aa_wp[j]]]$label == 1) ]
    aa_tilename_pos_test <- Partition_1_chunk[[aa_wp[j]]]$tile_name[(Partition_1_chunk[[aa_wp[j]]]$dataset == "test") & (Partition_1_chunk[[aa_wp[j]]]$label == 1) ]
    aa_tilename_pos_valid <- Partition_1_chunk[[aa_wp[j]]]$tile_name[(Partition_1_chunk[[aa_wp[j]]]$dataset == "valid") & (Partition_1_chunk[[aa_wp[j]]]$label == 1) ]
    
    
    aa_tilename_pos_train2<- as.numeric(unlist(lapply(strsplit(aa_tilename_pos_train, "_"), "[[", 2)))
    aa_tilename_pos_test2<- as.numeric(unlist(lapply(strsplit(aa_tilename_pos_test, "_"), "[[", 2)))
    aa_tilename_pos_valid2<- as.numeric(unlist(lapply(strsplit(aa_tilename_pos_valid, "_"), "[[", 2)))
    
    aa_tilename_neg_train <- Partition_1_chunk[[aa_wp[j]]]$tile_name[(Partition_1_chunk[[aa_wp[j]]]$dataset == "train") & (Partition_1_chunk[[aa_wp[j]]]$label == 0) ]
    aa_tilename_neg_test <- Partition_1_chunk[[aa_wp[j]]]$tile_name[(Partition_1_chunk[[aa_wp[j]]]$dataset == "test") & (Partition_1_chunk[[aa_wp[j]]]$label == 0) ]
    aa_tilename_neg_valid <- Partition_1_chunk[[aa_wp[j]]]$tile_name[(Partition_1_chunk[[aa_wp[j]]]$dataset == "valid") & (Partition_1_chunk[[aa_wp[j]]]$label == 0) ]
    
    
    aa_tilename_neg_train2<- as.numeric(unlist(lapply(strsplit(aa_tilename_neg_train, "_"), "[[", 2)))
    aa_tilename_neg_test2<- as.numeric(unlist(lapply(strsplit(aa_tilename_neg_test, "_"), "[[", 2)))
    aa_tilename_neg_valid2<- as.numeric(unlist(lapply(strsplit(aa_tilename_neg_valid, "_"), "[[", 2)))
    
    
    aachr_vecs3[[i]][aa_tilename_pos_test2] <- 1.8
    aachr_vecs3[[i]][aa_tilename_pos_valid2] <- 1.6
    aachr_vecs3[[i]][aa_tilename_pos_train2] <- 1.4
    
    aachr_vecs3[[i]][aa_tilename_neg_test2] <- 0.8
    aachr_vecs3[[i]][aa_tilename_neg_valid2] <- 0.6
    aachr_vecs3[[i]][aa_tilename_neg_train2] <- 0.4
    
    
    aachr_vecs4[[i]][aa_tilename_pos_train2] <- 5
    aachr_vecs4[[i]][aa_tilename_pos_test2] <- 2
    aachr_vecs4[[i]][aa_tilename_pos_valid2] <- 3
    
    aachr_vecs4[[i]][aa_tilename_neg_train2] <- 5
    aachr_vecs4[[i]][aa_tilename_neg_test2] <- 2
    aachr_vecs4[[i]][aa_tilename_neg_valid2] <- 3
    
  }
}

png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/chr_pos_neg_chosen_chunk_partitioning.png",       
    width = 40*300,        # 5 x 300 pixels
    height = 35*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10) 
par(mfrow = c(16,1), mar = c(2,4,3,4))
for(i in 1:length(aachr_vecs)){
  print(i)
  plot(aachr_vecs[[i]], main = paste(names(aachr_vecs)[i], aachr_size2$V2[i]),
       las = 2, xlab = "genomic position", xaxt = "n", ylim = c(0,3.1), pch = 19, cex = 0.3)
  points(aachr_vecs2[[i]], pch = 19, col = 4, cex = 0.3)
  points(aachr_vecs3[[i]], pch = 8, col = aachr_vecs4[[i]], cex = 0.3)
}
dev.off()

#################################################################################################################################################################
#################################################################################################################################################################
# Done with two modes of partitioning : random, chunk
# random: Partition_1_dfs
# chunk: Partition_1_chunk_dfs

# create the feature, label matrix for the whole dataset
load("~/Documents/Shayan/BioInf/lncRNA/ChIPATLAS_features.RData")
load("~/Documents/Shayan/BioInf/lncRNA/RBP_features.RData")

Triplex_feaures_partition1
Partition_1_5mer_mat
ChIPATLAS_features_partition1 <- ChIPATLAS_features[match(Partition_1_dfs$tile_name, rownames(ChIPATLAS_features)),]
RBP_features_partition1 <- RBP_features[match(Partition_1_dfs$tile_name, rownames(RBP_features)),]
ChIPATLAS_features_partition1[is.na(ChIPATLAS_features_partition1)] <- 0
RBP_features_partition1[is.na(RBP_features_partition1)] <- 0


mESC_CHIPATLAS_lncRNA_features
kmer_features_lncRNA28_normalized
mESC_RBP_lncRNA_features

# expression RNAseq tile
mESC_Renlab_RNASeq_partition1  <- numeric(nrow(Partition_1_dfs))
names(mESC_Renlab_RNASeq_partition1) <- Partition_1_dfs$tile_name
aamESC_Renlab_tiled_lisRNAseq <-  mESC_Renlab_tiled_list$`GSM723776_RenLab-RNA-Seq-mESC-ZY28.bed`[mESC_Renlab_tiled_list$`GSM723776_RenLab-RNA-Seq-mESC-ZY28.bed`$tile %in% Partition_1_dfs$tile_name,]
mESC_Renlab_RNASeq_partition1[match(aamESC_Renlab_tiled_lisRNAseq$tile, Partition_1_dfs$tile_name)] <- aamESC_Renlab_tiled_lisRNAseq$`aadf_agg$score`

# expression cage tile
aamESC_CAGE_mm9_tiled <- mESC_CAGE_mm9_tiled[mESC_CAGE_mm9_tiled$tile %in% Partition_1_dfs$tile_name,]
mESC_CAGE_mm9_tiled_partition1 <- numeric(nrow(Partition_1_dfs))
mESC_CAGE_mm9_tiled_partition1[match(aamESC_CAGE_mm9_tiled$tile, Partition_1_dfs$tile_name)] <- aamESC_CAGE_mm9_tiled$`aadf_agg$CAGE_score`
names(mESC_CAGE_mm9_tiled_partition1) <- Partition_1_dfs$tile_name

# expression RNAseq lncRNA
mESC_RNAseq_lncRNA_features
mESC_RNAseq_lncRNA_features_full <- numeric(length = nrow(lncRNA_chosen_gt1k_uniqTiles))
names(mESC_RNAseq_lncRNA_features_full) <- lncRNA_chosen_gt1k_uniqTiles$gene_name
mESC_RNAseq_lncRNA_features_full[match(names(mESC_RNAseq_lncRNA_features), names(mESC_RNAseq_lncRNA_features_full))] <- mESC_RNAseq_lncRNA_features



# expression cage lncRNA
mESC_CAGE_lncRNA_features
mESC_CAGE_lncRNA_features_full <- numeric(length = nrow(lncRNA_chosen_gt1k_uniqTiles))
names(mESC_CAGE_lncRNA_features_full) <- lncRNA_chosen_gt1k_uniqTiles$gene_name
mESC_CAGE_lncRNA_features_full[match(names(mESC_CAGE_lncRNA_features), names(mESC_CAGE_lncRNA_features_full))] <- mESC_CAGE_lncRNA_features



Partition_1_dfs$tile_name <- levels(Partition_1_dfs$tile_name)[as.numeric(Partition_1_dfs$tile_name)]
Partition_1_chunk_dfs$tile_name <- levels(Partition_1_chunk_dfs$tile_name)[as.numeric(Partition_1_chunk_dfs$tile_name)]

# 
Partition_1_feature_mat_tile <- cbind(ChIPATLAS_features_partition1, # remember to add repeats
                                      Partition_1_5mer_mat,
                                      RBP_features_partition1,
                                      mESC_CAGE_mm9_tiled_partition1,
                                      mESC_CAGE_mm9_tiled_partition1)

Partition_1_feature_mat_owner <- cbind(mESC_CHIPATLAS_lncRNA_features[match(Partition_1_dfs$owner, rownames(mESC_CHIPATLAS_lncRNA_features)),],
                                       kmer_features_lncRNA28_normalized[match(Partition_1_dfs$owner, rownames(kmer_features_lncRNA28_normalized)),],
                                       mESC_RBP_lncRNA_features[match(Partition_1_dfs$owner, rownames(mESC_RBP_lncRNA_features)),],
                                       mESC_RNAseq_lncRNA_features_full[match(Partition_1_dfs$owner, names(mESC_RNAseq_lncRNA_features_full))],
                                       mESC_CAGE_lncRNA_features_full[match(Partition_1_dfs$owner, names(mESC_CAGE_lncRNA_features_full))])

Partition_1_feature_mat_pair <- cbind(distance_feature_partition_1,
                                      Triplex_feaures_partition1[match(Partition_1_dfs$tile_name , Triplex_feaures_partition1$tile_name), 3])
colnames(Partition_1_feature_mat_pair) <- c("distance", "triplex")

# how to remove - from names
aatst <- colnames(ChIPATLAS_features_partition1)
aatst2 <- gsub(pattern = "-", replacement = "", x = aatst)




##########################
# repeating the random and chunk partitionings for Partition_2

###############
# make pairs of datapoints (tile, owner, label), for positive and negatives 


aa_df_list <- list()
for(i in 1:length(Partition_2)){
  aa_df_pos <- data.frame(tile_name = rownames(territory_assigned_interaction_3Mfilered)[Partition_2[[i]]$positive],
                          owner = rep(names(Partition_2)[i], length(Partition_2[[i]]$positive)), 
                          label = rep(1, length(Partition_2[[i]]$positive)), 
                          dataset = sample(c("train", "valid", "test"), 
                                           size = length(Partition_2[[i]]$positive),
                                           replace = T, 
                                           prob = c(0.7, 0.2, 0.1)))
  aa_df_neg <- data.frame(tile_name = rownames(territory_assigned_interaction_3Mfilered)[Partition_2[[i]]$negative],
                          owner = rep(names(Partition_2)[i], length(Partition_2[[i]]$negative)), 
                          label = rep(0, length(Partition_2[[i]]$negative)),
                          dataset = sample(c("train", "valid", "test"), 
                                           size = length(Partition_2[[i]]$negative),
                                           replace = T, 
                                           prob = c(0.7, 0.2, 0.1))
  )
  
  aa_df_list[[i]] <- rbind(aa_df_pos, aa_df_neg)
  
}

names(aa_df_list) <- names(Partition_2)

Partition_2_dfs <- do.call(rbind, aa_df_list)
Partition_2_df_list <- aa_df_list
for(i in 1:length(Partition_2_df_list)){
  Partition_2_df_list[[i]]$tile_name <- as.character(levels(Partition_2_df_list[[i]]$tile_name)[as.numeric(Partition_2_df_list[[i]]$tile_name)])
}



aachr_vecs3 <- list()
for(i in 1:nrow(aachr_size2)){
  aachr_vecs3[[i]] <-   array(dim =  ceiling(aachr_size2$V2[i]/1000))
  names(aachr_vecs3[[i]]) <- c(1:length(aachr_vecs3[[i]]))
  #aachr_vecs[[i]][aa_tilenamesp2[which(aa_tilenamesp1 %in% aachr_size2$V1[i])]] <- 1
}
names(aachr_vecs3) <- aachr_size2$V1

aachr_vecs4 <- list()
for(i in 1:nrow(aachr_size2)){
  aachr_vecs4[[i]] <-   array(dim =  ceiling(aachr_size2$V2[i]/1000))
  names(aachr_vecs4[[i]]) <- c(1:length(aachr_vecs4[[i]]))
  #aachr_vecs[[i]][aa_tilenamesp2[which(aa_tilenamesp1 %in% aachr_size2$V1[i])]] <- 1
}
names(aachr_vecs4) <- aachr_size2$V1


for(i in 1:length(aachr_vecs3)){
  print(i)
  aaw_lnc <- lncRNA_chosen_gt1k_uniqTiles$gene_name[lncRNA_chosen_gt1k_uniqTiles$chromosome_name  == names(aachr_vecs3)[i]]
  aa_wp <- which(names(Partition_2_df_list) %in% aaw_lnc)
  for(j in 1:length(aa_wp)){
    aa_tilename_pos_train <- Partition_2_df_list[[aa_wp[j]]]$tile_name[(Partition_2_df_list[[aa_wp[j]]]$dataset == "train") & (Partition_2_df_list[[aa_wp[j]]]$label == 1) ]
    aa_tilename_pos_test <- Partition_2_df_list[[aa_wp[j]]]$tile_name[(Partition_2_df_list[[aa_wp[j]]]$dataset == "test") & (Partition_2_df_list[[aa_wp[j]]]$label == 1) ]
    aa_tilename_pos_valid <- Partition_2_df_list[[aa_wp[j]]]$tile_name[(Partition_2_df_list[[aa_wp[j]]]$dataset == "valid") & (Partition_2_df_list[[aa_wp[j]]]$label == 1) ]
    
    
    aa_tilename_pos_train2<- as.numeric(unlist(lapply(strsplit(aa_tilename_pos_train, "_"), "[[", 2)))
    aa_tilename_pos_test2<- as.numeric(unlist(lapply(strsplit(aa_tilename_pos_test, "_"), "[[", 2)))
    aa_tilename_pos_valid2<- as.numeric(unlist(lapply(strsplit(aa_tilename_pos_valid, "_"), "[[", 2)))
    
    aa_tilename_neg_train <- Partition_2_df_list[[aa_wp[j]]]$tile_name[(Partition_2_df_list[[aa_wp[j]]]$dataset == "train") & (Partition_2_df_list[[aa_wp[j]]]$label == 0) ]
    aa_tilename_neg_test <- Partition_2_df_list[[aa_wp[j]]]$tile_name[(Partition_2_df_list[[aa_wp[j]]]$dataset == "test") & (Partition_2_df_list[[aa_wp[j]]]$label == 0) ]
    aa_tilename_neg_valid <- Partition_2_df_list[[aa_wp[j]]]$tile_name[(Partition_2_df_list[[aa_wp[j]]]$dataset == "valid") & (Partition_2_df_list[[aa_wp[j]]]$label == 0) ]
    
    
    aa_tilename_neg_train2<- as.numeric(unlist(lapply(strsplit(aa_tilename_neg_train, "_"), "[[", 2)))
    aa_tilename_neg_test2<- as.numeric(unlist(lapply(strsplit(aa_tilename_neg_test, "_"), "[[", 2)))
    aa_tilename_neg_valid2<- as.numeric(unlist(lapply(strsplit(aa_tilename_neg_valid, "_"), "[[", 2)))
    
    
    aachr_vecs3[[i]][aa_tilename_pos_test2] <- 1.8 + rnorm(n = length(aa_tilename_pos_test2), mean = 0, sd = 0.01)
    aachr_vecs3[[i]][aa_tilename_pos_valid2] <- 1.6 + rnorm(n = length(aa_tilename_pos_valid2), mean = 0, sd = 0.01)
    aachr_vecs3[[i]][aa_tilename_pos_train2] <- 1.4 + rnorm(n = length(aa_tilename_pos_train2), mean = 0, sd = 0.01)
    
    aachr_vecs3[[i]][aa_tilename_neg_test2] <- 0.8 + rnorm(n = length(aa_tilename_neg_test2), mean = 0, sd = 0.01)
    aachr_vecs3[[i]][aa_tilename_neg_valid2] <- 0.6 + rnorm(n = length(aa_tilename_neg_valid2), mean = 0, sd = 0.01)
    aachr_vecs3[[i]][aa_tilename_neg_train2] <- 0.4 + rnorm(n = length(aa_tilename_neg_train2), mean = 0, sd = 0.01)
    
    aachr_vecs4[[i]][aa_tilename_pos_train2] <- 5
    aachr_vecs4[[i]][aa_tilename_pos_test2] <- 2
    aachr_vecs4[[i]][aa_tilename_pos_valid2] <- 3
    
    aachr_vecs4[[i]][aa_tilename_neg_train2] <- 5
    aachr_vecs4[[i]][aa_tilename_neg_test2] <- 2
    aachr_vecs4[[i]][aa_tilename_neg_valid2] <- 3
    
  }
}

png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/chr_pos_neg_chosen_partition2_random_partitioning.png",       
    width = 40*300,        # 5 x 300 pixels
    height = 35*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10) 
par(mfrow = c(16,1), mar = c(3,4,3,4))
for(i in 1:length(aachr_vecs)){
  print(i)
  plot(aachr_vecs[[i]], main = paste(names(aachr_vecs)[i], aachr_size2$V2[i]),
       las = 2, xlab = "genomic position", xaxt = "n", ylim = c(0,3.1), pch = 19, cex = 0.3)
  points(aachr_vecs2[[i]], pch = 19, col = 4, cex = 0.3)
  points(aachr_vecs3[[i]], pch = 8, col = aachr_vecs4[[i]], cex = 0.3)
}
dev.off()
##########################################################################################################################################
write.table(Partition_2_dfs$tile_name, 
            row.names = F, 
            col.names = F,
            quote = F, 
            file = "~/Documents/Shayan/BioInf/lncRNA/partition_2_tile_list.txt")
##########################################################################################################################################
##########################################################################################################################################
# Chunk partitioning
Partition_2_chunk <- list()
aa_df_list2 <- list()
aa_nu_start <- 5
aa_test_frac <- 0.1
aa_Valid_frac <- 0.2
set.seed(seed = 123456)
for(i in 1:length(Partition_2)){
  # choosing test pos
  print(i)
  aa_nu_test_pos <- floor(length(Partition_2[[i]]$positive)*aa_test_frac)
  aa_nu_vali_pos <- floor(length(Partition_2[[i]]$positive)*aa_Valid_frac)
  
  aa_nu_test_neg <- floor(length(Partition_2[[i]]$negative)*aa_test_frac)
  aa_nu_vali_neg <- floor(length(Partition_2[[i]]$negative)*aa_Valid_frac)
  
  aa_distcheck_test <- 0
  aa_pos_ind_univ <- c(1:(length(Partition_2[[i]]$positive) - ceiling(aa_nu_test_pos/aa_nu_start)))
  while(aa_distcheck_test <= ceiling(aa_nu_test_pos/aa_nu_start)){ # making sure the chosen test positives do not overlap
    aa_test_pos_st <- sample(x = aa_pos_ind_univ,
                             size = aa_nu_start,replace = F)
    aa_distcheck_test <- min(dist(aa_test_pos_st))
  }
  aa_test_pos_index <- unlist(lapply(aa_test_pos_st, function(x) return(c(x:(x + floor(aa_nu_test_pos/aa_nu_start))))))
  stopifnot(sum(duplicated(aa_test_pos_index)) == 0)
  
  # choosing positives for validation set
  
  aa_pos_ind_univ_valid <- setdiff(c(1:(length(Partition_2[[i]]$positive) - ceiling(aa_nu_vali_pos/aa_nu_start))), aa_test_pos_index)
  aa_distcheck_valid <- 0
  aa_check_ovl <- 1
  while((aa_distcheck_valid <= ceiling(aa_nu_vali_pos/aa_nu_start)) | aa_check_ovl > 0 ){ # making sure the chosen valid positives do not overlap
    aa_vali_pos_st <- sample(x = aa_pos_ind_univ_valid,
                             size = aa_nu_start,replace = F)
    aa_distcheck_valid <- min(dist(aa_vali_pos_st))
    aa_valid_pos_index <- unlist(lapply(aa_vali_pos_st, function(x) return(c(x:(x + floor(aa_nu_vali_pos/aa_nu_start))))))
    
    aa_check_ovl <- length(intersect(aa_test_pos_index,aa_valid_pos_index))
    
  }
  stopifnot(length(intersect(aa_test_pos_index,aa_valid_pos_index)) == 0)
  stopifnot(sum(duplicated(aa_valid_pos_index)) == 0)
  # positives for training set:
  aa_train_pos_index <- setdiff(c(1:length(Partition_2[[i]]$positive)), c(aa_test_pos_index,aa_valid_pos_index))
  
  #negatives for test set
  aa_neg_univ_ind <- c(1:length(Partition_2[[i]]$negative))
  aa_neg_univ_pos_test_dis <- numeric(length(aa_neg_univ_ind))
  for(aa_neg in 1:length(aa_neg_univ_ind)){
    # calc distance for  negatives to test positives
    aa_neg_univ_pos_test_dis[aa_neg] <- min(abs(Partition_2[[i]]$negative[aa_neg_univ_ind[aa_neg]] - Partition_2[[i]]$positive[aa_test_pos_index]))
  }
  aa_exp <- exp(1) + 0.1
  aa_prob_sum <- 0
  while(aa_prob_sum < aa_nu_test_neg){
    aa_exp <- aa_exp - 0.1
    aa_prob <- (1/(aa_exp^aa_neg_univ_pos_test_dis))
    aa_prob_sum <- sum(aa_prob > 0)
    stopifnot(aa_exp > 1)
  }
  aa_test_neg_ind <- sample(x = aa_neg_univ_ind,
                            size = aa_nu_test_neg,
                            prob = aa_prob,
                            replace = F)
  aa_neg_univ_ind_valid <- setdiff(aa_neg_univ_ind, aa_test_neg_ind)
  aa_neg_univ_pos_valid_dis <- numeric(length(aa_neg_univ_ind_valid))
  for(aa_neg in 1:length(aa_neg_univ_ind_valid)){
    # calc distance for remaining negatives to valid positives
    aa_neg_univ_pos_valid_dis[aa_neg] <- min(abs(Partition_2[[i]]$negative[aa_neg_univ_ind_valid[aa_neg]] - Partition_2[[i]]$positive[aa_valid_pos_index]))
  }
  aa_exp <- exp(1) + 0.1
  aa_prob_sum <- 0
  while(aa_prob_sum < aa_nu_vali_neg){
    aa_exp <- aa_exp - 0.1
    aa_prob <- (1/(aa_exp^aa_neg_univ_pos_valid_dis))
    aa_prob_sum <- sum(aa_prob > 0)
    stopifnot(aa_exp > 1)
  }
  
  
  aa_valid_neg_ind <- sample(x = aa_neg_univ_ind_valid,
                             size = aa_nu_vali_neg,
                             prob = aa_prob,
                             replace = F)
  aa_train_neg_ind <- setdiff(aa_neg_univ_ind, 
                              c(aa_test_neg_ind, aa_valid_neg_ind))
  # forming the dataframes
  aa_Sets <- c("train", "valid", "test")
  aa_pos_frame_set <- character(length = length(Partition_2[[i]]$positive))
  aa_neg_frame_set <- character(length = length(Partition_2[[i]]$negative))
  
  aa_pos_frame_set[aa_test_pos_index] <- aa_Sets[3]
  aa_pos_frame_set[aa_valid_pos_index] <- aa_Sets[2]
  aa_pos_frame_set[aa_train_pos_index] <- aa_Sets[1]
  aa_neg_frame_set[aa_test_neg_ind] <- aa_Sets[3]
  aa_neg_frame_set[aa_valid_neg_ind] <- aa_Sets[2]
  aa_neg_frame_set[aa_train_neg_ind] <- aa_Sets[1]
  
  
  
  
  aa_df_pos <- data.frame(tile_name = rownames(territory_assigned_interaction_3Mfilered)[Partition_2[[i]]$positive],
                          owner = rep(names(Partition_2)[i], length(Partition_2[[i]]$positive)), 
                          label = rep(1, length(Partition_2[[i]]$positive)), 
                          dataset = aa_pos_frame_set)
  aa_df_neg <- data.frame(tile_name = rownames(territory_assigned_interaction_3Mfilered)[Partition_2[[i]]$negative],
                          owner = rep(names(Partition_2)[i], length(Partition_2[[i]]$negative)), 
                          label = rep(0, length(Partition_2[[i]]$negative)),
                          dataset = aa_neg_frame_set)
  
  aa_df_list2[[i]] <- rbind(aa_df_pos, aa_df_neg)
  
}
Partition_2_chunk <- aa_df_list2
names(Partition_2_chunk) <- names(Partition_2)
Partition_2_chunk_dfs <- do.call(rbind, Partition_2_chunk)
table(Partition_2_chunk$Malat1$dataset)

for(i in 1:length(Partition_2_chunk)){
  Partition_2_chunk[[i]]$tile_name <- as.character(levels(Partition_2_chunk[[i]]$tile_name)[as.numeric(Partition_2_chunk[[i]]$tile_name)])
}

# check if tiles are sorted correctly by their numeric part
aatts <- rownames(territory_assigned_interaction_3Mfilered)[Partition_2[[i]]$negative]

aanu <- as.numeric(unlist(lapply(strsplit(aatts, "_"), "[[", 2)))
identical(aanu, sort(aanu)) # True -> hence I can use the indecis for distance




aachr_vecs3 <- list()
for(i in 1:nrow(aachr_size2)){
  aachr_vecs3[[i]] <-   array(dim =  ceiling(aachr_size2$V2[i]/1000))
  names(aachr_vecs3[[i]]) <- c(1:length(aachr_vecs3[[i]]))
  #aachr_vecs[[i]][aa_tilenamesp2[which(aa_tilenamesp1 %in% aachr_size2$V1[i])]] <- 1
}
names(aachr_vecs3) <- aachr_size2$V1

aachr_vecs4 <- list()
for(i in 1:nrow(aachr_size2)){
  aachr_vecs4[[i]] <-   array(dim =  ceiling(aachr_size2$V2[i]/1000))
  names(aachr_vecs4[[i]]) <- c(1:length(aachr_vecs4[[i]]))
  #aachr_vecs[[i]][aa_tilenamesp2[which(aa_tilenamesp1 %in% aachr_size2$V1[i])]] <- 1
}
names(aachr_vecs4) <- aachr_size2$V1


for(i in 1:length(aachr_vecs3)){
  print(i)
  aaw_lnc <- lncRNA_chosen_gt1k_uniqTiles$gene_name[lncRNA_chosen_gt1k_uniqTiles$chromosome_name  == names(aachr_vecs3)[i]]
  aa_wp <- which(names(Partition_2_chunk) %in% aaw_lnc)
  for(j in 1:length(aa_wp)){
    aa_tilename_pos_train <- Partition_2_chunk[[aa_wp[j]]]$tile_name[(Partition_2_chunk[[aa_wp[j]]]$dataset == "train") & (Partition_2_chunk[[aa_wp[j]]]$label == 1) ]
    aa_tilename_pos_test <- Partition_2_chunk[[aa_wp[j]]]$tile_name[(Partition_2_chunk[[aa_wp[j]]]$dataset == "test") & (Partition_2_chunk[[aa_wp[j]]]$label == 1) ]
    aa_tilename_pos_valid <- Partition_2_chunk[[aa_wp[j]]]$tile_name[(Partition_2_chunk[[aa_wp[j]]]$dataset == "valid") & (Partition_2_chunk[[aa_wp[j]]]$label == 1) ]
    
    
    aa_tilename_pos_train2<- as.numeric(unlist(lapply(strsplit(aa_tilename_pos_train, "_"), "[[", 2)))
    aa_tilename_pos_test2<- as.numeric(unlist(lapply(strsplit(aa_tilename_pos_test, "_"), "[[", 2)))
    aa_tilename_pos_valid2<- as.numeric(unlist(lapply(strsplit(aa_tilename_pos_valid, "_"), "[[", 2)))
    
    aa_tilename_neg_train <- Partition_2_chunk[[aa_wp[j]]]$tile_name[(Partition_2_chunk[[aa_wp[j]]]$dataset == "train") & (Partition_2_chunk[[aa_wp[j]]]$label == 0) ]
    aa_tilename_neg_test <- Partition_2_chunk[[aa_wp[j]]]$tile_name[(Partition_2_chunk[[aa_wp[j]]]$dataset == "test") & (Partition_2_chunk[[aa_wp[j]]]$label == 0) ]
    aa_tilename_neg_valid <- Partition_2_chunk[[aa_wp[j]]]$tile_name[(Partition_2_chunk[[aa_wp[j]]]$dataset == "valid") & (Partition_2_chunk[[aa_wp[j]]]$label == 0) ]
    
    
    aa_tilename_neg_train2<- as.numeric(unlist(lapply(strsplit(aa_tilename_neg_train, "_"), "[[", 2)))
    aa_tilename_neg_test2<- as.numeric(unlist(lapply(strsplit(aa_tilename_neg_test, "_"), "[[", 2)))
    aa_tilename_neg_valid2<- as.numeric(unlist(lapply(strsplit(aa_tilename_neg_valid, "_"), "[[", 2)))
    
    
    aachr_vecs3[[i]][aa_tilename_pos_test2] <- 1.8 + rnorm(n = length(aa_tilename_pos_test2), mean = 0, sd = 0.01)
    aachr_vecs3[[i]][aa_tilename_pos_valid2] <- 1.6 + rnorm(n = length(aa_tilename_pos_valid2), mean = 0, sd = 0.01)
    aachr_vecs3[[i]][aa_tilename_pos_train2] <- 1.4 + rnorm(n = length(aa_tilename_pos_train2), mean = 0, sd = 0.01)
    
    aachr_vecs3[[i]][aa_tilename_neg_test2] <- 0.8 + rnorm(n = length(aa_tilename_neg_test2), mean = 0, sd = 0.01)
    aachr_vecs3[[i]][aa_tilename_neg_valid2] <- 0.6 + rnorm(n = length(aa_tilename_neg_valid2), mean = 0, sd = 0.01)
    aachr_vecs3[[i]][aa_tilename_neg_train2] <- 0.4 + rnorm(n = length(aa_tilename_neg_train2), mean = 0, sd = 0.01)
    
    
    aachr_vecs4[[i]][aa_tilename_pos_train2] <- 5
    aachr_vecs4[[i]][aa_tilename_pos_test2] <- 2
    aachr_vecs4[[i]][aa_tilename_pos_valid2] <- 3
    
    aachr_vecs4[[i]][aa_tilename_neg_train2] <- 5
    aachr_vecs4[[i]][aa_tilename_neg_test2] <- 2
    aachr_vecs4[[i]][aa_tilename_neg_valid2] <- 3
    
  }
}

png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/chr_pos_neg_chosen_patition_2_chunk_partitioning.png",       
    width = 40*300,        # 5 x 300 pixels
    height = 35*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10) 
par(mfrow = c(16,1), mar = c(2,4,3,4))
for(i in 1:length(aachr_vecs)){
  print(i)
  plot(aachr_vecs[[i]], main = paste(names(aachr_vecs)[i], aachr_size2$V2[i]),
       las = 2, xlab = "genomic position", xaxt = "n", ylim = c(0,3.1), pch = 19, cex = 0.3)
  points(aachr_vecs2[[i]], pch = 19, col = 4, cex = 0.3)
  points(aachr_vecs3[[i]], pch = 8, col = aachr_vecs4[[i]], cex = 0.3)
}
dev.off()
save(list = c("Partition_2_dfs","Partition_2_chunk_dfs"), file= "~/Documents/Shayan/BioInf/lncRNA/Partition_2_chunk_random.RData")

################################################################################################################
# partition two feature gathering

Triplex_feaures_partition2
Tile_kmer_features_partition2
ChIPATLAS_features_partition2 <- ChIPATLAS_features[match(Partition_2_dfs$tile_name, rownames(ChIPATLAS_features)),]
RBP_features_partition2 <- RBP_features[match(Partition_2_dfs$tile_name, rownames(RBP_features)),]
ChIPATLAS_features_partition2[is.na(ChIPATLAS_features_partition2)] <- 0
RBP_features_partition2[is.na(RBP_features_partition2)] <- 0


mESC_CHIPATLAS_lncRNA_features
kmer_features_lncRNA28_normalized
mESC_RBP_lncRNA_features

# expression RNAseq tile
mESC_Renlab_RNASeq_partition2  <- numeric(nrow(Partition_2_dfs))
names(mESC_Renlab_RNASeq_partition2) <- Partition_2_dfs$tile_name
aamESC_Renlab_tiled_lisRNAseq <-  mESC_Renlab_tiled_list$`GSM723776_RenLab-RNA-Seq-mESC-ZY28.bed`[mESC_Renlab_tiled_list$`GSM723776_RenLab-RNA-Seq-mESC-ZY28.bed`$tile %in% Partition_2_dfs$tile_name,]
mESC_Renlab_RNASeq_partition2[match(aamESC_Renlab_tiled_lisRNAseq$tile, Partition_2_dfs$tile_name)] <- aamESC_Renlab_tiled_lisRNAseq$`aadf_agg$score`

# expression cage tile
aamESC_CAGE_mm9_tiled <- mESC_CAGE_mm9_tiled[mESC_CAGE_mm9_tiled$tile %in% Partition_2_dfs$tile_name,]
mESC_CAGE_mm9_tiled_partition2 <- numeric(nrow(Partition_2_dfs))
mESC_CAGE_mm9_tiled_partition2[match(aamESC_CAGE_mm9_tiled$tile, Partition_2_dfs$tile_name)] <- aamESC_CAGE_mm9_tiled$`aadf_agg$CAGE_score`
names(mESC_CAGE_mm9_tiled_partition2) <- Partition_2_dfs$tile_name

save(list = c("mESC_CAGE_mm9_tiled",
              "mESC_Renlab_tiled_list"),
     file = "~/Documents/Shayan/BioInf/lncRNA/CAGE_RNAseq_tiled.RData")
# expression RNAseq lncRNA
mESC_RNAseq_lncRNA_features
mESC_RNAseq_lncRNA_features_full <- numeric(length = nrow(lncRNA_chosen_gt1k_uniqTiles))
names(mESC_RNAseq_lncRNA_features_full) <- lncRNA_chosen_gt1k_uniqTiles$gene_name
mESC_RNAseq_lncRNA_features_full[match(names(mESC_RNAseq_lncRNA_features), names(mESC_RNAseq_lncRNA_features_full))] <- mESC_RNAseq_lncRNA_features



# expression cage lncRNA
mESC_CAGE_lncRNA_features
mESC_CAGE_lncRNA_features_full <- numeric(length = nrow(lncRNA_chosen_gt1k_uniqTiles))
names(mESC_CAGE_lncRNA_features_full) <- lncRNA_chosen_gt1k_uniqTiles$gene_name
mESC_CAGE_lncRNA_features_full[match(names(mESC_CAGE_lncRNA_features), names(mESC_CAGE_lncRNA_features_full))] <- mESC_CAGE_lncRNA_features

save(list = c("mESC_CAGE_lncRNA_features_full",
              "mESC_RNAseq_lncRNA_features_full", 
              "mESC_RBP_lncRNA_features", 
              "kmer_features_lncRNA28_normalized",
              "mESC_CHIPATLAS_lncRNA_features"),
     file = "~/Documents/Shayan/BioInf/lncRNA/lncRNA_feaures.RData")


Partition_2_dfs$tile_name <- levels(Partition_2_dfs$tile_name)[as.numeric(Partition_2_dfs$tile_name)]
Partition_2_chunk_dfs$tile_name <- levels(Partition_2_chunk_dfs$tile_name)[as.numeric(Partition_2_chunk_dfs$tile_name)]
# 
Partition_2_feature_mat_tile <- cbind(ChIPATLAS_features_partition2, # remember to add repeats
                                      Tile_kmer_features_partition2,
                                      RBP_features_partition2,
                                      mESC_Renlab_RNASeq_partition2,
                                      mESC_CAGE_mm9_tiled_partition2)
colnames(Partition_2_feature_mat_tile)[ncol(Partition_2_feature_mat_tile) - 1] <- "RNAseq"
colnames(Partition_2_feature_mat_tile)[ncol(Partition_2_feature_mat_tile)] <- "CAGE"
  
Partition_2_feature_mat_owner <- cbind(mESC_CHIPATLAS_lncRNA_features[match(Partition_2_dfs$owner, rownames(mESC_CHIPATLAS_lncRNA_features)),],
                                       kmer_features_lncRNA28_normalized[match(Partition_2_dfs$owner, rownames(kmer_features_lncRNA28_normalized)),],
                                       mESC_RBP_lncRNA_features[match(Partition_2_dfs$owner, rownames(mESC_RBP_lncRNA_features)),],
                                       mESC_RNAseq_lncRNA_features_full[match(Partition_2_dfs$owner, names(mESC_RNAseq_lncRNA_features_full))],
                                       mESC_CAGE_lncRNA_features_full[match(Partition_2_dfs$owner, names(mESC_CAGE_lncRNA_features_full))])
colnames(Partition_2_feature_mat_owner)[ncol(Partition_2_feature_mat_owner) - 1] <- "RNAseq"
colnames(Partition_2_feature_mat_owner)[ncol(Partition_2_feature_mat_owner)] <- "CAGE"
  
  
colnames(Partition_2_feature_mat_owner) <- paste(colnames(Partition_2_feature_mat_owner), "lncRNA", sep = "__")
Partition_2_feature_mat_pair <- cbind(distance_feature_partition_2,
                                      Triplex_feaures_partition2[match(Partition_2_dfs$tile_name , Triplex_feaures_partition2$tile_name), 3])
colnames(Partition_2_feature_mat_pair) <- c("distance", "triplex")

# how to remove - from names
aatst <- colnames(ChIPATLAS_features_partition2)
aatst2 <- gsub(pattern = "-", replacement = "", x = aatst)



#########################################################################################################
#########################################################################################################
# looking at distance from lncRNA feature for partition one and two, for each of traning, valid, and test points in each setting (rand, chunk)


(distance_feature_partition_1)
length(distance_feature_partition_2)

aadata_p1 <- cbind(Partition_1_dfs,Partition_1_chunk_dfs$dataset, distance_feature_partition_1)
aadata_p2 <- cbind(Partition_2_dfs,Partition_2_chunk_dfs$dataset, distance_feature_partition_2)
colnames(aadata_p1) <- c("tile_name", "owner",     "label" ,  "random_partition" , "chunk_partition", "distance" )
colnames(aadata_p2) <- c("tile_name", "owner",     "label" ,  "random_partition" , "chunk_partition", "distance" )
aadata_p1$label[aadata_p1$label == 1] <- "pos"
aadata_p1$label[aadata_p1$label == 0] <- "neg"
aadata_p2$label[aadata_p2$label == 1] <- "pos"
aadata_p2$label[aadata_p2$label == 0] <- "neg"
aadata_p1$label <- factor(aadata_p1$label, levels = c("neg", "pos"))
aadata_p2$label <- factor(aadata_p2$label, levels = c("neg", "pos"))
# aatttg_1_rand <- summarySE(aadata_p1, measurevar=c("distance"), groupvars=c("label", "random_partition"), na.rm=T)
# aatttg_1_chun <- summarySE(aadata_p1, measurevar=c("distance"), groupvars=c("label", "chunk_partition"), na.rm=T)
# aatttg_2_rand <- summarySE(aadata_p2, measurevar=c("distance"), groupvars=c("label", "random_partition"), na.rm=T)
# aatttg_2_chun <- summarySE(aadata_p2, measurevar=c("distance"), groupvars=c("label", "chunk_partition"), na.rm=T)
ggplot(aadata_p1, aes(x=random_partition, y=distance, fill=label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot() +
  ggtitle("partition_1_random distance to owner origin")+
  ylab("Dist (bp)")+
  xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank())

ggplot(aadata_p1, aes(x=chunk_partition, y=distance, fill=label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot() +
  ggtitle("partition_1_chunk distance to owner origin")+
  ylab("Dist (bp)")+
  xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank()) 

ggplot(aadata_p1, aes(x=owner, y=distance, fill=label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot() +
  geom_vline(xintercept = seq(1.5,28.5,1), color ="blue", size = 1.5) +
  ggtitle("partition_1 distance to owner origin by lncRNA")+
  ylab("Dist (bp)")+
  xlab("lncRNA") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90)) 

ggplot(aadata_p1, aes(x=label, y=distance+1)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot() +
  ggtitle("partition 1 distance to owner origin")+
  ylab("Dist (bp)")+
  xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank())

ggplot(aadata_p2, aes(x=random_partition, y=distance, fill=label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot() +
  ggtitle("partition_2_random distance to owner origin")+
  ylab("Dist (bp)")+
  xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank())

ggplot(aadata_p2, aes(x=chunk_partition, y=distance, fill=label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot() +
  ggtitle("partition_2_chunk distance to owner origin")+
  ylab("Dist (bp)")+
  xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank()) 

ggplot(aadata_p2, aes(x=owner, y=distance, fill=label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot() +
  geom_vline(xintercept = seq(1.5,28.5,1), color ="blue", size = 1.5) +
  ggtitle("partition_2 distance to owner origin by lncRNA")+
  ylab("Dist (bp)")+
  xlab("lncRNA") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90)) 

ggplot(aadata_p2, aes(x=label, y=distance+1)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot() +
  ggtitle("partition 2 distance to owner origin")+
  ylab("Dist (bp)")+
  xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank())

ggplot(aadata_p2, aes(x=distance)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_histogram() 
  # geom_vline(xintercept = seq(1.5,28.5,1), color ="blue", size = 1.5) +
  # ggtitle("partition_2 distance to owner origin by lncRNA")+
  # ylab("Dist (bp)")+
  # xlab("lncRNA") + 
  # theme(plot.title = element_text(hjust = 0.5),
  #       panel.grid.major.x = element_blank(),
  #       axis.text.x = element_text(angle = 90)) 


# compute the same plots as above for the full dataset to see the original distribution of data points
aa <- hist(log10(distance_feature_partition_1[Partition_1_dfs$label == 0]))

#MM9_1kb_tiled_owner_3Mfilter_noNA <- MM9_1kb_tiled_owner_3Mfilter[!is.na(MM9_1kb_tiled_owner_3Mfilter)]
aa_distance <- numeric(length = length(MM9_1kb_tiled_owner_3Mfilter))
aalnc <- unique(MM9_1kb_tiled_owner_3Mfilter)
aalnc <- aalnc[!is.na(aalnc)]
for(i in 1:length(aalnc)){
  print(i)
  aa_cur_qu <- MM9_1kb_tiled_GR_filtered[match(names(MM9_1kb_tiled_owner_3Mfilter)[MM9_1kb_tiled_owner_3Mfilter %in% aalnc[i]],
                                               MM9_1kb_tiled_GR_filtered$tile)]
  aa_cur_su <- lncRNA_chosen_gt1k_uniqTiles_GR[lncRNA_chosen_gt1k_uniqTiles_GR$gene_name == aalnc[i]]
  aadist <- GenomicRanges::distance(x = aa_cur_qu, y = aa_cur_su, select = "all")
  aa_distance[match(aa_cur_qu$tile, names(MM9_1kb_tiled_owner_3Mfilter))] <- aadist
}
aa_distance[is.na(aa_distance)] <- max(aa_distance, na.rm = T) + 1000000

distance_feature_allTiles <- aa_distance
names(distance_feature_allTiles) <- names(MM9_1kb_tiled_owner_3Mfilter)
distance_feature_allTiles_nozero <- distance_feature_allTiles
distance_feature_allTiles_nozero[distance_feature_allTiles == 0] <- NA
distance_feature_allTiles_nozero_cis <- distance_feature_allTiles_nozero
distance_feature_allTiles_nozero_cis[distance_feature_allTiles_nozero_cis == max(distance_feature_allTiles_nozero_cis)] <- NA

aadata_p3 <- data.frame(tile_name = names(distance_feature_allTiles),
                        owner =MM9_1kb_tiled_owner_3Mfilter,
                        label = MM9_1kb_tiled_owner_labels_binary_3Mfilter  , 
                        distance = distance_feature_allTiles_nozero_cis)
aadata_p3$label[aadata_p3$label == 1] <- "pos"
aadata_p3$label[aadata_p3$label == 0] <- "neg"
aadata_p3$label <- factor(aadata_p3$label, levels = c("neg", "pos"))

ggplot(aadata_p3, aes(x=owner, y=distance+1, fill=label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot() +
  geom_vline(xintercept = seq(1.5,28.5,1), color ="blue", size = 1.5) +
  ggtitle("All points distance to owner origin by lncRNA")+
  ylab("Dist (bp)")+
  xlab("lncRNA") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90)) 

ggplot(aadata_p3, aes(x=label, y=distance+1)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  #ggtitle("all points (2116029) distance to owner origin")+
  ggtitle("distance to lncRNA origin\nfor all datapoints")+
  ylab("Dist (bp)")+
  xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank())

aahist <- hist(distance_feature_allTiles, breaks = 100, freq = F)

qqplot(x = log10(distance_feature_allTiles_nozero_cis[MM9_1kb_tiled_owner_labels_binary_3Mfilter == 0] + 1),
       y = log10(distance_feature_allTiles_nozero_cis[MM9_1kb_tiled_owner_labels_binary_3Mfilter == 1] + 1),
       xlab = "Negative", ylab = "Positive", 
       xlim = range(log10(distance_feature_allTiles_nozero_cis+1), na.rm = T),
       ylim = range(log10(distance_feature_allTiles_nozero_cis+1), na.rm = T), 
       col = 3, pch = 19, cex = 0.6)
abline(a = 0, b = 1, col=2)
par(new= T)
qqplot(x = log10(distance_feature_partition_6_cis[Partition_6_dfs$label == 0] + 1),
       y = log10(distance_feature_partition_6_cis[Partition_6_dfs$label == 1] + 1),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       #xlim = range(distance_feature_partition_6_cis, na.rm = T),
       #ylim = range(distance_feature_partition_6_cis, na.rm = T),
       col = 4, pch = 19, cex = 0.6 )

png(filename = "~/Desktop/lncRNA_paper/Figs/qqplot_distance.png",    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 11) 
qqplot(x = (distance_feature_allTiles_nozero_cis[MM9_1kb_tiled_owner_labels_binary_3Mfilter == 0]),
       y = (distance_feature_allTiles_nozero_cis[MM9_1kb_tiled_owner_labels_binary_3Mfilter == 1] ),
       xlab = "Negative", ylab = "Positive", 
       xlim = range((distance_feature_allTiles_nozero_cis), na.rm = T),
       ylim = range((distance_feature_allTiles_nozero_cis), na.rm = T), 
       col = 3, pch = 19, cex = 0.6, main = "Q-Q plot of distance to lncRNA origin")
abline(a = 0, b = 1, col=2)
par(new= T)
qqplot(x = (distance_feature_partition_6_cis[Partition_6_dfs$label == 0] ),
       y = (distance_feature_partition_6_cis[Partition_6_dfs$label == 1] ),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       #xlim = range(distance_feature_partition_6_cis, na.rm = T),
       #ylim = range(distance_feature_partition_6_cis, na.rm = T),
       col = 4, pch = 19, cex = 0.6 )
legend("topleft", legend = c("All data-points", "selected data-points"), fill = c(3,4))
dev.off()
#########################################################################################################
# making a third partition by controlling for distance from lncRNA to negatives to be similar to lncRNA distance to positives.

Partition_3 <- list()
aa_ovl <- findOverlaps(query = MM9_1kb_tiled_GR_filtered, subject = lncRNA_chosen_gt1k_uniqTiles_GR)
aa_thresh <- 3000
neg_to_pos_ratio <- 10
aaallhist <- hist(distance_feature_allTiles, breaks = 1000, freq = F)
set.seed(seed = 86786)
for(i in 1:ncol(territory_assigned_interaction_label_3Mfiltered)){
  aa_w_pos <- which(territory_assigned_interaction_label_3Mfiltered[, i] == 1)
  aa_w_pos <- setdiff(aa_w_pos,unique(aa_ovl@from))
  aa_w_Terr <- which(territory_assigned_interaction_3Mfilered[, i] == 1)
  aa_w_neg <- setdiff(aa_w_Terr,aa_w_pos)
  aa_w_neg <- setdiff(aa_w_neg,unique(aa_ovl@from))
  aa_w_neg <- aa_w_neg[rownames(territory_assigned_interaction_label_3Mfiltered)[aa_w_neg] %in% names(MM9_1kb_tiled_owner_3Mfilter_noNA)]
  if(length(aa_w_pos) > aa_thresh){
    aamyneg_dist_hist <- hist(distance_feature_allTiles[match(rownames(territory_assigned_interaction_label_3Mfiltered)[aa_w_neg],
                                                         names(distance_feature_allTiles))],
                              breaks = aaallhist$breaks, probability = T)

    aamypos_dist <- distance_feature_allTiles[match(rownames(territory_assigned_interaction_label_3Mfiltered)[aa_w_pos],
                                                    names(distance_feature_allTiles))]
    aa_prob_sum <- 0
    aa_myk <- 10
    while(aa_prob_sum < aa_thresh){
      aa_density_smooth <- zoo::rollmean(aamyneg_dist_hist$counts, k = aa_myk, fill = "extend")
      aa_myprob_smooth <- aa_density_smooth[findInterval(x = aamypos_dist, vec =  aaallhist$breaks)]
      aa_prob_sum <- sum(aa_myprob_smooth > 0)
      aa_myk <- aa_myk + 10
      stopifnot(aa_myk < (length(aa_myprob)/2))
    }
    aa_myprob_smooth[aa_myprob_smooth <0] <- 0

    aa_w_pos_sampled <- sort(sample(x = aa_w_pos, size = aa_thresh, replace = F, prob = aa_myprob_smooth))
  }else{
    aa_w_pos_sampled <- sort(aa_w_pos)
  }
  if(length(aa_w_neg) > (neg_to_pos_ratio * length(aa_w_pos_sampled))){
    aamypos_dist_hist <- hist(distance_feature_allTiles[match(rownames(territory_assigned_interaction_label_3Mfiltered)[aa_w_pos_sampled],
                                                              names(distance_feature_allTiles))], breaks = aaallhist$breaks, probability = T)
    
    aamyneg_dist <- distance_feature_allTiles[match(rownames(territory_assigned_interaction_label_3Mfiltered)[aa_w_neg],
                                                    names(distance_feature_allTiles))]
    aa_prob_sum <- 0
    aa_myk <- 10
    while(aa_prob_sum < (neg_to_pos_ratio * length(aa_w_pos_sampled))){
      aa_density_smooth <- zoo::rollmean(aamypos_dist_hist$counts, k = aa_myk, fill = "extend")
      aa_myprob_smooth <- aa_density_smooth[findInterval(x = aamyneg_dist, vec =  aaallhist$breaks)]
      aa_prob_sum <- sum(aa_myprob_smooth > 0)
      aa_myk <- aa_myk + 10
      stopifnot(aa_myk < (length(aa_myprob)/2))
    }
    aa_myprob_smooth[aa_myprob_smooth <0] <- 0
    aa_w_neg_sampled <- sort(sample(x = aa_w_neg,
                                    size = (neg_to_pos_ratio * length(aa_w_pos_sampled)),
                                    replace = F,
                                    prob = aa_myprob_smooth))
  }else{
    aa_w_neg_sampled <- aa_w_neg
  }
  Partition_3[[i]] <- list(positive=aa_w_pos_sampled, 
                           negative=aa_w_neg_sampled,
                           extra_positive = setdiff(aa_w_pos,aa_w_pos_sampled), 
                           extra_negative = setdiff(aa_w_neg,aa_w_neg_sampled))
  
}
names(Partition_3) <- colnames(territory_assigned_interaction_label_3Mfiltered)

aa_ps <- lapply(Partition_3, "[[", 1)
aa_ng <- lapply(Partition_3, "[[", 2)
aa_data_cnt <- rbind(unlist(lapply(aa_ps, length)),unlist(lapply(aa_ng, length)))
barplot(aa_data_cnt, las = 2, legend.text = c("+", "-")   ,
        main = "data points per lncRNA",
        args.legend=list(x = "topright", bty = "n", 
                         inset=c(0, 0),
                         cex = 1, 
                         y.intersp = 0.6,
                         x.intersp = 0.5,
                         text.width = 1))

aa_ps2 <- lapply(Partition_3, "[[", 3)
aa_ng2 <- lapply(Partition_3, "[[", 4)
aa_data_cnt2 <- rbind(unlist(lapply(aa_ps2, length)),unlist(lapply(aa_ng2, length)))
barplot(aa_data_cnt2, las = 2, legend.text = c("+", "-")   ,
        main = "extra data points per lncRNA",
        args.legend=list(x = "topright", bty = "n", 
                         inset=c(0, 0),
                         cex = 1, 
                         y.intersp = 0.6,
                         x.intersp = 0.5,
                         text.width = 1))


aachr_vecs2 <- list()
for(i in 1:nrow(aachr_size2)){
  aachr_vecs2[[i]] <-   array(dim =  ceiling(aachr_size2$V2[i]/1000))
  names(aachr_vecs2[[i]]) <- c(1:length(aachr_vecs2[[i]]))
  #aachr_vecs[[i]][aa_tilenamesp2[which(aa_tilenamesp1 %in% aachr_size2$V1[i])]] <- 1
}
names(aachr_vecs2) <- aachr_size2$V1

for(i in 1:length(aachr_vecs2)){
  aaw_lnc <- lncRNA_chosen_gt1k_uniqTiles$gene_name[lncRNA_chosen_gt1k_uniqTiles$chromosome_name  == names(aachr_vecs2)[i]]
  aa_wp <- which(names(Partition_3) %in% aaw_lnc)
  for(j in 1:length(aa_wp)){
    aa_tilename_pos <- rownames(territory_assigned_interaction_3Mfilered)[Partition_3[[aa_wp[j]]]$positive]
    aa_tilename_pos2<- as.numeric(unlist(lapply(strsplit(aa_tilename_pos, "_"), "[[", 2)))
    aa_tilename_neg <- rownames(territory_assigned_interaction_3Mfilered)[Partition_3[[aa_wp[j]]]$negative]
    aa_tilename_neg2<- as.numeric(unlist(lapply(strsplit(aa_tilename_neg, "_"), "[[", 2)))
    
    aachr_vecs2[[i]][aa_tilename_pos2] <- 1.25 + rnorm(n = length(aa_tilename_pos2), mean = 0, sd = 0.01)
    aachr_vecs2[[i]][aa_tilename_neg2] <- 0.25 + rnorm(n = length(aa_tilename_neg2), mean = 0, sd = 0.01)
    
    
  }
}

png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/chr_pos_neg_chosen_partition3.png",       
    width = 40*300,        # 5 x 300 pixels
    height = 30*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10) 
par(mfrow = c(16,1), mar = c(3,4,3,4))
for(i in 1:length(aachr_vecs)){
  print(i)
  plot(aachr_vecs[[i]], main = paste(names(aachr_vecs)[i], aachr_size2$V2[i]),
       las = 2, xlab = "genomic position", xaxt = "n", ylim = c(0,3.1), pch = 19, cex = 0.2)
  points(aachr_vecs2[[i]], pch = 8, col = 2, cex = 0.5)
}
dev.off()


aa_df_list <- list()
for(i in 1:length(Partition_3)){
  aa_df_pos <- data.frame(tile_name = rownames(territory_assigned_interaction_3Mfilered)[Partition_3[[i]]$positive],
                          owner = rep(names(Partition_3)[i], length(Partition_3[[i]]$positive)), 
                          label = rep(1, length(Partition_3[[i]]$positive)), 
                          dataset = sample(c("train", "valid", "test"), 
                                           size = length(Partition_3[[i]]$positive),
                                           replace = T, 
                                           prob = c(0.7, 0.2, 0.1)))
  aa_df_neg <- data.frame(tile_name = rownames(territory_assigned_interaction_3Mfilered)[Partition_3[[i]]$negative],
                          owner = rep(names(Partition_3)[i], length(Partition_3[[i]]$negative)), 
                          label = rep(0, length(Partition_3[[i]]$negative)),
                          dataset = sample(c("train", "valid", "test"), 
                                           size = length(Partition_3[[i]]$negative),
                                           replace = T, 
                                           prob = c(0.7, 0.2, 0.1))
  )
  
  aa_df_list[[i]] <- rbind(aa_df_pos, aa_df_neg)
  
}

names(aa_df_list) <- names(Partition_3)

Partition_3_dfs <- do.call(rbind, aa_df_list)
Partition_3_df_list <- aa_df_list
for(i in 1:length(Partition_3_df_list)){
  Partition_3_df_list[[i]]$tile_name <- as.character(levels(Partition_3_df_list[[i]]$tile_name)[as.numeric(Partition_3_df_list[[i]]$tile_name)])
}
Partition_3_dfs$tile_name <- as.character(levels(Partition_3_dfs$tile_name)[as.numeric(Partition_3_dfs$tile_name)])

aadata_p4 <- cbind(Partition_3_dfs, 
                   distance_feature_allTiles[match(Partition_3_dfs$tile_name, names(distance_feature_allTiles))])

colnames(aadata_p4) <- c("tile_name", "owner",     "label" ,  "dataset"  ,"distance" )
aadata_p4$label[aadata_p4$label == 1] <- "pos"
aadata_p4$label[aadata_p4$label == 0] <- "neg"
aadata_p4$label <- factor(aadata_p4$label, levels = c("neg", "pos"))



ggplot(aadata_p4, aes(x=owner, y=distance, fill=label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot() +
  geom_vline(xintercept = seq(1.5,28.5,1), color ="blue", size = 1.5) +
  ggtitle("Partition3 points distance to owner origin by lncRNA")+
  ylab("Dist (bp)")+
  xlab("lncRNA") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90)) 

ggplot(aadata_p4, aes(x=dataset, y=distance, fill=label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot() +
  ggtitle("partition_3_random distance to owner origin")+
  ylab("Dist (bp)")+
  xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank())

ggplot(aadata_p4, aes(x=label, y=distance+1)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot() +
  ggtitle("partition 3 distance to owner origin")+
  ylab("Dist (bp)")+
  xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank())

#########################################################################################################
# Partition_4 --> use the same strategy as partition 3 but choose cap 2000 pos and 4 times that negatives

Partition_4 <- list()
aa_ovl <- findOverlaps(query = MM9_1kb_tiled_GR_filtered, subject = lncRNA_chosen_gt1k_uniqTiles_GR)
aa_thresh <- 2000
neg_to_pos_ratio <- 4
nu_closest_neg_thresh <- 100 # number of negatives with the same closest positive
aaallhist <- hist(distance_feature_allTiles, breaks = 5000, freq = F)
set.seed(seed = 57668)
for(i in 1:ncol(territory_assigned_interaction_label_3Mfiltered)){
  print(i)
  aa_w_pos <- which(territory_assigned_interaction_label_3Mfiltered[, i] == 1)
  aa_w_pos <- setdiff(aa_w_pos,unique(aa_ovl@from))
  aa_w_Terr <- which(territory_assigned_interaction_3Mfilered[, i] == 1)
  aa_w_neg <- setdiff(aa_w_Terr,aa_w_pos)
  aa_w_neg <- setdiff(aa_w_neg,unique(aa_ovl@from))
  aa_w_neg <- aa_w_neg[rownames(territory_assigned_interaction_label_3Mfiltered)[aa_w_neg] %in% names(MM9_1kb_tiled_owner_3Mfilter_noNA)]
  
  # # compute the closest positive for each negative of this lncRNA
  # # assign each negative to its closest positive
  # aa_neg_pos_df <- data.frame(neg_index = aa_w_neg, 
  #                             closest_pos = numeric(length = length(aa_w_neg)))
  # for(aamuneg in 1:length(aa_w_neg)){
  #   aa_neg_pos_df$closest_pos[aamuneg] <- aa_w_pos[which.min(abs(aa_w_neg[aamuneg] - aa_w_pos))]
  # }
  # # allow each positive to have maximum of 100 negative points who are closest to it
  # aatab <- table(aa_neg_pos_df$closest_pos)
  # aa_w_pos_ABOVElimit <- as.numeric(names(aatab)[which(aatab > nu_closest_neg_thresh)])
  # aa_toberemoved <- numeric(0)
  # for(aacur_lim in 1:length(aa_w_pos_ABOVElimit)){
  #   aa_curlimneg <- aa_neg_pos_df$neg_index[aa_neg_pos_df$closest_pos %in% aa_w_pos_ABOVElimit[aacur_lim]]
  #   aarem <- aa_curlimneg[sort(abs(aa_curlimneg - aa_w_pos_ABOVElimit[aacur_lim]), index.return = T)$ix[(nu_closest_neg_thresh+1):length(aa_curlimneg)]]
  #   aa_toberemoved <- c(aa_toberemoved, aarem)
  # }
  # print("length(aa_w_neg)")
  # print(length(aa_w_neg))
  # print("length(aa_toberemoved)")
  # print(length(aa_toberemoved))
  # aa_w_neg <- setdiff(aa_w_neg,aa_toberemoved)

  
  
  if(length(aa_w_pos) > aa_thresh){
    aamyneg_dist_hist <- hist(distance_feature_allTiles[match(rownames(territory_assigned_interaction_label_3Mfiltered)[aa_w_neg],
                                                              names(distance_feature_allTiles))],
                              breaks = aaallhist$breaks, probability = T)
    
    aamypos_dist <- distance_feature_allTiles[match(rownames(territory_assigned_interaction_label_3Mfiltered)[aa_w_pos],
                                                    names(distance_feature_allTiles))]
    aa_prob_sum <- 0
    aa_myk <- 10
    while(aa_prob_sum < aa_thresh){
      print("pos sampling")
      print(aa_myk)
      aa_density_smooth <- zoo::rollmean(aamyneg_dist_hist$counts, k = aa_myk, fill = "extend")
      aa_myprob_smooth <- aa_density_smooth[findInterval(x = aamypos_dist, vec =  aaallhist$breaks)]
      aa_prob_sum <- sum(aa_myprob_smooth > 0)
      aa_myk <- aa_myk + 10
      stopifnot(aa_myk < (length(aa_myprob)/2))
    }
    aa_myprob_smooth[aa_myprob_smooth <0] <- 0
    
    aa_w_pos_sampled <- sort(sample(x = aa_w_pos, size = aa_thresh, replace = F, prob = aa_myprob_smooth))
  }else{
    aa_w_pos_sampled <- sort(aa_w_pos)
  }
  if(length(aa_w_neg) > (neg_to_pos_ratio * length(aa_w_pos_sampled))){
    aamypos_dist_hist <- hist(distance_feature_allTiles[match(rownames(territory_assigned_interaction_label_3Mfiltered)[aa_w_pos_sampled],
                                                              names(distance_feature_allTiles))], breaks = aaallhist$breaks, probability = T)
    
    aamyneg_dist <- distance_feature_allTiles[match(rownames(territory_assigned_interaction_label_3Mfiltered)[aa_w_neg],
                                                    names(distance_feature_allTiles))]
    aa_prob_sum <- 0
    aa_myk <- 10
    while(aa_prob_sum < (neg_to_pos_ratio * length(aa_w_pos_sampled))){
      print("neg sampling")
      print(aa_myk)
      aa_density_smooth <- zoo::rollmean(aamypos_dist_hist$counts, k = aa_myk, fill = "extend")
      aa_myprob_smooth <- aa_density_smooth[findInterval(x = aamyneg_dist, vec =  aaallhist$breaks)]
      aa_prob_sum <- sum(aa_myprob_smooth > 0)
      aa_myk <- aa_myk + 10
      stopifnot(aa_myk < (length(aa_myprob)/2))
    }
    aa_myprob_smooth[aa_myprob_smooth <0] <- 0
    aa_w_neg_sampled <- sort(sample(x = aa_w_neg,
                                    size = (neg_to_pos_ratio * length(aa_w_pos_sampled)),
                                    replace = F,
                                    prob = aa_myprob_smooth))
  }else{
    aa_w_neg_sampled <- aa_w_neg
  }
  Partition_4[[i]] <- list(positive=aa_w_pos_sampled, 
                           negative=aa_w_neg_sampled,
                           extra_positive = setdiff(aa_w_pos,aa_w_pos_sampled), 
                           extra_negative = setdiff(aa_w_neg,aa_w_neg_sampled))
  
}

names(Partition_4) <- colnames(territory_assigned_interaction_label_3Mfiltered)

aa_ps <- lapply(Partition_4, "[[", 1)
aa_ng <- lapply(Partition_4, "[[", 2)
aa_data_cnt <- rbind(unlist(lapply(aa_ps, length)),unlist(lapply(aa_ng, length)))
barplot(aa_data_cnt, las = 2, legend.text = c("+", "-")   ,
        main = "data points per lncRNA",
        args.legend=list(x = "topright", bty = "n", 
                         inset=c(0, 0),
                         cex = 1, 
                         y.intersp = 0.6,
                         x.intersp = 0.5,
                         text.width = 1))

aa_ps2 <- lapply(Partition_4, "[[", 3)
aa_ng2 <- lapply(Partition_4, "[[", 4)
aa_data_cnt2 <- rbind(unlist(lapply(aa_ps2, length)),unlist(lapply(aa_ng2, length)))
barplot(aa_data_cnt2, las = 2, legend.text = c("+", "-")   ,
        main = "extra data points per lncRNA",
        args.legend=list(x = "topright", bty = "n", 
                         inset=c(0, 0),
                         cex = 1, 
                         y.intersp = 0.6,
                         x.intersp = 0.5,
                         text.width = 1))


aachr_vecs2 <- list()
for(i in 1:nrow(aachr_size2)){
  aachr_vecs2[[i]] <-   array(dim =  ceiling(aachr_size2$V2[i]/1000))
  names(aachr_vecs2[[i]]) <- c(1:length(aachr_vecs2[[i]]))
  #aachr_vecs[[i]][aa_tilenamesp2[which(aa_tilenamesp1 %in% aachr_size2$V1[i])]] <- 1
}
names(aachr_vecs2) <- aachr_size2$V1

for(i in 1:length(aachr_vecs2)){
  aaw_lnc <- lncRNA_chosen_gt1k_uniqTiles$gene_name[lncRNA_chosen_gt1k_uniqTiles$chromosome_name  == names(aachr_vecs2)[i]]
  aa_wp <- which(names(Partition_4) %in% aaw_lnc)
  for(j in 1:length(aa_wp)){
    aa_tilename_pos <- rownames(territory_assigned_interaction_3Mfilered)[Partition_4[[aa_wp[j]]]$positive]
    aa_tilename_pos2<- as.numeric(unlist(lapply(strsplit(aa_tilename_pos, "_"), "[[", 2)))
    aa_tilename_neg <- rownames(territory_assigned_interaction_3Mfilered)[Partition_4[[aa_wp[j]]]$negative]
    aa_tilename_neg2<- as.numeric(unlist(lapply(strsplit(aa_tilename_neg, "_"), "[[", 2)))
    
    aachr_vecs2[[i]][aa_tilename_pos2] <- 1.25 + rnorm(n = length(aa_tilename_pos2), mean = 0, sd = 0.01)
    aachr_vecs2[[i]][aa_tilename_neg2] <- 0.25 + rnorm(n = length(aa_tilename_neg2), mean = 0, sd = 0.01)
    
    
  }
}

png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/chr_pos_neg_chosen_partition4.png",       
    width = 40*300,        # 5 x 300 pixels
    height = 30*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10) 
par(mfrow = c(16,1), mar = c(3,4,3,4))
for(i in 1:length(aachr_vecs)){
  print(i)
  plot(aachr_vecs[[i]], main = paste(names(aachr_vecs)[i], aachr_size2$V2[i]),
       las = 2, xlab = "genomic position", xaxt = "n", ylim = c(0,3.1), pch = 19, cex = 0.2)
  points(aachr_vecs2[[i]], pch = 8, col = 2, cex = 0.5)
}
dev.off()


aa_df_list <- list()
for(i in 1:length(Partition_4)){
  aa_df_pos <- data.frame(tile_name = rownames(territory_assigned_interaction_3Mfilered)[Partition_4[[i]]$positive],
                          owner = rep(names(Partition_4)[i], length(Partition_4[[i]]$positive)), 
                          label = rep(1, length(Partition_4[[i]]$positive)), 
                          dataset = sample(c("train", "valid", "test"), 
                                           size = length(Partition_4[[i]]$positive),
                                           replace = T, 
                                           prob = c(0.7, 0.2, 0.1)))
  aa_df_neg <- data.frame(tile_name = rownames(territory_assigned_interaction_3Mfilered)[Partition_4[[i]]$negative],
                          owner = rep(names(Partition_4)[i], length(Partition_4[[i]]$negative)), 
                          label = rep(0, length(Partition_4[[i]]$negative)),
                          dataset = sample(c("train", "valid", "test"), 
                                           size = length(Partition_4[[i]]$negative),
                                           replace = T, 
                                           prob = c(0.7, 0.2, 0.1))
  )
  
  aa_df_list[[i]] <- rbind(aa_df_pos, aa_df_neg)
  
}

names(aa_df_list) <- names(Partition_4)

Partition_4_dfs <- do.call(rbind, aa_df_list)
Partition_4_df_list <- aa_df_list
for(i in 1:length(Partition_4_df_list)){
  Partition_4_df_list[[i]]$tile_name <- as.character(levels(Partition_4_df_list[[i]]$tile_name)[as.numeric(Partition_4_df_list[[i]]$tile_name)])
}
Partition_4_dfs$tile_name <- as.character(levels(Partition_4_dfs$tile_name)[as.numeric(Partition_4_dfs$tile_name)])

aadata_p5 <- cbind(Partition_4_dfs,
                   Partition_4_chunk_dfs$dataset,
                   distance_feature_partition_4)

colnames(aadata_p5) <- c("tile_name", "owner",     "label" ,  "dataset_random",  "dataset_chunk" ,"distance" )
aadata_p5$label[aadata_p5$label == 1] <- "pos"
aadata_p5$label[aadata_p5$label == 0] <- "neg"
aadata_p5$label <- factor(aadata_p5$label, levels = c("neg", "pos"))



ggplot(aadata_p5, aes(x=owner, y=distance, fill=label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_violin() +
  geom_vline(xintercept = seq(1.5,28.5,1), color ="blue", size = 1.5) +
  ggtitle("Partition4 points distance to owner origin by lncRNA")+
  ylab("Dist (bp)")+
  xlab("lncRNA") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90)) 

ggplot(aadata_p5, aes(x=dataset_random, y=distance, fill=label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)
              #,scale = "count"
              ) +
  ggtitle("partition_4_random distance to owner origin")+
  ylab("Dist (bp)")+
  xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank())

ggplot(aadata_p5, aes(x=dataset_chunk, y=distance, fill=label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)
              #,scale = "count"
              ) +
  ggtitle("partition_4_chunk distance to owner origin")+
  ylab("Dist (bp)")+
  xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank())

ggplot(aadata_p5, aes(x=label, y=distance+1)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot() +
  ggtitle("partition 4 distance to owner origin")+
  ylab("Dist (bp)")+
  xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank())
################################################################################################################
aatst <- read.table("~/Documents/Shayan/BioInf/lncRNA/partition_4_tile_list.txt", stringsAsFactors = F)

write.table(Partition_4_dfs$tile_name, 
            row.names = F, 
            col.names = F,
            quote = F, 
            file = "~/Documents/Shayan/BioInf/lncRNA/partition_4_tile_list.txt")
################################################################################################################
# partition 4
# Chunk partitioning
Partition_4_chunk <- list()
aa_df_list2 <- list()
aa_nu_start <- 40
aa_test_frac <- 0.1
aa_Valid_frac <- 0.2
set.seed(seed = 546456)
for(i in 1:length(Partition_4)){
  # choosing test pos
  print(i)
  aa_nu_pos <- length(Partition_4[[i]]$positive)
  aa_nu_test_pos <- floor(aa_nu_pos*aa_test_frac)
  aa_nu_vali_pos <- floor(aa_nu_pos*aa_Valid_frac)
  aa_nu_train_pos <- aa_nu_pos - (aa_nu_test_pos + aa_nu_vali_pos)
  
  aa_nu_neg <- length(Partition_4[[i]]$negative)
  aa_nu_test_neg <- floor(aa_nu_neg*aa_test_frac)
  aa_nu_vali_neg <- floor(aa_nu_neg*aa_Valid_frac)
  aa_nu_train_neg <- aa_nu_neg - (aa_nu_test_neg + aa_nu_vali_neg)
  
  # assign each negative to its closest positive
  aa_neg_pos_df <- data.frame(neg_index = Partition_4[[i]]$negative, 
                              closest_pos = numeric(length = aa_nu_neg))
  for(aamuneg in 1:aa_nu_neg){
    aa_neg_pos_df$closest_pos[aamuneg] <- Partition_4[[i]]$positive[which.min(abs(Partition_4[[i]]$negative[aamuneg] - Partition_4[[i]]$positive))]
  }
  
  # divide the positives into chunks based on index
  aa_pos_chunk <- list()
  aa_poschunk_size <- floor(aa_nu_pos/aa_nu_start)
  aapos_assignment <- cbind(Partition_4[[i]]$positive,
                            cut(c(1:aa_nu_pos),
                                labels = F,
                                breaks = aa_nu_start))
  aaneg_assignment <- cbind(Partition_4[[i]]$negative, 
                            aapos_assignment[match(aa_neg_pos_df$closest_pos,
                                                   aapos_assignment[, 1]), 2])
  # choosing parts for each group
  aa_ne_prop_error <- 0.5
  aa_rep <- 1
  aa_ne_prop_error_thr <- 0.02
  aathr_step <- 0.04
  while(aa_ne_prop_error > aa_ne_prop_error_thr){
    aatest_parts <- sample(x = c(1:aa_nu_start), size = floor(aa_nu_start * aa_test_frac))
    aavalid_parts <- sample(x =setdiff(c(1:aa_nu_start), aatest_parts), size = floor(aa_nu_start * aa_Valid_frac))
    aatrain_parts <- setdiff(setdiff(c(1:aa_nu_start), aatest_parts), aavalid_parts)
    
    
    # forming the dataframes
    aa_Sets <- c("train", "valid", "test")
    aa_pos_frame_set <- character(length = length(Partition_4[[i]]$positive))
    aa_neg_frame_set <- character(length = length(Partition_4[[i]]$negative))
    
    aa_pos_frame_set[which(aapos_assignment[, 2] %in% aatest_parts)] <- aa_Sets[3]
    aa_pos_frame_set[which(aapos_assignment[, 2] %in% aavalid_parts)] <- aa_Sets[2]
    aa_pos_frame_set[which(aapos_assignment[, 2] %in% aatrain_parts)] <- aa_Sets[1]
    aa_neg_frame_set[which(aaneg_assignment[, 2] %in% aatest_parts)] <- aa_Sets[3]
    aa_neg_frame_set[which(aaneg_assignment[, 2] %in% aavalid_parts)] <- aa_Sets[2]
    aa_neg_frame_set[which(aaneg_assignment[, 2] %in% aatrain_parts)] <- aa_Sets[1]
    
    aa_df_pos <- data.frame(tile_name = rownames(territory_assigned_interaction_3Mfilered)[Partition_4[[i]]$positive],
                            owner = rep(names(Partition_4)[i], length(Partition_4[[i]]$positive)), 
                            label = rep(1, length(Partition_4[[i]]$positive)), 
                            dataset = aa_pos_frame_set)
    aa_df_neg <- data.frame(tile_name = rownames(territory_assigned_interaction_3Mfilered)[Partition_4[[i]]$negative],
                            owner = rep(names(Partition_4)[i], length(Partition_4[[i]]$negative)), 
                            label = rep(0, length(Partition_4[[i]]$negative)),
                            dataset = aa_neg_frame_set)
    aa_neg_tab <- table(aa_df_neg$dataset)/sum(table(aa_df_neg$dataset))
    if(length(aa_neg_tab) < 3){
      aa_ne_prop_error <- 1
    }else{
      aa_ne_prop_error <- sqrt(mean((aa_neg_tab[c(1,3)] - c(aa_test_frac, aa_Valid_frac))^2))
    }
    aa_rep <- aa_rep + 1
    if(aa_rep == 1000){
      print(paste0("changing thresh to ", aa_ne_prop_error_thr + aathr_step))
      aa_rep <- 1
      aa_ne_prop_error_thr <- aa_ne_prop_error_thr + aathr_step
    }
  }

  
  aa_df_list2[[i]] <- rbind(aa_df_pos, aa_df_neg)
  print(names(Partition_4)[i])
  print("neg")
  print(table(aa_df_neg$dataset))
  print(table(aa_df_neg$dataset)/sum(table(aa_df_neg$dataset)))
  print("pos")
  print(table(aa_df_pos$dataset))
  print(table(aa_df_pos$dataset)/sum(table(aa_df_pos$dataset)))
  
}
Partition_4_chunk <- aa_df_list2
names(Partition_4_chunk) <- names(Partition_4)
Partition_4_chunk_dfs <- do.call(rbind, Partition_4_chunk)
table(Partition_4_chunk$Malat1$dataset)

for(i in 1:length(Partition_4_chunk)){
  Partition_4_chunk[[i]]$tile_name <- as.character(levels(Partition_4_chunk[[i]]$tile_name)[as.numeric(Partition_4_chunk[[i]]$tile_name)])
}
Partition_4_chunk_dfs$tile_name <- levels(Partition_4_chunk_dfs$tile_name)[as.numeric(Partition_4_chunk_dfs$tile_name)]
#########################################################################################################
# gather features fro partition 4
save(list = c("Partition_4_dfs", "Partition_4_chunk_dfs"), file= "~/Documents/Shayan/BioInf/lncRNA/Partition_4_chunk_random.RData")

Triplex_feaures_partition4
Tile_kmer_features_partition4
ChIPATLAS_features_partition4 <- ChIPATLAS_features[match(Partition_4_dfs$tile_name, rownames(ChIPATLAS_features)),]
RBP_features_partition4 <- RBP_features[match(Partition_4_dfs$tile_name, rownames(RBP_features)),]
ChIPATLAS_features_partition4[is.na(ChIPATLAS_features_partition4)] <- 0
RBP_features_partition4[is.na(RBP_features_partition4)] <- 0


mESC_CHIPATLAS_lncRNA_features
kmer_features_lncRNA28_normalized
mESC_RBP_lncRNA_features

# expression RNAseq tile
mESC_Renlab_RNASeq_partition4  <- numeric(nrow(Partition_4_dfs))
names(mESC_Renlab_RNASeq_partition4) <- Partition_4_dfs$tile_name
aamESC_Renlab_tiled_lisRNAseq <-  mESC_Renlab_tiled_list$`GSM723776_RenLab-RNA-Seq-mESC-ZY28.bed`[mESC_Renlab_tiled_list$`GSM723776_RenLab-RNA-Seq-mESC-ZY28.bed`$tile %in% Partition_4_dfs$tile_name,]
mESC_Renlab_RNASeq_partition4[match(aamESC_Renlab_tiled_lisRNAseq$tile, Partition_4_dfs$tile_name)] <- aamESC_Renlab_tiled_lisRNAseq$`aadf_agg$score`

# expression cage tile
aamESC_CAGE_mm9_tiled <- mESC_CAGE_mm9_tiled[mESC_CAGE_mm9_tiled$tile %in% Partition_4_dfs$tile_name,]
mESC_CAGE_mm9_tiled_partition4 <- numeric(nrow(Partition_4_dfs))
mESC_CAGE_mm9_tiled_partition4[match(aamESC_CAGE_mm9_tiled$tile, Partition_4_dfs$tile_name)] <- aamESC_CAGE_mm9_tiled$`aadf_agg$CAGE_score`
names(mESC_CAGE_mm9_tiled_partition4) <- Partition_4_dfs$tile_name

# save(list = c("mESC_CAGE_mm9_tiled",
#               "mESC_Renlab_tiled_list"),
#      file = "~/Documents/Shayan/BioInf/lncRNA/CAGE_RNAseq_tiled.RData")
# expression RNAseq lncRNA
# mESC_RNAseq_lncRNA_features
# mESC_RNAseq_lncRNA_features_full <- numeric(length = nrow(lncRNA_chosen_gt1k_uniqTiles))
# names(mESC_RNAseq_lncRNA_features_full) <- lncRNA_chosen_gt1k_uniqTiles$gene_name
# mESC_RNAseq_lncRNA_features_full[match(names(mESC_RNAseq_lncRNA_features), names(mESC_RNAseq_lncRNA_features_full))] <- mESC_RNAseq_lncRNA_features



# expression cage lncRNA
# mESC_CAGE_lncRNA_features
# mESC_CAGE_lncRNA_features_full <- numeric(length = nrow(lncRNA_chosen_gt1k_uniqTiles))
# names(mESC_CAGE_lncRNA_features_full) <- lncRNA_chosen_gt1k_uniqTiles$gene_name
# mESC_CAGE_lncRNA_features_full[match(names(mESC_CAGE_lncRNA_features), names(mESC_CAGE_lncRNA_features_full))] <- mESC_CAGE_lncRNA_features

# save(list = c("mESC_CAGE_lncRNA_features_full",
#               "mESC_RNAseq_lncRNA_features_full", 
#               "mESC_RBP_lncRNA_features", 
#               "kmer_features_lncRNA28_normalized",
#               "mESC_CHIPATLAS_lncRNA_features"),
#      file = "~/Documents/Shayan/BioInf/lncRNA/lncRNA_feaures.RData")


# Partition_4_dfs$tile_name <- levels(Partition_2_dfs$tile_name)[as.numeric(Partition_2_dfs$tile_name)]
# Partition_4_chunk_dfs$tile_name <- levels(Partition_2_chunk_dfs$tile_name)[as.numeric(Partition_2_chunk_dfs$tile_name)]
# 
Partition_4_feature_mat_tile <- cbind(ChIPATLAS_features_partition4, # remember to add repeats
                                      Tile_kmer_features_partition4,
                                      RBP_features_partition4,
                                      mESC_Renlab_RNASeq_partition4,
                                      mESC_CAGE_mm9_tiled_partition4)
colnames(Partition_4_feature_mat_tile)[ncol(Partition_4_feature_mat_tile) - 1] <- "RNAseq"
colnames(Partition_4_feature_mat_tile)[ncol(Partition_4_feature_mat_tile)] <- "CAGE"

Partition_4_feature_mat_owner <- cbind(mESC_CHIPATLAS_lncRNA_features[match(Partition_4_dfs$owner, rownames(mESC_CHIPATLAS_lncRNA_features)),],
                                       kmer_features_lncRNA28_normalized[match(Partition_4_dfs$owner, rownames(kmer_features_lncRNA28_normalized)),],
                                       mESC_RBP_lncRNA_features[match(Partition_4_dfs$owner, rownames(mESC_RBP_lncRNA_features)),],
                                       mESC_RNAseq_lncRNA_features_full[match(Partition_4_dfs$owner, names(mESC_RNAseq_lncRNA_features_full))],
                                       mESC_CAGE_lncRNA_features_full[match(Partition_4_dfs$owner, names(mESC_CAGE_lncRNA_features_full))])
colnames(Partition_4_feature_mat_owner)[ncol(Partition_4_feature_mat_owner) - 1] <- "RNAseq"
colnames(Partition_4_feature_mat_owner)[ncol(Partition_4_feature_mat_owner)] <- "CAGE"


colnames(Partition_4_feature_mat_owner) <- paste(colnames(Partition_4_feature_mat_owner), "lncRNA", sep = "__")

#load("~/Documents/Shayan/BioInf/lncRNA/Tile_owner_pair_RBP_ChIP_pair.RData")
Partition_4_feature_mat_pair <- cbind(distance_feature_partition_4,
                                      Triplex_feaures_partition4[match(Partition_4_dfs$tile_name , Triplex_feaures_partition4$tile_name), 3],
                                      RBP_tile_owner_pair_feature[match(Partition_4_dfs$tile_name , names(RBP_tile_owner_pair_feature))], 
                                      ChIP_tile_owner_pair_feature[match(Partition_4_dfs$tile_name ,  names(ChIP_tile_owner_pair_feature))],
                                      Chromatin_tile_owner_pair_feature[match(Partition_4_dfs$tile_name ,  names(Chromatin_tile_owner_pair_feature))])
colnames(Partition_4_feature_mat_pair) <- c("distance", "triplex", "RBP_pair", "ChIP_pair", "Chromatin_pair")

# how to remove - from names
aatst <- colnames(ChIPATLAS_features_partition4)
aatst2 <- gsub(pattern = "-", replacement = "", x = aatst)

aa_all_feat <- cbind(Partition_4_feature_mat_tile, Partition_4_feature_mat_owner, Partition_4_feature_mat_pair, Partition_4_dfs$label)
colnames(aa_all_feat)[ncol(aa_all_feat)] <- "label"

aa_zerovar <- apply(aa_all_feat, MARGIN = 2, FUN = var)
Zero_variance_columns <- colnames(aa_all_feat)[aa_zerovar == 0]
aa_all_feat <- aa_all_feat[, -which(aa_zerovar == 0)]

aadescrCor <- cor(aa_all_feat[, 1:(ncol(aa_all_feat) - 1)])
aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .90)
Highly_correlating_columns <- colnames(aa_all_feat)[aahighlyCorDescr]
aa_all_feat <- aa_all_feat[,-aahighlyCorDescr]

aa_all_feat <- as.data.frame(aa_all_feat)
aa_all_feat$label[aa_all_feat$label == 0] <- "neg"
aa_all_feat$label[aa_all_feat$label == 1] <- "pos"

aa_all_feat$label <- factor(aa_all_feat$label , levels = c("neg", "pos"))

aa_name_dic_p4 <- cbind(colnames(aa_all_feat), paste0("feature_", c(1:ncol(aa_all_feat))))
aa_name_dic_p4[ncol(aa_all_feat), 2] <- aa_name_dic_p4[ncol(aa_all_feat), 1] 
colnames(aa_all_feat) <- aa_name_dic_p4[, 2]

Partition_4_modified_dataset <- aa_all_feat

#########################################################################################################
# Partition_5 --> use the same strategy as partition 4 but filtere negatives such that there are not more than nu_closest_neg_thresh negatives whose closest positive is the same

Partition_5 <- list()
aa_ovl <- findOverlaps(query = MM9_1kb_tiled_GR_filtered, subject = lncRNA_chosen_gt1k_uniqTiles_GR)
aa_thresh <- 2000
neg_to_pos_ratio <- 4
nu_closest_neg_thresh <- 100 # number of negatives with the same closest positive
aaallhist <- hist(distance_feature_allTiles, breaks = 5000, freq = F)
set.seed(seed = 57668)
for(i in 1:ncol(territory_assigned_interaction_label_3Mfiltered)){
  print(i)
  aa_w_pos <- which(territory_assigned_interaction_label_3Mfiltered[, i] == 1)
  aa_w_pos <- setdiff(aa_w_pos,unique(aa_ovl@from))
  aa_w_Terr <- which(territory_assigned_interaction_3Mfilered[, i] == 1)
  aa_w_neg <- setdiff(aa_w_Terr,aa_w_pos)
  aa_w_neg <- setdiff(aa_w_neg,unique(aa_ovl@from))
  aa_w_neg <- aa_w_neg[rownames(territory_assigned_interaction_label_3Mfiltered)[aa_w_neg] %in% names(MM9_1kb_tiled_owner_3Mfilter_noNA)]
  
  # compute the closest positive for each negative of this lncRNA
  # assign each negative to its closest positive
  aa_neg_pos_df <- data.frame(neg_index = aa_w_neg, 
                              closest_pos = numeric(length = length(aa_w_neg)))
  for(aamuneg in 1:length(aa_w_neg)){
    aa_neg_pos_df$closest_pos[aamuneg] <- aa_w_pos[which.min(abs(aa_w_neg[aamuneg] - aa_w_pos))]
  }
  # allow each positive to have maximum of 100 negative points who are closest to it
  aatab <- table(aa_neg_pos_df$closest_pos)
  aa_w_pos_ABOVElimit <- as.numeric(names(aatab)[which(aatab > nu_closest_neg_thresh)])
  aa_toberemoved <- numeric(0)
  for(aacur_lim in 1:length(aa_w_pos_ABOVElimit)){
    aa_curlimneg <- aa_neg_pos_df$neg_index[aa_neg_pos_df$closest_pos %in% aa_w_pos_ABOVElimit[aacur_lim]]
    aarem <- aa_curlimneg[sort(abs(aa_curlimneg - aa_w_pos_ABOVElimit[aacur_lim]), index.return = T)$ix[(nu_closest_neg_thresh+1):length(aa_curlimneg)]]
    aa_toberemoved <- c(aa_toberemoved, aarem)
  }
  print("length(aa_w_neg)")
  print(length(aa_w_neg))
  print("length(aa_toberemoved)")
  print(length(aa_toberemoved))
  aa_w_neg <- setdiff(aa_w_neg,aa_toberemoved)
  
  
  
  if(length(aa_w_pos) > aa_thresh){
    aamyneg_dist_hist <- hist(distance_feature_allTiles[match(rownames(territory_assigned_interaction_label_3Mfiltered)[aa_w_neg],
                                                              names(distance_feature_allTiles))],
                              breaks = aaallhist$breaks, probability = T)
    
    aamypos_dist <- distance_feature_allTiles[match(rownames(territory_assigned_interaction_label_3Mfiltered)[aa_w_pos],
                                                    names(distance_feature_allTiles))]
    aa_prob_sum <- 0
    aa_myk <- 10
    while(aa_prob_sum < aa_thresh){
      print("pos sampling")
      print(aa_myk)
      aa_density_smooth <- zoo::rollmean(aamyneg_dist_hist$counts, k = aa_myk, fill = "extend")
      aa_myprob_smooth <- aa_density_smooth[findInterval(x = aamypos_dist, vec =  aaallhist$breaks)]
      aa_prob_sum <- sum(aa_myprob_smooth > 0)
      aa_myk <- aa_myk + 10
      stopifnot(aa_myk < (length(aa_myprob)/2))
    }
    aa_myprob_smooth[aa_myprob_smooth <0] <- 0
    
    aa_w_pos_sampled <- sort(sample(x = aa_w_pos, size = aa_thresh, replace = F, prob = aa_myprob_smooth))
  }else{
    aa_w_pos_sampled <- sort(aa_w_pos)
  }
  if(length(aa_w_neg) > (neg_to_pos_ratio * length(aa_w_pos_sampled))){
    aamypos_dist_hist <- hist(distance_feature_allTiles[match(rownames(territory_assigned_interaction_label_3Mfiltered)[aa_w_pos_sampled],
                                                              names(distance_feature_allTiles))], breaks = aaallhist$breaks, probability = T)
    
    aamyneg_dist <- distance_feature_allTiles[match(rownames(territory_assigned_interaction_label_3Mfiltered)[aa_w_neg],
                                                    names(distance_feature_allTiles))]
    aa_prob_sum <- 0
    aa_myk <- 10
    while(aa_prob_sum < (neg_to_pos_ratio * length(aa_w_pos_sampled))){
      print("neg sampling")
      print(aa_myk)
      aa_density_smooth <- zoo::rollmean(aamypos_dist_hist$counts, k = aa_myk, fill = "extend")
      aa_myprob_smooth <- aa_density_smooth[findInterval(x = aamyneg_dist, vec =  aaallhist$breaks)]
      aa_prob_sum <- sum(aa_myprob_smooth > 0)
      aa_myk <- aa_myk + 10
      stopifnot(aa_myk < (length(aa_myprob)/2))
    }
    aa_myprob_smooth[aa_myprob_smooth <0] <- 0
    aa_w_neg_sampled <- sort(sample(x = aa_w_neg,
                                    size = (neg_to_pos_ratio * length(aa_w_pos_sampled)),
                                    replace = F,
                                    prob = aa_myprob_smooth))
  }else{
    aa_w_neg_sampled <- aa_w_neg
  }
  Partition_5[[i]] <- list(positive=aa_w_pos_sampled, 
                           negative=aa_w_neg_sampled,
                           extra_positive = setdiff(aa_w_pos,aa_w_pos_sampled), 
                           extra_negative = setdiff(aa_w_neg,aa_w_neg_sampled))
  
}

names(Partition_5) <- colnames(territory_assigned_interaction_label_3Mfiltered)

aa_ps <- lapply(Partition_5, "[[", 1)
aa_ng <- lapply(Partition_5, "[[", 2)
aa_data_cnt <- rbind(unlist(lapply(aa_ps, length)),unlist(lapply(aa_ng, length)))
barplot(aa_data_cnt, las = 2, legend.text = c("+", "-")   ,
        main = "data points per lncRNA partition 5",
        args.legend=list(x = "topright", bty = "n", 
                         inset=c(0, 0),
                         cex = 1, 
                         y.intersp = 0.6,
                         x.intersp = 0.5,
                         text.width = 1))

aa_ps2 <- lapply(Partition_5, "[[", 3)
aa_ng2 <- lapply(Partition_5, "[[", 4)
aa_data_cnt2 <- rbind(unlist(lapply(aa_ps2, length)),unlist(lapply(aa_ng2, length)))
barplot(aa_data_cnt2, las = 2, legend.text = c("+", "-")   ,
        main = "extra data points per lncRNA",
        args.legend=list(x = "topright", bty = "n", 
                         inset=c(0, 0),
                         cex = 1, 
                         y.intersp = 0.6,
                         x.intersp = 0.5,
                         text.width = 1))


aachr_vecs2 <- list()
for(i in 1:nrow(aachr_size2)){
  aachr_vecs2[[i]] <-   array(dim =  ceiling(aachr_size2$V2[i]/1000))
  names(aachr_vecs2[[i]]) <- c(1:length(aachr_vecs2[[i]]))
  #aachr_vecs[[i]][aa_tilenamesp2[which(aa_tilenamesp1 %in% aachr_size2$V1[i])]] <- 1
}
names(aachr_vecs2) <- aachr_size2$V1

for(i in 1:length(aachr_vecs2)){
  aaw_lnc <- lncRNA_chosen_gt1k_uniqTiles$gene_name[lncRNA_chosen_gt1k_uniqTiles$chromosome_name  == names(aachr_vecs2)[i]]
  aa_wp <- which(names(Partition_5) %in% aaw_lnc)
  for(j in 1:length(aa_wp)){
    aa_tilename_pos <- rownames(territory_assigned_interaction_3Mfilered)[Partition_5[[aa_wp[j]]]$positive]
    aa_tilename_pos2<- as.numeric(unlist(lapply(strsplit(aa_tilename_pos, "_"), "[[", 2)))
    aa_tilename_neg <- rownames(territory_assigned_interaction_3Mfilered)[Partition_5[[aa_wp[j]]]$negative]
    aa_tilename_neg2<- as.numeric(unlist(lapply(strsplit(aa_tilename_neg, "_"), "[[", 2)))
    
    aachr_vecs2[[i]][aa_tilename_pos2] <- 1.25 + rnorm(n = length(aa_tilename_pos2), mean = 0, sd = 0.01)
    aachr_vecs2[[i]][aa_tilename_neg2] <- 0.25 + rnorm(n = length(aa_tilename_neg2), mean = 0, sd = 0.01)
    
    
  }
}

png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/chr_pos_neg_chosen_partition5.png",       
    width = 40*300,        # 5 x 300 pixels
    height = 30*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10) 
par(mfrow = c(16,1), mar = c(3,4,3,4))
for(i in 1:length(aachr_vecs)){
  print(i)
  plot(aachr_vecs[[i]], main = paste(names(aachr_vecs)[i], aachr_size2$V2[i]),
       las = 2, xlab = "genomic position", xaxt = "n", ylim = c(0,3.1), pch = 19, cex = 0.2)
  points(aachr_vecs2[[i]], pch = 8, col = 2, cex = 0.5)
}
dev.off()


aa_df_list <- list()
for(i in 1:length(Partition_5)){
  aa_df_pos <- data.frame(tile_name = rownames(territory_assigned_interaction_3Mfilered)[Partition_5[[i]]$positive],
                          owner = rep(names(Partition_5)[i], length(Partition_5[[i]]$positive)), 
                          label = rep(1, length(Partition_5[[i]]$positive)), 
                          dataset = sample(c("train", "valid", "test"), 
                                           size = length(Partition_5[[i]]$positive),
                                           replace = T, 
                                           prob = c(0.7, 0.2, 0.1)))
  aa_df_neg <- data.frame(tile_name = rownames(territory_assigned_interaction_3Mfilered)[Partition_5[[i]]$negative],
                          owner = rep(names(Partition_5)[i], length(Partition_5[[i]]$negative)), 
                          label = rep(0, length(Partition_5[[i]]$negative)),
                          dataset = sample(c("train", "valid", "test"), 
                                           size = length(Partition_5[[i]]$negative),
                                           replace = T, 
                                           prob = c(0.7, 0.2, 0.1))
  )
  
  aa_df_list[[i]] <- rbind(aa_df_pos, aa_df_neg)
  
}

names(aa_df_list) <- names(Partition_5)

Partition_5_dfs <- do.call(rbind, aa_df_list)
Partition_5_df_list <- aa_df_list
for(i in 1:length(Partition_5_df_list)){
  Partition_5_df_list[[i]]$tile_name <- as.character(levels(Partition_5_df_list[[i]]$tile_name)[as.numeric(Partition_5_df_list[[i]]$tile_name)])
}
Partition_5_dfs$tile_name <- as.character(levels(Partition_5_dfs$tile_name)[as.numeric(Partition_5_dfs$tile_name)])

aadata_p6 <- cbind(Partition_5_dfs, 
                   Partition_5_chunk_dfs$dataset,
                   distance_feature_partition_5)

colnames(aadata_p6) <- c("tile_name", "owner",     "label" ,  "dataset_random", "dataset_chunk","distance" )
aadata_p6$label[aadata_p6$label == 1] <- "pos"
aadata_p6$label[aadata_p6$label == 0] <- "neg"
aadata_p6$label <- factor(aadata_p6$label, levels = c("neg", "pos"))



ggplot(aadata_p6, aes(x=owner, y=distance, fill=label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot() +
  geom_vline(xintercept = seq(1.5,28.5,1), color ="blue", size = 1.5) +
  ggtitle("Partition 5 points distance to owner origin by lncRNA")+
  ylab("Dist (bp)")+
  xlab("lncRNA") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90)) 

ggplot(aadata_p6, aes(x=dataset_random, y=distance, fill=label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),scale = "count") +
  ggtitle("partition_5_random distance to owner origin")+
  ylab("Dist (bp)")+
  xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank())

ggplot(aadata_p6, aes(x=dataset_chunk, y=distance, fill=label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75) ,scale = "count") +
  ggtitle("partition_5_chunk distance to owner origin")+
  ylab("Dist (bp)")+
  xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank())

ggplot(aadata_p6, aes(x=label, y=distance+1)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), scale =  "count") +
  ggtitle("partition 5 distance to owner origin")+
  ylab("Dist (bp)")+
  xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank())

ggplot(aadata_p3, aes(x=label, y=distance+1)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), scale =  "count") +
  ggtitle("all points distance to lncRNA origin")+
  ylab("Dist (bp)")+
  xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank())
################################################################################################################
write.table(Partition_5_dfs$tile_name, 
            row.names = F, 
            col.names = F,
            quote = F, 
            file = "~/Documents/Shayan/BioInf/lncRNA/partition_5_tile_list.txt")
################################################################################################################
# partition 5
# Chunk partitioning
Partition_5_chunk <- list()
aa_df_list2 <- list()
aa_nu_start <- 40
aa_test_frac <- 0.1
aa_Valid_frac <- 0.2
set.seed(seed = 78543)
for(i in 1:length(Partition_5)){
  # choosing test pos
  print(i)
  aa_nu_pos <- length(Partition_5[[i]]$positive)
  aa_nu_test_pos <- floor(aa_nu_pos*aa_test_frac)
  aa_nu_vali_pos <- floor(aa_nu_pos*aa_Valid_frac)
  aa_nu_train_pos <- aa_nu_pos - (aa_nu_test_pos + aa_nu_vali_pos)
  
  aa_nu_neg <- length(Partition_5[[i]]$negative)
  aa_nu_test_neg <- floor(aa_nu_neg*aa_test_frac)
  aa_nu_vali_neg <- floor(aa_nu_neg*aa_Valid_frac)
  aa_nu_train_neg <- aa_nu_neg - (aa_nu_test_neg + aa_nu_vali_neg)
  
  # assign each negative to its closest positive
  aa_neg_pos_df <- data.frame(neg_index = Partition_5[[i]]$negative, 
                              closest_pos = numeric(length = aa_nu_neg))
  for(aamuneg in 1:aa_nu_neg){
    aa_neg_pos_df$closest_pos[aamuneg] <- Partition_5[[i]]$positive[which.min(abs(Partition_5[[i]]$negative[aamuneg] - Partition_5[[i]]$positive))]
  }
  
  # divide the positives into chunks based on index
  aa_pos_chunk <- list()
  aa_poschunk_size <- floor(aa_nu_pos/aa_nu_start)
  aapos_assignment <- cbind(Partition_5[[i]]$positive,
                            cut(c(1:aa_nu_pos),
                                labels = F,
                                breaks = aa_nu_start))
  aaneg_assignment <- cbind(Partition_5[[i]]$negative, 
                            aapos_assignment[match(aa_neg_pos_df$closest_pos,
                                                   aapos_assignment[, 1]), 2])
  # choosing parts for each group
  aa_ne_prop_error <- 0.5
  aa_rep <- 1
  aa_ne_prop_error_thr <- 0.02
  aathr_step <- 0.04
  while(aa_ne_prop_error > aa_ne_prop_error_thr){
    aatest_parts <- sample(x = c(1:aa_nu_start), size = floor(aa_nu_start * aa_test_frac))
    aavalid_parts <- sample(x =setdiff(c(1:aa_nu_start), aatest_parts), size = floor(aa_nu_start * aa_Valid_frac))
    aatrain_parts <- setdiff(setdiff(c(1:aa_nu_start), aatest_parts), aavalid_parts)
    
    
    # forming the dataframes
    aa_Sets <- c("train", "valid", "test")
    aa_pos_frame_set <- character(length = length(Partition_5[[i]]$positive))
    aa_neg_frame_set <- character(length = length(Partition_5[[i]]$negative))
    
    aa_pos_frame_set[which(aapos_assignment[, 2] %in% aatest_parts)] <- aa_Sets[3]
    aa_pos_frame_set[which(aapos_assignment[, 2] %in% aavalid_parts)] <- aa_Sets[2]
    aa_pos_frame_set[which(aapos_assignment[, 2] %in% aatrain_parts)] <- aa_Sets[1]
    aa_neg_frame_set[which(aaneg_assignment[, 2] %in% aatest_parts)] <- aa_Sets[3]
    aa_neg_frame_set[which(aaneg_assignment[, 2] %in% aavalid_parts)] <- aa_Sets[2]
    aa_neg_frame_set[which(aaneg_assignment[, 2] %in% aatrain_parts)] <- aa_Sets[1]
    
    aa_df_pos <- data.frame(tile_name = rownames(territory_assigned_interaction_3Mfilered)[Partition_5[[i]]$positive],
                            owner = rep(names(Partition_5)[i], length(Partition_5[[i]]$positive)), 
                            label = rep(1, length(Partition_5[[i]]$positive)), 
                            dataset = aa_pos_frame_set)
    aa_df_neg <- data.frame(tile_name = rownames(territory_assigned_interaction_3Mfilered)[Partition_5[[i]]$negative],
                            owner = rep(names(Partition_5)[i], length(Partition_5[[i]]$negative)), 
                            label = rep(0, length(Partition_5[[i]]$negative)),
                            dataset = aa_neg_frame_set)
    aa_neg_tab <- table(aa_df_neg$dataset)/sum(table(aa_df_neg$dataset))
    if(length(aa_neg_tab) == 3){
      aa_ne_prop_error <- sqrt(mean((aa_neg_tab[c(1,3)] - c(aa_test_frac, aa_Valid_frac))^2))
    }
    
    aa_rep <- aa_rep + 1
    if(aa_rep == 1000){
      print(paste0("changing thresh to ", aa_ne_prop_error_thr + aathr_step))
      aa_rep <- 1
      aa_ne_prop_error_thr <- aa_ne_prop_error_thr + aathr_step
    }

  }
  
  
  aa_df_list2[[i]] <- rbind(aa_df_pos, aa_df_neg)
  print(names(Partition_5)[i])
  print("neg")
  print(table(aa_df_neg$dataset))
  print(table(aa_df_neg$dataset)/sum(table(aa_df_neg$dataset)))
  print("pos")
  print(table(aa_df_pos$dataset))
  print(table(aa_df_pos$dataset)/sum(table(aa_df_pos$dataset)))
  
}
Partition_5_chunk <- aa_df_list2
names(Partition_5_chunk) <- names(Partition_5)
Partition_5_chunk_dfs <- do.call(rbind, Partition_5_chunk)
table(Partition_5_chunk$Malat1$dataset)

for(i in 1:length(Partition_5_chunk)){
  Partition_5_chunk[[i]]$tile_name <- as.character(levels(Partition_5_chunk[[i]]$tile_name)[as.numeric(Partition_5_chunk[[i]]$tile_name)])
}
Partition_5_chunk_dfs$tile_name <- levels(Partition_5_chunk_dfs$tile_name)[as.numeric(Partition_5_chunk_dfs$tile_name)]
#########################################################################################################
# gather features fro partition 5
save(list = c("Partition_5_dfs", "Partition_5_chunk_dfs"), file= "~/Documents/Shayan/BioInf/lncRNA/Partition_5_chunk_random.RData")

Triplex_feaures_partition5
Tile_kmer_features_partition5
ChIPATLAS_features_partition5 <- ChIPATLAS_features[match(Partition_5_dfs$tile_name, rownames(ChIPATLAS_features)),]
RBP_features_partition5 <- RBP_features[match(Partition_5_dfs$tile_name, rownames(RBP_features)),]
ChIPATLAS_features_partition5[is.na(ChIPATLAS_features_partition5)] <- 0
RBP_features_partition5[is.na(RBP_features_partition5)] <- 0


mESC_CHIPATLAS_lncRNA_features
kmer_features_lncRNA28_normalized
mESC_RBP_lncRNA_features

# expression RNAseq tile
mESC_Renlab_RNASeq_partition5  <- numeric(nrow(Partition_5_dfs))
names(mESC_Renlab_RNASeq_partition5) <- Partition_5_dfs$tile_name
aamESC_Renlab_tiled_lisRNAseq <-  mESC_Renlab_tiled_list$`GSM723776_RenLab-RNA-Seq-mESC-ZY28.bed`[mESC_Renlab_tiled_list$`GSM723776_RenLab-RNA-Seq-mESC-ZY28.bed`$tile %in% Partition_5_dfs$tile_name,]
mESC_Renlab_RNASeq_partition5[match(aamESC_Renlab_tiled_lisRNAseq$tile, Partition_5_dfs$tile_name)] <- aamESC_Renlab_tiled_lisRNAseq$`aadf_agg$score`

# expression cage tile
aamESC_CAGE_mm9_tiled <- mESC_CAGE_mm9_tiled[mESC_CAGE_mm9_tiled$tile %in% Partition_5_dfs$tile_name,]
mESC_CAGE_mm9_tiled_partition5 <- numeric(nrow(Partition_5_dfs))
mESC_CAGE_mm9_tiled_partition5[match(aamESC_CAGE_mm9_tiled$tile, Partition_5_dfs$tile_name)] <- aamESC_CAGE_mm9_tiled$`aadf_agg$CAGE_score`
names(mESC_CAGE_mm9_tiled_partition5) <- Partition_5_dfs$tile_name

# save(list = c("mESC_CAGE_mm9_tiled",
#               "mESC_Renlab_tiled_list"),
#      file = "~/Documents/Shayan/BioInf/lncRNA/CAGE_RNAseq_tiled.RData")
# expression RNAseq lncRNA
# mESC_RNAseq_lncRNA_features
# mESC_RNAseq_lncRNA_features_full <- numeric(length = nrow(lncRNA_chosen_gt1k_uniqTiles))
# names(mESC_RNAseq_lncRNA_features_full) <- lncRNA_chosen_gt1k_uniqTiles$gene_name
# mESC_RNAseq_lncRNA_features_full[match(names(mESC_RNAseq_lncRNA_features), names(mESC_RNAseq_lncRNA_features_full))] <- mESC_RNAseq_lncRNA_features



# expression cage lncRNA
# mESC_CAGE_lncRNA_features
# mESC_CAGE_lncRNA_features_full <- numeric(length = nrow(lncRNA_chosen_gt1k_uniqTiles))
# names(mESC_CAGE_lncRNA_features_full) <- lncRNA_chosen_gt1k_uniqTiles$gene_name
# mESC_CAGE_lncRNA_features_full[match(names(mESC_CAGE_lncRNA_features), names(mESC_CAGE_lncRNA_features_full))] <- mESC_CAGE_lncRNA_features

mESC_CHIPATLAS_lncRNA_features[is.na(mESC_CHIPATLAS_lncRNA_features)] <- 0
mESC_RBP_lncRNA_features[is.na(mESC_RBP_lncRNA_features)] <- 0
save(list = c("mESC_CAGE_lncRNA_features_full",
              "mESC_RNAseq_lncRNA_features_full",
              "mESC_RBP_lncRNA_features",
              "kmer_features_lncRNA28_normalized",
              "mESC_CHIPATLAS_lncRNA_features"),
     file = "~/Documents/Shayan/BioInf/lncRNA/lncRNA_feaures.RData")


# Partition_4_dfs$tile_name <- levels(Partition_2_dfs$tile_name)[as.numeric(Partition_2_dfs$tile_name)]
# Partition_4_chunk_dfs$tile_name <- levels(Partition_2_chunk_dfs$tile_name)[as.numeric(Partition_2_chunk_dfs$tile_name)]
# 
Partition_5_feature_mat_tile <- cbind(ChIPATLAS_features_partition5, # remember to add repeats
                                      Tile_kmer_features_partition5,
                                      RBP_features_partition5,
                                      mESC_Renlab_RNASeq_partition5,
                                      mESC_CAGE_mm9_tiled_partition5)
colnames(Partition_5_feature_mat_tile)[ncol(Partition_5_feature_mat_tile) - 1] <- "RNAseq"
colnames(Partition_5_feature_mat_tile)[ncol(Partition_5_feature_mat_tile)] <- "CAGE"

Partition_5_feature_mat_owner <- cbind(mESC_CHIPATLAS_lncRNA_features[match(Partition_5_dfs$owner, rownames(mESC_CHIPATLAS_lncRNA_features)),],
                                       kmer_features_lncRNA28_normalized[match(Partition_5_dfs$owner, rownames(kmer_features_lncRNA28_normalized)),],
                                       mESC_RBP_lncRNA_features[match(Partition_5_dfs$owner, rownames(mESC_RBP_lncRNA_features)),],
                                       mESC_RNAseq_lncRNA_features_full[match(Partition_5_dfs$owner, names(mESC_RNAseq_lncRNA_features_full))],
                                       mESC_CAGE_lncRNA_features_full[match(Partition_5_dfs$owner, names(mESC_CAGE_lncRNA_features_full))])
colnames(Partition_5_feature_mat_owner)[ncol(Partition_5_feature_mat_owner) - 1] <- "RNAseq"
colnames(Partition_5_feature_mat_owner)[ncol(Partition_5_feature_mat_owner)] <- "CAGE"


colnames(Partition_5_feature_mat_owner) <- paste(colnames(Partition_5_feature_mat_owner), "lncRNA", sep = "__")

#load("~/Documents/Shayan/BioInf/lncRNA/Tile_owner_pair_RBP_ChIP_pair.RData")
Partition_5_feature_mat_pair <- cbind(distance_feature_partition_5,
                                      Triplex_feaures_partition5[match(Partition_5_dfs$tile_name , Triplex_feaures_partition5$tile_name), 3],
                                      RBP_tile_owner_pair_feature[match(Partition_5_dfs$tile_name , names(RBP_tile_owner_pair_feature))], 
                                      ChIP_tile_owner_pair_feature[match(Partition_5_dfs$tile_name ,  names(ChIP_tile_owner_pair_feature))],
                                      Chromatin_tile_owner_pair_feature[match(Partition_5_dfs$tile_name ,  names(Chromatin_tile_owner_pair_feature))])
colnames(Partition_5_feature_mat_pair) <- c("distance", "triplex", "RBP_pair", "ChIP_pair", "Chromatin_pair")


aa_all_feat <- cbind(Partition_5_feature_mat_tile, Partition_5_feature_mat_owner, Partition_5_feature_mat_pair, Partition_5_dfs$label)
colnames(aa_all_feat)[ncol(aa_all_feat)] <- "label"

aa_zerovar <- apply(aa_all_feat, MARGIN = 2, FUN = var)
Zero_variance_columns <- colnames(aa_all_feat)[aa_zerovar == 0]
aa_all_feat <- aa_all_feat[, -which(aa_zerovar == 0)]

aadescrCor <- cor(aa_all_feat[, 1:(ncol(aa_all_feat) - 1)])
aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .90)
Highly_correlating_columns <- colnames(aa_all_feat)[aahighlyCorDescr]
aa_all_feat <- aa_all_feat[,-aahighlyCorDescr]


aa_all_feat <- as.data.frame(aa_all_feat)
aa_all_feat$label[aa_all_feat$label == 0] <- "neg"
aa_all_feat$label[aa_all_feat$label == 1] <- "pos"

aa_all_feat$label <- factor(aa_all_feat$label , levels = c("neg", "pos"))

#aa_name_dic_p5 <- cbind(colnames(aa_all_feat), paste0("feature_", c(1:ncol(aa_all_feat))))
aa_name_dic_p5 <- cbind(colnames(aa_all_feat), paste0("feature_", c(1:ncol(aa_all_feat))))
aa_name_dic_p5[ncol(aa_all_feat), 2] <- aa_name_dic_p5[ncol(aa_all_feat), 1] 
colnames(aa_all_feat) <- aa_name_dic_p5[, 2]

Partition_5_modified_dataset <- aa_all_feat
####################################################################################################################################
# create partitions per lncRNA --> partition 6, running separate RFs for each lncRNA
# index_first_3000_perchr
# lncRNA_interaction_distance_matrix_full
# lncRNA_interaction_label_matrix
# colnames(lncRNA_interaction_label_matrix) <- lncRNA_chosen_gt1k_uniqTiles$gene_name[match(colnames(lncRNA_interaction_label_matrix),
#                                                                                           lncRNA_chosen_gt1k_uniqTiles$geneID)]
Partition_6 <- list()
aa_ovl <- findOverlaps(query = MM9_1kb_tiled_GR_filtered, subject = lncRNA_chosen_gt1k_uniqTiles_GR)
aa_thresh <- 200000
neg_to_pos_ratio <- 5

aartoremv <- rownames(lncRNA_interaction_distance_matrix_full)[index_first_3000_perchr]
lncRNA_interaction_distance_matrix_full_3MF <- lncRNA_interaction_distance_matrix_full[-index_first_3000_perchr,]
#lncRNA_interaction_label_matrix_3MF <- lncRNA_interaction_label_matrix[!rownames(lncRNA_interaction_label_matrix) %in% aartoremv,]
nu_closest_neg_thresh <- 200 # number of negatives with the same closest positive
#aaallhist <- hist(distance_feature_allTiles, breaks = 5000, freq = F)
set.seed(seed = 68521)
for(i in 1:ncol(lncRNA_interaction_label_matrix_3MF)){
  if(i == 1){ # for malat1 considering the non cis interactions
    aa_w_pos_name <- rownames(lncRNA_interaction_label_matrix_3MF)[which(lncRNA_interaction_label_matrix_3MF[, i] == 1)]
    aa_w_pos <- which(MM9_1kb_tiled_GR_filtered$tile %in% aa_w_pos_name)
    aa_w_pos <- setdiff(aa_w_pos, unique(aa_ovl@from))
    aa_w_neg <- unique(unlist((lapply(aa_w_pos, FUN = function(x) c((x -25):(x + 25))))))
    aa_w_neg <- setdiff(aa_w_neg,aa_w_pos)
    aa_w_neg <- setdiff(aa_w_neg,unique(aa_ovl@from))
    aa_w_neg <- aa_w_neg[aa_w_neg > 0]
    aa_w_neg <- aa_w_neg[aa_w_neg < length(MM9_1kb_tiled_GR_filtered)]
  }else{
    aa_w_pos_name <- rownames(lncRNA_interaction_label_matrix_3MF)[which(lncRNA_interaction_label_matrix_3MF[, i] == 1)]
    aa_w_pos <- which(MM9_1kb_tiled_GR_filtered$tile %in% aa_w_pos_name)
    aa_w_pos <- setdiff(aa_w_pos, unique(aa_ovl@from))
    aa_w_Terr <- which(!is.na(lncRNA_interaction_distance_matrix_full_3MF[, i]))
    aa_w_neg <- setdiff(aa_w_Terr,aa_w_pos)
    aa_w_neg <- setdiff(aa_w_neg,unique(aa_ovl@from))
    
  }
  print(i)
  #aa_w_neg <- aa_w_neg[rownames(territory_assigned_interaction_label_3Mfiltered)[aa_w_neg] %in% names(MM9_1kb_tiled_owner_3Mfilter_noNA)]
  
  # compute the closest positive for each negative of this lncRNA
  # assign each negative to its closest positive
  aa_neg_pos_df <- data.frame(neg_index = aa_w_neg, 
                              closest_pos = numeric(length = length(aa_w_neg)))
  for(aamuneg in 1:length(aa_w_neg)){
    aa_neg_pos_df$closest_pos[aamuneg] <- aa_w_pos[which.min(abs(aa_w_neg[aamuneg] - aa_w_pos))]
  }
  # allow each positive to have maximum of 100 negative points who are closest to it
  aatab <- table(aa_neg_pos_df$closest_pos)
  aa_w_pos_ABOVElimit <- as.numeric(names(aatab)[which(aatab > nu_closest_neg_thresh)])
  aa_toberemoved <- numeric(0)
  for(aacur_lim in 1:length(aa_w_pos_ABOVElimit)){
    aa_curlimneg <- aa_neg_pos_df$neg_index[aa_neg_pos_df$closest_pos %in% aa_w_pos_ABOVElimit[aacur_lim]]
    aarem <- aa_curlimneg[sort(abs(aa_curlimneg - aa_w_pos_ABOVElimit[aacur_lim]), index.return = T)$ix[(nu_closest_neg_thresh+1):length(aa_curlimneg)]]
    aa_toberemoved <- c(aa_toberemoved, aarem)
  }
  print("length(aa_w_neg)")
  print(length(aa_w_neg))
  print("length(aa_toberemoved)")
  print(length(aa_toberemoved))
  aa_w_neg <- setdiff(aa_w_neg,aa_toberemoved)
  
  
  
  if(length(aa_w_pos) > aa_thresh){
    stop("the threshold is too small, losing positive points")
    # aamyneg_dist_hist <- hist(distance_feature_allTiles[match(rownames(territory_assigned_interaction_label_3Mfiltered)[aa_w_neg],
    #                                                           names(distance_feature_allTiles))],
    #                           breaks = aaallhist$breaks, probability = T)
    # 
    # aamypos_dist <- distance_feature_allTiles[match(rownames(territory_assigned_interaction_label_3Mfiltered)[aa_w_pos],
    #                                                 names(distance_feature_allTiles))]
    # aa_prob_sum <- 0
    # aa_myk <- 10
    # while(aa_prob_sum < aa_thresh){
    #   print("pos sampling")
    #   print(aa_myk)
    #   aa_density_smooth <- zoo::rollmean(aamyneg_dist_hist$counts, k = aa_myk, fill = "extend")
    #   aa_myprob_smooth <- aa_density_smooth[findInterval(x = aamypos_dist, vec =  aaallhist$breaks)]
    #   aa_prob_sum <- sum(aa_myprob_smooth > 0)
    #   aa_myk <- aa_myk + 10
    #   stopifnot(aa_myk < (length(aa_myprob)/2))
    # }
    # aa_myprob_smooth[aa_myprob_smooth <0] <- 0
    # 
    # aa_w_pos_sampled <- sort(sample(x = aa_w_pos, size = aa_thresh, replace = F, prob = aa_myprob_smooth))
  }else{ # in this case there will be no sampling of positive points as all of them will be useed
    aa_w_pos_sampled <- sort(aa_w_pos)
  }
  if(length(aa_w_pos_sampled) > 100000){
    neg_to_pos_ratio <- 2
  }else{
    neg_to_pos_ratio <- 5
  }
  if(length(aa_w_neg) > (neg_to_pos_ratio * length(aa_w_pos_sampled))){
    if(i > 1){
      aa_allhist <- hist(lncRNA_interaction_distance_matrix_full[, i], breaks = 1000, probability = T)
      
      aamypos_dist_hist <- hist(lncRNA_interaction_distance_matrix_full[aa_w_pos_sampled, i], breaks = aa_allhist$breaks, probability = T)
      
      aamyneg_dist <- hist(lncRNA_interaction_distance_matrix_full[aa_w_neg, i], breaks = aa_allhist$breaks, probability = T)
    

    aa_prob_sum <- 0
    aa_myk <- 10
    while(aa_prob_sum < (neg_to_pos_ratio * length(aa_w_pos_sampled))){
      print("neg sampling")
      print(aa_myk)
      aa_density_smooth <- zoo::rollmean(aamypos_dist_hist$counts, k = aa_myk, fill = "extend")
      aa_myprob_smooth <- aa_density_smooth[findInterval(x = aamyneg_dist, vec =  aa_allhist$breaks)]
      aa_prob_sum <- sum(aa_myprob_smooth > 0)
      aa_myk <- aa_myk + 10
      stopifnot(aa_myk < (length(aa_myprob)/2))
    }
    aa_myprob_smooth[aa_myprob_smooth <0] <- 0
    aa_w_neg_sampled <- sort(sample(x = aa_w_neg,
                                    size = (neg_to_pos_ratio * length(aa_w_pos_sampled)),
                                    replace = F,
                                    prob = aa_myprob_smooth))
    }else{
      print("There are too many negatives for malat, think of a way to cut them")
      aa_w_neg_sampled <- aa_w_neg
    }
  }else{
    aa_w_neg_sampled <- aa_w_neg
  }
  Partition_6[[i]] <- list(positive=aa_w_pos_sampled, 
                           negative=aa_w_neg_sampled,
                           extra_positive = setdiff(aa_w_pos,aa_w_pos_sampled), 
                           extra_negative = setdiff(aa_w_neg,aa_w_neg_sampled))
  
}

names(Partition_6) <- colnames(lncRNA_interaction_distance_matrix_full)

aa_ps <- lapply(Partition_6, "[[", 1)
aa_ng <- lapply(Partition_6, "[[", 2)
aa_data_cnt <- rbind(unlist(lapply(aa_ps, length)),unlist(lapply(aa_ng, length)))
barplot((aa_data_cnt), las = 2, legend.text = c("+", "-")   ,
        main = "data points per lncRNA partition 6",
        args.legend=list(x = "topright", bty = "n", 
                         inset=c(0, 0),
                         cex = 1, 
                         y.intersp = 0.6,
                         x.intersp = 0.5,
                         text.width = 1))

aa_ps2 <- lapply(Partition_6, "[[", 3)
aa_ng2 <- lapply(Partition_6, "[[", 4)
aa_data_cnt2 <- rbind(unlist(lapply(aa_ps2, length)),unlist(lapply(aa_ng2, length)))
barplot(aa_data_cnt2, las = 2, legend.text = c("+", "-")   ,
        main = "extra data points per lncRNA P6",
        args.legend=list(x = "topright", bty = "n", 
                         inset=c(0, 0),
                         cex = 1, 
                         y.intersp = 0.6,
                         x.intersp = 0.5,
                         text.width = 1))


aachr_vecs2 <- list()
for(i in 1:nrow(aachr_size2)){
  aachr_vecs2[[i]] <-   array(dim =  ceiling(aachr_size2$V2[i]/1000))
  names(aachr_vecs2[[i]]) <- c(1:length(aachr_vecs2[[i]]))
  #aachr_vecs[[i]][aa_tilenamesp2[which(aa_tilenamesp1 %in% aachr_size2$V1[i])]] <- 1
}
names(aachr_vecs2) <- aachr_size2$V1

for(i in 1:length(aachr_vecs2)){
  aaw_lnc <- lncRNA_chosen_gt1k_uniqTiles$gene_name[lncRNA_chosen_gt1k_uniqTiles$chromosome_name  == names(aachr_vecs2)[i]]
  aa_wp <- which(names(Partition_6) %in% aaw_lnc)
  for(j in 1:length(aa_wp)){
    aa_tilename_pos <- rownames(lncRNA_interaction_distance_matrix_full_3MF)[Partition_6[[aa_wp[j]]]$positive]
    aa_tilename_pos2<- as.numeric(unlist(lapply(strsplit(aa_tilename_pos, "_"), "[[", 2)))
    aa_tilename_neg <- rownames(lncRNA_interaction_distance_matrix_full_3MF)[Partition_6[[aa_wp[j]]]$negative]
    aa_tilename_neg2<- as.numeric(unlist(lapply(strsplit(aa_tilename_neg, "_"), "[[", 2)))
    
    aachr_vecs2[[i]][aa_tilename_pos2] <- 1.25 + rnorm(n = length(aa_tilename_pos2), mean = 0, sd = 0.01)
    aachr_vecs2[[i]][aa_tilename_neg2] <- 0.25 + rnorm(n = length(aa_tilename_neg2), mean = 0, sd = 0.01)
    
    
  }
}

png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/chr_pos_neg_chosen_partition6.png",       
    width = 40*300,        # 5 x 300 pixels
    height = 30*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10) 
par(mfrow = c(16,1), mar = c(3,4,3,4))
for(i in 1:length(aachr_vecs)){
  print(i)
  plot(aachr_vecs[[i]], main = paste(names(aachr_vecs)[i], aachr_size2$V2[i]),
       las = 2, xlab = "genomic position", xaxt = "n", ylim = c(0,3.1), pch = 19, cex = 0.2)
  points(aachr_vecs2[[i]], pch = 8, col = 2, cex = 0.5)
}
dev.off()


aa_df_list <- list()
for(i in 1:length(Partition_6)){
  aa_df_pos <- data.frame(tile_name = rownames(lncRNA_interaction_distance_matrix_full_3MF)[Partition_6[[i]]$positive],
                          owner = rep(names(Partition_6)[i], length(Partition_6[[i]]$positive)), 
                          label = rep(1, length(Partition_6[[i]]$positive)), 
                          dataset = sample(c("train", "valid", "test"), 
                                           size = length(Partition_6[[i]]$positive),
                                           replace = T, 
                                           prob = c(0.7, 0.2, 0.1)))
  aa_df_neg <- data.frame(tile_name = rownames(lncRNA_interaction_distance_matrix_full_3MF)[Partition_6[[i]]$negative],
                          owner = rep(names(Partition_6)[i], length(Partition_6[[i]]$negative)), 
                          label = rep(0, length(Partition_6[[i]]$negative)),
                          dataset = sample(c("train", "valid", "test"), 
                                           size = length(Partition_6[[i]]$negative),
                                           replace = T, 
                                           prob = c(0.7, 0.2, 0.1))
  )
  
  aa_df_list[[i]] <- rbind(aa_df_pos, aa_df_neg)
  
}

names(aa_df_list) <- names(Partition_6)

Partition_6_dfs <- do.call(rbind, aa_df_list)
Partition_6_df_list <- aa_df_list
for(i in 1:length(Partition_6_df_list)){
  Partition_6_df_list[[i]]$tile_name <- as.character(levels(Partition_6_df_list[[i]]$tile_name)[as.numeric(Partition_6_df_list[[i]]$tile_name)])
}
for(i in 1:length(Partition_6_df_list)){
  print(i)
  print(sum(duplicated(Partition_6_df_list[[i]]$tile_name)))
}
Partition_6_dfs$tile_name <- as.character(levels(Partition_6_dfs$tile_name)[as.numeric(Partition_6_dfs$tile_name)])
######
aadata_p6 <- cbind(Partition_6_dfs, 
                   Partition_6_chunk_dfs$dataset,
                   distance_feature_partition_6_cis)

colnames(aadata_p6) <- c("tile_name", "owner",     "label" ,  "dataset_random", "dataset_chunk","distance" )
aadata_p6$label[aadata_p6$label == 1] <- "pos"
aadata_p6$label[aadata_p6$label == 0] <- "neg"
aadata_p6$label <- factor(aadata_p6$label, levels = c("neg", "pos"))



ggplot(aadata_p6, aes(x=owner, y=distance + 1, fill=label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot() +
  geom_vline(xintercept = seq(1.5,28.5,1), color ="blue", size = 1.5) +
  ggtitle("Partition 6 points distance to owner origin by lncRNA")+
  ylab("Dist (bp)")+
  xlab("lncRNA") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90)) 

ggplot(aadata_p6, aes(x=dataset_random, y=distance+1, fill=label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),scale = "count") +
  ggtitle("partition_6_random distance to owner origin")+
  ylab("Dist (bp)")+
  xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank())

ggplot(aadata_p6, aes(x=dataset_chunk, y=distance+1, fill=label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75) ,scale = "count") +
  ggtitle("partition_6_chunk distance to owner origin")+
  ylab("Dist (bp)")+
  xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank())

ggplot(aadata_p6, aes(x=label, y=distance+1)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)
              #, scale =  "count"
              ) +
  ggtitle("distance to lncRNA origin\nfor selected data points")+
  ylab("Dist (bp)")+
  xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank())
################################################################################################################
write.table(Partition_6_dfs$tile_name, 
            row.names = F, 
            col.names = F,
            quote = F, 
            file = "~/Documents/Shayan/BioInf/lncRNA/partition_6_tile_list.txt")
################################################################################################################
# partition 6
# Chunk partitioning
Partition_6_chunk <- list()
aa_df_list2 <- list()
aa_nu_start <- 40
aa_test_frac <- 0.1
aa_Valid_frac <- 0.2
set.seed(seed = 86721)
for(i in 1:length(Partition_6)){
  # choosing test pos
  print(i)
  aa_nu_pos <- length(Partition_6[[i]]$positive)
  aa_nu_test_pos <- floor(aa_nu_pos*aa_test_frac)
  aa_nu_vali_pos <- floor(aa_nu_pos*aa_Valid_frac)
  aa_nu_train_pos <- aa_nu_pos - (aa_nu_test_pos + aa_nu_vali_pos)
  
  aa_nu_neg <- length(Partition_6[[i]]$negative)
  aa_nu_test_neg <- floor(aa_nu_neg*aa_test_frac)
  aa_nu_vali_neg <- floor(aa_nu_neg*aa_Valid_frac)
  aa_nu_train_neg <- aa_nu_neg - (aa_nu_test_neg + aa_nu_vali_neg)
  
  # assign each negative to its closest positive
  aa_neg_pos_df <- data.frame(neg_index = Partition_6[[i]]$negative, 
                              closest_pos = numeric(length = aa_nu_neg))
  for(aamuneg in 1:aa_nu_neg){
    aa_neg_pos_df$closest_pos[aamuneg] <- Partition_6[[i]]$positive[which.min(abs(Partition_6[[i]]$negative[aamuneg] - Partition_6[[i]]$positive))]
  }
  
  # divide the positives into chunks based on index
  aa_pos_chunk <- list()
  aa_poschunk_size <- floor(aa_nu_pos/aa_nu_start)
  aapos_assignment <- cbind(Partition_6[[i]]$positive,
                            cut(c(1:aa_nu_pos),
                                labels = F,
                                breaks = aa_nu_start))
  aaneg_assignment <- cbind(Partition_6[[i]]$negative, 
                            aapos_assignment[match(aa_neg_pos_df$closest_pos,
                                                   aapos_assignment[, 1]), 2])
  # choosing parts for each group
  aa_ne_prop_error <- 0.5
  aa_rep <- 1
  aa_ne_prop_error_thr <- 0.02
  aathr_step <- 0.04
  while(aa_ne_prop_error > aa_ne_prop_error_thr){
    aatest_parts <- sample(x = c(1:aa_nu_start), size = floor(aa_nu_start * aa_test_frac))
    aavalid_parts <- sample(x =setdiff(c(1:aa_nu_start), aatest_parts), size = floor(aa_nu_start * aa_Valid_frac))
    aatrain_parts <- setdiff(setdiff(c(1:aa_nu_start), aatest_parts), aavalid_parts)
    
    
    # forming the dataframes
    aa_Sets <- c("train", "valid", "test")
    aa_pos_frame_set <- character(length = length(Partition_6[[i]]$positive))
    aa_neg_frame_set <- character(length = length(Partition_6[[i]]$negative))
    
    aa_pos_frame_set[which(aapos_assignment[, 2] %in% aatest_parts)] <- aa_Sets[3]
    aa_pos_frame_set[which(aapos_assignment[, 2] %in% aavalid_parts)] <- aa_Sets[2]
    aa_pos_frame_set[which(aapos_assignment[, 2] %in% aatrain_parts)] <- aa_Sets[1]
    aa_neg_frame_set[which(aaneg_assignment[, 2] %in% aatest_parts)] <- aa_Sets[3]
    aa_neg_frame_set[which(aaneg_assignment[, 2] %in% aavalid_parts)] <- aa_Sets[2]
    aa_neg_frame_set[which(aaneg_assignment[, 2] %in% aatrain_parts)] <- aa_Sets[1]
    
    aa_df_pos <- data.frame(tile_name = rownames(lncRNA_interaction_distance_matrix_full_3MF)[Partition_6[[i]]$positive],
                            owner = rep(names(Partition_6)[i], length(Partition_6[[i]]$positive)), 
                            label = rep(1, length(Partition_6[[i]]$positive)), 
                            dataset = aa_pos_frame_set)
    aa_df_neg <- data.frame(tile_name = rownames(lncRNA_interaction_distance_matrix_full_3MF)[Partition_6[[i]]$negative],
                            owner = rep(names(Partition_6)[i], length(Partition_6[[i]]$negative)), 
                            label = rep(0, length(Partition_6[[i]]$negative)),
                            dataset = aa_neg_frame_set)
    aa_neg_tab <- table(aa_df_neg$dataset)/sum(table(aa_df_neg$dataset))
    if(length(aa_neg_tab) == 3){
      aa_ne_prop_error <- sqrt(mean((aa_neg_tab[c(1,3)] - c(aa_test_frac, aa_Valid_frac))^2))
    }
    
    aa_rep <- aa_rep + 1
    if(aa_rep == 1000){
      print(paste0("changing thresh to ", aa_ne_prop_error_thr + aathr_step))
      aa_rep <- 1
      aa_ne_prop_error_thr <- aa_ne_prop_error_thr + aathr_step
    }
    
  }
  
  
  aa_df_list2[[i]] <- rbind(aa_df_pos, aa_df_neg)
  print(names(Partition_6)[i])
  print("neg")
  print(table(aa_df_neg$dataset))
  print(table(aa_df_neg$dataset)/sum(table(aa_df_neg$dataset)))
  print("pos")
  print(table(aa_df_pos$dataset))
  print(table(aa_df_pos$dataset)/sum(table(aa_df_pos$dataset)))
  
}
Partition_6_chunk <- aa_df_list2
names(Partition_6_chunk) <- names(Partition_6)
Partition_6_chunk_dfs <- do.call(rbind, Partition_6_chunk)
table(Partition_6_chunk$Malat1$dataset)

for(i in 1:length(Partition_6_chunk)){
  Partition_6_chunk[[i]]$tile_name <- as.character(levels(Partition_6_chunk[[i]]$tile_name)[as.numeric(Partition_6_chunk[[i]]$tile_name)])
}
Partition_6_chunk_dfs$tile_name <- levels(Partition_6_chunk_dfs$tile_name)[as.numeric(Partition_6_chunk_dfs$tile_name)]
#########################################################################################################
# gather features fro partition 6
save(list = c("Partition_6_dfs", "Partition_6_chunk_dfs"),
     file= "~/Documents/Shayan/BioInf/lncRNA/Partition_6_chunk_random.RData")

Triplex_feaures_partition6
Tile_kmer_features_partition6
ChIPATLAS_features_partition6 <- ChIPATLAS_features[match(Partition_6_dfs$tile_name, rownames(ChIPATLAS_features)),]
RBP_features_partition6 <- RBP_features[match(Partition_6_dfs$tile_name, rownames(RBP_features)),]
ChIPATLAS_features_partition6[is.na(ChIPATLAS_features_partition6)] <- 0
RBP_features_partition6[is.na(RBP_features_partition6)] <- 0


# mESC_CHIPATLAS_lncRNA_features
# kmer_features_lncRNA28_normalized
# mESC_RBP_lncRNA_features

# expression RNAseq tile
mESC_Renlab_RNASeq_partition6  <- numeric(nrow(Partition_6_dfs))
names(mESC_Renlab_RNASeq_partition6) <- Partition_6_dfs$tile_name
aamESC_Renlab_tiled_lisRNAseq <-  mESC_Renlab_tiled_list$`GSM723776_RenLab-RNA-Seq-mESC-ZY28.bed`[mESC_Renlab_tiled_list$`GSM723776_RenLab-RNA-Seq-mESC-ZY28.bed`$tile %in% Partition_6_dfs$tile_name,]
mESC_Renlab_RNASeq_partition6[match(aamESC_Renlab_tiled_lisRNAseq$tile, Partition_6_dfs$tile_name)] <- aamESC_Renlab_tiled_lisRNAseq$`aadf_agg$score`

# expression cage tile
aamESC_CAGE_mm9_tiled <- mESC_CAGE_mm9_tiled[mESC_CAGE_mm9_tiled$tile %in% Partition_6_dfs$tile_name,]
mESC_CAGE_mm9_tiled_partition6 <- numeric(nrow(Partition_6_dfs))
mESC_CAGE_mm9_tiled_partition6[match(aamESC_CAGE_mm9_tiled$tile, Partition_6_dfs$tile_name)] <- aamESC_CAGE_mm9_tiled$`aadf_agg$CAGE_score`
names(mESC_CAGE_mm9_tiled_partition6) <- Partition_6_dfs$tile_name

# save(list = c("mESC_CAGE_mm9_tiled",
#               "mESC_Renlab_tiled_list"),
#      file = "~/Documents/Shayan/BioInf/lncRNA/CAGE_RNAseq_tiled.RData")
# expression RNAseq lncRNA
# mESC_RNAseq_lncRNA_features
# mESC_RNAseq_lncRNA_features_full <- numeric(length = nrow(lncRNA_chosen_gt1k_uniqTiles))
# names(mESC_RNAseq_lncRNA_features_full) <- lncRNA_chosen_gt1k_uniqTiles$gene_name
# mESC_RNAseq_lncRNA_features_full[match(names(mESC_RNAseq_lncRNA_features), names(mESC_RNAseq_lncRNA_features_full))] <- mESC_RNAseq_lncRNA_features



# expression cage lncRNA
# mESC_CAGE_lncRNA_features
# mESC_CAGE_lncRNA_features_full <- numeric(length = nrow(lncRNA_chosen_gt1k_uniqTiles))
# names(mESC_CAGE_lncRNA_features_full) <- lncRNA_chosen_gt1k_uniqTiles$gene_name
# mESC_CAGE_lncRNA_features_full[match(names(mESC_CAGE_lncRNA_features), names(mESC_CAGE_lncRNA_features_full))] <- mESC_CAGE_lncRNA_features
# 
# mESC_CHIPATLAS_lncRNA_features[is.na(mESC_CHIPATLAS_lncRNA_features)] <- 0
# mESC_RBP_lncRNA_features[is.na(mESC_RBP_lncRNA_features)] <- 0
# save(list = c("mESC_CAGE_lncRNA_features_full",
#               "mESC_RNAseq_lncRNA_features_full",
#               "mESC_RBP_lncRNA_features",
#               "kmer_features_lncRNA28_normalized",
#               "mESC_CHIPATLAS_lncRNA_features"),
#      file = "~/Documents/Shayan/BioInf/lncRNA/lncRNA_feaures.RData")


# Partition_4_dfs$tile_name <- levels(Partition_2_dfs$tile_name)[as.numeric(Partition_2_dfs$tile_name)]
# Partition_4_chunk_dfs$tile_name <- levels(Partition_2_chunk_dfs$tile_name)[as.numeric(Partition_2_chunk_dfs$tile_name)]
# 
Partition_6_feature_mat_tile <- cbind(ChIPATLAS_features_partition6, # remember to add repeats
                                      Tile_kmer_features_partition6,
                                      RBP_features_partition6,
                                      mESC_Renlab_RNASeq_partition6,
                                      mESC_CAGE_mm9_tiled_partition6)
colnames(Partition_6_feature_mat_tile)[ncol(Partition_6_feature_mat_tile) - 1] <- "RNAseq"
colnames(Partition_6_feature_mat_tile)[ncol(Partition_6_feature_mat_tile)] <- "CAGE"

# Partition_5_feature_mat_owner <- cbind(mESC_CHIPATLAS_lncRNA_features[match(Partition_5_dfs$owner, rownames(mESC_CHIPATLAS_lncRNA_features)),],
#                                        kmer_features_lncRNA28_normalized[match(Partition_5_dfs$owner, rownames(kmer_features_lncRNA28_normalized)),],
#                                        mESC_RBP_lncRNA_features[match(Partition_5_dfs$owner, rownames(mESC_RBP_lncRNA_features)),],
#                                        mESC_RNAseq_lncRNA_features_full[match(Partition_5_dfs$owner, names(mESC_RNAseq_lncRNA_features_full))],
#                                        mESC_CAGE_lncRNA_features_full[match(Partition_5_dfs$owner, names(mESC_CAGE_lncRNA_features_full))])
# colnames(Partition_5_feature_mat_owner)[ncol(Partition_5_feature_mat_owner) - 1] <- "RNAseq"
# colnames(Partition_5_feature_mat_owner)[ncol(Partition_5_feature_mat_owner)] <- "CAGE"


# colnames(Partition_5_feature_mat_owner) <- paste(colnames(Partition_5_feature_mat_owner), "lncRNA", sep = "__")

#load("~/Documents/Shayan/BioInf/lncRNA/Tile_owner_pair_RBP_ChIP_pair.RData")
Partition_6_feature_mat_pair <- cbind(#distance_feature_partition_6,
                                      Triplex_feaures_partition6[match(Partition_6_dfs$tile_name , Triplex_feaures_partition6$tile_name), 3],
                                      RBP_tile_owner_pair_feature[match(Partition_6_dfs$tile_name , names(RBP_tile_owner_pair_feature))], 
                                      ChIP_tile_owner_pair_feature[match(Partition_6_dfs$tile_name ,  names(ChIP_tile_owner_pair_feature))],
                                      Chromatin_tile_owner_pair_feature[match(Partition_6_dfs$tile_name ,  names(Chromatin_tile_owner_pair_feature))])
colnames(Partition_6_feature_mat_pair) <- c(#"distance",
                                            "triplex", "RBP_pair", "ChIP_pair", "Chromatin_pair")


aa_all_feat <- cbind(Partition_6_feature_mat_tile,
                     #Partition_5_feature_mat_owner,
                     Partition_6_feature_mat_pair,
                     Partition_6_dfs$label)
colnames(aa_all_feat)[ncol(aa_all_feat)] <- "label"
#Partition_6_dfs$owner <- levels(Partition_6_dfs$owner)[as.numeric(Partition_6_dfs$owner )]
aa_un_own <- unique(Partition_6_dfs$owner)
Zero_variance_columns_list <- list()
name_dic_list <- list()
Highly_correlating_columns_list <- list()
for(i in 1:length(aa_un_own)){
  print(i)
  aacurf  <- aa_all_feat[Partition_6_dfs$owner %in% aa_un_own[i],]
  aacurf[is.na(aacurf)] <- 0
  print("computing variance ...")
  aa_zerovar <- apply(aacurf, MARGIN = 2, FUN = var)
  Zero_variance_columns_list[[i]] <- colnames(aacurf)[aa_zerovar == 0]
  aacurf <- aacurf[, -which(aa_zerovar == 0)]
  
  aadescrCor <- cor(aacurf[, 1:(ncol(aacurf) - 1)])
  print("computing correlation ...")
  aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .90)
  Highly_correlating_columns_list[[i]] <- colnames(aacurf)[aahighlyCorDescr]
  aacurf <- aacurf[,-aahighlyCorDescr]
  aacurf <- as.data.frame(aacurf)
  aacurf$label[aacurf$label == 0] <- "neg"
  aacurf$label[aacurf$label == 1] <- "pos"
  aacurf$label <- factor(aacurf$label , levels = c("neg", "pos"))
  name_dic_list[[i]]<- cbind(colnames(aacurf), paste0("feature_", c(1:ncol(aacurf))))
  name_dic_list[[i]][ncol(aacurf), 2] <- name_dic_list[[i]][ncol(aacurf), 1] 
  colnames(aacurf) <-name_dic_list[[i]][, 2]
  my_Dataset <- aacurf
  my_partition_rand  <- Partition_6_dfs[Partition_6_dfs$owner %in% aa_un_own[i],]
  my_partition_chunk <- Partition_6_chunk_dfs[Partition_6_chunk_dfs$owner %in% aa_un_own[i],]
  my_name_dic <- name_dic_list[[i]]
  print("saving  ...")
  save(list = c("my_Dataset", "my_partition_rand", "my_partition_chunk", "my_name_dic"),
       file = paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_modified_dataset_",
                     aa_un_own[i], ".RData"))
}
save(list = c("Zero_variance_columns_list", "Highly_correlating_columns_list"),
     file = "Partition_6_removedColumns_namedic.RData")

# write jobs ro run RF per lncRNA
aa_un_own <- unique(Partition_6_dfs$owner)
for(i in 1:length(aa_un_own)){
  cat(c("Rscript --vanilla RF_load_run_p6.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_",
        aa_un_own[i],
        ".RData rand ", sample(c(100:10000), 1), "\n"),sep = "", file = "run_rf_p6.job", append = T)
  cat(c("Rscript --vanilla RF_load_run_p6.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_",
        aa_un_own[i],
        ".RData chunk ", sample(c(100:10000), 1), "\n"),sep = "", file = "run_rf_p6.job", append = T)
}




# create datasets with the new feature set
# read the previous ones and rewrite the to add the three sub triplex features
aafiles <- list.files("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets")
aaname <- unlist(lapply(strsplit(unlist(lapply(strsplit(aafiles, "\\."), "[[", 1)), "_"), "[[", 5))
for(i in 1:length(aafiles)){
  print(i)
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets/", aafiles[i]))
  aa_cur_lab <- my_Dataset[, ncol(my_Dataset)]
  aa_m_trp <- Triplex_feaures_partition6_extended_merge[Triplex_feaures_partition6_extended_merge$owner %in% aaname[i],]
  aa_new_feat <- cbind(my_Dataset[, 1:(ncol(my_Dataset) - 1)], aa_m_trp[match(rownames(my_Dataset), aa_m_trp$tile_name), c(6:8)])
  stopifnot(sum(is.na(aa_new_feat)) == 0)
  aa_wh_trp <- which(my_name_dic[, 1] == "triplex")
  aa_new_feat[, aa_wh_trp] <- aa_m_trp[match(rownames(my_Dataset), aa_m_trp$tile_name), 5]
  aanamedic <- my_name_dic[1:(ncol(my_Dataset) - 1),]
  aanamedic <- rbind(aanamedic, cbind(colnames(Triplex_feaures_partition6_extended_merge)[c(6:8)],
                                      paste0("feature_", c(ncol(my_Dataset):(ncol(my_Dataset) + 2)))))
  colnames(aa_new_feat) <- aanamedic[, 2]
  my_name_dic <- rbind(aanamedic, my_name_dic[nrow(my_name_dic),])
  
  my_Dataset <- cbind(aa_new_feat, aa_cur_lab)
  colnames(my_Dataset)[ncol(my_Dataset)] <- "label"
  save(list = c("my_Dataset", "my_name_dic", "my_partition_chunk", "my_partition_rand"), 
       file = paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/", aafiles[i]))
}


#############
#running MALAT1 without expression features: RNA-seq, CAGE, RNApol2


load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Malat1.RData")
View(my_name_dic)
# Write MALAT1 jobs without expression features
cat(c("Rscript --vanilla RF_load_run_p6_impurity.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Malat1.RData rand ",
      sample(c(100:10000), 1), " 229 1484 1485 ", "\n"),
    sep = "", file = "run_rf_p6_malat1_noExp.job", append = T)
cat(c("Rscript --vanilla RF_load_run_p6_impurity.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Malat1.RData chunk ",
      sample(c(100:10000), 1), " 229 1484 1485 ", "\n"),
    sep = "", file = "run_rf_p6_malat1_noExp.job", append = T)

#############
#running MALAT1 without expression and chromatin features



load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Malat1.RData")
View(my_name_dic)
aaexp <-c(229, 1484, 1485 )
aachrom <- grep(pattern = '^H[0-9]{1,2}', x = (my_name_dic[, 1]))
# Write MALAT1 jobs without expression features
cat(c("Rscript --vanilla RF_load_run_p6_impurity.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Malat1.RData rand ",
      sample(c(100:10000), 1), aaexp, aachrom, "\n"),
    sep = " ", file = "run_rf_p6_malat1_noExp_noChrom.job", append = T)
cat(c("Rscript --vanilla RF_load_run_p6_impurity.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Malat1.RData chunk ",
      sample(c(100:10000), 1), aaexp, aachrom, "\n"),
    sep = " ", file = "run_rf_p6_malat1_noExp_noChrom.job", append = T)
#############
#running MALAT1 without expression and chromatin features and no chrom pair feaure

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Malat1.RData")
View(my_name_dic)
aaexp <-c(229, 1484, 1485 )
aachrom <- grep(pattern = '^H[0-9]{1,2}', x = my_name_dic[, 1])
aachrom_pair <- grep(pattern = 'Chromatin_pair', x = my_name_dic[, 1])
# Write MALAT1 jobs without expression features
cat(c("Rscript --vanilla RF_load_run_p6_impurity.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Malat1.RData rand ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair, "\n"),
    sep = " ", file = "run_rf_p6_malat1_noExp_noChromPair.job", append = F)
cat(c("Rscript --vanilla RF_load_run_p6_impurity.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Malat1.RData chunk ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair, "\n"),
    sep = " ", file = "run_rf_p6_malat1_noExp_noChromPair.job", append = T)
#############
#running MALAT1 without expression and chromatin features and no chip and no DNase features

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Malat1.RData")
View(my_name_dic)
aaexp <-c(229, 1484, 1485 )
aachrom <- grep(pattern = '^H[0-9]{1,2}', x = my_name_dic[, 1])
aachrom_pair <- grep(pattern = 'Chromatin_pair', x = my_name_dic[, 1])
aabed <- grep(pattern = '*.bed', x =my_name_dic[, 1])
aabed <- setdiff(aabed, aachrom)
aa_chipPair <- grep(pattern = 'ChIP_pair', x = my_name_dic[, 1])
#aaDNase <- grep(pattern = 'DNase-Seq', x = my_name_dic[, 1])
# Write MALAT1 jobs without expression features
cat(c("Rscript --vanilla RF_load_run_p6_impurity.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Malat1.RData rand ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair ,"\n"),
    sep = " ", file = "run_rf_p6_malat1_noExp_noChrom_noBed.job", append = F)
cat(c("Rscript --vanilla RF_load_run_p6_impurity.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Malat1.RData chunk ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair ,"\n"),
    sep = " ", file = "run_rf_p6_malat1_noExp_noChrom_noBed.job", append = T)
#############
#running MALAT1 without expression and chromatin features and no chip features and no RBP features



load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Malat1.RData")
View(my_name_dic)
aaexp <-c(229, 1484, 1485 )
aachrom <- grep(pattern = '^H[0-9]{1,2}', x = my_name_dic[, 1])
aachrom_pair <- grep(pattern = 'Chromatin_pair', x = my_name_dic[, 1])
aabed <- grep(pattern = '*.bed', x =my_name_dic[, 1])
aabed <- setdiff(aabed, aachrom)
aa_chipPair <- grep(pattern = 'ChIP_pair', x = my_name_dic[, 1])
aarbp <- grep(pattern = '_mm9', x =my_name_dic[, 1])
aa_rbpPair <- grep(pattern = 'RBP_pair', x = my_name_dic[, 1])
# Write MALAT1 jobs without expression features
cat(c("Rscript --vanilla RF_load_run_p6_impurity.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Malat1.RData rand ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair,aarbp,aa_rbpPair, "\n"),
    sep = " ", file = "run_rf_p6_malat1_noExp_noChrom_noBed_noRBP.job", append = F)
cat(c("Rscript --vanilla RF_load_run_p6_impurity.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Malat1.RData chunk ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair,aarbp,aa_rbpPair, "\n"),
    sep = " ", file = "run_rf_p6_malat1_noExp_noChrom_noBed_noRBP.job", append = T)

#############
#running MALAT1 without expression and chromatin features and no chip features and no RBP features and no triplex features



load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Malat1.RData")
View(my_name_dic)
aaexp <-c(229, 1484, 1485 )
aachrom <- grep(pattern = '^H[0-9]{1,2}', x = my_name_dic[, 1])
aachrom_pair <- grep(pattern = 'Chromatin_pair', x = my_name_dic[, 1])
aabed <- grep(pattern = '*.bed', x =my_name_dic[, 1])
aabed <- setdiff(aabed, aachrom)
aa_chipPair <- grep(pattern = 'ChIP_pair', x = my_name_dic[, 1])
aarbp <- grep(pattern = '_mm9', x =my_name_dic[, 1])
aa_rbpPair <- grep(pattern = 'RBP_pair', x = my_name_dic[, 1])
aa_triplex <- grep(pattern = 'triplex*', x = my_name_dic[, 1])
# Write MALAT1 jobs without expression features
cat(c("Rscript --vanilla RF_load_run_p6_impurity.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Malat1.RData rand ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair,aarbp,aa_rbpPair,aa_triplex, "\n"),
    sep = " ", file = "run_rf_p6_malat1_noExp_noChrom_noBed_noRBP_noTriplex.job", append = F)
cat(c("Rscript --vanilla RF_load_run_p6_impurity.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Malat1.RData chunk ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair,aarbp,aa_rbpPair,aa_triplex, "\n"),
    sep = " ", file = "run_rf_p6_malat1_noExp_noChrom_noBed_noRBP_noTriplex.job", append = T)

#############
#running for GM14820, kcnq1ot1, Trerf1, Neat1 using 2 feature sets each: one where all cell context is removed, another with only kmers



load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Neat1.RData")
#View(my_name_dic)
aaexp_names <- c("RNA_polymerase_II__315.bed", "RNAseq" ,"CAGE")
aaexp <- which(my_name_dic[,1] %in% aaexp_names )
aachrom <- grep(pattern = '^H[0-9]{1,2}', x = my_name_dic[, 1])
aachrom_pair <- grep(pattern = 'Chromatin_pair', x = my_name_dic[, 1])
aabed <- grep(pattern = '*.bed', x =my_name_dic[, 1])
aabed <- setdiff(aabed, aachrom)
aa_chipPair <- grep(pattern = 'ChIP_pair', x = my_name_dic[, 1])
aarbp <- grep(pattern = '_mm9', x =my_name_dic[, 1])
aa_rbpPair <- grep(pattern = 'RBP_pair', x = my_name_dic[, 1])
aa_triplex <- grep(pattern = 'triplex*', x = my_name_dic[, 1])
# Write NEAT1 jobs without cell-context features
cat(c("Rscript --vanilla RF_load_run_p6_impurity_2.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Neat1.RData /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Learned_models_p6_NEAT1/ /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/performance_plots_p6_NEAT1/ rand ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair, "\n"),
    sep = " ", file = "run_rf_p6_Neat1_Gm14820_kcnq1ot1_Trerf1_featureExp.job", append = F)
cat(c("Rscript --vanilla RF_load_run_p6_impurity_2.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Neat1.RData /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Learned_models_p6_NEAT1/ /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/performance_plots_p6_NEAT1/ chunk ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair, "\n"),
    sep = " ", file = "run_rf_p6_Neat1_Gm14820_kcnq1ot1_Trerf1_featureExp.job", append = T)
# Write NEAT1 jobs with only kmers
cat(c("Rscript --vanilla RF_load_run_p6_impurity_2.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Neat1.RData /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Learned_models_p6_NEAT1/ /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/performance_plots_p6_NEAT1/ rand ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair,aarbp,aa_rbpPair,aa_triplex, "\n"),
    sep = " ", file = "run_rf_p6_Neat1_Gm14820_kcnq1ot1_Trerf1_featureExp.job", append = T)
cat(c("Rscript --vanilla RF_load_run_p6_impurity_2.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Neat1.RData /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Learned_models_p6_NEAT1/ /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/performance_plots_p6_NEAT1/ chunk ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair,aarbp,aa_rbpPair,aa_triplex, "\n"),
    sep = " ", file = "run_rf_p6_Neat1_Gm14820_kcnq1ot1_Trerf1_featureExp.job", append = T)
#######
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Gm14820.RData")
#View(my_name_dic)
aaexp_names <- c("RNA_polymerase_II__315.bed", "RNAseq" ,"CAGE")
aaexp <- which(my_name_dic[,1] %in% aaexp_names )
aachrom <- grep(pattern = '^H[0-9]{1,2}', x = my_name_dic[, 1])
aachrom_pair <- grep(pattern = 'Chromatin_pair', x = my_name_dic[, 1])
aabed <- grep(pattern = '*.bed', x =my_name_dic[, 1])
aabed <- setdiff(aabed, aachrom)
aa_chipPair <- grep(pattern = 'ChIP_pair', x = my_name_dic[, 1])
aarbp <- grep(pattern = '_mm9', x =my_name_dic[, 1])
aa_rbpPair <- grep(pattern = 'RBP_pair', x = my_name_dic[, 1])
aa_triplex <- grep(pattern = 'triplex*', x = my_name_dic[, 1])
# Write Gm14820 jobs without cell-context features
cat(c("Rscript --vanilla RF_load_run_p6_impurity_2.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Gm14820.RData /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Learned_models_p6_GM14820/ /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/performance_plots_p6_GM14820/ rand ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair, "\n"),
    sep = " ", file = "run_rf_p6_Neat1_Gm14820_kcnq1ot1_Trerf1_featureExp.job", append = T)
cat(c("Rscript --vanilla RF_load_run_p6_impurity_2.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Gm14820.RData /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Learned_models_p6_GM14820/ /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/performance_plots_p6_GM14820/ chunk ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair, "\n"),
    sep = " ", file = "run_rf_p6_Neat1_Gm14820_kcnq1ot1_Trerf1_featureExp.job", append = T)
# Write Gm14820 jobs with only kmers
cat(c("Rscript --vanilla RF_load_run_p6_impurity_2.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Gm14820.RData /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Learned_models_p6_GM14820/ /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/performance_plots_p6_GM14820/ rand ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair,aarbp,aa_rbpPair,aa_triplex, "\n"),
    sep = " ", file = "run_rf_p6_Neat1_Gm14820_kcnq1ot1_Trerf1_featureExp.job", append = T)
cat(c("Rscript --vanilla RF_load_run_p6_impurity_2.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Gm14820.RData /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Learned_models_p6_GM14820/ /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/performance_plots_p6_GM14820/ chunk ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair,aarbp,aa_rbpPair,aa_triplex, "\n"),
    sep = " ", file = "run_rf_p6_Neat1_Gm14820_kcnq1ot1_Trerf1_featureExp.job", append = T)

#######
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Kcnq1ot1.RData")
#View(my_name_dic)
aaexp_names <- c("RNA_polymerase_II__315.bed", "RNAseq" ,"CAGE")
aaexp <- which(my_name_dic[,1] %in% aaexp_names )
aachrom <- grep(pattern = '^H[0-9]{1,2}', x = my_name_dic[, 1])
aachrom_pair <- grep(pattern = 'Chromatin_pair', x = my_name_dic[, 1])
aabed <- grep(pattern = '*.bed', x =my_name_dic[, 1])
aabed <- setdiff(aabed, aachrom)
aa_chipPair <- grep(pattern = 'ChIP_pair', x = my_name_dic[, 1])
aarbp <- grep(pattern = '_mm9', x =my_name_dic[, 1])
aa_rbpPair <- grep(pattern = 'RBP_pair', x = my_name_dic[, 1])
aa_triplex <- grep(pattern = 'triplex*', x = my_name_dic[, 1])
# Write Kcnq1ot1 jobs without cell-context features
cat(c("Rscript --vanilla RF_load_run_p6_impurity_2.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Kcnq1ot1.RData /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Learned_models_p6_KCNQ1OT1/ /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/performance_plots_p6_KCNQ1OT1/ rand ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair, "\n"),
    sep = " ", file = "run_rf_p6_Neat1_Gm14820_kcnq1ot1_Trerf1_featureExp.job", append = T)
cat(c("Rscript --vanilla RF_load_run_p6_impurity_2.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Kcnq1ot1.RData /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Learned_models_p6_KCNQ1OT1/ /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/performance_plots_p6_KCNQ1OT1/ chunk ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair, "\n"),
    sep = " ", file = "run_rf_p6_Neat1_Gm14820_kcnq1ot1_Trerf1_featureExp.job", append = T)
# Write Kcnq1ot1 jobs with only kmers
cat(c("Rscript --vanilla RF_load_run_p6_impurity_2.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Kcnq1ot1.RData /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Learned_models_p6_KCNQ1OT1/ /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/performance_plots_p6_KCNQ1OT1/ rand ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair,aarbp,aa_rbpPair,aa_triplex, "\n"),
    sep = " ", file = "run_rf_p6_Neat1_Gm14820_kcnq1ot1_Trerf1_featureExp.job", append = T)
cat(c("Rscript --vanilla RF_load_run_p6_impurity_2.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Kcnq1ot1.RData /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Learned_models_p6_KCNQ1OT1/ /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/performance_plots_p6_KCNQ1OT1/ chunk ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair,aarbp,aa_rbpPair,aa_triplex, "\n"),
    sep = " ", file = "run_rf_p6_Neat1_Gm14820_kcnq1ot1_Trerf1_featureExp.job", append = T)

#######
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_Trerf1.RData")
#View(my_name_dic)
aaexp_names <- c("RNA_polymerase_II__315.bed", "RNAseq" ,"CAGE")
aaexp <- which(my_name_dic[,1] %in% aaexp_names )
aachrom <- grep(pattern = '^H[0-9]{1,2}', x = my_name_dic[, 1])
aachrom_pair <- grep(pattern = 'Chromatin_pair', x = my_name_dic[, 1])
aabed <- grep(pattern = '*.bed', x =my_name_dic[, 1])
aabed <- setdiff(aabed, aachrom)
aa_chipPair <- grep(pattern = 'ChIP_pair', x = my_name_dic[, 1])
aarbp <- grep(pattern = '_mm9', x =my_name_dic[, 1])
aa_rbpPair <- grep(pattern = 'RBP_pair', x = my_name_dic[, 1])
aa_triplex <- grep(pattern = 'triplex*', x = my_name_dic[, 1])
# Write Trerf1 jobs without cell-context features
cat(c("Rscript --vanilla RF_load_run_p6_impurity_2.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Trerf1.RData /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Learned_models_p6_TRERF1/ /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/performance_plots_p6_TRERF1/ rand ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair, "\n"),
    sep = " ", file = "run_rf_p6_Neat1_Gm14820_kcnq1ot1_Trerf1_featureExp.job", append = T)
cat(c("Rscript --vanilla RF_load_run_p6_impurity_2.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Trerf1.RData /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Learned_models_p6_TRERF1/ /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/performance_plots_p6_TRERF1/ chunk ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair, "\n"),
    sep = " ", file = "run_rf_p6_Neat1_Gm14820_kcnq1ot1_Trerf1_featureExp.job", append = T)
# Write Trerf1 jobs with only kmers
cat(c("Rscript --vanilla RF_load_run_p6_impurity_2.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Trerf1.RData /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Learned_models_p6_TRERF1/ /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/performance_plots_p6_TRERF1/ rand ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair,aarbp,aa_rbpPair,aa_triplex, "\n"),
    sep = " ", file = "run_rf_p6_Neat1_Gm14820_kcnq1ot1_Trerf1_featureExp.job", append = T)
cat(c("Rscript --vanilla RF_load_run_p6_impurity_2.R /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_modified_dataset_Trerf1.RData /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Learned_models_p6_TRERF1/ /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/performance_plots_p6_TRERF1/ chunk ",
      sample(c(100:10000), 1), aaexp, aachrom,aachrom_pair,aabed,aa_chipPair,aarbp,aa_rbpPair,aa_triplex, "\n"),
    sep = " ", file = "run_rf_p6_Neat1_Gm14820_kcnq1ot1_Trerf1_featureExp.job", append = T)


###################################################################################################################################################
###################################################################################################################################################
# Examining transcription signal at predicted binding site of RBPs
# for each lncRNA draw a violin plot with 224 violins. One corresponding to the expression (RNA-seq, CAGE, pol2) of all the examples related to that
#  lncRNA (partition6 pos and neg) and the remaining 223 correspond to the expression of tiles which carry at least one binding site for particular RBPs.

# also do this for all lncRNAs combined together, see how the signals persist.

# load features
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_RBP_scanned.RData")
load("~/Documents/Shayan/BioInf/lncRNA/ChIPATLAS_features_partition6.RData")
load("~/Documents/Shayan/BioInf/lncRNA/partition_6_expression_features.RData")

#mESC_Renlab_RNASeq_partition6
#mESC_CAGE_mm9_tiled_partition6

aaexp_names <- c("RNA_polymerase_II__315.bed", "RNAseq" ,"CAGE")
aafiles <- list.files("~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_dataset/")
aafilelnc <- unlist(lapply(strsplit(unlist(lapply(strsplit(aafiles, "\\."), "[[", 1)), "_"), "[[", 5))

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_U1snRNA.RData")
aa_lnc_feat_list <- list()
aa_lnc_nuRBPsites <- list()

aa_lnc_all_featues <- list()
#mESC_CAGE_mm9_tiled
#mESC_Renlab_tiled_list
for(i in 1:length(aafiles)){
  print(aafilelnc[i])

  
  if(i > 0){
    load(paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_dataset/", aafiles[i]))
    #aaexp <- match(aaexp_names, my_name_dic[,1])
    
    #aarbp <- grep(pattern = '_mm9', x =my_name_dic[, 1])
    aa_nu_tiles <- numeric(length = ncol(Partition_6_RBP_scanned))
    names(aa_nu_tiles) <- colnames(Partition_6_RBP_scanned)
    
    aa_alltiles <- rownames(my_Dataset)
    aamch <- match(aa_alltiles, rownames(Partition_6_RBP_scanned))
    aadata_p3 <- data.frame(tile_name = character(0),
                            label =  factor(levels = c("neg", "pos")),
                            RNA_polymerase_II__315.bed = numeric(0),
                            RNAseq= numeric(0),
                            CAGE= numeric(0))
    aa_tmp_RBP <- Partition_6_RBP_scanned[aamch,]
    aa_tmp_pol2 <- partition_6_expression_features[aamch, "RNA_polymerase_II__315.bed"]
    aa_tmp_pol2[is.na(aa_tmp_pol2)] <- 0
    #aaarcm <- match(aa_alltiles, names(mESC_Renlab_RNASeq_partition6))
    aa_tmp_RNAseq <-  partition_6_expression_features[aamch, "RNAseq"]
    aa_tmp_cage <-  partition_6_expression_features[aamch, "CAGE"]
    
    aadata_p4 <- data.frame(tile_name = aa_alltiles,
                            label =  my_Dataset$label,
                            RNA_polymerase_II__315.bed = aa_tmp_pol2,
                            RNAseq= aa_tmp_RNAseq,
                            CAGE= aa_tmp_cage)
  }




  
  aaaxc <- 1
  # if(i == 19){
  #   aaaxc <- 132
  # }
  for(j in aaaxc:ncol(Partition_6_RBP_scanned)){
    print(j)

    aa_cirind <- which(aa_tmp_RBP[, j] > 0)
    aa_haveTiles <- aa_alltiles[aa_cirind]
    aa_nu_tiles[j] <- length(aa_haveTiles)
    if(length(aa_haveTiles) > 0){
    
      aa_cur <- data.frame(tile_name = aa_haveTiles, 
                           label=my_Dataset$label[aa_cirind],
                           RNA_polymerase_II__315.bed=aa_tmp_pol2[aa_cirind],
                           RNAseq=aa_tmp_RNAseq[aa_cirind],
                           CAGE=aa_tmp_cage[aa_cirind], 
                           RBP = rep(colnames(Partition_6_RBP_scanned)[j], length(aa_haveTiles)))
      rownames(aa_cur) <- NULL
      aadata_p3 <- rbind(aadata_p3, aa_cur)
      
    }
  }
  aa_lnc_nuRBPsites[[i]] <- aa_nu_tiles
  aa_lnc_feat_list[[i]] <- aadata_p3
  aa_lnc_all_featues[[i]] <- aadata_p4
}
names(aa_lnc_nuRBPsites) <- aafilelnc
names(aa_lnc_feat_list) <- aafilelnc
names(aa_lnc_all_featues) <- aafilelnc

RBP_binding_site_expression <- aa_lnc_feat_list
save(list = c("RBP_binding_site_expression"), file = "~/Documents/Shayan/BioInf/lncRNA/RBP_binding_site_expression_UPDATAED.RData")

aa_folder_p <- "~/Documents/Shayan/BioInf/lncRNA/RBP_and_Transcription_UPDATE/"

load("~/Documents/Shayan/BioInf/lncRNA/RBP_binding_site_expression_UPDATAED.RData")
aa_lnc_feat_list <- RBP_binding_site_expression
aa_gplot_RNAseq_violin_list <- list()
aa_gplot_CAGE_violin_list <- list()
aa_gplot_POL2_violin_list <- list()
aa_gplot_RNAseq_box_list <- list()
aa_gplot_CAGE_box_list <- list()
aa_gplot_POL2_box_list <- list()

for(i in 1:length(aafilelnc)){
  print(i)
 # aa_lnc_feat_list[[i]]$RBP <- factor(aa_lnc_feat_list[[i]]$RBP, levels = colnames(RBP_features))
  # dawing for RNAseq
  # png(filename = paste0(aa_folder_p, aafilelnc[i], "_perRBP_RNAseq_violin.png"),    
  #     width = 40*300,       
  #     height = 5*300,
  #     res = 300,            
  #     pointsize = 10)
 # rownames(aa_lnc_feat_list[[i]]) <- NULL
  aa_gplot_RNAseq_violin_list[[i]] <- ggplot(aa_lnc_feat_list[[i]], aes(x=RBP, y=RNAseq+ 1, fill = label)) + 
    scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                  labels = trans_format("log10", math_format(10^.x))) +
    
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), scale = "area") +
    geom_vline(xintercept = seq(1.5,223.5,1), color ="blue", size = 0.75) +
    ggtitle(paste0("RNA-seq at RBP predicted sites for ", aafilelnc[i]))+
    ylab("RNAseq signal")+
    xlab("RBP") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 90))
  #dev.off()
  
  # png(filename = paste0(aa_folder_p, aafilelnc[i], "_perRBP_RNAseq_boxplot.png"),    
  #     width = 40*300,       
  #     height = 5*300,
  #     res = 300,            
  #     pointsize = 10)
 # rownames(aa_lnc_feat_list[[i]]) <- NULL
  aa_gplot_RNAseq_box_list[[i]] <- ggplot(aa_lnc_feat_list[[i]], aes(x=RBP, y=RNAseq+ 1, fill = label)) + 
    scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                  labels = trans_format("log10", math_format(10^.x))) +
    
    geom_boxplot() +
    geom_vline(xintercept = seq(1.5,223.5,1), color ="blue", size = 0.75) +
    ggtitle(paste0("RNA-seq at RBP predicted sites for ", aafilelnc[i]))+
    ylab("RNAseq signal")+
    xlab("RBP") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 90))
  #dev.off()
  
  # dawing for CAGE
  # png(filename = paste0(aa_folder_p, aafilelnc[i], "_perRBP_CAGE_violin.png"),    
  #     width = 40*300,        
  #     height = 5*300,
  #     res = 300,            
  #     pointsize = 10)
  aa_gplot_CAGE_violin_list[[i]] <- ggplot(aa_lnc_feat_list[[i]], aes(x=RBP, y=CAGE + 1, fill = label)) + 
    scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_vline(xintercept = seq(1.5,223.5,1), color ="blue", size = 0.75) +
    ggtitle(paste0("CAGE at RBP predicted sites for ", aafilelnc[i]))+
    ylab("CAGE signal")+
    xlab("RBP") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 90))
  #dev.off()
  
  # png(filename = paste0(aa_folder_p, aafilelnc[i], "_perRBP_CAGE_boxplot.png"),    
  #     width = 40*300,        
  #     height = 5*300,
  #     res = 300,            
  #     pointsize = 10)
  aa_gplot_CAGE_box_list[[i]] <- ggplot(aa_lnc_feat_list[[i]], aes(x=RBP, y=CAGE + 1, fill = label)) + 
    scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() +
    geom_vline(xintercept = seq(1.5,223.5,1), color ="blue", size = 0.75) +
    ggtitle(paste0("CAGE at RBP predicted sites for ", aafilelnc[i]))+
    ylab("CAGE signal")+
    xlab("RBP") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 90))
  #dev.off()
  # drawing for pol2
  # png(filename = paste0(aa_folder_p, aafilelnc[i], "_perRBP_pol2_violin.png"),   
  #     width = 40*300,       
  #     height = 5*300,
  #     res = 300,           
  #     pointsize = 10)
  
  aa_gplot_POL2_violin_list[[i]] <- ggplot(aa_lnc_feat_list[[i]], aes(x=RBP, y=RNA_polymerase_II__315.bed + 1, fill = label)) + 
    scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_vline(xintercept = seq(1.5,223.5,1), color ="blue", size = 0.75) +
    ggtitle(paste0("pol2 at RBP predicted sites for ", aafilelnc[i]))+
    ylab("pol2 signal")+
    xlab("RBP") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 90))
  #dev.off()
  
  # png(filename = paste0(aa_folder_p, aafilelnc[i], "_perRBP_pol2_boxplot.png"),   
  #     width = 40*300,       
  #     height = 5*300,
  #     res = 300,           
  #     pointsize = 10)
  
  aa_gplot_POL2_box_list[[i]] <- ggplot(aa_lnc_feat_list[[i]], aes(x=RBP, y=RNA_polymerase_II__315.bed + 1, fill = label)) + 
    scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() +
    geom_vline(xintercept = seq(1.5,223.5,1), color ="blue", size = 0.75) +
    ggtitle(paste0("pol2 at RBP predicted sites for ", aafilelnc[i]))+
    ylab("pol2 signal")+
    xlab("RBP") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 90))
  # dev.off()
  
}

for(i in 27:length(aafilelnc)){
  print(i)
  # png(filename = paste0(aa_folder_p, aafilelnc[i], "_perRBP_pol2_boxplot.png"),
  #     width = 40*300,
  #     height = 5*300,
  #     res = 300,
  #     pointsize = 10)
  # print(aa_gplot_POL2_box_list[[i]])
  # dev.off()
  
  png(filename = paste0(aa_folder_p, aafilelnc[i], "_perRBP_pol2_violin.png"),
      width = 40*300,
      height = 5*300,
      res = 300,
      pointsize = 10)
  print(aa_gplot_POL2_violin_list[[i]])
  dev.off()
  
  png(filename = paste0(aa_folder_p, aafilelnc[i], "_perRBP_RNAseq_boxplot.png"),
      width = 40*300,
      height = 5*300,
      res = 300,
      pointsize = 10)
  print(aa_gplot_RNAseq_box_list[[i]])
  dev.off()
  
  png(filename = paste0(aa_folder_p, aafilelnc[i], "_perRBP_RNAseq_violin.png"),
      width = 40*300,
      height = 5*300,
      res = 300,
      pointsize = 10)
  print(aa_gplot_RNAseq_violin_list[[i]])
  dev.off()
  
  png(filename = paste0(aa_folder_p, aafilelnc[i], "_perRBP_CAGE_boxplot.png"),
      width = 40*300,
      height = 5*300,
      res = 300,
      pointsize = 10)
  print(aa_gplot_CAGE_box_list[[i]])
  dev.off()
  
  png(filename = paste0(aa_folder_p, aafilelnc[i], "_perRBP_CAGE_violin.png"),
      width = 40*300,
      height = 5*300,
      res = 300,
      pointsize = 10)
  print(aa_gplot_CAGE_violin_list[[i]])
  dev.off()
}

png(filename = paste0(aa_folder_p, "nu_BSperRBP_perlncRNA.png"),
    width = 40*300,
    height = 60*300,
    res = 300,
    pointsize = 12)
par(mfrow = c(28, 1), mar = c(1,4,4,1))
for(i in 1:length(aafilelnc)){
  if(i == length(aafilelnc)){
    par( mar = c(12,4,4,1))
    barplot(aa_lnc_nuRBPsites[[i]], las = 2, main = aafilelnc[i])
  }else{
    barplot(aa_lnc_nuRBPsites[[i]], xaxt = "n", main = aafilelnc[i])
  }
  
}
dev.off()


aa_tot_Count <- as.matrix(table(partition_6_expression_features$lncRNA, partition_6_expression_features$label))

aa_pos_rbp_cnt <- matrix(nrow = 28, ncol = 223)
aa_neg_rbp_cnt <- matrix(nrow = 28, ncol = 223)
colnames(aa_pos_rbp_cnt) <- colnames(RBP_features)
colnames(aa_neg_rbp_cnt) <- colnames(RBP_features)
rownames(aa_pos_rbp_cnt) <- names(RBP_binding_site_expression)
rownames(aa_neg_rbp_cnt) <- names(RBP_binding_site_expression)

aa_pos_rbp_pct <- matrix(nrow = 28, ncol = 223)
aa_neg_rbp_pct <- matrix(nrow = 28, ncol = 223)
colnames(aa_pos_rbp_pct) <- colnames(RBP_features)
colnames(aa_neg_rbp_pct) <- colnames(RBP_features)
rownames(aa_pos_rbp_pct) <- names(RBP_binding_site_expression)
rownames(aa_neg_rbp_pct) <- names(RBP_binding_site_expression)

aa_all_rbp_cnt <- matrix(nrow = 28, ncol = 223)
aa_all_rbp_pct <- matrix(nrow = 28, ncol = 223)
colnames(aa_all_rbp_pct) <- colnames(RBP_features)
rownames(aa_all_rbp_pct) <- names(RBP_binding_site_expression)
colnames(aa_all_rbp_cnt) <- colnames(RBP_features)
rownames(aa_all_rbp_cnt) <- names(RBP_binding_site_expression)

for(i in 1:length(RBP_binding_site_expression)){
  print(i)
  for(j in 1:ncol(RBP_features)){
    aa_pos_rbp_cnt[i, j] <- sum(RBP_binding_site_expression[[i]][RBP_binding_site_expression[[i]]$RBP == colnames(RBP_features)[j], 2] == "pos")
    aa_neg_rbp_cnt[i, j] <- sum(RBP_binding_site_expression[[i]][RBP_binding_site_expression[[i]]$RBP == colnames(RBP_features)[j], 2] == "neg")
    aa_pos_rbp_pct[i, j] <- aa_pos_rbp_cnt[i, j]/aa_tot_Count[i, 2]
    aa_neg_rbp_pct[i, j] <- aa_neg_rbp_cnt[i, j]/aa_tot_Count[i, 1]
    aa_all_rbp_cnt[i, j] <- aa_pos_rbp_cnt[i, j] + aa_neg_rbp_cnt[i, j]
    aa_all_rbp_pct[i, j] <- (aa_all_rbp_cnt[i, j])/(sum(aa_tot_Count[i, ]))
  }
}


png(filename = paste0(aa_folder_p, "percent_BSperRBP_perlncRNA.png"),
    width = 40*300,
    height = 60*300,
    res = 300,
    pointsize = 12)
par(mfrow = c(28, 1), mar = c(2,4,4,1))
for(i in 1:length(aafilelnc)){
  if(i == length(aafilelnc)){
    par( mar = c(10,4,4,1))
    barplot(aa_all_rbp_pct[i,], las = 2, main = aafilelnc[i])
  }else{
    barplot(aa_all_rbp_pct[i,], xaxt = "n", main = aafilelnc[i])
  }
  
}
dev.off()

png(filename = paste0(aa_folder_p, "percent_BSperRBP_perlncRNA_perlabel.png"),
    width = 40*300,
    height = 60*300,
    res = 300,
    pointsize = 12)
aa_perlabel_rbp_pct <- matrix(nrow= nrow(aa_all_rbp_pct), ncol = ncol(aa_all_rbp_pct)*2)
colnames(aa_perlabel_rbp_pct) <- character(length = ncol(aa_perlabel_rbp_pct))
colnames(aa_perlabel_rbp_pct)[seq(1, 445, 2)] <- paste0(colnames(aa_all_rbp_pct), "_neg")
colnames(aa_perlabel_rbp_pct)[seq(2, 446, 2)] <- paste0(colnames(aa_all_rbp_pct), "_pos")
rownames(aa_perlabel_rbp_pct) <- rownames(aa_all_rbp_pct)
aa_perlabel_rbp_pct[, seq(1, 445, 2)] <- aa_neg_rbp_pct
aa_perlabel_rbp_pct[, seq(2, 446, 2)] <- aa_pos_rbp_pct


par(mfrow = c(28, 1), mar = c(2,4,4,1))
for(i in 1:length(aafilelnc)){
  if(i == length(aafilelnc)){
    par( mar = c(10,4,4,1))
    barplot(aa_perlabel_rbp_pct[i,], las = 2, main = aafilelnc[i], col = rep(c("red", "green"), 223))
    abline(v = seq(2.45, 536, 2.4), col = "blue")
  }else{
    barplot(aa_perlabel_rbp_pct[i,], xaxt = "n", main = aafilelnc[i], col = rep(c("red", "green"), 223))
    abline(v = seq(2.45, 536, 2.4), col = "blue")
  }
  
}
dev.off()


for(i in 1:length(aafilelnc)){
  aa_lnc_all_featues[[i]]$lncRNA <- rep(aafilelnc[i], nrow(aa_lnc_all_featues[[i]]))
}
aa_lnc_all_featues_all <- do.call(rbind, aa_lnc_all_featues)
aa_lnc_all_featues_all$lncRNA <- factor(aa_lnc_all_featues_all$lncRNA, levels = aafilelnc)

partition_6_expression_features <- aa_lnc_all_featues_all
aatst <- rownames(partition_6_expression_features)
aatst2 <- unlist(lapply(strsplit(aatst, "\\."), "[[", 2))
aamch <- match(Partition_6_dfs$tile_name, aatst2)
partition_6_expression_features <- partition_6_expression_features[aamch, ]
#rownames(partition_6_expression_features) <- aatst2[aamch]
save(list = c("partition_6_expression_features"), file = "~/Documents/Shayan/BioInf/lncRNA/partition_6_expression_features.RData")


png(filename = paste0(aa_folder_p, "RNAseq_by_lncRNA_by_label.png"),    # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10)

ggplot(aa_lnc_all_featues_all, aes(x=lncRNA, y=RNAseq + 1, fill = label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  #ggtitle("all points (2116029) distance to owner origin")+
  #ggtitle(paste0("pol2 at RBP predicted sites for ", aafilelnc[i]))+
  ylab("RNAseq signal")+
  xlab("lncRNA") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90))
dev.off()

png(filename = paste0(aa_folder_p, "CAGE_by_lncRNA_by_label.png"),    # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10)

ggplot(aa_lnc_all_featues_all, aes(x=lncRNA, y=CAGE+1, fill = label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  #ggtitle("all points (2116029) distance to owner origin")+
  #ggtitle(paste0("pol2 at RBP predicted sites for ", aafilelnc[i]))+
  ylab("CAGE signal")+
  xlab("lncRNA") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90))
dev.off()

png(filename = paste0(aa_folder_p, "pol2_by_lncRNA_by_label.png"),    # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10)

ggplot(aa_lnc_all_featues_all, aes(x=lncRNA, y=RNA_polymerase_II__315.bed+1, fill = label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  #ggtitle("all points (2116029) distance to owner origin")+
  #ggtitle(paste0("pol2 at RBP predicted sites for ", aafilelnc[i]))+
  ylab("pol2 signal")+
  xlab("lncRNA") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90))
dev.off()

png(filename = paste0(aa_folder_p, "RNAseq_by_lncRNA_by_label_boxplot.png"),    # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10)

ggplot(aa_lnc_all_featues_all, aes(x=lncRNA, y=RNAseq + 1, fill = label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot() +
  #ggtitle("all points (2116029) distance to owner origin")+
  #ggtitle(paste0("pol2 at RBP predicted sites for ", aafilelnc[i]))+
  ylab("RNAseq signal")+
  xlab("lncRNA") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90))
dev.off()

png(filename = paste0(aa_folder_p, "CAGE_by_lncRNA_by_label_boxplot.png"),    # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10)

ggplot(aa_lnc_all_featues_all, aes(x=lncRNA, y=CAGE+1, fill = label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot() +
  #ggtitle("all points (2116029) distance to owner origin")+
  #ggtitle(paste0("pol2 at RBP predicted sites for ", aafilelnc[i]))+
  ylab("CAGE signal")+
  xlab("lncRNA") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90))
dev.off()

png(filename = paste0(aa_folder_p, "pol2_by_lncRNA_by_label_boxplot.png"),    # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10)

ggplot(aa_lnc_all_featues_all, aes(x=lncRNA, y=RNA_polymerase_II__315.bed+1, fill = label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot() +
  #ggtitle("all points (2116029) distance to owner origin")+
  #ggtitle(paste0("pol2 at RBP predicted sites for ", aafilelnc[i]))+
  ylab("pol2 signal")+
  xlab("lncRNA") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90))
dev.off()

#####
# counting the number of RBPs with predicted sites on each tile


aa_nu_rbp_perlnc_list <- list()
for(i in 1:length(aafilelnc)){
  print(i)
  load(paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_dataset/", aafiles[i]))
  aamtch <- match(rownames(my_Dataset), rownames(Partition_6_RBP_scanned))
  aa_nu_rbp_perlnc_list[[i]] <-  data.frame(tile_name = rownames(my_Dataset), nu_BS = rowSums(Partition_6_RBP_scanned[aamtch,] > 0, na.rm = T), label = my_Dataset$label)
}
names(aa_nu_rbp_perlnc_list) <- aafilelnc


png(filename = paste0(aa_folder_p, "nu_RBP_with_BS_perlnc_pie.png"),    # create PNG for the heat map        
    width = 6*300,        # 5 x 300 pixels
    height = 56*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10)

par(mfrow = c(28,3), mar = c(3,3,3,3))
for(i in 1:length(aa_nu_rbp_perlnc_list)){
  print(i)
#  aa_perc <- quantile(x = aa_nu_rbp_perlnc_list[[i]]$nu_BS, seq(0.01,1, 0.01))
  pie(table(aa_nu_rbp_perlnc_list[[i]]$nu_BS), main = paste0(aafilelnc[i], "_all"))
  pie(table(aa_nu_rbp_perlnc_list[[i]]$nu_BS[aa_nu_rbp_perlnc_list[[i]]$label == "neg"]), main = paste0(aafilelnc[i], "_neg"))
  pie(table(aa_nu_rbp_perlnc_list[[i]]$nu_BS[aa_nu_rbp_perlnc_list[[i]]$label == "pos"]), main = paste0(aafilelnc[i], "_pos"))
}
dev.off()


png(filename = paste0(aa_folder_p, "nu_RBP_with_BS_perlnc_comulative.png"),    # create PNG for the heat map        
    width = 10*300,        # 5 x 300 pixels
    height = 20*300,
    res = 300,            # 300 pixels per inch
    pointsize = 10)

par(mfrow = c(7,4), mar = c(3,3,3,3))
for(i in 1:length(aa_nu_rbp_perlnc_list)){
  print(i)
  aa_perc <- quantile(x = aa_nu_rbp_perlnc_list[[i]]$nu_BS, seq(0.01,1, 0.01))
  aa_perc1 <- quantile(x = aa_nu_rbp_perlnc_list[[i]]$nu_BS[aa_nu_rbp_perlnc_list[[i]]$label == "neg"],
                       seq(0.01,1, 0.01))
  aa_perc2 <- quantile(x = aa_nu_rbp_perlnc_list[[i]]$nu_BS[aa_nu_rbp_perlnc_list[[i]]$label == "pos"],
                       seq(0.01,1, 0.01))
  plot(aa_perc, type = "l", ylab = "# uniq RBP site", main = aafilelnc[i])
  lines(aa_perc1, col = 2)
  lines(aa_perc2, col = 3)
}
dev.off()


#partition_6_expression_features

aa_nu_rbp_perlnc_all <- do.call(rbind, aa_nu_rbp_perlnc_list)
load("partition_6_expression_features.RData")
all(aa_nu_rbp_perlnc_all$tile_name == partition_6_expression_features$tile_name)

aa_df <- cbind(aa_nu_rbp_perlnc_all[match(partition_6_expression_features$tile_name,aa_nu_rbp_perlnc_all$tile_name),], partition_6_expression_features$RNA_polymerase_II__315.bed,partition_6_expression_features$RNAseq, partition_6_expression_features$CAGE)
names(aa_df) <- c(names(aa_nu_rbp_perlnc_all), "Pol2", "RNAseq", "CAGE")
aa_df_names <- unlist(lapply(strsplit(rownames(aa_df), "\\."), "[[", 1))
aa_df$lncRNA <- aa_df_names
aa_df$lncRNA <-factor(aa_df$lncRNA )

plot(aa_df$nu_BS, log10(aa_df$RNAseq+ 1))

aa_df$RBPYN <- aa_df$nu_BS > 0

boxplot(log10(aa_df$RNAseq+ 1) ~ aa_df$RBPYN)

png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/RNASeq_byPredictedRBP_byLabel2.png",    # create PNG for the heat map        
    width = 4*300,        # 5 x 300 pixels
    height = 4*300,
    res = 300,            # 300 pixels per inch
    pointsize = 12) 

ggplot(aa_df, aes(x=RBPYN, y=RNAseq+1, fill = label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  #ggtitle("all points (2116029) distance to owner origin")+
  ggtitle("RNA_Seq at RBP vs noRBP")+
  ylab("count")+
  #xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank())

dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/Pol2_byPredictedRBP_byLabel2.png",    # create PNG for the heat map        
    width = 4*300,        # 5 x 300 pixels
    height = 4*300,
    res = 300,            # 300 pixels per inch
    pointsize = 12) 

ggplot(aa_df, aes(x=RBPYN, y=Pol2+1, fill = label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  #ggtitle("all points (2116029) distance to owner origin")+
  ggtitle("RNA_Seq at RBP vs noRBP")+
  ylab("count")+
  #xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank())

dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/CAGE_byPredictedRBP_byLabel2.png",    # create PNG for the heat map        
    width = 4*300,        # 5 x 300 pixels
    height = 4*300,
    res = 300,            # 300 pixels per inch
    pointsize = 12) 

ggplot(aa_df, aes(x=RBPYN, y=CAGE+1, fill = label)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  #ggtitle("all points (2116029) distance to owner origin")+
  ggtitle("RNA_Seq at RBP vs noRBP")+
  ylab("count")+
  #xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank())

dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/RNASeq_byPredictedRBP_byLabel_bylncRNA2.png",    # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 12) 

ggplot(aa_df, aes(x=lncRNA, y=RNAseq+1, fill = RBPYN)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  #ggtitle("all points (2116029) distance to owner origin")+
  ggtitle("RNA_Seq at RBP vs noRBP")+
  ylab("count")+
  #xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90))

dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/Pol2_byPredictedRBP_byLabel_bylncRNA2.png",    # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 12) 

ggplot(aa_df, aes(x=lncRNA, y=Pol2+1, fill = RBPYN)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  #ggtitle("all points (2116029) distance to owner origin")+
  ggtitle("pol2 at RBP vs noRBP")+
  ylab("count")+
  #xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90))

dev.off()

png(filename = "~/Documents/Shayan/BioInf/lncRNA/plots/CAGE_byPredictedRBP_byLabel2.png",    # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 12) 

ggplot(aa_df, aes(x=lncRNA, y=CAGE+1, fill = RBPYN)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  #ggtitle("all points (2116029) distance to owner origin")+
  ggtitle("CAGE at RBP vs noRBP")+
  ylab("count")+
  #xlab("partition") + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90))

dev.off()
################################################################################################################################################################################################################
#########################################################################################################################################################################
# create new feature matrices for partition6

head(Partition_6_dfs)
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_dfs.RData")
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_chunk_random.RData")
#load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_RBP_scanned.RData")
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_repeat_features.RData")
load("~/Documents/Shayan/BioInf/lncRNA/partition_6_expression_features.RData")
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_U1snRNA.RData")
load("~/Documents/Shayan/BioInf/lncRNA/Triplex_feaures_partition6_extended.RData")
load("~/Documents/Shayan/BioInf/lncRNA/Tile_kmer_features_partition6.RData")
load("~/Documents/Shayan/BioInf/lncRNA/ChIPATLAS_features_partition6.RData")

ChIPATLAS_features_partition6 <- ChIPATLAS_features[match(Partition_6_dfs$tile_name, rownames(ChIPATLAS_features)),]
ChIPATLAS_features_partition6 <- ChIPATLAS_features_partition6[,-c(grep("RNA_polymerase_II__315.bed", colnames(ChIPATLAS_features_partition6)))]
save(list=c("ChIPATLAS_features_partition6"), file = "~/Documents/Shayan/BioInf/lncRNA/ChIPATLAS_features_partition6.RData")

Partition_6_TF_scanned <- read.table("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/All_TF_profile_P6.txt", header = T)
Partition_6_RBP_scanned <- read.table("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/All_RBP_profile_P6.txt", header = T)

colnames(Tile_kmer_features_partition6) <- paste0(colnames(Tile_kmer_features_partition6), "__5mer")
colnames(Partition_6_repeat_features) <-  paste0(colnames(Partition_6_repeat_features), "__rep")
colnames(Partition_6_TF_scanned) <-  paste0(colnames(Partition_6_TF_scanned), "__TFscan")
colnames(Partition_6_RBP_scanned) <-  paste0(colnames(Partition_6_RBP_scanned), "__RBPscan")
#colnames(ChIPATLAS_features_partition6) <-  paste0(colnames(ChIPATLAS_features_partition6), "__Chip")
#colnames(partition_6_expression_features) <-  paste0(colnames(partition_6_expression_features), "__RBPscan")
Partition_6_U1snRNA_mat <- matrix(nrow = length(Partition_6_U1snRNA), ncol = 1)
colnames(Partition_6_U1snRNA_mat) <- "U1snRNA_site__RBPscan"
Partition_6_U1snRNA_mat[, 1] <-  Partition_6_U1snRNA


Partition_6_feature_mat_tile_2 <- cbind(Tile_kmer_features_partition6,
                                        Partition_6_repeat_features, 
                                        Partition_6_TF_scanned, 
                                        Partition_6_RBP_scanned,
                                        Partition_6_U1snRNA_mat,
                                        ChIPATLAS_features_partition6,
                                        partition_6_expression_features[,(3:5)])

colnames(Partition_6_feature_mat_tile_2)[which(duplicated(colnames(Partition_6_feature_mat_tile_2)))] <- paste0(colnames(Partition_6_feature_mat_tile_2)[which(duplicated(colnames(Partition_6_feature_mat_tile_2)))], "_v2")
# make the pair features: RBP_pair, TF_pair, repeat_pair, chip_pair, chromatin_pair
aa_lncNames <- read.table("~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/mES/MM9_lncRNA_28.bed", stringsAsFactors = F)$V7

aa_lncRNA_RBP <- read.table("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/RBP_pwm_lncRNA_all.txt", header = T)
rownames(aa_lncRNA_RBP) <- aa_lncNames
aa_RBP_scanPair <- matrix(nrow = nrow(Partition_6_U1snRNA_mat), ncol = 1)
colnames(aa_RBP_scanPair) <- "RBPscan_pair"
aaPartition_6_RBP_scanned <- Partition_6_RBP_scanned[, match(colnames(aa_lncRNA_RBP), colnames(Partition_6_RBP_scanned))]
for(i in 1:nrow(aa_lncRNA_RBP)){
  print(i)
  aa_tmp <- add_pair_features(lncRNA_feat = aa_lncRNA_RBP[i, ], 
                              tile_feat = aaPartition_6_RBP_scanned[Partition_6_dfs$owner %in% aa_lncNames[i],],
                              my_network_list = numeric(0), presence_thr=0)
  aa_RBP_scanPair[Partition_6_dfs$owner %in% aa_lncNames[i], 1] <- aa_tmp
}
Partition_6_RBP_scanPair <- aa_RBP_scanPair
save(list = c("Partition_6_RBP_scanPair"), file = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_RBP_scanPair.RData")


aa_lncRNA_TF <- read.table("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/TF_pwm_lncRNA_all.txt", header = T)
rownames(aa_lncRNA_TF) <- aa_lncNames
aa_TF_scanPair <- matrix(nrow = nrow(aa_RBP_scanPair), ncol = 1)
colnames(aa_TF_scanPair) <- "TFscan_pair"
aacnam <- unlist(lapply(strsplit(colnames(Partition_6_TF_scanned), split = "__TFscan"), "[[", 1))
aaPartition_6_TF_scanned <- Partition_6_TF_scanned[, match(colnames(aa_lncRNA_TF), aacnam)]
for(i in 1:nrow(aa_lncRNA_TF)){
  print(i)
  aa_tmp <- add_pair_features(lncRNA_feat = aa_lncRNA_TF[i, ], 
                              tile_feat = aaPartition_6_TF_scanned[Partition_6_dfs$owner %in% aa_lncNames[i],],
                              my_network_list = numeric(0), presence_thr=0)
  aa_TF_scanPair[Partition_6_dfs$owner %in% aa_lncNames[i], 1] <- aa_tmp
}
Partition_6_TF_scanPair <- aa_TF_scanPair
save(list = c("Partition_6_TF_scanPair"), file = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_TF_scanPair.RData")

aa_lncRNA_repeat <- read.table("~/Documents/Shayan/BioInf/lncRNA/mESC_repeat_lncRNA_features.txt", header = T)
rownames(aa_lncRNA_repeat) <- aa_lncNames
aa_repeat_Pair <- matrix(nrow =  nrow(Partition_6_U1snRNA_mat), ncol = 1)
colnames(aa_repeat_Pair) <- "Repeat_pair"
aacnam <- unlist(lapply(strsplit(colnames(Partition_6_repeat_features), split = "__rep"), "[[", 1))
aaPartition_6_repeat <- Partition_6_repeat_features[, match(colnames(aa_lncRNA_repeat), aacnam)]
for(i in 1:nrow(aa_lncRNA_repeat)){
  print(i)
  aa_tmp <- add_pair_features(lncRNA_feat = aa_lncRNA_repeat[i, ], 
                              tile_feat = aaPartition_6_repeat[Partition_6_dfs$owner %in% aa_lncNames[i],],
                              my_network_list = numeric(0), presence_thr=0)
  aa_repeat_Pair[Partition_6_dfs$owner %in% aa_lncNames[i], 1] <- aa_tmp
}
Partition_6_Repeat_pair <- aa_repeat_Pair
save(list = c("Partition_6_Repeat_pair"), file = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Repeat_pair.RData")


load("~/Documents/Shayan/BioInf/lncRNA/Tile_owner_pair_RBP_ChIP_pair.RData")


Triplex_feaures_partition6_extended_merge
# colnames(Partition_5_feature_mat_owner) <- paste(colnames(Partition_5_feature_mat_owner), "lncRNA", sep = "__")

Partition_6_feature_mat_pair_2 <- cbind(Partition_6_RBP_scanPair, 
                                        Partition_6_TF_scanPair,
                                        Partition_6_Repeat_pair,
                                        Triplex_feaures_partition6_extended_merge[match(Partition_6_dfs$tile_name , 
                                                                                        Triplex_feaures_partition6_extended_merge$tile_name), c(5:8)],
                                        ChIP_tile_owner_pair_feature[match(Partition_6_dfs$tile_name ,  names(ChIP_tile_owner_pair_feature))],
                                        Chromatin_tile_owner_pair_feature[match(Partition_6_dfs$tile_name ,  names(Chromatin_tile_owner_pair_feature))])
colnames(Partition_6_feature_mat_pair_2) <- c("RBP_pair","TFscan_pair", "repeat_pair", "triplex_total","triplex_GA","triplex_TC","triplex_GT","ChIP_pair", "Chromatin_pair")


# aa_all_feat <- cbind(Partition_6_feature_mat_tile_2,
#                      Partition_6_feature_mat_pair_2,
#                      Partition_6_dfs$label)
# colnames(aa_all_feat)[ncol(aa_all_feat)] <- "label"
# #Partition_6_dfs$owner <- levels(Partition_6_dfs$owner)[as.numeric(Partition_6_dfs$owner )]
aa_un_own <- unique(Partition_6_dfs$owner)
#Zero_variance_columns_list <- list()
name_dic_list <- list()
#Highly_correlating_columns_list <- list()
for(i in 1:length(aa_un_own)){
  print(i)
  aacurf  <- cbind(Partition_6_feature_mat_tile_2[Partition_6_dfs$owner %in% aa_un_own[i],],
                   Partition_6_feature_mat_pair_2[Partition_6_dfs$owner %in% aa_un_own[i],],
                   Partition_6_dfs$label[Partition_6_dfs$owner %in% aa_un_own[i]])
  colnames(aacurf)[ncol(aacurf)] <- "label"
  aacurf[is.na(aacurf)] <- 0
#  print("computing variance ...")
#  aa_zerovar <- apply(aacurf, MARGIN = 2, FUN = var)
#  Zero_variance_columns_list[[i]] <- colnames(aacurf)[aa_zerovar == 0]
#  aacurf <- aacurf[, -which(aa_zerovar == 0)]
  
  # aadescrCor <- cor(aacurf[, 1:(ncol(aacurf) - 1)])
  # print("computing correlation ...")
  # aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .90)
  # Highly_correlating_columns_list[[i]] <- colnames(aacurf)[aahighlyCorDescr]
  # aacurf <- aacurf[,-aahighlyCorDescr]
  # aacurf <- as.data.frame(aacurf)
  aacurf$label[aacurf$label == 0] <- "neg"
  aacurf$label[aacurf$label == 1] <- "pos"
  aacurf$label <- factor(aacurf$label , levels = c("neg", "pos"))
  name_dic_list[[i]]<- cbind(colnames(aacurf), paste0("feature_", c(1:ncol(aacurf))))
  name_dic_list[[i]][ncol(aacurf), 2] <- name_dic_list[[i]][ncol(aacurf), 1] 
  colnames(aacurf) <-name_dic_list[[i]][, 2]
  my_Dataset <- aacurf
  my_partition_rand  <- Partition_6_dfs[Partition_6_dfs$owner %in% aa_un_own[i],]
  my_partition_chunk <- Partition_6_chunk_dfs[Partition_6_chunk_dfs$owner %in% aa_un_own[i],]
  my_name_dic <- name_dic_list[[i]]
  print("saving  ...")
  save(list = c("my_Dataset", "my_partition_rand", "my_partition_chunk", "my_name_dic"),
       file = paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_dataset/Partition_6_Dec_dataset_",
                     aa_un_own[i], ".RData"))
}
# save(list = c("Zero_variance_columns_list"#, 
#               #"Highly_correlating_columns_list"
#               ),
#      file = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_removedColumns_namedic.RData")



aa_filter_list <- list()
aa_feature_remove_list <- list()
aa_colnn <- c(colnames(Partition_6_feature_mat_tile_2), colnames(Partition_6_feature_mat_pair_2))
# 1) RBP + RBP_pair  "__RBPscan" "RBP_pair" 
aafeat1 <- grep(pattern = "__RBPscan", x = aa_colnn)
aafeat2 <- grep(pattern = "RBP_pair", x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[1]] <- aa_remove
aa_filter_list[[1]] <- numeric(0)
# 2) RBP + RBP_pair  "__RBPscan" "RBP_pair"  expression filtered
aafeat1 <- grep(pattern = "__RBPscan", x = aa_colnn)
aafeat2 <- grep(pattern = "RBP_pair", x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[2]] <- aa_remove
aa_filter_list[[2]] <- c(aafeat1, aafeat2)
# 3) TF + TF_pair  "__TFscan" "TFscan_pair" 
aafeat1 <- grep(pattern = "__TFscan", x = aa_colnn)
aafeat2 <- grep(pattern = "TFscan_pair", x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[3]] <- aa_remove
aa_filter_list[[3]] <- numeric(0)
# 4) TF + TF_pair  "__TFscan" "TFscan_pair"   expression filtered
aafeat1 <- grep(pattern = "__TFscan", x = aa_colnn)
aafeat2 <- grep(pattern = "TFscan_pair", x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[4]] <- aa_remove
aa_filter_list[[4]] <- c(aafeat1, aafeat2)
# 5) repeat + repeat_pair  "__rep" "repeat_pair" 
aafeat1 <- grep(pattern = "__rep", x = aa_colnn)
aafeat2 <- grep(pattern = "repeat_pair", x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[5]] <- aa_remove
aa_filter_list[[5]] <- numeric(0)
# 6) repeat + repeat_pair  "__rep" "repeat_pair"  expression filtered
aafeat1 <- grep(pattern = "__rep", x = aa_colnn)
aafeat2 <- grep(pattern = "repeat_pair", x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[6]] <- aa_remove
aa_filter_list[[6]] <- c(aafeat1, aafeat2)
# 7) RBP + RBP_pair + TF + TF_pair  repeat + repeat_pair "__rep" "repeat_pair" 
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), grep(pattern = "__TFscan", x = aa_colnn), grep(pattern = "__rep", x = aa_colnn))
aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn), grep(pattern = "TFscan_pair", x = aa_colnn), grep(pattern = "repeat_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[7]] <- aa_remove
aa_filter_list[[7]] <- numeric(0)
# 8) RBP + RBP_pair + TF + TF_pair  repeat + repeat_pair "__rep" "repeat_pair" expression filtered
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), grep(pattern = "__TFscan", x = aa_colnn), grep(pattern = "__rep", x = aa_colnn))
aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn), grep(pattern = "TFscan_pair", x = aa_colnn), grep(pattern = "repeat_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[8]] <- aa_remove
aa_filter_list[[8]] <- c(grep(pattern = "__RBPscan", x = aa_colnn),
                         grep(pattern = "__rep", x = aa_colnn),
                         grep(pattern = "repeat_pair", x = aa_colnn),
                         grep(pattern = "RBP_pair", x = aa_colnn))
# 9) kmer_only
aafeat1 <- c(grep(pattern = "__5mer", x = aa_colnn))
#aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn), grep(pattern = "TFscan_pair", x = aa_colnn), grep(pattern = "repeat_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), aafeat1)
aa_feature_remove_list[[9]] <- aa_remove
aa_filter_list[[9]] <- numeric(0)

# 10) RBP + RBP_pair + TF + TF_pair  repeat + repeat_pair "__rep" "repeat_pair"  + triplex
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), grep(pattern = "__TFscan", x = aa_colnn), grep(pattern = "__rep", x = aa_colnn))
aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn), grep(pattern = "TFscan_pair", x = aa_colnn), grep(pattern = "repeat_pair", x = aa_colnn), grep(pattern = "triplex", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[10]] <- aa_remove
aa_filter_list[[10]] <- numeric(0)

# 11) RBP + RBP_pair + TF + TF_pair  repeat + repeat_pair "__rep" "repeat_pair" expression filtered + triplex
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), grep(pattern = "__TFscan", x = aa_colnn), grep(pattern = "__rep", x = aa_colnn))
aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn), grep(pattern = "TFscan_pair", x = aa_colnn), grep(pattern = "repeat_pair", x = aa_colnn), grep(pattern = "triplex", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[11]] <- aa_remove
aa_filter_list[[11]] <- c(grep(pattern = "__RBPscan", x = aa_colnn),
                         grep(pattern = "__rep", x = aa_colnn),
                         grep(pattern = "repeat_pair", x = aa_colnn),
                         grep(pattern = "RBP_pair", x = aa_colnn))
# 12) chromatin + chromtin_pair
aafeat1 <- grep(pattern = '^H[0-9]{1,2}', x = aa_colnn)
aafeat2 <- grep(pattern = "Chromatin_pair", x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[12]] <- aa_remove
aa_filter_list[[12]] <- numeric(0)

# 13) Chip + ChIP_pair
aachrom <- grep(pattern = '^H[0-9]{1,2}', x = aa_colnn)
aabed <- grep(pattern = '*.bed', x = aa_colnn)
aabed <- setdiff(aabed, aachrom)
aafeat1 <- aabed
aafeat2 <- grep(pattern = "ChIP_pair", x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[13]] <- aa_remove
aa_filter_list[[13]] <- numeric(0)
# 14) expression
aaexp_names <- c("RNA_polymerase_II__315.bed", "RNAseq" ,"CAGE")
aafeat1 <- which(aa_colnn %in%  aaexp_names)
aa_remove <- setdiff(c(1:length(aa_colnn)), aafeat1)
aa_feature_remove_list[[14]] <- aa_remove
aa_filter_list[[14]] <- numeric(0)
# 15) chromatin + chromtin_pair + Chip + ChIP_pair + expression
aafeat1 <- c(grep(pattern = '*.bed', x = aa_colnn),  which(aa_colnn %in%  aaexp_names))
aafeat2 <- c(grep(pattern = "ChIP_pair", x = aa_colnn), grep(pattern = "Chromatin_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[15]] <- aa_remove
aa_filter_list[[15]] <- numeric(0)
# 16) (10) + (15)
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), grep(pattern = "__TFscan", x = aa_colnn), grep(pattern = "__rep", x = aa_colnn),
             grep(pattern = '*.bed', x = aa_colnn),which(aa_colnn %in%  aaexp_names))
aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn), grep(pattern = "TFscan_pair", x = aa_colnn),
             grep(pattern = "repeat_pair", x = aa_colnn), grep(pattern = "triplex", x = aa_colnn),
             grep(pattern = "ChIP_pair", x = aa_colnn), grep(pattern = "Chromatin_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[16]] <- aa_remove
aa_filter_list[[16]] <- numeric(0)
# 17) (11) + (15)
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), grep(pattern = "__TFscan", x = aa_colnn), grep(pattern = "__rep", x = aa_colnn),
             grep(pattern = '*.bed', x = aa_colnn),which(aa_colnn %in%  aaexp_names))
aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn), grep(pattern = "TFscan_pair", x = aa_colnn),
             grep(pattern = "repeat_pair", x = aa_colnn), grep(pattern = "triplex", x = aa_colnn),
             grep(pattern = "ChIP_pair", x = aa_colnn), grep(pattern = "Chromatin_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[17]] <- aa_remove
aa_filter_list[[17]] <-  c(grep(pattern = "__RBPscan", x = aa_colnn),
                           grep(pattern = "__rep", x = aa_colnn),
                           grep(pattern = "repeat_pair", x = aa_colnn),
                           grep(pattern = "RBP_pair", x = aa_colnn))
# 18) (16) + kmer
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), grep(pattern = "__TFscan", x = aa_colnn), grep(pattern = "__rep", x = aa_colnn),
             grep(pattern = '*.bed', x = aa_colnn),which(aa_colnn %in%  aaexp_names), grep(pattern = "__5mer", x = aa_colnn))
aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn), grep(pattern = "TFscan_pair", x = aa_colnn),
             grep(pattern = "repeat_pair", x = aa_colnn), grep(pattern = "triplex", x = aa_colnn),
             grep(pattern = "ChIP_pair", x = aa_colnn), grep(pattern = "Chromatin_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[18]] <- aa_remove
aa_filter_list[[18]] <- numeric(0)

# 19) (17) + kmer
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), grep(pattern = "__TFscan", x = aa_colnn), grep(pattern = "__rep", x = aa_colnn),
             grep(pattern = '*.bed', x = aa_colnn),which(aa_colnn %in%  aaexp_names), grep(pattern = "__5mer", x = aa_colnn))
aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn), grep(pattern = "TFscan_pair", x = aa_colnn),
             grep(pattern = "repeat_pair", x = aa_colnn), grep(pattern = "triplex", x = aa_colnn),
             grep(pattern = "ChIP_pair", x = aa_colnn), grep(pattern = "Chromatin_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[19]] <- aa_remove
aa_filter_list[[19]] <-  c(grep(pattern = "__RBPscan", x = aa_colnn),
                           grep(pattern = "__rep", x = aa_colnn),
                           grep(pattern = "repeat_pair", x = aa_colnn),
                           grep(pattern = "RBP_pair", x = aa_colnn))
names(aa_filter_list) <- c("RBPscanned", "RBPscannedTX", "TFscanned", "TFscannedTX", "Repeat", "RepeatTX", 
                           "RBPTFRepeat", "RBPTFRepeatTXRBPrep", "Kmer", "RBPTFRepeatTriplex", "RBPTFRepeatTriplexTXRBPrep", "Chromatin",
                           "ChIP", "Transcription", "ChromChIPexp", "RBPTFRepeatChromChIPexp", "RBPTFRepeatTXRBPrepChromChIPexp", "RBPTFRepeatChromChIPexpKmer",
                           "RBPTFRepeatTXRBPrepChromChIPexpKmer")

# write filter_list files + jobs
aa_un_own <- unique(Partition_6_dfs$owner)
for(j in 1:length(aa_feature_remove_list)){
  if(length(aa_filter_list[[j]]) == 0){
    aamyfilt <- " no_filter"
  }else{
    write.table(aa_filter_list[[j]], quote = F, row.names = F, col.names = F, 
                file = paste0("~/Documents/Shayan/BioInf/lncRNA/partition_6_Dec_TOBE_exp_FILTERED/my_filter_", j,".txt"), append = F)
    aamyfilt <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_Dec_TOBE_exp_FILTERED/my_filter_", j,".txt")
  }
  for(i in 1:length(aa_un_own)){
    cat(c("Rscript --vanilla RF_load_run_p6_Dec_feature.R",
          paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_Dec_dataset_",
                 aa_un_own[i],".RData "), 
          paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_Dec_Results/Learned_models/",aa_un_own[i],"/"),
          paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_Dec_Results/performance_plots/",aa_un_own[i],"/"),
          " rand ",
          sample(c(100:10000), 1),
          names(aa_filter_list)[j], 
          aamyfilt,
          aa_feature_remove_list[[j]],
          "\n"),
        sep = " ", file = "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full.job", append = T)
    
    cat(c("Rscript --vanilla RF_load_run_p6_Dec_feature.R",
          paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_Dec_dataset_",
                 aa_un_own[i],".RData "), 
          paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_Dec_Results/Learned_models/",aa_un_own[i],"/"),
          paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_Dec_Results/performance_plots/",aa_un_own[i],"/"),
          " chunk ",
          sample(c(100:10000), 1),
          names(aa_filter_list)[j],
          aamyfilt,
          aa_feature_remove_list[[j]],
          "\n"),
        sep = " ", file = "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full.job", append = T)
  }
  
}

##########
# Write new feature list jobs

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_dataset/Partition_6_Dec_dataset_Firre.RData")
aa_filter_list2 <- list()
aa_feature_remove_list2 <- list()
aa_colnn <- my_name_dic[1:(nrow(my_name_dic) - 1),1]  #c(colnames(Partition_6_feature_mat_tile_2), colnames(Partition_6_feature_mat_pair_2))



# 1) kmer_and chromatin
aafeat1 <- c(grep(pattern = "__5mer", x = aa_colnn))
aafeat2 <- c(grep(pattern = '^H[0-9]{1,2}', x = aa_colnn),
             grep(pattern = "Chromatin_pair", x = aa_colnn))
#aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn), grep(pattern = "TFscan_pair", x = aa_colnn), grep(pattern = "repeat_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list2[[1]] <- aa_remove
aa_filter_list2[[1]] <- numeric(0)


# 2) kmer_and chromatin and chip
aafeat1 <- c(grep(pattern = "__5mer", x = aa_colnn))

aafeat2 <- c(grep(pattern = '*.bed', x = aa_colnn),
             grep(pattern = 'ChIP_pair', x = aa_colnn),
             grep(pattern = "Chromatin_pair", x = aa_colnn))
aa_pol2 <- grep("RNA_polymerase_II__315.bed", x = aa_colnn)
#aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn), grep(pattern = "TFscan_pair", x = aa_colnn), grep(pattern = "repeat_pair", x = aa_colnn))
aa_remove <- c(setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2)), aa_pol2)
aa_feature_remove_list2[[2]] <- aa_remove
aa_filter_list2[[2]] <- numeric(0)


# 3) kmer_and chromatin and chip and exp
aaexp_names <- c("RNA_polymerase_II__315.bed", "RNAseq" ,"CAGE")

aafeat1 <- c(grep(pattern = "__5mer", x = aa_colnn), which(aa_colnn %in%  aaexp_names))

aafeat2 <- c(grep(pattern = '*.bed', x = aa_colnn),
             grep(pattern = 'ChIP_pair', x = aa_colnn),
             grep(pattern = "Chromatin_pair", x = aa_colnn))
#aa_pol2 <- grep("RNA_polymerase_II__315.bed", x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list2[[3]] <- aa_remove
aa_filter_list2[[3]] <- numeric(0)

# 4) kmer_and triplex
aafeat1 <- c(grep(pattern = "__5mer", x = aa_colnn))
aafeat2 <- c(grep(pattern = 'triplex', x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list2[[4]] <- aa_remove
aa_filter_list2[[4]] <- numeric(0)

# 5) kmer_and triplex and pairs
aafeat1 <- c(grep(pattern = "__5mer", x = aa_colnn))
aafeat2 <- c(grep(pattern = 'triplex', x = aa_colnn),
             grep("_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list2[[5]] <- aa_remove
aa_filter_list2[[5]] <- numeric(0)

#6) triplex only
aafeat2 <- c(grep(pattern = 'triplex', x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)),  aafeat2)
aa_feature_remove_list2[[6]] <- aa_remove
aa_filter_list2[[6]] <- numeric(0)

#7) pairs only
aafeat2 <- c(grep("_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)),  aafeat2)
aa_feature_remove_list2[[7]] <- aa_remove
aa_filter_list2[[7]] <- numeric(0)


# 8) triplex and pairs
#aafeat1 <- c(grep(pattern = "__5mer", x = aa_colnn))
aafeat2 <- c(grep(pattern = 'triplex', x = aa_colnn),
             grep("_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)),  aafeat2)
aa_feature_remove_list2[[8]] <- aa_remove
aa_filter_list2[[8]] <- numeric(0)



names(aa_feature_remove_list2) <- c("KmerChrom", "KmerChromChIP", "KmerChromChipExp", "KmerTriplex", 
                                    "KmerTriplexPairs", "Triplex", "Pairs", "TriplexPairs")
names(aa_filter_list2) <- c("KmerChrom", "KmerChromChIP", "KmerChromChipExp", "KmerTriplex", 
                            "KmerTriplexPairs", "Triplex", "Pairs", "TriplexPairs")
# write filter_list files + jobs
aa_un_own <- unique(Partition_6_dfs$owner)
for(j in 1:length(aa_feature_remove_list2)){
  if(length(aa_filter_list2[[j]]) == 0){
    aamyfilt <- " no_filter"
  }else{
    write.table(aa_filter_list2[[j]], quote = F, row.names = F, col.names = F, 
                file = paste0("~/Documents/Shayan/BioInf/lncRNA/partition_6_Dec_TOBE_exp_FILTERED2/my_filter_", j,".txt"), append = F)
    aamyfilt <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_Dec_TOBE_exp_FILTERED2/my_filter_", j,".txt")
  }
  for(i in 1:length(aa_un_own)){
    cat(c("Rscript --vanilla RF_load_run_p6_Dec_feature.R",
          paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_Dec_dataset_",
                 aa_un_own[i],".RData "), 
          paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_Dec_Results/Learned_models/",aa_un_own[i],"/"),
          paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_Dec_Results/performance_plots/",aa_un_own[i],"/"),
          " rand ",
          sample(c(100:10000), 1),
          names(aa_filter_list2)[j], 
          aamyfilt,
          aa_feature_remove_list2[[j]],
          "\n"),
        sep = " ", file = "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full2.job", append = T)
    
    cat(c("Rscript --vanilla RF_load_run_p6_Dec_feature.R",
          paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_Dec_dataset_",
                 aa_un_own[i],".RData "), 
          paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_Dec_Results/Learned_models/",aa_un_own[i],"/"),
          paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_Dec_Results/performance_plots/",aa_un_own[i],"/"),
          " chunk ",
          sample(c(100:10000), 1),
          names(aa_filter_list2)[j],
          aamyfilt,
          aa_feature_remove_list2[[j]],
          "\n"),
        sep = " ", file = "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full2.job", append = T)
  }
  
}



######## run all expression filtered jobs using the new expression features (GROseq, PROseq)
aa_job1 <- readLines("~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full.job")
aa_job1Fil <- aa_job1[grep("filter_", aa_job1)]
aa_job1Fil2 <- gsub(pattern = "TX", replacement = "TXGRO", x = aa_job1Fil)
writeLines(aa_job1Fil2, "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full_TXGRO.job")

######## run randomized expression filtered jobs 
aa_job1 <- readLines("~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full.job")
aa_job1Fil <- aa_job1[grep("filter_", aa_job1)]
aa_job1Fil2 <- gsub(pattern = "TX", replacement = "TXRand", x = aa_job1Fil)
aa_job1Fil2 <- gsub(pattern = "RF_load_run_p6_Dec_feature.R", replacement = "RF_load_run_p6_Dec_feature_shuffle.R", x = aa_job1Fil2)

writeLines(aa_job1Fil2, "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full_TXrandomized.job")

######## run only using GRO-PRO data
aa_job1 <- readLines("~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full.job")
aa_job1Fil <- aa_job1[grep("RBPTFRepeatChromChIPexpKmer", aa_job1)]
aa_job1Fil2 <- gsub(pattern = "RBPTFRepeatChromChIPexpKmer", replacement = "GROPRO", x = aa_job1Fil)
aa_job1Fil2 <- gsub(pattern = "RF_load_run_p6_Dec_feature.R", replacement = "RF_load_run_p6_Dec_feature_GROPRO.R", x = aa_job1Fil2)

aa_job1Fil2[1]
writeLines(aa_job1Fil2, "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full_ONLYGROPRO.job")


######## Add  GRO-PRO to some models
aa_job1 <- readLines("~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full.job")
aa_job1Fil1 <- aa_job1[grep("\\bRBPTFRepeat\\b", aa_job1)]
aa_job1Fil2 <- aa_job1[grep("\\bRBPTFRepeatChromChIPexp\\b", aa_job1)]
aa_job1Fil3 <- aa_job1[grep("\\bKmerChromChipExp\\b", aa_job1)]
aa_job1Fil4 <- aa_job1[grep("\\bKmer\\b", aa_job1)]
aa_job1Fil5 <- aa_job1[grep("\\bChromChIPexp\\b", aa_job1)]

aa_job1Fil1 <- gsub(pattern = "RF_load_run_p6_Dec_feature.R", replacement = "RF_load_run_p6_Dec_feature_ADDGROPRO.R", x = aa_job1Fil1)
aa_job1Fil1 <-  gsub(pattern = "RBPTFRepeat", replacement = "RBPTFRepeatAddGROPRO", x = aa_job1Fil1)

aa_job1Fil2 <- gsub(pattern = "RF_load_run_p6_Dec_feature.R", replacement = "RF_load_run_p6_Dec_feature_ADDGROPRO.R", x = aa_job1Fil2)
aa_job1Fil2 <-  gsub(pattern = "RBPTFRepeatChromChIPexp", replacement = "RBPTFRepeatChromChIPexpAddGROPRO", x = aa_job1Fil2)

aa_job1Fil3 <- gsub(pattern = "RF_load_run_p6_Dec_feature.R", replacement = "RF_load_run_p6_Dec_feature_ADDGROPRO.R", x = aa_job1Fil3)
aa_job1Fil3 <-  gsub(pattern = "KmerChromChipExp", replacement = "KmerChromChipExpAddGROPRO", x = aa_job1Fil3)

aa_job1Fil4 <- gsub(pattern = "RF_load_run_p6_Dec_feature.R", replacement = "RF_load_run_p6_Dec_feature_ADDGROPRO.R", x = aa_job1Fil4)
aa_job1Fil4 <-  gsub(pattern = "Kmer", replacement = "KmerAddGROPRO", x = aa_job1Fil4)

aa_job1Fil5 <- gsub(pattern = "RF_load_run_p6_Dec_feature.R", replacement = "RF_load_run_p6_Dec_feature_ADDGROPRO.R", x = aa_job1Fil5)
aa_job1Fil5 <-  gsub(pattern = "ChromChIPexp", replacement = "ChromChIPexpAddGROPRO", x = aa_job1Fil5)


aa_job1Fil_AddGROPRO <- c(aa_job1Fil1, aa_job1Fil2, aa_job1Fil3, aa_job1Fil4, aa_job1Fil5)

writeLines(aa_job1Fil_AddGROPRO, "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full_AddedGROPRO.job")


############################################################################################################
############################################################################################################
# add partitions to both random and chunk settings, such that you make a cross validation system

# aadata_p6 <- cbind(Partition_6_dfs, 
#                    Partition_6_chunk_dfs$dataset,
#                    distance_feature_partition_6)
par(mar = c(10,4,4,4))
barplot(t(table(aadata_p6$owner, aadata_p6$label)[2:28,]), las = 2, legend.text = c("neg", "pos"))

# colnames(aadata_p6) <- c("tile_name", "owner",     "label" ,  "dataset_random", "dataset_chunk","distance" )
# aadata_p6$label[aadata_p6$label == 1] <- "pos"
# aadata_p6$label[aadata_p6$label == 0] <- "neg"
# aadata_p6$label <- factor(aadata_p6$label, levels = c("neg", "pos"))

# create 5 random CV partitions
#Partition_6_dfs$owner <- levels(Partition_6_dfs$owner)[as.numeric(Partition_6_dfs$owner)]
aa_lnc <- unique(Partition_6_dfs$owner)
Partition_6_random_cv <- list()


length(Partition_6_dfs_GR)
set.seed(87634)
for(cur_lnc in 1:length(aa_lnc)){
  print(cur_lnc)
  aaw <- which(Partition_6_dfs$owner %in% aa_lnc[cur_lnc])
  # check overlap with lncRNA's own gene
  aa_ovl <- findOverlaps(query = Partition_6_dfs_GR[aaw],
                         subject = lncRNA_chosen_gt1k_uniqTiles_GR[which(lncRNA_chosen_gt1k_uniqTiles_GR$gene_name == aa_lnc[cur_lnc])])
  if(length(aa_ovl@from) > 0){
    aaw <- aaw[-unique(aa_ovl@from)]
  }
  
  
  #aa_tmp_tile <- Partition_6_dfs$tile_name[]
  
  Partition_6_random_cv[[cur_lnc]] <- data.frame(tile_name = Partition_6_dfs$tile_name[aaw],
                                                 label = Partition_6_dfs$label[aaw],
                                                 owner = Partition_6_dfs$owner[aaw],
                                                 CV1 = numeric(length(aaw)),
                                                 CV2 = numeric(length(aaw)),
                                                 CV3 = numeric(length(aaw)),
                                                 CV4 = numeric(length(aaw)),
                                                 CV5 = numeric(length(aaw)))
  
  aa_allpos <- which(Partition_6_random_cv[[cur_lnc]]$label == 1)
  aa_allneg <- which(Partition_6_random_cv[[cur_lnc]]$label == 0)
  aa_pos_shuff <- sample(aa_allpos, size = length(aa_allpos), replace = F)
  aa_neg_shuff <- sample(aa_allneg, size = length(aa_allneg), replace = F)
  
  aa_pos_shuff_part <- cut(c(1:length(aa_pos_shuff)), breaks = 5, labels = paste0("CV", c(1:5)))
  aa_neg_shuff_part <- cut(c(1:length(aa_neg_shuff)), breaks = 5, labels = paste0("CV", c(1:5)))
  
  for(cur_cv in c(1:5)){
    aa_test_pos <- aa_pos_shuff[which(aa_pos_shuff_part %in% paste0("CV", cur_cv))]
    aa_train_pos <- setdiff(aa_allpos, aa_test_pos)
    aa_test_neg <- aa_neg_shuff[which(aa_neg_shuff_part %in% paste0("CV", cur_cv))]
    aa_train_neg <- setdiff(aa_allneg, aa_test_neg)
    Partition_6_random_cv[[cur_lnc]][c(aa_test_pos, aa_test_neg),(3+cur_cv)] <- 0
    Partition_6_random_cv[[cur_lnc]][c(aa_train_pos, aa_train_neg),(3+cur_cv)] <- 1
  }
}
names(Partition_6_random_cv) <- aa_lnc

for(i in 1:length(aa_lnc)){
  print(table(Partition_6_random_cv[[i]]$label, Partition_6_random_cv[[i]]$CV1))
}

Partition_6_chunk_cv <- list()
aa_lnc <- unique(Partition_6_dfs$owner)
#aa_nu_start <- 40
aa_test_frac <- 0.2
aa_Valid_frac <- 0.2
set.seed(seed = 3246)
for(cur_lnc in 12:length(aa_lnc)){
  print(cur_lnc)
  print(aa_lnc[cur_lnc])
  aaw <- which(Partition_6_dfs$owner %in% aa_lnc[cur_lnc])
  # check overlap with lncRNA's own gene
  aa_ovl <- findOverlaps(query = Partition_6_dfs_GR[aaw],
                         subject = lncRNA_chosen_gt1k_uniqTiles_GR[which(lncRNA_chosen_gt1k_uniqTiles_GR$gene_name == aa_lnc[cur_lnc])])
  if(length(aa_ovl@from) > 0){
    aaw <- aaw[-unique(aa_ovl@from)]
  }
  
  Partition_6_chunk_cv[[cur_lnc]] <- data.frame(tile_name = Partition_6_dfs$tile_name[aaw],
                                                 label = Partition_6_dfs$label[aaw],
                                                 owner = Partition_6_dfs$owner[aaw],
                                                 CV1 = numeric(length(aaw)),
                                                 CV2 = numeric(length(aaw)),
                                                 CV3 = numeric(length(aaw)),
                                                 CV4 = numeric(length(aaw)),
                                                 CV5 = numeric(length(aaw)))
  aa_nu_pos <- sum(Partition_6_chunk_cv[[cur_lnc]]$label == 1)
  aa_nu_test_pos <- floor(aa_nu_pos*aa_test_frac)
  aa_nu_vali_pos <- floor(aa_nu_pos*aa_Valid_frac)
  aa_nu_train_pos <- aa_nu_pos - (aa_nu_test_pos + aa_nu_vali_pos)
  
  aa_nu_neg <- sum(Partition_6_chunk_cv[[cur_lnc]]$label == 0)
  aa_nu_test_neg <- floor(aa_nu_neg*aa_test_frac)
  aa_nu_vali_neg <- floor(aa_nu_neg*aa_Valid_frac)
  aa_nu_train_neg <- aa_nu_neg - (aa_nu_test_neg + aa_nu_vali_neg)
  aanegname <- Partition_6_chunk_cv[[cur_lnc]]$tile_name[Partition_6_chunk_cv[[cur_lnc]]$label == 0]
  aaposname <- Partition_6_chunk_cv[[cur_lnc]]$tile_name[Partition_6_chunk_cv[[cur_lnc]]$label == 1]
  
  aa_postile_num <- as.numeric(unlist(lapply(strsplit(aaposname, "_"), "[[", 2)))
  aa_negtile_num <- as.numeric(unlist(lapply(strsplit(aanegname, "_"), "[[", 2)))
  aa_postile_chr <- (unlist(lapply(strsplit(aaposname, "_"), "[[", 1)))
  aa_negtile_chr <- (unlist(lapply(strsplit(aanegname, "_"), "[[", 1)))
  
  aa_neg_pos_df <- data.frame(neg_index = aanegname, 
                              closest_pos = character(length = aa_nu_neg))
  aaunlen <- length(unique(aa_negtile_chr))
  
  if(aaunlen > 1){
    
    for(i in 1:aa_nu_neg){
      
      aawposchr <- which(aa_postile_chr %in% aa_negtile_chr[i])
      aatmppos_name <- aaposname[aawposchr]
      aatmp_postile_num <- aa_postile_num[aawposchr]
      aa_neg_pos_df[i, 2] <- aatmppos_name[which.min(abs(aa_negtile_num[i] - aatmp_postile_num))]
    }
  }else{
    for(i in 1:aa_nu_neg){
    aa_neg_pos_df[i, 2] <- aaposname[which.min(abs(aa_negtile_num[i] - aa_postile_num))]
    }
  }

 
  # divide the positives into chunks based on index
  aa_pos_chunk <- list()
  aa_nu_start <- ceiling(sqrt(aa_nu_pos))
  aa_poschunk_size <- floor(sqrt(aa_nu_pos))
  aapos_assignment <- cbind(aaposname,
                            cut(c(1:aa_nu_pos),
                                labels = F,
                                breaks = ceiling(sqrt(aa_nu_pos))))
  
  aaneg_assignment <- cbind(aanegname, 
                            aapos_assignment[match(aa_neg_pos_df$closest_pos,
                                                   aapos_assignment[, 1]), 2])
  print("done pos neg assign")
  
  # choosing parts for each group
  aa_ne_prop_error <- 0.5
  aa_rep <- 1
  aa_ne_prop_error_thr <- 0.02
  aathr_step <- 0.04
  while(aa_ne_prop_error > aa_ne_prop_error_thr){
    # print("started cv")
    aashuffle <- sample(c(1:aa_nu_start), size = aa_nu_start, replace = F)
    # print("aashuffle")
    # print(aashuffle)
    aatest_parts_cut <- cut(x = c(1:aa_nu_start), breaks = 5, labels = paste0("CV", c(1:5)))
    # aatest_parts <- sample(x = c(1:aa_nu_start), size = floor(aa_nu_start * aa_test_frac))
    # 
    # aavalid_parts <- sample(x =setdiff(c(1:aa_nu_start), aatest_parts), size = floor(aa_nu_start * aa_Valid_frac))
    # aatrain_parts <- setdiff(setdiff(c(1:aa_nu_start), aatest_parts), aavalid_parts)
    # forming the dataframes
    aa_cv_check <- numeric(5)
    for(cur_cv in 1:5){
      aa_test_chunk <- aashuffle[which(aatest_parts_cut %in% paste0("CV", cur_cv))]
      aa_test_pos <- aapos_assignment[aapos_assignment[, 2] %in% aa_test_chunk, 1]
      aa_test_neg <- aaneg_assignment[aaneg_assignment[, 2] %in% aa_test_chunk, 1]
      aa_train_all <- setdiff(c(aanegname, aaposname), c(aa_test_pos, aa_test_neg))
      
      Partition_6_chunk_cv[[cur_lnc]][match(c(aa_test_pos, aa_test_neg), Partition_6_chunk_cv[[cur_lnc]]$tile_name),(3+cur_cv)] <- 0
      Partition_6_chunk_cv[[cur_lnc]][match(aa_train_all, Partition_6_chunk_cv[[cur_lnc]]$tile_name),(3+cur_cv)] <- 1
      
      aa_neg_tab <- table(Partition_6_chunk_cv[[cur_lnc]][match(aanegname, Partition_6_chunk_cv[[cur_lnc]]$tile_name), (3+cur_cv)])/aa_nu_neg
      # print("aa_neg_tab")
      # print(aa_neg_tab)
      if(length(aa_neg_tab) == 2){
        aa_cv_check[cur_cv] <- sqrt(mean((aa_neg_tab - c((aa_test_frac), (1- aa_test_frac)))^2))
      }else{
        aa_cv_check[cur_cv] <- 10
      }
    }
    # print("aa_cv_check")
    # print(aa_cv_check)
    aa_ne_prop_error <- mean(aa_cv_check)
    # print("aa_ne_prop_error")
    # print(aa_ne_prop_error)
    # print("aa_ne_prop_error_thr")
    # print(aa_ne_prop_error_thr)
    aa_rep <- aa_rep + 1
    print(aa_rep)
    if(aa_rep == 300){
      print(paste0("changing thresh to ", aa_ne_prop_error_thr + aathr_step))
      aa_rep <- 1
      aa_ne_prop_error_thr <- aa_ne_prop_error_thr + aathr_step
    }
  }

}

names(Partition_6_chunk_cv) <- aa_lnc

for(i in 1:length(aa_lnc)){
  print(table(Partition_6_chunk_cv[[i]]$label, Partition_6_chunk_cv[[i]]$CV5))
}
save(list = c("Partition_6_random_cv", "Partition_6_chunk_cv"), file = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_random_chunk_cv_updated.RData")
Partition_6_random_cv_df <- do.call(rbind, Partition_6_random_cv)
Partition_6_chunk_cv_df <- do.call(rbind, Partition_6_chunk_cv)
Partition_6_chunk_cv_df$owner <- levels(Partition_6_chunk_cv_df$owner)[as.numeric(Partition_6_chunk_cv_df$owner)]
colnames(Partition_6_random_cv_df)[4:8] <- paste0("RCV", c(1:5))
colnames(Partition_6_chunk_cv_df)[4:8] <- paste0("CCV", c(1:5))

Partition_6_random_chunk_cv_df <- cbind(Partition_6_random_cv_df, Partition_6_chunk_cv_df[, c(4:8)])
#Partition_6_random_chunk_cv_df <- Partition_6_random_chunk_cv_df[match(Partition_6_dfs$tile_name, Partition_6_random_chunk_cv_df$tile_name),]

save(list = c("Partition_6_random_chunk_cv_df"), file = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_random_chunk_cv_df_updated.RData")



#########################################################################################################
#########################################################################################################
Partition_6_random_cv_df
Partition_6_chunk_cv_df
sum(Partition_6_chunk_cv_df$tile_name == Partition_6_dfs$tile_name)
sum(Partition_6_chunk_cv_df$owner == Partition_6_dfs$owner)
sum(Partition_6_chunk_cv_df$label == Partition_6_dfs$label)
sum(Partition_6_random_chunk_cv_df$tile_name == Partition_6_dfs$tile_name)
sum(Partition_6_random_chunk_cv_df$owner == Partition_6_dfs$owner)
sum(Partition_6_random_chunk_cv_df$label == Partition_6_dfs$label)


################################################################################################################################################################################################################
#########################################################################################################################################################################
# constructing feature sets to use in the cross-validation setting (partition6 cross validation) . improvements: 

# 
# add methylation 
# add dinucleotide_freq
# use kmer_ReverseCompl instead of the full kmer mat (512 kmer features instead of 1024)
# add the network-driven pair features (hetero)
# add distance 

## DNase-Seq__46.bed RNA_polymerase_II__315.bed
# mESC_ATAC_Seq_tiled


head(Partition_6_dfs)
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_random_chunk_cv_df_updated.RData")
#load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_RBP_scanned.RData")
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_repeat_features.RData")
load("~/Documents/Shayan/BioInf/lncRNA/partition_6_expression_features.RData")
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_U1snRNA.RData")
load("~/Documents/Shayan/BioInf/lncRNA/Triplex_feaures_partition6_extended.RData")
load("~/Documents/Shayan/BioInf/lncRNA/Tile_kmer_features_partition6_revComp.RData")
load("~/Documents/Shayan/BioInf/lncRNA/ChIPATLAS_features_partition6.RData")

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_accessiblity.RData")
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_methylation.RData")
load("~/Documents/Shayan/BioInf/lncRNA/distance_to_owner_partition6.RData")
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_dinuc_freq.RData")

sum(is.na(distance_feature_partition_6))
max(distance_feature_partition_6, na.rm = T) # 106540001
distance_feature_partition_6[is.na(distance_feature_partition_6)] <- max(distance_feature_partition_6, na.rm = T) + 1000
save(list= c("distance_feature_partition_6"), file = "~/Documents/Shayan/BioInf/lncRNA/distance_feature_partition_6.RData")

head(Partition_6_random_chunk_cv_df)
sum( Partition_6_random_chunk_cv_df2$tile_name == Partition_6_dfs$tile_name)
aadup1 <- Partition_6_dfs$tile_name[which(duplicated(Partition_6_dfs$tile_name))]
aadup2 <- Partition_6_random_chunk_cv_df$tile_name[which(duplicated(Partition_6_random_chunk_cv_df$tile_name))]

table(Partition_6_random_chunk_cv_df$owner[Partition_6_random_chunk_cv_df$tile_name %in% aadup1])
table(Partition_6_dfs$owner[Partition_6_dfs$tile_name %in% aadup1])

Partition_6_random_chunk_cv_df2[Partition_6_random_chunk_cv_df2$tile_name == aadup[2], ]
table(Partition_6_random_chunk_cv_df2$owner[which(duplicated(Partition_6_random_chunk_cv_df2$tile_name))])

which(duplicated(Partition_6_random_chunk_cv_df2$tile_name))
ChIPATLAS_features_partition6 <- ChIPATLAS_features[match(Partition_6_dfs$tile_name, rownames(ChIPATLAS_features)),]
ChIPATLAS_features_partition6 <- ChIPATLAS_features_partition6[,-c(grep("RNA_polymerase_II__315.bed", colnames(ChIPATLAS_features_partition6)))]
ChIPATLAS_features_partition6 <- ChIPATLAS_features_partition6[,-c(grep("DNase-Seq__46.bed", colnames(ChIPATLAS_features_partition6)))]

ChIPATLAS_features_partition6[is.na(ChIPATLAS_features_partition6)] <- 0
save(list=c("ChIPATLAS_features_partition6"), file = "~/Documents/Shayan/BioInf/lncRNA/ChIPATLAS_features_partition6.RData")

Partition_6_TF_scanned <- read.table("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/All_TF_profile_P6.txt", header = T)
Partition_6_RBP_scanned <- read.table("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/All_RBP_profile_P6.txt", header = T)

colnames(Tile_kmer_features_partition6_revComp) <- paste0(colnames(Tile_kmer_features_partition6_revComp), "__5mer")
colnames(Partition_6_repeat_features) <-  paste0(colnames(Partition_6_repeat_features), "__rep")
colnames(Partition_6_TF_scanned) <-  paste0(colnames(Partition_6_TF_scanned), "__TFscan")
colnames(Partition_6_RBP_scanned) <-  paste0(colnames(Partition_6_RBP_scanned), "__RBPscan")
aachrom <- grep(pattern = '^H[0-9]{1,2}', x = (colnames(ChIPATLAS_features_partition6)))
aachip <- c(1:ncol(ChIPATLAS_features_partition6))[-aachrom]
partition_6_chromatin_features <-  ChIPATLAS_features_partition6[, aachrom]
partition_6_chipOnly_features <-  ChIPATLAS_features_partition6[, aachip]
save(list=c("partition_6_chromatin_features"), file = "~/Documents/Shayan/BioInf/lncRNA/partition_6_chromatin_features.RData")
save(list=c("partition_6_chipOnly_features"), file = "~/Documents/Shayan/BioInf/lncRNA/partition_6_chipOnly_features.RData")

colnames(partition_6_chromatin_features) <-  paste0(colnames(partition_6_chromatin_features), "__HISTONE")
colnames(partition_6_chipOnly_features) <-  paste0(colnames(partition_6_chipOnly_features), "__ChIP")

#colnames(partition_6_expression_features) <-  paste0(colnames(partition_6_expression_features), "__RBPscan")
Partition_6_U1snRNA_mat <- matrix(nrow = length(Partition_6_U1snRNA), ncol = 1)
colnames(Partition_6_U1snRNA_mat) <- "U1snRNA_site__RBPscan"
Partition_6_U1snRNA_mat[, 1] <-  Partition_6_U1snRNA

Partition_6_methylation_mat <- matrix(nrow = length(Partition_6_methylation), ncol = 1)
Partition_6_methylation_mat[, 1] <- Partition_6_methylation
colnames(Partition_6_methylation_mat) <- "Methylation__WGBS"

colnames(Partition_6_accessiblity) <- paste0(colnames(Partition_6_accessiblity), "__AXSB")

partition_6_expression_features_only <- partition_6_expression_features[, c(3,4,5,7,8)]
colnames(partition_6_expression_features_only) <- paste0(colnames(partition_6_expression_features_only), "__TRNSCRP")

distance_feature_partition_6_mat <- matrix(nrow = length(distance_feature_partition_6), ncol = 1)
distance_feature_partition_6_mat[, 1] <- distance_feature_partition_6
colnames(distance_feature_partition_6_mat) <- "Genomic_distance"

colnames(Partition_6_dinuc_freq) <- paste0(colnames(Partition_6_dinuc_freq), "__DiFreq")
Partition_6_feature_mat_tile_2 <- cbind(distance_feature_partition_6_mat,
                                        Partition_6_dinuc_freq,
                                        Tile_kmer_features_partition6_revComp,
                                        Partition_6_repeat_features, 
                                        Partition_6_TF_scanned, 
                                        Partition_6_RBP_scanned,
                                        Partition_6_U1snRNA_mat,
                                        Partition_6_methylation_mat,
                                        Partition_6_accessiblity,
                                        partition_6_chromatin_features,
                                        partition_6_chipOnly_features,
                                        partition_6_expression_features_only)

colnames(Partition_6_feature_mat_tile_2)[which(duplicated(colnames(Partition_6_feature_mat_tile_2)))] <- paste0(colnames(Partition_6_feature_mat_tile_2)[which(duplicated(colnames(Partition_6_feature_mat_tile_2)))], "_v2")
# make the pair features: RBP_pair, TF_pair, repeat_pair, chip_pair, chromatin_pair
aa_lncNames <- read.table("~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/mES/MM9_lncRNA_28.bed", stringsAsFactors = F)$V7

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_RBP_scanPair.RData")
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_TF_scanPair.RData")
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_Repeat_pair.RData")
load("~/Documents/Shayan/BioInf/lncRNA/Tile_owner_pair_RBP_ChIP_pair.RData")
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_network_pairs.RData")

Triplex_feaures_partition6_extended_merge
# colnames(Partition_5_feature_mat_owner) <- paste(colnames(Partition_5_feature_mat_owner), "lncRNA", sep = "__")

Partition_6_feature_mat_pair_2 <- cbind(Partition_6_RBP_scanPair, 
                                        Partition_6_TF_scanPair,
                                        Partition_6_Repeat_pair,
                                        Partition_6_network_pairs,
                                        Triplex_feaures_partition6_extended_merge[match(Partition_6_random_chunk_cv_df$tile_name , 
                                                                                        Triplex_feaures_partition6_extended_merge$tile_name), c(5:8)],
                                        ChIP_tile_owner_pair_feature[match(Partition_6_random_chunk_cv_df$tile_name ,  names(ChIP_tile_owner_pair_feature))],
                                        Chromatin_tile_owner_pair_feature[match(Partition_6_random_chunk_cv_df$tile_name ,  names(Chromatin_tile_owner_pair_feature))])
colnames(Partition_6_feature_mat_pair_2) <- c("RBP_pair","TFscan_pair", "repeat_pair", "PPI_pair", "triplex_total","triplex_GA","triplex_TC","triplex_GT","ChIP_pair", "Chromatin_pair")


# aa_all_feat <- cbind(Partition_6_feature_mat_tile_2,
#                      Partition_6_feature_mat_pair_2,
#                      Partition_6_dfs$label)
# colnames(aa_all_feat)[ncol(aa_all_feat)] <- "label"
# #Partition_6_dfs$owner <- levels(Partition_6_dfs$owner)[as.numeric(Partition_6_dfs$owner )]
#rm(list=setdiff(ls(), c("Partition_6_random_chunk_cv_df", "Partition_6_feature_mat_tile_2", "Partition_6_feature_mat_pair_2")))
aa_un_own <- unique(Partition_6_random_chunk_cv_df$owner)
#Zero_variance_columns_list <- list()
name_dic_list <- list()
#Highly_correlating_columns_list <- list()
for(i in 1:length(aa_un_own)){
  print(i)
  aaw <- which(Partition_6_random_chunk_cv_df$owner %in% aa_un_own[i])
  stopifnot(sum(duplicated(Partition_6_random_chunk_cv_df$tile_name[aaw])) == 0)
  aacurf  <- cbind(Partition_6_feature_mat_tile_2[aaw,],
                   Partition_6_feature_mat_pair_2[aaw,],
                   Partition_6_random_chunk_cv_df$label[aaw])
  colnames(aacurf)[ncol(aacurf)] <- "label"
  aacurf[is.na(aacurf)] <- 0
  rownames(aacurf) <- Partition_6_random_chunk_cv_df$tile_name[aaw]
  #  print("computing variance ...")
  #  aa_zerovar <- apply(aacurf, MARGIN = 2, FUN = var)
  #  Zero_variance_columns_list[[i]] <- colnames(aacurf)[aa_zerovar == 0]
  #  aacurf <- aacurf[, -which(aa_zerovar == 0)]
  
  # aadescrCor <- cor(aacurf[, 1:(ncol(aacurf) - 1)])
  # print("computing correlation ...")
  # aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .90)
  # Highly_correlating_columns_list[[i]] <- colnames(aacurf)[aahighlyCorDescr]
  # aacurf <- aacurf[,-aahighlyCorDescr]
  # aacurf <- as.data.frame(aacurf)
  aacurf$label[aacurf$label == 0] <- "neg"
  aacurf$label[aacurf$label == 1] <- "pos"
  aacurf$label <- factor(aacurf$label , levels = c("neg", "pos"))
  name_dic_list[[i]]<- cbind(colnames(aacurf), paste0("feature_", c(1:ncol(aacurf))))
  name_dic_list[[i]][ncol(aacurf), 2] <- name_dic_list[[i]][ncol(aacurf), 1] 
  colnames(aacurf) <-name_dic_list[[i]][, 2]
  my_Dataset <- aacurf
  my_partition <- Partition_6_random_chunk_cv_df[aaw,]
  my_name_dic <- name_dic_list[[i]]
  print("saving  ...")
  save(list = c("my_Dataset", "my_partition", "my_name_dic"),
       file = paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/Partition_6_CVfirst_dataset_",
                     aa_un_own[i], ".RData"))
}

# for(i in 1:length(aa_un_own)){
  # load(paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/Partition_6_CVfirst_dataset_",
  #             aa_un_own[i], ".RData"))
#   aaw <- which(Partition_6_random_chunk_cv_df$owner %in% aa_un_own[i])
#   rownames(my_Dataset) <- Partition_6_random_chunk_cv_df$tile_name[aaw]
#   save(list = c("my_Dataset", "my_partition", "my_name_dic"),
#        file = paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/Partition_6_CVfirst_dataset_",
#                      aa_un_own[i], ".RData"))
# }
# filter and jobs




aa_filter_list <- list()
aa_feature_remove_list <- list()
aa_colnn <- c(colnames(Partition_6_feature_mat_tile_2), colnames(Partition_6_feature_mat_pair_2))
# 1) RBP + RBP_pair  "__RBPscan" "RBP_pair" 
aafeat1 <- grep(pattern = "__RBPscan", x = aa_colnn)
aafeat2 <- grep(pattern = "RBP_pair", x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[1]] <- aa_remove
aa_filter_list[[1]] <- numeric(0)
# 2) RBP + RBP_pair  "__RBPscan" "RBP_pair"  expression filtered
aafeat1 <- grep(pattern = "__RBPscan", x = aa_colnn)
aafeat2 <- grep(pattern = "RBP_pair", x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[2]] <- aa_remove
aa_filter_list[[2]] <- c(aafeat1, aafeat2)
# 3) TF + TF_pair  "__TFscan" "TFscan_pair" 
aafeat1 <- grep(pattern = "__TFscan", x = aa_colnn)
aafeat2 <- grep(pattern = "TFscan_pair", x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[3]] <- aa_remove
aa_filter_list[[3]] <- numeric(0)
# 4) TF + TF_pair  "__TFscan" "TFscan_pair"   expression filtered
aafeat1 <- grep(pattern = "__TFscan", x = aa_colnn)
aafeat2 <- grep(pattern = "TFscan_pair", x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[4]] <- aa_remove
aa_filter_list[[4]] <- c(aafeat1, aafeat2)
# 5) repeat + repeat_pair  "__rep" "repeat_pair" 
aafeat1 <- grep(pattern = "__rep", x = aa_colnn)
aafeat2 <- grep(pattern = "repeat_pair", x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[5]] <- aa_remove
aa_filter_list[[5]] <- numeric(0)
# 6) repeat + repeat_pair  "__rep" "repeat_pair"  expression filtered
aafeat1 <- grep(pattern = "__rep", x = aa_colnn)
aafeat2 <- grep(pattern = "repeat_pair", x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[6]] <- aa_remove
aa_filter_list[[6]] <- c(aafeat1, aafeat2)
# 7) RBP + RBP_pair + TF + TF_pair  repeat + repeat_pair "__rep" "repeat_pair" 
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), grep(pattern = "__TFscan", x = aa_colnn), grep(pattern = "__rep", x = aa_colnn))
aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn), grep(pattern = "TFscan_pair", x = aa_colnn), grep(pattern = "repeat_pair", x = aa_colnn), grep(pattern = "PPI_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[7]] <- aa_remove
aa_filter_list[[7]] <- numeric(0)
# 8) RBP + RBP_pair + TF + TF_pair  repeat + repeat_pair "__rep" "repeat_pair" expression filtered
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), grep(pattern = "__TFscan", x = aa_colnn), grep(pattern = "__rep", x = aa_colnn))
aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn), grep(pattern = "TFscan_pair", x = aa_colnn), grep(pattern = "repeat_pair", x = aa_colnn), grep(pattern = "PPI_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[8]] <- aa_remove
aa_filter_list[[8]] <- c(grep(pattern = "__RBPscan", x = aa_colnn),
                         grep(pattern = "RBP_pair", x = aa_colnn))
# 9) kmer_only
aafeat1 <- c(grep(pattern = "__5mer", x = aa_colnn))
#aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn), grep(pattern = "TFscan_pair", x = aa_colnn), grep(pattern = "repeat_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), aafeat1)
aa_feature_remove_list[[9]] <- aa_remove
aa_filter_list[[9]] <- numeric(0)

# 10) RBP + RBP_pair + TF + TF_pair  repeat + repeat_pair "__rep" "repeat_pair"  + triplex
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), grep(pattern = "__TFscan", x = aa_colnn), grep(pattern = "__rep", x = aa_colnn))
aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn), 
             grep(pattern = "TFscan_pair", x = aa_colnn),
             grep(pattern = "repeat_pair", x = aa_colnn), 
             grep(pattern = "PPI_pair", x = aa_colnn), 
             grep(pattern = "triplex", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[10]] <- aa_remove
aa_filter_list[[10]] <- numeric(0)

# 11) RBP + RBP_pair + TF + TF_pair  repeat + repeat_pair "__rep" "repeat_pair" expression filtered + triplex
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), grep(pattern = "__TFscan", x = aa_colnn), grep(pattern = "__rep", x = aa_colnn))
aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn), 
             grep(pattern = "TFscan_pair", x = aa_colnn),
             grep(pattern = "repeat_pair", x = aa_colnn),
             grep(pattern = "PPI_pair", x = aa_colnn), 
             grep(pattern = "triplex", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[11]] <- aa_remove
aa_filter_list[[11]] <- c(grep(pattern = "__RBPscan", x = aa_colnn),
                          grep(pattern = "RBP_pair", x = aa_colnn))
# 12) chromatin + chromtin_pair
aafeat1 <- grep(pattern = 'HISTONE', x = aa_colnn)
aafeat2 <- grep(pattern = "Chromatin_pair", x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[12]] <- aa_remove
aa_filter_list[[12]] <- numeric(0)

# 13) Chip + ChIP_pair
aachip <- grep(pattern = 'ChIP', x = aa_colnn)
#aabed <- grep(pattern = '*.bed', x = aa_colnn)
#aabed <- setdiff(aabed, aachrom)
aafeat1 <- aachip
aafeat2 <- grep(pattern = "ChIP_pair", x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[13]] <- aa_remove
aa_filter_list[[13]] <- numeric(0)
# 14) expression
aafeat1 <- grep(pattern = 'TRNSCRP', x = aa_colnn)
# aaexp_names <- c("RNA_polymerase_II__315.bed", "RNAseq" ,"CAGE")
# aafeat1 <- which(aa_colnn %in%  aaexp_names)
aa_remove <- setdiff(c(1:length(aa_colnn)), aafeat1)
aa_feature_remove_list[[14]] <- aa_remove
aa_filter_list[[14]] <- numeric(0)

#15 accessibility
aafeat1 <- grep(pattern = '_AXSB', x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), aafeat1)
aa_feature_remove_list[[15]] <- aa_remove
aa_filter_list[[15]] <- numeric(0)
#16 methylation
aafeat1 <- grep(pattern = 'Methylation__WGBS', x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), aafeat1)
aa_feature_remove_list[[16]] <- aa_remove
aa_filter_list[[16]] <- numeric(0)

#17 distance
aafeat1 <- grep(pattern = 'Genomic_distance', x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), aafeat1)
aa_feature_remove_list[[17]] <- aa_remove
aa_filter_list[[17]] <- numeric(0)

#18 dinuc_freq
aafeat1 <- grep(pattern = '__DiFreq', x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), aafeat1)
aa_feature_remove_list[[18]] <- aa_remove
aa_filter_list[[18]] <- numeric(0)


#19 dinuc_freq + distance
aafeat1 <- grep(pattern = '__DiFreq', x = aa_colnn)
aafeat2 <- grep(pattern = 'Genomic_distance', x = aa_colnn)
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[19]] <- aa_remove
aa_filter_list[[19]] <- numeric(0)

#20 dinuc_freq + distance + Methylation__WGBS + accessibility
aafeat1 <- grep(pattern = '__DiFreq', x = aa_colnn)
aafeat2 <- c(grep(pattern = 'Genomic_distance', x = aa_colnn) ,
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn) )
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[20]] <- aa_remove
aa_filter_list[[20]] <- numeric(0)

#21 dinuc_freq + distance + Methylation__WGBS + accessibility + histone
# aafeat1 <- grep(pattern = 'HISTONE', x = aa_colnn)
# aafeat2 <- grep(pattern = "Chromatin_pair", x = aa_colnn)

aafeat1 <- grep(pattern = '__DiFreq', x = aa_colnn)
aafeat2 <- c(grep(pattern = 'Genomic_distance', x = aa_colnn) ,
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn), 
             grep(pattern = 'HISTONE', x = aa_colnn), 
             grep(pattern = "Chromatin_pair", x = aa_colnn) )
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[21]] <- aa_remove
aa_filter_list[[21]] <- numeric(0)

#22 kmer_freq + distance + Methylation__WGBS + accessibility + histone
# aafeat1 <- grep(pattern = 'HISTONE', x = aa_colnn)
# aafeat2 <- grep(pattern = "Chromatin_pair", x = aa_colnn)

aafeat1 <- grep(pattern = '__5mer', x = aa_colnn)
aafeat2 <- c(grep(pattern = 'Genomic_distance', x = aa_colnn) ,
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn), 
             grep(pattern = 'HISTONE', x = aa_colnn), 
             grep(pattern = "Chromatin_pair", x = aa_colnn) )
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[22]] <- aa_remove
aa_filter_list[[22]] <- numeric(0)

#23 accessibility + methylation
aafeat2 <- c(#grep(pattern = 'Genomic_distance', x = aa_colnn) ,
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn)
             #,grep(pattern = 'HISTONE', x = aa_colnn), 
             #grep(pattern = "Chromatin_pair", x = aa_colnn)
             )
aa_remove <- setdiff(c(1:length(aa_colnn)),aafeat2)
aa_feature_remove_list[[23]] <- aa_remove
aa_filter_list[[23]] <- numeric(0)
#24 accessibility + methylation + chromatin + chromtin_pair 
aafeat2 <- c(#grep(pattern = 'Genomic_distance', x = aa_colnn) ,
  grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
  grep(pattern = '_AXSB', x = aa_colnn)
  ,grep(pattern = 'HISTONE', x = aa_colnn), 
  grep(pattern = "Chromatin_pair", x = aa_colnn)
)
aa_remove <- setdiff(c(1:length(aa_colnn)),aafeat2)
aa_feature_remove_list[[24]] <- aa_remove
aa_filter_list[[24]] <- numeric(0)

#25 accessibility + methylation + chromatin + chromtin_pair + Chip + ChIP_pair 

aafeat2 <- c(#grep(pattern = 'Genomic_distance', x = aa_colnn) ,
  grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
  grep(pattern = '_AXSB', x = aa_colnn)
  ,grep(pattern = 'HISTONE', x = aa_colnn), 
  grep(pattern = "Chromatin_pair", x = aa_colnn), 
  grep(pattern = 'ChIP', x = aa_colnn), 
  grep(pattern = "ChIP_pair", x = aa_colnn) 
  #,grep(pattern = "PPI_pair", x = aa_colnn)
)
aa_remove <- setdiff(c(1:length(aa_colnn)),aafeat2)
aa_feature_remove_list[[25]] <- aa_remove
aa_filter_list[[25]] <- numeric(0)




# 26) accessibility + methylation + chromatin + chromtin_pair + Chip + ChIP_pair + expression 
#aafeat1 <- c(grep(pattern = '*.bed', x = aa_colnn),  which(aa_colnn %in%  aaexp_names))
aafeat2 <- c(#grep(pattern = 'Genomic_distance', x = aa_colnn) ,
  grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
  grep(pattern = '_AXSB', x = aa_colnn)
  ,grep(pattern = 'HISTONE', x = aa_colnn), 
  grep(pattern = "Chromatin_pair", x = aa_colnn), 
  grep(pattern = 'ChIP', x = aa_colnn), 
  grep(pattern = "ChIP_pair", x = aa_colnn), 
  #grep(pattern = "PPI_pair", x = aa_colnn), 
  grep(pattern = "TRNSCRP", x = aa_colnn)
)
aa_remove <- setdiff(c(1:length(aa_colnn)),aafeat2)
aa_feature_remove_list[[26]] <- aa_remove
aa_filter_list[[26]] <- numeric(0)

# 27) RBP-TF-Repeat-triplex-histone-accessibility-methylation
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), 
             grep(pattern = "__TFscan", x = aa_colnn),
             grep(pattern = "__rep", x = aa_colnn),
             grep(pattern = "triplex", x = aa_colnn),
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn)
             ,grep(pattern = 'HISTONE', x = aa_colnn), 
             grep(pattern = "Chromatin_pair", x = aa_colnn), 
             # grep(pattern = 'ChIP', x = aa_colnn), 
             # grep(pattern = "ChIP_pair", x = aa_colnn), 
             grep(pattern = "PPI_pair", x = aa_colnn)
             )
             #grep(pattern = '*.bed', x = aa_colnn),
        #     which(aa_colnn %in%  aaexp_names))
aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn),
             grep(pattern = "TFscan_pair", x = aa_colnn),
             grep(pattern = "repeat_pair", x = aa_colnn))
             #grep(pattern = "triplex", x = aa_colnn),
             #grep(pattern = "ChIP_pair", x = aa_colnn),
             #grep(pattern = "Chromatin_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[27]] <- aa_remove
aa_filter_list[[27]] <- numeric(0)



# 28) (27) + expression filter
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), 
             grep(pattern = "__TFscan", x = aa_colnn),
             grep(pattern = "__rep", x = aa_colnn),
             grep(pattern = "triplex", x = aa_colnn),
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn)
             ,grep(pattern = 'HISTONE', x = aa_colnn), 
             grep(pattern = "Chromatin_pair", x = aa_colnn), 
             # grep(pattern = 'ChIP', x = aa_colnn), 
             # grep(pattern = "ChIP_pair", x = aa_colnn), 
             grep(pattern = "PPI_pair", x = aa_colnn)
)

aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn),
             grep(pattern = "TFscan_pair", x = aa_colnn),
             grep(pattern = "repeat_pair", x = aa_colnn))

aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[28]] <- aa_remove
aa_filter_list[[28]] <-  c(grep(pattern = "__RBPscan", x = aa_colnn),
                           grep(pattern = "RBP_pair", x = aa_colnn))
#29) 27 + chip + chippair
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), 
             grep(pattern = "__TFscan", x = aa_colnn),
             grep(pattern = "__rep", x = aa_colnn),
             grep(pattern = "triplex", x = aa_colnn),
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn)
             ,grep(pattern = 'HISTONE', x = aa_colnn), 
             grep(pattern = "Chromatin_pair", x = aa_colnn), 
              grep(pattern = 'ChIP', x = aa_colnn), 
              grep(pattern = "ChIP_pair", x = aa_colnn), 
             grep(pattern = "PPI_pair", x = aa_colnn)
)

aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn),
             grep(pattern = "TFscan_pair", x = aa_colnn),
             grep(pattern = "repeat_pair", x = aa_colnn))

aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[29]] <- aa_remove
aa_filter_list[[29]] <- numeric(0)

#30) 29 + transcription
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), 
             grep(pattern = "__TFscan", x = aa_colnn),
             grep(pattern = "__rep", x = aa_colnn),
             grep(pattern = "triplex", x = aa_colnn),
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn)
             ,grep(pattern = 'HISTONE', x = aa_colnn), 
             grep(pattern = "Chromatin_pair", x = aa_colnn), 
             grep(pattern = 'ChIP', x = aa_colnn), 
             grep(pattern = "ChIP_pair", x = aa_colnn), 
             grep(pattern = "PPI_pair", x = aa_colnn),
             grep(pattern = "TRNSCRP", x = aa_colnn)
)

aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn),
             grep(pattern = "TFscan_pair", x = aa_colnn),
             grep(pattern = "repeat_pair", x = aa_colnn))

aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[30]] <- aa_remove
aa_filter_list[[30]] <- numeric(0)


#31)  29 + dinucletide
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), 
             grep(pattern = "__TFscan", x = aa_colnn),
             grep(pattern = "__rep", x = aa_colnn),
             grep(pattern = "triplex", x = aa_colnn),
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn)
             ,grep(pattern = 'HISTONE', x = aa_colnn), 
             grep(pattern = "Chromatin_pair", x = aa_colnn), 
             grep(pattern = 'ChIP', x = aa_colnn), 
             grep(pattern = "ChIP_pair", x = aa_colnn), 
             grep(pattern = "PPI_pair", x = aa_colnn),
             grep(pattern = '__DiFreq', x = aa_colnn)
)

aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn),
             grep(pattern = "TFscan_pair", x = aa_colnn),
             grep(pattern = "repeat_pair", x = aa_colnn))

aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[31]] <- aa_remove
aa_filter_list[[31]] <- numeric(0)

#32)  29 + 5mer
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), 
             grep(pattern = "__TFscan", x = aa_colnn),
             grep(pattern = "__rep", x = aa_colnn),
             grep(pattern = "triplex", x = aa_colnn),
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn)
             ,grep(pattern = 'HISTONE', x = aa_colnn), 
             grep(pattern = "Chromatin_pair", x = aa_colnn), 
             grep(pattern = 'ChIP', x = aa_colnn), 
             grep(pattern = "ChIP_pair", x = aa_colnn), 
             grep(pattern = "PPI_pair", x = aa_colnn),
             grep(pattern = '__5mer', x = aa_colnn)
)

aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn),
             grep(pattern = "TFscan_pair", x = aa_colnn),
             grep(pattern = "repeat_pair", x = aa_colnn))

aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[32]] <- aa_remove
aa_filter_list[[32]] <- numeric(0)


#33)  29 + dinucletide +TXfiltered
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), 
             grep(pattern = "__TFscan", x = aa_colnn),
             grep(pattern = "__rep", x = aa_colnn),
             grep(pattern = "triplex", x = aa_colnn),
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn)
             ,grep(pattern = 'HISTONE', x = aa_colnn), 
             grep(pattern = "Chromatin_pair", x = aa_colnn), 
             grep(pattern = 'ChIP', x = aa_colnn), 
             grep(pattern = "ChIP_pair", x = aa_colnn), 
             grep(pattern = "PPI_pair", x = aa_colnn),
             grep(pattern = '__DiFreq', x = aa_colnn)
)

aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn),
             grep(pattern = "TFscan_pair", x = aa_colnn),
             grep(pattern = "repeat_pair", x = aa_colnn))

aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[33]] <- aa_remove
aa_filter_list[[33]] <-  c(grep(pattern = "__RBPscan", x = aa_colnn),
                           grep(pattern = "RBP_pair", x = aa_colnn))


#34)  29 + 5mer +TXfiltered
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), 
             grep(pattern = "__TFscan", x = aa_colnn),
             grep(pattern = "__rep", x = aa_colnn),
             grep(pattern = "triplex", x = aa_colnn),
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn)
             ,grep(pattern = 'HISTONE', x = aa_colnn), 
             grep(pattern = "Chromatin_pair", x = aa_colnn), 
             grep(pattern = 'ChIP', x = aa_colnn), 
             grep(pattern = "ChIP_pair", x = aa_colnn), 
             grep(pattern = "PPI_pair", x = aa_colnn),
             grep(pattern = '__5mer', x = aa_colnn)
)

aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn),
             grep(pattern = "TFscan_pair", x = aa_colnn),
             grep(pattern = "repeat_pair", x = aa_colnn))

aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[34]] <- aa_remove
aa_filter_list[[34]] <-  c(grep(pattern = "__RBPscan", x = aa_colnn),
                           grep(pattern = "RBP_pair", x = aa_colnn))

#35)  29 + 5mer + dinucletide
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), 
             grep(pattern = "__TFscan", x = aa_colnn),
             grep(pattern = "__rep", x = aa_colnn),
             grep(pattern = "triplex", x = aa_colnn),
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn)
             ,grep(pattern = 'HISTONE', x = aa_colnn), 
             grep(pattern = "Chromatin_pair", x = aa_colnn), 
             grep(pattern = 'ChIP', x = aa_colnn), 
             grep(pattern = "ChIP_pair", x = aa_colnn), 
             grep(pattern = "PPI_pair", x = aa_colnn),
             grep(pattern = '__5mer', x = aa_colnn),
             grep(pattern = '__DiFreq', x = aa_colnn)
)

aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn),
             grep(pattern = "TFscan_pair", x = aa_colnn),
             grep(pattern = "repeat_pair", x = aa_colnn))

aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[35]] <- aa_remove
aa_filter_list[[35]] <- numeric(0)

#36)  29 + 5mer + dinucletide + distance
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), 
             grep(pattern = "__TFscan", x = aa_colnn),
             grep(pattern = "__rep", x = aa_colnn),
             grep(pattern = "triplex", x = aa_colnn),
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn)
             ,grep(pattern = 'HISTONE', x = aa_colnn), 
             grep(pattern = "Chromatin_pair", x = aa_colnn), 
             grep(pattern = 'ChIP', x = aa_colnn), 
             grep(pattern = "ChIP_pair", x = aa_colnn), 
             grep(pattern = "PPI_pair", x = aa_colnn),
             grep(pattern = '__5mer', x = aa_colnn),
             grep(pattern = '__DiFreq', x = aa_colnn),
             grep(pattern = 'Genomic_distance', x = aa_colnn)
)

aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn),
             grep(pattern = "TFscan_pair", x = aa_colnn),
             grep(pattern = "repeat_pair", x = aa_colnn))

aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[36]] <- aa_remove
aa_filter_list[[36]] <- numeric(0)

#37)  29 + 5mer + dinucletide + expression 
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), 
             grep(pattern = "__TFscan", x = aa_colnn),
             grep(pattern = "__rep", x = aa_colnn),
             grep(pattern = "triplex", x = aa_colnn),
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn)
             ,grep(pattern = 'HISTONE', x = aa_colnn), 
             grep(pattern = "Chromatin_pair", x = aa_colnn), 
             grep(pattern = 'ChIP', x = aa_colnn), 
             grep(pattern = "ChIP_pair", x = aa_colnn), 
             grep(pattern = "PPI_pair", x = aa_colnn),
             grep(pattern = '__5mer', x = aa_colnn),
             grep(pattern = '__DiFreq', x = aa_colnn),
             #grep(pattern = 'Genomic_distance', x = aa_colnn),
             grep(pattern = 'TRNSCRP', x = aa_colnn)
             
)

aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn),
             grep(pattern = "TFscan_pair", x = aa_colnn),
             grep(pattern = "repeat_pair", x = aa_colnn))

aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[37]] <- aa_remove
aa_filter_list[[37]] <- numeric(0)


#38)  29 + 5mer + dinucletide + expression + tx filter
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), 
             grep(pattern = "__TFscan", x = aa_colnn),
             grep(pattern = "__rep", x = aa_colnn),
             grep(pattern = "triplex", x = aa_colnn),
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn)
             ,grep(pattern = 'HISTONE', x = aa_colnn), 
             grep(pattern = "Chromatin_pair", x = aa_colnn), 
             grep(pattern = 'ChIP', x = aa_colnn), 
             grep(pattern = "ChIP_pair", x = aa_colnn), 
             grep(pattern = "PPI_pair", x = aa_colnn),
             grep(pattern = '__5mer', x = aa_colnn),
             grep(pattern = '__DiFreq', x = aa_colnn),
             #grep(pattern = 'Genomic_distance', x = aa_colnn),
             grep(pattern = 'TRNSCRP', x = aa_colnn)
             
)

aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn),
             grep(pattern = "TFscan_pair", x = aa_colnn),
             grep(pattern = "repeat_pair", x = aa_colnn))

aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[38]] <- aa_remove
aa_filter_list[[38]] <- c(grep(pattern = "__RBPscan", x = aa_colnn),
                          grep(pattern = "RBP_pair", x = aa_colnn))
#39)  29 + 5mer + dinucletide + expression + distance
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), 
             grep(pattern = "__TFscan", x = aa_colnn),
             grep(pattern = "__rep", x = aa_colnn),
             grep(pattern = "triplex", x = aa_colnn),
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn)
             ,grep(pattern = 'HISTONE', x = aa_colnn), 
             grep(pattern = "Chromatin_pair", x = aa_colnn), 
             grep(pattern = 'ChIP', x = aa_colnn), 
             grep(pattern = "ChIP_pair", x = aa_colnn), 
             grep(pattern = "PPI_pair", x = aa_colnn),
             grep(pattern = '__5mer', x = aa_colnn),
             grep(pattern = '__DiFreq', x = aa_colnn),
             grep(pattern = 'Genomic_distance', x = aa_colnn),
             grep(pattern = 'TRNSCRP', x = aa_colnn)
             
)

aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn),
             grep(pattern = "TFscan_pair", x = aa_colnn),
             grep(pattern = "repeat_pair", x = aa_colnn))

aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[39]] <- aa_remove
aa_filter_list[[39]] <- numeric(0)
#40)  29 + 5mer + dinucletide + expression + distance + txfilter
aafeat1 <- c(grep(pattern = "__RBPscan", x = aa_colnn), 
             grep(pattern = "__TFscan", x = aa_colnn),
             grep(pattern = "__rep", x = aa_colnn),
             grep(pattern = "triplex", x = aa_colnn),
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn)
             ,grep(pattern = 'HISTONE', x = aa_colnn), 
             grep(pattern = "Chromatin_pair", x = aa_colnn), 
             grep(pattern = 'ChIP', x = aa_colnn), 
             grep(pattern = "ChIP_pair", x = aa_colnn), 
             grep(pattern = "PPI_pair", x = aa_colnn),
             grep(pattern = '__5mer', x = aa_colnn),
             grep(pattern = '__DiFreq', x = aa_colnn),
             grep(pattern = 'Genomic_distance', x = aa_colnn),
             grep(pattern = 'TRNSCRP', x = aa_colnn)
             
)

aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn),
             grep(pattern = "TFscan_pair", x = aa_colnn),
             grep(pattern = "repeat_pair", x = aa_colnn))

aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[40]] <- aa_remove
aa_filter_list[[40]] <-c(grep(pattern = "__RBPscan", x = aa_colnn),
                         grep(pattern = "RBP_pair", x = aa_colnn))




# 41) kmer_and chromatin
aafeat1 <- c(grep(pattern = "__5mer", x = aa_colnn))

aafeat2 <- c(grep(pattern = 'HISTONE', x = aa_colnn),
             grep(pattern = "Chromatin_pair", x = aa_colnn))
#aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn), grep(pattern = "TFscan_pair", x = aa_colnn), grep(pattern = "repeat_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[41]] <- aa_remove
aa_filter_list[[41]] <- numeric(0)

# 42) kmer_and chromatin + methylation + accessibility
aafeat1 <- c(grep(pattern = "__5mer", x = aa_colnn))

aafeat2 <- c(grep(pattern = 'HISTONE', x = aa_colnn),
             grep(pattern = "Chromatin_pair", x = aa_colnn),
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn))
#aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn), grep(pattern = "TFscan_pair", x = aa_colnn), grep(pattern = "repeat_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[42]] <- aa_remove
aa_filter_list[[42]] <- numeric(0)

# 43) kmer_and chromatin + methylation + accessibility + chip
aafeat1 <- c(grep(pattern = "__5mer", x = aa_colnn))

aafeat2 <- c(grep(pattern = 'HISTONE', x = aa_colnn),
             grep(pattern = "Chromatin_pair", x = aa_colnn),
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn),       
             grep(pattern = 'ChIP', x = aa_colnn), 
             grep(pattern = "ChIP_pair", x = aa_colnn))
#aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn), grep(pattern = "TFscan_pair", x = aa_colnn), grep(pattern = "repeat_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[43]] <- aa_remove
aa_filter_list[[43]] <- numeric(0)
# 44) kmer_and chromatin + methylation + accessibility + chip +transcription
aafeat1 <- c(grep(pattern = "__5mer", x = aa_colnn))

aafeat2 <- c(grep(pattern = 'HISTONE', x = aa_colnn),
             grep(pattern = "Chromatin_pair", x = aa_colnn),
             grep(pattern = 'Methylation__WGBS', x = aa_colnn), 
             grep(pattern = '_AXSB', x = aa_colnn),       
             grep(pattern = 'ChIP', x = aa_colnn), 
             grep(pattern = "ChIP_pair", x = aa_colnn),
             grep(pattern = 'TRNSCRP', x = aa_colnn))
#aafeat2 <- c(grep(pattern = "RBP_pair", x = aa_colnn), grep(pattern = "TFscan_pair", x = aa_colnn), grep(pattern = "repeat_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[44]] <- aa_remove
aa_filter_list[[44]] <- numeric(0)


#45) kmer_and triplex
aafeat1 <- c(grep(pattern = "__5mer", x = aa_colnn))
aafeat2 <- c(grep(pattern = 'triplex', x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[45]] <- aa_remove
aa_filter_list[[45]] <- numeric(0)

# 46) kmer_and triplex and pairs
aafeat1 <- c(grep(pattern = "__5mer", x = aa_colnn))
aafeat2 <- c(grep(pattern = 'triplex', x = aa_colnn),
             grep("_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)), c(aafeat1, aafeat2))
aa_feature_remove_list[[46]] <- aa_remove
aa_filter_list[[46]] <- numeric(0)

#47) triplex only
aafeat2 <- c(grep(pattern = 'triplex', x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)),  aafeat2)
aa_feature_remove_list[[47]] <- aa_remove
aa_filter_list[[47]] <- numeric(0)

#48) pairs only
aafeat2 <- c(grep("_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)),  aafeat2)
aa_feature_remove_list[[48]] <- aa_remove
aa_filter_list[[48]] <- numeric(0)


# 49) triplex and pairs
#aafeat1 <- c(grep(pattern = "__5mer", x = aa_colnn))
aafeat2 <- c(grep(pattern = 'triplex', x = aa_colnn),
             grep("_pair", x = aa_colnn))
aa_remove <- setdiff(c(1:length(aa_colnn)),  aafeat2)
aa_feature_remove_list[[49]] <- aa_remove
aa_filter_list[[49]] <- numeric(0)



# dinuc_freq + distance
# dinuc_freq + distance + Methylation__WGBS + accessibility
# dinuc_freq + distance + Methylation__WGBS + accessibility + histone
# kmer_freq + distance + Methylation__WGBS + accessibility + histone
# accessibility + methylation
# accessibility + methylation + chromatin + chromtin_pair 
# accessibility + methylation + chromatin + chromtin_pair + Chip + ChIP_pair 
# accessibility + methylation + chromatin + chromtin_pair + Chip + ChIP_pair + expression 
# RBP-TF-Repeat-triplex-histone-accessibility-methylation
# RBP-TF-Repeat-triplex-histone-accessibility-methylation + expression filter
# RBP-TF-Repeat-triplex-histone-accessibility-methylation + chip + chippair
# RBP-TF-Repeat-triplex-histone-accessibility-methylation + chip + chippair  + expression 
# 
# RBP-TF-Repeat-triplex-histone-accessibility-methylation + chip + chippair  + dinucletide
# RBP-TF-Repeat-triplex-histone-accessibility-methylation + chip + chippair  + kmer
# RBP-TF-Repeat-triplex-histone-accessibility-methylation + chip + chippair  + dinucletide + TXfilter
# RBP-TF-Repeat-triplex-histone-accessibility-methylation + chip + chippair  + kmer + TXfilter
# RBP-TF-Repeat-triplex-histone-accessibility-methylation + chip + chippair  + dinucletide + kmer
# RBP-TF-Repeat-triplex-histone-accessibility-methylation + chip + chippair  + dinucletide + kmer + distance
# RBP-TF-Repeat-triplex-histone-accessibility-methylation + chip + chippair  + dinucletide + kmer + expression
# RBP-TF-Repeat-triplex-histone-accessibility-methylation + chip + chippair  + dinucletide + kmer + expression + TXfilter
# RBP-TF-Repeat-triplex-histone-accessibility-methylation + chip + chippair  + dinucletide + kmer + expression + distance
# RBP-TF-Repeat-triplex-histone-accessibility-methylation + chip + chippair  + dinucletide + kmer + expression + distance + TXfilter
# 
# kmer_and chromatin
# kmer_and chromatin + methylation + accessibility
# kmer_and chromatin + methylation + accessibility + chip
# kmer_and chromatin + methylation + accessibility + chip +transcription
# kmer_and triplex
# kmer_and triplex and pairs
# triplex only
# pairs only
# triplex and pairs


names(aa_filter_list) <- c("RBPscanned", "RBPscannedTX", "TFscanned", "TFscannedTX", "Repeat", "RepeatTX", 
                           "RBPTFRepeat", "RBPTFRepeatTXRBPrep", "Kmer", "RBPTFRepeatTriplex", "RBPTFRepeatTriplexTXRBPrep", "Chromatin",
                           "ChIP", "Transcription",
                           "Accessibility", "Methylation", "DistanceOnly", "DinucFreq",  "DinucDist", "DinucDistMethAccess", "DinucDistMethAccessChrom",
                           "kmerDistMethAccessChrom", "AccessMeth", "AccessMethChrom", "AccessMethChromChIP", "AccessMethChromChIPexp", 
                           "RBPTFRepeatTriplexChromAccessMeth", "RBPTFRepeatTriplexChromAccessMethTX", "RBPTFRepeatTriplexChromAccessMethChIP",
                           "RBPTFRepeatTriplexChromAccessMethChIPExp", "RBPTFRepeatTriplexChromAccessMethChIPDinuc" , "RBPTFRepeatTriplexChromAccessMethChIPkmer",
                           "RBPTFRepeatTriplexChromAccessMethChIPDinucTX", "RBPTFRepeatTriplexChromAccessMethChIPkmerTX", "RBPTFRepeatTriplexChromAccessMethChIPDinucKmer",
                           "RBPTFRepeatTriplexChromAccessMethChIPDinucKmerDistance","RBPTFRepeatTriplexChromAccessMethChIPDinucKmerExp", "RBPTFRepeatTriplexChromAccessMethChIPDinucKmerExpTX",
                           "RBPTFRepeatTriplexChromAccessMethChIPDinucKmerExpDist", "RBPTFRepeatTriplexChromAccessMethChIPDinucKmerExpDistTX",
                           "KmerChrom", "KmerChromMethAccess", "KmerChromMethAccessChIP", "KmerChromMethAccessChIPExp", "KmerTriplex","KmerTriplexPairs","Triplex", "Pairs", "TriplexPairs")
names(aa_feature_remove_list) <- names(aa_filter_list)
# write filter_list files + jobs
aa_un_own <- unique(Partition_6_random_chunk_cv_df$owner)
aa_all_cv <- c(paste0("RCV", c(1:5)), paste0("CCV", c(1:5)))
for(j in 1:length(aa_feature_remove_list)){
  if(length(aa_filter_list[[j]]) == 0){
    aamyfilt <- " no_filter"
  }else{
    write.table(aa_filter_list[[j]], quote = F, row.names = F, col.names = F, 
                file = paste0("~/Documents/Shayan/BioInf/lncRNA/partition_6_CVfirst_TOBE_exp_FILTERED/my_filter_", j,".txt"), append = F)
    aamyfilt <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_CVfirst_TOBE_exp_FILTERED/my_filter_", j,".txt")
  }
  
  for(i in 1:length(aa_un_own)){
    for(aacurcv in 1:length(aa_all_cv)){
      cat(c("Rscript --vanilla RF_load_run_p6_CVfirst_feature.R",
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_CVfirst_dataset_",
                   aa_un_own[i],".RData "), 
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CVfirst_Results/Learned_models/",aa_un_own[i],"/"),
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CVfirst_Results/performance_plots/",aa_un_own[i],"/"),
            aa_all_cv[aacurcv],
            
            sample(c(100:10000), 1),
            names(aa_filter_list)[j], 
            aamyfilt,
            aa_feature_remove_list[[j]],
            "\n"),
          sep = " ", file = "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_CVfirst_full.job", append = T)
    }

  }
  
}

##########
# Write new feature list jobs

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_Dec_dataset/Partition_6_Dec_dataset_Firre.RData")
aa_filter_list2 <- list()
aa_feature_remove_list2 <- list()
aa_colnn <- my_name_dic[1:(nrow(my_name_dic) - 1),1]  #c(colnames(Partition_6_feature_mat_tile_2), colnames(Partition_6_feature_mat_pair_2))




# write filter_list files + jobs
aa_un_own <- unique(Partition_6_dfs$owner)
for(j in 1:length(aa_feature_remove_list2)){
  if(length(aa_filter_list2[[j]]) == 0){
    aamyfilt <- " no_filter"
  }else{
    write.table(aa_filter_list2[[j]], quote = F, row.names = F, col.names = F, 
                file = paste0("~/Documents/Shayan/BioInf/lncRNA/partition_6_Dec_TOBE_exp_FILTERED2/my_filter_", j,".txt"), append = F)
    aamyfilt <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_Dec_TOBE_exp_FILTERED2/my_filter_", j,".txt")
  }
  for(i in 1:length(aa_un_own)){
    cat(c("Rscript --vanilla RF_load_run_p6_Dec_feature.R",
          paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_Dec_dataset_",
                 aa_un_own[i],".RData "), 
          paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_Dec_Results/Learned_models/",aa_un_own[i],"/"),
          paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_Dec_Results/performance_plots/",aa_un_own[i],"/"),
          " rand ",
          sample(c(100:10000), 1),
          names(aa_filter_list2)[j], 
          aamyfilt,
          aa_feature_remove_list2[[j]],
          "\n"),
        sep = " ", file = "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full2.job", append = T)
    
    cat(c("Rscript --vanilla RF_load_run_p6_Dec_feature.R",
          paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_Dec_dataset_",
                 aa_un_own[i],".RData "), 
          paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_Dec_Results/Learned_models/",aa_un_own[i],"/"),
          paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_Dec_Results/performance_plots/",aa_un_own[i],"/"),
          " chunk ",
          sample(c(100:10000), 1),
          names(aa_filter_list2)[j],
          aamyfilt,
          aa_feature_remove_list2[[j]],
          "\n"),
        sep = " ", file = "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_dec_full2.job", append = T)
  }
  
}


#########################################################################################################
# investigate the duplicates
load("Partition_6_random_chunk_cv_df.RData")
load("Partition_6_chunk_random.RData")
head(Partition_6_dfs)
head(Partition_6_random_chunk_cv_df)
sum(Partition_6_random_chunk_cv_df$owner == Partition_6_dfs$owner)
sum(Partition_6_random_chunk_cv_df$label == Partition_6_dfs$label)
sum(Partition_6_random_chunk_cv_df$tile_name == Partition_6_dfs$tile_name)


sum(duplicated(Partition_6_dfs$tile_name))
sum(duplicated(Partition_6_random_chunk_cv_df$tile_name))

for(i in 1:length(Partition_6_chunk_cv)){
  print(names(Partition_6_chunk_cv)[i])
  print(sum(! duplicated(Partition_6_chunk_cv[[i]]$tile_name)))
  print(sum(! duplicated(Partition_6_random_cv[[i]]$tile_name)))
}


#########################################################################################################
# write job files for RF_svd analysis

# should define three variables for each scenario: my_feat_family_list my_family_transform my_family_trans_len

# families:
all_families <- list(dist = c('Genomic_distance'),
                     kmerFq = c('__DiFreq',"__5mer"),
                     motSc = c("__RBPscan", "__TFscan"),
                     Rep = c( "__rep"),
                     pairs = c("RBP_pair", "TFscan_pair","PPI_pair", "repeat_pair"),
                     triplx = c("triplex"), 
                     Chrom = c('HISTONE', "Chromatin_pair"),
                     Meth = c('Methylation__WGBS'),
                     AX = c('_AXSB'), 
                     ChIP = c('ChIP', "ChIP_pair"),
                     trnsp = c('TRNSCRP'))

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/Partition_6_CVfirst_dataset_Firre.RData")
aa_colnn <- my_name_dic[1:(nrow(my_name_dic) - 1),1]  #c(colnames(Partition_6_feature_mat_tile_2), colnames(Partition_6_feature_mat_pair_2))
my_feat_family_list_all <- list()
# generate the set of family combinations that want to test
# first sequence families, then add epigenome, chip and expression on top of it
# all combinations (with at least one family) containing members of distance kmerFreq motifScanRepPair triplex
aa_alph <- c(1,2,3,4,5, 6)
aa1_list <- list()
for(i in 1:length(aa_alph)){
  aatmp <- generate.kmers(.k = i, .alphabet = aa_alph)
  aatmp2 <- aatmp[unlist(lapply(strsplit(aatmp, ""), function(x) sum(duplicated(x)) == 0))]
  aatmp3 <- lapply(strsplit(aatmp2, ""), sort)
  aatmp2 <- aatmp2[!duplicated(aatmp3)]
  # remove ones with the same element but only different permutation
  
  aa1_list[[i]] <- lapply(strsplit(aatmp2, ""), as.numeric)
}
aa1_list_all <- unlist(aa1_list, recursive = F)
aa1_list_all_names <- character(length = length(aa1_list_all))
for(i in 1:length(aa1_list_all)){
  aa1_list_all_names[i] <- paste(names(all_families)[aa1_list_all[[i]]], collapse = "_")
}
names(aa1_list_all) <- aa1_list_all_names



for(i in 1:length(aa1_list_all)){
  my_feat_family_list_all[[i]] <- list()
  aanm <- character(0)
  for(j in 1:length(aa1_list_all[[i]])){
    aatmp <- numeric(0)
    aanm <- c(aanm, names(all_families)[[aa1_list_all[[i]][j]]])
    for(k in 1:length(all_families[[aa1_list_all[[i]][j]]])){
      aatmp <- c(aatmp, grep(pattern = all_families[[aa1_list_all[[i]][j]]][k], x = aa_colnn))
      
    }
    my_feat_family_list_all[[i]][[j]] <- aatmp
  }
  names(my_feat_family_list_all[[i]]) <- aanm
}
names(my_feat_family_list_all) <- aa1_list_all_names
#my_feat_family_list_all
# all combinations containing union of distance kmerFreq motifScanRepPair triplex as one family (sequence family), plus at least one member of the following families: epigenome ChIP transcription
all_families2 <- list(sequ = c('__DiFreq',"__5mer","__RBPscan", "__TFscan","__rep","RBP_pair", "TFscan_pair","PPI_pair", "repeat_pair","triplex"),
                      dist = c('Genomic_distance'),
                      Chrom = c('HISTONE', "Chromatin_pair"),
                      Meth = c('Methylation__WGBS'),
                      AX = c('_AXSB'), 
                      ChIP = c('ChIP', "ChIP_pair"),
                      trnsp = c('TRNSCRP'))


aa_alph2 <- c(2,3,4,5,6,7)
aa1_list2 <- list()
for(i in 1:length(aa_alph2)){
  aatmp <- generate.kmers(.seq = "1",.k = i, .alphabet = aa_alph2)
  aatmp2 <- aatmp[unlist(lapply(strsplit(aatmp, ""), function(x) sum(duplicated(x)) == 0))]
  aatmp3 <- lapply(strsplit(aatmp2, ""), sort)
  aatmp2 <- aatmp2[!duplicated(aatmp3)]
  # remove ones with the same element but only different permutation
  
  aa1_list2[[i]] <- lapply(strsplit(aatmp2, ""), as.numeric)
}
aa1_list_all2 <- unlist(aa1_list2, recursive = F)
aa1_list_all_names2 <- character(length = length(aa1_list_all2))
for(i in 1:length(aa1_list_all2)){
  aa1_list_all_names2[i] <- paste(names(all_families2)[aa1_list_all2[[i]]], collapse = "_")
}
names(aa1_list_all2) <- aa1_list_all_names2



aac <- length(my_feat_family_list_all)
for(i in 1:length(aa1_list_all2)){
  my_feat_family_list_all[[aac + i]] <- list()
  aanm <- character(0)
  for(j in 1:length(aa1_list_all2[[i]])){
    aatmp <- numeric(0)
    aanm <- c(aanm, names(all_families2)[[aa1_list_all2[[i]][j]]])
    for(k in 1:length(all_families2[[aa1_list_all2[[i]][j]]])){
      aatmp <- c(aatmp, grep(pattern = all_families2[[aa1_list_all2[[i]][j]]][k], x = aa_colnn))
      
    }
    my_feat_family_list_all[[aac + i]][[j]] <- aatmp
  }
  names(my_feat_family_list_all[[aac + i]]) <- aanm
}
names(my_feat_family_list_all)[(aac+1):length(my_feat_family_list_all)] <- aa1_list_all_names2

# make the other two required vectors for each exp: my_family_transform my_family_trans_len
aaal_fam <- union(names(all_families), names(all_families2))
aa_logical_len_dic <- data.frame(family =aaal_fam,
                                 transf= c(F,T,T,T,T,T,T,F,T,T,T,T),
                                 ful_len=numeric(length(aaal_fam)),
                                 tran_len = numeric(length(aaal_fam)))
View(aa_logical_len_dic)
all_families_both <- list(dist = c('Genomic_distance'),
                          kmerFq = c('__DiFreq',"__5mer"),
                          motSc = c("__RBPscan", "__TFscan"),
                          Rep = c( "__rep"),
                          pairs = c("RBP_pair", "TFscan_pair","PPI_pair", "repeat_pair"),
                          triplx = c("triplex"), 
                          sequ = c('__DiFreq',"__5mer","__RBPscan", "__TFscan","__rep","RBP_pair", "TFscan_pair","PPI_pair", "repeat_pair","triplex"),
                          Chrom = c('HISTONE', "Chromatin_pair"),
                          Meth = c('Methylation__WGBS'),
                          AX = c('_AXSB'), 
                          ChIP = c('ChIP', "ChIP_pair"),
                          trnsp = c('TRNSCRP'))


all_families_both_ind <- list()
for(i in 1:length(all_families_both)){
  tmp <- numeric(0)
  for(j in 1:length(all_families_both[[i]])){
    tmp <- c(tmp, grep(pattern = all_families_both[[i]][j], x = aa_colnn))
  }
  all_families_both_ind[[i]] <- tmp
}
names(all_families_both_ind) <- names(all_families_both)
aassf <- unlist(lapply(all_families_both_ind, length))
aa_logical_len_dic$ful_len <- aassf[match(aa_logical_len_dic$family, names(aassf))]
aassf2 <- unlist(lapply(aassf, function(x) min(x, ceiling(2*sqrt(x)))))
aa_logical_len_dic$tran_len <- aassf2[match(aa_logical_len_dic$family, names(aassf2))]


length(my_feat_family_list_all)
sample(names(my_feat_family_list_all), 20)


aadress <- "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Input1/"
for(i in 1:length(my_feat_family_list_all)){
  aa_fam <- unlist(strsplit(names(my_feat_family_list_all)[i], "_"))
  my_family_trans_len <- aa_logical_len_dic$tran_len[match(aa_fam, aa_logical_len_dic$family)] 
  my_family_transform <- aa_logical_len_dic$transf[match(aa_fam, aa_logical_len_dic$family)] 
  my_feat_family_list <- my_feat_family_list_all[[i]]
  
  aafile <- paste0(aadress, "famSVD_inp1_",names(my_feat_family_list_all)[i], ".RData")
  save(list = c("my_family_trans_len", "my_family_transform", "my_feat_family_list"), file = aafile)
}

my_family_trans_len
my_feat_family_list
my_family_transform

aa_un_own <- unique(Partition_6_random_chunk_cv_df$owner)
aa_all_cv <- c(paste0("RCV", c(1:5)), paste0("CCV", c(1:5)))
for(i in 1:length(my_feat_family_list_all)){
  aamyfilt <- " no_filter"
  for(j in 1:length(aa_un_own)){
    for(aacurcv in 1:length(aa_all_cv)){
      cat(c("Rscript --vanilla RF_load_run_p6_CV_svd_feature.R",
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_CVfirst_dataset_",
                   aa_un_own[j],".RData "), 
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_svd1_Results/Learned_models/",aa_un_own[j],"/"),
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_svd1_Results/performance_plots/",aa_un_own[j],"/"),
            aa_all_cv[aacurcv],
            
            sample(c(100:10000), 1),
            names(my_feat_family_list_all)[i], 
            aamyfilt,
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_svd1_input/famSVD_inp1_",names(my_feat_family_list_all)[i],".RData"),
            "\n"),
          sep = " ", file = "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_CV_svd1.job", append = T)
    }
    
  }
  
}












#########################################################################################################
# write job files for RF analysis while decreasing the resolution of distance --> to avoid overfiting in models


# should define three variables for each scenario: my_feat_family_list my_family_transform my_family_trans_len

# families:
all_families <- list(dist = c('Genomic_distance'),
                     kmerFq = c('__DiFreq',"__5mer"),
                     motSc = c("__RBPscan", "__TFscan"),
                     Rep = c( "__rep"),
                     pairs = c("RBP_pair", "TFscan_pair","PPI_pair", "repeat_pair"),
                     triplx = c("triplex"), 
                     Chrom = c('HISTONE', "Chromatin_pair"),
                     Meth = c('Methylation__WGBS'),
                     AX = c('_AXSB'), 
                     ChIP = c('ChIP', "ChIP_pair"),
                     trnsp = c('TRNSCRP'))

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/Partition_6_CVfirst_dataset_Firre.RData")
aa_colnn <- my_name_dic[1:(nrow(my_name_dic) - 1),1]  #c(colnames(Partition_6_feature_mat_tile_2), colnames(Partition_6_feature_mat_pair_2))
my_feat_family_list_all_2 <- list()
# generate the set of family combinations that want to test
# first sequence families, then add epigenome, chip and expression on top of it
# all combinations (with at least one family) containing members of distance kmerFreq motifScanRepPair triplex
aa_alph <- c(1,2,3,4,5, 6)
aa1_list <- list()
for(i in 1:length(aa_alph)){
  aatmp <- generate.kmers(.k = i, .alphabet = aa_alph)
  aatmp2 <- aatmp[unlist(lapply(strsplit(aatmp, ""), function(x) sum(duplicated(x)) == 0))]
  aatmp3 <- lapply(strsplit(aatmp2, ""), sort)
  aatmp2 <- aatmp2[!duplicated(aatmp3)]
  # remove ones with the same element but only different permutation
  
  aa1_list[[i]] <- lapply(strsplit(aatmp2, ""), as.numeric)
}
aa1_list_all <- unlist(aa1_list, recursive = F)
aa1_list_all_names <- character(length = length(aa1_list_all))
for(i in 1:length(aa1_list_all)){
  aa1_list_all_names[i] <- paste(names(all_families)[aa1_list_all[[i]]], collapse = "_")
}
names(aa1_list_all) <- aa1_list_all_names



for(i in 1:length(aa1_list_all)){
  my_feat_family_list_all_2[[i]] <- list()
  aanm <- character(0)
  for(j in 1:length(aa1_list_all[[i]])){
    aatmp <- numeric(0)
    aanm <- c(aanm, names(all_families)[[aa1_list_all[[i]][j]]])
    for(k in 1:length(all_families[[aa1_list_all[[i]][j]]])){
      aatmp <- c(aatmp, grep(pattern = all_families[[aa1_list_all[[i]][j]]][k], x = aa_colnn))
      
    }
    my_feat_family_list_all_2[[i]][[j]] <- aatmp
  }
  names(my_feat_family_list_all_2[[i]]) <- aanm
}
names(my_feat_family_list_all_2) <- aa1_list_all_names
#my_feat_family_list_all
# all combinations containing union of distance kmerFreq motifScanRepPair triplex as one family (sequence family), plus at least one member of the following families: epigenome ChIP transcription
all_families2 <- list(sequ = c('__DiFreq',"__5mer","__RBPscan", "__TFscan","__rep","RBP_pair", "TFscan_pair","PPI_pair", "repeat_pair","triplex"),
                      dist = c('Genomic_distance'),
                      Chrom = c('HISTONE', "Chromatin_pair"),
                      Meth = c('Methylation__WGBS'),
                      AX = c('_AXSB'), 
                      ChIP = c('ChIP', "ChIP_pair"),
                      trnsp = c('TRNSCRP'))


aa_alph2 <- c(3,4,5,6,7)
aa1_list2 <- list()
for(i in 1:length(aa_alph2)){
  aatmp <- generate.kmers(.seq = "12",.k = i, .alphabet = aa_alph2)
  aatmp2 <- aatmp[unlist(lapply(strsplit(aatmp, ""), function(x) sum(duplicated(x)) == 0))]
  aatmp3 <- lapply(strsplit(aatmp2, ""), sort)
  aatmp2 <- aatmp2[!duplicated(aatmp3)]
  # remove ones with the same element but only different permutation
  
  aa1_list2[[i]] <- lapply(strsplit(aatmp2, ""), as.numeric)
}
aa1_list_all2 <- unlist(aa1_list2, recursive = F)
aa1_list_all_names2 <- character(length = length(aa1_list_all2))
for(i in 1:length(aa1_list_all2)){
  aa1_list_all_names2[i] <- paste(names(all_families2)[aa1_list_all2[[i]]], collapse = "_")
}
names(aa1_list_all2) <- aa1_list_all_names2

aa1_list_all3 <- list()
aa1_list_all3[2:(length(aa1_list_all2) + 1)] <-  aa1_list_all2
aa1_list_all3[[1]] <- c(1,2)
names(aa1_list_all3)[1] <- "sequ_dist"
names(aa1_list_all3)[2:(length(aa1_list_all2) + 1)] <- names(aa1_list_all2)

aac <- length(my_feat_family_list_all_2)
for(i in 1:length(aa1_list_all3)){
  my_feat_family_list_all_2[[aac + i]] <- list()
  aanm <- character(0)
  for(j in 1:length(aa1_list_all3[[i]])){
    aatmp <- numeric(0)
    aanm <- c(aanm, names(all_families2)[[aa1_list_all3[[i]][j]]])
    for(k in 1:length(all_families2[[aa1_list_all3[[i]][j]]])){
      aatmp <- c(aatmp, grep(pattern = all_families2[[aa1_list_all3[[i]][j]]][k], x = aa_colnn))
      
    }
    my_feat_family_list_all_2[[aac + i]][[j]] <- aatmp
  }
  names(my_feat_family_list_all_2[[aac + i]]) <- aanm
}
names(my_feat_family_list_all_2)[(aac+1):length(my_feat_family_list_all_2)] <- names(aa1_list_all3)

# make the other two required vectors for each exp: my_family_transform my_family_trans_len
aaal_fam <- union(names(all_families), names(all_families2))
aa_logical_len_dic <- data.frame(family =aaal_fam,
                                 transf= c(F,F,F,F,F,F,F,F,F,F,F,T),
                                 ful_len=numeric(length(aaal_fam)),
                                 tran_len = numeric(length(aaal_fam)))
View(aa_logical_len_dic)
all_families_both <- list(dist = c('Genomic_distance'),
                          kmerFq = c('__DiFreq',"__5mer"),
                          motSc = c("__RBPscan", "__TFscan"),
                          Rep = c( "__rep"),
                          pairs = c("RBP_pair", "TFscan_pair","PPI_pair", "repeat_pair"),
                          triplx = c("triplex"), 
                          sequ = c('__DiFreq',"__5mer","__RBPscan", "__TFscan","__rep","RBP_pair", "TFscan_pair","PPI_pair", "repeat_pair","triplex"),
                          Chrom = c('HISTONE', "Chromatin_pair"),
                          Meth = c('Methylation__WGBS'),
                          AX = c('_AXSB'), 
                          ChIP = c('ChIP', "ChIP_pair"),
                          trnsp = c('TRNSCRP'))


all_families_both_ind <- list()
for(i in 1:length(all_families_both)){
  tmp <- numeric(0)
  for(j in 1:length(all_families_both[[i]])){
    tmp <- c(tmp, grep(pattern = all_families_both[[i]][j], x = aa_colnn))
  }
  all_families_both_ind[[i]] <- tmp
}
names(all_families_both_ind) <- names(all_families_both)
aassf <- unlist(lapply(all_families_both_ind, length))
aa_logical_len_dic$ful_len <- aassf[match(aa_logical_len_dic$family, names(aassf))]
aassf2 <- unlist(lapply(aassf, function(x) min(x, ceiling(2*sqrt(x)))))
aa_logical_len_dic$tran_len <- aassf2[match(aa_logical_len_dic$family, names(aassf2))]


length(my_feat_family_list_all_2)
sample(names(my_feat_family_list_all_2), 20)


aadress <- "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Input2/"
for(i in 1:length(my_feat_family_list_all_2)){
  aa_fam <- unlist(strsplit(names(my_feat_family_list_all_2)[i], "_"))
  my_family_trans_len <- aa_logical_len_dic$tran_len[match(aa_fam, aa_logical_len_dic$family)] 
  my_family_transform <- aa_logical_len_dic$transf[match(aa_fam, aa_logical_len_dic$family)] 
  my_feat_family_list <- my_feat_family_list_all_2[[i]]
  
  aafile <- paste0(aadress, "famSVD_inp2_",names(my_feat_family_list_all_2)[i], ".RData")
  save(list = c("my_family_trans_len", "my_family_transform", "my_feat_family_list"), file = aafile)
}

my_family_trans_len
my_feat_family_list
my_family_transform

aa_un_own <- unique(Partition_6_random_chunk_cv_df$owner)
aa_all_cv <- c(paste0("RCV", c(1:5)), paste0("CCV", c(1:5)))
for(i in 1:length(my_feat_family_list_all_2)){
  aamyfilt <- " no_filter"
  for(j in 1:length(aa_un_own)){
    for(aacurcv in 1:length(aa_all_cv)){
      cat(c("Rscript --vanilla RF_load_run_p6_CV_svd2_lowres_feature.R",
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_CVfirst_dataset_",
                   aa_un_own[j],".RData "), 
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_svd2_lowres_Results/Learned_models/",aa_un_own[j],"/"),
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_svd2_lowres_Results/performance_plots/",aa_un_own[j],"/"),
            aa_all_cv[aacurcv],
            
            sample(c(100:10000), 1),
            names(my_feat_family_list_all_2)[i], 
            aamyfilt,
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_svd2_input/famSVD_inp2_",names(my_feat_family_list_all_2)[i],".RData"),
            "\n"),
          sep = " ", file = "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_CV_svd2_lowres.job", append = T)
    }
    
  }
  
}


nrow(Partition_6_random_chunk_cv_df)
length(distance_feature_partition_6)
all(Partition_6_random_chunk_cv_df$tile_name == names(distance_feature_partition_6))

aa_job1 <- readLines("~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_CV_svd2_lowres.job")
aa_job1Fil <- aa_job1[grep("dist", aa_job1)]
aa_job1Fil2 <- gsub(pattern = "RF_load_run_p6_CV_svd2_lowres_feature.R", replacement = "RF_load_run_p6_CV_svd2_lowres10MG_feature.R", x = aa_job1Fil)
aa_job1Fil3 <- gsub(pattern = "Partition_6_CV_svd2_lowres_Results", replacement = "Partition_6_CV_svd2_lowres10MG_Results", x = aa_job1Fil2)

writeLines(aa_job1Fil3, "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_CV_svd2_lowres10MG.job")

################# write jobs for the updated cv partitioning (fixed same)
aa_job1 <- readLines("~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_CV_svd2_lowres.job")
#aa_job1Fil <- aa_job1[grep("dist", aa_job1)]
aa_job1Fil2 <- gsub(pattern = "RF_load_run_p6_CV_svd2_lowres_feature.R", replacement = "RF_load_run_p6_CV2_svd_lowres1MG_feature.R", x = aa_job1)
aa_job1Fil3 <- gsub(pattern = "Partition_6_CV_svd2_lowres_Results", replacement = "Partition_6_CV2_svd_lowres1MG", x = aa_job1Fil2)

writeLines(aa_job1Fil3, "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_CV2_svd_lowres1MG.job")


################# write jobs for the updated cv partitioning (fixed same) --> fixed models not to use all features
aa_job1 <- readLines("~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_CV2_svd_lowres1MG.job")
aa_job1Fil <- aa_job1[-grep("sequ_dist", aa_job1)]
#aa_job1Fil2 <- gsub(pattern = "RF_load_run_p6_CV_svd2_lowres_feature.R", replacement = "RF_load_run_p6_CV2_svd_lowres1MG_feature.R", x = aa_job1Fil)
#aa_job1Fil3 <- gsub(pattern = "Partition_6_CV_svd2_lowres_Results", replacement = "Partition_6_CV2_svd_lowres1MG", x = aa_job1Fil2)

writeLines(aa_job1Fil, "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_CV2_svd_lowres1MG_fix.job")

################# write jobs with sequence SVD but not distance
#my_feat_family_list_all
# all combinations containing union of distance kmerFreq motifScanRepPair triplex as one family (sequence family), plus at least one member of the following families: epigenome ChIP transcription
all_families2 <- list(sequ = c('__DiFreq',"__5mer","__RBPscan", "__TFscan","__rep","RBP_pair", "TFscan_pair","PPI_pair", "repeat_pair","triplex"),
                      dist = c('Genomic_distance'),
                      Chrom = c('HISTONE', "Chromatin_pair"),
                      Meth = c('Methylation__WGBS'),
                      AX = c('_AXSB'), 
                      ChIP = c('ChIP', "ChIP_pair"),
                      trnsp = c('TRNSCRP'))


aa_alph2 <- c(3,4,5,6,7)
aa1_list4 <- list()
for(i in 1:length(aa_alph2)){
  aatmp <- generate.kmers(.seq = "1",.k = i, .alphabet = aa_alph2)
  aatmp2 <- aatmp[unlist(lapply(strsplit(aatmp, ""), function(x) sum(duplicated(x)) == 0))]
  aatmp3 <- lapply(strsplit(aatmp2, ""), sort)
  aatmp2 <- aatmp2[!duplicated(aatmp3)]
  # remove ones with the same element but only different permutation
  
  aa1_list4[[i]] <- lapply(strsplit(aatmp2, ""), as.numeric)
}
aa1_list_all4 <- unlist(aa1_list4, recursive = F)
aa1_list_all_names4 <- character(length = length(aa1_list_all4))
for(i in 1:length(aa1_list_all4)){
  aa1_list_all_names4[i] <- paste(names(all_families2)[aa1_list_all4[[i]]], collapse = "_")
}
names(aa1_list_all4) <- aa1_list_all_names4

aa1_list_all5 <- list()
aa1_list_all5[2:(length(aa1_list_all4) + 1)] <-  aa1_list_all4
aa1_list_all5[[1]] <- c(1)
names(aa1_list_all5)[1] <- "sequ"
names(aa1_list_all5)[2:(length(aa1_list_all4) + 1)] <- names(aa1_list_all4)

my_feat_family_list_all_4 <- my_feat_family_list_all_2
aac <- length(my_feat_family_list_all_4)
for(i in 1:length(aa1_list_all5)){
  my_feat_family_list_all_4[[aac + i]] <- list()
  aanm <- character(0)
  for(j in 1:length(aa1_list_all5[[i]])){
    aatmp <- numeric(0)
    aanm <- c(aanm, names(all_families2)[[aa1_list_all5[[i]][j]]])
    for(k in 1:length(all_families2[[aa1_list_all5[[i]][j]]])){
      aatmp <- c(aatmp, grep(pattern = all_families2[[aa1_list_all5[[i]][j]]][k], x = aa_colnn))
      
    }
    my_feat_family_list_all_4[[aac + i]][[j]] <- aatmp
  }
  names(my_feat_family_list_all_4[[aac + i]]) <- aanm
}
names(my_feat_family_list_all_4)[(aac+1):length(my_feat_family_list_all_4)] <- names(aa1_list_all5)
names(my_feat_family_list_all_4)


################# write jobs with transcriptional filtering for any experiment that uses RBP
# index of RBP features in the dataset
my_feat_family_list_all_5 <- my_feat_family_list_all_4
aa_filter_list3 <- list()
aass1 <- grep(pattern = "__RBPscan", x = aa_colnn)
aass2 <- grep(pattern = "RBP_pair", x = aa_colnn)
aa_filter_list3[[1]] <- union(aass1, aass2)

aamyi1 <- grep(pattern = "motSc",x = names(my_feat_family_list_all_4))
aamyi2 <- grep(pattern = "pairs",x = names(my_feat_family_list_all_4))
aamyi3 <- grep(pattern = "sequ",x = names(my_feat_family_list_all_4))
aamyi4 <- grep(pattern = "trnsp",x = names(my_feat_family_list_all_4))


aamyi_all1 <- sort(union(union(aamyi1, aamyi2), intersect(aamyi3,aamyi4)))
aackl <- length(my_feat_family_list_all_4)
my_feat_family_list_all_5_filter <- list()
for(i in 1:length(aamyi_all1)){
  my_feat_family_list_all_5[[aackl + i]] <- my_feat_family_list_all_4[[aamyi_all1[i]]]
  names(my_feat_family_list_all_5)[aackl + i] <- paste0(names(my_feat_family_list_all_4)[aamyi_all1[i]], "_TXRB")
  my_feat_family_list_all_5_filter[[aackl + i]] <- aa_filter_list3[[1]] 
}


################# write jobs with transcriptional filtering for any experiment that uses repeats
# index of repeat features in the dataset
aass1r <- grep(pattern = "__rep", x = aa_colnn)
aass2r <- grep(pattern = "repeat_pair", x = aa_colnn)
aa_filter_list3[[2]] <- union(aass1r, aass2r)

write.table(aa_filter_list3[[1]], quote = F, row.names = F, col.names = F, 
            file = paste0("~/Documents/Shayan/BioInf/lncRNA/partition_6_CV2_filt/my_filter_RBP.txt"), append = F)
write.table(aa_filter_list3[[2]], quote = F, row.names = F, col.names = F, 
            file = paste0("~/Documents/Shayan/BioInf/lncRNA/partition_6_CV2_filt/my_filter_repeat.txt"), append = F)

aamyi1r <- grep(pattern = "Rep",x = names(my_feat_family_list_all_4))
aamyi2r <- grep(pattern = "pairs",x = names(my_feat_family_list_all_4))
aamyi3r <- grep(pattern = "sequ",x = names(my_feat_family_list_all_4))
aamyi4r <- grep(pattern = "trnsp",x = names(my_feat_family_list_all_4))

aamyi_all1r <- sort(union(union(aamyi1r, aamyi2r), intersect(aamyi3r, aamyi4r)))


aackl <- length(my_feat_family_list_all_5)

for(i in 1:length(aamyi_all1r)){
  my_feat_family_list_all_5[[aackl + i]] <- my_feat_family_list_all_4[[aamyi_all1r[i]]]
  names(my_feat_family_list_all_5)[aackl + i] <- paste0(names(my_feat_family_list_all_4)[aamyi_all1r[i]], "_TXRE")
  my_feat_family_list_all_5_filter[[aackl + i]] <- aa_filter_list3[[2]] 
}

# if(length(aa_filter_list2[[j]]) == 0){
#   aamyfilt <- " no_filter"
# }else{
#   write.table(aa_filter_list2[[j]], quote = F, row.names = F, col.names = F, 
#               file = paste0("~/Documents/Shayan/BioInf/lncRNA/partition_6_CV2_filt/my_filter_", j,".txt"), append = F)
#   aamyfilt <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_Dec_TOBE_exp_FILTERED2/my_filter_", j,".txt")
# }


aadress <- "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Input2/"
for(i in 96:length(my_feat_family_list_all_5)){
  aa_fam <- unlist(strsplit(names(my_feat_family_list_all_5)[i], "_"))
  aa_fam <- setdiff(aa_fam, c("TXRE","TXRB"))
  my_family_trans_len <- aa_logical_len_dic$tran_len[match(aa_fam, aa_logical_len_dic$family)] 
  my_family_transform <- aa_logical_len_dic$transf[match(aa_fam, aa_logical_len_dic$family)] 
  my_feat_family_list <- my_feat_family_list_all_5[[i]]
  
  aafile <- paste0(aadress, "famSVD_inp2_",names(my_feat_family_list_all_5)[i], ".RData")
  save(list = c("my_family_trans_len", "my_family_transform", "my_feat_family_list"), file = aafile)
}

# my_family_trans_len
# my_feat_family_list
# my_family_transform
# 
aa_un_own <- unique(Partition_6_random_chunk_cv_df$owner)
aa_all_cv <- c(paste0("RCV", c(1:5)), paste0("CCV", c(1:5)))
for(i in 96:length(my_feat_family_list_all_5)){
  if(length(my_feat_family_list_all_5_filter[[i]]) == 0){
    aamyfilt <- " no_filter"
  }else if(length(my_feat_family_list_all_5_filter[[i]]) == 80){
    aamyfilt <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_CV2_filt/my_filter_RBP.txt")
  }else if(length(my_feat_family_list_all_5_filter[[i]]) == 60){
    aamyfilt <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_CV2_filt/my_filter_repeat.txt")
  }
  
  for(j in 1:length(aa_un_own)){
    for(aacurcv in 1:length(aa_all_cv)){
      cat(c("Rscript --vanilla RF_load_run_p6_CV2_svd_lowres1MG_feature.R",
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_CVfirst_dataset_",
                   aa_un_own[j],".RData "), 
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_lowres1MG/Learned_models/",aa_un_own[j],"/"),
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_lowres1MG/performance_plots/",aa_un_own[j],"/"),
            aa_all_cv[aacurcv],
            
            sample(c(100:10000), 1),
            names(my_feat_family_list_all_5)[i], 
            aamyfilt,
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_svd2_input/famSVD_inp2_",names(my_feat_family_list_all_5)[i],".RData"),
            "\n"),
          sep = " ", file = "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_CV2_svdTX_lowres1MG.job", append = T)
    }
    
  }
  
}



# write jobs for running the random forests for top200 features, with no_filter, RBP_filter, repeat_filter, or both
aa_filter_list4 <- list()

aass1 <- grep(pattern = "__RBPscan", x = aa_colnn)
aass2 <- grep(pattern = "RBP_pair", x = aa_colnn)
aa_filter_list4[[1]] <- union(aass1, aass2)

aass1r <- grep(pattern = "__rep", x = aa_colnn)
aass2r <- grep(pattern = "repeat_pair", x = aa_colnn)
aa_filter_list4[[2]] <- union(aass1r, aass2r)

aa_filter_list4[[3]] <- union(aa_filter_list4[[1]], aa_filter_list4[[2]])

write.table(aa_filter_list4[[1]], quote = F, row.names = F, col.names = F, 
            file = paste0("~/Documents/Shayan/BioInf/lncRNA/partition_6_CV2_filt/my_filter_RBP.txt"), append = F)
write.table(aa_filter_list4[[2]], quote = F, row.names = F, col.names = F, 
            file = paste0("~/Documents/Shayan/BioInf/lncRNA/partition_6_CV2_filt/my_filter_repeat.txt"), append = F)
write.table(aa_filter_list4[[3]], quote = F, row.names = F, col.names = F, 
            file = paste0("~/Documents/Shayan/BioInf/lncRNA/partition_6_CV2_filt/my_filter_both.txt"), append = F)

aa_filt_dic <- c("no_filter",
                 "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_CV2_filt/my_filter_RBP.txt",
                 "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_CV2_filt/my_filter_repeat.txt",
                 "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_CV2_filt/my_filter_both.txt")
aa_filt_dic_names <- c("noFilter", "RBPFilter", "RepFilter","BothFilter")

aa_un_own <- unique(Partition_6_random_chunk_cv_df$owner)
aa_all_cv <- c(paste0("RCV", c(1:5)), paste0("CCV", c(1:5)))
for(i in 1:length(aa_filt_dic)){
  aamyfilt <- aa_filt_dic[i]
  
  for(j in 1:length(aa_un_own)){
    for(aacurcv in 1:length(aa_all_cv)){
      cat(c("Rscript --vanilla RF_load_run_p6_CV2_200feat.R",
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_CVfirst_dataset_",
                   aa_un_own[j],".RData "), 
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_200feat/Learned_models/",aa_un_own[j],"/"),
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_200feat/performance_plots/",aa_un_own[j],"/"),
            aa_all_cv[aacurcv],
            
            sample(c(100:10000), 1),
            paste0("top200Features_", aa_filt_dic_names[i]), 
            aamyfilt,
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/top200_feature_index/index_",aa_un_own[j],"_",substr(aa_all_cv[aacurcv],1,1),".RData"),
            "\n"),
          sep = " ", file = "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_CV2_top200.job", append = T)
    }
    
  }
  
}



########################################################################################################################################################################
########################################################################################################################################################################
# write jobs for running all model combinations after fixing RBP and TFscan features (problem regarding duplicates) --> do not include distance in any of them


# families:
all_families <- list(dist = c('Genomic_distance'),
                     kmerFq = c('__DiFreq',"__5mer"),
                     motSc = c("__RBPscan", "__TFscan"),
                     Rep = c( "__rep"),
                     pairs = c("RBP_pair", "TFscan_pair","PPI_pair", "repeat_pair"),
                     triplx = c("triplex"), 
                     Chrom = c('HISTONE', "Chromatin_pair"),
                     Meth = c('Methylation__WGBS'),
                     AX = c('_AXSB'), 
                     ChIP = c('ChIP', "ChIP_pair"),
                     trnsp = c('TRNSCRP'))

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/Partition_6_CVfirst_dataset_Firre.RData")
aa_colnn <- my_name_dic[1:(nrow(my_name_dic) - 1),1]  #c(colnames(Partition_6_feature_mat_tile_2), colnames(Partition_6_feature_mat_pair_2))
my_feat_family_list_all_6 <- list()
# generate the set of family combinations that want to test
# first sequence families, then add epigenome, chip and expression on top of it
# all combinations (with at least one family) containing members of distance kmerFreq motifScanRepPair triplex
aa_alph <- c(2,3,4,5, 6)
aa1_list <- list()
for(i in 1:length(aa_alph)){
  aatmp <- generate.kmers(.k = i, .alphabet = aa_alph)
  aatmp2 <- aatmp[unlist(lapply(strsplit(aatmp, ""), function(x) sum(duplicated(x)) == 0))]
  aatmp3 <- lapply(strsplit(aatmp2, ""), sort)
  aatmp2 <- aatmp2[!duplicated(aatmp3)]
  # remove ones with the same element but only different permutation
  
  aa1_list[[i]] <- lapply(strsplit(aatmp2, ""), as.numeric)
}
aa1_list_all <- unlist(aa1_list, recursive = F)
aa1_list_all_names <- character(length = length(aa1_list_all))
for(i in 1:length(aa1_list_all)){
  aa1_list_all_names[i] <- paste(names(all_families)[aa1_list_all[[i]]], collapse = "_")
}
names(aa1_list_all) <- aa1_list_all_names



for(i in 1:length(aa1_list_all)){
  my_feat_family_list_all_6[[i]] <- list()
  aanm <- character(0)
  for(j in 1:length(aa1_list_all[[i]])){
    aatmp <- numeric(0)
    aanm <- c(aanm, names(all_families)[[aa1_list_all[[i]][j]]])
    for(k in 1:length(all_families[[aa1_list_all[[i]][j]]])){
      aatmp <- c(aatmp, grep(pattern = all_families[[aa1_list_all[[i]][j]]][k], x = aa_colnn))
      
    }
    my_feat_family_list_all_6[[i]][[j]] <- aatmp
  }
  names(my_feat_family_list_all_6[[i]]) <- aanm
}
names(my_feat_family_list_all_6) <- aa1_list_all_names



all_families2 <- list(sequ = c('__DiFreq',"__5mer","__RBPscan", "__TFscan","__rep","RBP_pair", "TFscan_pair","PPI_pair", "repeat_pair","triplex"),
                      dist = c('Genomic_distance'),
                      Chrom = c('HISTONE', "Chromatin_pair"),
                      Meth = c('Methylation__WGBS'),
                      AX = c('_AXSB'), 
                      ChIP = c('ChIP', "ChIP_pair"),
                      trnsp = c('TRNSCRP'))


aa_alph2 <- c(3,4,5,6,7)
aa1_list4 <- list()
for(i in 1:length(aa_alph2)){
  aatmp <- generate.kmers(.seq = "1",.k = i, .alphabet = aa_alph2)
  aatmp2 <- aatmp[unlist(lapply(strsplit(aatmp, ""), function(x) sum(duplicated(x)) == 0))]
  aatmp3 <- lapply(strsplit(aatmp2, ""), sort)
  aatmp2 <- aatmp2[!duplicated(aatmp3)]
  # remove ones with the same element but only different permutation
  
  aa1_list4[[i]] <- lapply(strsplit(aatmp2, ""), as.numeric)
}
aa1_list_all4 <- unlist(aa1_list4, recursive = F)
aa1_list_all_names4 <- character(length = length(aa1_list_all4))
for(i in 1:length(aa1_list_all4)){
  aa1_list_all_names4[i] <- paste(names(all_families2)[aa1_list_all4[[i]]], collapse = "_")
}
names(aa1_list_all4) <- aa1_list_all_names4

aa1_list_all5 <- list()
aa1_list_all5[2:(length(aa1_list_all4) + 1)] <-  aa1_list_all4
aa1_list_all5[[1]] <- c(1)
names(aa1_list_all5)[1] <- "sequ"
names(aa1_list_all5)[2:(length(aa1_list_all4) + 1)] <- names(aa1_list_all4)

#my_feat_family_list_all_4 <- my_feat_family_list_all_2
aac <- length(my_feat_family_list_all_6)
for(i in 1:length(aa1_list_all5)){
  my_feat_family_list_all_6[[aac + i]] <- list()
  aanm <- character(0)
  for(j in 1:length(aa1_list_all5[[i]])){
    aatmp <- numeric(0)
    aanm <- c(aanm, names(all_families2)[[aa1_list_all5[[i]][j]]])
    for(k in 1:length(all_families2[[aa1_list_all5[[i]][j]]])){
      aatmp <- c(aatmp, grep(pattern = all_families2[[aa1_list_all5[[i]][j]]][k], x = aa_colnn))
      
    }
    my_feat_family_list_all_6[[aac + i]][[j]] <- aatmp
  }
  names(my_feat_family_list_all_6[[aac + i]]) <- aanm
}
names(my_feat_family_list_all_6)[(aac+1):length(my_feat_family_list_all_6)] <- names(aa1_list_all5)
names(my_feat_family_list_all_6)

aaal_fam <- union(names(all_families), names(all_families2))
aa_logical_len_dic <- data.frame(family =aaal_fam,
                                 transf= c(F,F,F,F,F,F,F,F,F,F,F,T),
                                 ful_len=numeric(length(aaal_fam)),
                                 tran_len = numeric(length(aaal_fam)))
all_families_both <- list(dist = c('Genomic_distance'),
                          kmerFq = c('__DiFreq',"__5mer"),
                          motSc = c("__RBPscan", "__TFscan"),
                          Rep = c( "__rep"),
                          pairs = c("RBP_pair", "TFscan_pair","PPI_pair", "repeat_pair"),
                          triplx = c("triplex"), 
                          sequ = c('__DiFreq',"__5mer","__RBPscan", "__TFscan","__rep","RBP_pair", "TFscan_pair","PPI_pair", "repeat_pair","triplex"),
                          Chrom = c('HISTONE', "Chromatin_pair"),
                          Meth = c('Methylation__WGBS'),
                          AX = c('_AXSB'), 
                          ChIP = c('ChIP', "ChIP_pair"),
                          trnsp = c('TRNSCRP'))
load("Partition_6_CVfirst_dataset/Partition_6_CVfirst_dataset_2410003L11Rik.RData")
aa_colnn <- my_name_dic[1:(nrow(my_name_dic) - 1),1]  #c(colnames(Partition_6_feature_mat_tile_2), colnames(Partition_6_feature_mat_pair_2))
all_families_both_ind <- list()
for(i in 1:length(all_families_both)){
  tmp <- numeric(0)
  for(j in 1:length(all_families_both[[i]])){
    tmp <- c(tmp, grep(pattern = all_families_both[[i]][j], x = aa_colnn))
  }
  all_families_both_ind[[i]] <- tmp
}
names(all_families_both_ind) <- names(all_families_both)
# aassf <- unlist(lapply(all_families_both_ind, length))
# aa_logical_len_dic$ful_len <- aassf[match(aa_logical_len_dic$family, names(aassf))]
# aassf2 <- unlist(lapply(aassf, function(x) min(x, ceiling(2*sqrt(x)))))
# aa_logical_len_dic$tran_len <- aassf2[match(aa_logical_len_dic$family, names(aassf2))]

aassf <- unlist(lapply(all_families_both_ind, length))
aa_logical_len_dic$ful_len <- aassf[match(aa_logical_len_dic$family, names(aassf))]
aassf2 <- unlist(lapply(aassf, function(x) min(x, ceiling(4*sqrt(x)))))
aa_logical_len_dic$tran_len <- aassf2[match(aa_logical_len_dic$family, names(aassf2))]

### adding TX filtering jobs for RBP and Repeats
#my_feat_family_list_all_5 <- my_feat_family_list_all_6
aa_filter_list3 <- list()
aass1 <- grep(pattern = "__RBPscan", x = aa_colnn)
aass2 <- grep(pattern = "RBP_pair", x = aa_colnn)
aa_filter_list3[[1]] <- union(aass1, aass2)

aamyi1 <- grep(pattern = "motSc",x = names(my_feat_family_list_all_6))
aamyi2 <- grep(pattern = "pairs",x = names(my_feat_family_list_all_6))
aamyi3 <- grep(pattern = "sequ",x = names(my_feat_family_list_all_6))
aamyi4 <- grep(pattern = "trnsp",x = names(my_feat_family_list_all_6))


aamyi_all1 <- sort(union(union(aamyi1, aamyi2), intersect(aamyi3,aamyi4)))
aackl <- length(my_feat_family_list_all_6)
my_feat_family_list_all_6_filter <- list()
for(i in 1:length(aamyi_all1)){
  my_feat_family_list_all_6[[aackl + i]] <- my_feat_family_list_all_6[[aamyi_all1[i]]]
  names(my_feat_family_list_all_6)[aackl + i] <- paste0(names(my_feat_family_list_all_6)[aamyi_all1[i]], "_TXRB")
  my_feat_family_list_all_6_filter[[aackl + i]] <- aa_filter_list3[[1]] 
}
########
aass1r <- grep(pattern = "__rep", x = aa_colnn)
aass2r <- grep(pattern = "repeat_pair", x = aa_colnn)
aa_filter_list3[[2]] <- union(aass1r, aass2r)

write.table(aa_filter_list3[[1]], quote = F, row.names = F, col.names = F, 
            file = paste0("~/Documents/Shayan/BioInf/lncRNA/partition_6_CV2_filt/my_filter_RBP.txt"), append = F)
write.table(aa_filter_list3[[2]], quote = F, row.names = F, col.names = F, 
            file = paste0("~/Documents/Shayan/BioInf/lncRNA/partition_6_CV2_filt/my_filter_repeat.txt"), append = F)

aamyi1r <- grep(pattern = "Rep",x = names(my_feat_family_list_all_6)[1:63])
aamyi2r <- grep(pattern = "pairs",x = names(my_feat_family_list_all_6)[1:63])
aamyi3r <- grep(pattern = "sequ",x = names(my_feat_family_list_all_6)[1:63])
aamyi4r <- grep(pattern = "trnsp",x = names(my_feat_family_list_all_6)[1:63])

aamyi_all1r <- sort(union(union(aamyi1r, aamyi2r), intersect(aamyi3r, aamyi4r)))


aackl <- length(my_feat_family_list_all_6)

for(i in 1:length(aamyi_all1r)){
  my_feat_family_list_all_6[[aackl + i]] <- my_feat_family_list_all_6[[aamyi_all1r[i]]]
  names(my_feat_family_list_all_6)[aackl + i] <- paste0(names(my_feat_family_list_all_6)[aamyi_all1r[i]], "_TXRE")
  my_feat_family_list_all_6_filter[[aackl + i]] <- aa_filter_list3[[2]] 
}


########


aadress <- "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Input3/"
for(i in 1:length(my_feat_family_list_all_6)){
  aa_fam <- unlist(strsplit(names(my_feat_family_list_all_6)[i], "_"))
  aa_fam <- setdiff(aa_fam, c("TXRE","TXRB"))
  my_family_trans_len <- aa_logical_len_dic$tran_len[match(aa_fam, aa_logical_len_dic$family)] 
  my_family_transform <- aa_logical_len_dic$transf[match(aa_fam, aa_logical_len_dic$family)] 
  my_feat_family_list <- my_feat_family_list_all_6[[i]]
  
  aafile <- paste0(aadress, "famSVD_inp3_",names(my_feat_family_list_all_6)[i], ".RData")
  save(list = c("my_family_trans_len", "my_family_transform", "my_feat_family_list"), file = aafile)
}

# my_family_trans_len
# my_feat_family_list
# my_family_transform
# 
aa_un_own <- unique(Partition_6_random_chunk_cv_df$owner)
aa_all_cv <- c(paste0("RCV", c(1:5)), paste0("CCV", c(1:5)))
for(i in 1:length(my_feat_family_list_all_6)){
  if(length(my_feat_family_list_all_6_filter[[i]]) == 0){
    aamyfilt <- " no_filter"
  }else if(length(my_feat_family_list_all_6_filter[[i]]) == 80){
    aamyfilt <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_CV2_filt/my_filter_RBP.txt")
  }else if(length(my_feat_family_list_all_6_filter[[i]]) == 60){
    aamyfilt <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_CV2_filt/my_filter_repeat.txt")
  }
  
  for(j in 1:length(aa_un_own)){
    for(aacurcv in 1:length(aa_all_cv)){
      cat(c("Rscript --vanilla RF_load_run_p6_CV2_svd_NoDinuc_Nodist_feature.R",
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_CVfirst_dataset_",
                   aa_un_own[j],".RData "), 
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_NOD/Learned_models/",aa_un_own[j],"/"),
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_NOD/performance_plots/",aa_un_own[j],"/"),
            aa_all_cv[aacurcv],
            
            sample(c(100:10000), 1),
            names(my_feat_family_list_all_6)[i], 
            aamyfilt,
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_svd2_input/famSVD_inp3_",names(my_feat_family_list_all_6)[i],".RData"),
            "\n"),
          sep = " ", file = "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_CV2_svdTX_NoD.job", append = T)
    }
    
  }
  
}

#################### # rerun filtered jobs after fixing the filtered indecis
aa_job1 <- readLines("~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_CV2_svdTX_NoD.job")
aa_job1Fil <- aa_job1[grep("my_filter_", aa_job1)]
#aa_job1Fil2 <- gsub(pattern = "TX", replacement = "TXGRO", x = aa_job1Fil)
writeLines(aa_job1Fil, "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_CV2_svdTX_NoD_filtered.job")

##########
# write jobs for running the random forests for top200 features, with no_filter, RBP_filter, repeat_filter, or both
aa_filter_list4 <- list()
load("Partition_6_CVfirst_dataset/Partition_6_CVfirst_dataset_2410003L11Rik.RData")
aa_colnn <- my_name_dic[1:(nrow(my_name_dic) - 1),1]  #c(colnames(Partition_6_feature_mat_tile_2), colnames(Partition_6_feature_mat_pair_2))

aass1 <- grep(pattern = "__RBPscan", x = aa_colnn)
aass2 <- grep(pattern = "RBP_pair", x = aa_colnn)
aa_filter_list4[[1]] <- union(aass1, aass2)

aass1r <- grep(pattern = "__rep", x = aa_colnn)
aass2r <- grep(pattern = "repeat_pair", x = aa_colnn)
aa_filter_list4[[2]] <- union(aass1r, aass2r)

aa_filter_list4[[3]] <- union(aa_filter_list4[[1]], aa_filter_list4[[2]])

write.table(aa_filter_list4[[1]], quote = F, row.names = F, col.names = F, 
            file = paste0("~/Documents/Shayan/BioInf/lncRNA/partition_6_CV2_filt/my_filter_RBP.txt"), append = F)
write.table(aa_filter_list4[[2]], quote = F, row.names = F, col.names = F, 
            file = paste0("~/Documents/Shayan/BioInf/lncRNA/partition_6_CV2_filt/my_filter_repeat.txt"), append = F)
write.table(aa_filter_list4[[3]], quote = F, row.names = F, col.names = F, 
            file = paste0("~/Documents/Shayan/BioInf/lncRNA/partition_6_CV2_filt/my_filter_both.txt"), append = F)

aa_filt_dic <- c("no_filter",
                 "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_CV2_filt/my_filter_RBP.txt",
                 "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_CV2_filt/my_filter_repeat.txt",
                 "/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_CV2_filt/my_filter_both.txt")
aa_filt_dic_names <- c("noFilter", "RBPFilter", "RepFilter","BothFilter")

aa_un_own <- unique(Partition_6_random_chunk_cv_df$owner)
aa_all_cv <- c(paste0("RCV", c(1:5)), paste0("CCV", c(1:5)))
for(i in 1:length(aa_filt_dic)){
  aamyfilt <- aa_filt_dic[i]
  
  for(j in 1:length(aa_un_own)){
    for(aacurcv in 1:length(aa_all_cv)){
      cat(c("Rscript --vanilla RF_load_run_p6_CV2_200feat.R",
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_CVfirst_dataset_",
                   aa_un_own[j],".RData "), 
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_200feat_rerun1/Learned_models/",aa_un_own[j],"/"),
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_200feat_rerun1/performance_plots/",aa_un_own[j],"/"),
            aa_all_cv[aacurcv],
            
            sample(c(100:10000), 1),
            paste0("top200Features_", aa_filt_dic_names[i]), 
            aamyfilt,
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/top200_feature_index_rerun1/index_",aa_un_own[j],"_",substr(aa_all_cv[aacurcv],1,1),".RData"),
            "\n"),
          sep = " ", file = "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_CV2_top200_rerun1.job", append = T)
    }
    
  }
  
}






#######################################################################################################################################
# run RF models with individual families and all pairs of families to compute complementarit of information vs redundancy


# families:
all_families <- list( kmerFq = c("__5mer"),
                     motSc = c("__RBPscan", "__TFscan"),
                     Rep = c( "__rep"),
                     pairs = c("RBP_pair", "TFscan_pair","PPI_pair", "repeat_pair"),
                     triplx = c("triplex"), 
                     Chrom = c('HISTONE', "Chromatin_pair"),
                     Meth = c('Methylation__WGBS'),
                     AX = c('_AXSB'), 
                     ChIP = c('ChIP', "ChIP_pair"),
                     trnsp = c('TRNSCRP'))

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/Partition_6_CVfirst_dataset_Firre.RData")
aa_colnn <- my_name_dic[1:(nrow(my_name_dic) - 1),1]  #c(colnames(Partition_6_feature_mat_tile_2), colnames(Partition_6_feature_mat_pair_2))
library(tcR)
aak <- generate.kmers(.alphabet = paste0("_",c(1:10), "_"), .k = 2)

my_feat_family_list_all_pair <- list()
aanmn <- character(0)
aaccnt <- 0
for(i in 1:length(all_families)){
  
  for(j in i:length(all_families)){
    aaccnt <- aaccnt + 1
    my_feat_family_list_all_pair[[aaccnt]] <- grep(paste(unlist(all_families[union(i,j)]),collapse = "|"), aa_colnn)
    aanmn <- c(aanmn, paste(unlist(all_families[union(i,j)]),collapse = "_a_"))
  }
}
names(my_feat_family_list_all_pair) <- paste0("feat_",aanmn)

# aa_logical_len_dic <- data.frame(family =aaal_fam,
#                                  transf= c(F,F,F,F,F,F,F,F,F,F,F,T),
#                                  ful_len=numeric(length(aaal_fam)),
#                                  tran_len = numeric(length(aaal_fam)))
# 
# aassf <- unlist(lapply(all_families_both_ind, length))
# aa_logical_len_dic$ful_len <- aassf[match(aa_logical_len_dic$family, names(aassf))]
# aassf2 <- unlist(lapply(aassf, function(x) min(x, ceiling(4*sqrt(x)))))
# aa_logical_len_dic$tran_len <- aassf2[match(aa_logical_len_dic$family, names(aassf2))]


aadress <- "~/Documents/Shayan/BioInf/lncRNA/Family_svd/Input_pair/"
for(i in 1:length(my_feat_family_list_all_pair)){
  
  my_family_trans_len <- 10
  my_family_transform <- F
  my_feat_family_list <- my_feat_family_list_all_pair[i]
  
  aafile <- paste0(aadress, "famSVD_inpPair_",names(my_feat_family_list_all_pair)[i], ".RData")
  save(list = c("my_family_trans_len", "my_family_transform", "my_feat_family_list"), file = aafile)
}

aa_un_own <- unique(Partition_6_random_chunk_cv_df$owner)
aa_all_cv <- c(paste0("RCV", c(1:5)), paste0("CCV", c(1:5)))
for(i in 1:length(my_feat_family_list_all_pair)){
  aamyfilt <- " no_filter"
  for(j in 1:length(aa_un_own)){
    for(aacurcv in 1:length(aa_all_cv)){
      cat(c("Rscript --vanilla RF_load_run_p6_CV2_svd_NoDinuc_Nodist_feature.R",
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/input_RData/Partition_6_CVfirst_dataset_",
                   aa_un_own[j],".RData "), 
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_pair/Learned_models/",aa_un_own[j],"/"),
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_pair/performance_plots/",aa_un_own[j],"/"),
            aa_all_cv[aacurcv],
            
            sample(c(100:10000), 1),
            names(my_feat_family_list_all_pair)[i], 
            aamyfilt,
            paste0("/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV_pair_input/famSVD_inpPair_",names(my_feat_family_list_all_pair)[i],".RData"),
            "\n"),
          sep = " ", file = "~/Documents/Shayan/BioInf/lncRNA/run_rf_p6_CV_pair.job", append = T)
    }
    
  }
  
}

