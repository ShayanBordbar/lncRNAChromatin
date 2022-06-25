# lncRNA feature construction
#BiocManager::install("BSgenome.Mmusculus.UCSC.mm9")
library(BSgenome.Mmusculus.UCSC.mm9)



Distance_based_territory_assignment

MM9_1kb_tiled
MM9_1kb_tiled_owner
MM9_1kb_tiled_owner_labels
MM9_1kb_tiled_owner_labels_binary
MM9_1kb_tiled_TADnu


# lncRNA df
lncRNA_chosen_gt1k_uniqTiles

####################################################################################################
####################################################################################################
############################    create feaures representing chip tracks   ##########################
####################################################################################################
####################################################################################################
# CAGE
aaov <- findOverlaps(MM9_1kb_tiled_GR, mESC_CAGE_mm9_GR)
aadf <- data.frame(tile_index = aaov@from,
                   atac_index = aaov@to, 
                   CAGE_score = mESC_CAGE_mm9_GR$aa_max[aaov@to])
aadf_agg <- aggregate(aadf[c("CAGE_score")],
                      by = aadf[c("tile_index")],
                      FUN = max)

mESC_CAGE_mm9_tiled <- cbind(MM9_1kb_tiled[aadf_agg$tile_index,], aadf_agg$CAGE_score)
mESC_CAGE_mm9_tiledGR <- makeGRangesFromDataFrame(mESC_CAGE_mm9_tiled, keep.extra.columns = T)

mESC_CAGE_mm9_tiled$tile <- as.character(levels(mESC_CAGE_mm9_tiled$tile)[as.numeric(mESC_CAGE_mm9_tiled$tile)])
####################################################################################################
# CAGE for lncRNA

#mESC_CAGE_lncRNA_features <- numeric(length = length(lncRNA_chosen_gt1k_uniqTiles_GR))
#names(mESC_CAGE_lncRNA_features) <- lncRNA_chosen_gt1k_uniqTiles_GR$gene_name
aaov <- findOverlaps(lncRNA_chosen_gt1k_uniqTiles_GR, mESC_CAGE_mm9_GR)
aadf <- data.frame(lncRNA_index = aaov@from,
                   atac_index = aaov@to, 
                   CAGE_score = mESC_CAGE_mm9_GR$aa_max[aaov@to])
aadf_agg <- aggregate(aadf[c("CAGE_score")],
                      by = aadf[c("lncRNA_index")],
                      FUN = max)

aa_my_ind <- match(aadf_agg$lncRNA_index, c(1:length(lncRNA_chosen_gt1k_uniqTiles_GR)))
  
mESC_CAGE_lncRNA_features <- aadf_agg$CAGE_score
names(mESC_CAGE_lncRNA_features) <- lncRNA_chosen_gt1k_uniqTiles_GR$gene_name[aadf_agg$lncRNA_index]



####################################################################################################
####################################################################################################
# ATAC-seq
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
####################################################################################################
Partition_6_atac_Seq <- numeric(length = nrow(Partition_6_random_chunk_cv_df))
names(Partition_6_atac_Seq) <- Partition_6_random_chunk_cv_df$tile_name
aamESC_ATAC_Seq_tiled <- mESC_ATAC_Seq_tiled[mESC_ATAC_Seq_tiled$tile %in% names(Partition_6_atac_Seq), ]
Partition_6_atac_Seq[match(aamESC_ATAC_Seq_tiled$tile, names(Partition_6_atac_Seq))] <- aamESC_ATAC_Seq_tiled$`aadf_agg$atac_score`
save(list = c("Partition_6_atac_Seq"), file = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_atac_Seq.RData")

Partition_6_accessiblity <- data.frame(ATAC = Partition_6_atac_Seq, DNase = ChIPATLAS_features_partition6[,"DNase-Seq__46.bed"])
Partition_6_accessiblity[is.na(Partition_6_accessiblity)] <- 0
save(list = c("Partition_6_accessiblity"), file = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_accessiblity.RData")


####################################################################################################
####################################################################################################
#  chromatin marks:
#  "H3K27me3-mouse" "H3K9me3-mouse"  "H3K4me3-mouse"  "H3K4me1-mouse"  "H3K36me3-mouse"
#  "H3K9ac-mouse"   "H3K27ac-mouse" 

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
####################################################################################################
####################################################################################################
# data from Ren Lab:
# "GSM723015_RenLab-CTCF-mESC.bed"         "GSM723016_RenLab-H3K4me1-mESC.bed"     
# "GSM723017_RenLab-H3K4me3-mESC.bed"      "GSM723018_RenLab-P300-mESC.bed"        
# "GSM723019_RenLab-Pol2-mESC.bed"         "GSM723776_RenLab-RNA-Seq-mESC-ZY28.bed"

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
####################################################################################################
# RNA seq for lncRNA


aaRNAseqGR <- mESC_Renlab_tiled_list_GR$`GSM723776_RenLab-RNA-Seq-mESC-ZY28.bed`

aaov <- findOverlaps(lncRNA_chosen_gt1k_uniqTiles_GR, aaRNAseqGR)
aadf <- data.frame(lncRNA_index = aaov@from,
                   atac_index = aaov@to, 
                   RNA_score = aaRNAseqGR$`aadf_agg$score`[aaov@to])
aadf_agg <- aggregate(aadf[c("RNA_score")],
                      by = aadf[c("lncRNA_index")],
                      FUN = max)

mESC_RNAseq_lncRNA_features <- aadf_agg$RNA_score
names(mESC_RNAseq_lncRNA_features) <- lncRNA_chosen_gt1k_uniqTiles_GR$gene_name[aadf_agg$lncRNA_index]



####################################################################################################
####################################################################################################
# DNAse seq data
aaov <- findOverlaps(MM9_1kb_tiled_GR, mESC_DNAse_GR)
aadf <- data.frame(tile_index = aaov@from,
                   dnase_index = aaov@to, 
                   dnase_score = mESC_DNAse_GR$V7[aaov@to])
aadf_agg <- aggregate(aadf[c("dnase_score")],
                      by = aadf[c("tile_index")],
                      FUN = max)

mESC_DNAse_mm9_tiled <- cbind(MM9_1kb_tiled[aadf_agg$tile_index,], aadf_agg$dnase_score)
mESC_DNAse_mm9_tiledGR <- makeGRangesFromDataFrame(mESC_DNAse_mm9_tiled, keep.extra.columns = T)

####################################################################################################
####################################################################################################
# repeat elements
aaov <- findOverlaps(MM9_1kb_tiled_GR, MM9_repeats_GR)
aadf <- data.frame(tile_index = aaov@from,
                   repeat_index = aaov@to, 
                   repeat_name = MM9_repeats_GR$V5[aaov@to])


aatab <- table(MM9_repeats_GR$V4)
sort(aatab, decreasing = T)[1:20]
pie(sort(aatab, decreasing = T))
aadf_agg <- aggregate(aadf[c("repeat_name")],
                      by = aadf[c("tile_index")],
                      FUN = c)
aal <- unlist(lapply(aadf_agg$repeat_name, length))

mESC_repeat_mm9_tiled <- cbind(MM9_1kb_tiled[aadf_agg$tile_index,], aal)
mESC_repeat_mm9_tiledGR <- makeGRangesFromDataFrame(mESC_repeat_mm9_tiled, keep.extra.columns = T)
####################################################################################################
####################################################################################################

# Chip ATLAS:
aa_atlas_exp <- read.delim("~/Documents/Shayan/BioInf/lncRNA/ChipATLAS/experimentList.tab.txt", stringsAsFactors = F, header = F, sep = "\t")
aa_antigen <- read.delim("~/Documents/Shayan/BioInf/lncRNA/ChipATLAS/antigenList.tab", stringsAsFactors = F, header = T, sep = "\t")
aa_celltype <- read.delim("~/Documents/Shayan/BioInf/lncRNA/ChipATLAS/celltypeList.tab", stringsAsFactors = F, header = T, sep = "\t")
#aa_atlas_Files <- read.delim("~/Downloads/fileList.tab",stringsAsFactors = F)
#table(aa_atlas_Files$Adult.1)
table(aa_atlas_exp$V2)
# create a dataframe containing the meta data information of the Chip ATLAS, with three columns: SRX ID, cell type, antibody

table(aa_antigen$Genome)
table(aa_celltype$Genome)
#table(aa_atlas_exp$V1)
sum(aa_celltype$V3 %in% "")

View(aa_celltype[(aa_celltype$Cell_type_class %in% "Pluripotent stem cell") & (aa_celltype$Genome %in% "mm9"),])
aa_my_Exp_id <- aa_celltype$ID[(aa_celltype$Cell_type %in% "Embryonic Stem Cells")& (aa_celltype$Genome %in% "mm9")]
aa_my_Exp_id <- unlist(strsplit(aa_my_Exp_id, ","))

View(aa_atlas_exp[(aa_atlas_exp$V3 %in% c("Histone", "TFs and others", "RNA polymerase","DNase-seq")) & 
                    (aa_atlas_exp$V2 %in% "mm9") & (aa_atlas_exp$V6 %in% "Embryonic Stem Cells"),])

ChIP_ATLAS_exps <- aa_atlas_exp[(aa_atlas_exp$V3 %in% c("Histone", "TFs and others", "RNA polymerase","DNase-seq")) & 
                                  (aa_atlas_exp$V2 %in% "mm9") & (aa_atlas_exp$V6 %in% "Embryonic Stem Cells"),]
colnames(ChIP_ATLAS_exps)[1:6] <- c("ID","Genome","Antigen_class","Antigen", "Cell_type_class", "Cell_type")

ChIP_ATLAS_exps <- ChIP_ATLAS_exps[!(ChIP_ATLAS_exps$Antigen %in% c("Epitope tags","Biotin","GFP")),]
pie(table(ChIP_ATLAS_exps$Antigen_class), main = "ChIP-ATLAS mm9\nEmbryonic Stem Cell\ndataset")
pie(sort(table(ChIP_ATLAS_exps$Antigen), decreasing = T), main = "ChIP-ATLAS mm9\nEmbryonic Stem Cell\ndataset")
#View(ChIP_ATLAS_exps[ChIP_ATLAS_exps$Antigen %in% "Biotin",])


write.table(x = ChIP_ATLAS_exps$ID, file = "~/Documents/Shayan/BioInf/lncRNA/ChipATLAS/mm9_mESC_exp_names.txt", 
            row.names = F, col.names = F, quote = F, append = F)


sum(aa_my_Exp_id %in% aa_atlas_exp$V1)
sum((aa_my_Exp_id %in% aa_atlas_exp$V1) & !(aa_my_Exp_id %in% ChIP_ATLAS_exps$V1))
aaw <- which((aa_my_Exp_id %in% aa_atlas_exp$V1) & !(aa_my_Exp_id %in% ChIP_ATLAS_exps$V1))
View(aa_atlas_exp[aa_atlas_exp$V1 %in% aa_my_Exp_id[aaw],])


# write a grep command to filter the full file for these names: "~/Documents/Shayan/BioInf/lncRNA/ChipATLAS/mm9_mESC_exp_names.txt"
# awk 'NR==FNR{a[$2];next} ($1 in a)' file1.txt file2.txt
# awk 'NR==FNR{A[$1];next}$4 in A'  mm9_mESC_exp_names.txt allPeaks_light.mm9.05.bed > allPeaks_light_mESC.mm9.05.bed
cat(c("#!/bin/bash\n"), file = "~/Documents/Shayan/BioInf/lncRNA/ChipATLAS/filter_atlas.sh", append = F)
aaf1 <- "~/Documents/Shayan/BioInf/lncRNA/ChipATLAS/allPeaks_light.mm9.05.bed"
aaf2 <- "~/Documents/Shayan/BioInf/lncRNA/ChipATLAS/mm9_mESC_exp_names.txt"
aaout <- "~/Documents/Shayan/BioInf/lncRNA/ChipATLAS/allPeaks_light_mesc.mm9.05.bed"
cat(c("awk 'NR==FNR{A[$1];next}$4 in A'",aaf2, aaf1, ">",aaout,"\n"), 
    file = "~/Documents/Shayan/BioInf/lncRNA/ChipATLAS/filter_atlas.sh", sep = " ", append = T)

# write grep commands to make antigen specific bed files
aa_ant_agg <- aggregate(ChIP_ATLAS_exps[c("ID")],
                        by = ChIP_ATLAS_exps[c("Antigen")],
                        FUN = c)
aa_id_name <- paste0("/shared-mounts/sinhas/tabebor2/mESC_chip/Antigen_ids/",aa_ant_agg$Antigen, "__ID.txt")
for(i in 1:nrow(aa_ant_agg)){
  write.table(x = aa_ant_agg$ID[[i]], file = aa_id_name[i], 
              row.names = F, col.names = F, quote = F, append = F)
}
aaln <- unlist(lapply(aa_ant_agg$ID, length))
aa_bed_name <- paste0("/shared-mounts/sinhas/tabebor2/mESC_chip/Antigen_specific_beds/",aa_ant_agg$Antigen, "__",aaln,".bed")
cat(c("#!/bin/bash\n"), 
    file = "~/Documents/Shayan/BioInf/lncRNA/ChipATLAS/separate_by_ag.sh",
    append = F)
for(i in 1:nrow(aa_ant_agg)){
  aaf1 <- "/shared-mounts/sinhas/tabebor2/mESC_chip/allPeaks_light_mESC.mm9.05.bed"
  cat(c("awk 'NR==FNR{A[$1];next}$4 in A'",aa_id_name[i], aaf1, ">",aa_bed_name[i],"\n"), 
      file = "~/Documents/Shayan/BioInf/lncRNA/ChipATLAS/separate_by_ag.sh", 
      sep = " ", append = T)
}


aa_mm9_catlas_files <- list.files("~/Documents/Shayan/BioInf/lncRNA/ChipATLAS/Antigen_specific_beds.nosync", pattern = ".bed", full.names = T)
aa_mm9_catlas_name2 <- list.files("~/Documents/Shayan/BioInf/lncRNA/ChipATLAS/Antigen_specific_beds.nosync", pattern = ".bed", full.names = F)
aanameee1 <- unlist(lapply(strsplit(aa_mm9_catlas_name2, "__"), "[[", 1))
aanameee2 <- unlist(lapply(strsplit(aa_mm9_catlas_name2, "__"), "[[", 2))
aanameee2 <- unlist(lapply(strsplit(aanameee2, "\\."), "[[", 1))
aanameee2 <- as.numeric(aanameee2)

# read bed files
aa_catlas <- list()
for(i in 1:length(aa_mm9_catlas_files)){
  print(i)
  print(aanameee1[i])
  aa_catlas[[i]] <- read.table(aa_mm9_catlas_files[i], stringsAsFactors = F)
}
mESC_CHIPATLAS_GR_list <- lapply(aa_catlas,
                           makeGRangesFromDataFrame,
                           keep.extra.columns = T,
                           seqnames.field = "V1",
                           start.field = "V2",
                           end.field = "V3")
names(mESC_CHIPATLAS_GR_list) <- aa_mm9_catlas_name2
mESC_ChIP_ATLAS_tiled_list <- list()
mESC_ChIP_ATLAS_tiled_list_GR <- list()

for(i in 1:length(mESC_CHIPATLAS_GR_list)){
  print(i)
  print(aanameee1[i])
  
  aaov <- findOverlaps(MM9_1kb_tiled_GR, mESC_CHIPATLAS_GR_list[[i]])
  print("overlap done")
  aadf <- data.frame(tile_index = aaov@from,
                     chip_binding_index = aaov@to, 
                     chip_score = mESC_CHIPATLAS_GR_list[[i]]$V5[aaov@to],
                     exp_name = mESC_CHIPATLAS_GR_list[[i]]$V4[aaov@to])
  aadf$exp_name <- as.numeric(aadf$exp_name )
  aadf_agg1 <- aggregate(aadf[c("chip_score")],
                        by = aadf[c("tile_index")],
                        FUN = max)
  print("first agg done")
  aadf_agg_int <- aggregate(aadf[c("exp_name")],
                        by = aadf[c("tile_index")],
                        FUN = c)
  print("second agg done")
  aadf_agg_2 <- lapply(aadf_agg_int$exp_name, unique)
  print("uniq done")
  aadf_agg_2 <-  unlist(lapply(aadf_agg_2, length))
  print("length done")
  #aa_nu_int <- unlist(lapply(aadf_agg$interaction_index, length))
  mESC_ChIP_ATLAS_tiled_list[[i]] <- cbind(MM9_1kb_tiled[aadf_agg1$tile_index,], aadf_agg1$chip_score, aadf_agg_2, aadf_agg_2/aanameee2[i])
  colnames(mESC_ChIP_ATLAS_tiled_list[[i]])[c(5,6,7)]<- c("score_max", "nu_uniq_exp", "percent_uniq_exp")
  mESC_ChIP_ATLAS_tiled_list_GR[[i]] <- makeGRangesFromDataFrame(mESC_ChIP_ATLAS_tiled_list[[i]], keep.extra.columns = T)
  
}
names(mESC_ChIP_ATLAS_tiled_list) <- names(mESC_CHIPATLAS_GR_list)
names(mESC_ChIP_ATLAS_tiled_list_GR) <- names(mESC_CHIPATLAS_GR_list)

ChIPATLAS_features <- matrix(nrow = length(aaMM9_1kb_tiled_GR_filtered),
                       ncol = length(mESC_ChIP_ATLAS_tiled_list))

rownames(ChIPATLAS_features) <- MM9_1kb_tiled_GR_filtered$tile
colnames(ChIPATLAS_features) <- names(mESC_ChIP_ATLAS_tiled_list)
for(i in 1:length(mESC_ChIP_ATLAS_tiled_list)){
  print(i)
  aatmp <- mESC_ChIP_ATLAS_tiled_list[[i]]
  aatmp <- aatmp[aatmp$tile %in% rownames(ChIPATLAS_features),]
  aatmp <- aatmp[(aatmp$nu_uniq_exp >= 3 | aatmp$percent_uniq_exp >= 0.5),]
  ChIPATLAS_features[match(aatmp$tile, rownames(ChIPATLAS_features)),i] <- aatmp$score_max
}
save(list= c("ChIPATLAS_features"), file = "~/Documents/Shayan/BioInf/lncRNA/ChIPATLAS_features.RData")

mESC_CHIPATLAS_lncRNA_features <- matrix(nrow = length(lncRNA_chosen_gt1k_uniqTiles_GR), 
                                         ncol = length(mESC_ChIP_ATLAS_tiled_list))
rownames(mESC_CHIPATLAS_lncRNA_features) <- lncRNA_chosen_gt1k_uniqTiles_GR$gene_name
colnames(mESC_CHIPATLAS_lncRNA_features) <- names(mESC_ChIP_ATLAS_tiled_list)
length(mESC_CHIPATLAS_GR_list)
for(i in 1:length(mESC_CHIPATLAS_GR_list)){
  print(i)
  aaov <- findOverlaps(lncRNA_chosen_gt1k_uniqTiles_GR,
                       mESC_CHIPATLAS_GR_list[[i]])
  if(length(aaov) > 0){
    print("overlap done")
    aadf <- data.frame(lnc_index = aaov@from,
                       chip_binding_index = aaov@to, 
                       chip_score = mESC_CHIPATLAS_GR_list[[i]]$V5[aaov@to],
                       exp_name = mESC_CHIPATLAS_GR_list[[i]]$V4[aaov@to])
    aadf$exp_name <- as.numeric(aadf$exp_name)
    aadf_agg1 <- aggregate(aadf[c("chip_score")],
                           by = aadf[c("lnc_index")],
                           FUN = max, simplify = F)
    print("first agg done")
    aadf_agg_int <- aggregate(aadf[c("exp_name")],
                              by = aadf[c("lnc_index")],
                              FUN = c, simplify = F)
    print("second agg done")
    aadf_agg_2 <- lapply(aadf_agg_int$exp_name, unique)
    print("uniq done")
    aadf_agg_2 <-  unlist(lapply(aadf_agg_2, length))
    print("length done")
    aa_my_ind <- match(aadf_agg1$lnc_index, 
                       c(1:length(lncRNA_chosen_gt1k_uniqTiles_GR)))
    aadf_agg_3 <- aadf_agg_2/aanameee2[i]
    if(all(aadf_agg_int$lnc_index == aadf_agg1$lnc_index)){
      aadf_agg1$chip_score[(aadf_agg_2 < 3) & (aadf_agg_3 < 0.5)] <- NA
      mESC_CHIPATLAS_lncRNA_features[aa_my_ind,i] <- as.numeric(aadf_agg1$chip_score)
    }else{
      print("stop stop")
    }
  }
}
aatst <- mESC_CHIPATLAS_lncRNA_features
aatst[is.na(aatst)] <- 0
aatst <- aatst + 1

library(gplots)
z <- zClust(x=log10(aatst), scale="none", method="average")

require(RColorBrewer)
aacols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(7)

heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=z$Colv)

heatmap.2(log10(aatst), Rowv = T, Colv = T, trace = "none")
####################################################################################################
####################################################################################################
# getting the kmers
start(MM9_1kb_tiled_GR) <- start(MM9_1kb_tiled_GR) + 1
end(MM9_1kb_tiled_GR) <- end(MM9_1kb_tiled_GR) + 1

aa_tile_Seq <- getSeq(x = BSgenome.Mmusculus.UCSC.mm9, names = MM9_1kb_tiled_GR[3001])
aa_tile_Seq_kmer <- get.kmers(as.character(aa_tile_Seq), .k = 5)
aa_tile_Seq_kmer

which(MM9_1kb_tiled_chrname == "chr2")[1:5]
aa_tile_Seq <- getSeq(x = BSgenome.Mmusculus.UCSC.mm9, names = MM9_1kb_tiled_GR[200196])
aa_tile_Seq_kmer <- get.kmers(as.character(aa_tile_Seq), .k = 5)
aa_tile_Seq_kmer

which(MM9_1kb_tiled_chrname == "chrX")[1:5]
aa_tile_Seq <- getSeq(x = BSgenome.Mmusculus.UCSC.mm9, names = MM9_1kb_tiled_GR[2472354 + 2999])
aa_tile_Seq_kmer <- get.kmers(as.character(aa_tile_Seq), .k = 5)
aa_tile_Seq_kmer

options(scipen=999)
write.table(lncRNA_chosen_gt1k_uniqTiles_GR, 
            file="~/Documents/Shayan/BioInf/lncRNA/lncRNA_28_mm9.bed",
            quote=F, sep="\t", row.names=F, col.names=F)
options(scipen=0)
aa_lncRNA_Seq <- getSeq(x = BSgenome.Mmusculus.UCSC.mm9, names = lncRNA_chosen_gt1k_uniqTiles_GR)
names(aa_lncRNA_Seq) <- lncRNA_chosen_gt1k_uniqTiles_GR$gene_name
lncRNA_sequence_28 <- aa_lncRNA_Seq

par(mar = c(8,4,4,4))
barplot(width(lncRNA_sequence_28), las = 2, names.arg = names(lncRNA_sequence_28), main = "lncRNA width")

ShortRead::writeFasta(object = aa_lncRNA_Seq, 
                      file = "~/Documents/Shayan/BioInf/lncRNA/lncRNA_28_mm9_seq.fasta")


# Tile_kmer_features <- matrix(nrow = length(MM9_1kb_tiled_GR_filtered), ncol = (1024))
# rownames(Tile_kmer_features) <- MM9_1kb_tiled_GR_filtered$tile
# colnames(Tile_kmer_features) 

save(list = c("MM9_1kb_tiled_GR_filtered"), file = "~/Documents/Shayan/BioInf/lncRNA/MM9_1kb_tiled_GR_filtered.RData")
# aa_myseq <- getSeq(x = BSgenome.Mmusculus.UCSC.mm9, names = MM9_1kb_tiled_GR_filtered )
# names(aa_myseq) <- MM9_1kb_tiled_GR_filtered$tile
# aa_kmer_list <- lapply(as.character(aa_myseq), get.kmers, .k=5)
 all_5mers <- sort(generate.kmers(.k = 5))
map_to_all5mer <- function(x, my_5mers = all_5mers){
  aa <- numeric(length = length(my_5mers))
  aa[match(x[[1]], my_5mers)] <- x[[2]]
  return(aa)
}
# aa_kmer_list2 <- lapply(aa_kmer_list, map_to_all5mer, all_5mers)
# aa_kmer_all <- do.call(rbind, aa_kmer_list2)
# Tile_kmer_features <- aa_kmer_all
# rownames(Tile_kmer_features) <- MM9_1kb_tiled_GR_filtered$tile
# colnames(Tile_kmer_features) <- all_5mers
# save(list = c("Tile_kmer_features"), 
#      file = "~/Documents/Shayan/BioInf/lncRNA/Tile_kmer_features.RData")

save(list = c("Tile_kmer_features_partition2"),
     file = "~/Documents/Shayan/BioInf/lncRNA/Tile_kmer_features_partition2.RData")

# using Jellyfish
# writing each tile in a different fasta file --> breaking previously written files
aa_file_list <- list.files("/Users/Shayan/Documents/Shayan/BioInf/lncRNA/MM9_tiles.nosync/")

#cat("#!/bin/bash\n", file = "break_to_tiles_mm9.sh", append = F)
for(i in 1:length(aa_file_list)){
  cat(c("/shared-mounts/sinhas/tabebor2/lncRNA/break_fasta.sh <", paste0("../Seq/", aa_file_list[i], "\n")),
      file = "break_to_tiles_mm9.job", append = (i != 1), sep = " " )
}

# write jellyfish jobs

for(i in 1:length(MM9_1kb_tiled_GR_filtered)){
  if((i %% 1000) == 0){
    print(i)
    }
  cat(c("jellyfish count -m 5 -t 1 --text  -s 1024 -o", 
        paste0("kmer_by_Tile/", MM9_1kb_tiled_GR_filtered$tile[i], ".5mer"), 
        paste0("Seq_by_Tile/", MM9_1kb_tiled_GR_filtered$tile[i], ".fa\n")),
      file = "Jellyfish_mm9.job", append = !(i == 1), sep = " " )
}
# read the lncRNA results for all tiles of partition1: Partition_1_5mer_mat.RData
load("~/Documents/Shayan/BioInf/lncRNA/Partition_1_5mer_mat.RData")
colnames(Partition_1_5mer_mat) <- all_5mers



#### Get kmers for the lncRNAs


aa_myseq_lnc <- getSeq(x = BSgenome.Mmusculus.UCSC.mm9, names = lncRNA_chosen_gt1k_uniqTiles_GR )
names(aa_myseq_lnc) <- lncRNA_chosen_gt1k_uniqTiles_GR$gene_name
aa_kmer_list_lnc <- lapply(as.character(aa_myseq_lnc), get.kmers, .k=5)
aa_kmer_list2_lnc <- lapply(aa_kmer_list_lnc, map_to_all5mer, all_5mers)
aa_kmer_all_lnc <- do.call(rbind, aa_kmer_list2_lnc)
aa_kmer_all_lnc_norm <- aa_kmer_all_lnc/width(lncRNA_chosen_gt1k_uniqTiles_GR)
boxplot.matrix(t(aa_kmer_all_lnc/width(lncRNA_chosen_gt1k_uniqTiles_GR)), las = 2)
kmer_features_lncRNA28_normalized <- aa_kmer_all_lnc_norm
colnames(kmer_features_lncRNA28_normalized) <- all_5mers

aawlkmer <- sort(colSums(kmer_features_lncRNA28_normalized), decreasing = T, index.return=T)
aamalats <- sort(kmer_features_lncRNA28_normalized[1,], decreasing = T, index.return=T)
which(colnames(kmer_features_lncRNA28_normalized)[aamalats$ix] == "GCTGG")
which(colnames(kmer_features_lncRNA28_normalized)[aamalats$ix] == "CCAGC")
which(colnames(kmer_features_lncRNA28_normalized)[aamalats$ix] == "CCAGC")




####################################################################################################
####################################################################################################
# TRIPLEX scores

# rgt-TDF get_TTS -i ../../../Data_sets/Cis_RNA_interaction/Cis_RNA28_interacting_positive.bed  -tts MEG3_TTSs.bed -r FENDRR.fasta  -organism mm9

# write FASTA files for each lncRNA separately
aa_fasta_name <- paste0("~/Documents/Shayan/BioInf/lncRNA/TDF/", names(lncRNA_sequence_28), ".fasta")
aa_out_name <- paste0("~/Documents/Shayan/BioInf/lncRNA/TDF/", names(lncRNA_sequence_28), "__triplex.bed")
cat(c("#!/bin/bash\n"), 
    file = "~/Documents/Shayan/BioInf/lncRNA/TDF/TDF_per_lncRNA_run.sh",
    append = F)
for(i in 1:length(lncRNA_sequence_28)){
  ShortRead::writeFasta(object = lncRNA_sequence_28[i], 
                        file = aa_fasta_name[i])
  cat(c("rgt-TDF get_TTS -i ~/Documents/Shayan/BioInf/lncRNA/mm9_1kb_tiles.bed -tts",
        aa_out_name[i], " -r ", aa_fasta_name[i], " -organism mm9\n"),
      file = "~/Documents/Shayan/BioInf/lncRNA/TDF/TDF_per_lncRNA_run.sh",
      append = T, sep = " ")
}
####################################################################################################
# TDF doesn't work well. installed triplexator on hal and will use that
# write fasta file of mm9 genome, named by tiles. remove the first 3000 tiles for each chromosome.

start(MM9_1kb_tiled_GR) <- start(MM9_1kb_tiled_GR) + 1
end(MM9_1kb_tiled_GR) <- end(MM9_1kb_tiled_GR) + 1

aa_tile_Seq <- getSeq(x = BSgenome.Mmusculus.UCSC.mm9, names = MM9_1kb_tiled_GR[3000])
aa_tile_Seq_kmer <- get.kmers(as.character(aa_tile_Seq), .k = 5)
aa_tile_Seq_kmer

aaMM9_1kb_tiled_GR_filtered <- MM9_1kb_tiled_GR
aarem <- numeric(0)
aa_unq_chr <- unique(MM9_1kb_tiled_chrname)
for(i in 1:length(aa_unq_chr)){
  aarem <-  c(aarem, which(MM9_1kb_tiled_chrname == aa_unq_chr[i])[1:min(3000, sum(MM9_1kb_tiled_chrname == aa_unq_chr[i]))])
}
aaMM9_1kb_tiled_GR_filtered <- aaMM9_1kb_tiled_GR_filtered[-c(aarem)]
index_first_3000_perchr <- aarem
MM9_1kb_tiled_GR_filtered <- aaMM9_1kb_tiled_GR_filtered

aa_tile_Seq <- getSeq(x = BSgenome.Mmusculus.UCSC.mm9, names = aaMM9_1kb_tiled_GR_filtered[1:1000])
names(aa_tile_Seq) <- aaMM9_1kb_tiled_GR_filtered$tile[1:1000]

aal <- ceiling(length(aaMM9_1kb_tiled_GR_filtered) / 1000)
aacnt <- 1
for(i in 1:aal){
  aast <-((i - 1)* 1000) + 1
  aaed <- min((i * 1000), length(aaMM9_1kb_tiled_GR_filtered))
  aa_tile_Seq <- getSeq(x = BSgenome.Mmusculus.UCSC.mm9, names = aaMM9_1kb_tiled_GR_filtered[aast:aaed ])
  names(aa_tile_Seq) <- aaMM9_1kb_tiled_GR_filtered$tile[aast:aaed]
  ShortRead::writeFasta(object = aa_tile_Seq, 
                        file = paste0("~/Documents/Shayan/BioInf/lncRNA/MM9_tiles.nosync/mm9_seg_", i, ".fasta"))
  for(j in 1:length(lncRNA_sequence_28))
  cat(c("/shared-mounts/sinhas/tabebor2/lncRNA/Triplexator/triplexator/triplexator/bin/triplexator ",
        "-ss",paste0("/shared-mounts/sinhas/tabebor2/lncRNA/lncRNA_28/", names(lncRNA_sequence_28)[j], ".fasta"),
        "-ds", paste0("/shared-mounts/sinhas/tabebor2/lncRNA/Seq/mm9_seg_", i, ".fasta"),
        "-l 10 -e 20 -g 40 -of 2 rm 0", 
        "-o", paste0("mm9_seg_", i, "__", names(lncRNA_sequence_28)[j]),
        "-od /shared-mounts/sinhas/tabebor2/lncRNA/Triple_output\n"), 
      file = "~/Documents/Shayan/BioInf/lncRNA/triplexator.job",
      sep = " ",
      append = !((i == 1) & (j == 1)))
}

aa_job <- read.table("~/Documents/Shayan/BioInf/lncRNA/job_names.txt", stringsAsFactors = F)
aa_job <- aa_job$V1

cat(c("#!/bin/bash\n"), 
    file = "~/Documents/Shayan/BioInf/lncRNA/trip_job_TO_submit.sh",
    append = F)
for(i in 1:length(aa_job)){
  cat(c("/shared-mounts/sinhas/tabebor2/lncRNA/create_jobs_submit.pl ", 
        aa_job[i], 1000,
        paste0("tmp_",aa_job[i]), " > ",
        paste0(aa_job[i], ".submit\n")), 
      file = "~/Documents/Shayan/BioInf/lncRNA/trip_job_TO_submit.sh",
      sep = " ",
      append = T)
}

# make them executable
cat(c("#!/bin/bash\n"), 
    file = "~/Documents/Shayan/BioInf/lncRNA/submit_exec.sh",
    append = F)
for(i in 1:length(aa_job)){
  cat(c("chmod +x", 
        paste0(aa_job[i], ".submit\n")), 
      file = "~/Documents/Shayan/BioInf/lncRNA/submit_exec.sh",
      sep = " ",
      append = T)
}
#execute them
cat(c("#!/bin/bash\n"), 
    file = "~/Documents/Shayan/BioInf/lncRNA/submit_the_jobs.sh",
    append = F)
for(i in 1:length(aa_job)){
  cat(c(paste0("./",aa_job[i], ".submit\n")), 
      file = "~/Documents/Shayan/BioInf/lncRNA/submit_the_jobs.sh",
      sep = " ",
      append = T)
}


# read results from hal

# gather summaries per lncRNA
# cat Triple_output/*__2410003L11Rik.summary > all__2410003L11Rik.summary
cat("#!/bin/bash\n", file = "~/Documents/Shayan/BioInf/lncRNA/Gather_Triplex.sh", append = F)
for(i in 1:ncol(territory_assigned_interaction_3Mfilered)){
  cat(c(paste0("cat Triple_output/*__", colnames(territory_assigned_interaction_3Mfilered)[i],
               ".summary"), " > "  , paste0("all__", colnames(territory_assigned_interaction_3Mfilered)[i], ".summary\n")), 
      file = "~/Documents/Shayan/BioInf/lncRNA/Gather_Triplex.sh", sep= " " ,append = T)
  
}
# filter only the chosen tiles for each lncRNA
for(i in 1:length(Partition_1)){
  write.table(Partition_1_dfs$tile_name[Partition_1_dfs$owner %in% names(Partition_1)[i]],
              row.names = F,
              col.names = F, 
              quote = F,
              file = paste0("~/Documents/Shayan/BioInf/lncRNA/partition_1_tile_", names(Partition_1)[i], ".txt"))
}

cat(c("#!/bin/bash\n"), file = "~/Documents/Shayan/BioInf/lncRNA/filter_triplex.sh", append = F)
for(i in 1:length(Partition_1)){
  aaf1 <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/Triple_summary_all/all__", names(Partition_1)[i], ".summary")
  aaf2 <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/Partition_1/partition_1_tile_", names(Partition_1)[i], ".txt")
  aaout<- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/Partition_1/triplex_part1_",names(Partition_1)[i], ".txt")
  cat(c("awk 'NR==FNR{A[$1];next}$1 in A'",aaf2, aaf1, "|  cut -f 1,2,3 ", ">",aaout,"\n"), 
      file = "~/Documents/Shayan/BioInf/lncRNA/filter_triplex.sh", sep = " ", append = T)
  
}

aafnames_full <- list.files(path = "~/Documents/Shayan/BioInf/lncRNA/Partition_1", pattern = "triplex_*", full.names = T)
aafnames <- list.files(path = "~/Documents/Shayan/BioInf/lncRNA/Partition_1", pattern = "triplex_*", full.names = F)
aafnames2 <- unlist(lapply(strsplit(unlist(lapply(strsplit(aafnames, "\\."), "[[", 1)), "_"), "[[", 3))
aa_trip_df <- list()
for(i in 1:length(aafnames_full)){
  aa_trip_df[[i]] <- read.table(aafnames_full[i], stringsAsFactors = F)
}
names(aa_trip_df) <- aafnames2
Triplex_feaures_partition1 <- do.call(rbind, aa_trip_df)
colnames(Triplex_feaures_partition1) <- c("tile_name", "owner", "triplexScore")

aaez <- setdiff(Partition_1_dfs$tile_name, Triplex_feaures_partition1$tile_name)
aaez2 <- Partition_1_dfs[Partition_1_dfs$tile_name %in% aaez, 1:2]
aaez2 <-  cbind(aaez2, rep(0, nrow(aaez2)))
colnames(aaez2)[3] <- "triplexScore"
Triplex_feaures_partition1 <- rbind(Triplex_feaures_partition1, aaez2)
### reading triplexes for partition 2
aafnames_full <- list.files(path = "~/Documents/Shayan/BioInf/lncRNA/Partition_2", pattern = "triplex_*", full.names = T)
aafnames <- list.files(path = "~/Documents/Shayan/BioInf/lncRNA/Partition_2", pattern = "triplex_*", full.names = F)
aafnames2 <- unlist(lapply(strsplit(unlist(lapply(strsplit(aafnames, "\\."), "[[", 1)), "_"), "[[", 3))
aa_trip_df <- list()
for(i in 1:length(aafnames_full)){
  aa_trip_df[[i]] <- read.table(aafnames_full[i], stringsAsFactors = F)
}
names(aa_trip_df) <- aafnames2
Triplex_feaures_partition2 <- do.call(rbind, aa_trip_df)
colnames(Triplex_feaures_partition2) <- c("tile_name", "owner", "triplexScore")
aaez <- setdiff(Partition_2_dfs$tile_name, Triplex_feaures_partition2$tile_name)
aaez2 <- Partition_2_dfs[Partition_2_dfs$tile_name %in% aaez, 1:2]
aaez2 <-  cbind(aaez2, rep(0, nrow(aaez2)))
colnames(aaez2)[3] <- "triplexScore"
Triplex_feaures_partition2 <- rbind(Triplex_feaures_partition2, aaez2)

####################################################################################################
# read triplex features for partition2

# filter only the chosen tiles for each lncRNA
for(i in 1:length(Partition_2)){
  write.table(Partition_2_dfs$tile_name[Partition_2_dfs$owner %in% names(Partition_2)[i]],
              row.names = F,
              col.names = F, 
              quote = F,
              file = paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_2/partition_2_tile_", names(Partition_1)[i], ".txt"))
}

cat(c("#!/bin/bash\n"), file = "~/Documents/Shayan/BioInf/lncRNA/filter_triplex_p2.sh", append = F)
for(i in 1:length(Partition_2)){
  aaf1 <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/Triple_summary_all/all__", names(Partition_2)[i], ".summary")
  aaf2 <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/Partition_2/partition_2_tile_", names(Partition_2)[i], ".txt")
  aaout<- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/Partition_2/triplex_part2_",names(Partition_2)[i], ".txt")
  cat(c("awk 'NR==FNR{A[$1];next}$1 in A'",aaf2, aaf1, "|  cut -f 1,2,3 ", ">",aaout,"\n"), 
      file = "~/Documents/Shayan/BioInf/lncRNA/filter_triplex_p2.sh", sep = " ", append = T)
  
}
# read results
aafnames_full <- list.files(path = "~/Documents/Shayan/BioInf/lncRNA/Partition_2", pattern = "triplex_*", full.names = T)
aafnames <- list.files(path = "~/Documents/Shayan/BioInf/lncRNA/Partition_2", pattern = "triplex_*", full.names = F)
aafnames2 <- unlist(lapply(strsplit(unlist(lapply(strsplit(aafnames, "\\."), "[[", 1)), "_"), "[[", 3))
aa_trip_df <- list()
for(i in 1:length(aafnames_full)){
  aa_trip_df[[i]] <- read.table(aafnames_full[i], stringsAsFactors = F)
}
names(aa_trip_df) <- aafnames2
Triplex_feaures_partition2 <- do.call(rbind, aa_trip_df)
colnames(Triplex_feaures_partition2) <- c("tile_name", "owner", "triplexScore")

aaez <- setdiff(Partition_2_dfs$tile_name, Triplex_feaures_partition2$tile_name)
aaez2 <- Partition_2_dfs[Partition_2_dfs$tile_name %in% aaez, 1:2]
aaez2 <-  cbind(aaez2, rep(0, nrow(aaez2)))
colnames(aaez2)[3] <- "triplexScore"
Triplex_feaures_partition2 <- rbind(Triplex_feaures_partition2, aaez2)
save(list = c("Triplex_feaures_partition2"), file = "~/Documents/Shayan/BioInf/lncRNA/Triplex_feaures_partition2.RData")

####################################################################################################
# read triplex features for partition4

# filter only the chosen tiles for each lncRNA
for(i in 1:length(Partition_4)){
  write.table(Partition_4_dfs$tile_name[Partition_4_dfs$owner %in% names(Partition_4)[i]],
              row.names = F,
              col.names = F, 
              quote = F,
              file = paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_4/partition_4_tile_", names(Partition_4)[i], ".txt"))
}

cat(c("#!/bin/bash\n"), file = "~/Documents/Shayan/BioInf/lncRNA/filter_triplex_p4.sh", append = F)
for(i in 1:length(Partition_4)){
  aaf1 <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/Triple_summary_all/all__", names(Partition_4)[i], ".summary")
  aaf2 <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/Partition_4/partition_4_tile_", names(Partition_4)[i], ".txt")
  aaout<- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/Partition_4/triplex_part4_",names(Partition_4)[i], ".txt")
  cat(c("awk 'NR==FNR{A[$1];next}$1 in A'",aaf2, aaf1, "|  cut -f 1,2,3 ", ">",aaout,"\n"), 
      file = "~/Documents/Shayan/BioInf/lncRNA/filter_triplex_p4.sh", sep = " ", append = T)
  
}
# read results
aafnames_full <- list.files(path = "~/Documents/Shayan/BioInf/lncRNA/Partition_4", pattern = "triplex_*", full.names = T)
aafnames <- list.files(path = "~/Documents/Shayan/BioInf/lncRNA/Partition_4", pattern = "triplex_*", full.names = F)
aafnames2 <- unlist(lapply(strsplit(unlist(lapply(strsplit(aafnames, "\\."), "[[", 1)), "_"), "[[", 3))
aa_trip_df <- list()
for(i in 1:length(aafnames_full)){
  aa_trip_df[[i]] <- read.table(aafnames_full[i], stringsAsFactors = F)
}
names(aa_trip_df) <- aafnames2
Triplex_feaures_partition4 <- do.call(rbind, aa_trip_df)
colnames(Triplex_feaures_partition4) <- c("tile_name", "owner", "triplexScore")

aaez <- setdiff(Partition_4_dfs$tile_name, Triplex_feaures_partition4$tile_name)
aaez2 <- Partition_4_dfs[Partition_4_dfs$tile_name %in% aaez, 1:2]
aaez2 <-  cbind(aaez2, rep(0, nrow(aaez2)))
colnames(aaez2)[3] <- "triplexScore"
Triplex_feaures_partition4 <- rbind(Triplex_feaures_partition4, aaez2)
Triplex_feaures_partition4 <- Triplex_feaures_partition4[match(Partition_4_dfs$tile_name, Triplex_feaures_partition4$tile_name ),]
save(list = c("Triplex_feaures_partition4"), file = "~/Documents/Shayan/BioInf/lncRNA/Triplex_feaures_partition4.RData")
####################################################################################################
# read triplex features for partition5

# filter only the chosen tiles for each lncRNA
for(i in 1:length(Partition_5)){
  write.table(Partition_5_dfs$tile_name[Partition_5_dfs$owner %in% names(Partition_5)[i]],
              row.names = F,
              col.names = F, 
              quote = F,
              file = paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_5/partition_5_tile_", names(Partition_5)[i], ".txt"))
}

cat(c("#!/bin/bash\n"), file = "~/Documents/Shayan/BioInf/lncRNA/filter_triplex_p5.sh", append = F)
for(i in 1:length(Partition_5)){
  aaf1 <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/Triple_summary_all/all__", names(Partition_5)[i], ".summary")
  aaf2 <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/Partition_5/partition_5_tile_", names(Partition_5)[i], ".txt")
  aaout<- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/Partition_5/triplex_part5_",names(Partition_5)[i], ".txt")
  cat(c("awk 'NR==FNR{A[$1];next}$1 in A'",aaf2, aaf1, "|  cut -f 1,2,3 ", ">",aaout,"\n"), 
      file = "~/Documents/Shayan/BioInf/lncRNA/filter_triplex_p5.sh", sep = " ", append = T)
  
}
# read results
aafnames_full <- list.files(path = "~/Documents/Shayan/BioInf/lncRNA/Partition_5", pattern = "triplex_*", full.names = T)
aafnames <- list.files(path = "~/Documents/Shayan/BioInf/lncRNA/Partition_5", pattern = "triplex_*", full.names = F)
aafnames2 <- unlist(lapply(strsplit(unlist(lapply(strsplit(aafnames, "\\."), "[[", 1)), "_"), "[[", 3))
aa_trip_df <- list()
for(i in 1:length(aafnames_full)){
  aa_trip_df[[i]] <- read.table(aafnames_full[i], stringsAsFactors = F)
}
names(aa_trip_df) <- aafnames2
Triplex_feaures_partition5 <- do.call(rbind, aa_trip_df)
colnames(Triplex_feaures_partition5) <- c("tile_name", "owner", "triplexScore")

aaez <- setdiff(Partition_5_dfs$tile_name, Triplex_feaures_partition5$tile_name)
aaez2 <- Partition_5_dfs[Partition_5_dfs$tile_name %in% aaez, 1:2]
aaez2 <-  cbind(aaez2, rep(0, nrow(aaez2)))
colnames(aaez2)[3] <- "triplexScore"
Triplex_feaures_partition5 <- rbind(Triplex_feaures_partition5, aaez2)
Triplex_feaures_partition5 <- Triplex_feaures_partition5[match(Partition_5_dfs$tile_name, Triplex_feaures_partition5$tile_name ),]
save(list = c("Triplex_feaures_partition5"), file = "~/Documents/Shayan/BioInf/lncRNA/Triplex_feaures_partition5.RData")

####################################################################################################
# read triplex features for partition6

# filter only the chosen tiles for each lncRNA
#Partition_6_dfs$owner <- levels(Partition_6_dfs$owner)[as.numeric(Partition_6_dfs$owner)]
aaunq <- unique(Partition_6_dfs$owner)
for(i in 1:length(aaunq)){
  write.table(Partition_6_dfs$tile_name[Partition_6_dfs$owner %in% aaunq[i]],
              row.names = F,
              col.names = F, 
              quote = F,
              file = paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6/partition_6_tile_", aaunq[i], ".txt"))
}

cat(c("#!/bin/bash\n"), file = "~/Documents/Shayan/BioInf/lncRNA/filter_triplex_p6.sh", append = F)
for(i in 1:length(aaunq)){
  aaf1 <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/Triple_summary_all/all__", aaunq[i], ".summary")
  aaf2 <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/Partition_6/partition_6_tile_", aaunq[i], ".txt")
  aaout<- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/Partition_6/triplex_part6_",aaunq[i], ".txt")
  cat(c("awk 'NR==FNR{A[$1];next}$1 in A'",aaf2, aaf1, "|  cut -f 1,2,3 ", ">",aaout,"\n"), 
      file = "~/Documents/Shayan/BioInf/lncRNA/filter_triplex_p6.sh", sep = " ", append = T)
  
}
# read results
aafnames_full <- list.files(path = "~/Documents/Shayan/BioInf/lncRNA/Partition_6", pattern = "triplex_*", full.names = T)
aafnames <- list.files(path = "~/Documents/Shayan/BioInf/lncRNA/Partition_6", pattern = "triplex_*", full.names = F)
aafnames2 <- unlist(lapply(strsplit(unlist(lapply(strsplit(aafnames, "\\."), "[[", 1)), "_"), "[[", 3))
aa_trip_df <- list()
for(i in 1:length(aafnames_full)){
  aa_trip_df[[i]] <- read.table(aafnames_full[i], stringsAsFactors = F)
}
names(aa_trip_df) <- aafnames2
Triplex_feaures_partition6 <- do.call(rbind, aa_trip_df)
colnames(Triplex_feaures_partition6) <- c("tile_name", "owner", "triplexScore")

aaez <- setdiff(Partition_6_dfs$tile_name, Triplex_feaures_partition6$tile_name)
aaez2 <- Partition_6_dfs[Partition_6_dfs$tile_name %in% aaez, 1:2]
aaez2 <-  cbind(aaez2, rep(0, nrow(aaez2)))
colnames(aaez2)[3] <- "triplexScore"
Triplex_feaures_partition6 <- rbind(Triplex_feaures_partition6, aaez2)
Triplex_feaures_partition6 <- Triplex_feaures_partition6[match(Partition_6_dfs$tile_name, Triplex_feaures_partition6$tile_name ),]
save(list = c("Triplex_feaures_partition6"), file = "~/Documents/Shayan/BioInf/lncRNA/Triplex_feaures_partition6.RData")
####### read more features
aanems <- levels(Partition_6_chunk_dfs$owner)
cat(c("#!/bin/bash\n"), file = "~/Documents/Shayan/BioInf/lncRNA/filter_triplex_p6_three.sh", append = F)
for(i in 1:length(aanems)){
  aaf1 <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/Triple_summary_all/all__", aanems[i], ".summary")
  aaf2 <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/Partition_6/partition_6_tile_", aanems[i], ".txt")
  aaout<- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/Partition_6/triplex_part6_",aanems[i], "__Three.txt")
  cat(c("awk 'NR==FNR{A[$1];next}$1 in A'",aaf2, aaf1, "|  cut -f 1,2,3,5,7,9 ", ">",aaout,"\n"), 
      file = "~/Documents/Shayan/BioInf/lncRNA/filter_triplex_p6_three.sh", sep = " ", append = T)
  
}
# read results
aafnames_full <- list.files(path = "~/Documents/Shayan/BioInf/lncRNA/Partition_6", pattern = "*__Three.txt", full.names = T)
aafnames <- list.files(path = "~/Documents/Shayan/BioInf/lncRNA/Partition_6", pattern = "*__Three.txt", full.names = F)
aafnames2 <- unlist(lapply(strsplit(unlist(lapply(strsplit(aafnames, "\\."), "[[", 1)), "_"), "[[", 3))
aa_trip_df <- list()
for(i in 1:length(aafnames_full)){
  aa_trip_df[[i]] <- read.table(aafnames_full[i], stringsAsFactors = F)
}
#names(aa_trip_df) <- aafnames2
#Partition_6_chunk_dfs$owner <- levels(Partition_6_chunk_dfs$owner)[as.numeric(Partition_6_chunk_dfs$owner)]
Triplex_feaures_partition6_extended <- do.call(rbind, aa_trip_df)
colnames(Triplex_feaures_partition6_extended) <- c("tile_name", "owner", "triplex_total", "triplex_GA", "triplex_TC", "triplex_GT")

aaez <- setdiff(Partition_6_dfs$tile_name, Triplex_feaures_partition6_extended$tile_name)
aaez2 <- Partition_6_dfs[Partition_6_dfs$tile_name %in% aaez, 1:2]
aaez2 <-  cbind(aaez2, rep(0, nrow(aaez2)), rep(0, nrow(aaez2)), rep(0, nrow(aaez2)), rep(0, nrow(aaez2)))
colnames(aaez2)[3:6] <- c("triplex_total", "triplex_GA", "triplex_TC", "triplex_GT")
Triplex_feaures_partition6_extended <- rbind(Triplex_feaures_partition6_extended, aaez2)
#Triplex_feaures_partition6_extended <- Triplex_feaures_partition6_extended[match(Partition_6_dfs$tile_name, Triplex_feaures_partition6_extended$tile_name ),]
Triplex_feaures_partition6_extended_merge <- left_join(x = Partition_6_dfs, y = Triplex_feaures_partition6_extended, by = c("tile_name", "owner"))
Triplex_feaures_partition6_extended_merge[is.na(Triplex_feaures_partition6_extended_merge)] <- 0
table(Triplex_feaures_partition6_extended_merge$owner)
table(Partition_6_dfs$owner)
save(list = c("Triplex_feaures_partition6_extended_merge"), 
     file = "~/Documents/Shayan/BioInf/lncRNA/Triplex_feaures_partition6_extended.RData")

####################################################################################################


####################################################################################################
####################################################################################################
# RBP predictions

aa_fname <- list.files("~/Documents/Shayan/BioInf/lncRNA/oRNAment/MM/", pattern = ".bed.gz$")
aa_fname1 <- unlist(lapply(strsplit(aa_fname, "\\."), "[[", 1))
aafiles_full <- paste0("~/Documents/Shayan/BioInf/lncRNA/oRNAment/MM/",aa_fname)
aafiles_full_2 <- paste0("~/Documents/Shayan/BioInf/lncRNA/oRNAment/MM/",aa_fname1, ".bed" )
#system("mkdir ~/Documents/Shayan/BioInf/lncRNA/oRNAment/processed_mm9/")
cat(c("#!/bin/bash\n"), 
    file = "~/Documents/Shayan/BioInf/lncRNA/oRNAment/liftover_to_mm9.sh",
    append = F)
for(i in 1:length(aafiles_full)){
  cat(c("gunzip ", aafiles_full[i], "\n"), 
      file = "~/Documents/Shayan/BioInf/lncRNA/oRNAment/liftover_to_mm9.sh",
      sep = " ",
      append = T)
  cat(c("~/Documents/Shayan/BioInf/liftOver -bedPlus=4", 
        aafiles_full_2[i],
        " ~/Documents/Shayan/BioInf/Liftover_chain_files/mm10ToMm9.over.chain ", 
        paste0("~/Documents/Shayan/BioInf/lncRNA/oRNAment/processed_mm9/",
               aa_fname1[i],
               "_mm9.bed"),
        paste0("~/Documents/Shayan/BioInf/lncRNA/oRNAment/processed_mm9/unconverted_", aa_fname1[i]),
        "\n"),
      file = "~/Documents/Shayan/BioInf/lncRNA/oRNAment/liftover_to_mm9.sh",
      sep = " ", append = T)
  cat(c("gzip ", aafiles_full_2[i], "\n"),
      file = "~/Documents/Shayan/BioInf/lncRNA/oRNAment/liftover_to_mm9.sh",
      sep = " ", append = T)
}


aa_mm9_rbp_files <- list.files("~/Documents/Shayan/BioInf/lncRNA/oRNAment/processed_mm9.nosync", pattern = ".bed", full.names = T)
aa_mm9_rbp_name2 <- list.files("~/Documents/Shayan/BioInf/lncRNA/oRNAment/processed_mm9.nosync", pattern = ".bed", full.names = F)
aaname2 <- unlist(lapply(strsplit(aa_mm9_rbp_name2, "\\."), "[[", 1))

# read bed files
aa_rbp <- list()
for(i in 1:length(aa_mm9_rbp_files)){
  aa_rbp[[i]] <- read.table(aa_mm9_rbp_files[i], stringsAsFactors = F)
}
mESC_RBP_GR_list <- lapply(aa_rbp,
                                      makeGRangesFromDataFrame,
                                      keep.extra.columns = T,
                                      seqnames.field = "V1",
                                      start.field = "V2",
                                      end.field = "V3")
names(mESC_RBP_GR_list) <- aaname2
mESC_RBP_tiled_list <- list()
mESC_RBP_tiled_list_GR <- list()

for(i in 1:length(mESC_RBP_GR_list)){
  aaov <- findOverlaps(MM9_1kb_tiled_GR, mESC_RBP_GR_list[[i]])
  aadf <- data.frame(tile_index = aaov@from,
                     rbp_binding_index = aaov@to, 
                     interaction_index = c(1:length(aaov@to)))
  aadf_agg <- aggregate(aadf[c("interaction_index")],
                        by = aadf[c("tile_index")],
                        FUN = c)
  aa_nu_int <- unlist(lapply(aadf_agg$interaction_index, length))
  mESC_RBP_tiled_list[[i]] <- cbind(MM9_1kb_tiled[aadf_agg$tile_index,], aa_nu_int)
  mESC_RBP_tiled_list_GR[[i]] <- makeGRangesFromDataFrame(mESC_RBP_tiled_list[[i]], 
                                                          keep.extra.columns = T)
  
}
names(mESC_RBP_tiled_list) <- names(mESC_RBP_GR_list)
names(mESC_RBP_tiled_list_GR) <- names(mESC_RBP_GR_list)

summary(mESC_RBP_tiled_list$A1CF_mm9$aa_nu_int)


RBP_features <- matrix(nrow = length(aaMM9_1kb_tiled_GR_filtered),
                       ncol = length(mESC_RBP_tiled_list))

rownames(RBP_features) <- MM9_1kb_tiled_GR_filtered$tile
colnames(RBP_features) <- names(mESC_RBP_tiled_list)
for(i in 1:length(mESC_RBP_tiled_list)){
  print(i)
  aatmp <- mESC_RBP_tiled_list[[i]]
  aatmp <- aatmp[aatmp$tile %in% rownames(RBP_features),]
  RBP_features[match(aatmp$tile, rownames(RBP_features)),i] <- aatmp$aa_nu_int
}
save(list= c("RBP_features"), file = "~/Documents/Shayan/BioInf/lncRNA/RBP_features.RData")

# for each RBP get the mean number of predicted binding sites in + vs - examples
aa_pos_mean <- numeric(ncol(RBP_features))
aa_neg_mean <- numeric(ncol(RBP_features))
names(aa_pos_mean) <- colnames(RBP_features)
names(aa_neg_mean) <- colnames(RBP_features)

aa_pos_sd <- numeric(ncol(RBP_features))
aa_neg_sd <- numeric(ncol(RBP_features))
names(aa_pos_sd) <- colnames(RBP_features)
names(aa_neg_sd) <- colnames(RBP_features)


aaps <- which(MM9_1kb_tiled_owner_labels_binary_3Mfilter == 1)
aang <- which(MM9_1kb_tiled_owner_labels_binary_3Mfilter == 0)
for(i in 1:ncol(RBP_features)){
  print(i)
  aatmp <- RBP_features[, i]
  aatmp[is.na(aatmp)] <- 0
  aa_pos_mean[i] <- mean(aatmp[aaps])
  aa_pos_sd[i] <- sd(aatmp[aaps])
  aa_neg_mean[i] <- mean(aatmp[aang])
  aa_neg_sd[i] <- sd(aatmp[aang])
}

barplot(rbind(aa_pos_mean, aa_neg_mean), beside = T, las= 2)
hist(aa_pos_mean-aa_neg_mean , main = " (mean count in + examples) - (mean count in - examples) ")
sum((aa_pos_mean   - aa_pos_sd )> aa_neg_mean)


mESC_RBP_lncRNA_features <- matrix(nrow = length(lncRNA_chosen_gt1k_uniqTiles_GR), ncol = length(mESC_RBP_GR_list))
rownames(mESC_RBP_lncRNA_features) <- lncRNA_chosen_gt1k_uniqTiles_GR$gene_name
colnames(mESC_RBP_lncRNA_features) <- names(mESC_RBP_tiled_list)

for(i in 1:length(mESC_RBP_GR_list)){
  print(i)
  aaov <- findOverlaps(lncRNA_chosen_gt1k_uniqTiles_GR, mESC_RBP_GR_list[[i]])
  aadf <- data.frame(lncRNA_index = aaov@from,
                     rbp_binding_index = aaov@to, 
                     interaction_index = c(1:length(aaov@to)))
  aadf_agg <- aggregate(aadf[c("interaction_index")],
                        by = aadf[c("lncRNA_index")],
                        FUN = c, simplify=F)
  aa_nu_int <- unlist(lapply(aadf_agg$interaction_index, length))
  aa_my_ind <- match(aadf_agg$lncRNA_index, c(1:length(lncRNA_chosen_gt1k_uniqTiles_GR)))
  
  mESC_RBP_lncRNA_features[aa_my_ind,i] <- aa_nu_int
  
}
boxplot.matrix(t(mESC_RBP_lncRNA_features), las = 2, outline= F)
save(list = c("mESC_RBP_lncRNA_features"), 
     file ="~/Documents/Shayan/BioInf/lncRNA/mESC_RBP_lncRNA_features.RData")
##############################################################################################################################################################################################
# recalculating RBP features from scanned RBP files 

aa_mm9_rbp_files <- list.files("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/Mouse_RBP_beds", pattern = "*.bed", full.names = T)
aa_mm9_rbp_name2 <- list.files("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/Mouse_RBP_beds/", pattern = "*.bed", full.names = F)
aaname2 <- unlist(lapply(strsplit(aa_mm9_rbp_name2, "_mm9"), "[[", 1))

# read bed files
aa_p6_gr <- MM9_1kb_tiled_GR_filtered[match(Partition_6_dfs$tile_name, MM9_1kb_tiled_GR_filtered$tile)]
aa_p6_RBP <- matrix(0L, nrow = length(aa_p6_gr), ncol = length(aaname2))
rownames(aa_p6_RBP) <- aa_p6_gr$tile
colnames(aa_p6_RBP) <- aaname2
for(i in 28:length(aa_mm9_rbp_files)){
  print(i)
  aa_tmp <- read.table(aa_mm9_rbp_files[i], stringsAsFactors = F)
  aa_tGR <- makeGRangesFromDataFrame(aa_tmp,
                                     seqnames.field = "V1",
                                      start.field = "V2",
                                      end.field = "V3")
  aaov <- findOverlaps(aa_p6_gr, aa_tGR)
  aadf <- data.frame(tile_index = aaov@from,
                     rbp_binding_index = aaov@to, 
                     interaction_index = c(1:length(aaov@to)))
  aadf_agg <- aggregate(aadf[c("interaction_index")],
                        by = aadf[c("tile_index")],
                        FUN = length)
  aa_p6_RBP[match(aa_p6_gr$tile[aadf_agg$tile_index], rownames(aa_p6_RBP)), i] <- aadf_agg$interaction_index
}
Partition_6_RBP_scanned <- aa_p6_RBP
save(list = c("Partition_6_RBP_scanned"), file = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_RBP_scanned.RData")

##############################################################################################################################################################################################
##############################################################################################################################################################################################
save(list=c("Partition_1_dfs", "MM9_1kb_tiled", "lncRNA_chosen_gt1k_uniqTiles"), file = "~/Documents/Shayan/BioInf/lncRNA/calc_dist.RData")

# distance
# aa_distance <- numeric(length = nrow(Partition_1_dfs))
# for(i in 1:nrow(Partition_1_dfs)){
#   print(i)
#   aastdis <- abs(MM9_1kb_tiled$start[MM9_1kb_tiled$tile == Partition_1_dfs$tile_name[i]] - lncRNA_chosen_gt1k_uniqTiles$start_position[lncRNA_chosen_gt1k_uniqTiles$gene_name == Partition_1_dfs$owner[i]])
#   aaendis <- abs(MM9_1kb_tiled$start[MM9_1kb_tiled$tile == Partition_1_dfs$tile_name[i]] - lncRNA_chosen_gt1k_uniqTiles$end_position[lncRNA_chosen_gt1k_uniqTiles$gene_name == Partition_1_dfs$owner[i]])
#   aa_distance[i] <- min(aastdis, aaendis)
# }
# distance_feature_partition_1 <- aa_distance
# 


aa_distance <- numeric(length = nrow(Partition_1_dfs))
aalnc <- levels((Partition_1_dfs$owner))
for(i in 1:length(aalnc)){
  aa_cur_qu <- MM9_1kb_tiled_GR_filtered[match(Partition_1_dfs$tile_name[Partition_1_dfs$owner == aalnc[i]], MM9_1kb_tiled_GR_filtered$tile)]
  aa_cur_su <- lncRNA_chosen_gt1k_uniqTiles_GR[lncRNA_chosen_gt1k_uniqTiles_GR$gene_name == aalnc[i]]
  aadist <- GenomicRanges::distance(x = aa_cur_qu, y = aa_cur_su, select = "all")
  aa_distance[match(aa_cur_qu$tile, Partition_1_dfs$tile_name)] <- aadist
}

distance_feature_partition_1 <- aa_distance
names(distance_feature_partition_1) <- Partition_1_dfs$tile_name

aa_distance <- numeric(length = nrow(Partition_2_dfs))
aalnc <- levels((Partition_2_dfs$owner))
for(i in 1:length(aalnc)){
  aa_cur_qu <- MM9_1kb_tiled_GR_filtered[match(Partition_2_dfs$tile_name[Partition_2_dfs$owner == aalnc[i]], MM9_1kb_tiled_GR_filtered$tile)]
  aa_cur_su <- lncRNA_chosen_gt1k_uniqTiles_GR[lncRNA_chosen_gt1k_uniqTiles_GR$gene_name == aalnc[i]]
  aadist <- GenomicRanges::distance(x = aa_cur_qu, y = aa_cur_su, select = "all")
  aa_distance[match(aa_cur_qu$tile, Partition_2_dfs$tile_name)] <- aadist
}

distance_feature_partition_2 <- aa_distance
names(distance_feature_partition_2) <- Partition_2_dfs$tile_name
save(list  = c("distance_feature_partition_2"), 
     file = "~/Documents/Shayan/BioInf/lncRNA/distance_to_closest_lncRNA_partition2.RData")

aa_distance <- numeric(length = nrow(Partition_4_dfs))
aalnc <- levels((Partition_4_dfs$owner))
for(i in 1:length(aalnc)){
  aa_cur_qu <- MM9_1kb_tiled_GR_filtered[match(Partition_4_dfs$tile_name[Partition_4_dfs$owner == aalnc[i]], MM9_1kb_tiled_GR_filtered$tile)]
  aa_cur_su <- lncRNA_chosen_gt1k_uniqTiles_GR[lncRNA_chosen_gt1k_uniqTiles_GR$gene_name == aalnc[i]]
  aadist <- GenomicRanges::distance(x = aa_cur_qu, y = aa_cur_su, select = "all")
  aa_distance[match(aa_cur_qu$tile, Partition_4_dfs$tile_name)] <- aadist
}

distance_feature_partition_4 <- aa_distance
names(distance_feature_partition_4) <- Partition_4_dfs$tile_name
save(list  = c("distance_feature_partition_4"), 
     file = "~/Documents/Shayan/BioInf/lncRNA/distance_to_closest_lncRNA_partition4.RData")

aa_distance <- numeric(length = nrow(Partition_5_dfs))
aalnc <- levels((Partition_5_dfs$owner))
for(i in 1:length(aalnc)){
  aa_cur_qu <- MM9_1kb_tiled_GR_filtered[match(Partition_5_dfs$tile_name[Partition_5_dfs$owner == aalnc[i]], MM9_1kb_tiled_GR_filtered$tile)]
  aa_cur_su <- lncRNA_chosen_gt1k_uniqTiles_GR[lncRNA_chosen_gt1k_uniqTiles_GR$gene_name == aalnc[i]]
  aadist <- GenomicRanges::distance(x = aa_cur_qu, y = aa_cur_su, select = "all")
  aa_distance[match(aa_cur_qu$tile, Partition_5_dfs$tile_name)] <- aadist
}

distance_feature_partition_5 <- aa_distance
names(distance_feature_partition_5) <- Partition_5_dfs$tile_name
save(list  = c("distance_feature_partition_5"), 
     file = "~/Documents/Shayan/BioInf/lncRNA/distance_to_closest_lncRNA_partition5.RData")
##########
aa_distance <- numeric(length = nrow(Partition_6_dfs))
aalnc <- levels((Partition_6_dfs$owner))
aamydf_list <- list()
for(i in 1:length(aalnc)){
  print(aalnc[i])
  aamytil <- Partition_6_dfs$tile_name[Partition_6_dfs$owner == aalnc[i]]
  aa_cur_qu <- MM9_1kb_tiled_GR_filtered[match(aamytil, MM9_1kb_tiled_GR_filtered$tile)]
  aa_cur_su <- lncRNA_chosen_gt1k_uniqTiles_GR[lncRNA_chosen_gt1k_uniqTiles_GR$gene_name == aalnc[i]]
  aadist <- GenomicRanges::distance(x = aa_cur_qu, y = aa_cur_su, select = "all")
  print(sum(is.na(aadist)))
  aamydf_list[[i]] <- data.frame(tile_name = aamytil, 
                       owner = rep(aalnc[i], sum(Partition_6_dfs$owner == aalnc[i])),
                       distance = aadist[match(aa_cur_qu$tile,aamytil)])
 # aa_distance[match(aa_cur_qu$tile, Partition_6_dfs$tile_name)] <- aadist
}

aamydf_list_all <- do.call(rbind, aamydf_list)
aadistall <- join(Partition_6_dfs, aamydf_list_all)
all(aadistall$tile_name == Partition_6_dfs$tile_name)
aadistall$distance[is.na(aadistall$distance)] <- max(aadistall$distance, na.rm = T) + 1000000
distance_feature_partition_6_df <- aadistall
distance_feature_partition_6 <- aadistall$distance
names(distance_feature_partition_6) <- Partition_6_dfs$tile_name
distance_feature_partition_6_cis <- distance_feature_partition_6
distance_feature_partition_6_cis[distance_feature_partition_6_cis == 107540001] <- NA
#distance_feature_partition_6[is.na(distance_feature_partition_6)] <- max(distance_feature_partition_6, na.rm = T) + 1000000
save(list  = c("distance_feature_partition_6"), 
     file = "~/Documents/Shayan/BioInf/lncRNA/distance_to_owner_partition6.RData")
par(mfrow = c(1,1), mar = c(10,4,4,4))
boxplot(distance_feature_partition_6 ~ Partition_6_random_chunk_cv_df$owner, las = 2)
boxplot(log10(distance_feature_partition_6) ~ Partition_6_random_chunk_cv_df$owner, las = 2)
sum(aadistall$distance == 0, na.rm = T)
# add lower resolution version of distance
distance_feature_partition_6_lowRes <- distance_feature_partition_6
distance_feature_partition_6_lowRes <- round(distance_feature_partition_6_lowRes/100000)
boxplot(distance_feature_partition_6_lowRes ~ Partition_6_random_chunk_cv_df$owner, las = 2)

distance_feature_partition_6_lowRes_df <- aadistall
distance_feature_partition_6_lowRes_df$distance <- round(distance_feature_partition_6_lowRes_df$distance/100000)
distance_feature_partition_6_lowRes_df$highres <- distance_feature_partition_6
distance_feature_partition_6_lowRes_df$lowerres1Mg <- round(distance_feature_partition_6/1000000)
distance_feature_partition_6_lowRes_df$lowerres10Mg <- round(distance_feature_partition_6/10000000)


aadistance_feature_partition_6_lowRes_df <- distance_feature_partition_6_lowRes_df
aadistance_feature_partition_6_lowRes_df$owner <- factor(aadistance_feature_partition_6_lowRes_df$owner)
aadistance_feature_partition_6_lowRes_df$label <- factor(aadistance_feature_partition_6_lowRes_df$label)
aadistance_feature_partition_6_lowRes_df$highres <- distance_feature_partition_6
aadistance_feature_partition_6_lowRes_df$lowerres1Mg <- round(distance_feature_partition_6/1000000)
aadistance_feature_partition_6_lowRes_df$lowerres10Mg <- round(distance_feature_partition_6/10000000)

save(list  = c("distance_feature_partition_6_lowRes"), 
     file = "~/Documents/Shayan/BioInf/lncRNA/distance_feature_partition_6_lowRes.RData")
save(list  = c("distance_feature_partition_6_lowRes_df"), 
     file = "~/Documents/Shayan/BioInf/lncRNA/distance_feature_partition_6_lowRes_df.RData")

all(Partition_6_dfs_GR$tile == distance_feature_partition_6_lowRes_df$tile_name)
distance_feature_partition_6_lowRes_dfGR <- Partition_6_dfs_GR
distance_feature_partition_6_lowRes_dfGR$highres <- distance_feature_partition_6_lowRes_df$highres
distance_feature_partition_6_lowRes_dfGR$lowerres1Mg <- distance_feature_partition_6_lowRes_df$lowerres1Mg
distance_feature_partition_6_lowRes_dfGR$lowerres10Mg <- distance_feature_partition_6_lowRes_df$lowerres10Mg
save(list  = c("distance_feature_partition_6_lowRes_dfGR"), 
     file = "~/Documents/Shayan/BioInf/lncRNA/distance_feature_partition_6_lowRes_dfGRs.RData")

ggplot(aadistance_feature_partition_6_lowRes_df, aes(x=owner, y=distance, fill=label)) + 
  geom_boxplot() +
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  # ylab("delta auprc chunk")+
  # xlab("lncRNA") + 
  # #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))
ggplot(aadistance_feature_partition_6_lowRes_df, aes(x=owner, y=highres, fill=label)) + 
  geom_boxplot() +
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  # ylab("delta auprc chunk")+
  # xlab("lncRNA") + 
  # #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

ggplot(aadistance_feature_partition_6_lowRes_df, aes(x=owner, y=lowerres1Mg, fill=label)) + 
  geom_boxplot() +
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  # ylab("delta auprc chunk")+
  # xlab("lncRNA") + 
  # #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

ggplot(aadistance_feature_partition_6_lowRes_df, aes(x=owner, y=lowerres10Mg, fill=label)) + 
  geom_boxplot() +
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  # ylab("delta auprc chunk")+
  # xlab("lncRNA") + 
  # #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

ggplot(aadistance_feature_partition_6_lowRes_df, aes(x=owner, y=round(log10(highres)), fill=label)) + 
  geom_boxplot() +
  #ggtitle() + # \npaired t-test p-value: ",
  #format(aattest$p.value, scientific = T, digits = 3))) +
  # ylab("delta auprc chunk")+
  # xlab("lncRNA") + 
  # #  ylim(0, 0.2)+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 90),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

##############################

aanam <- unique(unlist(lapply(strsplit(sort(unique(MM9_1kb_tiled_chrname)), "_"), "[[", 1)))
cat(aanam, sep = " ")
cat(c("#!/bin/bash\n"), sep = " ", file = "move_kmers.sh", append = F)
for(i in length(aanam):1){
  aa_in <- paste0("'",aanam[i],"_*","'")
    #paste0("/shared-mounts/sinhas/tabebor2/lncRNA/kmer_by_Tile/", aanam[i],"_*")
  aa_out <- paste0("/shared-mounts/sinhas/tabebor2/lncRNA/tile_kmer_byChr/",  aanam[i], "/")
  cat(c("find /shared-mounts/sinhas/tabebor2/lncRNA/kmer_by_Tile -maxdepth 1 -type f -name",  aa_in ,"| xargs", "mv -t", aa_out, '"{}"\n'), sep = " ", file = "move_kmers.sh", append = T)
  cat(c("echo ", aanam[i], "done\n"), sep = " ", file = "move_kmers.sh", append = T)
}

# find /shared-mounts/sinhas/tabebor2/lncRNA/kmer_by_Tile -maxdepth 1 -type f -name '??????????_a1ac*.5mer'  -exec mv -t destination "{}" +

#rewriting Jellyfish Jobs to store results in directories for each chromosome
cat(c("#!/bin/bash\n"), sep = " ", file = "modify_Jellyjob.sh", append = F)
for(i in 1:length(aanam)){
  aa_in <- paste0("'",aanam[i],"_*","'")
  cat(c(paste0("sed -i 's:kmer_by_Tile/", aanam[i], "_:kmer_by_chr_byTile/", aanam[i], "/", aanam[i], "_:g'"), "Jellyfish_mm9.job\n"),
      sep = " ", 
      file = "modify_Jellyjob.sh", append = T)
  cat(c("echo ", aanam[i], "done\n"), sep = " ", file = "modify_Jellyjob.sh", append = T)
}

#sed -i 's:kmer_by_Tile/chr1_:kmer_by_chr_byTile/chr1/chr1_:g' aahead

#############################################################################
# creating a feature that counts the number of RBPs common to a tile and its owner
load("RBP_features.RData")
#RBP_features
#3mESC_RBP_lncRNA_features
all(colnames(RBP_features) == colnames(mESC_RBP_lncRNA_features))
aa_rbp_pair_intersect_list <- list()
aa_RBP_features_noNA <- RBP_features[match(names(MM9_1kb_tiled_owner_3Mfilter_noNA), rownames(RBP_features)),]



# creating a feature that counts the number of ChIPs (from chip atlas) common to a tile and its owner
load("ChIPATLAS_features.RData")
#ChIPATLAS_features
#mESC_CHIPATLAS_lncRNA_features
all(colnames(ChIPATLAS_features) == colnames(mESC_CHIPATLAS_lncRNA_features))
aa_chrom <- grep(pattern = '^H[0-9]{1,2}', x = colnames(ChIPATLAS_features))
aa_chip_other <- setdiff(c(1:ncol(ChIPATLAS_features)), aa_chrom)

aa_chrom_features_tile  <- ChIPATLAS_features[, aa_chrom]
aa_chip_features_tile  <- ChIPATLAS_features[, aa_chip_other]
aa_chrom_features_tile_noNA <- aa_chrom_features_tile[match(names(MM9_1kb_tiled_owner_3Mfilter_noNA), rownames(aa_chrom_features_tile)),]
aa_chip_features_tile_noNA <- aa_chip_features_tile[match(names(MM9_1kb_tiled_owner_3Mfilter_noNA), rownames(aa_chip_features_tile)),]

aa_chrom_features_lncRNA  <- mESC_CHIPATLAS_lncRNA_features[, aa_chrom]
aa_chip_features_lncRNA  <- mESC_CHIPATLAS_lncRNA_features[, aa_chip_other]

aa_chrom_pair_intersect_list <- list()
aa_chip_pair_intersect_list <- list()

for(i in 1:length(MM9_1kb_tiled_owner_3Mfilter_noNA)){
  print(i)
  aa_rbp_pair_intersect_list[[i]] <- colnames(aa_RBP_features_noNA)[intersect(which(aa_RBP_features_noNA[i, ] > 0),
                                                                              which(mESC_RBP_lncRNA_features[rownames(mESC_RBP_lncRNA_features) %in% MM9_1kb_tiled_owner_3Mfilter_noNA[i],] > 0))]
  aa_chrom_pair_intersect_list[[i]] <- colnames(aa_chrom_features_lncRNA)[intersect(which(aa_chrom_features_tile_noNA[i, ] > 0),
                                                                              which(aa_chrom_features_lncRNA[rownames(aa_chrom_features_lncRNA) %in% MM9_1kb_tiled_owner_3Mfilter_noNA[i],] > 0))]
  aa_chip_pair_intersect_list[[i]] <- colnames(aa_chip_features_tile_noNA)[intersect(which(aa_chip_features_tile_noNA[i, ] > 0),
                                                                              which(aa_chip_features_lncRNA[rownames(aa_chip_features_lncRNA) %in% MM9_1kb_tiled_owner_3Mfilter_noNA[i],] > 0))]
  
  
}
RBP_tile_owner_pair_feature <- unlist(lapply(aa_rbp_pair_intersect_list, length))
ChIP_tile_owner_pair_feature <- unlist(lapply(aa_chip_pair_intersect_list, length))
Chromatin_tile_owner_pair_feature <- unlist(lapply(aa_chrom_pair_intersect_list, length))
names(RBP_tile_owner_pair_feature) <- names(MM9_1kb_tiled_owner_3Mfilter_noNA)
names(ChIP_tile_owner_pair_feature) <- names(MM9_1kb_tiled_owner_3Mfilter_noNA)
names(Chromatin_tile_owner_pair_feature) <- names(MM9_1kb_tiled_owner_3Mfilter_noNA)

save(list = c("RBP_tile_owner_pair_feature", "ChIP_tile_owner_pair_feature", "Chromatin_tile_owner_pair_feature"),
     file = "~/Documents/Shayan/BioInf/lncRNA/Tile_owner_pair_RBP_ChIP_pair.RData")

boxplot(RBP_tile_owner_pair_feature~MM9_1kb_tiled_owner_labels_binary_3Mfilter_noNA)
boxplot(ChIP_tile_owner_pair_feature~MM9_1kb_tiled_owner_labels_binary_3Mfilter_noNA)
boxplot(Chromatin_tile_owner_pair_feature~MM9_1kb_tiled_owner_labels_binary_3Mfilter_noNA)

t.test(RBP_tile_owner_pair_feature[MM9_1kb_tiled_owner_labels_binary_3Mfilter_noNA == 0], RBP_tile_owner_pair_feature[MM9_1kb_tiled_owner_labels_binary_3Mfilter_noNA == 1])
t.test(ChIP_tile_owner_pair_feature[MM9_1kb_tiled_owner_labels_binary_3Mfilter_noNA == 0], ChIP_tile_owner_pair_feature[MM9_1kb_tiled_owner_labels_binary_3Mfilter_noNA == 1])
t.test(Chromatin_tile_owner_pair_feature[MM9_1kb_tiled_owner_labels_binary_3Mfilter_noNA == 0], Chromatin_tile_owner_pair_feature[MM9_1kb_tiled_owner_labels_binary_3Mfilter_noNA == 1])


aa_perlnc_rbp_pair <- list()
aa_perlnc_chrom_pair <- list()
aa_perlnc_chip_pair <- list()
for(i in 1:nrow(aa_chip_features_lncRNA)){
  print(i)
  aa_perlnc_rbp_pair[[i]] <- sort(table(unlist(aa_rbp_pair_intersect_list[MM9_1kb_tiled_owner_3Mfilter_noNA %in% rownames(aa_chip_features_lncRNA)[i]])), decreasing = T)
  aa_perlnc_chrom_pair[[i]] <- sort(table(unlist(aa_chrom_pair_intersect_list[MM9_1kb_tiled_owner_3Mfilter_noNA %in% rownames(aa_chip_features_lncRNA)[i]])), decreasing = T)
  aa_perlnc_chip_pair[[i]] <- sort(table(unlist(aa_chip_pair_intersect_list[MM9_1kb_tiled_owner_3Mfilter_noNA %in% rownames(aa_chip_features_lncRNA)[i]])), decreasing = T)
  
}
names(aa_perlnc_rbp_pair) <- rownames(aa_chip_features_lncRNA)
names(aa_perlnc_chrom_pair) <- rownames(aa_chip_features_lncRNA)
names(aa_perlnc_chip_pair) <- rownames(aa_chip_features_lncRNA)

#############################################################################
# creating a feature that counts the number of repeat elements per tile




#############################################################################
# Check the coverage of RBP binding sites for a bunch of RBPs
# write scripts using bedtools genomecov

cat(c("#!/bin/bash\n"), file = "~/Documents/Shayan/BioInf/lncRNA/oRNAment/get_RBP_genomeCoverage.sh", append = F)
aa_files <- list.files("~/Documents/Shayan/BioInf/lncRNA/oRNAment/processed_mm9/")
aa_files <- aa_files[-grep(pattern = "unconverted*", x = aa_files)]
aa_files_rbp <- unlist(lapply(strsplit(aa_files, "_mm9"), "[[", 1))
for(i in 1:length(aa_files)){
  aa_inp <- paste0("<((sort -k 1,1 ~/Documents/Shayan/BioInf/lncRNA/oRNAment/processed_mm9/", aa_files[i], "))")
  cat(c("bedtools genomecov -i ", aa_inp,
         " -g ~/Documents/Shayan/BioInf/lncRNA/mm9_chr_sizes.txt > ~/Documents/Shayan/BioInf/lncRNA/oRNAment/coverage/",
        aa_files_rbp[i], "_coverage.txt\n" ), 
      sep = "", file = "~/Documents/Shayan/BioInf/lncRNA/oRNAment/get_RBP_genomeCoverage.sh", append = T)
}


############################
#Add u1snRNA binding site info
Partition_6_dfs$tile_name[1:5]
Partition_6_dfs_GR <- MM9_1kb_tiled_GR[match(Partition_6_dfs$tile_name, MM9_1kb_tiled$tile),]
save(list = c("Partition_6_dfs_GR"), file = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_dfs_GR.RData")
aa_p6_seq <- getSeq(x = BSgenome.Mmusculus.UCSC.mm9, names = Partition_6_dfs_GR)

# aa_fw <- grep(pattern = "TACTTACC", aa_p6_seq)
# aa_rv <- grep(pattern = "GGTAAGTA", aa_p6_seq)

aa_fw <- vcountPattern(pattern = "TACTTACC", subject = aa_p6_seq,
                       max.mismatch=1, min.mismatch=0,
                       with.indels=FALSE, fixed=TRUE,
                       algorithm="auto")
aa_rv <- vcountPattern(pattern = "GGTAAGTA", subject = aa_p6_seq,
                       max.mismatch=1, min.mismatch=0,
                       with.indels=FALSE, fixed=TRUE,
                       algorithm="auto")
aa_both <- apply(rbind(aa_fw, aa_rv), 2, max)

Partition_6_U1snRNA <- aa_both
names(Partition_6_U1snRNA) <- Partition_6_dfs_GR$tile

save(list = c("Partition_6_U1snRNA"), file = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_U1snRNA.RData")
################################################################################################################
#creating TF motif features
# creating RBP motif features
# motifs for RBP and TFs for both mouse and human from Cisbp database
# choose TFs for mouse
# TF pwms mouse
aa_tf_mouse_info <- read.delim(file = "~/Documents/Shayan/BioInf/PWMs/Mus_musculus_2020_11_24_7-38_am/TF_Information.txt", header = T,sep = "\t" )
aa_tf_mouse_info2 <- aa_tf_mouse_info[aa_tf_mouse_info$TF_Status == "D",]
TF_mouse_info <- aa_tf_mouse_info2
aaunq <- unique(TF_mouse_info$TF_Name)
aa_new_unq_name <- character(length = length(TF_mouse_info$TF_Name))
for(i in 1:length(aaunq)){
  aa_which <- which(TF_mouse_info$TF_Name == aaunq[i])
  aa_new_unq_name[aa_which] <- paste0(aaunq[i], "__", c(1:length(aa_which)))
}
TF_mouse_info$unique_name <- aa_new_unq_name
TF_mouse_list <- list()
for(i in 1:nrow(TF_mouse_info)){
  TF_mouse_list[[i]] <- t(read.delim(paste0("~/Documents/Shayan/BioInf/PWMs/Mus_musculus_2020_11_24_7-38_am/pwms_all_motifs/",
                                            TF_mouse_info$Motif_ID[i], ".txt"),
                                     header = T, row.names = 1))
}
names(TF_mouse_list) <- TF_mouse_info$unique_name
aancol <- unlist(lapply(TF_mouse_list, ncol))
TF_mouse_list <- TF_mouse_list[-which(aancol == 0)]
save(list = c("TF_mouse_list"), file = "~/Documents/Shayan/BioInf/PWMs/TF_mouse_list_pwm.RData")
my_pwm_list <- TF_mouse_list
save(list = c("my_pwm_list"), file = "~/Documents/Shayan/BioInf/PWMs/my_pwm_list_mouse_TF.RData")

my_pwm_list <- TF_mouse_list[1:2]
save(list = c("my_pwm_list"), file = "~/Documents/Shayan/BioInf/PWMs/my_pwm_list_mouse_TF_test.RData")

for(i in 1:length(TF_mouse_list)){
  cat(c("> ", names(TF_mouse_list)[i], "\n"), sep = "", append = F,
      file = paste0("~/Documents/Shayan/BioInf/PWMs/TF_mouse_PFM/", names(TF_mouse_list)[i], ".lpm"))
  
  write.table(t(TF_mouse_list[[i]]), append = T, quote = F, row.names = F, col.names = F,
              file = paste0("~/Documents/Shayan/BioInf/PWMs/TF_mouse_PFM/", names(TF_mouse_list)[i], ".lpm"))
}

#write job file for motif conversion
cat(c("#!/bin/bash\n"), sep = "", append = F, file = "~/Documents/Shayan/BioInf/PWMs/lpm_to_mat_mouse_TF.sh")
for(i in 1:length(TF_mouse_list)){
  cat(c("/shared-mounts/sinhas/tabebor2/pwmscan/perl_tools/lpmconvert.pl ", "-o /shared-mounts/sinhas/tabebor2/pwmscan/pwm_TF_mouse/", names(TF_mouse_list)[i], ".mat ", "/shared-mounts/sinhas/tabebor2/pwmscan/pwm_TF_mouse/",
        names(TF_mouse_list)[i], ".lpm\n"),
      sep = "", append = T, file = "~/Documents/Shayan/BioInf/PWMs/lpm_to_mat_mouse_TF.sh" )
}
#write job file for scanning mm9 genome ucsc style using motifScanner

for(i in 1:length(TF_mouse_list)){
  cat(c("/shared-mounts/sinhas/tabebor2/pwmscan/pwm_mscan_wrapper_ucsc ", 
        "-m /shared-mounts/sinhas/tabebor2/pwmscan/pwm_TF_mouse/", names(TF_mouse_list)[i], 
        ".mat -e 0.0001 -d /shared-mounts/sinhas/tabebor2/pwmscan/genomedb -s mm9 > /shared-mounts/sinhas/tabebor2/pwmscan/bed_TF_mouse/", names(TF_mouse_list)[i], "____mm9_hits.bed\n"),
      sep = "", append = T, file = "~/Documents/Shayan/BioInf/PWMs/scan_mm9_TF.job" )
}

#pwm_mscan_wrapper_ucsc
#pwm_mscan_wrapper -m <matrix-file> -e <p-value> -d <genome-root-dir> -s <assembly[hg19|mm9|..]> [options] -w 
################################################################################################################
# RBP pwms mouse
aa_rbp_mouse_info <- read.delim(file = "~/Documents/Shayan/BioInf/PWMs/Mus_musculus_2020_11_24_8-02_am_RBP/RBP_Information_all_motifs.txt", header = T,sep = "\t" )
aa_rbp_mouse_info2 <- aa_rbp_mouse_info[(aa_rbp_mouse_info$Motif_ID != "."),]
aa_rbp_mouse_info2 <- aa_rbp_mouse_info2[!duplicated(aa_rbp_mouse_info2$Motif_ID),]
aa_file_list <- list.files("~/Documents/Shayan/BioInf/PWMs/Mus_musculus_2020_11_24_8-02_am_RBP/pwms_all_motifs", full.names = T)
aanames <- list.files("~/Documents/Shayan/BioInf/PWMs/Mus_musculus_2020_11_24_8-02_am_RBP/pwms_all_motifs", full.names = F)
aanames2 <- unlist(lapply(strsplit(aanames, ".txt"), "[[", 1))
aax <- numeric(length = length(aanames2))

for(i in 1:length(aa_file_list)){
  aax[i] <- file.size(aa_file_list[i])
  
}
names(aax) <- aanames2
aanames22 <- aanames2[aax == 0]
aa_rbp_mouse_info2 <- aa_rbp_mouse_info2[!(aa_rbp_mouse_info2$Motif_ID %in% aanames22),]

RBP_mouse_info <- aa_rbp_mouse_info2
aaunq <- unique(RBP_mouse_info$RBP_Name)
aa_new_unq_name <- character(length = length(RBP_mouse_info$RBP_Name))
for(i in 1:length(aaunq)){
  aa_which <- which(RBP_mouse_info$RBP_Name == aaunq[i])
  aa_new_unq_name[aa_which] <- paste0(aaunq[i], "__", c(1:length(aa_which)))
}
RBP_mouse_info$unique_name <- aa_new_unq_name

RBP_mouse_list <- list()
for(i in 1:nrow(RBP_mouse_info)){
  
  RBP_mouse_list[[i]] <- t(read.delim(file = paste0("~/Documents/Shayan/BioInf/PWMs/Mus_musculus_2020_11_24_8-02_am_RBP/pwms_all_motifs/",
                                             RBP_mouse_info$Motif_ID[i], ".txt"),
                                     header = T, row.names = 1))
}
names(RBP_mouse_list) <- RBP_mouse_info$unique_name 

aancol <- unlist(lapply(RBP_mouse_list, ncol))
table(aancol)
#RBP_mouse_list <- RBP_mouse_list[-which(aancol == 0)]

save(list = c("RBP_mouse_list"), file = "~/Documents/Shayan/BioInf/PWMs/RBP_mouse_list_pwm.RData")
my_pwm_list <- RBP_mouse_list
save(list = c("my_pwm_list"), file = "~/Documents/Shayan/BioInf/PWMs/my_pwm_list_mouse_RBP.RData")

for(i in 1:length(RBP_mouse_list)){
  cat(c("> ", names(RBP_mouse_list)[i], "\n"), sep = "", append = F,
      file = paste0("~/Documents/Shayan/BioInf/PWMs/RBP_mouse_PFM/", names(RBP_mouse_list)[i], ".lpm"))
  
  write.table(t(RBP_mouse_list[[i]]), append = T, quote = F, row.names = F, col.names = F,
              file = paste0("~/Documents/Shayan/BioInf/PWMs/RBP_mouse_PFM/", names(RBP_mouse_list)[i], ".lpm"))
}

#write job file for motif conversion
cat(c("#!/bin/bash\n"), sep = "", append = F, file = "~/Documents/Shayan/BioInf/PWMs/lpm_to_mat_mouse_RBP.sh")
for(i in 1:length(RBP_mouse_list)){
  cat(c("/shared-mounts/sinhas/tabebor2/pwmscan/perl_tools/lpmconvert.pl ", "-o /shared-mounts/sinhas/tabebor2/pwmscan/pwm_RBP_mouse/", names(RBP_mouse_list)[i], ".mat ", "/shared-mounts/sinhas/tabebor2/pwmscan/pwm_RBP_mouse/",
        names(RBP_mouse_list)[i], ".lpm\n"),
      sep = "", append = T, file = "~/Documents/Shayan/BioInf/PWMs/lpm_to_mat_mouse_RBP.sh" )
}
#write job file for scanning mm9 genome ucsc style using motifScanner

for(i in 1:length(RBP_mouse_list)){
  cat(c("/shared-mounts/sinhas/tabebor2/pwmscan/pwm_mscan_wrapper_ucsc ", 
        "-m /shared-mounts/sinhas/tabebor2/pwmscan/pwm_RBP_mouse/", names(RBP_mouse_list)[i], 
        ".mat -e 0.0001 -d /shared-mounts/sinhas/tabebor2/pwmscan/genomedb -s mm9 > /shared-mounts/sinhas/tabebor2/pwmscan/bed_RBP_mouse/", names(RBP_mouse_list)[i], "____mm9_hits.bed\n"),
      sep = "", append = T, file = "~/Documents/Shayan/BioInf/PWMs/scan_mm9_RBP.job" )
}


################################################################################################################
################################################################################################################
# TF pwms human
aa_tf_human_info <- read.delim(file = "~/Documents/Shayan/BioInf/PWMs/Homo_sapiens_2020_11_24_7-39_am/TF_Information.txt", header = T,sep = "\t" )
aa_tf_human_info2 <- aa_tf_human_info[aa_tf_human_info$TF_Status == "D",]

TF_human_info <- aa_tf_human_info2
aaunq <- unique(TF_human_info$TF_Name)
aa_new_unq_name <- character(length = length(TF_human_info$TF_Name))
for(i in 1:length(aaunq)){
  aa_which <- which(TF_human_info$TF_Name == aaunq[i])
  aa_new_unq_name[aa_which] <- paste0(aaunq[i], "__", c(1:length(aa_which)))
}
TF_human_info$unique_name <- aa_new_unq_name
TF_human_list <- list()
for(i in 1:nrow(TF_human_info)){
  print(i)
  TF_human_list[[i]] <- t(read.delim(paste0("~/Documents/Shayan/BioInf/PWMs/Homo_sapiens_2020_11_24_7-39_am/pwms_all_motifs/",
                                            TF_human_info$Motif_ID[i], ".txt"),
                                     header = T, row.names = 1))
}
names(TF_human_list) <- TF_human_info$unique_name

aancol <- unlist(lapply(TF_human_list, ncol))
table(aancol)
TF_human_list <- TF_human_list[-which(aancol == 0)]
TF_human_info <- TF_human_info[aancol > 0,]


save(list = c("TF_human_list"), file = "~/Documents/Shayan/BioInf/PWMs/TF_human_list_pwm.RData")
my_pwm_list <- TF_human_list
save(list = c("my_pwm_list"), file = "~/Documents/Shayan/BioInf/PWMs/my_pwm_list_human_TF.RData")

for(i in 1:length(TF_human_list)){
  cat(c("> ", names(TF_human_list)[i], "\n"), sep = "", append = F,
      file = paste0("~/Documents/Shayan/BioInf/PWMs/TF_human_PFM/", names(TF_human_list)[i], ".lpm"))
  
  write.table(t(TF_human_list[[i]]), append = T, quote = F, row.names = F, col.names = F,
              file = paste0("~/Documents/Shayan/BioInf/PWMs/TF_human_PFM/", names(TF_human_list)[i], ".lpm"))
}

#write job file for motif conversion
cat(c("#!/bin/bash\n"), sep = "", append = F, file = "~/Documents/Shayan/BioInf/PWMs/lpm_to_mat_human_TF.sh")
for(i in 1:length(TF_human_list)){
  cat(c("/shared-mounts/sinhas/tabebor2/pwmscan/perl_tools/lpmconvert.pl ", "-o /shared-mounts/sinhas/tabebor2/pwmscan/pwm_TF_human/", names(TF_human_list)[i], ".mat ", "/shared-mounts/sinhas/tabebor2/pwmscan/pwm_TF_human/",
        names(TF_human_list)[i], ".lpm\n"),
      sep = "", append = T, file = "~/Documents/Shayan/BioInf/PWMs/lpm_to_mat_human_TF.sh" )
}
#write job file for scanning mm9 genome ucsc style using motifScanner

for(i in 1:length(TF_human_list)){
  cat(c("/shared-mounts/sinhas/tabebor2/pwmscan/pwm_mscan_wrapper ", 
        "-m /shared-mounts/sinhas/tabebor2/pwmscan/pwm_TF_human/", names(TF_human_list)[i], 
        ".mat -e 0.0001 -d /shared-mounts/sinhas/tabebor2/pwmscan/genomedb -s hg19 > /shared-mounts/sinhas/tabebor2/pwmscan/bed_TF_human/", names(TF_human_list)[i], "____hg19_hits.bed\n"),
      sep = "", append = T, file = "~/Documents/Shayan/BioInf/PWMs/scan_hg19_TF.job" )
}

################################################################################################################
# RBP pwms human
aa_rbp_human_info <- read.delim(file = "~/Documents/Shayan/BioInf/PWMs/Homo_sapiens_2020_11_24_8-03_am_RBP/RBP_Information_all_motifs.txt", header = T,sep = "\t" )
aa_rbp_human_info <- aa_rbp_human_info[(aa_rbp_human_info$Motif_ID != "."),]
aa_rbp_human_info <- aa_rbp_human_info[!duplicated(aa_rbp_human_info$Motif_ID),]
aa_file_list <- list.files("~/Documents/Shayan/BioInf/PWMs/Homo_sapiens_2020_11_24_8-03_am_RBP/pwms_all_motifs", full.names = T)
aanames <- list.files("~/Documents/Shayan/BioInf/PWMs/Homo_sapiens_2020_11_24_8-03_am_RBP/pwms_all_motifs", full.names = F)
aanames2 <- unlist(lapply(strsplit(aanames, ".txt"), "[[", 1))
aax <- numeric(length = length(aanames2))

for(i in 1:length(aa_file_list)){
  aax[i] <- file.size(aa_file_list[i])
  
}
names(aax) <- aanames2
aanames22 <- aanames2[aax == 0]
aa_rbp_human_info <- aa_rbp_human_info[!(aa_rbp_human_info$Motif_ID %in% aanames22),]

RBP_human_info <- aa_rbp_human_info
aaunq <- unique(RBP_human_info$RBP_Name)
aa_new_unq_name <- character(length = length(RBP_human_info$RBP_Name))
for(i in 1:length(aaunq)){
  aa_which <- which(RBP_human_info$RBP_Name == aaunq[i])
  aa_new_unq_name[aa_which] <- paste0(aaunq[i], "__", c(1:length(aa_which)))
}
RBP_human_info$unique_name <- aa_new_unq_name

RBP_human_list <- list()
for(i in 1:nrow(RBP_human_info)){
  
  RBP_human_list[[i]] <- t(read.delim(file = paste0("~/Documents/Shayan/BioInf/PWMs/Homo_sapiens_2020_11_24_8-03_am_RBP/pwms_all_motifs/",
                                                    RBP_human_info$Motif_ID[i], ".txt"),
                                      header = T, row.names = 1))
}
names(RBP_human_list) <- RBP_human_info$unique_name 

aancol <- unlist(lapply(RBP_human_list, ncol))
table(aancol)
#RBP_human_list <- RBP_human_list[-which(aancol == 0)]

save(list = c("RBP_human_list"), file = "~/Documents/Shayan/BioInf/PWMs/RBP_human_list_pwm.RData")
my_pwm_list <- RBP_human_list
save(list = c("my_pwm_list"), file = "~/Documents/Shayan/BioInf/PWMs/my_pwm_list_human_RBP.RData")

for(i in 1:length(RBP_human_list)){
  cat(c("> ", names(RBP_human_list)[i], "\n"), sep = "", append = F,
      file = paste0("~/Documents/Shayan/BioInf/PWMs/RBP_human_PFM/", names(RBP_human_list)[i], ".lpm"))
  
  write.table(t(RBP_human_list[[i]]), append = T, quote = F, row.names = F, col.names = F,
              file = paste0("~/Documents/Shayan/BioInf/PWMs/RBP_human_PFM/", names(RBP_human_list)[i], ".lpm"))
}

#write job file for motif conversion
cat(c("#!/bin/bash\n"), sep = "", append = F, file = "~/Documents/Shayan/BioInf/PWMs/lpm_to_mat_human_RBP.sh")
for(i in 1:length(RBP_human_list)){
  cat(c("/shared-mounts/sinhas/tabebor2/pwmscan/perl_tools/lpmconvert.pl ", "-o /shared-mounts/sinhas/tabebor2/pwmscan/pwm_RBP_human/", names(RBP_human_list)[i], ".mat ", "/shared-mounts/sinhas/tabebor2/pwmscan/pwm_RBP_human/",
        names(RBP_human_list)[i], ".lpm\n"),
      sep = "", append = T, file = "~/Documents/Shayan/BioInf/PWMs/lpm_to_mat_human_RBP.sh" )
}
#write job file for scanning mm9 genome ucsc style using motifScanner

for(i in 1:length(RBP_human_list)){
  cat(c("/shared-mounts/sinhas/tabebor2/pwmscan/pwm_mscan_wrapper ", 
        "-m /shared-mounts/sinhas/tabebor2/pwmscan/pwm_RBP_human/", names(RBP_human_list)[i], 
        ".mat -e 0.0001 -d /shared-mounts/sinhas/tabebor2/pwmscan/genomedb -s hg19 > /shared-mounts/sinhas/tabebor2/pwmscan/bed_RBP_human/", names(RBP_human_list)[i], "____hg19_hits.bed\n"),
      sep = "", append = T, file = "~/Documents/Shayan/BioInf/PWMs/scan_hg19_RBP.job" )
}
######################################################################################################
######################################################################################################
# write scan job for mouse TF
# aa_seq_list <- list.files("~/Documents/Shayan/BioInf/lncRNA/MM9_tiles", full.names = F)
# aa_pre <- "/shared-mounts/sinhas/tabebor2/lncRNA/Seq/"
# 
# for(i in 1:length(aa_seq_list)){
#   cat(c("Rscript --vanilla /shared-mounts/sinhas/tabebor2/lncRNA/motif_scan_run.R",
#         paste0(aa_pre, aa_seq_list[i]),
#         "DNA", 
#         "/shared-mounts/sinhas/tabebor2/lncRNA/PWMs/my_pwm_list_mouse_TF.RData",
#         "/shared-mounts/sinhas/tabebor2/lncRNA/TFscore_by_tile_mm9\n"),
#       file = "~/Documents/Shayan/BioInf/lncRNA/PWM_scan_TF_mm9.job", 
#       append = !(i == 1), 
#       sep = " ")
# }
# 
# # write scan job for mouse TF --> single
# aa_seq_list <- names(MM9_1kb_tiled_owner_3Mfilter)
# aax2 <- unlist(lapply(strsplit(aa_seq_list, "_"), "[[", 1))
# 
# aa_pre <- "Seq_by_Tile/"
# 
# for(i in 1:length(aa_seq_list)){
#   cat(c("Rscript --vanilla motif_scan_run_single.R",
#         paste0(aa_pre,  aa_seq_list[i],".fa"),
#         "DNA", 
#         "PWMs/my_pwm_list_mouse_TF.RData",
#         "TFscore_by_tile_mm9\n"),
#       file = "~/Documents/Shayan/BioInf/lncRNA/PWM_scan_TF_mm9_single.job", 
#       append = !(i == 1), 
#       sep = " ")
# }
# ######################################################################################################
# # write scan job for mouse RBP
# aa_seq_list <- list.files("~/Documents/Shayan/BioInf/lncRNA/MM9_tiles", full.names = F)
# aa_pre <- "/shared-mounts/sinhas/tabebor2/lncRNA/Seq/"
# 
# for(i in 1:length(aa_seq_list)){
#   cat(c("Rscript --vanilla /shared-mounts/sinhas/tabebor2/lncRNA/motif_scan_run.R",
#         paste0(aa_pre, aa_seq_list[i]),
#         "RNA", 
#         "/shared-mounts/sinhas/tabebor2/lncRNA/PWMs/my_pwm_list_mouse_RBP.RData",
#         "/shared-mounts/sinhas/tabebor2/lncRNA/RBPscore_by_tile_mm9\n"),
#       file = "~/Documents/Shayan/BioInf/lncRNA/PWM_scan_RBP_mm9.job", 
#       append = !(i == 1), 
#       sep = " ")
# }
# 
# 
# # write scan job for mouse RBP --> single
# aa_seq_list <- names(MM9_1kb_tiled_owner_3Mfilter)
# aax2 <- unlist(lapply(strsplit(aa_seq_list, "_"), "[[", 1))
# 
# aa_pre <- "Seq_by_Tile/"
# 
# for(i in 1:length(aa_seq_list)){
#   cat(c("Rscript --vanilla motif_scan_run_single.R",
#         paste0(aa_pre, aa_seq_list[i],".fa"),
#         "RNA", 
#         "PWMs/my_pwm_list_mouse_RBP.RData",
#         "RBPscore_by_tile_mm9\n"),
#       file = "~/Documents/Shayan/BioInf/lncRNA/PWM_scan_RBP_mm9_single.job", 
#       append = !(i == 1), 
#       sep = " ")
# }
###################################################################################################################################################
# read the name of bed files, and write job files to merge the overlapping
aa_human_TF <- read.table("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/human_bed_files_TF.txt", )
aa_human_RBP <- read.table("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/human_bed_files_RBP.txt")
aa_mouse_TF <- read.table("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/mouse_bed_files_TF.txt")
aa_mouse_RBP <- read.table("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/mouse_bed_files_RBP.txt")

aahum <- unlist(lapply(strsplit(aa_human_TF$V1, "____"), "[[", 1))
aahum <- unlist(lapply(strsplit(aahum, "__"), "[[", 1))
aahum_unq <- unique(aahum)
cat(c("#!/bin/bash\n"),
    file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_human_TF.sh", 
    append = F)
for (i in 1:length(aahum_unq)){
  aa_wch <- which(aahum %in% aahum_unq[i])
  cat(c("cat", aa_human_TF$V1[aa_wch], " > ", paste0(aahum_unq[i], "_hg19.bed\n")),
      sep = " ",
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_human_TF.sh", 
      append = T)
  cat(c("sort -k1,1 -k2,2n -o", paste0(aahum_unq[i], "_hg19.bed"), paste0(aahum_unq[i], "_hg19.bed\n")), 
      sep = " ",
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_human_TF.sh", 
      append = T)
  
  cat(c("/shared-mounts/sinhas/tabebor2/bedtools merge -d 1 -c 7,8 -o collapse,collapse -delim '|' -i", paste0(aahum_unq[i], "_hg19.bed"), " > ", paste0(aahum_unq[i], "_hg19_TF_merged.bed\n") ), 
      sep = " ",
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_human_TF.sh", 
      append = T)
  cat(c("rm", paste0(aahum_unq[i], "_hg19.bed\n") ), 
      sep = " ",
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_human_TF.sh", 
      append = T)
}

aahum <- unlist(lapply(strsplit(aa_human_RBP$V1, "____"), "[[", 1))
aahum <- unlist(lapply(strsplit(aahum, "__"), "[[", 1))
aahum_unq <- unique(aahum)
cat(c("#!/bin/bash\n"),
    file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_human_RBP.sh", 
    append = F)
for (i in 1:length(aahum_unq)){
  aa_wch <- which(aahum %in% aahum_unq[i])
  cat(c("cat", aa_human_RBP$V1[aa_wch], " > ", paste0(aahum_unq[i], "_hg19.bed\n")),
      sep = " ",
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_human_RBP.sh", 
      append = T)
  cat(c("sort -k1,1 -k2,2n -o", paste0(aahum_unq[i], "_hg19.bed"), paste0(aahum_unq[i], "_hg19.bed\n")), 
      sep = " ",
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_human_RBP.sh", 
      append = T)
  
  cat(c("/shared-mounts/sinhas/tabebor2/bedtools merge -d 1 -c 7,8 -o collapse,collapse -delim '|' -i", paste0(aahum_unq[i], "_hg19.bed"), " > ", paste0(aahum_unq[i], "_hg19_RBP_merged.bed\n") ), 
      sep = " ",
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_human_RBP.sh", 
      append = T)
  cat(c("rm", paste0(aahum_unq[i], "_hg19.bed\n") ), 
      sep = " ",
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_human_RBP.sh", 
      append = T)
}


aahum <- unlist(lapply(strsplit(aa_mouse_TF$V1, "____"), "[[", 1))
aahum <- unlist(lapply(strsplit(aahum, "__"), "[[", 1))
aahum_unq <- unique(aahum)
cat(c("#!/bin/bash\n"),
    file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_mouse_TF.sh", 
    append = F)
for (i in 1:length(aahum_unq)){
  aa_wch <- which(aahum %in% aahum_unq[i])
  cat(c("cat", aa_mouse_TF$V1[aa_wch], " > ", paste0(aahum_unq[i], "_mm9.bed\n")),
      sep = " ",
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_mouse_TF.sh", 
      append = T)
  cat(c("sort -k1,1 -k2,2n -o", paste0(aahum_unq[i], "_mm9.bed"), paste0(aahum_unq[i], "_mm9.bed\n")), 
      sep = " ",
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_mouse_TF.sh", 
      append = T)
  
  cat(c("/shared-mounts/sinhas/tabebor2/bedtools merge -d 1 -c 7,8 -o collapse,collapse -delim '|' -i", paste0(aahum_unq[i], "_mm9.bed"), " > ", paste0(aahum_unq[i], "_mm9_TF_merged.bed\n") ), 
      sep = " ",
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_mouse_TF.sh", 
      append = T)
  cat(c("rm", paste0(aahum_unq[i], "_mm9.bed\n") ), 
      sep = " ",
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_mouse_TF.sh", 
      append = T)
}


aahum <- unlist(lapply(strsplit(aa_mouse_RBP$V1, "____"), "[[", 1))
aahum <- unlist(lapply(strsplit(aahum, "__"), "[[", 1))
aahum_unq <- unique(aahum)
cat(c("#!/bin/bash\n"),
    file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_mouse_RBP.sh", 
    append = F)
for (i in 1:length(aahum_unq)){
  aa_wch <- which(aahum %in% aahum_unq[i])
  cat(c("cat", aa_mouse_RBP$V1[aa_wch], " > ", paste0(aahum_unq[i], "_mm9.bed\n")),
      sep = " ",
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_mouse_RBP.sh",
      append = T)
  cat(c("sort -k1,1 -k2,2n -o", paste0(aahum_unq[i], "_mm9.bed"), paste0(aahum_unq[i], "_mm9.bed\n")),
      sep = " ",
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_mouse_RBP.sh",
      append = T)

  cat(c("/shared-mounts/sinhas/tabebor2/bedtools merge -d 1 -c 7,8 -o collapse,collapse -delim '|' -i", paste0(aahum_unq[i], "_mm9.bed"), " > ", paste0(aahum_unq[i], "_mm9_TF_merged.bed\n") ),
      sep = " ",
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_mouse_RBP.sh",
      append = T)
  cat(c("rm", paste0(aahum_unq[i], "_mm9.bed\n") ), 
      sep = " ",
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_mouse_RBP.sh", 
      append = T)
}


# write one executable file for each TF and submit to hal
aahum <- unlist(lapply(strsplit(aa_mouse_RBP$V1, "____"), "[[", 1))
aahum <- unlist(lapply(strsplit(aahum, "__"), "[[", 1))
aahum_unq <- unique(aahum)
# cat(c("#!/bin/bash\n"),
#     file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/concat_mouse_RBP.sh", 
#     append = F)
for (i in 63:length(aahum_unq)){
  aa_wch <- which(aahum %in% aahum_unq[i])
  cat(c("#!/bin/bash\n"),
      file = paste0("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/single_jobs_RBP_mouse/concat_mouse_RBP_", aahum_unq[i], ".sh"), 
      append = F)
  cat(c("cat", aa_mouse_RBP$V1[aa_wch], " > ", paste0(aahum_unq[i], "_mm9.bed\n")),
      sep = " ",
      file =  paste0("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/single_jobs_RBP_mouse/concat_mouse_RBP_", aahum_unq[i], ".sh"),
      append = T)
  cat(c("sort -k1,1 -k2,2n -o", paste0(aahum_unq[i], "_mm9.bed"), paste0(aahum_unq[i], "_mm9.bed\n")),
      sep = " ",
      file =  paste0("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/single_jobs_RBP_mouse/concat_mouse_RBP_", aahum_unq[i], ".sh"),
      append = T)
  
  cat(c("/shared-mounts/sinhas/tabebor2/bedtools merge -d 1 -c 7,8 -o collapse,collapse -delim '|' -i", paste0(aahum_unq[i], "_mm9.bed"), " > ", paste0(aahum_unq[i], "_mm9_TF_merged.bed\n") ),
      sep = " ",
      file =  paste0("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/single_jobs_RBP_mouse/concat_mouse_RBP_", aahum_unq[i], ".sh"),
      append = T)
  cat(c("rm", paste0(aahum_unq[i], "_mm9.bed\n") ), 
      sep = " ",
      file =  paste0("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/single_jobs_RBP_mouse/concat_mouse_RBP_", aahum_unq[i], ".sh"), 
      append = T)
}
for(i in 63:length(aahum_unq)){
  cat(c(paste0("/shared-mounts/sinhas/tabebor2/pwmscan/bed_RBP_mouse/concat_mouse_RBP_", aahum_unq[i], ".sh\n")), 
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/run_concat_mouse_RBP.sh", append = !(i==1), sep = "")
}
for(i in 63:length(aahum_unq)){
  cat(c(paste0("chmod +x /shared-mounts/sinhas/tabebor2/pwmscan/bed_RBP_mouse/concat_mouse_RBP_", aahum_unq[i], ".sh\n")), 
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/run_concat_mouse_RBP_exec.sh", append = !(i==1), sep = "")
}

aahum <- unlist(lapply(strsplit(aa_mouse_TF$V1, "____"), "[[", 1))
aahum <- unlist(lapply(strsplit(aahum, "__"), "[[", 1))
aahum_unq <- unique(aahum)

for (i in 1:length(aahum_unq)){
  cat(c("#!/bin/bash\n"),
      file = paste0("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/single_jobs_TF_mouse/concat_mouse_TF_", aahum_unq[i], ".sh"), 
      append = F)
  aa_wch <- which(aahum %in% aahum_unq[i])
  cat(c("cat", aa_mouse_TF$V1[aa_wch], " > ", paste0(aahum_unq[i], "_mm9.bed\n")),
      sep = " ",
      file = paste0("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/single_jobs_TF_mouse/concat_mouse_TF_", aahum_unq[i], ".sh"),
      append = T)
  cat(c("sort -k1,1 -k2,2n -o", paste0(aahum_unq[i], "_mm9.bed"), paste0(aahum_unq[i], "_mm9.bed\n")), 
      sep = " ",
      file = paste0("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/single_jobs_TF_mouse/concat_mouse_TF_", aahum_unq[i], ".sh"),
      append = T)
  
  cat(c("/shared-mounts/sinhas/tabebor2/bedtools merge -d 1 -c 7,8 -o collapse,collapse -delim '|' -i", paste0(aahum_unq[i], "_mm9.bed"), " > ", paste0(aahum_unq[i], "_mm9_TF_merged.bed\n") ), 
      sep = " ",
      file = paste0("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/single_jobs_TF_mouse/concat_mouse_TF_", aahum_unq[i], ".sh"), 
      append = T)
  cat(c("rm", paste0(aahum_unq[i], "_mm9.bed\n") ), 
      sep = " ",
      file = paste0("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/single_jobs_TF_mouse/concat_mouse_TF_", aahum_unq[i], ".sh"), 
      append = T)
}
for(i in 1:length(aahum_unq)){
  cat(c(paste0("/shared-mounts/sinhas/tabebor2/pwmscan/bed_TF_mouse/concat_mouse_TF_", aahum_unq[i], ".sh\n")), 
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/run_concat_mouse_TF.sh", append = !(i==1), sep = "")
}
for(i in 1:length(aahum_unq)){
  cat(c(paste0("chmod +x /shared-mounts/sinhas/tabebor2/pwmscan/bed_TF_mouse/concat_mouse_TF_", aahum_unq[i], ".sh\n")), 
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/run_concat_mouse_TF_exec.sh", append = !(i==1), sep = "")
}


aahum <- unlist(lapply(strsplit(aa_human_RBP$V1, "____"), "[[", 1))
aahum <- unlist(lapply(strsplit(aahum, "__"), "[[", 1))
aahum_unq <- unique(aahum)

for (i in 74:length(aahum_unq)){
  aa_wch <- which(aahum %in% aahum_unq[i])
  cat(c("#!/bin/bash\n"),
      file = paste0("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/single_jobs_RBP_human/concat_human_RBP_", aahum_unq[i], ".sh"),
      append = F)
  cat(c("cat", aa_human_RBP$V1[aa_wch], " > ", paste0(aahum_unq[i], "_hg19.bed\n")),
      sep = " ",
      file = paste0("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/single_jobs_RBP_human/concat_human_RBP_", aahum_unq[i], ".sh"),
      append = T)
  cat(c("sort -k1,1 -k2,2n -o", paste0(aahum_unq[i], "_hg19.bed"), paste0(aahum_unq[i], "_hg19.bed\n")), 
      sep = " ",
      file = paste0("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/single_jobs_RBP_human/concat_human_RBP_", aahum_unq[i], ".sh"),
      append = T)
  
  cat(c("/shared-mounts/sinhas/tabebor2/bedtools merge -d 1 -c 7,8 -o collapse,collapse -delim '|' -i", paste0(aahum_unq[i], "_hg19.bed"), " > ", paste0(aahum_unq[i], "_hg19_RBP_merged.bed\n") ), 
      sep = " ",
      file = paste0("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/single_jobs_RBP_human/concat_human_RBP_", aahum_unq[i], ".sh"),
      append = T)
  cat(c("rm", paste0(aahum_unq[i], "_hg19.bed\n") ), 
      sep = " ",
      file = paste0("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/single_jobs_RBP_human/concat_human_RBP_", aahum_unq[i], ".sh"),
      append = T)
}
for(i in 74:length(aahum_unq)){
  cat(c(paste0("/shared-mounts/sinhas/tabebor2/pwmscan/bed_RBP_human/concat_human_RBP_", aahum_unq[i], ".sh\n")), 
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/run_concat_human_RBP.sh", append = !(i==1), sep = "")
}

for(i in 74:length(aahum_unq)){
  cat(c(paste0("chmod +x /shared-mounts/sinhas/tabebor2/pwmscan/bed_RBP_human/concat_human_RBP_", aahum_unq[i], ".sh\n")), 
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/run_concat_human_RBP_exec.sh", append = !(i==1), sep = "")
}

aahum <- unlist(lapply(strsplit(aa_human_TF$V1, "____"), "[[", 1))
aahum <- unlist(lapply(strsplit(aahum, "__"), "[[", 1))
aahum_unq <- unique(aahum)

for (i in 1:length(aahum_unq)){
  cat(c("#!/bin/bash\n"),
      file = paste0("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/single_jobs_TF_human/concat_human_TF_", aahum_unq[i], ".sh"), 
      append = F)
  aa_wch <- which(aahum %in% aahum_unq[i])
  cat(c("cat", aa_human_TF$V1[aa_wch], " > ", paste0(aahum_unq[i], "_hg19.bed\n")),
      sep = " ",
      file = paste0("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/single_jobs_TF_human/concat_human_TF_", aahum_unq[i], ".sh"),
      append = T)
  cat(c("sort -k1,1 -k2,2n -o", paste0(aahum_unq[i], "_hg19.bed"), paste0(aahum_unq[i], "_hg19.bed\n")), 
      sep = " ",
      file = paste0("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/single_jobs_TF_human/concat_human_TF_", aahum_unq[i], ".sh"),
      append = T)
  
  cat(c("/shared-mounts/sinhas/tabebor2/bedtools merge -d 1 -c 7,8 -o collapse,collapse -delim '|' -i", paste0(aahum_unq[i], "_hg19.bed"), " > ", paste0(aahum_unq[i], "_hg19_TF_merged.bed\n") ), 
      sep = " ",
      file = paste0("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/single_jobs_TF_human/concat_human_TF_", aahum_unq[i], ".sh"), 
      append = T)
  cat(c("rm", paste0(aahum_unq[i], "_hg19.bed\n") ), 
      sep = " ",
      file = paste0("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/single_jobs_TF_human/concat_human_TF_", aahum_unq[i], ".sh"), 
      append = T)
}
for(i in 1:length(aahum_unq)){
  cat(c(paste0("/shared-mounts/sinhas/tabebor2/pwmscan/bed_TF_human/concat_human_TF_", aahum_unq[i], ".sh\n")), 
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/run_concat_human_TF.sh", append = !(i==1), sep = "")
}

for(i in 1:length(aahum_unq)){
  cat(c(paste0("chmod +x /shared-mounts/sinhas/tabebor2/pwmscan/bed_TF_human/concat_human_TF_", aahum_unq[i], ".sh\n")), 
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/run_concat_human_TF_exec.sh", append = !(i==1), sep = "")
}
#/shared-mounts/sinhas/tabebor2/pwmscan/bin/filterOverlaps -l10
#sort -s -c -k1,1 -k2,2n -k6,6 <BED file>.

# write jobs to get the number of sites for each TF on each p6 tile

aahum <- unlist(lapply(strsplit(aa_mouse_TF$V1, "____"), "[[", 1))
aahum <- unlist(lapply(strsplit(aahum, "__"), "[[", 1))
aahum_unq <- unique(aahum)


for (i in 1:length(aahum_unq)){
  cat(c("Rscript --vanilla read_bed_ovl_partition6.R /shared-mounts/sinhas/tabebor2/lncRNA/partition_6_tile_list.txt ", paste0("/shared-mounts/sinhas/tabebor2/pwmscan/bed_TF_mouse/" , aahum_unq[i], "_mm9_TF_merged.bed ", "/shared-mounts/sinhas/tabebor2/lncRNA/TFpwmFeatures_mm9_partition6\n")),
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/mouse_TF_overlap_p6.job", 
      append = !(i==1), sep = "")
}


aahum <- unlist(lapply(strsplit(aa_mouse_RBP$V1, "____"), "[[", 1))
aahum <- unlist(lapply(strsplit(aahum, "__"), "[[", 1))
aahum_unq <- unique(aahum)


for (i in 1:length(aahum_unq)){
  cat(c("Rscript --vanilla read_bed_ovl_partition6.R /shared-mounts/sinhas/tabebor2/lncRNA/partition_6_tile_list.txt ", paste0("/shared-mounts/sinhas/tabebor2/pwmscan/bed_RBP_mouse/" , aahum_unq[i], "_mm9_TF_merged.bed ", "/shared-mounts/sinhas/tabebor2/lncRNA/RBPpwmFeatures_mm9_partition6\n")),
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/mouse_RBP_overlap_p6.job", 
      append = !(i==1), sep = "")
}


###################################################################################################################################################
# adding repeat element features

Partition_6_repeat_features <- matrix(nrow = nrow(Partition_6_dfs), 
                                      ncol = length(MM9_repeats_list_tiled))
rownames(Partition_6_repeat_features) <- Partition_6_dfs$tile_name
colnames(Partition_6_repeat_features) <- names(MM9_repeats_list_tiled)

for(i in 1:length(MM9_repeats_list_tiled)){
  aamch <- match(Partition_6_dfs$tile_name, MM9_repeats_list_tiled[[i]]$tile_name)
  Partition_6_repeat_features[ , i] <- MM9_repeats_list_tiled[[i]]$nu_hits[aamch]
}
Partition_6_repeat_features[is.na(Partition_6_repeat_features)]  <- 0

save(list= c("Partition_6_repeat_features"), file = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_repeat_features.RData")

###################################################################################################################################################
# get repeat element features for lncRNAs

aalnc <- read.table("~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/mES/MM9_lncRNA_28.bed")
aalncGR <- makeGRangesFromDataFrame(aalnc, keep.extra.columns = T, seqnames.field = "V1", start.field = "V2", end.field = "V3")
mESC_repeat_lncRNA_features <- matrix(nrow = length(aalncGR), ncol = length(MM9_repeats_list_GR))
rownames(mESC_repeat_lncRNA_features) <- aalncGR$V7
colnames(mESC_repeat_lncRNA_features) <- names(MM9_repeats_list_GR)

for(i in 1:length(MM9_repeats_list_GR)){
  print(i)
  aaov <- findOverlaps(aalncGR, MM9_repeats_list_GR[[i]])
  if(length(aaov@from ) > 0){
    aadf <- data.frame(lncRNA_index = aaov@from,
                       repeat_index = aaov@to, 
                       interaction_index = c(1:length(aaov@to)))
    aadf_agg <- aggregate(aadf[c("interaction_index")],
                          by = aadf[c("lncRNA_index")],
                          FUN = length)
    #aa_nu_int <- unlist(lapply(aadf_agg$interaction_index, length))
    aa_my_ind <- match(aadf_agg$lncRNA_index, c(1:length(aalncGR)))
    
    mESC_repeat_lncRNA_features[aa_my_ind,i] <- aadf_agg$interaction_index
  }

  
}
mESC_repeat_lncRNA_features[is.na(mESC_repeat_lncRNA_features)] <- 0
save(list = c("mESC_repeat_lncRNA_features"), file = "~/Documents/Shayan/BioInf/lncRNA/mESC_repeat_lncRNA_features.RData")

# get repeat element features for TF scanned sites
aa_mouse_TF <- read.table("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/mouse_bed_files_TF.txt")
aa_mouse_RBP <- read.table("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/mouse_bed_files_RBP.txt")

aahum <- unlist(lapply(strsplit(aa_mouse_TF$V1, "____"), "[[", 1))
aahum <- unlist(lapply(strsplit(aahum, "__"), "[[", 1))
aahum_unq <- unique(aahum)


for (i in 1:length(aahum_unq)){
  cat(c("Rscript --vanilla read_bed_ovl_lncRNA.R /shared-mounts/sinhas/tabebor2/lncRNA/MM9_lncRNA_28.bed ", paste0("/shared-mounts/sinhas/tabebor2/pwmscan/bed_TF_mouse/" , aahum_unq[i], "_mm9_TF_merged.bed ", "/shared-mounts/sinhas/tabebor2/lncRNA/TFpwmFeatures_mm9_lncRNA\n")),
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/mouse_TF_overlap_lncRNA.job", 
      append = !(i==1), sep = "")
}

# get repeat element features for RBP scanned sites


aahum <- unlist(lapply(strsplit(aa_mouse_RBP$V1, "____"), "[[", 1))
aahum <- unlist(lapply(strsplit(aahum, "__"), "[[", 1))
aahum_unq <- unique(aahum)


for (i in 1:length(aahum_unq)){
  cat(c("Rscript --vanilla read_bed_ovl_lncRNA.R /shared-mounts/sinhas/tabebor2/lncRNA/MM9_lncRNA_28.bed ", paste0("/shared-mounts/sinhas/tabebor2/pwmscan/bed_RBP_mouse/" , aahum_unq[i], "_mm9_TF_merged.bed ", "/shared-mounts/sinhas/tabebor2/lncRNA/RBPpwmFeatures_mm9_lncRNA\n")),
      file = "~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/mouse_RBP_overlap_lncRNA.job", 
      append = !(i==1), sep = "")
}

#########################################
# Read GRO-Seq and PROseq data on mESC
library(rtracklayer)
#mESC_GROseq_list <- list()
"~/Documents/Shayan/BioInf/lncRNA/Partition_6_chunk_random.RData"
Partition_6_dfs_GR <- MM9_1kb_tiled_GR[match(Partition_6_dfs$tile_name, MM9_1kb_tiled$tile)]

#mESC_GROseq_list[[1]] 
aatst <- import(format = "BigWig",
                con = "~/Documents/Shayan/BioInf/lncRNA/GROSeq/2/GSE73441_GRO-SEQ-mESC-shControl-NO_DOX.plus.bw",
                which = Partition_6_dfs_GR)
aatst2 <- import(format = "BigWig",
                con = "~/Documents/Shayan/BioInf/lncRNA/GROSeq/2/GSE73441_GRO-SEQ-mESC-shControl-NO_DOX.minus.bw",
                which = Partition_6_dfs_GR)

aatstpro <- import(format = "BigWig",
                con = "~/Documents/Shayan/BioInf/lncRNA/GROSeq/3/GSE130691_mm9_mESC_PROseq_WT_N3_dedup_F.bw",
                which = Partition_6_dfs_GR)
aatst2pro <- import(format = "BigWig",
                 con = "~/Documents/Shayan/BioInf/lncRNA/GROSeq/3/GSE130691_mm9_mESC_PROseq_WT_N3_dedup_R.bw",
                 which = Partition_6_dfs_GR)

aatstpro_all <- c(aatstpro, aatst2pro)
aatst2$score <- -1 * aatst2$score
aatst2_all <- c(aatst, aatst2)


aaov1 <- findOverlaps(Partition_6_dfs_GR, aatstpro_all)
aadf1 <- data.frame(tile_index = aaov1@from,
                   RL_index = aaov1@to, 
                   score = aatstpro_all$score[aaov1@to])

aadf_agg1 <- aggregate(aadf1[c("score")],
                        by = aadf1[c("tile_index")],
                        FUN = sum)

aaov2 <- findOverlaps(Partition_6_dfs_GR, aatst2_all)
aadf2 <- data.frame(tile_index = aaov2@from,
                    RL_index = aaov2@to, 
                    score = aatst2_all$score[aaov2@to])

aadf_agg2 <- aggregate(aadf2[c("score")],
                       by = aadf2[c("tile_index")],
                       FUN = sum)

Partition_6_GROseq <- numeric(length(Partition_6_dfs_GR))
Partition_6_PROseq <- numeric(length(Partition_6_dfs_GR))
Partition_6_GROseq[aadf_agg2$tile_index] <- aadf_agg2$score
Partition_6_PROseq[aadf_agg1$tile_index] <- aadf_agg1$score


summary(Partition_6_PROseq)
summary(Partition_6_GROseq)

save(list = c("Partition_6_PROseq", "Partition_6_GROseq"), 
     file = "~/Documents/Shayan/BioInf/lncRNA/partition_6_GROseqPROseq_features.RData")

load("~/Documents/Shayan/BioInf/lncRNA/partition_6_expression_features.RData")
partition_6_expression_features$GROseq <- Partition_6_GROseq
partition_6_expression_features$PROseq <- Partition_6_PROseq
save(list = c("partition_6_expression_features"), 
     file = "~/Documents/Shayan/BioInf/lncRNA/partition_6_expression_features.RData")


######################################################################
# add methylation feature

aatst <- import(format = "bedGRaph",
                con = "~/Documents/Shayan/BioInf/lncRNA/Methylation_mESC/GSM3752614_WGBS_WT_Rep1_mm9.bedgraph"
                ,which = Partition_6_dfs_GR
                )
aatst2 <- import(format = "bedGRaph",
                 con = "~/Documents/Shayan/BioInf/lncRNA/Methylation_mESC/GSM4558210_WGBS_WT_Rep2_mm9.bedgraph"
                 ,which = Partition_6_dfs_GR
                 )

aaov1 <- findOverlaps(Partition_6_dfs_GR, aatst)
aadf1 <- data.frame(tile_index = aaov1@from,
                    RL_index = aaov1@to, 
                    score = aatst$score[aaov1@to])

aadf_agg1 <- aggregate(aadf1[c("score")],
                       by = aadf1[c("tile_index")],
                       FUN = sum)

aaov2 <- findOverlaps(Partition_6_dfs_GR, aatst2)
aadf2 <- data.frame(tile_index = aaov2@from,
                    RL_index = aaov2@to, 
                    score = aatst2$score[aaov2@to])

aadf_agg2 <- aggregate(aadf2[c("score")],
                       by = aadf2[c("tile_index")],
                       FUN = sum)

aaPartition_6_me1 <- numeric(length(Partition_6_dfs_GR))
aaPartition_6_me2 <- numeric(length(Partition_6_dfs_GR))
aaPartition_6_me1[aadf_agg1$tile_index] <- aadf_agg1$score
aaPartition_6_me2[aadf_agg2$tile_index] <- aadf_agg2$score
Partition_6_methylation <- (aaPartition_6_me1 + aaPartition_6_me2)/2
summary(Partition_6_methylation)
names(Partition_6_methylation) <- Partition_6_dfs_GR$tile
cor(aaPartition_6_me1, aaPartition_6_me2)
hist(Partition_6_methylation)
hist(aaPartition_6_me1)
hist(aaPartition_6_me2)
save(list = c("Partition_6_methylation"), 
     file = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_methylation.RData")

######################################################################
# create network-guided pair features for RBP and TFs
aa_ppi <- read.table("~/Documents/Shayan/BioInf/lncRNA/Protein_protein_int/ppi.txt", 
                     header = T,
                     stringsAsFactors = F)
aa_ppi <- aa_ppi[aa_ppi$cellType %in% c("MESC", "MESC-J1"),]
aa_all_ppi <- tolower(unique(union(aa_ppi$protein1Name, aa_ppi$protein2Name)))
aa_all_ppi[which(aa_all_ppi %in% tolower(colnames(Partition_6_TF_scanned)))]
aa_all_ppi[which(aa_all_ppi %in% tolower(colnames(Partition_6_RBP_scanned)))]
aa_all_ppi_filt <- union(aa_all_ppi[which(aa_all_ppi %in% tolower(colnames(Partition_6_TF_scanned)))],
                         aa_all_ppi[which(aa_all_ppi %in% tolower(colnames(Partition_6_RBP_scanned)))])

aa_ppi_filt <- aa_ppi[((tolower(aa_ppi$protein1Name) %in% aa_all_ppi_filt) & (tolower(aa_ppi$protein2Name) %in% aa_all_ppi_filt)),]
aa_ppi_filt <- aa_ppi_filt[! duplicated(aa_ppi_filt[,c("protein1Name", "protein2Name")]),]
aa_ppi_filt2 <- aa_ppi_filt[, c("protein1Name", "protein2Name")]
write.table(aa_ppi_filt2, file = "~/Documents/Shayan/BioInf/lncRNA/Protein_protein_int/PPI_relevant_mesc.txt", quote = F, row.names = F, col.names = F, append = F)
PPI_relevant_mesc <- aa_ppi_filt2
##############################
all_prot_predicted_site <- cbind(Partition_6_TF_scanned, Partition_6_RBP_scanned)
PPI_relevant_mesc <- read.table("~/Documents/Shayan/BioInf/lncRNA/Protein_protein_int/PPI_relevant_mesc.txt", header = F, stringsAsFactors = F)

aa_lncRNA_RBP <- read.table("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/RBP_pwm_lncRNA_all.txt", header = T)
aa_lncRNA_TF <- read.table("~/Documents/Shayan/BioInf/lncRNA/Motif_scan_cisbp/TF_pwm_lncRNA_all.txt", header = T)
aa_lnc_all <- cbind(aa_lncRNA_TF, aa_lncRNA_RBP)
aaw <- which(! (colnames(all_prot_predicted_site) %in%  colnames(aa_lnc_all)))
if(length(aaw) > 0){
  all_prot_predicted_site <- all_prot_predicted_site[, -aaw]
}

aa_lnc_all <- aa_lnc_all[,match(colnames(all_prot_predicted_site), colnames(aa_lnc_all))]

aam1 <- match(tolower(PPI_relevant_mesc[, 1]), tolower(colnames(all_prot_predicted_site)))
aam2 <- match(tolower(PPI_relevant_mesc[, 2]), tolower(colnames(all_prot_predicted_site)))
aamy_net1 <- as.matrix(cbind(aam1, aam2))
#aamy_net2 <- as.matrix(cbind(c(1:ncol(all_prot_predicted_site)), c(1:ncol(all_prot_predicted_site))))
aa_netlist <- list(heteroDimer = aamy_net1
                   # , homoDimer = aamy_net2
)


#Partition_6_network_pairs <- add_pair_features_network(lncRNA_feat = , tile_feat, my_network_list = numeric(0), presence_thr=0)
aa_lncNames <- read.table("~/Documents/Shayan/BioInf/lncRNA/TAD_annotations/mES/MM9_lncRNA_28.bed", stringsAsFactors = F)$V7


load("Partition_6_chunk_random.RData")
Partition_6_network_pairs <- numeric(length = nrow(Partition_6_dfs))
source("~/Documents/Shayan/BioInf/lncRNA/pair_feat_cons.R")
#
for(i in 1:nrow(aa_lnc_all)){
  print(i)
  aa_tmp_tile <- which(Partition_6_dfs$owner %in% aa_lncNames[i])
  aa_tmp <- add_pair_features_network(lncRNA_feat = aa_lnc_all[i, ],
                              tile_feat = all_prot_predicted_site[aa_tmp_tile,],
                              my_network_list = aa_netlist,
                              presence_thr=0)
  Partition_6_network_pairs[aa_tmp_tile] <- aa_tmp
}
save(list = c("Partition_6_network_pairs"), file = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_network_pairs.RData")

############################################################
# get dinucleotide freq matrix
load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_dfs_GR.RData")
aa_p6_seq <- getSeq(x = BSgenome.Mmusculus.UCSC.mm9, names = Partition_6_dfs_GR)
aa_p6_dinuc <- dinucleotideFrequency(aa_p6_seq,)
Partition_6_dinuc_freq <- aa_p6_dinuc
save(list = c("Partition_6_dinuc_freq"), file = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_dinuc_freq.RData")


############################################################
#Get CG skewnesss
#get AT skewness
#
Partition_6_nuc_freq <- oligonucleotideFrequency(aa_p6_seq, width = 1)
Partition_6_nuc_GC_skew <-( Partition_6_nuc_freq[,3]-Partition_6_nuc_freq[,2])/( Partition_6_nuc_freq[,3]+Partition_6_nuc_freq[,2] + 1)
Partition_6_nuc_TA_skew <-( Partition_6_nuc_freq[,4]-Partition_6_nuc_freq[,1])/( Partition_6_nuc_freq[,4]+Partition_6_nuc_freq[,1] + 1)

save(list = c("Partition_6_nuc_freq", "Partition_6_nuc_GC_skew", "Partition_6_nuc_TA_skew"), file = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_nuc_freq.RData")



############################################################
# sum over reverse compelments
aacol <- DNAStringSet(colnames(Tile_kmer_features_partition6))
aacol_rev <- reverseComplement(aacol)
tmp_Tile_kmer_features_partition6 <- Tile_kmer_features_partition6
aa_kmer_revComp <- matrix(nrow = nrow(tmp_Tile_kmer_features_partition6), ncol = 0)
aaname <- character(0)
for(i in 1:length(aacol)){
  print(i)
  if(ncol(tmp_Tile_kmer_features_partition6) > 0){
    aa_tmp <- unique(as.character(c(aacol[i], aacol_rev[i])))
    if(all(aa_tmp %in% colnames(tmp_Tile_kmer_features_partition6))){
      aaname <- c(aaname, paste(as.character(c(aacol[i], aacol_rev[i])), collapse = "__"))
      if(length(aa_tmp) == 2){
        aa_kmer_revComp <- cbind(aa_kmer_revComp, (tmp_Tile_kmer_features_partition6[, aa_tmp[1]] + tmp_Tile_kmer_features_partition6[,aa_tmp[2]]))
      }else if(length(aa_tmp) == 1){
        aa_kmer_revComp <- cbind(aa_kmer_revComp, tmp_Tile_kmer_features_partition6[, aa_tmp[1]] )
      }
      tmp_Tile_kmer_features_partition6 <- tmp_Tile_kmer_features_partition6[, -c(which(colnames(tmp_Tile_kmer_features_partition6) %in% aa_tmp))]
      
    }
    
  }
  
  
}


colnames(aa_kmer_revComp) <- aaname
Tile_kmer_features_partition6_revComp <- aa_kmer_revComp
save(list = c("Tile_kmer_features_partition6_revComp"),
     file = "~/Documents/Shayan/BioInf/lncRNA/Tile_kmer_features_partition6_revComp.RData")


#########################################################################################################################################################
#########################################################################################################################################################
# hg19 to hg38 for motifscnas


aa_fname <- list.files("/shared-mounts/sinhas/tabebor2/pwmscan/bed_TF_human", pattern = "TF_merged.bed$")
aa_fname1 <- unlist(lapply(strsplit(aa_fname, "_"), "[[", 1))
aafiles_full <- paste0("/shared-mounts/sinhas/tabebor2/pwmscan/bed_TF_human/",aa_fname)
aafiles_full_2 <- paste0("/shared-mounts/sinhas/tabebor2/pwmscan/bed_TF_human_hg38/",aa_fname1, "_TFhg38.bed" )
#system("mkdir ~/Documents/Shayan/BioInf/lncRNA/oRNAment/processed_mm9/")
# cat(c("#!/bin/bash\n"), 
#     file = "/shared-mounts/sinhas/tabebor2/pwmscan/liftover_TF_to_hg38.sh",
#     append = F)
for(i in 1:length(aafiles_full)){
  cat(c("/shared-mounts/sinhas/tabebor2/liftOver -bedPlus=3", 
        aafiles_full[i],
        "/shared-mounts/sinhas/tabebor2/liftOver_chain/hg19ToHg38.over.chain", 
        aafiles_full_2[i],
        paste0("/shared-mounts/sinhas/tabebor2/pwmscan/unmapped_TF_hg19_to_hg38/", aa_fname1[i], ".unmapped\n")),
      file = "/shared-mounts/sinhas/tabebor2/pwmscan/liftover_TF_to_hg38.job",
      sep = " ", append = !(i==1))
}

aa_fname <- list.files("/shared-mounts/sinhas/tabebor2/pwmscan/bed_RBP_human", pattern = "RBP_merged.bed$")
aa_fname1 <- unlist(lapply(strsplit(aa_fname, "_"), "[[", 1))
aafiles_full <- paste0("/shared-mounts/sinhas/tabebor2/pwmscan/bed_RBP_human/",aa_fname)
aafiles_full_2 <- paste0("/shared-mounts/sinhas/tabebor2/pwmscan/bed_RBP_human_hg38/",aa_fname1, "_RBPhg38.bed" )
#system("mkdir ~/Documents/Shayan/BioInf/lncRNA/oRNAment/processed_mm9/")
# cat(c("#!/bin/bash\n"), 
#     file = "/shared-mounts/sinhas/tabebor2/pwmscan/liftover_RBP_to_hg38.sh",
#     append = F)
for(i in 1:length(aafiles_full)){
  cat(c("/shared-mounts/sinhas/tabebor2/liftOver -bedPlus=3", 
        aafiles_full[i],
        "/shared-mounts/sinhas/tabebor2/liftOver_chain/hg19ToHg38.over.chain", 
        aafiles_full_2[i],
        paste0("/shared-mounts/sinhas/tabebor2/pwmscan/unmapped_RBP_hg19_to_hg38/", aa_fname1[i], ".unmapped\n")),
      file = "/shared-mounts/sinhas/tabebor2/pwmscan/liftover_RBP_to_hg38.job",
      sep = " ", append = !(i==1))
}
##########################################################################################################
# Adding R-loop feature


# wigToBigWig(x = "~/Documents/Shayan/BioInf/lncRNA/Rloop_mesc/GSM1720620_E14_DRIP.wig", seqinfo = Seqinfo(genome="GRCm38"),
#             dest = paste("~/Documents/Shayan/BioInf/lncRNA/Rloop_mesc/GSM1720620_E14_DRIP", "bw", sep = "."),
#             clip = FALSE)

aatst <- import(format = "wig",
                con = "~/Documents/Shayan/BioInf/lncRNA/Rloop_mesc/GSM1720620_E14_DRIP.wig"
                ,which = Partition_6_dfs_GR
)

aaov1 <- findOverlaps(Partition_6_dfs_GR, aatst)
aadf1 <- data.frame(tile_index = aaov1@from,
                    RL_index = aaov1@to, 
                    score = aatst$score[aaov1@to])

aadf_agg1 <- aggregate(aadf1[c("score")],
                       by = aadf1[c("tile_index")],
                       FUN = sum)
aadf_agg2 <- aggregate(aadf1[c("score")],
                       by = aadf1[c("tile_index")],
                       FUN = max)
aaPartition_6_rloop_sum <- numeric(length(Partition_6_dfs_GR))
aaPartition_6_rloop_sum[aadf_agg1$tile_index] <- aadf_agg1$score
Partition_6_RLoop_sum <- aaPartition_6_rloop_sum

aaPartition_6_rloop_max <- numeric(length(Partition_6_dfs_GR))
aaPartition_6_rloop_max[aadf_agg2$tile_index] <- aadf_agg2$score
Partition_6_RLoop_max <- aaPartition_6_rloop_max

#summary(Partition_6_RLoop)
names(Partition_6_RLoop_max) <- Partition_6_dfs_GR$tile
names(Partition_6_RLoop_sum) <- Partition_6_dfs_GR$tile

#cor(aaPartition_6_rloop, aaPartition_6_me2)
hist(Partition_6_RLoop_max)
save(list = c("Partition_6_RLoop_max","Partition_6_RLoop_sum"), 
     file = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_RLoop.RData")
##########################################################################################################

# add conservation scores for all k-mers in every lncRNA
aalncall <- read.delim("~/Documents/Shayan/BioInf/lncRNA/lncRNA_28_mm9.bed", header = F)
aalncall <- makeGRangesFromDataFrame(aalncall,keep.extra.columns = T,
                                     ignore.strand = T,
                                     seqnames.field = "V1",
                                     start.field = "V2",
                                     end.field = "V3")
aachrn <- as.character(seqnames(aalncall))
aachrn_un <- unique(aachrn)
# convert wig to big wig
cat("#!/bin/bash\n", file = "~/Documents/Shayan/BioInf/mm9_phastcons/convertwig2bigwig.sh", append = F)
for(i in 1:length(aachrn_un)){
  aafile <- paste0("~/Documents/Shayan/BioInf/mm9_phastcons/",aachrn_un[i],".phyloP30way.wigFix.gz" )
  cat(c("~/Documents/Shayan/BioInf/wigToBigWig ", 
        paste0("<(gzcat ", aafile,") "), 
        "~/Documents/Shayan/BioInf/lncRNA/mm9_chr_sizes.txt ",
        paste0("~/Documents/Shayan/BioInf/mm9_phastcons/",aachrn_un[i],".phyloP30way.bw \n")),
      sep = "", append = T,
      file = "~/Documents/Shayan/BioInf/mm9_phastcons/convertwig2bigwig.sh")
}
aaphastlnc <- list()

for (i in 1:length(aalncall)){
  aafile <- paste0("~/Documents/Shayan/BioInf/mm9_phastcons/",aachrn[i],".phyloP30way.bw" )
  aaphastlnc[[i]] <- import(format = "bigwig",
                                     con = aafile,
                                     which = aalncall[i]
  )
}
names(aaphastlnc) <- aalncall$V7
lncRNA_phastCons_score_mm9 <- aaphastlnc
save(list = c("lncRNA_phastCons_score_mm9"), 
     file = "~/Documents/Shayan/BioInf/lncRNA/lncRNA_phastCons_score_mm9.RData")
aa_lncRNA_Seq <- getSeq(x = BSgenome.Mmusculus.UCSC.mm9, names = aalncall)
names(aa_lncRNA_Seq) <- aalncall$V7

#########################3
# count kmers in repeat elements
library(ShortRead)
library(Biostrings)
mouse_repeat_dfam <- readFasta("~/Documents/Shayan/BioInf/lncRNA/Repeat_elements/families.fa")


mouse_repeat_dfam_5mer <- oligonucleotideFrequency(mouse_repeat_dfam@sread, 
                                          width = 5,
                                          step = 1,
                                          as.prob = T,
                                          with.labels = T)
rownames(mouse_repeat_dfam_5mer) <- as.character(mouse_repeat_dfam@id)
save(list= "mouse_repeat_dfam_5mer", file = "~/Documents/Shayan/BioInf/lncRNA/mouse_repeat_dfam_5mer.RData")

load("~/Documents/Shayan/BioInf/lncRNA/lncRNA_kmer_df.RData")
mouse_repeat_dfam_5mer_all <- colMeans(mouse_repeat_dfam_5mer)
lncRNA_kmer_df$repeat_density <- numeric(nrow(lncRNA_kmer_df))
lncRNA_kmer_df$repeat_density <- mouse_repeat_dfam_5mer_all[match(lncRNA_kmer_df$k1, names(mouse_repeat_dfam_5mer_all))] + mouse_repeat_dfam_5mer_all[match(lncRNA_kmer_df$k2, names(mouse_repeat_dfam_5mer_all))]


