args <- commandArgs(trailingOnly = TRUE)
# args[1] is the address to txt file containing the coordinates of lncRNAs (or basically any feature of interest)
# args[2] is the address to bed file of interest
# args[3] is the address to write the outputs to
library(GenomicRanges)
#load("/shared-mounts/sinhas/tabebor2/lncRNA/MM9_1kb_tiled_GR_filtered.RData")
#load(args[1])
#load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_chunk_random.RData")
#lncRNA_coord <- read.table(args[1] , stringsAsFactors = F, header = F)
aa_name <- unlist(strsplit(args[2], "\\/"))
aa_name <- aa_name[length(aa_name)]
aaname2 <- unlist(strsplit(aa_name, "_mm9"))[1]
print(aaname2)
lncRNA_coord <- read.table(args[1] , stringsAsFactors = F, header = F)
aalncGR <- makeGRangesFromDataFrame(lncRNA_coord, keep.extra.columns = T, seqnames.field = "V1", start.field = "V2", end.field = "V3")

# read bed files

#aa_p6_gr <- MM9_1kb_tiled_GR_filtered[match(Partition_6_dfs$V1, MM9_1kb_tiled_GR_filtered$tile)]
aa_p6_RBP <- matrix(0L, nrow = length(aalncGR), ncol = 1)
#rownames(aa_p6_RBP) <- aalncGR$tile
colnames(aa_p6_RBP) <- aaname2
aa_tmp <- read.table(args[2], stringsAsFactors = F)
aa_tGR <- makeGRangesFromDataFrame(aa_tmp,
                                   seqnames.field = "V1",
                                   start.field = "V2",
                                   end.field = "V3")
print("before_ovl")
aaov <- findOverlaps(aalncGR, aa_tGR)
aadf <- data.frame(tile_index = aaov@from,
                   rbp_binding_index = aaov@to, 
                   interaction_index = c(1:length(aaov@to)))
aadf_agg <- aggregate(aadf[c("interaction_index")],
                      by = aadf[c("tile_index")],
                      FUN = length)
print("after agg")
print(aadf_agg)
if(length(aadf_agg$tile_index) > 0){
  aa_p6_RBP[aadf_agg$tile_index, 1] <- aadf_agg$interaction_index
}

write.table(aa_p6_RBP, file = paste0(args[3],"/", aaname2, "__lncRNAProfile.txt"), append = F, quote = F,row.names = F,  col.names = T)





