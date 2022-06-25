args <- commandArgs(trailingOnly = TRUE)
# First argument is the address of the DNA segment
# Second argument is the type "DNA" or "RNA"
# Third argument is the address to motiflist RData file
# Fourth argument is the parent directoy for output
# 
require(ShortRead)
require(PWMEnrich)
source("~/Documents/Shayan/BioInf/lncRNA/tile_motif_LLR_calc_functions.R")
load(args[3])
aa_seq_all <- readFasta(dirPath = args[1])
aa_seq <- as.character(aa_seq_all@sread)
names(aa_seq) <- as.character(aa_seq_all@id)

aafn <- unlist(strsplit(args[1], "\\/"))
aafn <- aafn[length(aafn)]
aafn <- unlist(strsplit(aafn, ".fa"))[1]

MotifScore_multiseq(seq_vector = aa_seq, 
                    seq_type_all = args[2], 
                    motifList_all = my_pwm_list, 
                    output_parent_folder=args[4],
                    output_filename = aafn,
                    bg = c(0.25, 0.25, 0.25, 0.25))