library(irlba)
setwd("~/Documents/Shayan/BioInf/lncRNA")

load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/Partition_6_CVfirst_dataset_Neat1.RData")
nrow(my_Dataset)
ncol(my_Dataset)
colnames(my_Dataset)[1:10]
my_name_dic[1:5,]


aac <- c(grep(pattern = "DiFreq", x = my_name_dic[,1]),
         grep(pattern = "5mer", x = my_name_dic[,1]) )
aa_mydkmer <- my_Dataset[, aac]
aa_mydkmer_trcv1 <- aa_mydkmer[my_partition$RCV1 == 1,]
aa_mydkmer_tscv1 <- aa_mydkmer[my_partition$RCV1 == 0,]

# aac <- c(grep(pattern = "ChIP", x = my_name_dic[,1]))
# aa_mydchip <- my_Dataset[, aac]
# aa_mydchip_tccv1 <- aa_mydchip[my_partition$CCV1 == 1,]


aa_mydkmer_rcv1_irlba <- irlba(A = t(aa_mydkmer_trcv1), nv = 30, maxit = 600)
dim(aa_mydkmer_rcv1_irlba$u)
dim(aa_mydkmer_rcv1_irlba$v)
plot(aa_mydkmer_rcv1_irlba$d, type = "l")

sigma_inverse <- 1/aa_mydkmer_rcv1_irlba$d
u_transpose <- as.matrix(t(aa_mydkmer_rcv1_irlba$u))
#aa_neat1t1 <- t(as.matrix(aa_mydkmer_trcv1[1,]))

aa_neat1_pr <- t(sigma_inverse * (u_transpose %*% t(aa_mydkmer_trcv1)))
aa_neat1_prt <- t(sigma_inverse * (u_transpose %*% t(aa_mydkmer_tscv1)))



hist(aa_neat1_pr - aa_mydkmer_rcv1_irlba$v[1,])

aad <-aa_mydkmer_rcv1_irlba$d/sum(aa_mydkmer_rcv1_irlba$d)
plot(cumsum(aad), type = "l")
which(cumsum(aad) >= 0.95)[1]



Sys.setenv('R_MAX_VSIZE'=32000000000)
aatime <- proc.time()
aa_mydkmer_rcv1_irlba <- irlba(A = t(aa_mydkmer_trcv1), nv = 300, maxit = 600)
#aa_mydkmer_rcv1_irlba <- mclapply(t(aa_mydkmer_trcv1), irlba, mc.cores = detectCores(), nv = 10, maxit = 600)
proc.time() - aatime
plot(cumsum(aa_mydkmer_rcv1_irlba$d/sum(aa_mydkmer_rcv1_irlba$d)), type = "l")
plot(aa_mydkmer_rcv1_irlba$d, type = "l")

which(cumsum(aa_mydkmer_rcv1_irlba$d/sum(aa_mydkmer_rcv1_irlba$d)) >= 0.95)[1]


