Filter_RBP <- function(inp_dataset, 
                       inp_nameDic,
                       feature_index=numeric(0),
                       old_imp=F,randomized=F,
                       partition_Exp_features_address,
                       exp_feature_name = c("RNA_polymerase_II__315.bed", "RNAseq" ,"CAGE", "GROseq", "PROseq"),
                       rand_seed=numeric(0),
                       exp_thr = c(0,0,0,0,0), scale_by=F){
  # /shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/partition_6_expression_features.RData
  # this function takes in dataset and its name_dic, identifies and filters RBP features using the expression data
  #load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_U1snRNA.RData")
  # randomized: logical, if True instead of using the actual expression profiles it uses a shuffled version of the expression profiles.
  # scale_by: if True, instead of setting some features to zero, scales the features by log2(max(GROseq, PROseq) + 1)
  load(partition_Exp_features_address)
  aa_cur_tile <- rownames(inp_dataset)
  aamyma1 <- match(aa_cur_tile, partition_6_expression_features$tile_name)
  aamyma2 <- match(exp_feature_name, colnames(partition_6_expression_features))
  stopifnot(sum(is.na(aamyma1)) == 0,
            sum(is.na(aamyma2)) == 0)
  aa_cur_exp <- partition_6_expression_features[aamyma1, aamyma2]
  if(randomized){
    if(length(rand_seed) == 1){
      set.seed(seed = rand_seed)
    }
    aa_smp <- sample(x = c(1:nrow(aa_cur_exp)), size = nrow(aa_cur_exp), replace = F)
    aa_cur_exp <- aa_cur_exp[aa_smp,]
  }
  if(isTRUE(old_imp)){
    aarbp <- grep(pattern = '_mm9', x =inp_nameDic[, 1])
    aa_rbpPair <- grep(pattern = 'RBP_pair', x = inp_nameDic[, 1])
    aa_u1 <- match("U1snRNA_site", colnames(inp_dataset))
    aa_features <- c(aarbp, aa_rbpPair, aa_u1)
  }else if(length(feature_index) == 0){
    print("length of feature_index has to be greater than one in the new implimentation")
  }else{
    aa_features <- feature_index
  }
  
  if(! scale_by){
    which_filter <- list()
    for(i in 1:ncol(aa_cur_exp)){
      print(i)
      which_filter[[i]] <- which(aa_cur_exp[, i] <= exp_thr[i])
      #print(length(which_filter[[i]]))
    }
    which_filter_int <- Reduce(intersect, which_filter)
    print(length(which_filter_int))
    
    
    if(length(which_filter_int) > 0){
      inp_dataset[which_filter_int, aa_features] <- 0
    }
  }else{
    scale_vec <- log2(apply(aa_cur_exp[, c("GROseq", "PROseq")], 1, max) + 1)
    for(i in 1:length(aa_features)){
      inp_dataset[, aa_features[i]] <- inp_dataset[, aa_features[i]] * scale_vec
    }
  }

  
  
  return(inp_dataset)
}

# load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_datasets_corrected/Partition_6_modified_dataset_2410003L11Rik.RData")
# load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_U1snRNA.RData")
# my_Dataset$U1snRNA_site <-  Partition_6_U1snRNA[match(rownames(my_Dataset), names(Partition_6_U1snRNA))]
# my_Dataset <- my_Dataset[, c(c(1:(ncol(my_Dataset)-2)), (ncol(my_Dataset)), (ncol(my_Dataset) - 1))]
# my_name_dic <- rbind(my_name_dic, c("U1snRNA_site", "U1snRNA_site"))
# 
# aa_tst <- Filter_RBP(inp_dataset = my_Dataset,
#                      inp_nameDic = my_name_dic,
#                      partition_Exp_features_address = "~/Documents/Shayan/BioInf/lncRNA/partition_6_expression_features.RData")


# ransvd = function(A, k=10, p=5) {
#   n = nrow(A)
#   y = A %*% matrix(rnorm(n * (k+p)), nrow=n)
#   q = qr.Q(qr(y))
#   b = t(q) %*% A
#   svd = svd(b)
#   list(u=q %*% svd$u, d=svd$d, v=svd$v)
# }


# load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/Partition_6_CVfirst_dataset_Malat1.RData")
# aax <- as.matrix(my_Dataset[, 1:(ncol(my_Dataset) - 1)])
# 
# ptm <- proc.time()
# aaxsvd <- svd(x = aax, nu = 10, nv = ncol(aax))
# proc.time() - ptm
# aatime <- system.time()
# aatime2 <- system.time(aaxsvd2 <- ransvd(aax))



feature_transformer <- function(dataset,
                                partitioning_vec,
                                name_Dic,
                                feat_family_list,
                                family_transform, 
                                family_trans_len,
                                might_exist_here="/shared-mounts/sinhas/tabebor2/lncRNA/RandomForest/Partition_6_CV2_svd_NOD/precomputed_svd",
                                all_train=F, 
                                lnRNA_name,
                                partition_name){
  # dataset: is the full dataset for a lncRNA --> columns are features, rows are examples
  # partitioning_vec: is a numeric vecror with length equal to nrow(dataset), training examples are 1. test examples are 0.
  # name_Dic: name mapping df
  #feat_family_list: list with one entry for each family of features, entry is a numeric vector containing the index of features (column index) in that family
  #family_transform: logical with one entry for each family of features, if T the family will be transformed
  #family_trans_len: numeric with one entry for each family of features, indicating the number of dominant orthogonal axis to keep
  #all_train: logical, if True performs svd on the full dataset, not just a training subset --> when this is
  # provided then partitioning_vec should be a vector of length equal to nrow(dataset) , where all entries are 1
  #might_exist_here: is the address where a precomputed svd for the family might exist
  # lnRNA_name : name of the lncRNA to look for precomputed svd
  # partition_name : name of the partitioning to look for precomputed svd
  require(irlba)
  stopifnot(length(feat_family_list) == length(names(feat_family_list)),
            length(feat_family_list) == length(family_transform),
            length(family_trans_len) == length(feat_family_list), 
            "label" %in% colnames(dataset), 
            ((all_train & all(partitioning_vec == 1)) | ((!all_train) & (!all(partitioning_vec == 1)))))
  aatrainind <- which(partitioning_vec %in% 1)
  train_orig <- dataset[aatrainind,]
  label_train <- train_orig[, which(colnames(train_orig) == "label")]
  if(! all_train){
    test_orig <- dataset[-aatrainind,]
    label_test <- test_orig[, which(colnames(train_orig) == "label")]
  }
 
  rm(dataset)
  
  
  #print("before var")
  aa_zerovar <- apply(train_orig, MARGIN = 2, FUN = var)
  aawzv <- which(aa_zerovar == 0)
  
  if(length(aawzv) > 0){
    for(i in 1:length(feat_family_list)){
      feat_family_list[[i]] <- setdiff(feat_family_list[[i]] , aawzv)
      if(length(feat_family_list[[i]]) == 0){
        print(paste0("removed all members because of zero variance: ", names(feat_family_list)[i]))
      }
    }
  }
  train_trans_list <- list()
  if(! all_train){
    test_trans_list <- list()
    
  }
  irlba_list <- list()
  for(cur_fam in 1:length(feat_family_list)){
    if(length(feat_family_list[[cur_fam]]) > 0){
      cur_tr <- train_orig[, feat_family_list[[cur_fam]]]
      if(family_transform[cur_fam]){ # should transform
        # check if a precomputed one exists
        aa_pre_fname <- paste0(might_exist_here, "/", names(feat_family_list)[cur_fam], "_", lnRNA_name, "_", partition_name, ".RData")
        if(file.exists(aa_pre_fname)){
          load(aa_pre_fname)
          irlba_list[[cur_fam]] <- my_yoyo
        }else{
          if(family_trans_len[cur_fam] == ncol(cur_tr)){
            irlba_list[[cur_fam]] <- svd(x = t(cur_tr), 
                                         nv = family_trans_len[cur_fam],
                                         nu = family_trans_len[cur_fam])
          }else{
            irlba_list[[cur_fam]] <- irlba(A = t(cur_tr), 
                                           nv = family_trans_len[cur_fam],
                                           maxit = 1000)
          }
          my_yoyo <- irlba_list[[cur_fam]]
          save(list = c("my_yoyo"), file = aa_pre_fname)
        }

        
        rm(cur_tr)
        new_tr <- irlba_list[[cur_fam]]$v
        aacol <- paste0("TR_",names(feat_family_list)[cur_fam],"_",c(1:ncol(new_tr)))
        colnames(new_tr) <- aacol
        if(! all_train){
          sigma_inverse <- 1/irlba_list[[cur_fam]]$d
          u_transpose <- as.matrix(t(irlba_list[[cur_fam]]$u))
          new_ts <- t(sigma_inverse * (u_transpose %*% t(test_orig[, feat_family_list[[cur_fam]]])))
          colnames(new_ts) <- aacol
        }
        

      }else{ # should not transform
        if(length(feat_family_list[[cur_fam]]) == 1){
          new_tr <- matrix(nrow = nrow(train_orig), ncol = length(feat_family_list[[cur_fam]]))
          colnames(new_tr) <- colnames(train_orig)[feat_family_list[[cur_fam]]]
          new_tr[,1:ncol(new_tr)] <- cur_tr
          rm(cur_tr)
          if(! all_train){
            new_ts <- matrix(nrow = nrow(test_orig), ncol = length(feat_family_list[[cur_fam]]))
            colnames(new_ts) <- colnames(test_orig)[feat_family_list[[cur_fam]]]
            new_ts[,1:ncol(new_ts)]<- test_orig[, feat_family_list[[cur_fam]]]
          }

        }else{
          new_tr <- cur_tr
          if(! all_train){
            new_ts <- test_orig[, feat_family_list[[cur_fam]]]
          }
        }
      }
      
      
    }else{ # if there are no features in this family (after removing zero variance ones)
      new_tr <- matrix(nrow = nrow(train_orig), ncol = 0)
      if(! all_train){
        new_ts <- matrix(nrow = nrow(test_orig), ncol = 0)
        }
      
    }
    
    train_trans_list[[cur_fam]] <- new_tr
    rm(new_tr)
    if(! all_train){
      test_trans_list[[cur_fam]] <- new_ts
      rm(new_ts)
    }
  }
  #print("before comb")
  train_trans <- as.data.frame(do.call(cbind, train_trans_list))
  train_trans$label <- label_train
  stopifnot(sum(is.na(train_trans)) == 0)
  if(! all_train){
    test_trans <- as.data.frame(do.call(cbind, test_trans_list))
    test_trans$label <- label_test
    stopifnot(sum(is.na(test_trans)) == 0)
  }

  
  aanewname <- setdiff(colnames(train_trans), name_Dic[, 2])
  if(length(aanewname) > 0){
    name_Dic_new <- rbind(name_Dic, cbind(aanewname, aanewname))
  }else{
    name_Dic_new <- name_Dic
  }
  
  if(! all_train){
    my_res <- list(train_new = train_trans, 
                   test_new = test_trans,
                   name_dic_new = name_Dic_new,
                   irlba_List =irlba_list)
  }else{
    my_res <- list(train_new = train_trans, 
                   name_dic_new = name_Dic_new,
                   irlba_List =irlba_list)
  }
  return(my_res)
}


#load("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/Partition_6_CVfirst_dataset_2410003L11Rik.RData")
#load("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Input1/famSVD_inp1_dist.RData")
# load("~/Documents/Shayan/BioInf/lncRNA/Family_svd/Input1/famSVD_inp1_sequ_dist_Chrom_Meth_AX_ChIP_trnsp.RData")
#transformed_dataset <- feature_transformer(dataset = my_Dataset,
#                                           partitioning_vec = my_partition$CCV2,
#                                           name_Dic = my_name_dic,
#                                           feat_family_list = my_feat_family_list,
#                                           family_transform = my_family_transform,
#                                           family_trans_len = my_family_trans_len)
# sum(is.na(transformed_dataset$train_new))
# View(transformed_dataset)
# aatst <- distance_feature_partition_6[match(rownames(my_Dataset), names(distance_feature_partition_6))]
# aanm <- sort(unique(Partition_6_chunk_cv_df$owner))
# aa_file <- list.files(path = "~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset", full.names = T)
# sum(names(distance_feature_partition_6) == Partition_6_chunk_cv_df$tile_name)
# length(aa_file)
# for(i in 1:1){
#   load(aa_file[i])
#   aax <- distance_feature_partition_6[Partition_6_chunk_cv_df$owner %in% aanm[i]]
#   print(aanm[i])
#   print(sum(Partition_6_chunk_cv_df$owner %in% aanm[i]))
#   print(summary(aax))
#   print(sum(my_Dataset[, 1] == 0))
#   my_Dataset[, 1] <- aax
#   # save(list = c("my_Dataset", "my_partition", "my_name_dic"),
#   #      file = paste0("~/Documents/Shayan/BioInf/lncRNA/Partition_6_CVfirst_dataset/Partition_6_CVfirst_dataset_",
#   #                    aanm[i], ".RData"))
#   
# } 
 
 