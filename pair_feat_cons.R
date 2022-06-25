#create pair features based on a provided network
add_pair_features_network <- function(lncRNA_feat, tile_feat, my_network_list = numeric(0), presence_thr=0){
  # lncRNA_feat is the feature set for lncRNA of interest --> numeric
  # tile_feat is the feature set for the tiles (matrix with ncol(tile_feat) == length(lncRNA_feat))
  # my_network_list is a set of edges over which the aggregation takes place, if nothing is porvided (and if name name_pattern is not provided) homo pairing of all feature will be done
  #  each entry of my_network_list is an n*2 matrix, where n is the number of interactions and each column shows the index of interacting protein in the provided feature set
  # note that each interaction should be shown with two entries ((a , b) and (b , a))
  stopifnot(ncol(tile_feat) == length(lncRNA_feat)
            ,all(names(lncRNA_feat) == colnames(tile_feat))
            )
  
  if(length(my_network_list) == 0){
    my_network_list = list()
    my_network_list[[1]] <- cbind(c(1: length(lncRNA_feat)),
                                  c(1: length(lncRNA_feat)))
    names(my_network_list) <- "All_Homo_pairs"
  }else{
    stopifnot(all(unlist(lapply(my_network_list, is.matrix))))
  }
  pair_feat <- matrix(nrow = nrow(tile_feat), ncol = length(my_network_list))
  colnames(pair_feat) <- names(my_network_list)
  #pair_feat <- numeric(length = nrow(tile_feat))
  for(cur_tile in 1:nrow(tile_feat)){
    print(cur_tile)
    for(cur_net in 1:length(my_network_list)){
    pair_feat[cur_tile, cur_net] <- sum((lncRNA_feat[my_network_list[[cur_net]][, 1]] > presence_thr) &
                                          (tile_feat[cur_tile, my_network_list[[cur_net]][, 2]] > presence_thr))
    # pair_feat[cur_tile] <- sum((lncRNA_feat > presence_thr) &
    #                              (tile_feat[cur_tile,] > presence_thr))
    }
  }
  return(pair_feat)
}
