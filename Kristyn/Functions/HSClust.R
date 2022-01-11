library(sClust)

euclidean_distance = function(X1, X2){
  return(sum((X1 - X2)^2))
}

getDistanceMatrix = function(X, distance_measure){
  p = ncol(X)
  M = matrix(0, p, p)
  for (j in 1:(p - 1)){
    for (k in (j + 1):p){
      val = distance_measure(X[, j], X[, k])
      M[j, k] = val
      M[k, j] = val
    }
  }
  max_M = max(M)
  S = M / max_M
  return(S)
}
getSimilarityMatrix = function(
  X = NULL, similarity_measure = NULL, distance_measure = NULL, 
  unnormalized_similarity_matrix = NULL, 
  normalized_distance_matrix = NULL, unnormalized_distance_matrix = NULL
){
  if(sum(c(!is.null(similarity_measure), !is.null(normalized_distance_matrix),
      !is.null(unnormalized_similarity_matrix), 
      !is.null(unnormalized_distance_matrix))) != 1){
    stop("getSimilarityMatrix: must give either a similarity measure, a distance measure, a normalized distance matrix, an unnormalized similarity matrix, OR an unormalized distance matrix!!")
  }
  if(!is.null(normalized_distance_matrix) | !is.null(distance_measure)){
    if(!is.null(distance_measure)){
      if(is.null(X)){
        stop("getSimilarityMatrix: need X to calculate similarity matrix!!")
      }
      normalized_distance_matrix = getDistanceMatrix(
        X = X, distance_measure = distance_measure)
    }
    return(1 - normalized_distance_matrix)
  } else if(!is.null(unnormalized_distance_matrix)){
    return(1 - (unnormalized_distance_matrix / 
                  max(unnormalized_distance_matrix)))
  } else if(!is.null(unnormalized_similarity_matrix)){
    return(unnormalized_similarity_matrix / max(unnormalized_similarity_matrix))
  } else{ # a similarity measure is given
    if(is.null(X)){
      stop("getSimilarityMatrix: need X to calculate similarity matrix!!")
    }
    p = ncol(X)
    M = matrix(0, p, p)
    for (j in 1:(p - 1)){
      for (k in (j + 1):p){
        val = similarity_measure(X[, j], X[, k])
        M[j, k] = val
        M[k, j] = val
      }
    }
    max_M = max(M)
    S = M / max_M
    S = S + diag(p)
    return(S)
  }
}

getNormalizedLaplacian = function(S){
  # S is similarity matrix (called A in TooManyCells)
  if(nrow(S) != ncol(S)) stop("S isn't pxp!")
  if(!isSymmetric(S)) stop("S isn't symmetric!")
  p = nrow(S)
  Dmat_diag = as.numeric(S %*% rep(1, p)) # assumed to be the matrix with degrees of nodes
  Dmat = diag(Dmat_diag)
  negsqrtD = diag(1 / sqrt(Dmat_diag))
  AdjMat = S # assume S is the adjacency matrix?
  normAdjMat = negsqrtD  %*% AdjMat %*% negsqrtD
  normLaplacMat = diag(p) - normAdjMat
  return(normLaplacMat)
}


spectral.clustering <- function(
  W, # similarity matrix
  n_eig = 2
){
  L = graph.laplacian(W)          # 2. compute graph laplacian
  ei = eigen(L, symmetric = TRUE) # 3. Compute the eigenvectors and values of L
  # we will use k-means to cluster the data
  # using the leading eigenvalues in absolute values
  ei$vectors <- ei$vectors[
    , base::order(abs(ei$values), decreasing=TRUE), drop = FALSE]
  if(nrow(W) == 2){
    obj <- kmeans(
      ei$vectors[, 2:n_eig, drop = FALSE], centers = n_eig, nstart = 100, 
      algorithm = "Lloyd")
  } else{
    obj <- kmeans(
      ei$vectors[, 2:n_eig, drop = FALSE], centers = n_eig, nstart = 100)
  }
  # if (n_eig==2){
  #   cl <- 2*(obj$cluster - 1) - 1
  # } else {
    cl <- obj$cluster
  # }
  names(cl) <- rownames(W)
  # return the cluster membership
  return(cl)
}

graph.laplacian <- function(
  W, 
  normalized = TRUE, 
  zeta=0.01
){
  stopifnot(nrow(W) == ncol(W))

  n = nrow(W)    # number of vertices
  # We perturb the network by adding some links with low edge weights
  W <- W + zeta * mean(colSums(W))/n * tcrossprod(rep(1,n))
  g <- colSums(W) # degrees of vertices

  if(normalized){
    D_half = diag(1 / sqrt(g) )
    return(D_half %*% W %*% D_half )
  } else {
    return(W)
  }
}



# binary tree code is inspired by sClust::hierClust
# https://mawenzi.univ-littoral.fr/sClust/software/
# https://rdrr.io/cran/sClust/src/R/recursClust.R
HSClust <- function(
  W, # similarity matrix
  levelMax = NULL, # p for full binary partition
  force_levelMax = TRUE, 
  stopping_rule = NULL, 
  #   NULL means go with force_levelMax
  #   "TooManyCells", "newmangirmanmodularity", "ngmod", "tmc", "ngm"
  method = "ShiMalik" # ShiMalik or kmeans
){
  p = nrow(W)
  if(is.null(levelMax)) levelMax = p
  if(is.null(stopping_rule)){
    stopping_rule = "none"
  } else{
    if(stopping_rule %in% 
       c("TooManyCells", "newmangirmanmodularity", "ngmod", "tmc", "ngm")){
      force_levelMax = FALSE
    }
  }
  
  stop <- FALSE
  # level of the clustering, i.e. how many cluster partitions have already been 
  #   made, plus 1
  level <- 1
  # clusterToCut = set of clusters to be further partitioned into more clusters
  clusterToCut <- 1
  # the clusters
  cl <- matrix(1, nrow = p, ncol = levelMax)
  
  while(!stop){
    newCluster <- c()
    # get the previous cluster assignments (base = 1, ..., 1)
    cl[,level + 1] <- cl[, level]
    # iterate over clusterToCut to further partition the clusters in that set
    for(i in 1:length(clusterToCut)){
      x = clusterToCut[i]
      indices = which(cl[, level] == x)
      Wprime = W[indices, indices]
      if(length(indices) > 1){
        # spectral clustering to partition into 2 clusters!
        if(method %in% c("ShiMalik", "shimalik", "sm", "SM")){
          invisible(capture.output(results <- ShiMalikSC(
            W = Wprime, flagDiagZero = FALSE, verbose = FALSE)))
          # groups <- results$cluster
          # # the above code might be choosing second largest instead of second 
          # #   smallest eigenvalue
          # NL = getNormalizedLaplacian(Wprime)
          # ei = eigen(NL, symmetric = TRUE)
          # second_smallest_idx = order(ei$values, decreasing = FALSE)[2]
          # groups = ifelse(ei$vectors[, second_smallest_idx] < 0, 1, 2)
          second_smallest_idx = order(results$eigenVal, decreasing = FALSE)[2]
          groups = ifelse(results$eigenVect[, second_smallest_idx] < 0, 1, 2)
        } else if(method %in% c("kmeans", "KMeans", "k", "K")){
          groups = spectral.clustering(W = Wprime, n_eig = 2)
        } else{
          stop("invalid method argument")
        }
        # Changing the cluster if necessary
        if(is.null(groups)){
          stop("HSClust(): groups is null for some reason! Clustering didn't happen correctly.")
        }
        if(length(unique(groups)) > 1){
          # cluster assignments for each observation
          if(level == 1){
            cluster_assignments <- as.character(groups)
          } else{
            cluster_assignments <- paste0(cl[indices, level], ".", groups)
          }
          cl[indices, level + 1] <- cluster_assignments
          # clusters to be further split
          newCluster <- c(newCluster, unique(cluster_assignments))
          if(stopping_rule %in% 
             c("TooManyCells", "newmangirmanmodularity", "ngmod", "tmc", "ngm")
             ){
            # stopping rule to be implemented here! --------------------------------
          }
        } else {
          # artificially split clusters if spectral clustering doesn't split
          #   and force_levelMax == TRUE
          if(length(groups) >= 2 && force_levelMax){
            groups = c(rep(1, length(groups) - 1), 2)
            if(length(groups) == p){
              cl[indices, level + 1] <- as.character(groups)
            } else{
              cl[indices, level + 1] <- paste0(cl[indices, level], ".", groups)
            }
            newCluster <- c(
              newCluster, unique(cl[indices, level + 1]))
          }
        }
      }
    }
    clusterToCut <- newCluster
    stop <- (level + 1) >= levelMax
    level <- level + 1
  }
  
  out <- list(
    cluster = cl[,level],
    allLevels = cl,
    nbLevels <- level
  )
}

sbp.fromHSClust = function(levels_matrix, row_names = NULL){
  p = nrow(levels_matrix)
  levels_matrix_old = levels_matrix # levels_matrix = levels_matrix_old
  levels_matrix = levels_matrix[, -1, drop = FALSE]
  sbp = matrix(ifelse(levels_matrix[, 1] == "1", 1, -1)) # "1" = 1, "2" = -1
  max_level_split_len_prev = 1
  if(ncol(levels_matrix) >= 2){
    for(j in 2:ncol(levels_matrix)){
      col.tmp = levels_matrix[, j] # get jth new cluster assignments (string)
      # check if any values are repeats from the previous column
      #   -- then these values aren't in this col's contrast (assign value 0)
      sbp.tmp0 = rep(NA, nrow(levels_matrix))
      col.prev.tmp = levels_matrix[, j - 1]
      sbp.tmp0[col.prev.tmp == col.tmp] = 0
      # make a matrix (nrow = p) for each of the unique columns in col.tmp
      #   that are the same (max) length -- two columns will form one contrast
      #   in the sbp matrix
      levels_all = unique(col.tmp)
      levels_split_lens = sapply(strsplit(levels_all, split = "\\."), length)
      max_level_split_len = max(levels_split_lens)
      if(max_level_split_len == max_level_split_len_prev){
        break
      }
      max_level_split_len_prev = max_level_split_len
      levels = levels_all[levels_split_lens == max_level_split_len]
      num_levels = length(levels)
      sbp_submat0 = matrix(rep(sbp.tmp0, num_levels), ncol = num_levels)
      # check each element in levels (there are p of them) to assign its
      #   corresponding contrast value to 1 or -1
      for(i in 1:p){
        # print(i)
        if(is.na(sbp.tmp0[i])){
          level.tmp = col.tmp[i]
          level_idx = which(levels == level.tmp)
          if(length(level_idx) == 1){
            sbp_submat0[i, level_idx] = ifelse(
              strsplit(level.tmp, split = "\\.")[[1]][j] == "1", 1, -1)
          } else if(!(length(level_idx) %in% c(0, 1))){
            stop("invalid clustering -- one covariate is in multiple clusters!")
          }
        }
      }
      # combine the columns that have the same previous level 
      #   (e.g. 1.1.1 & 1.1.2, 1.2.1, 1.2.2, 2.1.1, 2.1.2, 2.2.1, 2.2.2)
      # note: we can use matrix() bc all level-splits have
      #   length == max_level_split_len
      levels_split_mat = matrix(
        unlist(strsplit(levels, split = "\\.")), nrow = num_levels,
        byrow = TRUE)
      # combine columns that have same next_last_split_idx (should be 2 cols)
      #   because they define one contrast / split
      # get all unique previous levels (all those before the last level), 
      #   to loop through, e.g. c(1, 1), c(1, 2), c(2, 1), c(2, 2), corresp. to
      #   1.1, 1.2, 2.1, 2.2 respectively
      levels_prev = unique(
        levels_split_mat[, -ncol(levels_split_mat), drop = FALSE])
      sbp_submat = matrix(NA, nrow = p, ncol = nrow(levels_prev)) # ncol = ncol(sbp_submat0) / 2)
      for(l in 1:nrow(levels_prev)){
        # current previous levels to mathc
        prev_cur = levels_prev[l, ] # e.g. c(1, 1)
        # levels_prev cols with previous levels that match prev_cur --
        #   there should only be two of them
        col_idxs = which(apply(
          levels_split_mat[, -ncol(levels_split_mat), drop = FALSE], 1, 
          function(row) isTRUE(all.equal(unname(row), unname(prev_cur)))))
        if(length(col_idxs) != 2){
          print(col_idxs)
          stop("the number of columns for this split should be equal to 2!!")
        }
        sbp_submat0_subset_pt1 = sbp_submat0[, col_idxs[1]]
        sbp_submat0_subset_pt2 = sbp_submat0[, col_idxs[2]]
        sbp_col.tmp = sbp_submat0_subset_pt1
        sbp_col.tmp[!is.na(sbp_submat0_subset_pt2)] = sbp_submat0_subset_pt2[
          !is.na(sbp_submat0_subset_pt2)]
        sbp_submat[, l] = sbp_col.tmp
      }
      sbp = cbind(sbp, sbp_submat)
    }
  }
  sbp[is.na(sbp)] = 0
  if(!is.null(row_names)){
    rownames(sbp) = row_names
  } else{
    rownames(sbp) = paste("s", 1:nrow(sbp), sep = "")
  }
  colnames(sbp) = paste("z", 1:ncol(sbp), sep = "")
  return(sbp)
}







# plotting functions

getEdgesFromSBP = function(sbp){
  # make sure the rows (leaf nodes) and & columns (inner nodes) are named
  if(is.null(rownames(sbp))){
    rownames(sbp) = paste("s", 1:nrow(sbp), sep = "")
  }
  if(is.null(colnames(sbp))){
    colnames(sbp) = paste("z", 1:ncol(sbp), sep = "")
  }
  leafs = rownames(sbp)
  inners = colnames(sbp)
  
  # make edges data frame
  edges.df = data.frame(from = character(), to = character())
  for(j in 1:ncol(sbp)){
    from.tmp = rep(inners[j], 2) # because binary tree
    # left child - "1"'s
    if(sum(sbp[, j] == 1) == 1){ # if one "1", it's a leaf node
      left.tmp = leafs[which(sbp[, j] == 1)]
    } else{ # if multiple "1"'s, it's an inner node
      is.left.tmp = rep(FALSE, ncol(sbp))
      for(k in 1:ncol(sbp)){
        # there is only one column with the same locations of non-zeroes as 
        #   there are "1"'s in sbp[, j] -- that column has the left child node
        if(isTRUE(all.equal((sbp[, k] != 0), (sbp[, j] == 1)))){
          is.left.tmp[k] = TRUE
        }
      }
      left.tmp = inners[is.left.tmp]
    }
    # right child - "-1"'s
    if(sum(sbp[, j] == -1) == 1){ 
      right.tmp = leafs[which(sbp[, j] == -1)]
    } else{ 
      is.right.tmp = rep(FALSE, ncol(sbp))
      for(k in 1:ncol(sbp)){
        if(isTRUE(all.equal((sbp[, k] != 0), (sbp[, j] == -1)))){
          is.right.tmp[k] = TRUE
        }
      }
      right.tmp = inners[is.right.tmp]
    }
    to.tmp = c(left.tmp, right.tmp) # left & right children -- swap for ggraph
    if(length(to.tmp) == 2 & length(from.tmp) == 2){
      edges.df = rbind(edges.df, data.frame(from = from.tmp, to = to.tmp))
    }
  }
  return(edges.df)
}

# still need to figure out how to label!
plotSBP = function(
  sbp = NULL, edges = NULL, title = NULL, 
  nodes_types = NULL, 
  # a data frame of the nodes/vertices (column = "name") labeled as 
  #   "balance", "covariate", "significant covariate" (column = "type")
  text_size = 2
){
  # get edges data frame, if one isn't given
  if(is.null(sbp) & is.null(edges)){
    stop("plotSBP: provide either sbp or edges arguments!!")
  }
  if(is.null(edges) & !is.null(sbp)){
    edges = getEdgesFromSBP(sbp)
  }
  
  # make the graph input for ggraph plot
  if(!is.null(nodes_types)){
    # check nodes_types data frame
    if(!all(unname(nodes_types[, 1]) %in% unique(unlist(edges))) ||
       !all(unique(unlist(edges)) %in% unname(nodes_types[, 1]))){
      stop("names of nodes in nodes_types data frame don't match the row names of sbp matrix")
    } 
    mygraph <- as_tbl_graph(graph_from_data_frame(
      d = edges, 
      vertices = nodes_types
      ),
      directed = TRUE)
  } else{
    mygraph <- as_tbl_graph(graph_from_data_frame(
      d = edges, vertices = data.frame(nodes = unique(unlist(edges)))),
      directed = TRUE)
  }
  
  # make the plot
  plt = ggraph(mygraph, layout = 'igraph', algorithm = 'tree') +
    geom_edge_diagonal() +
    geom_node_point()
  if(!is.null(title)){
    plt = plt + ggtitle(title)
  }
  if(!is.null(nodes_types)){
    plt = plt + 
      geom_node_label(aes(label = name, fill = type), size = text_size)
  } else{
    plt = plt + 
      geom_node_label(aes(label = name), size = text_size)
  }
  plt = plt + theme_void()
  return(plt)
}
