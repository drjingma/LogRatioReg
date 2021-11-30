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

# getNormalizedLaplacian = function(S){
#   # S is similarity matrix (called A in TooManyCells)
#   if(nrow(S) != ncol(S)) stop("S isn't pxp!")
#   if(!isSymmetric(S)) stop("S isn't symmetric!")
#   p = nrow(S)
#   Dmat_diag = as.numeric(S %*% rep(1, p)) # assumed to be the matrix with degrees of nodes
#   Dmat = diag(Dmat_diag)
#   negsqrtD = diag(1 / sqrt(Dmat_diag))
#   AdjMat = S # assume S is the adjacency matrix?
#   normAdjMat = negsqrtD  %*% AdjMat %*% negsqrtD
#   normLaplacMat = diag(p) - normAdjMat
#   return(normLaplacMat)
# }


# hierarchical spectral clustering recursively using sClust
# https://mawenzi.univ-littoral.fr/sClust/software/
# https://rdrr.io/cran/sClust/src/R/recursClust.R
# hierClust from library sClust
HSClust <- function(
  # dataFrame, 
  W, # similarity matrix
  levelMax = NULL, # p - 1
  # clustFunction = ShiMalikSC, 
  # similarity = TRUE, 
  flagDiagZero = FALSE,
  # biparted = FALSE, 
  minPoint = 1
){
  similarity = FALSE
  biparted = TRUE
  if(is.null(levelMax)) levelMax = nrow(W)
  
  stop <- FALSE
  level <- 1
  clusterToCut <- 1
  cl <- matrix(1, nrow = nrow(W), ncol = levelMax)
  
  while(!stop){
    newCluster <- c()
    cl[,level+1] <- cl[,level]
    sapply(clusterToCut, FUN = function(x){
      indices = which(cl[,level]==x)
      Wprime = W[indices, indices]
      if(length(indices) > minPoint){
        invisible(capture.output(results <- ShiMalikSC(
          W = Wprime, flagDiagZero = flagDiagZero, verbose = FALSE)))
        groups <- results$cluster
        #Changing the cluster if necessary
        if(!is.null(groups) && length(unique(groups))>1){
          if(level == 1 ){
            cl[indices,level+1] <<- paste0(groups)
            newCluster <<- c(newCluster, unique(paste0(groups)))
          }else{
            cl[indices,level+1] <<- paste0(cl[indices,level],".",groups)
            newCluster <<- c(newCluster, unique(paste0(cl[indices,level],".",groups)))
          }
        }
      }
    })
    
    clusterToCut <- newCluster
    stop <- (level+1) >= levelMax
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
  levels_matrix = levels_matrix[, -1, drop = FALSE]
  sbp = matrix(ifelse(levels_matrix[, 1] == "1", 1, -1)) # "1" = 1, "2" = -1
  if(ncol(levels_matrix) >= 2){
    for(j in 2:ncol(levels_matrix)){
      col.tmp = levels_matrix[, j]
      sbp.tmp0 = rep(NA, nrow(levels_matrix))
      # check if any values are repeats from the previous column
      col.prev.tmp = levels_matrix[, j - 1]
      sbp.tmp0[col.prev.tmp == col.tmp] = 0
      # separate strings by period
      levels = strsplit(col.tmp, split = "\\.")
      sbp.tmp = sbp.tmp0
      for(i in 1:nrow(levels_matrix)){
        levels.tmp = levels[[i]]
        # check if previous levels have changed
        if(i > 1){
          levels.prev.tmp = levels[[i - 1]]
          if((length(levels.tmp[-length(levels.tmp)]) != 
             length(levels.prev.tmp[-length(levels.prev.tmp)])) || 
              (!all(levels.tmp[-length(levels.tmp)] == 
                   levels.prev.tmp[-length(levels.prev.tmp)]))){
            if(!all(na.omit(sbp.tmp) == 0)){
              sbp = cbind(sbp, sbp.tmp)
            }
            sbp.tmp = sbp.tmp0
          } 
        }
        # assign 1, -1 for last level
        if(is.na(sbp.tmp0[i])){
          sbp.tmp[i] = ifelse(levels.tmp[j] == "1", 1, -1)
        }
        if(i == nrow(levels_matrix)){
          if(!all(na.omit(sbp.tmp) == 0)){
            sbp = cbind(sbp, sbp.tmp)
          }
        }
      }
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
    edges.df = rbind(edges.df, data.frame(from = from.tmp, to = to.tmp))
  }
  return(edges.df)
}

# still need to figure out how to label!
plotSBP = function(sbp = NULL, edges = NULL){
  if(is.null(sbp) & is.null(edges)){
    stop("plotSBP: provide either sbp or edges arguments!!")
  }
  if(is.null(edges) & !is.null(sbp)){
    edges = getEdgesFromSBP(sbp)
  }
  mygraph <- graph_from_data_frame(
    d = edges, vertices = data.frame(labels = unique(unlist(edges))),
    directed = TRUE)
  plt = ggraph(mygraph, layout = 'igraph', algorithm = 'tree') +
    geom_edge_diagonal() +
    geom_node_point() +
    # geom_node_label(aes(label = vertices)) + #c(edges$from[1], edges$to))) +
    theme_void()
  return(plt)
}
