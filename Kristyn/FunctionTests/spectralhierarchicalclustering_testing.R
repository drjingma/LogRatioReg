# spectral decomposition
library(matlib)

# mat <- matrix(c(1, 2, 3, 2, 7, 4, 3, 4, 0), nrow=3)
# e <- eigen(mat)
# o <- order(e$values, decreasing=FALSE)
# e$values[o]

# hierarchical spectral clustering
library(matlib)

X = matrix(c(1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 4, 4, 5, 4, 4, 4, 5, 4, 4, 4.5, 3, 3, 3, 3, 2, 2, 2, 1), nrow = 4)
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
  X = NULL, similarity_measure = NULL, distance_matrix = NULL
){
  if((is.null(similarity_measure) & is.null(distance_matrix)) | 
     (!is.null(similarity_measure) & !is.null(distance_matrix))){
    stop("getSimilarityMatrix: must give either a similarity measure or a distance matrix, but not both!!")
  }
  if(!is.null(distance_matrix)){
    return(1 - distance_matrix)
  } else{
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

# get2ClustersOld = function(S){
#   normLaplacMat = getNormalizedLaplacian(S)
#   eigendecomp = eigen(normLaplacMat)
#   # find eigenvector corresponding to second smallest eigenvalue (Fiedler value)
#   #   https://towardsdatascience.com/spectral-clustering-aba2640c0d5b
#   eigenorder = order(eigendecomp$values, decreasing = FALSE) # small -> large
#   eigenvector_2ndsmallest = eigendecomp$vectors[, eigenorder[2]]
#   cluster_assignments = eigenvector_2ndsmallest >= 0
#   return(cluster_assignments)
# }
# hierSpectralClust = function(X, similarity_measure){
#   p = ncol(X)
#   sbp = matrix(NA, nrow = p, ncol = p - 1)
#   S0 = similarity_measure(X)
#   sbp[, 1] = get2Clusters(S0)
#   left_cluster0 = X[, sbp[, i]]
#   sbp[, 2] = get2Clusters(S0[ sbp[, i], sbp[, i]])
#   right_cluster0 = X[, -sbp[, i]]
#   sbp[, 3] = get2Clusters(S0[ -sbp[, i], -sbp[, i]])
#   
#   for(i in 1:(p - 1)){
#     S.tmp = similarity_measure(X)
#     res = get2Clusters(S.tmp)
#     sbp[, i] = 
#   }
# }

# # recursive function example
# factorial = function(k){
#   if(k == 1){
#     return(1)
#   } else{
#     return(k * factorial(k - 1))
#   }
# }
# factorial(5)
# 
# # roll out into for loop
# k = 5
# factorialk = 1
# for(i in 1:k){
#   factorialk = factorialk * i
# }

# get2Clusters = function(S_full, all_indices, subset_indices, clustering_list){
#   if(is.null(subset_indices)) subset_indices = all_indices
#   if(length(subset_indices) < 2) return(clustering_list)
#   S_subset = S_full[ subset_indices, subset_indices]
#   normLaplacMat = getNormalizedLaplacian(S_subset)
#   eigendecomp = eigen(normLaplacMat)
#   # find eigenvector corresponding to second smallest eigenvalue (Fiedler value)
#   #   https://towardsdatascience.com/spectral-clustering-aba2640c0d5b
#   eigenorder = order(eigendecomp$values, decreasing = FALSE) # small -> large
#   eigenvector_2ndsmallest = eigendecomp$vectors[, eigenorder[2]]
#   cluster_assignments = eigenvector_2ndsmallest >= 0
#   
#   # find the parent
#   parent_idx = NA
#   for(i in 1:length(clustering_list)){
#     if(clustering_list[[i]]$value = subset_indices) parent_idx = i
#   }
#   # add the left and right children to the clustering_list
#   left_child_idx = length(clustering_list) + 1
#   clustering_list[[left_child_idx]] = list(
#     index = left_child_idx, 
#     position = "left",
#     parent_index = parent_idx,
#     value = subset_indices[cluster_assignments]
#   )
#   right_child_idx = length(clustering_list) + 1
#   clustering_list[[right_child_idx]] = list(
#     index = right_child_idx, 
#     position = "right",
#     parent_index = parent_idx,
#     value = subset_indices[-cluster_assignments]
#   )
#   # find the next left and right children for each left and right child above
#   clustering_list[[length(clustering_list) + 1]] = list(
#     index = length(clustering_list) + 1, 
#     position = "left",
#     parent_index = left_child_idx,
#     value = subset_indices[cluster_assignments]
#     
#   )
#   clustering_list[[length(clustering_list) + 1]] = list(
#     index = length(clustering_list) + 1, 
#     position = "right",
#     parent_index = right_child_idx,
#     value = subset_indices[-cluster_assignments]
#   )
#   return(clustering_list)
# }
# 
# hierSpectralClust = function(X, similarity_measure){
#   p = ncol(X)
#   sbp = matrix(NA, nrow = p, ncol = p - 1)
#   S_full = getSimilarityMatrix(
#     distance_matrix = getDistanceMatrix(
#       X = X, distance_measure = euclidean_distance))
#   clustering_list = list()
#   clustering_list[[1]] = list(
#     index = 0,
#     position = "root", 
#     parent_index = NA,
#     value = 1:p
#   )
#   sbp[, 1] = get2Clusters(S_full, all_indices = 1:p, subset_indices = NULL, clustering_list = NULL)
#   
#   
#   left_cluster0 = X[, sbp[, i]]
#   sbp[, 2] = get2Clusters(S0[ sbp[, i], sbp[, i]])
#   right_cluster0 = X[, -sbp[, i]]
#   sbp[, 3] = get2Clusters(S0[ -sbp[, i], -sbp[, i]])
# 
#   for(i in 1:(p - 1)){
#     S.tmp = similarity_measure(X)
#     res = get2Clusters(S.tmp)
#     sbp[, i] =
#   }
# }

# hierarchical spectral clustering recursively using sClust
library(sClust)
#' @title Perform a multi level clustering
#' @author Emilie Poisson Caillault and Erwan Vincent
#' @description The function, for a given dataFrame, will separate the data using the input clustering method in several levels.
#' @param dataFrame The dataFrame.
#' @param levelMax The maximum depth level.
#' @param clustFunction the clustering function to apply on data.
#' @param similarity if True, will use the similarity matrix for the clustering function.
#' @param vois number of points that will be selected for the similarity computation. 
#' @param flagDiagZero if True, Put zero on the similarity matrix W.
#' @param biparted if True, the function will not automatically choose the number of clusters to compute.
#' @param method The method that will be used. "default" to let the function choose the most suitable method. "PEV" for the Principal EigenValue method. "GAP" for the GAP method.
#' @param tolerence The tolerance allowed for the Principal EigenValue method.
#' @param threshold The threshold to select the dominant eigenvalue for the GAP method.
#' @param minPoint The minimum number of points required to compute a cluster.
#' @param verbose To output the verbose in the terminal. 
#' @param ... additional arguments for the clustering function.
#' @return returns a list containing the following elements:
#' \itemize{
#'  \item{cluster: }{vector that contain the result of the last level}
#'  \item{allLevels: }{dataframe containing the clustering results of each levels}
#'  \item{nbLevels: }{the number of computed levels}
#' }
#' @examples
#' ### Example 1: 2 disks of the same size
#' n<-100 ; r1<-1
#' x<-(runif(n)-0.5)*2;
#' y<-(runif(n)-0.5)*2
#' keep1<-which((x*2+y*2)<(r1*2))
#' disk1<-data.frame(x+3*r1,y)[keep1,]
#' disk2 <-data.frame(x-3*r1,y)[keep1,]
#' sameTwoDisks <- rbind(disk1,disk2)
#' res <- recursClust(scale(sameTwoDisks),levelMax=3, clustFunction =ShiMalikSC,
#'                    similarity = TRUE, vois = 7, flagDiagZero = FALSE,
#'                    biparted = TRUE, verbose = TRUE)
#' plot(sameTwoDisks, col = as.factor(res$cluster))
#' 
#' ### Example 2: Speed and Stopping Distances of Cars
#' res <- recursClust(scale(iris[,-5]),levelMax=4, clustFunction = spectralPAM,
#'                    similarity = TRUE, vois = 7, flagDiagZero = FALSE,
#'                    biparted = FALSE, method = "PEV", tolerence =  0.99,
#'                    threshold = 0.9, verbose = TRUE)
#' plot(iris, col = as.factor(res$cluster))
HSClust <- function(
  # dataFrame, 
  W, 
  levelMax = 2, 
  # clustFunction = ShiMalikSC, 
  # similarity = TRUE, 
  flagDiagZero = FALSE,
  # biparted = FALSE, 
  minPoint = 1,
  verbose = TRUE
){
  similarity = FALSE
  biparted = TRUE
  
  stop <- FALSE
  level <- 1
  clusterToCut <- 1
  cl <- matrix(1, nrow = nrow(W), ncol = levelMax)
  
  while(!stop){
    newCluster <- c()
    cl[,level+1] <- cl[,level]
    message(level)
    sapply(clusterToCut, FUN = function(x){
      print(x)
      indices = which(cl[,level]==x)
      Wprime = W[indices, indices]
      if(length(indices) > minPoint){
        results <- ShiMalikSC(
          W = Wprime, flagDiagZero = flagDiagZero, verbose = verbose)
        groups <- results$cluster
        print(groups)
        #Changing the cluster if necessary
        if(!is.null(groups) && length(unique(groups))>1){
          if(level == 1 ){
            cl[indices,level+1] <<- paste0(groups)
            newCluster <<- c(newCluster, unique(paste0(groups)))
          }else{
            cl[indices,level+1] <<- paste0(cl[indices,level],".",groups)
            newCluster <<- c(newCluster, unique(paste0(cl[indices,level],".",groups)))
          }
          if(verbose){
            message("CLUSTERS VECTOR")
            print(cl[,level+1])
          }
          print(newCluster)
          
        }
      }
    })
    
    clusterToCut <- newCluster
    print(clusterToCut)
    stop <- (level+1) >= levelMax
    level <- level + 1
  }
  
  out <- list(
    cluster = cl[,level],
    allLevels = cl,
    nbLevels <- level
  )
}

p = ncol(X)
res = HSClust(
  W = getSimilarityMatrix(
    distance_matrix = getDistanceMatrix(
      X = X, distance_measure = euclidean_distance)), 
  levelMax = p - 1, 
  flagDiagZero = TRUE,
  verbose = TRUE
)

res$cluster
res$allLevels

levels_matrix = res$allLevels

sbp.fromHSClust = function(levels_matrix){
  levels_matrix = levels_matrix[, -1]
  sbp = ifelse(levels_matrix[, 1] == "1", 1, -1) # "1" = 1, "2" = -1
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
        if(!all(levels.tmp[-length(levels.tmp)] == 
                levels.prev.tmp[-length(levels.prev.tmp)])){
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
  sbp[is.na(sbp)] = 0
  rownames(sbp) = colnames(sbp) = NULL
  return(sbp)
}
sbp.fromHSClust(res$allLevels)

# res = recursClust(
#   dataFrame = getSimilarityMatrix(
#     distance_matrix = getDistanceMatrix(
#       X = X, distance_measure = euclidean_distance)), # t() to cluster the variables, not the observations
#   levelMax = p - 1, clustFunction = ShiMalikSC, 
#   similarity = FALSE, # if TRUE, use similarity matrix for clustering function
#   # W = getSimilarityMatrix(
#   #   distance_matrix = getDistanceMatrix(
#   #     X = X, distance_measure = euclidean_distance)),
#   # vois = 2, # (default) # pts that will be selected for similarity computation
#   biparted = TRUE, method = "default", minPoint = 1, 
#   verbose = TRUE
# )


