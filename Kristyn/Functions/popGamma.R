popGammajk = function(
  alpha1j, alpha1k, beta1, var_epsilon, var_epsilonj, var_epsilonk, U){
  varU = stats::var(U) #(1 / 12) * (0.5 - (-0.5))
  corrjk = ((alpha1j - alpha1k) * beta1 * varU) / 
    sqrt((beta1^2 * varU + var_epsilon) * 
           ((alpha1j - alpha1k)^2 * varU + var_epsilonj + var_epsilonk))
  return(abs(corrjk))
}
popGamma = function(
  alpha1, beta1, var_epsilon, var_epsilon2, U
){
  p = length(alpha1)
  if(length(var_epsilon2) == 1) var_epsilon2 = rep(var_epsilon2, p)
  
  rhoMat = matrix(0, p, p)
  for (j in 1:p){
    for (k in 1:p){
      if (k==j){next}
      else {
        rhoMat[j, k] = popGammajk(
          alpha1j = alpha1[j], alpha1k = alpha1[k], beta1 = beta1, 
          var_epsilon = var_epsilon, var_epsilonj = var_epsilon2[j], 
          var_epsilonk = var_epsilon2[k], U = U)
      }
    }
  }
  return(rhoMat)
}

slr_popGamma <- function(
  x, y, popGamma, subtractFrom1 = FALSE,
  num.clusters = 2, classification = FALSE, approx = TRUE, 
  amini.regularization = TRUE, 
  amini.regularization.parameter = 0.01, 
  highdegree.regularization.summary = "mean",
  highdegree.regularization = FALSE,
  include.leading.eigenvector = TRUE
){
  n = nrow(x)
  p = ncol(x)
  
  ## Compute pairwise correlation
  if(!subtractFrom1){
    rhoMat <- popGamma
  } else{
    rhoMat <- 1 - popGamma
  }
  out <- list()
  out$kernel <- rhoMat
  
  ## Split into active/inactive sets
  if(approx){
    rhoMat.svd <- svd(rhoMat)
    rhoMat_approx <- tcrossprod(
      rhoMat.svd$u[, 1], 
      rhoMat.svd$v[, 1]) *
      rhoMat.svd$d[1]
    rownames(rhoMat_approx) <- colnames(rhoMat_approx) <- rownames(rhoMat)
    clusters1 <- spectral.clustering2(
      rhoMat_approx, n_eig = num.clusters, 
      amini.regularization = amini.regularization,
      amini.regularization.parameter = amini.regularization.parameter,
      highdegree.regularization.summary = highdegree.regularization.summary,
      highdegree.regularization = highdegree.regularization,
      include.leading.eigenvector = include.leading.eigenvector)
  } else{
    clusters1 <- spectral.clustering2(
      rhoMat, n_eig = num.clusters, 
      amini.regularization = amini.regularization,
      amini.regularization.parameter = amini.regularization.parameter,
      highdegree.regularization.summary = highdegree.regularization.summary,
      highdegree.regularization = highdegree.regularization,
      include.leading.eigenvector = include.leading.eigenvector)
  }
  cluster.lengths = table(clusters1)
  clusters = names(cluster.lengths)
  out$num.clusters = num.clusters
  
  ## Find which set has the more predictive balance
  cors = rep(NA, num.clusters)
  Rsqs = rep(NA, num.clusters)
  sbp.ests = matrix(0, nrow = p, ncol = num.clusters)
  bal.ests = matrix(NA, nrow = n, ncol = num.clusters)
  rownames(sbp.ests) <- colnames(x)
  for(i in 1:num.clusters){ 
    cluster.label = as.numeric(clusters[i])
    index = which(clusters1 == cluster.label)
    if(length(index) >= 2){ # if a log-ratio can be made from the variables in this cluster
      ## Perform spectral clustering to get the numerator/denominator groups
      if(all(rhoMat[index, index] == 0)){
        sbp.ests[, i] = rep(NA, p)
      } else{
        subset = spectral.clustering2(
          rhoMat[index, index], n_eig = 2, reindex = TRUE, 
          amini.regularization = amini.regularization,
          amini.regularization.parameter = amini.regularization.parameter,
          highdegree.regularization.summary = highdegree.regularization.summary, 
          highdegree.regularization = highdegree.regularization,
          include.leading.eigenvector = include.leading.eigenvector)
        ## calculate balance and its correlation with y
        sbp.ests[match(names(subset), rownames(sbp.ests)), i] = subset
        bal.ests[, i] = balance::balance.fromSBP(
          x = x, y = sbp.ests[, i, drop = FALSE])
      }
    } else{ # otherwise, if there's just one variable in the cluster, set all other variables as -1
      sbp.ests[, i] = rep(NA, p)
    }
    cors[i] = stats::cor(bal.ests[, i], y)
    Rsqs[i] = cors[i]^2
  }
  out$cors = cors
  out$Rsqs = Rsqs
  # The correct active set should have largest correlation magnitude, i.e. Rsq.
  ## We refit the linear model on the balance from the set with the 
  ##    largest correlation magnitude.
  selected.cluster = which.max(Rsqs)
  out$index = sbp.ests[, selected.cluster]
  refit_data = data.frame(V1 = bal.ests[, selected.cluster], y = y)
  if(!classification){
    refit <- lm(y~V1, data = refit_data)
  } else{
    refit = stats::glm(
      y~V1, data = refit_data, family = binomial(link = "logit"))
  }
  out$model <- refit
  
  # return the full SBP vector (with entries for all p variables)
  full.sbp.est = sbp.ests[, selected.cluster, drop = FALSE]
  rownames(full.sbp.est) = colnames(x)
  out$sbp = full.sbp.est
  
  # return the slr object
  return(out)
}