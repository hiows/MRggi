#' @name FineMapping
#' @title fine-mapping
#'
#' @import susieR
#'
#' @param X genotype matrix
#' @param y expression matrix
#' @param L susie param
#' @param pip.thr posterior inclusion probabilities thr
#'
#' @return X.FineMap
#' @export

FineMapping = function(X, y, L = 10, pip.thr = 0){
  # X : genotype matrix
  # y : expression matrix
  
  X = sapply(1:ncol(X), function(i) as.double(X[,i]))
  
  X.FineMap = lapply(1:ncol(y), function(i){
    
    res = susie(X, y[,i])  # fine-mapping
    res.cs = res$sets$cs
    res.pip = res$pip
    res.pip.idx = which(res.pip > pip.thr)
    
    snp.idx = unlist(lapply(1:length(res.cs), function(x){
      temp1 = res.cs[[x]][res.cs[[x]] %in% res.pip.idx]
      
      if (length(temp1) > 0){
        temp1.pip = res.pip[temp1]
        temp1.max.pip.idx = which(temp1.pip == max(temp1.pip))
        temp1.snp.idx = temp1.max.pip.idx[which.max(temp1.max.pip.idx)]
        a = temp1[temp1.snp.idx]
        a
      }
    }))  # cs > pip.thr
    
    if (is.null(snp.idx)){
      snp.idx = max(which(res.pip==max(res.pip)))
    }  # cs < pip.thr
    
    if (is.null(res.cs)){
      if (max(res.pip) != 0){
        snp.idx = max(which(res.pip==max(res.pip)))
      } else if (max(res.pip) == 0){
        stat = univariate_regression(X, y[,i])
        z = stat$betahat/stat$sebetahat
        snp.idx = max(which(z==max(z)))
      }
    }  # no sets
    
    as.matrix(X[,snp.idx])
  }) ; names(X.FineMap) = paste0("S",1:length(X.FineMap))
  
  X.FineMap
  
  return(X.FineMap)
}