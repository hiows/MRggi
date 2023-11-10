#' @name createSimulData
#' @title create simulation data set
#'
#' @import mvtnorm
#'
#' @param N sample size
#' @param Bs1g1 effect size between S1 and g1
#' @param Bs2g2 effect size between S1 and g1
#' @param Bg1g2 gene-gene interaction between g1 and g2
#'
#' @return res
#' @export

createSimulData = function(N = 200, Bs1g1 = c(0.55, 0.4, 0.25), Bs2g2 = c(0.55, 0.4, 0.25), Bg1g2 = 0, ro = 0){
  
  createSimulData.genotype = function(N = 200, Bs1g1 = c(0.55, 0.4, 0.25), Bs2g2 = c(0.55, 0.4, 0.25)){
    
    geno1 = rep(NULL, N)
    for (i in 1:length(Bs1g1)){
      snp = rbinom(N,2,0.3)
      geno1 = cbind(geno1, snp)
    }
    colnames(geno1) = paste0("snp",1:length(Bs1g1))
    
    geno2 = rep(NULL, N)
    for (i in 1:length(Bs2g2)){
      snp = rbinom(N,2,0.3)
      geno2 = cbind(geno2, snp)
    }
    colnames(geno2) = paste0("snp",(1 + length(Bs1g1)):(length(Bs2g2) + length(Bs1g1)))
    
    X = list(S1 = geno1, S2 = geno2)  # Genotype input data
    
    scale1 = scale(geno1)
    scale2 = scale(geno2)
    Bs1g1.S1 = sweep(scale1, 2, Bs1g1, "*")  # Bs1g1*S1
    Bs1g1.S2 = sweep(scale2, 2, Bs2g2, "*")  # Bs2g2*S2
    
    return(list(X = X, Bs1g1.S1 = Bs1g1.S1, Bs1g1.S2 = Bs1g1.S2))
  }
  
  creatSimulData.expression = function(N = 200, Bs1g1.S1, Bs1g1.S2, Bg1g2 = 0, ro = 0){
    u = rmvnorm(N, mean = c(0, 0), sigma = matrix(c(1, ro, ro, 1), 2))
    u1 = u[,1]
    u2 = u[,2]
    
    G1 = rowSums(Bs1g1.S1) + rnorm(N,0,1)    # Expression data 1: G1 = Bs1g1*S1 + e
    G2 = rowSums(Bs1g1.S2) + rnorm(N,0,1) + Bg1g2*G1  # Expression data 2: G2 = Bs2g2*S2 + Bg1g2*G1 + e
    
    if (ro != 0){
      G1 = G1 + u1
      G2 = G2 + u2
    }
    
    y = cbind(G1,G2)
    colnames(y) = c("G1","G2")
    
    return(y)
  }
  
  geno = createSimulData.genotype(N = N, Bs1g1 = Bs1g1, Bs2g2 = Bs2g2)
  expr = creatSimulData.expression(N = N, Bs1g1.S1 = geno$Bs1g1.S1, Bs1g1.S2 = geno$Bs1g1.S2, Bg1g2 = Bg1g2, ro = ro)
  
  res = list(genoList = geno$X, exprMat = expr)
  
  return(res)
}
