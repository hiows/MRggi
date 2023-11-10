#' @name .multiple_regression
#' @title performing multiple regression
#'
#' @param y1 expression
#' @param X1 genotype
#'
#' @return beta

.multiple_regression = function(y1, X1){
  beta_fit = summary(lm(y1 ~ X1))
  beta_coef = coef(beta_fit)
  beta = beta_coef[,1][-1]
  
  return(beta)
}


#' @name .TSLS
#' @title performing two-stage least square
#'
#' @param X1 genotype for gene1
#' @param X2 genotype for gene2
#' @param y1 expression for gene1
#' @param y2 expression for gene2
#'
#' @return stat

.TSLS = function(X1, X2, y1, y2){
  
  # IV -> exposure
  Bzx_fit = summary(lm(y1 ~ X1))
  Bzx_coef = coef(Bzx_fit)
  Bzx = Bzx_coef[,1][-1]
  
  # IV -> outcome
  y2.resid = resid(lm(y2 ~ X2))
  Bzy_fit = summary(lm(y2.resid ~ X1))
  Bzy_coef = coef(Bzy_fit)
  Bzy = Bzy_coef[,1][-1]
  se_Bzy = Bzy_coef[,2][-1]
  IV = 1/(se_Bzy^2)
  
  # outcome -> exposure
  Bxy = sum(Bzx*Bzy*IV)/sum(Bzx*Bzx*IV)
  se_Bxy = sqrt(1/sum(Bzx*Bzx*IV))
  pval = 2*pnorm(abs(Bxy/se_Bxy), lower.tail = F)  # because of two-sided test
  
  stat = list(Bxy = Bxy, se_Bxy = se_Bxy, Pvalue = pval)
  
  return(stat)
}


#' @name MRggi
#' @title MR-GGI
#'
#' @param y gene expression; matrix
#' @param X cis-SNPs genotype for each gene; list
#' @param cor.thr threshold of the correlation between genes
#' @param Bsg.thr threshold of the effect size between cis-Snp set and its gene
#' @param p.adjust.method method for p value adjustment
#'
#' @return df
#' @export

MRggi = function(y, X, cor.thr = 0, p.adjust.method = "bonferroni"){
  
  # output
  g1 = c()
  g2 = c()
  GGcor = c()
  Bg1g2 = c()
  Bg2g1 = c()
  pval_Bg1g2 = c()
  pval_Bg2g1 = c()
  
  # data scale
  scale.y = scale(y)
  scale.X = lapply(X, scale)
  
  # filtering with threshold
  cor.y = cor(scale.y)
  corMat = matrix(0, nrow(cor.y), ncol(cor.y))
  corMat[which(abs(cor.y) > cor.thr)] = 1
  corMat[lower.tri(corMat, diag = T)] = 0
  calc.idx = which(corMat == 1, arr.ind = T)
  calc.idx = as.data.frame(calc.idx)
  
  # message
  message("gene-gene correlation : >", cor.thr)
  message("Possible network : ",nrow(calc.idx))
  
  for (k in 1:nrow(calc.idx)){  # Calculate by length temp
    
    # Get data index
    i = calc.idx[k,1]
    j = calc.idx[k,2]
    y1 = y[,i]
    y2 = y[,j]
    X1 = X[[i]]
    X2 = X[[j]]
    
    # Bzx
    Bs1g1 = .multiple_regression(y1, X1)
    Bs2g2 = .multiple_regression(y2, X2)
    
    Bg1g2.data = .TSLS(X1, X2, y1, y2)
    Bg2g1.data = .TSLS(X2, X1, y2, y1)
    
    # data append
    g1 = append(g1, colnames(y)[i])
    g2 = append(g2, colnames(y)[j])
    GGcor = append(GGcor, cor(y1, y2)) 
    Bg1g2 = append(Bg1g2, Bg1g2.data$Bxy)
    pval_Bg1g2 = append(pval_Bg1g2, Bg1g2.data$Pvalue)
    Bg2g1 = append(Bg2g1, Bg2g1.data$Bxy)
    pval_Bg2g1 = append(pval_Bg2g1, Bg2g1.data$Pvalue)
  }
  
  FDR_Bg1g2 = rep(0,length(g1))
  g1.lv = levels(as.factor(g1))
  for (i in g1.lv){
    idx = which(g1 == i)
    pval.idx = pval_Bg1g2[idx]
    FDR.idx = p.adjust(pval.idx, method = p.adjust.methods)
    FDR_Bg1g2[idx] = FDR.idx
  }
  FDR_Bg2g1 = rep(0,length(g2))
  g2.lv = levels(as.factor(g2))
  for (i in g2.lv){
    idx = which(g2 == i)
    pval.idx = pval_Bg2g1[idx]
    FDR.idx = p.adjust(pval.idx, method = p.adjust.methods)
    FDR_Bg2g1[idx] = FDR.idx
  }
  
  df = data.frame(g1 = g1,
                  g2 = g2,
                  GGcor = round(GGcor, 3),
                  Bg1g2 = round(Bg1g2, 3),
                  pval_Bg1g2 = round(pval_Bg1g2, 3),
                  FDR_Bg1g2 = round(FDR_Bg1g2, 3),
                  Bg2g1 = round(Bg2g1, 3),
                  pval_Bg2g1 = round(pval_Bg2g1, 3),
                  FDR_Bg2g1 = round(FDR_Bg2g1, 3))
  return(df)
}
