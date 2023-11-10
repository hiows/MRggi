#' @name visGRN
#' @title GRN visualization
#'
#' @import visNetwork
#' @import RColorBrewer
#' @import colorspace
#' @import igraph
#'
#' @param nodes  genes
#' @param edges  gene-gene interaction
#'
#' @export

visGRN <- function(nodes, edges){
  
  #if(!require(visNetwork)){
  #  warning("instll visNetwork r package")
  #  install.packages("visNetwork")
  #} else if (!require(RColorBrewer)){
  #  install.packages("RColorBrewer")
  #} else if (!require(colorspace)){
  #  install.packages("colorspace")
  #} else if (!require(igraph)){
  #  install.packages("igraph")
  #}
  
  #require(visNetwork, quietly = TRUE)
  #require(RColorBrewer, quietly = TRUE)
  #require(colorspace, quietly = TRUE)
  #require(igraph, quietly = TRUE)
  
  # nodes color
  pal.nodes = "Set3"
  pal.idx.nodes = rev(qualitative_hcl(n = nlevels(as.factor(nodes$group)), palette = pal.nodes))  # sequential palette
  clu.idx = levels(as.factor(nodes$group))
  
  color = rep(0,nrow(nodes))
  for (i in 1:nlevels(as.factor(nodes$group))){
    color[which(nodes$group==clu.idx[i])] = pal.idx.nodes[i]
  }  
  nodes$color = color
  
  # edges color
  pal.edges = "Blue-Red"
  pal.idx = diverging_hcl(n = 7, palette = pal.edges)
  Bgg.idx = sort(c(Inf, as.vector(summary(abs(edges$Bgg))[c("1st Qu.","Median","3rd Qu.")]),
                   -as.vector(summary(abs(edges$Bgg))[c("1st Qu.","Median","3rd Qu.")]), -Inf))
  temp = edges$Bgg
  
  edges$color = sapply(1:nrow(edges), function(x){
    for (i in 1:7){
      if (temp[x]>=Bgg.idx[i] & temp[x]<Bgg.idx[i+1]){col = pal.idx[i]}
    }
    col
  })
  
  lnodes <- data.frame(label = clu.idx, color = pal.idx.nodes)
  
  visNetwork(nodes, edges) %>%
    visNodes(size = 10, shape = "circle") %>%
    visEdges(arrows ="to", width = 2) %>%
    visOptions(manipulation = TRUE)%>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visLegend(addNodes = lnodes, useGroups = F)
}
