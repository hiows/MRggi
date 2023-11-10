#' @name convertnet
#' @title convert MRggi results to network data
#'
#' @import igraph
#'
#' @param res  MRggi results
#' @param fdr.thr  FDR threshold
#' @param c.size  minimum size of cluster
#'
#' @return network
#' @export

convertnet <- function(res, fdr.thr = 0.05, c.size = 10){
  
  #if (!require(igraph)){
  #  install.packages("igraph")
  #}
  #require(igraph, quietly = T)
  
  sig_res = res[res$FDR_Bg1g2 < fdr.thr | res$FDR_Bg2g1 < fdr.thr,]
  
  message("filtering < " ,fdrth,"...")
  if (nrow(sig_res)==0){
    stop("There is no network.")
  }
  
  sig_res$fdrth.g1g2 = sapply(1:length(sig_res$FDR_Bg1g2), function(x){ if (sig_res$FDR_Bg1g2[x]<fdrth){1} else{0} })
  sig_res$fdrth.g2g1 = sapply(1:length(sig_res$FDR_Bg2g1), function(x){ if (sig_res$FDR_Bg2g1[x]<fdrth){1} else{0} })
  
  temp1=sig_res[sig_res$fdrth.g1g2 == 1,c("g1","g2","Bg1g2")]
  temp2=sig_res[sig_res$fdrth.g2g1 == 1,c("g2","g1","Bg2g1")]
  colnames(temp1) = colnames(temp2) = c("from", "to", "Bgg")
  edges = rbind(temp1, temp2)
  edges[,c(1,2)] <- apply(edges[,c(1,2)],2,as.character)
  
  nodes <- data.frame(id = sort(unique(c(edges$from,edges$to))), label = sort(unique(c(edges$from,edges$to))))
  message("nodes: ",nrow(nodes)," / edges: ",nrow(edges))
  
  
  # nodes clustering
  graph <- graph_from_data_frame(edges, directed = FALSE)
  cluster <- cluster_louvain(graph)
  
  cluster_df <- as.data.frame(t(data.frame(as.list(membership(cluster)))))
  cluster_df$label <- rownames(cluster_df)
  cluster_df = cluster_df[order(cluster_df$label),]
  temp = as.data.frame(table(cluster_df$V1))
  temp.idx = temp$Var1[ifelse(temp$Freq>c.size, T, F)]
  cluster_df$temp = ifelse(cluster_df$V1 %in% temp.idx, cluster_df$V1, 0)
  
  nodes = nodes[order(nodes$id),]
  nodes$temp = cluster_df$temp
  temp.group.idx = as.factor(nodes$temp)
  nodes$group = sapply(1:nrow(nodes), function(x){
    for (i in 1:nlevels(temp.group.idx)){
      if (nodes$temp[x] == as.double(levels(temp.group.idx))[i]){
        a = paste("cluster", i-1)
      } else if (nodes$temp[x] == 0){
        a = "-"
      }
    }
    a
  })
  nodes = nodes[,-3]
  
  network = list(edges = edges, nodes = nodes)
  
  return(network)
}
