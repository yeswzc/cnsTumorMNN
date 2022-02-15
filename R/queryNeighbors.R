#@ x.ref: PCA of reference data
#@ x.query: PCA of query sample
#@ k: number of K neighbors to query
#@ required package RANN, plotly
#@ return a list of 2 tables and 1 figure list
#' 
#' @export
#' @import dplyr
#' @import RColorBrewer
#' @import RANN
knn.search = function(x.ref, x.query, y.ref = y.ref, k = 15){
  #library(dplyr)
  paired.color = RColorBrewer::brewer.pal(12,"Paired")
  nn2.search = RANN::nn2(x.ref, query = x.query, k = k, treetype = "kd")
  idx = nn2.search$nn.idx
  idx.label = as.character(y.ref[nn2.search$nn.idx])
  idx.dist = as.vector(nn2.search$nn.dists)
  idx.dist.qc = vector(length = length(idx.label))
  idx.dist.threshold = vector(length = length(idx.label))
  #
  candidate = unique(idx.label)
  for(x in 1:length(candidate)){
    class= candidate[x]
    #class = "GBM, G34"
    dist.threshold = min(as.vector( c(quantile(ref_distance.subclass[[class]][[1]], 0.75), 
                                      quantile(ref_distance.subclass[[class]][[2]], 0.01) )))
    #cat(class, dist.threshold,"\n")
    ii = which(idx.label == class)
    idx.dist.qc[ii] = ifelse(idx.dist[ii] < dist.threshold, "Y", "N")
    idx.dist.threshold[ii] = dist.threshold
  }
  
  res = data.frame(ID = row.names(x.ref)[idx], Cancer = idx.label, Distance = round(idx.dist,1), 
                   Threshold = round(idx.dist.threshold,1), PASS = idx.dist.qc, stringsAsFactors=F)
  ##figures
  pass.subclass = as.vector(unique(res$Cancer[which(res$PASS=="Y")]))
  if(length(pass.subclass) < 1){
      #print("Not have near neighbors")
      NN = ifelse(length(unique(res$Cancer))>3, 3, length(unique(res$Cancer))); 
      pass.subclass = unique(res$Cancer)[1:NN];
   }
   #print(pass.subclass)
  #
  figs = lapply(pass.subclass,function(class){
    #class = pass.subclass[1]
    #print(class)
    #print(class %in% names(ref_distance.subclass) )
    dens1 = density(ref_distance.subclass[[class]][[1]])
    dens2 = density(ref_distance.subclass[[class]][[2]])
    dist.thresholds = as.vector(c( quantile(ref_distance.subclass[[class]][[1]], 0.75), 
                                   quantile(ref_distance.subclass[[class]][[2]], 0.01)))
    
    
    pass.dist = res$Distance[res$Cancer==class & res$PASS=="Y"]
    lines = lapply(pass.dist, function(v){
      line <- list(type = "line", line = list(color = paired.color[4], dash="dot", opacity =0.3), 
                   xref = "x", yref = "y",text = class);
      line[c('x0','x1')] = v; line[c('y0','y1')] = c(0,0.15);
      line; })
    p = plotly::plot_ly(x = dens1$x, y = dens1$y, type = 'scatter', mode = 'lines', name = class, 
                        line = list(color = paired.color[2]),
                        height = 300, width = 300*length(pass.subclass)) %>%
      plotly::add_trace(x = dens2$x, y = dens2$y, name = 'others', mode = 'lines', line = list(color = paired.color[8])) %>%
#threshold 1
      plotly::add_trace(x = rep(dist.thresholds[1], 2), y = c(0, 0.12), text = c(NA, paste0("75% quantile in ", class)),
                         mode = 'lines',
                        textposition = "top right", textfont=list( color=paired.color[2]),
                        line = list(color = paired.color[2], dash="dash", opacity = 0.3)) %>%
#threshold 2
      plotly::add_trace(x = rep(dist.thresholds[2], 2), y = c(0, 0.10), 
                        text = c(NA,paste0("1% quantile of ", class," to other cancers")),
                        textposition = "bottom right", textfont=list( color=paired.color[8]),
                         mode = 'lines',
                        line = list(color = paired.color[8], dash="dash", opacity = 0.3, size = 5)) %>%
      plotly::layout(shapes = lines, showlegend = FALSE,
                     yaxis = list(range = c(0, 0.15)),
                     annotations = list(x = 0.5 , y = 1.0, text = class, showarrow = F, xref='paper', yref='paper'))
    p
  })
  #frequency table
  t = data.frame(table(idx.label))
  colnames(t) = c("Cancer", "Count")
  t2 = table(res$Cancer[which(res$PASS=="Y")])
  t$`N PASS` = as.vector(t2[match(t$Cancer, names(t2))])
  t$`N PASS`[which(is.na(t$`N PASS`))] = 0
  t = t[match(unique(idx.label), t$Cancer),] #keep order, closer distance - larger distance
  return(list(res, t, figs))  
  
}
