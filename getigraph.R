
library(igraph)
library(gplots)

## please note, this .RData is just an example of elastic net output, you may replace it with the correct output
load(file = '/export/home/xurren/WTCProject/Results/PCLPHQ9LRSgenes.RData')


## PCLgenes is a vector contains non-zero gene names of PCL. The same with PHQ9genes and LRSgenes
getigraph = function(PCLgenes, PHQ9genes, LRSgenes){
  phenotypes = c(rep("PCL", length(PCLgenes)), rep("PHQ9", length(PHQ9genes)), rep("LRS", length(LRSgenes)))
  nonzerogenes = c(PCLgenes, PHQ9genes, LRSgenes)
  el = cbind(phenotypes, nonzerogenes)
  graph = graph_from_edgelist(el, directed = F)
  
  myvenn = venn(list(PCL = PCLgenes,PHQ9 = PHQ9genes,LRS = LRSgenes), show.plot=F)
  intersect.list = attr(myvenn, "intersections")
  intersect.names = names(attr(myvenn, "intersections"))
  mycolors = c("green", "red", "yellow", "skyblue", "pink", "brown", "cyan")
  
  ## this is to set vertex sizes
  vsizes = rep(4, length(vertex_attr(graph)$name))
  vsizes[vertex_attr(graph)$name %in% c("PCL", "PHQ9", "LRS")] = 10
  
  ## this is to set vertex colors
  vcolors = rep(NA, length(vertex_attr(graph)$name))
  for(i in 1:length(intersect.names)){
    vcolors[vertex_attr(graph)$name %in% intersect.list[[i]] ] = mycolors[i]
  }
  
  vcolors[vertex_attr(graph)$name == "PCL"] = mycolors[which(intersect.names == "PCL")]     
  vcolors[vertex_attr(graph)$name == "PHQ9"] = mycolors[which(intersect.names == "PHQ9")]   
  vcolors[vertex_attr(graph)$name == "LRS"] = mycolors[which(intersect.names == "LRS")]   
  
  V(graph)$color = vcolors
  plot(graph, vertex.size = vsizes, vertex.color = vcolors, vertex.label.cex = 0.5)
  legend("topleft", legend=intersect.names, col=mycolors, pch=21, pt.bg = mycolors)
  
}

## Here is an example:
getigraph(PCLgenes, PHQ9genes, LRSgenes)



getigraph2 = function(PCLHgenes, PCLAgenes, PCLRgenes, PCLNgenes){
  phenotypes = c(rep("PCL_H", length(PCLHgenes)), rep("PCL_A", length(PCLAgenes)), rep("PCL_R", length(PCLRgenes)), rep("PCL_N", length(PCLNgenes)))
  nonzerogenes = c(PCLHgenes, PCLAgenes, PCLRgenes, PCLNgenes)
  el = cbind(phenotypes, nonzerogenes)
  graph = graph_from_edgelist(el, directed = F)
  
  myvenn = venn(list(PCL_H = PCLHgenes, PCL_A = PCLAgenes, PCL_R = PCLRgenes, PCL_N = PCLNgenes), show.plot=F)
  intersect.list = attr(myvenn, "intersections")
  intersect.names = names(attr(myvenn, "intersections"))
  mycolors = sample(colors(), length(intersect.names))
  
  ## this is to set vertex sizes
  vsizes = rep(4, length(vertex_attr(graph)$name))
  vsizes[vertex_attr(graph)$name %in% c("PCL_H", "PCL_A", "PCL_R", "PCL_N")] = 10
  
  ## this is to set vertex colors
  vcolors = rep(NA, length(vertex_attr(graph)$name))
  for(i in 1:length(intersect.names)){
    vcolors[vertex_attr(graph)$name %in% intersect.list[[i]] ] = mycolors[i]
  }
  
  vcolors[vertex_attr(graph)$name == "PCL_H"] = mycolors[which(intersect.names == "PCL_H")]     
  vcolors[vertex_attr(graph)$name == "PCL_A"] = mycolors[which(intersect.names == "PCL_A")]   
  vcolors[vertex_attr(graph)$name == "PCL_R"] = mycolors[which(intersect.names == "PCL_R")] 
  vcolors[vertex_attr(graph)$name == "PCL_N"] = mycolors[which(intersect.names == "PCL_N")] 
  
  V(graph)$color = vcolors
  plot(graph, vertex.size = vsizes, vertex.color = vcolors, vertex.label.cex = 0.5)
  legend("topleft", legend=intersect.names, col=mycolors, pch=21, pt.bg = mycolors)
  
}






