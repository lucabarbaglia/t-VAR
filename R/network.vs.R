#' network.vs
#' 
#' Function to obtain network of volatility spillovers
#' 
#' @param spill matrix of spillover effects. Dimension JxJ, where J is the number of time series.
#' @param main.label network names.  Default NULL.
#' @param names node names. Default is the column names of beta.
#' @param ww weight of the thickness of the arrows
#' @param colour arrows' colour. Default "#619CFF".
#' 
#' @return a network of spillover effects consisting of J nodes
#' 
#' @export



network.vs<-function(spill, main.label=NULL, names=NULL,  ww=5, colour="#619CFF"){
  
  # Dimension
  J<-ncol(spill)
  
  # Nodes labels
  if (!is.null(names)){
    colnames(spill)<-rownames(spill)<-names
  }
  
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  lab.locs <- radian.rescale(x=1:J, direction=-1, start=0) 
  
  GRAPH<-igraph::graph.adjacency(adjmatrix=t(spill), mode="directed", diag=FALSE, weighted=T)
  
  igraph::E(GRAPH)$width <- abs(igraph::E(GRAPH)$weight)
  igraph::E(GRAPH)$color<-colour
  
  loc<-igraph::layout.circle(GRAPH) # location of the nodes
  
  plot(margin=rep(10^-10,4), 
       GRAPH, 
       layout=igraph::layout.circle(GRAPH), 
       vertex.label.dist=4, 
       vertex.color="white", 
       vertex.frame.color="white",
       vertex.shape="circle", 
       vertex.label.color="black",
       vertex.label.degree=lab.locs,
       edge.curved=seq(-0.3, 0.3, length = igraph::ecount(GRAPH)),
       edge.arrow.width=1.3
  )
  
  if(!is.null(main.label)){
    title(main=main.label, cex.main=1)
  }
  
  
} # end function

