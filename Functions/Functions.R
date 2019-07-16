
#' multintersect
#'
#' Performs pair-wise overlaps between the elements of one or two lists
#'
#' @param ll A list of vectors to be compared
#' @param ll2 An optional second list of vectors
#' @param universe An optional vector of the universe (all terms), used to calculate overlap probabilities and enrichments. If NULL, the union of all lists will be used. Only elements of `ll` and `ll2` that are in the `universe` will be considered.
#' @param keyCol The values to use for the heatmap's colors. Possible values are: prob, enrichment, log2Enrichment, log10Prob, overlap, jaccard
#' @param keyWrite The values to be written in the cells. Possible values are: prob, enrichment, log2Enrichment, log10Prob, overlap, jaccard
#' @param addSetSize Logical; whether to add the set sizes to set names (default TRUE)
#' @param breakNames If not NULL (default), should be an number indicating the character length threshold above which set names should be split onto two lines.
#' @param margin Either a single number indicating the size of the left and bottom margins, or a list (l,r,b,t,pad) indicating the size of each margin.
#' @param returnTables Logical; whether to return tables in addition to the plot (default FALSE)
#'
#' @return Plots a heatmap and return a list of pair-wise overlap metrics
#'
#' @export
multintersect <- function(ll, ll2=NULL, universe=NULL, keyCol="log2Enrichment", keyWrite="overlap", addSetSize=TRUE, breakNames=NULL, margin=100, returnTables=FALSE){
  library(plotly)
  keyCol <- match.arg(keyCol, c("prob","enrichment","log2Enrichment","log10Prob","overlap","jaccard"))
  keyWrite <- match.arg(keyWrite, c("prob","enrichment","log2Enrichment","log10Prob","overlap","jaccard"))
  if(is.null(ll2)){
    symm <- TRUE
    ll2 <- ll
  }else{
    symm <- FALSE
  }
  if(is.null(universe)) universe <- unique(c(unlist(ll),unlist(ll2)))
  ll <- lapply(ll,y=universe,FUN=intersect)
  ll2 <- lapply(ll2,y=universe,FUN=intersect)
  ll <- ll[sapply(ll,FUN=length)>0]
  ll2 <- ll2[sapply(ll2,FUN=length)>0]
  if(addSetSize){
    names(ll) <- paste0(names(ll),"\n(",sapply(ll,FUN=length),")")
    names(ll2) <- paste0(names(ll2)," (",sapply(ll2,FUN=length),")")
  }
  if(!is.null(breakNames)){
    names(ll) <- breakStrings(names(ll),breakNames)
    names(ll2) <- breakStrings(names(ll2),breakNames)
  }
  n <- length(ll)
  j <- length(ll2)
  m <- matrix(0,nrow=n,ncol=j)
  colnames(m) <- names(ll2)
  rownames(m) <- names(ll)
  prob <- m
  enr <- m
  jacc <- m
  for(i in 1:n){
    for(j in 1:length(ll2)){
      m[i,j] <- length(intersect(ll[[i]],ll2[[j]]))
      if(names(ll)[i]==names(ll2)[j]){
        prob[i,j] <- NA
        enr[i,j] <- NA
        jacc[i,j] <- NA
      }else{
        getEnrichment <- function(sfor,sin,universe){
          sfor <- intersect(sfor,universe)
          sin <- intersect(sin, universe)
          expected <- length(sin)*(length(sfor)/length(universe))
          obs <- length(intersect(sfor,sin))
          return(obs/expected)
        }
        
        overlap.prob <- function(set1,set2,universe,lower=F){
          set1 <- as.character(set1)
          set2 <- as.character(set2)
          if(class(universe)=="character"){
            set1 <- intersect(set1,universe)
            set2 <- intersect(set2,universe)
            universe <- length(unique(universe))
          }
          set1 <- unique(set1)
          set2 <- unique(set2)
          ov <- sum(set1 %in% set2)
          phyper(max(0,ov-1), length(set1), universe-length(set1), length(set2), lower.tail=lower)
        }
        
        enr[i,j] <- getEnrichment(ll[[i]],ll2[[j]],universe)
        prob[i,j] <- overlap.prob(ll[[i]],ll2[[j]],universe,lower=enr[i,j]<1)
        jacc[i,j] <- m[i,j]/length(unique(c(ll[[i]],ll2[[j]])))
      }
    }
  }
  x <- switch(keyCol,
              prob=prob,
              enrichment=enr,
              log2Enrichment=log2(enr),
              log10Prob=-log10(prob),
              overlap=m,
              jaccard=jacc)
  y <- switch(keyWrite,
              prob=format(prob,digits=1),
              enrichment=format(enr,digits=2,trim=T,drop0trailing=T),
              log2Enrichment=round(log2(enr),1),
              log10Prob=round(-log10(prob),1),
              overlap=m,
              jaccard=round(jacc,3))
  y[is.infinite(y)] <- NA
  x[is.infinite(x) & x>0] <- max(x[which(!is.infinite(x))],na.rm=T)
  x[is.infinite(x) & x<0] <- min(x[which(!is.infinite(x))],na.rm=T)
  lab <- matrix(paste0( rep( gsub("\n"," ",row.names(x),fixed=T),ncol(x) ), "\n",
                        rep( gsub("\n"," ",colnames(x), fixed=T),each=nrow(x)), "\n",
                        "overlap: ", as.numeric(m), "\n",
                        as.character(format(enr,digits=2,trim=T,drop0trailing=T)),"-fold enrichment \n",
                        "p~",as.character(format(prob,digits=2))
  ),nrow=nrow(x),ncol=ncol(x))

  if(symm){
    for(i in 1:ncol(x)) lab[i,i] <- row.names(x)[i]
  }
  if(length(margin)==1){
    margin <- list(l=margin, r=5, b=margin, t=5, pad=4)
  }
  p <- plot_ly(x=colnames(x), y=row.names(x), z=x, type="heatmap", text=lab, hoverinfo = 'text') %>% colorbar(title = keyCol) %>%
    layout(xaxis = list(showgrid = FALSE), yaxis=list(showgrid = FALSE), margin = margin)
  p <- p %>% add_annotations(x = rep(colnames(x),each=nrow(x)), y = rep(row.names(x),ncol(x)), text = y, xref="x", yref="y", showarrow=FALSE)
  if(returnTables) return(list(overlaps=m, probability=prob, enrichment=enr, jaccard=jacc, plot=p))
  p
}

byheatmap <- function (x, scale = "none", ...) 
{
  require("pheatmap")
  scale <- match.arg(scale, c("none", "row", "column"))
  x <- x[which(apply(x, 1, FUN = function(y) {
    !all(is.na(y))
  })), which(apply(x, 2, FUN = function(y) {
    !all(is.na(y))
  }))]
  pheatmap(x, scale = scale, color = colorRampPalette(c("blue", 
                                                        "black", "yellow"))(29), border_color = NA, ...)
}