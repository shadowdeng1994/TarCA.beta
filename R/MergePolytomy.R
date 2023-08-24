#' Merge short edge into polytomy
#'
#' @param TTTree A tree file of class "phylo".
#' @param EEEpsilon Non-negative numeric, specifying the maximum edge length for an edge to be considered ``short'' and thus to be eliminated. Typically 0 or some small positive number.
#'
#' @return A list containing processed tree and distribution of the edge length.
#' @export
#'
#' @examples
#' library(TarCA)
#' set.seed(123)
#' tmp.tree <- rtree(100)
#' fun.MergePolytomy(tmp.tree,EEEpsilon = 0.1)
#'

fun.MergePolytomy <- function(TTTree,EEEpsilon=0.00001){
  tmp.out <- list()
  tmp.out[["tree"]] <- castor::merge_short_edges(TTTree,edge_length_epsilon = EEEpsilon) %>% .$tree
  tmp.out[["plot"]] <- fun.PlotEdgeLengthDistribution(TTTree)
  return(tmp.out)
}


#' Plot the distribution of edge length
#'
#' @param TTTree A tree file of class "phylo".
#'
fun.PlotEdgeLengthDistribution <- function(TTTree){
  data.frame("EdgeLength"=TTTree$edge.length) %>%
    ggplot2::ggplot(ggplot2::aes(x=EdgeLength))+
    ggplot2::geom_histogram()+
    ggplot2::scale_x_continuous(trans="log10")+
    ggplot2::theme_bw()
}
