#' Checking the input files for TCA. 
#'
#' @param Tree A tree file of class "phylo" with node labels (using node ID by default).
#' @param Ann A dataframe with TipLabel and TipAnn (tip labels on the tree file and corresponding cell annotations).
#'
#' @export
#'

check_TCA_input <- function(Tree,Ann){
    message("===> Checking input files.")
    if(is.null(Tree$node.label)){ message("**Error** Lacking node labels in the input phylogeny.");return("Break") }
    if(length(unique(Tree$node.label))!=Tree$Nnode){ message("**Error** Node labels in the input phylogeny are NOT unique.");return("Break") }
    if(is.na(match("TipLabel",colnames(Ann)))|is.na(match("TipAnn",colnames(Ann)))){ message("**Error** Both columns TipLabel and TipAnn are required.");return("Break") }
    if(any(is.na(match(Tree$tip.label,Ann$TipLabel)))){ message("**Error** Incongruous tip labels in the input phylogeny and annotation.");return("Break") }
    return("AllPass")
}