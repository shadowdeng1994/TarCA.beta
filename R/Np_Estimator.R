#' Estimate Np using cell phylogeny for each sub-population.
#'
#' @param Tree A tree file of class "phylo" with node labels (using node ID by default).
#' @param Ann A dataframe with TipLabel and TipAnn (tip labels on the tree file and corresponding cell annotations).
#' @param ExTree (optional) A ExTree file obtained from conv_ExTree.
#' @param ReturnExTree Output in ExTree format (default FALSE).
#'
#' @return A dataframe storing the inferred Np for each cell type. 
#'
#' @export
#'
#' @examples
#' library("TarCA.beta")
#' load(system.file("Exemplar","Exemplar_TCA.RData",package = "TarCA"))
#' tmp.res <- Np_Estimator(
#' Tree = ExemplarData_1$Tree,
#' Ann = ExemplarData_1$Ann
#' )

Np_Estimator <- function(Tree,Ann,ExTree=NULL,ReturnExTree=FALSE){
    if(is.null(ExTree)){ ExTree <- conv_ExTree(Tree,Ann) }
    if(is.null(ExTree[["AllDescendants"]])){ ExTree <- add_AllDescendants(ExTree) }
    if(is.null(ExTree[["MonoClades"]])){ ExTree <- add_MonoClades(ExTree) }
    
    message("===> Estimating Np.")
    
    tmp.CladeSizeTable <- 
    ExTree[["MonoClades"]] %>% 
    group_by(TipAnn,CladeSize) %>% summarise(Freq=n()) %>% 
    arrange(TipAnn,CladeSize)
    
    tmp.res <- 
    tmp.CladeSizeTable %>% 
    group_by %>% mutate(LLL=paste0(CladeSize," (",Freq,")")) %>% 
    group_by(TipAnn) %>% summarise(MonoInfo=paste(LLL,collapse = ", "),Total=sum(CladeSize*Freq),Np=Total*(Total-1)/sum(CladeSize*(CladeSize-1)*Freq))
    
    if(ReturnExTree){
        tmp.out <- ExTree
        tmp.out[["Np_Estimator"]] <- tmp.res
        return(tmp.out)
    }else{
        return(tmp.res)
    }
}