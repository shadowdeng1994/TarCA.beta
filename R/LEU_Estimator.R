#' Detection of lineage specific expression upregulation.
#'
#' @param Tree A tree file of class "phylo" with node labels (using node ID by default).
#' @param Ann A dataframe with TipLabel and TipAnn (tip labels on the tree file and corresponding cell annotations).
#' @param ExTree (optional) A ExTree file obtained from conv_ExTree.
#'
#' @return A list containing Tree, Ann, BiasNode, FilterBiasNode, BiasFig, PureNode, PureNode2Organ and EffN. The estimated Np are store in EffN.
#'
#' @export
#'
#' @examples
#' library("TarCA")
#' load(system.file("Exemplar","Exemplar_LEU.RData",package = "TarCA"))
#' tmp.res <- LEU_Estimator(
#' Tree = ExemplarData_2$Tree,
#' Ann = ExemplarData_2$Ann,
#' Fileout = NULL,
#' ReturnNp = TRUE
#' )

LEU_Estimator <- function(Tree,Ann,ExTree=NULL){
    if(is.null(ExTree)){ ExTree <- conv_ExTree(Tree,Ann,ForceFactor = FALSE) }
    if(is.null(ExTree[["AllDescendants"]])){ ExTree <- add_AllDescendants(ExTree) }
    if(is.null(ExTree[["ExpreBias"]])){ ExTree <- add_ExpreBias(ExTree) }
    if(is.null(ExTree[["FilterBiasParent"]])){ ExTree <- add_FilterBiasParent(ExTree) }
    
    message("===> Estimating Np.")
    
    tmp.CladeSizeTable <- 
    ExTree[["FilterBiasParent"]] %>% filter(Filter) %>% rownames %>% 
    lapply(function(nnn){ data.frame(NodeName=nnn,CladeSize=ExTree[["AllDescendants"]][[nnn]] %>% .$TipAnn %>% sum) }) %>% 
    bind_rows %>% 
    group_by(CladeSize) %>% summarise(Freq=n()) %>% 
    rbind(data.frame(CladeSize=1,Freq=sum(Ann$TipAnn)-sum(.$CladeSize*.$Freq)),.) %>% 
    group_by %>% mutate(TipAnn=TRUE) %>% arrange(TipAnn,CladeSize)

    tmp.CladeSizeTable %>% 
    group_by %>% mutate(LLL=paste0(CladeSize," (",Freq,")")) %>% 
    group_by(TipAnn) %>% summarise(MonoInfo=paste(LLL,collapse = ", "),Total=sum(CladeSize*Freq),Np=Total*(Total-1)/sum(CladeSize*(CladeSize-1)*Freq))
}