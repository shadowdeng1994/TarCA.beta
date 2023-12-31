#' Convert a tree file and a annotation dataframe to ExTree format. 
#'
#' @param Tree A tree file of class "phylo" with node labels (using node ID by default).
#' @param Ann A dataframe with TipLabel and TipAnn (tip labels on the tree file and corresponding cell annotations).
#' @param ForceFactor Force transform TipAnn into a factor (default TRUE).
#'
#' @return Return a ExTree file.
#'
#' @export
#'
#' @examples
#' library("TarCA.beta")
#' load(system.file("Exemplar","Exemplar_TCA.RData",package = "TarCA"))
#' tmp.ExTree <- conv_ExTree(Tree = ExemplarData_1$Tree,Ann = ExemplarData_1$Ann)

conv_ExTree <- function(Tree,Ann=NULL,ReplaceNode=FALSE,ForceFactor=TRUE){
    tmp.tree <- Tree
    tmp.ann <- Ann

    if(is.null(tmp.ann)){ tmp.ann <- tmp.tree$tip.label %>% data.frame(TipLabel=.,TipAnn=.) }
    if(is.null(tmp.tree$node.label)|ReplaceNode){ tmp.tree$node.label <- paste0("Node_",1:tmp.tree$Nnode) }

    if(check_TCA_input(tmp.tree,tmp.ann)!="AllPass"){ return() }
    
    message("===> Converting to ExTree.")
    
    TreeData <- tmp.tree %>% ggtree %>% .$data
    
    Name2Meta <- TreeData %>% 
    dplyr::mutate(ParentLabel=.$label[parent],isRoot=parent==node) %>% 
    dplyr::rename(TipLabel=label) %>% dplyr::left_join(tmp.ann,by="TipLabel") %>% tibble::column_to_rownames("TipLabel")
    if(ForceFactor & !is.factor(Name2Meta$TipAnn)){ Name2Meta <- Name2Meta %>% dplyr::mutate(TipAnn=factor(TipAnn)) }
    
    Name2Daughter <- Name2Meta %>% dplyr::filter(!isRoot) %>% split(f=.$ParentLabel,x=rownames(.))
    
    list(
        "Tree"=tmp.tree,
        "Ann"=tmp.ann,
        "NameList"=Name2Meta %>% rownames,
        "AnnList"=Name2Meta$TipAnn %>% unique %>% na.omit %>% sort,
        "Dim"=c(length(levels(Name2Meta$TipAnn)),nrow(Name2Meta)),
        "Name2Meta"=Name2Meta,
        "Name2Daughter"=Name2Daughter
    )
}