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
#' library("TarCA")
#' load(system.file("Exemplar","Exemplar_TCA.RData",package = "TarCA"))
#' tmp.ExTree <- conv_ExTree(Tree = ExemplarData_1$Tree,Ann = ExemplarData_1$Ann)

conv_ExTree <- function(Tree,Ann,ForceFactor=TRUE){
    if(check_TCA_input(Tree,Ann)!="AllPass"){ return() }
    
    message("===> Converting to ExTree.")
    
    TreeData <- Tree %>% ggtree %>% .$data
    
    Name2Meta <- TreeData %>% 
    mutate(ParentLabel=.$label[parent],isRoot=parent==node) %>% 
    rename(TipLabel=label) %>% left_join(Ann,by="TipLabel") %>% column_to_rownames("TipLabel")
    if(ForceFactor & !is.factor(Name2Meta$TipAnn)){ Name2Meta <- Name2Meta %>% mutate(TipAnn=factor(TipAnn)) }
    
    Name2Daughter <- Name2Meta %>% filter(!isRoot) %>% split(f=.$ParentLabel,x=rownames(.))
    
    list(
        "Tree"=Tree,
        "Ann"=Ann,
        "NameList"=Name2Meta %>% rownames,
        "AnnList"=Name2Meta$TipAnn %>% unique %>% na.omit %>% sort,
        "Dim"=c(length(levels(Name2Meta$TipAnn)),nrow(Name2Meta)),
        "Name2Meta"=Name2Meta,
        "Name2Daughter"=Name2Daughter
    )
}