#' Function for adding FilterBiasParent field to a ExTree file.
#'
#' @param ExTree A ExTree file obtained from conv_ExTree.
#'
#' @return Return a ExTree file.
#'
#' @export
#'

add_FilterBiasParent <- function(ExTree){
    ## this function is used for filtering the bias node on different hierarchy.
    ##      0.5-------(x)
    ##       |
    ##      0.8-------(v)
    ##       |
    ##      0.8-------(x)
    
    message("===> Adding FilterBiasParent")
    
    ExTree <- fun.EnsureNodeIndex(ExTree)
    ni <- ExTree[["NodeIndex"]]
    ## In a tree two clades overlap iff one contains the other; the DFS
    ## interval [tin, tout] makes that test O(1) instead of a tip-set intersect.
    subfun.isP2D <- function(Node_1,Node_2){
        i1 <- match(Node_1,ni$label)
        i2 <- match(Node_2,ni$label)
        ni$tin[i1]<=ni$tin[i2] && ni$tout[i2]<=ni$tout[i1]
    }
    
    tmp.biasparent <- ExTree$ExpreBias %>% filter(BiasParent) %>% arrange(-CladeSize)
    
    tmp.res <- tmp.biasparent %>% mutate(BeatByChild=NA,BeatByParent=NA,Filter=TRUE)
    if(nrow(tmp.biasparent)>1){
        tmp.unclearP2D <- 
        tmp.biasparent %>% rownames %>% list(.,.) %>% expand.grid %>% 
        mutate(Var1=Var1 %>% as.character,Var2=Var2 %>% as.character) %>% 
        mutate(N1=tmp.biasparent[Var1,"CladeSize"],N2=tmp.biasparent[Var2,"CladeSize"]) %>% 
        filter(N1>N2)
        if(nrow(tmp.unclearP2D)>0){
            tmp.isP2D <- 
            tmp.unclearP2D %>% 
            group_by(Var1,Var2) %>% mutate(isP2D=subfun.isP2D(Var1,Var2)) %>% 
            group_by %>% mutate(ParentBias=ExTree$ExpreBias[Var1,"Ratio"],ChildBias=ExTree$ExpreBias[Var2,"Ratio"]) %>% 
            filter(isP2D)
            if(nrow(tmp.isP2D)>0){
                tmp.res <- tmp.biasparent %>% rownames_to_column("NodeName") %>% 
                left_join(tmp.isP2D %>% group_by(NodeName=Var1) %>% summarise(BeatByChild=any(ChildBias>ParentBias)),by="NodeName") %>% 
                left_join(tmp.isP2D %>% group_by(NodeName=Var2) %>% summarise(BeatByParent=any(ChildBias<=ParentBias)),by="NodeName") %>% 
                mutate(BeatByChild=replace_na(BeatByChild,FALSE),BeatByParent=replace_na(BeatByParent,FALSE)) %>% 
                mutate(Filter=!BeatByChild & !BeatByParent) %>% 
                column_to_rownames("NodeName")
            }
        }
    }
    
    tmp.out <- ExTree
    tmp.out[["FilterBiasParent"]] <- tmp.res
    return(tmp.out)
}