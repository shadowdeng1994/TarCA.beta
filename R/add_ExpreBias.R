#' Title
#'
#' @param ExTree A ExTree file obtained from conv_ExTree.
#'
#' @return Return
#'
#' @export
#'
#' @examples
#' library("TarCA")
#' Example


add_ExpreBias <- function(ExTree){
    message("===> Adding ExpreBias")
    
    tmp.res <- 
    ExTree$Name2Meta %>% filter(!isTip) %>% rownames %>% 
    lapply(function(nnn){
        tmp.background <- sum(Ann$TipAnn)/nrow(Ann)
        tmp.modiann <- ExTree$AllDescendants[[nnn]] %>% mutate(TipAnn=as.logical(as.character(TipAnn)))
        tmp.size <- nrow(tmp.modiann)
        tmp.rate <- sum(tmp.modiann$TipAnn)/tmp.size
        
        data.frame(
            NodeName=nnn,
            CladeSize=tmp.size,
            Ratio=tmp.rate,
            p_chiq0=c(tmp.rate*tmp.size,(1-tmp.rate)*tmp.size) %>% chisq.test(.,p=c(tmp.background,1-tmp.background)) %>% .$p.value,
            OverBackground=tmp.rate>tmp.background
        )
    }) %>% bind_rows %>% 
    mutate(padj_chiq0=p.adjust(p_chiq0,"fdr")) %>% 
    mutate(BiasParent=OverBackground & padj_chiq0<0.05 & CladeSize>1) %>% 
    column_to_rownames("NodeName")
    
    tmp.out <- ExTree
    tmp.out[["ExpreBias"]] <- tmp.res
    return(tmp.out)
}