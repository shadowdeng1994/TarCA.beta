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


add_MonoClades <- function(ExTree){
    message("===> Adding MonoClades.")
    
    subfun.get_monophyletic_clade <- function(ExTree,NNNodeName){
        tmp.AnnVector <- ExTree$AllDescendants[[NNNodeName]]$TipAnn
        if(length(unique(tmp.AnnVector))>1){
            # keep going
            tmp.NextLayer <- 
            ExTree$Name2Daughter[[NNNodeName]] %>% 
            lapply(function(nnn){ subfun.get_monophyletic_clade(ExTree,nnn) })
            tmp.NextLayer %>% bind_rows
        }else{
            # return
            data.frame(NodeName=NNNodeName,TipAnn=unique(tmp.AnnVector),CladeSize=length(tmp.AnnVector))
        }
    }
    
    tmp.res <- 
    ExTree$Name2Meta %>% filter(isRoot) %>% rownames %>% 
    subfun.get_monophyletic_clade(ExTree,.) %>% 
    arrange(TipAnn,CladeSize)
    
    tmp.out <- ExTree
    tmp.out[["MonoClades"]] <- tmp.res
    return(tmp.out)
}