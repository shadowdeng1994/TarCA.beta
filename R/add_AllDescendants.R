#' Function for adding AllDescendants field to a ExTree file.
#'
#' @param ExTree A ExTree file obtained from conv_ExTree.
#'
#' @return Return a ExTree file.
#'
#' @export
#'


add_AllDescendants <- function(ExTree){
    message("===> Adding AllDescendants.")
    
    subfun.get_all_descendants <- function(ExTree,NNName){
        if(!ExTree$Name2Meta[NNName,"isTip"]){
            # keep going
            tmp.NextLayer <- 
            ExTree$Name2Daughter[[NNName]] %>% 
            lapply(function(nnn){ subfun.get_all_descendants(ExTree,nnn) })
            tmp.NextLayer %>% bind_rows %>% rbind(mutate(.,NodeName=NNName) %>% unique) %>% remove_rownames
        }else{
            # return
            data.frame(NodeName=NNName,TipLabel=NNName,TipAnn=ExTree$Name2Meta[NNName,"TipAnn"])
        }
    }
    
    tmp.res <- 
    ExTree$Name2Meta %>% filter(isRoot) %>% rownames %>% 
    subfun.get_all_descendants(ExTree,.) %>% 
    split(f=.$NodeName,x=select(.,-NodeName))

    tmp.out <- ExTree
    tmp.out[["AllDescendants"]] <- tmp.res
    return(tmp.out)
}