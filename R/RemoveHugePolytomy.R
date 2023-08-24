#' Remove huge polytomies from the phylogeny
#' @description Remove confusing tips from huge polytomy.
#'
#' @param TTTree A tree file of class "phylo".
#' @param NNNum The maximum size of polytomy node.
#'
#' @return The processed tree.
#' @export
#'
fun.RemoveHugePolytomies <- function(TTTree,NNNum){
  tmp.treedata <- fun.GetTreeData(TTTree)

  tmp.ConfusedParent <-
    tmp.treedata %>%
    filter(isTip) %>%
    group_by(parent) %>% summarise(Count=n()) %>%
    mutate(Rate=Count/sum(Count)) %>%
    arrange(-Count) %>%
    filter(Count>NNNum)

  if(nrow(tmp.ConfusedParent)>0){
    tmp.treedata %>%
      filter(isTip) %>%
      right_join(tmp.ConfusedParent) %>%
      select(label) %>%
      mutate(DDD=1) %>%
      right_join(TTTree$tip.label %>% tbl_df %>% rename(label=value)) %>%
      filter(is.na(DDD)) %>%
      .$label %>%
      castor::get_subtree_with_tips(TTTree,.) %>% .$subtree
  }else{
    return(TTTree)
  }
}
