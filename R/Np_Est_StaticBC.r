#' Estimate the relative Np for static barcoding system.
#'
#' @param BCTable A dataframe with TipLabel and StaticBC.
#' @param Ann A dataframe with TipLabel and TipAnn (tip labels on the tree file and corresponding cell annotations).
#'
#' @return A dataframe storing the inferred Np for each cell type.
#'
#' @export
#'
#' @examples
#' BCTable <- data.frame(TipLabel=c("Cell_1","Cell_2","Cell_3"), StaticBC=c("BC1","BC1","BC2"))
#' Ann <- data.frame(TipLabel=c("Cell_1","Cell_2","Cell_3"), TipAnn=c("O1","O1","O2"))
#' Np_Est_StaticBC(BCTable, Ann)
#'

Np_Est_StaticBC <- function(BCTable,Ann){
    message("===> Estimating Np.")

    tmp.CladeSizeTable <-
    Ann %>%
    left_join(BCTable) %>%
    group_by(TipAnn,StaticBC) %>% summarise(CladeSize=n()) %>%
    group_by(TipAnn,CladeSize) %>% summarise(Freq=n()) %>%
    arrange(TipAnn,CladeSize)

    tmp.res <-
    tmp.CladeSizeTable %>%
    group_by %>% dplyr::mutate(LLL=paste0(CladeSize," (",Freq,")")) %>%
    group_by(TipAnn) %>% summarise(MonoInfo=paste(LLL,collapse = ", "),Total=sum(CladeSize*Freq),Np=Total*(Total-1)/sum(CladeSize*(CladeSize-1)*Freq))

    return(tmp.res)
}
