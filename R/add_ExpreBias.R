#' Function for adding ExpreBias field during LEU estimation.
#'
#' @param ExTree A ExTree file obtained from conv_ExTree.
#'
#' @return Return a ExTree file.
#'
#' @export
#'

add_ExpreBias <- function(ExTree){
    message("===> Adding ExpreBias")

    ## Single bottom-up pass gives, for every node, the number of positive
    ## (TRUE) descendants and the clade size; the chi-square test is then
    ## evaluated vectorised instead of node by node.
    ExTree <- fun.EnsureNodeIndex(ExTree)
    ni <- ExTree[["NodeIndex"]]

    tipAnnRaw <- ExTree$Name2Meta$TipAnn[match(ni$label[ni$isTip],rownames(ExTree$Name2Meta))]
    pos <- fun.NodeSum(ExTree,as.numeric(as.logical(as.character(tipAnnRaw))))

    internalLabels <- rownames(ExTree$Name2Meta)[!ExTree$Name2Meta$isTip]
    ii <- match(internalLabels,ni$label)

    CladeSize <- ni$size[ii]
    tmp.background <- sum(ExTree$Ann$TipAnn)/nrow(ExTree$Ann)
    Ratio <- pos[ii]/CladeSize
    p_chiq0 <- pchisq((pos[ii]-CladeSize*tmp.background)^2/(CladeSize*tmp.background*(1-tmp.background)),df=1,lower.tail=FALSE)
    OverBackground <- Ratio>tmp.background
    padj_chiq0 <- p.adjust(p_chiq0,"fdr")
    BiasParent <- OverBackground & padj_chiq0<0.05 & CladeSize>1

    tmp.res <- data.frame(
        CladeSize=CladeSize,
        Ratio=Ratio,
        p_chiq0=p_chiq0,
        OverBackground=OverBackground,
        padj_chiq0=padj_chiq0,
        BiasParent=BiasParent,
        row.names=internalLabels
    )

    tmp.out <- ExTree
    tmp.out[["ExpreBias"]] <- tmp.res
    return(tmp.out)
}
