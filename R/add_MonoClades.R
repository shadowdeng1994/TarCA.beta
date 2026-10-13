#' Function for adding MonoClades field to a ExTree file.
#'
#' @param ExTree A ExTree file obtained from conv_ExTree.
#' @param Sensitivity A numeric in (0,1] giving the minimal fraction of the
#'   dominant annotation required for a clade to be called monophyletic. With
#'   the default 1 every tip in the clade must share one annotation; with e.g.
#'   0.9 a clade whose dominant annotation covers at least 90% of its tips is
#'   called monophyletic. Maximal monophyletic clades are returned, and each
#'   clade contributes the number of tips of its dominant annotation as
#'   `CladeSize`.
#'
#' @return Return a ExTree file.
#'
#' @export
#'

add_MonoClades <- function(ExTree,Sensitivity=1){
    message("===> Adding MonoClades.")

    ExTree <- fun.EnsureNodeIndex(ExTree)
    ni <- ExTree[["NodeIndex"]]
    N <- ni$N
    isTip <- ni$isTip
    children <- ni$children

    tipAnn <- ExTree$Name2Meta$TipAnn[match(ni$label[isTip],rownames(ExTree$Name2Meta))]
    if(is.factor(tipAnn)){
        lv <- levels(tipAnn)
        code <- as.integer(tipAnn)
        asFactor <- TRUE
    }else{
        lv <- sort(unique(tipAnn))
        code <- match(tipAnn,lv)
        asFactor <- FALSE
    }

    if(Sensitivity==1){
        ## O(N) streaming: a node is monophyletic iff all descendant tips share
        ## one annotation.
        all_same <- logical(N)
        common <- integer(N)
        all_same[isTip] <- TRUE
        common[isTip] <- code
        for(v in ni$order){
            if(!isTip[v]){
                kv <- children[[v]]
                if(length(kv)>0L && all(all_same[kv])){
                    c1 <- common[kv[1L]]
                    if(all(common[kv]==c1)){ all_same[v] <- TRUE; common[v] <- c1 }
                }
            }
        }
        mono <- all_same
        domType <- common
        domCount <- ni$size
    }else{
        ## Count matrix (N x K) accumulated bottom-up; a node is monophyletic
        ## when its dominant annotation reaches the requested fraction.
        K <- length(lv)
        Count <- matrix(0L,N,K)
        Count[cbind(which(isTip),code)] <- 1L
        for(v in ni$order){
            if(!isTip[v]){ Count[v,] <- colSums(Count[children[[v]],,drop=FALSE]) }
        }
        total <- rowSums(Count)
        domType <- max.col(Count,ties.method="first")
        domCount <- Count[cbind(seq_len(N),domType)]
        mono <- total>0L & (domCount/total)>=Sensitivity
    }

    ## keep maximal monophyletic clades: monophyletic node whose parent is not
    sel <- mono & (!mono[ni$parent] | ni$parent==seq_len(N))
    idx <- which(sel)

    if(asFactor){
        tipAnnOut <- factor(lv[domType[idx]],levels=lv)
    }else{
        tipAnnOut <- lv[domType[idx]]
    }

    tmp.res <- 
    data.frame(NodeName=ni$label[idx],TipAnn=tipAnnOut,CladeSize=domCount[idx]) %>% 
    arrange(TipAnn,CladeSize)
    attr(tmp.res,"Sensitivity") <- Sensitivity

    tmp.out <- ExTree
    tmp.out[["MonoClades"]] <- tmp.res
    return(tmp.out)
}
