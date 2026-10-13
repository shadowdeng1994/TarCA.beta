#' Build a lightweight node index for an ExTree (internal).
#'
#' @param Tree A tree file of class "phylo".
#'
#' @return A list with per-node label, parent, children, traversal order and
#'   subtree sizes / intervals, all indexed by node number (1..Ntip+Nnode).
#'
#' @noRd
#'
fun.NodeIndex <- function(Tree){
    ntip <- length(Tree$tip.label)
    nnode <- Tree$Nnode
    N <- ntip + nnode
    if(is.null(Tree$node.label)){ Tree$node.label <- paste0("Node_",seq_len(nnode)) }
    label <- c(Tree$tip.label,Tree$node.label)

    edge <- Tree$edge
    parent <- integer(N)
    parent[edge[,2L]] <- edge[,1L]
    root <- which(parent==0L)
    if(length(root)!=1L){ root <- edge[1L,1L] }
    parent[root] <- root

    kids <- split(edge[,2L],edge[,1L])
    children <- vector("list",N)
    children[as.integer(names(kids))] <- kids

    isTip <- seq_len(N)<=ntip

    ## iterative depth-first pre-order (parent before children, subtrees contiguous)
    stack <- integer(N)
    top <- 1L
    stack[1L] <- root
    preorder <- integer(N)
    idx <- 0L
    while(top>0L){
        v <- stack[top]
        top <- top-1L
        idx <- idx+1L
        preorder[idx] <- v
        kv <- children[[v]]
        nk <- length(kv)
        if(nk>0L){
            stack[(top+1L):(top+nk)] <- kv
            top <- top+nk
        }
    }

    ## reverse pre-order == children before parents, sufficient for bottom-up accumulation
    order <- preorder[N:1L]
    tin <- integer(N)
    tin[preorder] <- seq_len(N)

    size <- integer(N)
    size[isTip] <- 1L
    tout <- tin
    for(v in order){
        if(!isTip[v]){
            kv <- children[[v]]
            size[v] <- sum(size[kv])
            tout[v] <- max(tout[kv])
        }
    }

    list(label=label,isTip=isTip,parent=parent,children=children,root=root,
         preorder=preorder,order=order,tin=tin,tout=tout,size=size,
         ntip=ntip,nnode=nnode,N=N)
}

#' Ensure an ExTree carries a node index (internal).
#'
#' @param ExTree A ExTree file obtained from conv_ExTree.
#'
#' @return A ExTree file with a `NodeIndex` field.
#'
#' @noRd
#'
fun.EnsureNodeIndex <- function(ExTree){
    if(is.null(ExTree[["NodeIndex"]])){
        ExTree[["NodeIndex"]] <- fun.NodeIndex(ExTree[["Tree"]])
    }
    ExTree
}

#' Aggregate tip-level values up the tree (internal).
#'
#' @param ExTree A ExTree file carrying a `NodeIndex` field.
#' @param tipValues A numeric vector, one value per tip (in `Tree$tip.label`
#'   order).
#'
#' @return A numeric vector of subtree sums, indexed by node number.
#'
#' @noRd
#'
fun.NodeSum <- function(ExTree,tipValues){
    ni <- ExTree[["NodeIndex"]]
    isTip <- ni$isTip
    children <- ni$children
    out <- numeric(ni$N)
    out[isTip] <- tipValues
    for(v in ni$order){
        if(!isTip[v]){ out[v] <- sum(out[children[[v]]]) }
    }
    out
}
