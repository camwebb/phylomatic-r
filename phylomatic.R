library(ape)

phylomatic <- function(phylo, taxa) {

    ## 1. Create missing objects needed for C input 
    
    ## $node.label?
    if(length(phylo$node.label) == 0) {
        phylo$node.label <- vector("character", length(phylo$edge[,1]))
        phylo$node.label[] <- ""
        NL <- 0
    } else {
        NL <- 1
    }

    ## $tip.label?
    if(length(phylo$tip.label) == 0) {
        phylo$tip.label <-
            vector("character", length(phylo$edge[,1]) - phylo$Nnode)
        phylo$tip.label[] <- ""
        TL <- 0
    } else {
        TL <- 1
    }

    ## $edge.length?
    if(length(phylo$edge.length) == 0) {
        phylo$edge.length <- vector("double", length(phylo$edge[,1]))
        phylo$edge.length[] <- NaN
        ## any missing BLs in $edge.length are NaN.  NA could be used,
        ## and somehow converted to NULL in C (via stdio.h) but C
        ## deals with NA as NaN anyway, and NULL would need to be
        ## converted back to NaN in ape
        BL <- 0
    } else {
        BL <- 1
    }
    
    ## $root.edge?
    if(length(phylo$root.edge) == 0) {
        phylo$root.edge <- vector("double", 1)
        phylo$root.edge[] <- NaN
        RL <- 0
    } else {
        RL <- 1
    }

    dyn.load("phylomatic.so")

    x <- .C("phylomatic",
       NAOK = TRUE,
       ## INPUT = APE phylo
       ## dimension of $edge
       ape_nedge    = as.integer(length(phylo$edge[,1])),
       ## vectors of edge lookup table: to, from
       ape_to       = as.integer(phylo$edge[,1]),
       ape_from     = as.integer(phylo$edge[,2]),
       ## branch lengths
       ape_bl       = as.double(phylo$edge.length),
       ## root edge
       ape_rbl      = as.double(phylo$root.edge),
       ## number of tips
       ape_ntip     = as.integer(length(phylo$tip.label)),
       ## tip labels
       ape_tip_lab  = as.character(phylo$tip.label),
       ## inner node labels
       ape_node_lab = as.character(phylo$node.label),
       ## max label string length
       ape_max_lab_len = 
           as.integer(max(nchar(c(phylo$tip.label,phylo$node.label)))),

       ## OUTPUT ~= APE phylo
       ## Dimension arrays as size of input megatree
       ## $Nnode
       out_nnode    = as.integer(vector("integer", 1)),
       ## Two vectors of $edge lookup table: [,1], [,2]
       out_to       = as.integer(vector("integer", length(phylo$edge[,1]))),
       out_from     = as.integer(vector("integer", length(phylo$edge[,1]))),
       ## Branch lengths: $edge.length
       out_bl       = as.double(vector("double", length(phylo$edge[,1]))),
       ## Root branch length: $root.edge
       out_rbl      = as.double(vector("double", 1)),
       ## $tip.label
       out_tip_lab  = as.character(vector("character", length(phylo$edge[,1]))),
       ## $node.label
       out_node_lab = as.character(vector("character", length(phylo$edge[,1]))),
       ## N tips
       out_ntip     = as.integer(vector("integer", 1))
       )
    ## NB you can reference output sub-objects by out[[12]] or $out_from

    dyn.unload("phylomatic.so")

    ## All phylo have minimally $Nnode and $edge
    ape <- list(Nnode = x$out_nnode,
                edge = cbind(x$out_to[  1:(x$out_nnode + x$out_ntip -1)],
                             x$out_from[1:(x$out_nnode + x$out_ntip -1)]))

    ## Check for additional objects. In future, pass from C, not from original
    if (NL == 1) {
        ape <- c(ape, list(node.label = x$out_node_lab[1:x$out_nnode]))
    }
    if (TL == 1) {
        ape <- c(ape, list(tip.label = x$out_tip_lab[1:x$out_ntip]))
    }
    if (BL == 1) {
        ape <- c(ape, list(edge.length = x$out_bl))
    }
    if (RL == 1) {
        ape <- c(ape, list(root.edge = x$out_rbl))
    }

    class(ape) <- "phylo"
    return(ape)
}




