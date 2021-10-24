
phylomatic <- function(phylo, taxa) {

    ## edges?
    if(length(phylo$edge.length) == 0) {
        phylo$edge.length <- vector("double", length(phylo$edge[,1]))
        ## phylo$edge.length[] <- NA
        ## not needed, C deals with NA as NaN anyway 
    }
    ## recode NaNs to NAs
    ## phylo$edge.length[is.nan(phylo$edge.length)] <- NA
    ## root edge?
    if(length(phylo$root.edge) == 0) {
        phylo$root.edge <- vector("double", 1)
        ## phylo$root.edge[] <- NA
    }
    ## ## number of nodes for FY
    ## fn     <- length(phylo$edge[,1]) + 1

    dyn.load("phylomatic.so")

    ## TODO : null bls passed into C as 0.0 - need to inform if no bLs
    
    out <- .C("phylomatic",
       NAOK = TRUE,
       ## INPUT (APE tree)
       ## 1. dimension of $edge
       as.integer(length(phylo$edge[,1])),
       ## 2,3. vectors of edge lookup table: to, from
       as.integer(phylo$edge[,1]),
       as.integer(phylo$edge[,2]),
       ## 4. branch lengths
       as.double(phylo$edge.length),
       ## 5. root edge
       as.double(phylo$root.edge),
       ## 6. number of tips 
       as.integer(length(phylo$tip.label)),
       ## 7. tip labels
       as.character(phylo$tip.label),
       ## 8. inner node labels
       as.character(phylo$node.label),
       ## 9. max label string length
       as.integer(max(nchar(c(phylo$tip.label,phylo$node.label)))),

       ## OUTPUT (APE tree)
       ## Dimension arrays as size of input megatree
       ## 10. $Nnode
       as.integer(vector("integer", 1)),
       ## 11, 12. Two vectors of $edge lookup table: [,1], [,2]
       as.integer(vector("integer", length(phylo$edge[,1]))),
       as.integer(vector("integer", length(phylo$edge[,1]))),
       ## 13. Branch lengths: $edge.length
       as.double(vector("double", length(phylo$edge[,1]))),
       ## 14. Root branch length: $root.edge
       as.double(vector("double", 1)),
       ## 15. $tip.label
       as.character(vector("character", length(phylo$edge[,1]))),
       ## 16. $node.label
       as.character(vector("character", length(phylo$edge[,1])))
       )
    
    dyn.unload("phylomatic.so")

    ## trim $tip.label
    ape <- list(Nnode = out[[10]],
                edge = cbind(out[[11]], out[[12]]),
                tip.label = out[[15]][out[[15]] != ""],
                node.label = out[[16]][out[[16]] != ""])
    
    if (length(out[[13]][!is.nan(out[[13]])]) != 0) {
        out[[13]][is.nan(out[[13]])] <- NA
        ape <- c(ape, list(edge.length = out[[13]]))
    }
    if (!is.nan(out[[14]])) {
        ape <- c(ape, list(root.edge = out[[14]]))
    }

    class(ape) <- "phylo"
    return(ape)
}




