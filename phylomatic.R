
phylomatic <- function(phylo, taxa) {

    ## edges?
    if(length(phylo$edge.length) == 0) {
        phylo$edge.length <- vector("double", length(phylo$edge[,1]))
        phylo$edge.length[] <- NA
    }
    ## recode NaNs to NAs
    phylo$edge.length[is.nan(phylo$edge.length)] <- NA
    ## root edge?
    if(length(phylo$root.edge) == 0) {
        phylo$root.edge <- vector("double", 1)
        phylo$root.edge[] <- NA
    }
    ## ## number of nodes for FY
    ## fn     <- length(phylo$edge[,1]) + 1

    dyn.load("phylomatic.so")
    
    .C("phylomatic",
       NAOK = TRUE,
       ## INPUT
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
       as.integer(max(nchar(c(phylo$tip.label,phylo$node.label))))

       
       ## ## OUTPUT
       ## ## 9. Node identifier
       ## as.integer(vector("integer", fn)),
       ## ## 10. Node 'to'  
       ## as.integer(vector("integer", fn)),
       ## ## 11. Inner node?
       ## as.integer(vector("integer", fn)),
       ## ## 12. Branch length
       ## as.double(vector("double",fn)),
       ## ## 13. Label
       ## as.character(vector("character",fn))
       )
}



ape2fy <- function(t) {

    f <- data.frame(cbind(t$edge[,2], t$edge[,1]))
    colnames(f) <- c("id","to")
    
    ## edges
    if(length(t$edge.length) == 0) {
        f$bl <- vector("double", length(t$edge[,1]))
        t$bl[] <- NA
    } else {
        f$bl <- t$edge.length
        ## recode NaNs to NAs
        f$bl[is.nan(f$bl)] <- NA
    }

    ## is it a tip node? is the edge[,1] value <= number of tips?
    f$inner <- vector("integer", length(t$edge[,1]))
    f$inner <- 1
    f$inner[f$id <= length(t$tip.label)] <- 0

    f$at <- vector("character", length(t$edge[,1]))

    ## tips first
    b <- data.frame(a = 1:length(t$tip.label), b = t$tip.label)
    f$at <- b$b[match(t$edge[,2], b$a)] # or use merge()

    # the rest are filled with node labels, in order
    f$at[is.na(f$at)] <- t$node.label[-1]

    ## root node
    ## root edge?
    if(length(t$root.edge) == 0) {
        re <- NA
    } else {
        re <- t$root.edge
    }
    ## new row for root
    f <- rbind(f, data.frame(id = length(t$tip.label)+1, to = -1,
                             bl = re, inner = 1, at = t$node.label[1]))
    
    ## # sort?
    ## f <- f[order(f$id),]
 
    return(f)
}





## back to APE: t2 <- list(edge = cbind(c(4,5,5,4),c(5,1,2,3)), Nnode=2, node.label=c("root","inner"), tip.label=c("A","B","C"))
## class(t2) <- "phylo"



