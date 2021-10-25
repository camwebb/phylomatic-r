# Development notes

This is the first time I’ve tried to develop a C extension to R. Here
are a few notes, FYI.

## APE data structure

The `ape` `phylo` object is a bit tricky to work with, compared to my
`fy` format.  I left the internals of `phylomatic()` using the `fy`
structure, so I needed to convert ape-to-fy and fy-to-ape.

```
  <----------- APE -------------->     <--------- FY -------------->
                          .label       id   to  at      
                        $tip $node        
                                       0 -> -1 root
      4 ‘root’    4 <- 1   A           1 ->  0 A            0 ‘root’  
     /\ 5 ‘in’    4 <- 5      in       2 ->  0 in          /\ 2 ‘x’
    / /\          5 <- 2   B           3 ->  2 B          / /\
   1 2  3         5 <- 3   C           4 ->  2 C         1 3  4
  ‘A’‘B’’C’                                             ‘A’‘B’‘C’
```

## Debugging

See:

 * <https://www.stat.purdue.edu/~liu105/STAT598G_lab/lab6.pdf>
 * <https://www.maths.ed.ac.uk/~swood34/RCdebug/RCdebug.html>
 
### gdb

 * Compile with `MAKEFLAGS="CFLAGS=-g\ -O1" R CMD SHLIB foo.c`

```
R -d gdb
...
b phylomatic
y
run
... R commands to get to phylomatic
n n n n...
p i[1]
display i[1]
display i[2]
n n n n...
...until crash
```

### valgrind

`R -d valgrind`

Then watch for memory errors


