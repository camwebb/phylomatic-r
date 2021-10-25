#include <R.h>
#include <string.h> // for strcpy

void phylomatic(int *ape_nedge,  int *ape_to,  int *ape_from,           \
                double *ape_bl,  double *ape_rbl,                       \
                int *ape_ntip,  char **ape_tip_lab,                     \
                char **ape_node_lab, int *ape_max_lab_len,              \
                                                                        \
                int *out_nnode, int *out_to,  int *out_from,            \
                double *out_bl,  double *out_rbl,                       \
                char **out_tip_lab, char **out_node_lab, int *out_ntip)
{
  // NB, this APE -> FY conversion could go in the R code, but far easier
  // to use loops and arrays
  int i;
  int n = 0;  // node counter
  int j = 0;  // tip label counter
  int k = 1;  // node label counter (skip root)

  int ot = 0;  // out tip counter
  int on = 0;  // out node counter

  // BEWARE index of a2f[]: node values are 1 to nedge+1, [0] not used
  int    a2f[*ape_nedge +1 +1];
  // All Fy index begin with 0=root and run to ape_nedge 
  int    f2a[*ape_nedge+1];
  int    fy_to[*ape_nedge+1];
  int    fy_isinner[*ape_nedge+1];
  double fy_length[*ape_nedge+1];
  char   fy_lab[*ape_nedge+1][*ape_max_lab_len+1];

  // Start with root node
  fy_to[0] = -1; 
  fy_isinner[0] = 1;
  fy_length[0] = *ape_rbl;
  strcpy(fy_lab[0], ape_node_lab[0]);
  a2f[ape_to[0]] = 0;  // ape-to-fy lookup: edge is root
  // Rprintf("%2d %2d %2d %3.1f %8s\n",                             \
  //         n, fy_to[n], fy_isinner[n], fy_length[n], fy_lab[n] );
  
  // enumerate through ape's edge list
  for(i = 0; i < *ape_nedge; i++) {
    n++;
    // is it a tip node? is the edge[,1] value <= number of tips?
    if (ape_from[i] <= *ape_ntip) {
      fy_isinner[n] = 0;
      strcpy(fy_lab[n], ape_tip_lab[j]);
      j++;
    }
    else {
      fy_isinner[n] = 1;
      strcpy(fy_lab[n], ape_node_lab[k]);
      k++;
    }
    // copy $edge and branch lengths
    a2f[ape_from[i]] = n;
    fy_to[n]   = a2f[ape_to[i]];  // will always have been filled in prior loop
    fy_length[n]  = ape_bl[i];
    
    // Rprintf("%2d %2d %2d %3.1f %8s\n", \
    //   n, fy_to[n], fy_isinner[n], fy_length[n], fy_lab[n] );
  }

  // ----- fy2ape ----------------------------------------------------------
  // root info:
  *out_rbl = fy_length[0] ;
  strcpy(out_node_lab[0], fy_lab[0]);
  f2a[0] = *ape_ntip + 1;
  
  // for each member of fy, excluding the first
  for(i = 1; i < *ape_nedge+1; i++) {
    // tips
    if(fy_isinner[i] != 1) {
      ot++;
      out_from[i-1] = ot;
      f2a[i] = ot;
      strcpy(out_tip_lab[ot-1], fy_lab[i]);
    }
    else {
      on++;
      f2a[i] = on + *ape_ntip + 1;
      strcpy(out_node_lab[on], fy_lab[i]);
    }
    // Rprintf("%2d %2d %3.1f\n", f2a[fy_to[i]],f2a[i], fy_length[i]);

    out_from[i-1] = f2a[i];
    out_to[i-1]   = f2a[fy_to[i]];
    out_bl[i-1] = fy_length[i];
  }
  
  *out_ntip = ot;
  *out_nnode = on+1;
}

