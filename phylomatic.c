#include <R.h>
#include <stdio.h>  // for NULL
#include <string.h> // for strcpy

void phylomatic(int *ape_nedge,  int *ape_to,  int *ape_from,           \
                double *ape_length,  double *ape_rlength,               \
                int *ape_ntip,  char **ape_tip_lab,                     \
                char **ape_node_lab, int *ape_max_lab_len,              \
                                                                        \
                int *out_nnode, int *out_to,  int *out_from,            \
                double *out_length,  double *out_rlength,               \
                char **out_tip_lab, char **out_node_lab)
{
  // NB, this APE -> FY conversion could go in the R code, but far easier
  // to use loops and arrays
  int i;
  int n = 0;  // node counter
  int j = 0;  // tip label counter
  int k = 1;  // node label counter (skip root)

  int    a2f[*ape_nedge+1];
  int    fy_to[*ape_nedge+1];
  int    fy_isinner[*ape_nedge+1];
  double fy_length[*ape_nedge+1];
  char   fy_lab[*ape_nedge+1][*ape_max_lab_len+1];

  // Start with root node
  fy_to[0] = -1; 
  fy_isinner[0] = 1;
  fy_length[0] = *ape_rlength;
  strcpy(fy_lab[0], ape_node_lab[0]);
  a2f[ape_to[0]] = 0;  // ape-to-fy lookup: edge is root
  Rprintf("%2d %2d %2d %3.1f %8s\n", n, fy_to[n], fy_isinner[n], fy_length[n], fy_lab[n] );
  
  // enumerate through ape's edge list
  for(i = 0; i < *ape_nedge; i++) {
    n++;
    // is it a tip node? is the edge[,1] value <= number of tips?
    if (ape_from[i] <= *ape_ntip) {
      // Rprintf("terminal\n");
      fy_isinner[n] = 0;
      strcpy(fy_lab[n], ape_tip_lab[j]);
      j++;
    }
    else {
      // Rprintf("inner\n");
      fy_isinner[n] = 1;
      strcpy(fy_lab[n], ape_node_lab[k]);
      k++;
    }
    // copy $edge and branch lengths
    a2f[ape_from[i]] = n;
    fy_to[n]   = a2f[ape_to[i]];  // will always have been filled in prior loop 
    fy_length[n]  = ape_length[i];
    
    Rprintf("%2d %2d %2d %3.1f %8s\n", n, fy_to[n], fy_isinner[n], fy_length[n], fy_lab[n] );
  }

  // **fy2new**
  *out_nnode = *ape_nedge ;
  *out_rlength = fy_length[0] ;
  strcpy(out_node_lab[0], fy_lab[0]);

  int ot = 0;  // out tip counter
  int on = 0;  // out node counter
  // for each member of fy, excluding the first
  for(i = 1; i < *ape_nedge+1; i++) {
    // tips
    if(fy_isinner[i] != 1) {
      ot++;
      out_from[i-1] = ot;
      strcpy(out_tip_lab[ot-1], fy_lab[i]);
    }
    else {
      on++;
      out_from[i-1] = i +1 + *ape_ntip;
      strcpy(out_node_lab[on], fy_lab[i]);
    }
    out_to[i-1] = fy_to[i] +1+ *ape_ntip;
    out_length[i-1] = fy_length[i];
  }

  
}

//                                       labels            
//   APE TREE                    APE   tip  node        FY
//                                                
//       root = 4              4 <- 1   A             1 ->  4 A   
//     /\ inner = 5            4 <- 5      inner      5 ->  4 inner
//    / /\                     5 <- 2   B             2 ->  5 B
//   A  B C                    5 <- 3   C             3 ->  5 C
//   1  2 3                                           4 -> -1 root
