#include <R.h>
#include <stdio.h>  // for NULL
#include <string.h> // for strcpy

void phylomatic(int *ape_nedge,  int *ape_to,  int *ape_from,           \
                double *ape_length,  double *ape_rlength,               \
                int *ape_ntip,  char **ape_tip_lab,                     \
                char **ape_node_lab, int *ape_max_lab_len)
{
  // NB, this APE -> FY conversion could go in the R code, but far easier
  // to use loops and arrays
  int i;
  int j = 0;  // tips used
  int k = 1;  // inner used (skip root)

  int    fy_id[*ape_nedge+1];
  int    fy_to[*ape_nedge+1];
  int    fy_isinner[*ape_nedge+1];
  double fy_length[*ape_nedge+1];
  char   fy_lab[*ape_nedge+1][*ape_max_lab_len+1];

  // enumerate through ape's edge list 
  for(i = 0; i < *ape_nedge; i++) {
    Rprintf("node [%d] %d ", i, ape_from[i] );
    // is it a tip node? is the edge[,1] value <= number of tips?
    if (ape_from[i] <= *ape_ntip) {
      Rprintf("terminal\n");
      fy_isinner[i] = 0;
      strcpy(fy_lab[i], ape_tip_lab[j]);
      j++;
    }
    else {
      Rprintf("inner\n");
      fy_isinner[i] = 1;
      strcpy(fy_lab[i], ape_node_lab[k]);
      k++;
    }
    // copy $edge and branch lengths
    fy_id[i] = ape_from[i];
    fy_to[i]   = ape_to[i]; 
    fy_length[i]  = ape_length[i];
  }
  
  // add a root node at end
  Rprintf("node [%d] %d inner\n", *ape_nedge, *ape_ntip+1);

  fy_id[*ape_nedge] = *ape_ntip+1; // this is the temporary root id
  fy_to[*ape_nedge] = -1; 
  fy_isinner[*ape_nedge] = 1;
  fy_length[*ape_nedge] = *ape_rlength;
  strcpy(fy_lab[*ape_nedge], ape_node_lab[0]);

}

//                                       labels            
//   APE TREE                    APE   tip  node        FY
//                                                
//       root = 4              4 <- 1   A             4  1 A   
//     /\ inner = 5            4 <- 5      inner      4  5 inner
//    / /\                     5 <- 2   B             5  2 B
//   A  B C                    5 <- 3   C             5  2 C
//   1  2 3                                           4 -1 root
