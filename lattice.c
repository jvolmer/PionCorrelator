#include <stdio.h>
#include <stdlib.h>

int **lattice(int **nn, int ndim, int *lsize){
  /*
  Based on algorithm geom_pbs from B Bunk 12/2005, rev 4/2013

  Needs
      lattice dimension ndim
      lattice sizes     lsize[k], k=1..ndim

  Gives
      volume              nvol
      next-neighbor field nn[k][i], k=0..2*ndim, i=0..(nvol-1)
         next-neighbor of site i in direction k
	 for negative direction -k: ndim+k
	 nn[0][i] is reserved
  */
  int nvol;
  int  i, k;
  int *ibase, *ix; 

  ibase = (int *) malloc((ndim+2) * sizeof(int));
  ix = (int *) malloc((ndim+1) * sizeof(int));

  /* base for site indices */
  ibase[1] = 1;
  for (k=1; k<=ndim; k++)
    ibase[k+1] = ibase[k]*lsize[k];
  nvol = ibase[ndim+1];

  for (k=1; k<=ndim; k++) ix[k] = 0;   /* coordinates of starting site */

  for (i=0; i<nvol; i++){              /* loop over all sites */
    
    for (k=1; k<=ndim; k++){           /* loop over all neighbors e_k */
      
      nn[k][i] = i + ibase[k];            /* neighbor x + e_k */
      if (ix[k] == (lsize[k]-1)) {        /* at boundaries */
	nn[k][i] -= ibase[k+1];
	nn[0][i] = 1;                     /* for dirichlet */
      }

      nn[ndim+k][i] = i - ibase[k];       /* neighbor x - e_k */
      if (ix[k] == 0) {                   /* at boundaries */
	nn[ndim+k][i] += ibase[k+1];
	nn[0][i] = 1;                     /* for dirichlet */
      }

    }

    for (k=1; k<=ndim; k++){          /* coordinates of next site */
      ix[k]++;
      if (ix[k] < lsize[k]) break;
      ix[k] = 0;
    }
    
  }
  free(ibase); free(ix);
  return nn;
}
