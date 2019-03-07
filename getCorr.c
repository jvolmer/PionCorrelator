#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "getCorr.h"
#include "readNprintPhi.h"
#include "matrix.h"

#define PI 3.141592653589793238462643383279502884


void getCorr(char *file_prop, char *file_corr, int *Lsize, int Nxi){

  int m, confread, iconf;
  int mt, nt;
  int nxi;

  _Complex double ***phi, *sum;

  FILE *fout, *fcorr;

  int Lvol = Lsize[1]*Lsize[2]*Lsize[3]*Lsize[4];
  int T = Lsize[1];


  /* allocation */
  phi = Alloc3D(Lvol, 4, 3);
  sum = AllocComplexVector(T);

  fcorr = fopen(file_corr, "w");
  
  /* inverted configurations */
  for (iconf=1000; iconf<=1000; iconf++){ //1095; iconf=iconf+5){

    for (nxi=0; nxi<Nxi; nxi++){

      /* read produced solution for nxi and iconf */
      mt = readPhi(file_prop, 'p', phi, nxi, iconf, 1);

      /* calculate trace for fixed conf and source */
      sum = doSums(sum, phi, mt, nxi, iconf, Lsize);

      /* print correlator to file (because of stochastic sources:
                                   average over Nxi and L3) */
      printCorr(fcorr, sum, nxi, iconf, Lsize);

    } // end nxi


  } // end iconf

  free(sum);
  free(phi);

}


_Complex double *doSums(_Complex double *sum, 
			_Complex double ***phi, 
			int mt, int nxi, int iconf, int *Lsize){

  int n, nt, dt, beta, b;
  int T = Lsize[1];
  int Lvol = Lsize[1]*Lsize[2]*Lsize[3]*Lsize[4];
  _Complex double kahan_new[T], kahan_add[T], kahan_rest[T];

  fprintf(stderr, "conf %i.%i: ...Calculate phi...\n", iconf, nxi);

  for (nt=0; nt<T; nt++){             /* set added vars to zero */
    sum[nt] = 0.;
    kahan_rest[nt] = 0.;
  }
  
  for (n=0; n<Lvol; n++){             /* over all sink lattice points */
    
    nt = ( ( n % (Lsize[1]*Lsize[2]*Lsize[3]) /* get time of */
	     )   % (Lsize[1]*Lsize[2])          /* lattice point */
	   )     % Lsize[1];
    
    if (nt-mt >= 0) dt = nt - mt;     /* time difference of sink&source */
    else dt = nt-mt+T;
    
    for (beta=0; beta<4; beta++){     /* DC trace (sink) */
      for (b=0; b<3; b++){
	
	/* add G(n|m)_{alpha, a, beta, b}^2 with kahan summation (better 
	   floating point precision) */
	kahan_new[dt] = conj(phi[n][beta][b]) * phi[n][beta][b] - kahan_rest[dt];
	kahan_add[dt] = sum[dt] + kahan_new[dt];
	kahan_rest[dt] = (kahan_add[dt] - sum[dt]) - kahan_new[dt];
	sum[dt] = kahan_add[dt];
	
      } // end b
    } // end beta
    
  } // end n
  
  return sum;
  
}


void printCorr(FILE *fout, _Complex double *sum, int nxi, int iconf, int *Lsize){

  int T = Lsize[1];
  int L3 = Lsize[2]*Lsize[3]*Lsize[4];

  int nt;
  double norm;                        /* normalization of correlator */
  
  
  fprintf(stderr,
	  "conf %i.%i: ...Write correlator to file...\n",
	  iconf, nxi);
  
  norm = (double)L3;
  
  for (nt=0; nt<T; nt++){          /* print correlator/trace */

    fprintf(fout, "%i\t%i\t%i\t%.20f\t%.20f\n", iconf, nxi+1, nt,
  	    creal(sum[nt]/norm),
  	    cimag(sum[nt]/norm) );

    
  }
  


}
