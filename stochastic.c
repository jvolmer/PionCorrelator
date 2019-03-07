#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "stochastic.h"
#include "matrix.h"

_Complex double ***getStochasticSource(_Complex double ***xi, int mt, int Lvol, int *Lsize){

  int j, a, alpha;    /* counting vars */
  int jt;             /* time of j */
  int r;              /* random var */
  double sqrt2 = 1./ sqrt(2.);
  double si = 0., co=0.;

  xi = setZero3D(xi, Lvol);
	    
  for (j=0; j<Lvol; j++){          /* get stochastic noise vector */
    jt = ( ( j % (Lsize[1]*Lsize[2]*Lsize[3])
	     )   % (Lsize[1]*Lsize[2])
	   )     % Lsize[1];
    if (jt == mt) {                /* random time-slice source */
      for (alpha=0; alpha<4; alpha++){
	for (a=0; a<3; a++){

	  r = rand() % 4;

	  if      (r == 0){
	    si = sqrt(2);
	    co = sqrt(2);
	  }
	  else if (r == 1){
	    si = -sqrt(2);
	    co = sqrt(2);
	  }
	  else if (r == 2){
	    si = sqrt(2);
	    co = -sqrt(2);
	  }
	  else{
	    si = -sqrt(2);
	    co = -sqrt(2);
	  }

	  xi[j][alpha][a] = co + si * I;

	} // end a
      } // end alpha
    } // if jt == mt
  } //end j

  return xi;
}
