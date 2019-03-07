#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "operator.h"

#define PI 3.141592653589793238462643383279502884

/* Matrices (1 +- gamma_mu)_{alpha,beta}: 
   EinsMinusGamma[mu][alpha][beta],  */
    
  _Complex double EinsMinusGamma[4][4][4] = {
				   { {1 + 0*I, 0 + 0*I, 1 + 0*I, 0 + 0*I},
	 /*gamma4*/		     {0 + 0*I, 1 + 0*I, 0 + 0*I, 1 + 0*I},
				     {1 + 0*I, 0 + 0*I, 1 + 0*I, 0 + 0*I},
				     {0 + 0*I, 1 + 0*I, 0 + 0*I, 1 + 0*I}
				   },
				   { {1 + 0*I, 0 + 0*I, 0 + 0*I, 0 + 1*I},
	 /*gamma1*/		     {0 + 0*I, 1 + 0*I, 0 + 1*I, 0 + 0*I},
				     {0 + 0*I, 0 - 1*I, 1 + 0*I, 0 + 0*I},
				     {0 - 1*I, 0 + 0*I, 0 + 0*I, 1 + 0*I}
				   },
                                   { {1 + 0*I, 0 + 0*I, 0 + 0*I, 1 + 0*I},
				     {0 + 0*I, 1 + 0*I, -1 + 0*I, 0 + 0*I},
	 /*gamma2*/		     {0 + 0*I, -1 + 0*I, 1 + 0*I, 0 + 0*I},
				     {1 + 0*I, 0 + 0*I, 0 + 0*I, 1 + 0*I}
				   },
                                   { {1 + 0*I, 0 + 0*I, 0 + 1*I, 0 + 0*I},
	 /*gamma3*/		     {0 + 0*I, 1 + 0*I, 0 + 0*I, 0 - 1*I},
				     {0 - 1*I, 0 + 0*I, 1 + 0*I, 0 + 0*I},
				     {0 + 0*I, 0 + 1*I, 0 + 0*I, 1 + 0*I}
				   }
				 };

_Complex double EinsPlusGamma[4][4][4] = {
				   { {1 + 0*I, 0 + 0*I, -1 + 0*I, 0 + 0*I},
	 /*gamma4*/		     {0 + 0*I, 1 + 0*I, 0 + 0*I, -1 + 0*I},
				     {-1 + 0*I, 0 + 0*I, 1 + 0*I, 0 + 0*I},
				     {0 + 0*I, -1 + 0*I, 0 + 0*I, 1 + 0*I}
				   },
				   { {1 + 0*I, 0 + 0*I, 0 + 0*I, 0 - 1*I},
      	 /*gamma1*/		     {0 + 0*I, 1 + 0*I, 0 - 1*I, 0 + 0*I},
				     {0 + 0*I, 0 + 1*I, 1 + 0*I, 0 + 0*I},
				     {0 + 1*I, 0 + 0*I, 0 + 0*I, 1 + 0*I}
				   },
                                   { {1 + 0*I, 0 + 0*I, 0 + 0*I, -1 + 0*I},
				     {0 + 0*I, 1 + 0*I, 1 + 0*I, 0 + 0*I},
	 /*gamma2*/		     {0 + 0*I, 1 + 0*I, 1 + 0*I, 0 + 0*I},
				     {-1 + 0*I, 0 + 0*I, 0 + 0*I, 1 + 0*I}
				   },
                                   { {1 + 0*I, 0 + 0*I, 0 - 1*I, 0 + 0*I},
				     {0 + 0*I, 1 + 0*I, 0 + 0*I, 0 + 1*I},
	 /*gamma3*/		     {0 + 1*I, 0 + 0*I, 1 + 0*I, 0 + 0*I},
				     {0 + 0*I, 0 - 1*I, 0 + 0*I, 1 + 0*I}
				   }
				 };


_Complex double ***twistedMassToField(int ud,
				       _Complex double ***out,
				       _Complex double ***in, 
				       _Complex double ****U, 
				       int **nn, int T, int Lvol){

  /* CONSTANTS */
  double kappa = 0.160856;
  double mu_tw = 0.004;
  double C = 1./(2*kappa);

  int f, n, alpha, a;
  int beta, mu, b;
  _Complex double factor; /* diagonal factor */
  _Complex double anti;   /* anti-periodic b.c. in time */

  int flavorSign = 1; /* if f=0 */
  int pauliSign = 1;  /* if alpha <= 2 */

  if (ud == 0) flavorSign = 1;
  else flavorSign = -1;

  /* set DIAGONAL part */
#pragma omp parallel for private(n,alpha,a, pauliSign, factor), shared(out)
    for (n=0; n<Lvol; n++){              /* spacial loop */
      for (alpha=0; alpha<4; alpha++) {  /* lorentz loop */
	if (alpha >= 2) pauliSign = -1;  /* for alpha = 2,3 change sign */
	else pauliSign = 1;
	for (a=0; a<3; a++){             /* color loop */

	  factor = C + (pauliSign * flavorSign * mu_tw) * I;
	  out[n][alpha][a] = factor * in[n][alpha][a] ;

	  
	} // end a
      } // end alpha
    } // end n
  
  /* add OFF-DIAGONAL part */
#pragma omp parallel for private(n,alpha,a,mu,beta,b,anti), shared(out)
    for (n=0; n<Lvol; n++){              /* spacial loop */
      for (alpha=0; alpha<4; alpha++) {  /* lorentz loop */
	for (a=0; a<3; a++){             /* color loop */

	  for (mu=0; mu<4; mu++){         /* direction loop */
	    for (beta=0; beta<4; beta++){ /* lorentz loop */
	      for (b=0; b<3; b++){        /* color loop */
	  
		/* anti-p.b.c. in temporal direction */
		if (mu ==0)
		  anti = cos(PI/T) + I * sin(PI/T);
		else anti = 1;
		
		out[n][alpha][a] -= C * kappa *
		  ( EinsMinusGamma[mu][alpha][beta]
		    * anti*U[n][mu][a][b]
		    * in[nn[mu+1][n]][beta][b]
		    +  EinsPlusGamma[mu][alpha][beta]
		    * conj(anti*U[nn[4+mu+1][n]][mu][b][a])
		    * in[nn[4+mu+1][n]][beta][b]);

	      } // end b
	    } // end beta
	  } // end mu

	} // end a
      } // end alpha
    } // end n



  return out;

}

