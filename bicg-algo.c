#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "matrix.h"
#include "operator.h"
#include "bicg-algo.h"

void BiCG_finalize(_Complex double ***q, _Complex double ***v, _Complex double ***t, _Complex double ***r, _Complex double ***rtw, _Complex double ***p, _Complex double ***s){
  free(q); free(v); free(t);
  free(r); free(rtw);
  free(p); free(s);
}

void BiCG_error(int BiCG_out, int  threadnr, int iconf, int nxi){

  char errortext[100];
  
  if (BiCG_out == 0)
    sprintf(errortext,"BiStab CG converged");
  else if (BiCG_out == 1) 
    sprintf(errortext, "Test successful, input solves system :D");
  else if (BiCG_out == 2)
    sprintf(errortext, "BiStab CG is not working (rho = 0)!");
  else if (BiCG_out == 3)
    sprintf(errortext, "BiStab CG used maximal number of iterations!");
  else if (BiCG_out == 4)
    sprintf(errortext, "BiStab CG produced 'nan' for |r|^2!");
  else
    sprintf(errortext, "BiStab CG unknown error %i",
	    BiCG_out);

  fprintf(stderr,
	  "%i| conf %i.%i: %s\n",
	  threadnr, iconf, nxi, errortext);
}

int BiCG_Complex_Algo(int ud, int iconf, int nxi, _Complex double ***x, _Complex double ***xi, _Complex double ****U, int **nn, int T, int N, double tol, int maxIt) {

  _Complex double ***q;        /* 1. help vector q = A x */
  _Complex double ***v;        /* 2. help vector v = A p */
  _Complex double ***t;        /* 3. help vector t = A s*/

  _Complex double ***r;        /* carryover vector r = b - q */
  _Complex double ***rtw;      /* carryover vector at start */

  _Complex double ***p;        /* 1. translation vector */
  _Complex double ***s;        /* 2. translation vector */

  double s2, t2;                 /* vector products */
  _Complex double ts, rtwv;

  _Complex double alpha;         /* 1. coefficient */
  _Complex double omega;         /* 2. coefficient */
  _Complex double beta;          /* help coefficient */
  _Complex double rho;           /* parallelity coeff: 
				    checks sucess of method */
  _Complex double rho_old;       /* parallelity coeff of last iteration */
  double err, bnorm;             /* exit cond: err < tol*bnorm */

  int i,j,f,a,al;                /* counting vars */

  q   = Alloc3D(N,4,3);
  v   = Alloc3D(N,4,3);
  t   = Alloc3D(N,4,3);
  r   = Alloc3D(N,4,3);
  rtw = Alloc3D(N,4,3);
  p   = Alloc3D(N,4,3);
  s   = Alloc3D(N,4,3);

  q = twistedMassToField(ud, q, x, U, nn, T, N);     /* q = D x */

  for (j=0; j<N; j++) {
    for (al=0; al<4; al++){
      for (a=0; a<3; a++){
	  
	r[j][al][a] = xi[j][al][a] - q[j][al][a];    /* Rest r_0 = b - A x */
	  
	rtw[j][al][a] = r[j][al][a];      /* r_twidle = r_0 */
	p[j][al][a] = r[j][al][a];        /* p = r_0 */
      }
    }
  }

  
  rho = scalarProductComplex(rtw, r, N);    /* parallelity of actual and 
					       starting carryover vector 
					       rho = r_tw * r */

  bnorm = 1.;                               /* |b|^2 */

  for (i=0; i< maxIt; i++) {                // start iteration

    err = normSquared(r, N);                /* cg-error err = |r|^2 */

    fprintf(stderr,"conf %i.%i: %d, %e\n", 
	    iconf, nxi, i, err); /* output iteration, cg-error */

    if (err <= tol * bnorm){                /* sufficient exact */
      BiCG_finalize(q, v, t, r, rtw, p, s);   /* break */
      if (i==0) return(1);                    /* input x solves system */
      else return(0);                         /* normal exit */
    }
    else if (isnan(err)){
      BiCG_finalize(q, v, t, r, rtw, p, s); /* something went wrong */
      return(4);
    }

    v = twistedMassToField(ud, v, p, U, nn, T, N); /* v = D p */

    rtwv = scalarProductComplex(rtw, v, N); /* r_tw * v */
    alpha = rho / rtwv;                     /* first coefficient: 
					       alpha = rho / (rtw*v) */

    for (j=0; j<N; j++){
      for (al=0; al<4; al++){
	for (a=0; a<3; a++){
	  
	  s[j][al][a] = r[j][al][a] - alpha * v[j][al][a];
	                                    /* 2. translation vector:
					       s = r - alpha * v */
	}
      }
    }
    
    t = twistedMassToField(ud, t, s, U, nn, T, N); /* t = D s */
    t2 = normSquared(t, N);                 /* |t|^2 */
    ts = scalarProductComplex(t, s, N);     /* t * s */
    omega = ts / t2;                        /* second coefficient: 
					       omega = (t*s) / |t|^2*/

    for (j=0; j<N; j++){
      for (al=0; al<4; al++){
	for (a=0; a<3; a++){
	  
	  r[j][al][a] = s[j][al][a] - omega * t[j][al][a];
	                                    /* new carryover vector 
					       r = s - omega * t */
	  x[j][al][a] += alpha *p[j][al][a] + omega *s[j][al][a];
	                                    /* new and result vector 
					       x = x + alpha*p + omega*s */
	  
	}
      }
    }

    rho_old = rho;                          /* remember old rho */
    rho = scalarProductComplex(rtw, r, N);  /* rho = t_tw * r */

    if ( (fabs(creal(rho)) <= 1E-25) &&     /* algo fails if rho = 0 */
	 (fabs(cimag(rho)) <= 1E-25)){
      BiCG_finalize(q, v, t, r, rtw, p, s);
      return(2);
    }

    beta = (alpha * rho) / (omega * rho_old);/* help coefficient
						beta = (alpha * rho) / 
						       (omega * rho_old) */

    for (j=0; j<N; j++){
      for (al=0; al<4; al++){
	for (a=0; a<3; a++){
	  
	  p[j][al][a] = r[j][al][a] + beta *
	    (p[j][al][a] - omega * v[j][al][a]);
	                                      /* 1. translation vector
						 p = r + beta*(p-omega*v) */
	}
      }
    }
	  

  } // end iteration loop


  BiCG_finalize(q, v, t, r, rtw, p, s);       /* reached maxIt, */
  return(3);				      /* algo failed */

}
