#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <omp.h>
#include "readConfig.h"
#include "matrix.h"
#include "lattice.h"
#include "bicg-algo.h"
#include "operator.h"
#include "readNprintPhi.h"
#include "getCorr.h"
#include "stochastic.h"

#ifdef _OPENMP
   #include <omp.h>
   #define THREAD_NUM omp_get_thread_num()
   #define THREADS omp_get_num_procs()
#else
   #define THREAD_NUM 0
   #define THREADS 1
#endif

int main(int argc, char *argv[]){

  /* LATTICE */
  int Ldim;     /* LATTICE dimension */
  int Lvol;     /* LATTICE volume */
  int i;        /* counting var */
  int *Lsize;   /* LATTICE side lengths 1..Ldim+1 */
  int **nn;     /* LATTICE next neighbors
                           nn[k][i] next neighbor of i in kth dimension
                                    k: 1..Ldim         pos direction
				    k: Ldim+1..2Ldim+1 neg direction
				    i: 0..Lvol                           */

  /* CONFIGURATION */
  /* _Complex double ****U; */ /* CONFIGURATION
		      U[n][mu][a][b] link at n in mu direction,
		                     entry a,b in SU(3) matrix
				     n:  0..Lvol-1 lattice index
				     mu: 0..3      lorentz index
				     a:  0..2      color index
				     b:  0..2      color index    */
  
  /* PROPAGATOR */
  /* _Complex double ***phi; */ /* propagator
			     phi[n][alpha][a] 
			     f-propagator at n for alpha, a
				     n:     0..Lvol-1 lattice index
				     alpha: 0..3      lorentz index
				     a:     0..2      color index    */
  
  /* CG - VAR */
  double tol = 1E-16;
  int maxIt = 100000;

  /* IMPORTANT QUANTITIES */
  int ud = 0;                 /* D: 0 d-quark, 1 u-quark */
  int Nxi;                    /* number of sources */       

  char file_corr[100];         /* file to save correlator */
  char file_prop[100];       /* file of saved propagators (without number) */
  char file_source[100];     /* file of saved source (without number) */

  /* INPUT PARAMETERS */
  if (argc != 3){
    printf("./run_dis.x: Wrong number of input parameters\n\n");
    printf("Usage: ./run_dis.x <COMMAND> <Nsources>\n\n");
    printf("                   <COMMAND> = p(his) -> invert Dirac matrix with stochastic sources\n");
    printf("                   <COMMAND> = c(orr) -> calculate correlator with one-end trick (phi*phi)\n\n");
    printf("                             <Nsources> = 1 ... 9\n\n");
    return 1;
  }
  Nxi = argv[2][0] - '0';             /* get number of sources */
  if (Nxi == 0){
    printf("./run_dis.x: Wrong input parameter\n\n");
    printf("Usage: ./run_dis.x <COMMAND> <Nsources>\n\n");
    printf("                             <Nsources> = 1 ... 9\n\n");
    return 1;
  }

  /* LATTICE */
  Ldim = 4;                         /* lattice dim = 1 time + 3 space */
  Lsize = AllocIntVector(Ldim + 1); /* set lattice lenghts */
  Lsize[1] = 32;                    /* set temporal lattice length */
    Lsize[2] =                      /* set spacial lattice lengths */
    Lsize[3] =
    Lsize[4] = 16;
  Lvol = 1;                         /* compute lattice volume */
  for (i=1; i<Ldim+1; i++)
    Lvol *= Lsize[i];

  nn = AllocIntMatrix(2*Ldim + 1, Lvol); /* set next neighbor matrix */
  nn = lattice(nn, Ldim, Lsize);

  printf("THREADS = %i\n", THREADS);

  /* CALCULATE PROPAGATORS */
  if (*argv[1] == 'p'){

    fprintf(stderr,"\n----- Start: calculating propagators with %i sources -----\n", Nxi);

	int iconf;                        /* over used configurations */
	srand(time(0));                   /* random number seed */
	
	for (iconf=1000; iconf<=1000; iconf++){ //1095; iconf=iconf+5){
	  
	  /* COUNTING VARS */
	  int mt;                     /* random time-slice per conf */
	  int nxi;                    /* counting of sources */
	  
	  _Complex double ***phi;
	  _Complex double ***xi;     /* stochastic source */
	  _Complex double ****U;

	  int cg=10;                  /* returnvalue of bicg-algo */
	  char path[100];             /* path of config, prop */
	  
	  /* ALLOCATION */
	  phi = Alloc3D(Lvol, 4, 3);  /* inversion matrix */
	  xi = Alloc3D(Lvol, 4, 3);    /* stochastic noise vector */
	  U = Alloc4D(Lvol, 4, 3, 3);  /* configuration matrix */
	  

	  /* read configuration */
	  fprintf(stderr, 
		  "\nconf %i: ============ Configuration %i =============\n",
		  iconf, iconf);

	  snprintf(path, 100, "input/conf.%d.gauge", iconf);
	  U = readConfiguration(path, iconf, U, Lsize, THREAD_NUM);


	  /* random source time slice */
	  mt = rand() % Lsize[1];
	  

	  for (nxi=0; nxi<Nxi; nxi++){      /* loop over sources */

	    while ((cg != 0) && (cg != 1)) {  /* if cg didnt work do it again
						 with new source */
	      
	      fprintf(stderr,
		      "\nconf %i.%i: ------- CG-Algo for stochastic source %i -------\n",
		      iconf, nxi, nxi);
	      
	      /* get stochastic source */
	      xi = getStochasticSource(xi, mt, Lvol, Lsize);
	      
	      /* invert matrix with cg */
	      phi = setZero3D(phi, Lvol);   /* set phi to zero */
	      cg = BiCG_Complex_Algo(ud, iconf, nxi, phi, xi, U, nn,
				     Lsize[1], Lvol, tol, maxIt);
	      
	      /* evaluate cg outcome */
	      BiCG_error(cg, THREAD_NUM, iconf, nxi);
	      
	    } /* if cg != 0 or 1 do this configuration again */
	    
	    
	    /* End */
	    fprintf(stderr, 
		    "\nconf %i.%i: ------------- End ---------------\n", 
		    iconf, nxi);


	    /* print propagator to some file */
	    snprintf(path, 100, "output/propagator/prop_timeslice.%d.%d.%d", iconf, Nxi, nxi);
	    printPhi(path, 'p', phi, iconf, nxi, mt, Lvol);

	    /* print source to some file */
	    snprintf(path, 100, "output/source/source_timeslice.%d.%d.%d", iconf, Nxi, nxi);
	    printPhi(path, 's', xi, iconf, nxi, mt, Lvol);

	    cg = 10;                          /* cg should start again
						 normally for next source */

	    
	  } // end nxi

	   
	  free(U);
	  free(phi);
	  free(xi);
	  
	} // end iconf
	

    } // end calculating propagator


  /* CALCULATE CORRELATOR using the one-end trick */
  else if (*argv[1] == 'c'){ 

    snprintf(file_prop, 100, "output/propagator/prop_timeslice");
    snprintf(file_corr, 100, "output/correlator/corr_timeslice.%d", Nxi);
    /* read in propagators, calculate correlator, print it to file */
    getCorr(file_prop, file_corr, Lsize, Nxi);

    fprintf(stderr, "\n ====> Correlator printed to file %s\n", file_corr);

  }

  /* END */
  free(nn);
  free(Lsize);

  return 0;
}
