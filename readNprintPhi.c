#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "readNprintPhi.h"
#include "matrix.h"

void printPhi(char path[100], char ps, _Complex double ***phi, int iconf, int nxi, int mt, int Lvol){

  /* output FILE */
  FILE *fout;                 /* file to save result */
  char what[10];
  
  if (ps == 'p') strcpy(what, "phi");
  if (ps == 's') strcpy(what, "source");
  
  fprintf(stderr,
	  "conf %i.%i: ...Print %s to some file: %s...\n",
	  iconf, nxi, what, path);

  /* open and print */
  fout = fopen(path, "w");
  fprintf(fout, "conf = %i, Nsource = %i, mt = %i\n\n", iconf, nxi, mt);
  fPrint3D(fout, phi, Lvol, 4, 3);
  
  fclose(fout);

}


int readPhi(char file[100], char ps, _Complex double ***phi, int nxi, int iconf, int Nxi){

  char fileNumber[10];
  char fileName[100];
  FILE *fout;
  int confread, sourceread;
  int mt;
  char what[10];

  if (ps == 'p') strcpy(what, "phi");
  if (ps == 's') strcpy(what, "source");

    /* open propagator file */
    strcpy(fileName, file);
    snprintf(fileNumber, 10, ".%d.%d.%d", iconf, Nxi, nxi);
    strcat(fileName, fileNumber);
    /* get proper file name for output file: conf.nr.source */
    fout = fopen(fileName, "r");
    fprintf(stderr,
	    "conf %i.%i: ...Read %s from file: %s...\n", 
	    iconf, nxi, what, fileName);

    /* read propagator */
    phi = readmyOutput(fout, &confread, &sourceread,
		     &mt, phi); /* read prop */

    /* close propagator file */
    fclose(fout);

    /* check whether its the right file */
    if (confread != iconf){     /* if read vars do not fit 
				   to counting vars */
      fprintf(stderr,
	      "conf from file does not fit to loop var\n");
      fprintf(stderr, "conf_file = %i\n\n",
	      confread);
      fprintf(stderr, "conf_prog = %i\n\n",
	      iconf);
    }
    if (sourceread != nxi){     /* if read vars do not fit 
				   to counting vars */
      fprintf(stderr,
	      "source from file does not fit to loop var\n");
      fprintf(stderr, "source_file = %i\n\n",
	      sourceread);
      fprintf(stderr, "source_prog = %i\n\n",
	      nxi);
    }

  return mt;

}

/* read my propagator files */
/* to get prop G[n][alpha][a] */
/* additional: get time slice of source mt, sourcenr Nxi and confnr conf */
_Complex double ***readmyOutput(FILE *f,
				int *conf, int *nxi, int *m0,
				_Complex double ***out 
				){
  int n_read, alpha_read, a_read;
  /* int a0, alpha0, m0; */

  fscanf(f, "conf = %i, Nsource = %i, mt = %i\n\n", 
	 &(*conf), &(*nxi), &(*m0));
  fprintf(stderr, "conf %i.%i: mt = %i\n", *conf, *nxi, *m0);

  do{
    /* convention in file: z y x t */
    fscanf(f, "A[n=%i][alpha=%i][a=%i] = ", &n_read, &alpha_read, &a_read);
    fscanf(f, "%lf +i*%lf\n",
	   &(__real__ out[n_read][alpha_read][a_read]), 
	   &(__imag__ out[n_read][alpha_read][a_read])
	   );

  }while(!feof(f));

  return out;
}
