#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "readConfig.h"

_Complex double ****readConfiguration(char path[100], int iconf, _Complex double ****U, int *Lsize, int THREAD_NUM){

  /* input FILE */
  char filename[100];          /* filename configuration */
  FILE *f;                    /* file where configuration is saved in */

  snprintf(filename, 100, path, iconf);

  fprintf(stderr,
	  "%i| conf %i: ...Read configuration: %s...\n", 
	  THREAD_NUM, iconf, filename);

  /* open file and read from it */
  f = fopen(filename, "r");           /* file config. is written in */
  U = readConfChristian(f, U, Lsize);/* read config. from file */
  fclose(f);

  return U;

}

_Complex double ****readConfElena(FILE *f, _Complex double ****U, int *Lsize){
  int nold, mu;
  int n, x, y, z, t;
  int zVol   = Lsize[4];
  int yzVol  = Lsize[3]*Lsize[4];
  int xyzVol = Lsize[2]*yzVol;

  do{

    /* convention in file: z y x t */
    fscanf(f, "x = %i mu = %i\n", &nold, &mu);
    t = nold / xyzVol;
    x = (nold % xyzVol) / yzVol;
    y = ((nold % xyzVol) % yzVol) / zVol;
    z = ((nold % xyzVol) % yzVol) % zVol;
    n = t + Lsize[1]*x + Lsize[1]*Lsize[2]*y + Lsize[1]*Lsize[2]*Lsize[3]*z;

    fscanf(f, "[ %lf + %lf * i,  %lf + %lf * i, %lf + %lf * i;\n", 
	   &(__real__ U[n][mu][0][0]), &(__imag__ U[n][mu][0][0]), 
	   &(__real__ U[n][mu][0][1]), &(__imag__ U[n][mu][0][1]), 
	   &(__real__ U[n][mu][0][2]), &(__imag__ U[n][mu][0][2]));
    fscanf(f, "  %lf + %lf * i,  %lf + %lf * i, %lf + %lf * i;\n", 
	   &(__real__ U[n][mu][1][0]), &(__imag__ U[n][mu][1][0]), 
	   &(__real__ U[n][mu][1][1]), &(__imag__ U[n][mu][1][1]), 
	   &(__real__ U[n][mu][1][2]), &(__imag__ U[n][mu][1][2]));
    fscanf(f, "  %lf + %lf * i,  %lf + %lf * i, %lf + %lf * i  ]\n", 
	   &(__real__ U[n][mu][2][0]), &(__imag__ U[n][mu][2][0]), 
	   &(__real__ U[n][mu][2][1]), &(__imag__ U[n][mu][2][1]), 
	   &(__real__ U[n][mu][2][2]), &(__imag__ U[n][mu][2][2]));

    }while(!feof(f));
    return U;
}

_Complex double ****readConfChristian(FILE *f, _Complex double ****U, int *Lsize){
  int mu, i;
  int n, x, y, z, t;

  do{

    /* convention in file: z y x t */
    fscanf(f, "x: %i, y: %i, z: %i, t: %i\n", &x, &y, &z, &t);
    n = t + Lsize[1]*x + Lsize[1]*Lsize[2]*y + Lsize[1]*Lsize[2]*Lsize[3]*z;

    for(i=0; i<4; i++){
      fscanf(f, "mu: %i\n", &mu);

      fscanf(f, "%lf +i*%lf %lf +i*%lf %lf +i*%lf %lf +i*%lf %lf +i*%lf %lf +i*%lf %lf +i*%lf %lf +i*%lf %lf +i*%lf\n", 
	     &(__real__ U[n][mu][0][0]), &(__imag__ U[n][mu][0][0]), 
	     &(__real__ U[n][mu][0][1]), &(__imag__ U[n][mu][0][1]), 
	     &(__real__ U[n][mu][0][2]), &(__imag__ U[n][mu][0][2]),
	     &(__real__ U[n][mu][1][0]), &(__imag__ U[n][mu][1][0]), 
	     &(__real__ U[n][mu][1][1]), &(__imag__ U[n][mu][1][1]), 
	     &(__real__ U[n][mu][1][2]), &(__imag__ U[n][mu][1][2]), 
	     &(__real__ U[n][mu][2][0]), &(__imag__ U[n][mu][2][0]), 
	     &(__real__ U[n][mu][2][1]), &(__imag__ U[n][mu][2][1]), 
	     &(__real__ U[n][mu][2][2]), &(__imag__ U[n][mu][2][2]));
    }
  }while(!feof(f));
  return U;
}

/* _Complex double ****freeTheory(_Complex double ****U, int N){ */
/*   int a, c, d; */
/*   for (d=0; d<N; d++) { */
/*     for (c=0; c<4; c++) { */
/*       for (a=0; a<3; a++) { */
/* 	U[d][c][a][a] = 1.; */
/*       } */
/*     } */
/*   } */
/*   return U; */
/* } */
