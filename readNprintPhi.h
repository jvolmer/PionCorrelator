
void printPhi(char path[100], char ps, _Complex double ***phi, int iconf, int nxi, int mt, int Lvol);

/* read in all phis of one configuration (from all sources) */
int readPhi(char file[100], char ps, _Complex double ***phi, int nxi, int iconf, int Nxi);

/* read in output from me (from fprint4D) */
_Complex double ***readmyOutput(FILE *f,
				int *conf, int *nxi, int *m0,
				_Complex double ***out 
				);
