_Complex double *doSums(_Complex double *sum, 
			_Complex double ***phi, 
			int mt, int Nxi, int iconf, int *Lsize);

void getCorr(char *file_prop, char *file_corr, int *Lsize, int Nxi);

void printCorr(FILE *fout, _Complex double *sum, int Nxi, int iconf, int *Lsize);
