/* 
  cg-algo.h
*/


void BiCG_error(int BiCG_out, int  threadnr, int iconf, int nxi);
int BiCG_Complex_Algo(int ud, int iconf, int nxi, _Complex double ***x, _Complex double ***xi, _Complex double ****U, int **nn, int T, int N, double tol, int maxIt);
