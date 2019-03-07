#ifndef READ_CONFIG_H
#define READ_CONFIG_H

/* read configuration from file */
_Complex double ****readConfiguration(char path[100], int iconf, _Complex double ****U, int *Lsize, int THREAD_NUM);
_Complex double ****readConfElena(FILE *f, _Complex double ****U, int *Lsize);
_Complex double ****readConfChristian(FILE *f, _Complex double ****U, int *Lsize);
/* set configuration to 1 */
_Complex double ****freeTheory(_Complex double ****U, int N);

#endif
