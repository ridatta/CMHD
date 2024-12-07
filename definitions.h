// Definitions.h
#include "input.h"

double *S;
double **rho, **u, **v, **p;
double **internal_bound;
double ***U, ***Up, ***Fx, ***Fxp, ***Fy, ***Fyp;
double dx = 1.0/(double)nx; // Grid space

