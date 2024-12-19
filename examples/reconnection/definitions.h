// Definitions.h
#include "input.h"
#define mu0 1.0
#define nvar 6

double *S;
double **rho, **u, **v, **p, **Bx, **By;
double **internal_bound;
double ***U, ***Up, ***Fx, ***Fxp, ***Fy, ***Fyp;
double dx = 1.0/(double)nx; // Grid space

// Memory Allocation

extern void allocateArrays(){
	
	// Allocate Memory //
	S = (double *)malloc(nvar * sizeof(double *));

    // Function to allocate memory for 2D arrays
    rho = (double **)malloc(ny * sizeof(double *));
    u = (double **)malloc(ny * sizeof(double *));
    v = (double **)malloc(ny * sizeof(double *));
    p = (double **)malloc(ny * sizeof(double *));
    Bx = (double **)malloc(ny * sizeof(double *));
    By = (double **)malloc(ny * sizeof(double *));
    internal_bound = (double **)malloc(ny * sizeof(double *));

    for (int i = 0; i < ny; i++) {
        rho[i] = (double *)malloc(nx * sizeof(double));
        u[i] = (double *)malloc(nx * sizeof(double));
        v[i] = (double *)malloc(nx * sizeof(double));
        p[i] = (double *)malloc(nx * sizeof(double));
        Bx[i] = (double *)malloc(nx * sizeof(double));
        By[i] = (double *)malloc(nx * sizeof(double));
        internal_bound[i] = (double *)malloc(nx * sizeof(double));
    }

    // Allocate memory for 3D arrays (6 x ny x nx)

    // Function to allocate memory for 3D arrays
    U = (double ***)malloc(nvar * sizeof(double **));
    Up = (double ***)malloc(nvar * sizeof(double **));
    Fx = (double ***)malloc(nvar * sizeof(double **));
    Fxp = (double ***)malloc(nvar * sizeof(double **));
    Fy = (double ***)malloc(nvar * sizeof(double **));
    Fyp = (double ***)malloc(nvar * sizeof(double **));

    for (int idx = 0; idx < nvar; idx++) {
        U[idx] = (double **)malloc(ny * sizeof(double *));
        Up[idx] = (double **)malloc(ny * sizeof(double *));
        Fx[idx] = (double **)malloc(ny * sizeof(double *));
        Fxp[idx] = (double **)malloc(ny * sizeof(double *));
        Fy[idx] = (double **)malloc(ny * sizeof(double *));
        Fyp[idx] = (double **)malloc(ny * sizeof(double *));

        for (int i = 0; i < ny; i++) {
            U[idx][i] = (double *)malloc(nx * sizeof(double));
            Up[idx][i] = (double *)malloc(nx * sizeof(double));
            Fx[idx][i] = (double *)malloc(nx * sizeof(double));
            Fxp[idx][i] = (double *)malloc(nx * sizeof(double));
            Fy[idx][i] = (double *)malloc(nx * sizeof(double));
            Fyp[idx][i] = (double *)malloc(nx * sizeof(double));
        }
    }
    
}

extern void freeMemory(){
	// Free the allocated memory //
    for (int i = 0; i < ny; i++) {
        free(rho[i]);
        free(u[i]);
        free(v[i]);
        free(p[i]);
        free(Bx[i]);
        free(By[i]);
        free(internal_bound[i]);
    }

    free(rho);
    free(u);
    free(v);
    free(p);
    free(Bx);
    free(By);
    free(internal_bound);

    for (int idx = 0; idx < nvar; idx++) {
        for (int i = 0; i < ny; i++) {
            free(U[idx][i]);
            free(Up[idx][i]);
            free(Fx[idx][i]);
            free(Fxp[idx][i]);
            free(Fy[idx][i]);
            free(Fyp[idx][i]);
        }

        free(U[idx]);
        free(Up[idx]);
        free(Fx[idx]);
        free(Fxp[idx]);
        free(Fy[idx]);
        free(Fyp[idx]);
    }

    free(U);
    free(Up);
    free(Fx);
    free(Fxp);
    free(Fy);
    free(Fyp);
}
