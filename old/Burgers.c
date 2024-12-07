#include <stdio.h>
#include <stdlib.h>
#define nx 4096
void main(){
	double *x = malloc(sizeof(double) * nx);
	double *u = malloc(sizeof(double) * nx);
    double *up = malloc(sizeof(double) * nx);
    double *f = malloc(sizeof(double) * nx);
    double *fp = malloc(sizeof(double) * nx);
    double *temp;
	
	double dx = 1.0/(double)nx;
	// initial condition
	for (int i=0; i<nx; i++){
		if (i < nx / 2) {
			x[i] = i * dx;
			u[i] = 1.0;
			f[i] = 1.0 * u[i];
		}
			else {
				x[i] = i * dx;
				u[i] = 0.1;
				f[i] = 1.0 * u[i];
			}
	}
	
	double dt = 1.0e-4;
	double max_step = 0.4/dt;
	
	double C = 1.0 * dt / dx; // CFL No.
	double nu = 0.5 * (1-C*C)/6; // Artifical viscosity
	
	printf("Max step is %f\n", max_step);
	
	for (int num=0; num<(int)max_step; num++){
        
		// Predictor step
		for (int ii = 0; ii<nx; ii++){
			up[ii] = u[ii] - dt/dx * (f[ii+1] - f[ii]) + nu * (u[ii+1] + u[ii-1] - 2* u[ii]);
			fp[ii] = 1.0 * up[ii]; // Predicted fluxes
		}
        	
		// Corrector step
		for (int ii = 0; ii<nx; ii++){
			u[ii] = 0.5 * (u[ii] + up[ii]) - 0.5 * dt/dx * (fp[ii] - fp[ii-1]) + nu * (up[ii+1] + up[ii-1] - 2* up[ii]);
			f[ii] = 1.0 * u[ii]; // Update flux
		}
		
		// ANti-diffusion step
		for (int ii = 0; ii<nx; ii++){
			u[ii] = u[ii] -  nu * (up[ii+1] + up[ii-1] - 2* up[ii]);
			f[ii] = 1.0 * u[ii]; // Update flux
		}
	}
	
	
	// Write output
	FILE *fpt;
	fpt = fopen("output.csv", "w+");
	for(int i = 0; i < nx; i++){
     fprintf(fpt, "%f,%f\n", x[i],u[i]);
	}
	fclose(fpt);
}
