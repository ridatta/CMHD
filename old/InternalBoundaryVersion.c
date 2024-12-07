#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "input.c"

int main(){
	//~ double rho[ny][nx], u[ny][nx], v[ny][nx], p[ny][nx]; // Primitive variables
	//~ double U[4][ny][nx], Up[4][ny][nx]; // conserved quantities
	//~ double Fx[4][ny][nx], Fxp[4][ny][nx], Fy[4][ny][nx], Fyp[4][ny][nx]; // Fluxes
	//~ double Cs[ny][nx];
	//~ int internal_bound[ny][nx]; // Internal boundary location
	
	// Allocate Memory //
	
	 // Allocate memory for 2D arrays
    double **rho, **u, **v, **p;
    int **internal_bound;

    // Function to allocate memory for 2D arrays
    rho = (double **)malloc(ny * sizeof(double *));
    u = (double **)malloc(ny * sizeof(double *));
    v = (double **)malloc(ny * sizeof(double *));
    p = (double **)malloc(ny * sizeof(double *));
    internal_bound = (int **)malloc(ny * sizeof(int *));

    for (int i = 0; i < ny; i++) {
        rho[i] = (double *)malloc(nx * sizeof(double));
        u[i] = (double *)malloc(nx * sizeof(double));
        v[i] = (double *)malloc(nx * sizeof(double));
        p[i] = (double *)malloc(nx * sizeof(double));
        internal_bound[i] = (int *)malloc(nx * sizeof(int));
    }

    // Allocate memory for 3D arrays (4 x ny x nx)
    double ***U, ***Up, ***Fx, ***Fxp, ***Fy, ***Fyp;

    // Function to allocate memory for 3D arrays
    U = (double ***)malloc(4 * sizeof(double **));
    Up = (double ***)malloc(4 * sizeof(double **));
    Fx = (double ***)malloc(4 * sizeof(double **));
    Fxp = (double ***)malloc(4 * sizeof(double **));
    Fy = (double ***)malloc(4 * sizeof(double **));
    Fyp = (double ***)malloc(4 * sizeof(double **));

    for (int idx = 0; idx < 4; idx++) {
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
    
    // Declare Variables //
	
	int cntr = 0;
	double dx = 1.0/(double)nx; // Grid space
	double nu = visc_fac * (1.0 - CFL*CFL)/6.0; // Artificial viscosity
	
	
	// Initial Conditions //
	for (int i=0; i<nx; i++){
		for (int j=0;j<ny;j++){
				rho[j][i] = 1.0;
				u[j][i] = Ms * sqrt(gamma*1.0/1.0);
				v[j][i] = 0.0;
				p[j][i] = 1.0;
				if ((j*dx-1.0)-tan(-pi/4.0)*(i*dx-0.25)>0){
					u[j][i] = 0.5;
					rho[j][i] = 0.5;
					p[j][i] = 0.5;
				}
			
		}
	}
	
	
	// Calculate Conserved Variables
	for (int i=0; i<nx; i++){
		for (int j=0; j<ny; j++){
			U[0][j][i] = 1.0 * rho[j][i];
			U[1][j][i] =  rho[j][i] * u[j][i];
			U[2][j][i] =  rho[j][i] * v[j][i];
			U[3][j][i] =  p[j][i] / (gamma-1) + 0.5 * rho[j][i] * (u[j][i]*u[j][i] + v[j][i]*v[j][i]);
		}
	}
	
		for (int idx=0; idx < 4; idx++){ // Fill predicted Fluxes 
		for (int i=0; i < nx; i++){
			for (int j=0; j<ny; j++){
				Fxp[idx][j][i] = 1.0;
				Fyp[idx][j][i] = 1.0;
				Up[idx][j][i] =  1.0;
			}
		}
	}
	
	// Create Internal Boundary
	for (int j=0; j<ny; j++){
		for (int i=0; i<nx; i++){
			if ((j*dx-1.0)-tan(-pi/12)*(i*dx-0.25)>0 && i*dx<0.6){
				internal_bound[j][i] = 1;
			}
			else if (i*dx>=0.6 && j*dx > 1.0+tan(-pi/12)*(0.6-0.25)){
				internal_bound[j][i] = 1;
			}
				else {
				internal_bound[j][i] = 0;
				}
		}
	}
	
	//~ // Create Internal Boundary
	//~ for (int j=0; j<ny; j++){
		//~ for (int i=0; i<nx; i++){
			//~ if ((j*dx-1.0)-tan(-1.0*def_angle)*(i*dx-0.25)>0){
				//~ internal_bound[j][i] = 1;
			//~ }
				//~ else {
				//~ internal_bound[j][i] = 0;
				//~ }
		//~ }
	//~ }
	

	
	// Write output
	FILE *fpt;
	fpt = fopen("bound.txt", "w+");
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            fprintf(fpt,"%d ",internal_bound[j][i]);
        }
        fprintf(fpt,"\n");
    }
    fclose(fpt);
    fpt = fopen("rho-0.txt", "w+");
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            fprintf(fpt,"%f ",rho[j][i]);
        }
        fprintf(fpt,"\n");
    }
    fclose(fpt);
    
	
	
	// Main Time Loop //
	
	// Get Max. Velocity
	double maxV = 0.0;
	for (int i=0; i<nx; i++){
		for (int j=0; j<ny; j++){
			if (sqrt(gamma * p[j][i] / rho[j][i]) > maxV){
				maxV = sqrt(gamma * p[j][i] / rho[j][i]);
			}
			if (sqrt(u[j][i]*u[j][i] + v[j][i]*v[j][i]) > maxV){
				maxV = sqrt(u[j][i]*u[j][i] + v[j][i]*v[j][i]);
			}
		}
	}
	
	double dt = CFL * dx / maxV; // Get dt
	
	printf("Running McMormick Scheme\n\n");
	
	double t = 0.0;
	
	while (t < tend){
		// Calculate Fluxes
		for (int i=0; i < nx; i++){
			for (int j=0; j<ny; j++){
				// x-fluxes
				Fx[0][j][i] = rho[j][i] * u[j][i];
				Fx[1][j][i] = rho[j][i] * u[j][i] * u[j][i] + p[j][i];
				Fx[2][j][i] = rho[j][i] * v[j][i] * u[j][i];
				Fx[3][j][i] = u[j][i] * (p[j][i] / (gamma-1.0) + p[j][i] +  0.5 * rho[j][i] * (u[j][i]*u[j][i] + v[j][i]*v[j][i]));
				// y-fluxes
				Fy[0][j][i] = rho[j][i] * v[j][i];
				Fy[1][j][i] = rho[j][i] * u[j][i] * v[j][i];
				Fy[2][j][i] = rho[j][i] * v[j][i] * v[j][i] + p[j][i];
				Fy[3][j][i] = v[j][i] * (p[j][i] / (gamma-1.0) + p[j][i] +  0.5 * rho[j][i] * (u[j][i]*u[j][i] + v[j][i]*v[j][i]));
			}
		}
		
		double temp = 0.0;
		
		// Predictor Step //
		for (int idx=0; idx<4; idx++){
			for (int j=1; j<(ny-1); j++){
				for (int i=1; i<(nx-1); i++){
					if (internal_bound[j][i] == 1){
						Up[1][j][i] = 0.0; // If internal boundary set velocity to 0
						Up[2][j][i] = 0.0;
						
					}
					else {
						Up[idx][j][i] = (U[idx][j][i] 
						- dt/dx * (Fx[idx][j][i+1]-Fx[idx][j][i])
						- dt/dx * (Fy[idx][j+1][i]-Fy[idx][j][i])
						+ nu * (U[idx][j][i+1] + U[idx][j][i-1] - 2*U[idx][j][i])
						+ nu * (U[idx][j+1][i] + U[idx][j-1][i] - 2*U[idx][j][i]));
					}
				}
			}
		}
		
		// Boundary Conditions //
		for (int idx=0; idx<4; idx++){
			for (int j=0; j<ny; j++){
				Up[idx][j][0] = Up[idx][j][1]; // Left
				Up[idx][j][nx-1] = Up[idx][j][nx-2]; // Right
			}
			for (int i=0; i<nx; i++){
				Up[idx][0][i] = Up[idx][1][i]; // Bottom
				Up[idx][ny-1][i] = Up[idx][ny-2][i]; // Left
			}
		}
		
		// Recalculate Primitive Variables //
		for (int j = 0; j < ny; j++){
			for (int i = 0; i < nx; i++){
				rho[j][i] = 1.0 * Up[0][j][i];
				u[j][i] = Up[1][j][i] / Up[0][j][i];
				v[j][i] = Up[2][j][i] / Up[0][j][i];
				p[j][i] = (Up[3][j][i] - 0.5 * rho[j][i] * (u[j][i]*u[j][i]+v[j][i]*v[j][i])) * (gamma - 1.0);
			}
		}
		
		// Calculate Predicted Fluxes
		for (int i=0; i < nx; i++){
			for (int j=0; j<ny; j++){
				// x-fluxes
				Fxp[0][j][i] = rho[j][i] * u[j][i];
				Fxp[1][j][i] = rho[j][i] * u[j][i] * u[j][i] + p[j][i];
				Fxp[2][j][i] = rho[j][i] * v[j][i] * u[j][i];
				Fxp[3][j][i] = u[j][i] * (p[j][i] / (gamma-1.0) + p[j][i] +  0.5 * rho[j][i] * (u[j][i]*u[j][i] + v[j][i]*v[j][i]));
				// y-fluxes
				Fyp[0][j][i] = rho[j][i] * v[j][i];
				Fyp[1][j][i] = rho[j][i] * u[j][i] * v[j][i];
				Fyp[2][j][i] = rho[j][i] * v[j][i] * v[j][i] + p[j][i];
				Fyp[3][j][i] = v[j][i] * (p[j][i] / (gamma-1.0) + p[j][i] +  0.5 * rho[j][i] * (u[j][i]*u[j][i] + v[j][i]*v[j][i]));
			}
		}
		
		
		// Corrector Step //
		for (int idx=0; idx<4; idx++){
			for (int j=1; j<(ny-1); j++){
				for (int i=1; i<(nx-1); i++){
					if (internal_bound[j][i] == 1){
						U[1][j][i] = 0.0;
						U[2][j][i] = 0.0;
					}
					else{
						U[idx][j][i] = (0.5*(U[idx][j][i]+Up[idx][j][i])
						- 0.5*dt/dx * (Fxp[idx][j][i]-Fxp[idx][j][i-1])
						- 0.5*dt/dx * (Fyp[idx][j][i]-Fyp[idx][j-1][i])
						+nu*(Up[idx][j][i+1] + Up[idx][j][i-1] -2*Up[idx][j][i])
						+nu*(Up[idx][j+1][i] + Up[idx][j-1][i] -2*Up[idx][j][i]));
					}
				}
			}
		}
		
		// Boundary Conditions //
		for (int idx=0; idx<4; idx++){
			for (int j=0; j<ny; j++){
				U[idx][j][0] = U[idx][j][1]; // Left
				U[idx][j][nx-1] = U[idx][j][nx-2]; // Right
			}
			for (int i=0; i<nx; i++){
				U[idx][0][i] = U[idx][1][i]; // Bottom
				U[idx][ny-1][i] = U[idx][ny-2][i]; // Top
			}
		}
		
		// Recalculate Primitive Variables //
		for (int j = 0; j < ny; j++){
			for (int i = 0; i < nx; i++){
				rho[j][i] = 1.0 * U[0][j][i];
				u[j][i] = U[1][j][i] / U[0][j][i];
				v[j][i] = U[2][j][i] / U[0][j][i];
				p[j][i] = (U[3][j][i] - 0.5 * rho[j][i] * (u[j][i]*u[j][i]+v[j][i]*v[j][i])) * (gamma - 1.0);
			}
		}
		
		// Calculate and Update Time //
		

		// Get Max. Velocity
		maxV = 0.0;
		for (int i=0; i<nx; i++){
			for (int j=0; j<ny; j++){
				if (sqrt(gamma * p[j][i] / rho[j][i]) > maxV){
					maxV = sqrt(gamma * p[j][i] / rho[j][i]);
				}
				if (sqrt(u[j][i]*u[j][i] + v[j][i]*v[j][i]) > maxV){
					maxV = sqrt(u[j][i]*u[j][i] + v[j][i]*v[j][i]);
				}
			}
		}
		
		
		dt = CFL * dx / maxV; // Get dt
		t = t + dt;
		cntr++; // Update counter
		
		if (cntr % 25 == 0){
			printf("t = %f, step = %d , maxV = %f, dt = %f\n",t,cntr,maxV,dt);
		}
		
	}
	
	// Print //
	printf("Simulation ended at step %d at time %f \n",cntr,t);
	
	
	// Write output
	//~ FILE *fpt;
	fpt = fopen("rho.txt", "w+");
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            fprintf(fpt,"%f ",rho[j][i]);
        }
        fprintf(fpt,"\n");
    }
    fclose(fpt);
    
    fpt = fopen("p.txt", "w+");
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            fprintf(fpt,"%f ",p[j][i]);
        }
        fprintf(fpt,"\n");
    }
    fclose(fpt);
    
    fpt = fopen("u.txt", "w+");
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            fprintf(fpt,"%f ",u[j][i]);
        }
        fprintf(fpt,"\n");
    }
    fclose(fpt);
    
    fpt = fopen("v.txt", "w+");
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            fprintf(fpt,"%f ",v[j][i]);
        }
        fprintf(fpt,"\n");
    }
    fclose(fpt);
    
    // Free the allocated memory //
    for (int i = 0; i < ny; i++) {
        free(rho[i]);
        free(u[i]);
        free(v[i]);
        free(p[i]);
        free(internal_bound[i]);
    }

    free(rho);
    free(u);
    free(v);
    free(p);
    free(internal_bound);

    for (int idx = 0; idx < 4; idx++) {
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
    
    
    return 0;
}
	
