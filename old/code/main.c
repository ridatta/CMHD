#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "input.h"
#include "definitions.h"
#include <string.h>


void allocateArrays(){
	
	// Allocate Memory //
	S = (double *)malloc(4 * sizeof(double *));

    // Function to allocate memory for 2D arrays
    rho = (double **)malloc(ny * sizeof(double *));
    u = (double **)malloc(ny * sizeof(double *));
    v = (double **)malloc(ny * sizeof(double *));
    p = (double **)malloc(ny * sizeof(double *));
    internal_bound = (double **)malloc(ny * sizeof(double *));

    for (int i = 0; i < ny; i++) {
        rho[i] = (double *)malloc(nx * sizeof(double));
        u[i] = (double *)malloc(nx * sizeof(double));
        v[i] = (double *)malloc(nx * sizeof(double));
        p[i] = (double *)malloc(nx * sizeof(double));
        internal_bound[i] = (double *)malloc(nx * sizeof(double));
    }

    // Allocate memory for 3D arrays (4 x ny x nx)

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
    
}

void freeMemory(){
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
}

void initial_conditions(){
	// Initial Conditions //
	//~ for (int i=0; i<nx; i++){
		//~ for (int j=0;j<ny;j++){
				//~ rho[j][i] = 1.0;
				//~ u[j][i] = 2.5 * sqrt(gamma*1.0/1.0);
				//~ v[j][i] = 0.0;
				//~ p[j][i] = 1.0;
				//~ if ((j*dx-1.0)-tan(-pi/4.0)*(i*dx-0.25)>0){ // Oblique shock
					//~ u[j][i] = 0.5;
					//~ rho[j][i] = 0.5;
					//~ p[j][i] = 0.5;
				//~ }
			
		//~ }
	//~ }
	
	//~ // Initial Conditions // - RMI
	//~ for (int i=0; i<nx; i++){
		//~ for (int j=0;j<ny;j++){
				//~ if ((j*dx-1.25)>5e-2*cos(i*dx*2*pi/1.0)){ 
					//~ u[j][i] = 0.0;
					//~ rho[j][i] = 1.0;
					//~ p[j][i] = 0.1;
					//~ v[j][i] = -0.0;
				//~ }
				//~ else {
					//~ u[j][i] = 0.0;
					//~ rho[j][i] = 0.1;
					//~ p[j][i] = 0.1;
					//~ v[j][i] =-0.0;
				//~ }
				//~ if (j*dx<1.1){ 
					//~ u[j][i] = 0.0;
					//~ rho[j][i] = 0.25;
					//~ p[j][i] = 0.25;
					//~ v[j][i] = 2.0;
				//~ }		
		//~ }
	//~ }
	
	//~ // Initial Conditions // - RT
	//~ for (int i=0; i<nx; i++){
		//~ for (int j=0;j<ny;j++){
				//~ if (j>ny/2){ 
					//~ u[j][i] = 0.0;
					//~ rho[j][i] = 1.0;
					//~ p[j][i] = 20.0 + 1.0 * g * (j*dx);
					//~ v[j][i] = 0.0;
				//~ }
				//~ else {
					//~ u[j][i] = 0.0;
					//~ rho[j][i] = 0.1;
					//~ p[j][i] = 20.0 + 1.0 * g * (ny/2*dx) + 0.1 *  g * (j*dx) + 0.1 * -g * (ny/2*dx);
					//~ v[j][i] = 0.0;
				//~ }	
				//~ if (j==ny/2){
					//~ v[j][i] = 0.2 * sin(2 * 3.141 / (nx/2) * i);
				//~ }
		//~ }
	//~ }
	
	// Initial Conditions // - Kelvin Helmholtz
	for (int i=0; i<nx; i++){
		for (int j=0;j<ny;j++){
				if (j<=ny/4 || j>=ny*3/4){ 
					u[j][i] = 0.25;
					rho[j][i] = 1.0;
					p[j][i] = 1.0;
					v[j][i] = 0.0;
				}
				else {
					u[j][i] = -0.25;
					rho[j][i] = 0.5;
					p[j][i] = 1.0;
					v[j][i] = 0.0;
				}
				if (j==ny/4 || j==ny*3/4){
					v[j][i] = 0.2 * sin(2 * 3.141 / (nx/4) * i);
				}
			
		}
	}
	
	// Initial Conditions //
	//~ for (int i=0; i<nx; i++){
		//~ for (int j=0;j<ny;j++){
				//~ rho[j][i] = 1.0;
				//~ u[j][i] = 0.0 * sqrt(gamma*1.0/1.0);
				//~ v[j][i] = 0.0;
				//~ p[j][i] = 1.0;
		//~ }
	//~ }
	
	
	
	//~ // Create Internal Boundary - Circle
	//~ for (int j=0; j<ny; j++){
		//~ for (int i=0; i<nx; i++){
			//~ if ((i*dx-0.5)*(i*dx-0.5) + (j*dx-0.5)*(j*dx-0.5)<0.15*0.15){ // Circle
				//~ internal_bound[j][i] = 1.0;
			//~ }
				//~ else {
				//~ internal_bound[j][i] = 0.0;
				//~ }
		//~ }
	//~ }
	// Create Internal Boundary - Wedge
	for (int j=0; j<ny; j++){
		for (int i=0; i<nx; i++){
			if ((j*dx-1.0)-tan(-pi/12.0)*(i*dx-0.25)>0){ 
				internal_bound[j][i] = 0.0;
			}
				else {
				internal_bound[j][i] = 0.0;
				}
		}
	}
}


void writeOutput(double **arr, char *fname){
	// Write output
	FILE *fpt;
	fpt = fopen(fname, "w+");
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            fprintf(fpt,"%f ",arr[j][i]);
        }
        fprintf(fpt,"\n");
    }
    fclose(fpt);
}

double getMaxVelocity(double **u,double **v, double **p, double **rho){
	
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
	return maxV;
}


double getQrad(double density,double pressure){
	if (irad == 1){
		double A = -1.0;
		double alpha = 1.0, beta = 0.5;
		double T;
		T = pressure * 1.0 / density; // Temp 
		return A * density * pow(density,alpha) + pow(T,beta);
	}
	else {
		return 0.0;
	}
}

int main(){
	
	double nu = visc_fac * (1.0 - CFL*CFL)/6.0; // Artificial viscosity
	int cntr = 0;
	// Allocate memory 
	allocateArrays();
	// Initial conditions and internal boundary
	initial_conditions();

	
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
		
		// Source Terms
		for (int idx = 0; idx < 4; idx++){
			S[idx] = 0.0; 
		}
	
	// Write output
	writeOutput(internal_bound,"./output/bound-0.txt");
	writeOutput(rho,"./output/rho-0.txt");
	writeOutput(p,"./output/p-0.txt");
	writeOutput(u,"./output/u-0.txt");
	writeOutput(v,"./output/v-0.txt");
	
  
	// Main Time Loop //
	double maxV = getMaxVelocity(u,v,p,rho);// Get Max. Velocity	
	double dt = CFL * dx /  maxV;// Get dt
	double t = 0.0;
	
	
	printf("Running McMormick Scheme\n\n");
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
					if (internal_bound[j][i] == 1.0){
						Up[1][j][i] = 0.0; // If internal boundary set velocity to 0
						Up[2][j][i] = 0.0;
						
					}
					else {
						S[2] = 0.0; // Source term for y-momentum
						S[3] = getQrad(rho[j][i],p[j][i]); // Source term for energy eqn.
						
						Up[idx][j][i] = (U[idx][j][i] 
						- dt/dx * (Fx[idx][j][i+1]-Fx[idx][j][i])
						- dt/dx * (Fy[idx][j+1][i]-Fy[idx][j][i])
						+ nu * (U[idx][j][i+1] + U[idx][j][i-1] - 2*U[idx][j][i])
						+ nu * (U[idx][j+1][i] + U[idx][j-1][i] - 2*U[idx][j][i])
						+ S[idx]*dt);
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
				Up[idx][ny-1][i] = Up[idx][ny-2][i]; // Top
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
					if (internal_bound[j][i] == 1.0){
						U[1][j][i] = 0.0;
						U[2][j][i] = 0.0;
					}
					else{
						S[2] = 0.0; // Source term for y-momentum
						S[3] = getQrad(rho[j][i],p[j][i]); // Recalculate Raditaive cooling with upaded quantities
						
						U[idx][j][i] = (0.5*(U[idx][j][i]+Up[idx][j][i])
						- 0.5*dt/dx * (Fxp[idx][j][i]-Fxp[idx][j][i-1])
						- 0.5*dt/dx * (Fyp[idx][j][i]-Fyp[idx][j-1][i])
						+nu*(Up[idx][j][i+1] + Up[idx][j][i-1] -2*Up[idx][j][i])
						+nu*(Up[idx][j+1][i] + Up[idx][j-1][i] -2*Up[idx][j][i])
						+0.5*S[idx]*dt);
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
		maxV = getMaxVelocity(u,v,p,rho); // Get Max. Velocity
		dt = CFL * dx / maxV; // Get dt
		t = t + dt;
		cntr++; // Update counter
		
		
		//~ // OUTPUT //
		//~ FILE *fptr; // Pressure vs Time
		//~ if (cntr == 1){
			//~ fptr = fopen("./output/p_vs_t.txt", "w+");
			//~ }
		//~ else{
			//~ fptr = fopen("./output/p_vs_t.txt", "a");
		//~ }
		//~ fprintf(fptr,"%f,%f\n",t,p[ny/2][nx/2]);
		//~ fclose(fptr);
		//~ FILE *fptr;
		//~ if (cntr == 1){
			//~ fptr = fopen("./output/v_vs_t.txt", "w+"); // Total Energy vs Time
			//~ }
		//~ else{
			//~ fptr = fopen("./output/v_vs_t.txt", "a");
		//~ }
		//~ fprintf(fptr,"%f,%f\n",t,v[ny/2][nx/2]);
		//~ fclose(fptr);
		
		
		if (cntr % 100 == 0){
			printf("t = %f, step = %d , maxV = %f, dt = %f\n",t,cntr,maxV,dt);
			char fname[50];
			sprintf(fname,"./output/u-%d.txt",cntr);
			writeOutput(u,fname);
			sprintf(fname,"./output/v-%d.txt",cntr);
			writeOutput(v,fname);
			sprintf(fname,"./output/rho-%d.txt",cntr);
			writeOutput(rho,fname);
		}
		
	}
	
	// Print //
	printf("Simulation ended at step %d at time %f \n",cntr,t);
	
	// Write Output
	writeOutput(rho,"./output/rho.txt");
	writeOutput(u,"./output/u.txt");
	writeOutput(v,"./output/v.txt");
	writeOutput(p,"./output/p.txt");    
    
    return 0;
}
	
