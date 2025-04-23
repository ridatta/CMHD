#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "input.h"
#include "definitions.h"
#include "init_conditions.h"
#include <string.h>


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
		double A = -2.0;
		double alpha = 2.0, beta = 0.5;
		double T;
		T = pressure * 1.0 / density; // Temp 
		return A * density * pow(density,alpha) * pow(T,beta);
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
			U[3][j][i] =  p[j][i] / (gamma-1) + 0.5 * rho[j][i] * (u[j][i]*u[j][i] + v[j][i]*v[j][i]) + 0.5 * 1.0/mu0 * (Bx[j][i]*Bx[j][i] + By[j][i]*By[j][i]);
			U[4][j][i] = 1.0 * Bx[j][i];
			U[5][j][i] = 1.0 * By[j][i];
		}
	}
	
		for (int idx=0; idx < nvar; idx++){ // Fill predicted Fluxes 
			for (int i=0; i < nx; i++){
				for (int j=0; j<ny; j++){
					Fxp[idx][j][i] = 1.0;
					Fyp[idx][j][i] = 1.0;
					Up[idx][j][i] =  1.0;
				}
			}
		}
		
		// Source Terms
		for (int idx = 0; idx < nvar; idx++){
			S[idx] = 0.0; 
		}
	
	// Write output
	writeOutput(internal_bound,"./output/bound-0.txt");
	writeOutput(bnx,"./output/bnx-0.txt");
	writeOutput(bny,"./output/bny-0.txt");
	writeOutput(rho,"./output/rho-0.txt");
	writeOutput(p,"./output/p-0.txt");
	writeOutput(u,"./output/u-0.txt");
	writeOutput(v,"./output/v-0.txt");
	writeOutput(Bx,"./output/Bx-0.txt");
	writeOutput(By,"./output/By-0.txt");
	
  
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
				Fx[1][j][i] = rho[j][i] * u[j][i] * u[j][i] + p[j][i] + 0.5*(Bx[j][i]*Bx[j][i] + By[j][i]*By[j][i])/mu0 - Bx[j][i]*Bx[j][i]/mu0;
				Fx[2][j][i] = rho[j][i] * v[j][i] * u[j][i] - 1.0/mu0 * By[j][i] * Bx[j][i];
				Fx[3][j][i] = u[j][i] * (p[j][i] / (gamma-1.0) + p[j][i] +  0.5 * rho[j][i] * (u[j][i]*u[j][i] + v[j][i]*v[j][i]) + (Bx[j][i]*Bx[j][i] + By[j][i]*By[j][i])/mu0) - Bx[j][i] * (u[j][i]*Bx[j][i] + v[j][i] * By[j][i])/mu0;
				Fx[4][j][i] = 0.0;
				Fx[5][j][i] = u[j][i] * By[j][i] - Bx[j][i] * v[j][i];
				// y-fluxes
				Fy[0][j][i] = rho[j][i] * v[j][i];
				Fy[1][j][i] = rho[j][i] * u[j][i] * v[j][i] - 1.0/mu0 * Bx[j][i] * By[j][i];
				Fy[2][j][i] = rho[j][i] * v[j][i] * v[j][i] + p[j][i] + 0.5*(Bx[j][i]*Bx[j][i] + By[j][i]*By[j][i])/mu0 - By[j][i]*By[j][i]/mu0;;
				Fy[3][j][i] = v[j][i] * (p[j][i] / (gamma-1.0) + p[j][i] +  0.5 * rho[j][i] * (u[j][i]*u[j][i] + v[j][i]*v[j][i]) + (Bx[j][i]*Bx[j][i] + By[j][i]*By[j][i])/mu0) - By[j][i] * (u[j][i]*Bx[j][i] + v[j][i] * By[j][i])/mu0;
				Fy[4][j][i] = -1* u[j][i] * By[j][i] + Bx[j][i] * v[j][i];
				Fy[5][j][i] = 0.0;
			}
		}
		
		double temp = 0.0;
		double vn = 0.0;
		double ub = 0.0;
		double vb = 0.0;
		int locj = 1;
		int loci = 1;
		float thrs = 0.1;
		
		
		// Predictor Step //
		for (int idx=0; idx<nvar; idx++){
			for (int j=1; j<(ny-1); j++){
				for (int i=1; i<(nx-1); i++){
					if (internal_bound[j][i] == 1.0){
						
							if (fabs(fabs(bny[j][i]) - fabs(bnx[j][i])) <= thrs){
								locj = bny[j][i] / fabs(bny[j][i]);
								loci = bnx[j][i] / fabs(bnx[j][i]);
							}
							else if (fabs(bny[j][i])>fabs(bnx[j][i])) {
								locj = bny[j][i] / fabs(bny[j][i]);
								loci = 0;	
							}
							else {
								locj = 0;
								loci = bnx[j][i] / fabs(bnx[j][i]);	
							}
							
							Up[0][j][i] = Up[0][j+locj][i+loci];
							Up[3][j][i] = Up[3][j+locj][i+loci];
							vn = (Up[1][j][i]/Up[0][j][i]) * bnx[j][i] + (Up[2][j][i]/Up[0][j][i]) * bny[j][i];
							ub = (Up[1][j][i]/Up[0][j][i]) - 1.0 * vn * bnx[j][i];
							vb = (Up[2][j][i]/Up[0][j][i]) - 1.0 * vn * bny[j][i];
							Up[1][j][i] = Up[0][j][i] * ub * 1.0; // Velocity BC
							Up[2][j][i] = Up[0][j][i] * vb * 1.0; 
							Up[4][j][i] = 0.0;
							Up[5][j][i] = 0.0;
						
						
					}
					else {
						S[2] = 0.0; // Source term for y-momentum
						S[3] = getQrad(rho[j][i],p[j][i]); // Source term for energy eqn.
						
						Up[idx][j][i] = (U[idx][j][i] 
						- 0.5*dt/dx * (Fx[idx][j][i+1]-Fx[idx][j][i-1]) // Using CDD was originally FDD
						- 0.5*dt/dx * (Fy[idx][j+1][i]-Fy[idx][j-1][i])
						+ nu * (U[idx][j][i+1] + U[idx][j][i-1] - 2*U[idx][j][i])
						+ nu * (U[idx][j+1][i] + U[idx][j-1][i] - 2*U[idx][j][i])
						+ S[idx]*dt);
						//~ Up[idx][j][i] = (U[idx][j][i] 
						//~ - 1.0*dt/dx * (Fx[idx][j][i+1]-Fx[idx][j][i]) // FDD
						//~ - 1.0*dt/dx * (Fy[idx][j+1][i]-Fy[idx][j][i])
						//~ + nu * (U[idx][j][i+1] + U[idx][j][i-1] - 2*U[idx][j][i])
						//~ + nu * (U[idx][j+1][i] + U[idx][j-1][i] - 2*U[idx][j][i])
						//~ + S[idx]*dt);
					}
				}
			}
		}
		
		// Boundary Conditions //
		for (int idx=0; idx<nvar; idx++){
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
				Bx[j][i] = 1.0 * Up[4][j][i];
				By[j][i] = 1.0 * Up[5][j][i];
				u[j][i] = Up[1][j][i] / Up[0][j][i];
				v[j][i] = Up[2][j][i] / Up[0][j][i];
				p[j][i] = (Up[3][j][i] - 0.5 * rho[j][i] * (u[j][i]*u[j][i]+v[j][i]*v[j][i]) - 0.5 /mu0 * (Bx[j][i]*Bx[j][i] + By[j][i]*By[j][i])) * (gamma - 1.0);
			}
		}
		
		// Calculate Predicted Fluxes
		for (int i=0; i < nx; i++){
			for (int j=0; j<ny; j++){
				// x-fluxes
				Fxp[0][j][i] = rho[j][i] * u[j][i];
				Fxp[1][j][i] = rho[j][i] * u[j][i] * u[j][i] + p[j][i] + 0.5*(Bx[j][i]*Bx[j][i] + By[j][i]*By[j][i])/mu0 - Bx[j][i]*Bx[j][i]/mu0;
				Fxp[2][j][i] = rho[j][i] * v[j][i] * u[j][i] - 1.0/mu0 * By[j][i] * Bx[j][i];
				Fxp[3][j][i] = u[j][i] * (p[j][i] / (gamma-1.0) + p[j][i] +  0.5 * rho[j][i] * (u[j][i]*u[j][i] + v[j][i]*v[j][i]) + (Bx[j][i]*Bx[j][i] + By[j][i]*By[j][i])/mu0) - Bx[j][i] * (u[j][i]*Bx[j][i] + v[j][i] * By[j][i])/mu0;
				Fxp[4][j][i] = 0.0;
				Fxp[5][j][i] = u[j][i] * By[j][i] - Bx[j][i] * v[j][i];
				// y-fluxes
				Fyp[0][j][i] = rho[j][i] * v[j][i];
				Fyp[1][j][i] = rho[j][i] * u[j][i] * v[j][i] - 1.0/mu0 * Bx[j][i] * By[j][i];
				Fyp[2][j][i] = rho[j][i] * v[j][i] * v[j][i] + p[j][i] + 0.5*(Bx[j][i]*Bx[j][i] + By[j][i]*By[j][i])/mu0 - By[j][i]*By[j][i]/mu0;;
				Fyp[3][j][i] = v[j][i] * (p[j][i] / (gamma-1.0) + p[j][i] +  0.5 * rho[j][i] * (u[j][i]*u[j][i] + v[j][i]*v[j][i]) + (Bx[j][i]*Bx[j][i] + By[j][i]*By[j][i])/mu0) - By[j][i] * (u[j][i]*Bx[j][i] + v[j][i] * By[j][i])/mu0;
				Fyp[4][j][i] = -1* u[j][i] * By[j][i] + Bx[j][i] * v[j][i];
				Fyp[5][j][i] = 0.0;
			}
		}
		
		
		// Corrector Step //
		for (int idx=0; idx<nvar; idx++){
			for (int j=1; j<(ny-1); j++){
				for (int i=1; i<(nx-1); i++){
					if (internal_bound[j][i] == 1.0){
						
							if (fabs(fabs(bny[j][i]) - fabs(bnx[j][i])) <= thrs){
								locj = bny[j][i] / fabs(bny[j][i]);
								loci = bnx[j][i] / fabs(bnx[j][i]);
							}
							else if (fabs(bny[j][i])>fabs(bnx[j][i])) {
								locj = bny[j][i] / fabs(bny[j][i]);
								loci = 0;	
							}
							else {
								locj = 0;
								loci = bnx[j][i] / fabs(bnx[j][i]);	
							}
							
							 //~ printf("locj, loci = %d, %d \n",locj,loci);
						
							U[0][j][i] = U[0][j+locj][i+loci];
							U[3][j][i] = U[3][j+locj][i+loci];
							vn = (U[1][j][i]/U[0][j][i]) * bnx[j][i] + (U[2][j][i]/U[0][j][i]) * bny[j][i];
							ub = (U[1][j][i]/U[0][j][i]) - 1.0 * vn * bnx[j][i]; // No normal velocity at wall
							vb = (U[2][j][i]/U[0][j][i]) - 1.0 * vn * bny[j][i];
							U[1][j][i] = U[0][j][i] * ub * 1.0;  // Velocity BC
							U[2][j][i] = U[0][j][i] * vb * 1.0;  // Velocity BC
							U[4][j][i] = 0.0;
							U[5][j][i] = 0.0;
						
					}
					else{
						S[2] = 0.0; // Source term for y-momentum
						S[3] = getQrad(rho[j][i],p[j][i]); // Recalculate Raditaive cooling with upaded quantities
						
						U[idx][j][i] = (0.5*(U[idx][j][i]+Up[idx][j][i])
						- 0.25*dt/dx * (Fxp[idx][j][i+1]-Fxp[idx][j][i-1]) // Using  CDD was originally BDD
						- 0.25*dt/dx * (Fyp[idx][j+1][i]-Fyp[idx][j-1][i])
						+nu*(Up[idx][j][i+1] + Up[idx][j][i-1] -2*Up[idx][j][i])
						+nu*(Up[idx][j+1][i] + Up[idx][j-1][i] -2*Up[idx][j][i])
						+0.5*S[idx]*dt);
						//~ U[idx][j][i] = (0.5*(U[idx][j][i]+Up[idx][j][i])
						//~ - 0.5*dt/dx * (Fxp[idx][j][i]-Fxp[idx][j][i-1]) //  BDD
						//~ - 0.5*dt/dx * (Fyp[idx][j][i]-Fyp[idx][j-1][i])
						//~ +nu*(Up[idx][j][i+1] + Up[idx][j][i-1] -2*Up[idx][j][i])
						//~ +nu*(Up[idx][j+1][i] + Up[idx][j-1][i] -2*Up[idx][j][i])
						//~ +0.5*S[idx]*dt);
						
					}
				}
			}
		}
		
		// Boundary Conditions //
		for (int idx=0; idx<nvar; idx++){
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
				Bx[j][i] = 1.0 * U[4][j][i];
				By[j][i] = 1.0 * U[5][j][i];
				u[j][i] = U[1][j][i] / U[0][j][i];
				v[j][i] = U[2][j][i] / U[0][j][i];
				p[j][i] = (U[3][j][i] - 0.5 * rho[j][i] * (u[j][i]*u[j][i]+v[j][i]*v[j][i]) - 0.5 /mu0 * (Bx[j][i]*Bx[j][i] + By[j][i]*By[j][i])) * (gamma - 1.0);
			}
		}
		
		// Calculate and Update Time //
		maxV = getMaxVelocity(u,v,p,rho); // Get Max. Velocity
		dt = CFL * dx / maxV; // Get dt
		t = t + dt;
		cntr++; // Update counter
		
		
		// OUTPUT //
		FILE *fptr; // Pressure vs Time
		if (cntr == 1){
			fptr = fopen("./output/p_vs_t.txt", "w+");
			}
		else{
			fptr = fopen("./output/p_vs_t.txt", "a");
		}
		fprintf(fptr,"%f,%f\n",t,p[ny/2][nx/2]);
		fclose(fptr);
		//~ FILE *fptr;
		//~ if (cntr == 1){
			//~ fptr = fopen("./output/v_vs_t.txt", "w+"); // Total Energy vs Time
			//~ }
		//~ else{
			//~ fptr = fopen("./output/v_vs_t.txt", "a");
		//~ }
		//~ fprintf(fptr,"%f,%f\n",t,v[ny/2][nx/2]);
		//~ fclose(fptr);
		
		
		if (cntr % num_print == 0){
			printf("t = %f, step = %d , maxV = %f, dt = %f\n",t,cntr,maxV,dt);
			char fname[50];
			sprintf(fname,"./output/u-%d.txt",cntr);
			writeOutput(u,fname);
			sprintf(fname,"./output/v-%d.txt",cntr);
			writeOutput(v,fname);
			sprintf(fname,"./output/rho-%d.txt",cntr);
			writeOutput(rho,fname);
			sprintf(fname,"./output/p-%d.txt",cntr);
			writeOutput(p,fname);
			sprintf(fname,"./output/Bx-%d.txt",cntr);
			writeOutput(Bx,fname);
			sprintf(fname,"./output/By-%d.txt",cntr);
			writeOutput(By,fname);
		}
		
	}
	
	// Print //
	printf("Simulation ended at step %d at time %f \n",cntr,t);
	
	// Write Output
	writeOutput(rho,"./output/rho.txt");
	writeOutput(u,"./output/u.txt");
	writeOutput(v,"./output/v.txt");
	writeOutput(p,"./output/p.txt");    
	writeOutput(Bx,"./output/Bx.txt");
	writeOutput(By,"./output/By.txt"); 
    
    return 0;
}
	
