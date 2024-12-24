#include <stdlib.h>
#include <math.h>


extern void initial_conditions(){
	double phi = 5.0/180.0 * pi; // Deflection angle
	// Initial Conditions //
		for (int i=0; i<nx; i++){
		for (int j=0;j<ny;j++){
				rho[j][i] = 1.0;
				u[j][i] = 3.0 * sqrt(gamma*1.0/1.0) * cos(phi);
				v[j][i] = 0.0 * sqrt(gamma*1.0/1.0) * sin(phi);
				v[j][i] = -0.0 * sqrt(gamma*1.0/1.0) * sin(phi);
				p[j][i] = 1.0;
				Bx[j][i] = 0.0;
				By[j][i] = 1.0;
				
			
		}
	}


	// Create Internal Boundary - Wedge
	for (int j=0; j<ny; j++){
		for (int i=0; i<nx; i++){
			if ((i*dx-0.5)*(i*dx-0.5) + (j*dx-0.5)*(j*dx-0.5)<0.15*0.15){ 
				internal_bound[j][i] = 1.0;
			}
				else {
				internal_bound[j][i] = 0.0;
				}
		}
	}
}

