#include <stdlib.h>
#include <math.h>


extern void initial_conditions(){
	// Initial Conditions //
	double phi = 20.0/180.0 * pi;
		for (int i=0; i<nx; i++){
		for (int j=0;j<ny;j++){
				rho[j][i] = 1.0;
				u[j][i] = 1.5 * sqrt(gamma*1.0/1.0) * sin(phi);
				if (j*dx > 0.5){
					v[j][i] = -1.5 * sqrt(gamma*1.0/1.0) * cos(phi);
				}
				else {
					v[j][i] = 1.5 * sqrt(gamma*1.0/1.0) * cos(phi);
				}
				p[j][i] = 1.0;
				Bx[j][i] = 0.0;
				By[j][i] = 0.0;
				if ( fabs(j*dx-0.5) -tan(-pi/4.0)*(i*dx)>0){ // Oblique shock
					v[j][i] = 0.0;
				}
				
			
		}
	}


	// Create Internal Boundary - Wedge
	for (int j=0; j<ny; j++){
		for (int i=0; i<nx; i++){
			if ((j*dx-1.0)-tan(-pi/8.0)*(i*dx-0.25)>0){ // Oblique shock
				internal_bound[j][i] = 0.0;
			}
				else {
				internal_bound[j][i] = 0.0;
				}
		}
	}
}

