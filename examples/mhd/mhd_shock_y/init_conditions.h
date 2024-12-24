#include <stdlib.h>
#include <math.h>


extern void initial_conditions(){
	//~ // Initial Conditions //
	for (int i=0; i<nx; i++){
		for (int j=0;j<ny;j++){
				rho[j][i] = 1.0;
				if (j < ny*0.5){
					v[j][i] = 3.0 * sqrt(gamma*1.0/1.0);
				}
					else {
						v[j][i] = -3.0 * sqrt(gamma*1.0/1.0);
					}
				u[j][i] = 0.0;
				p[j][i] = 1.0;
				Bx[j][i] = 1.0;
				By[j][i] = 0.0;
		}
	}

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

