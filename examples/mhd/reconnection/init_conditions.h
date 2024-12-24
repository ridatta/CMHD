#include <stdlib.h>
#include <math.h>


extern void initial_conditions(){
	//~ // Initial Conditions //
	for (int i=0; i<nx; i++){
		for (int j=0;j<ny;j++){
				rho[j][i] = 1.0;
				u[j][i] = 0.1 * sin(2*pi*j*dx/0.5);
				v[j][i] = 0.0;
				p[j][i] = 0.1;
				Bx[j][i] = 0.0;
				if (i*dx <0.25){
					By[j][i] = 1.0;
				}
				else if (i*dx <0.75 & i*dx >= 0.25){
					By[j][i] = -1.0;
				}
				else {
					By[j][i] = 1.0;
				}
				
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

