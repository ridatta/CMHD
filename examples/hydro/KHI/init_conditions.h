#include <stdlib.h>
#include <math.h>


extern void initial_conditions(){
	//~ // Initial Conditions //
	for (int i=0; i<nx; i++){
		for (int j=0;j<ny;j++){
				rho[j][i] = 1.0;
				if (j*dx < 0.5){
					u[j][i] = 0.25 * sqrt(gamma*1.0/1.0);
					rho[j][i] = 0.1;
				}
					else {
						u[j][i] = -0.25 * sqrt(gamma*1.0/1.0);
						rho[j][i] = 1.0;
					}
				v[j][i] = 0.05 * sin(2*pi/0.25*i*dx);
				p[j][i] = 1.0;
				Bx[j][i] = 0.0;
				By[j][i] = 0.0 * tanh((i*dx-0.5)/0.25);
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
