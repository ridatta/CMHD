#include <stdlib.h>
#include <math.h>


extern void initial_conditions(){
	//~ // Initial Conditions //
	double dr = 3.5 * dx;
	for (int i=0; i<nx; i++){
		for (int j=0;j<ny;j++){
				rho[j][i] = 1.0;
				v[j][i] = 0.0;
				u[j][i] = 0.0;
				Bx[j][i] = 0.0;
				By[j][i] = 0.0;
				if ((i*dx-0.5)*(i*dx-0.5) + (j*dx-0.33)*(j*dx-0.33) <= dr*dr){
					p[j][i] = 3.0 * (gamma-1)*1.0/ ((2+1) * 3.141 * pow(dr,2));
				}
				else if ((i*dx-0.5)*(i*dx-0.5) + (j*dx-0.67)*(j*dx-0.67) <= dr*dr){
					p[j][i] = 3.0 * (gamma-1)*1.0/ ((2+1) * 3.141 * pow(dr,2));
				}
					else {
						p[j][i] = 1e-5;
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

