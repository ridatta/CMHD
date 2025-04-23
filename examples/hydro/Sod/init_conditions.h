#include <stdlib.h>
#include <math.h>

double f(double x, double y){ // Define the internal boundary 
	//~ return pow(y-0.5,2) + pow(x-0.5,2) - 0.15*0.15;
	//~ return fabs(y-0.5) - tan(pi/12.0) * (x-0.25);
	//~ return fabs(y+0.05) - tan(10.0/180.0*pi) * (x - 0.1);
	return 1.0;
}

void getNormalVector(double x, double y, double arr[2]){ // Retruns the normal vector at an internal boundary
	double eps = 1.0 * dx;
	double fx = (f(x+eps,y)-f(x-eps,y))/ (2 * eps);
	double fy = (f(x,y+eps)-f(x,y-eps))/ (2 * eps);
	double fp = sqrt(fx*fx + fy*fy);
	arr[0] = fx/fp;
	arr[1] = fy/fp;
}


extern void initial_conditions(){
	//~ // Initial Conditions //
	for (int i=0; i<nx; i++){
		for (int j=0;j<ny;j++){
				rho[j][i] = 1.0;
				if (i < nx*0.5){
					rho[j][i] = 1.0;
					p[j][i] = 1.0;
				}
					else {
						rho[j][i] = 0.125;
						p[j][i] = 0.1;
					}
				v[j][i] = 0.0;
				u[j][i] = 0.0;
				Bx[j][i] = 0.0;
				By[j][i] = 0.0;
		}
	}

	// Create Internal Boundary - Wedge
	for (int j=0; j<ny; j++){
		for (int i=0; i<nx; i++){
			if (f(i*dx,j*dx)<=0){ 
				internal_bound[j][i] = 1.0;
				double nvec[2];
				getNormalVector(i*dx, j*dx, nvec);
				bnx[j][i] = 1.0* nvec[0];
				bny[j][i] = 1.0* nvec[1];
			}
				else {
				internal_bound[j][i] = 0.0;
				bnx[j][i] = 0.0;
				bny[j][i] = 0.0;
				}
		}
	}
}

