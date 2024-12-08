#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

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
	
	//~ // Initial Conditions // - Kelvin Helmholtz
	//~ for (int i=0; i<nx; i++){
		//~ for (int j=0;j<ny;j++){
				//~ if (j<=ny/4 || j>=ny*3/4){ 
					//~ u[j][i] = 0.25;
					//~ rho[j][i] = 1.0;
					//~ p[j][i] = 1.0;
					//~ v[j][i] = 0.0;
				//~ }
				//~ else {
					//~ u[j][i] = -0.25;
					//~ rho[j][i] = 0.5;
					//~ p[j][i] = 1.0;
					//~ v[j][i] = 0.0;
				//~ }
				//~ if (j==ny/4 || j==ny*3/4){
					//~ v[j][i] = 0.2 * sin(2 * 3.141 / (nx/4) * i);
				//~ }
			
		//~ }
	//~ }
	
	//~ // Initial Conditions //
	for (int i=0; i<nx; i++){
		for (int j=0;j<ny;j++){
				rho[j][i] = 1.0;
				u[j][i] = 2.0 * sqrt(gamma*1.0/1.0);
				v[j][i] = 0.0;
				p[j][i] = 1.0;
				Bx[j][i] = 0.0;
				By[j][i] = 0.0;
		}
	}
	
	
	
	// Create Internal Boundary - Circle
	for (int j=0; j<ny; j++){
		for (int i=0; i<nx; i++){
			if ((i*dx-0.5)*(i*dx-0.5) + (j*dx-0.5)*(j*dx-0.5)<0.15*0.15){ // Circle
				internal_bound[j][i] = 1.0;
			}
				else {
				internal_bound[j][i] = 0.0;
				}
		}
	}
	//~ // Create Internal Boundary - Wedge
	//~ for (int j=0; j<ny; j++){
		//~ for (int i=0; i<nx; i++){
			//~ if ((j*dx-1.0)-tan(-pi/12.0)*(i*dx-0.25)>0){ 
				//~ internal_bound[j][i] = 0.0;
			//~ }
				//~ else {
				//~ internal_bound[j][i] = 0.0;
				//~ }
		//~ }
	//~ }
}
