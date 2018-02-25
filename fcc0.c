#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

time_t t;// I don't know what it means but it works for random number seed.
double r[1372][3]; //position array
double v[1372][3]; //velocity array
double f[1372][3]; //force array
int N_dim = 7; // 7 cells per dimension
int N_atoms = 1372; // total 1372 atoms
double ini_T = 2.5; // temp in reduced unit
double k_b = 1.0; // Bolzmann constant
double a = 5.26/3.4, m = 1.0; // cell length and mass
int tstep = 500;
double delta_t = 0.02;
double V2 = 0.0; //velocity square

void initpos(){ //position function

// initialize positions
	for(int i = 0; i < N_atoms; i++){
		r[i][0] = 0.0;
		r[i][1] = 0.0;
		r[i][2] = 0.0;
	}
	int n = 0;
	for (int x = 0; x < N_dim; x++){
		for (int y = 0; y < N_dim; y++){
			for (int z = 0; z < N_dim; z++){
				for (int p = 0; p < 4; p++){
					if (n < N_atoms && p%4 == 0) {
						r[n][0] = x * a;
						r[n][1] = y * a;
						r[n][2] = z * a;
						n = n + 1;
					}
					else if (n < N_atoms && p%4 == 1){
						r[n][0] = (x + 1.0/2.0)*a;
						r[n][1] = (y + 1.0/2.0)*a;
						r[n][2] = z * a;
						n = n + 1;
					}
					else if (n < N_atoms && p%4 == 2){
						r[n][0] = x * a;
						r[n][1] = (y + 1.0/2.0)*a;
						r[n][2] = (z + 1.0/2.0)*a;
						n = n + 1;
					}
					else if (n < N_atoms && p%4 == 3){
						r[n][0] = (x + 1.0/2.0)*a;
						r[n][1] = y * a;
						r[n][2] = (z + 1.0/2.0)*a;
						n = n + 1;
					}
				}
			}
		}
	}
}

void initvel(){ //initialize velocity vi = (r0,r1,r2)v0
	double v0 = sqrt(3.0 * ini_T);
	double r0, r1, r2, s2;
	for(int i = 0; i < N_atoms; i++){
		v[i][0] = 0.0;
		v[i][1] = 0.0;
		v[i][2] = 0.0;
	}
	srand((unsigned) time(&t));
	rand();
	for(int n=0; n < N_atoms; n++){
		do {
			r0 = 2.0 * ((double)rand() / (double)RAND_MAX) - 1.0;
			r1 = 2.0 * ((double)rand() / (double)RAND_MAX) - 1.0;
			s2 = r0*r0 + r1*r1; //s2 is actually s^2 as in book MDBasics
				}while (s2 > 1.0);
				
				r0 = 2.0*sqrt(1.0-s2)*r0;
				r1 = 2.0*sqrt(1.0-s2)*r1;
				r2 = 1 - 2.0*s2;
				//assign velocities to all atoms
				v[n][0] = v0 * r0;
				v[n][1] = v0 * r1;
				v[n][2] = v0 * r2;
				//calculate V total and V average
				V2 = V2 + v[n][0]*v[n][0]+v[n][1]*v[n][1]+v[n][2]*v[n][2];
	}
}

void MD_wi_pbc(){
	
	//calculate LJ and force. f_x = -dV/d|x| * x/|x| = 48/|r|**2 (r**-12 - 0.5r**-6) x in atomic units
	// 1 amu f = 8.2387225e-8 Newton
	
	double epsilon = 0.0104;
	double sigma = 3.4e-10;
	double boxL = a * N_dim;
	double r_c = boxL/2.0;
	double rxi, ryi, rzi, rxij, ryij, rzij, rijsq;
	for(int i = 0; i < N_atoms; i++){
		f[i][0] = 0.0;
		f[i][1] = 0.0;
		f[i][2] = 0.0;
	}
	//states at t=0
	for(int i = 0; i < N_atoms - 1; i++){
		rxi = r[i][0];
		ryi = r[i][1];
		rzi = r[i][2];
		//calculate force and LJ	
		for(int j = i + 1; j < N_atoms; j++){
			rxij = r[j][0]-rxi; 
			ryij = r[j][1]-ryi;
			rzij = r[j][2]-rzi;
			if (rxij > boxL/2.0){
				rxij -= boxL;
			}
			if (rxij < boxL/(-2.0)){
				rxij += boxL;
			}
			if (ryij > boxL/2.0){
				ryij -= boxL;
			}
			if (ryij < boxL/(-2.0)){
				ryij += boxL;
			}if (rzij > boxL/2.0){
				rzij -= boxL;
			}
			if (rzij < boxL/(-2.0)){
				rzij += boxL;
			}
			rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
			if(rijsq < r_c*r_c){
				// find force and acceleration 
				f[i][0] = f[i][0] - 48/(rijsq)*(pow(rijsq,-6)-0.5 * pow(rijsq,-3))*rxij;
				f[j][0] = f[j][0] + 48/(rijsq)*(pow(rijsq,-6)-0.5 * pow(rijsq,-3))*rxij;
				f[i][1] = f[i][1] - 48/(rijsq)*(pow(rijsq,-6)-0.5 * pow(rijsq,-3))*ryij;
				f[j][1] = f[j][1] + 48/(rijsq)*(pow(rijsq,-6)-0.5 * pow(rijsq,-3))*ryij;
				f[i][2] = f[i][2] - 48/(rijsq)*(pow(rijsq,-6)-0.5 * pow(rijsq,-3))*rzij;
				f[j][2] = f[j][2] + 48/(rijsq)*(pow(rijsq,-6)-0.5 * pow(rijsq,-3))*rzij;
			}
		}
	}
	//time steps, start simulation
	for(int t=0; t<tstep; t++){
		double LJ = 0.0;
		double K = 0.0;

		for(int i = 0; i < N_atoms; i++){	
		//update velocity v(t+delta_t/2)
			v[i][0] = v[i][0] + f[i][0]*delta_t/2.0;
			v[i][1] = v[i][1] + f[i][1]*delta_t/2.0;
			v[i][2] = v[i][2] + f[i][2]*delta_t/2.0;
			// K = K + (v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2])*0.5;
		//update position r(t+delta_t)
			r[i][0] = r[i][0] + v[i][0]*delta_t;
			if (r[i][0] < 0.0){
				r[i][0] += boxL;
			}
			if (r[i][0] >= boxL){
				r[i][0] -= boxL;
			}
			r[i][1] = r[i][1] + v[i][1]*delta_t;
			if (r[i][1] < 0.0){
				r[i][1] += boxL;
			}
			if (r[i][1] >= boxL){
				r[i][1] -= boxL;
			}
			r[i][2] = r[i][2] + v[i][2]*delta_t;
			if (r[i][2] < 0.0){
				r[i][2] += boxL;
			}
			if (r[i][2] >= boxL){
				r[i][2] -= boxL;
			}
		}
		//rezero force
		for(int i = 0; i < N_atoms; i++){
			f[i][0] -= f[i][0];
			f[i][1] -= f[i][1];
			f[i][2] -= f[i][2];
		}

		//update force a(t+delta_t)
		for(int i = 0; i < N_atoms - 1; i++){
			rxi = r[i][0];
			ryi = r[i][1];
			rzi = r[i][2];
			
			//calculate force and LJ	
			for(int j = i + 1; j < N_atoms; j++){
				rxij = r[j][0]-rxi; 
				ryij = r[j][1]-ryi;
				rzij = r[j][2]-rzi;
				if (rxij > boxL/2.0){
				rxij -= boxL;
				}
				if (rxij < boxL/(-2.0)){
				rxij += boxL;
				}
				if (ryij > boxL/2.0){
					ryij -= boxL;
				}
				if (ryij < boxL/(-2.0)){
					ryij += boxL;
				}
				if (rzij > boxL/2.0){
				rzij -= boxL;
				}
				if (rzij < boxL/(-2.0)){
					rzij += boxL;
				}
				rijsq = rxij*rxij + ryij*ryij + rzij*rzij;
				if(rijsq < r_c*r_c){
				//calculate force and acceleration 
					f[i][0] = f[i][0] - 48/rijsq*(pow(rijsq,-6)-0.5 * pow(rijsq,-3))*rxij;
					f[j][0] = f[j][0] + 48/rijsq*(pow(rijsq,-6)-0.5 * pow(rijsq,-3))*rxij;
					f[i][1] = f[i][1] - 48/rijsq*(pow(rijsq,-6)-0.5 * pow(rijsq,-3))*ryij;
					f[j][1] = f[j][1] + 48/rijsq*(pow(rijsq,-6)-0.5 * pow(rijsq,-3))*ryij;
					f[i][2] = f[i][2] - 48/rijsq*(pow(rijsq,-6)-0.5 * pow(rijsq,-3))*rzij;
					f[j][2] = f[j][2] + 48/rijsq*(pow(rijsq,-6)-0.5 * pow(rijsq,-3))*rzij;
					// LJ = LJ + pow(rijsq,-6)-pow(rijsq,-3);
				}
			}
		}
		
		//update velocity v(t+delta_t) = v(t+0.5delta_t) + 0.5a(t+delta_t)*delta_t
		for(int i = 0; i < N_atoms; i++){
			v[i][0] = v[i][0] + f[i][0]*delta_t/2.0;
			v[i][1] = v[i][1] + f[i][1]*delta_t/2.0;
			v[i][2] = v[i][2] + f[i][2]*delta_t/2.0;
		}
		// LJ = 4.0 * LJ;
	}	
}

int main(){
	initpos();
	initvel();
	MD_wi_pbc();
}
	
