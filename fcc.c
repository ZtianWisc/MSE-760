#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define tstep 2000
#define N_dim 7
#define N_atoms 1372 // total 4000 atoms

time_t t;// I don't know what it means but it works for random number seed.
double r[N_atoms][3]; // position array
double v[N_atoms][3]; // velocity array
double f[N_atoms][3]; // force array
double b[850]; //bin array
double ini_T = 240; // temp in reduced unit
double k_b = 0.000086173303;
double m = 1.0; // mass
double a = 1.68239;
double epsilon = 0.0104;
double sigma = 3.4e-10;
double delta_t = 0.005;
double V2 = 0.0; //velocity square
double rxi, ryi, rzi, rxij, ryij, rzij, rijsq, msd, vacf, LJ, LJold, LJnew;
int MCsteps = 500000;

void initpos(){ //position function
	double boxL = a * N_dim;
	double r_c = boxL/2.0;

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
	
	a = pow(4.0/0.8442,1.0/3.0); // to mimic liquid
	double boxL = a * N_dim;
	double r_c = boxL/2.0;
	double diff_coeff;
	double rxu[N_atoms], ryu[N_atoms], rzu[N_atoms];//store unfold positions
	double r0[N_atoms][3];// store initial position
	double v0[N_atoms][3];// store initial velocity
	for(int i = 0; i < N_atoms; i++){
		f[i][0] = 0.0;
		f[i][1] = 0.0;
		f[i][2] = 0.0;
		r0[i][0] = rxu[i] = r[i][0];
		r0[i][1] = ryu[i] = r[i][1];
		r0[i][2] = rzu[i] = r[i][2];
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
			}
			if (rzij > boxL/2.0){
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
	int c[N_atoms][3];// unfold coordinates count
	
	for (int i=0; i < N_atoms; i++){
		c[i][0]=c[i][1]=c[i][2]=0;
	}
	//time cycle
	for(int t=0; t<tstep; t++){
		double dx2=0.0, dy2=0.0, dz2=0.0;
		msd = 0.0;
		V2 = 0.0;

		for(int i = 0; i < N_atoms; i++){	
		//	update velocity v(t+delta_t/2)
			v[i][0] = v[i][0] + f[i][0]*delta_t/2.0;
			v[i][1] = v[i][1] + f[i][1]*delta_t/2.0;
			v[i][2] = v[i][2] + f[i][2]*delta_t/2.0;
			V2 = V2 + v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2];
			if(t == 100){
				v0[i][0] = v[i][0];
				v0[i][1] = v[i][1];
				v0[i][2] = v[i][2];
		// calculate vacf
			}
			if( t >= 100){
				vacf = vacf + (v[i][0]*v0[i][0]+v[i][1]*v0[i][1]+v[i][2]*v0[i][2])*delta_t;
			}
		// update position r(t+delta_t)
			r[i][0] = r[i][0] + v[i][0]*delta_t;
			if (r[i][0] < 0.0){
				r[i][0] += boxL;
				c[i][0]--;
			}
			if (r[i][0] >= boxL){
				r[i][0] -= boxL;
				c[i][0]++;
			}
			rxu[i] = r[i][0] + c[i][0]*boxL;
			
			r[i][1] = r[i][1] + v[i][1]*delta_t;
			if (r[i][1] < 0.0){
				r[i][1] += boxL;
				c[i][1]--;
			}
			if (r[i][1] >= boxL){
				r[i][1] -= boxL;
				c[i][1]++;
			}
			ryu[i] = r[i][1] + c[i][1]*boxL;
		
			r[i][2] = r[i][2] + v[i][2]*delta_t;
			if (r[i][2] < 0.0){
				r[i][2] += boxL;
				c[i][2]--;
			}
			if (r[i][2] >= boxL){
				r[i][2] -= boxL;
				c[i][2]++;
			}
			rzu[i] = r[i][2] + c[i][2]*boxL;
			
			dx2 = dx2 + pow(rxu[i] - r0[i][0],2);
			dy2 = dy2 + pow(ryu[i] - r0[i][1],2);
			dz2 = dz2 + pow(rzu[i] - r0[i][2],2);
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
				}
			}
		}
		
		//update velocity v(t+delta_t) = v(t+0.5delta_t) + 0.5a(t+delta_t)*delta_t
		
		/*for(int i = 0; i < N_atoms; i++){
			v[i][0] = v[i][0] + f[i][0]*delta_t/2.0;
			v[i][1] = v[i][1] + f[i][1]*delta_t/2.0;
			v[i][2] = v[i][2] + f[i][2]*delta_t/2.0;
		}
		msd = (dx2 + dy2 + dz2)/N_atoms;
		diff_coeff = msd/(6.0*delta_t*t);
		printf("%f, %f\n", t*delta_t, msd);
		if(t > 100){
		printf("%f, %f\n", t*delta_t, vacf/N_atoms/3.0);
		}*/
		printf("%f\n", V2/3.0/N_atoms);
	}
}

void gr(){
	a = pow(4.0/0.8442,1.0/3.0); // to mimic liquid
	double boxL = a * N_dim;
	double r_c = boxL/2.0;
	double ideal = 0.85; // ideal gas number density
	int bin;
	double dr = 0.01;//bin size
	int n_b = r_c/dr;
	
	//set up position
	for(int i = 0; i < N_atoms - 1; i++){
			rxi = r[i][0];
			ryi = r[i][1];
			rzi = r[i][2];
			
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
			if (rijsq < r_c*r_c){
				bin = (int)(sqrt(rijsq)/dr);
				b[bin] = b[bin] + 2;
			}
		}
	}
	for(int i = 0; i < n_b; i++){
		double V_b = (4.0/3.0)* M_PI * (pow(i+1,3)-pow(i,3))*pow(dr,3);
		double normalization = V_b * ideal * N_atoms;
		b[i] = b[i]/normalization;
		printf("%f, %f\n", i*dr, b[i]);
	}
}

void MC(){
	int naccept = 0;
	double boxL = a * N_dim;
	double deltmax = 0.1;
	double V = pow(a*N_dim,3);
	double P = 1.677661;
	LJ = -102.12894;
	srand((unsigned) time(&t));
	rand();
	for(int k = 0; k < MCsteps; k++){
		LJold = 0;
		LJnew = 0;
		double Pold = 0;
		double Pnew = 0;
		//rezero force
		for(int i = 0; i < N_atoms; i++){
		f[i][0] = 0.0;
		f[i][1] = 0.0;
		f[i][2] = 0.0;
		}
		int pick = 1372*((double)rand() / (double)RAND_MAX);
			
		for(int j = 0; j < N_atoms; j++){
			if(j != pick){
			rxij = r[j][0]-r[pick][0]; 
			ryij = r[j][1]-r[pick][1];
			rzij = r[j][2]-r[pick][2];
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
			LJold = LJold + (pow(rijsq,-6)-pow(rijsq,-3));
			f[pick][0] = f[pick][0] - 48/rijsq*(pow(rijsq,-6)-0.5 * pow(rijsq,-3))*rxij;
			f[pick][1] = f[pick][1] - 48/rijsq*(pow(rijsq,-6)-0.5 * pow(rijsq,-3))*ryij;
			f[pick][2] = f[pick][2] - 48/rijsq*(pow(rijsq,-6)-0.5 * pow(rijsq,-3))*rzij;
			}
		}
		Pold = Pold + (rxij*f[pick][0]+ryij*f[pick][1]+rzij*f[pick][2])/3.0/V; 
		
		//displace a random particle
		double x = (2.0 * ((double)rand() / (double)RAND_MAX) - 1.0)*deltmax;
		double y = (2.0 * ((double)rand() / (double)RAND_MAX) - 1.0)*deltmax;
		double z = (2.0 * ((double)rand() / (double)RAND_MAX) - 1.0)*deltmax;
		r[pick][0] = r[pick][0] + x;
		r[pick][1] = r[pick][1] + y;
		r[pick][2] = r[pick][2] + z;
		if (r[pick][0] < 0.0){
			r[pick][0] += boxL;
		}
		if (r[pick][0] >= boxL){
			r[pick][0] -= boxL;
		}
		if (r[pick][1] < 0.0){
			r[pick][1] += boxL;
		}
		if (r[pick][1] >= boxL){
			r[pick][1] -= boxL;
		}
		if (r[pick][2] < 0.0){
			r[pick][2] += boxL;
		}
		if (r[pick][2] >= boxL){
			r[pick][2] -= boxL;
		}
		//rezero force
		for(int i = 0; i < N_atoms; i++){
		f[i][0] = 0.0;
		f[i][1] = 0.0;
		f[i][2] = 0.0;
		}
			
		for(int j = 0;  j < N_atoms; j++){
			if(j != pick){
			rxij = r[j][0]-r[pick][0]; 
			ryij = r[j][1]-r[pick][1];
			rzij = r[j][2]-r[pick][2];
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
			LJnew = LJnew + (pow(rijsq,-6)-pow(rijsq,-3));
			f[pick][0] = f[pick][0] - 48/rijsq*(pow(rijsq,-6)-0.5 * pow(rijsq,-3))*rxij;
			f[pick][1] = f[pick][1] - 48/rijsq*(pow(rijsq,-6)-0.5 * pow(rijsq,-3))*ryij;
			f[pick][2] = f[pick][2] - 48/rijsq*(pow(rijsq,-6)-0.5 * pow(rijsq,-3))*rzij;
			}
		}
		Pnew = Pnew + (rxij*f[pick][0]+ryij*f[pick][1]+rzij*f[pick][2])/3.0/V;
		
		double acc_rate = exp(4.0*(LJold-LJnew)*epsilon/k_b/ini_T);
		if ((double)rand()/(double)RAND_MAX < acc_rate){
			naccept += 1;
			LJ = LJ + (LJnew - LJold)*4.0*epsilon;
			P = P + Pnew - Pold;
		}
		else {
			LJ = LJ;
			P = P;
			r[pick][0] = r[pick][0] - x;
			r[pick][1] = r[pick][1] - y;
			r[pick][2] = r[pick][2] - z;
			if (r[pick][0] < 0.0){
				r[pick][0] += boxL;
			}
			if (r[pick][0] >= boxL){
				r[pick][0] -= boxL;
			}
			if (r[pick][1] < 0.0){
				r[pick][1] += boxL;
			}
			if (r[pick][1] >= boxL){
				r[pick][1] -= boxL;
			}
			if (r[pick][2] < 0.0){
				r[pick][2] += boxL;
			}
			if (r[pick][2] >= boxL){
				r[pick][2] -= boxL;
			}
		}
	printf("%f\n", LJ);
	//printf("%d, %d, %f\n", k, naccept, P);
	}
}

int main(){
	initpos();
	MC();
}
	
