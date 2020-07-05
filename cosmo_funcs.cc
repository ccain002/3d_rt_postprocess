#include <math.h>
#include "healpix_funcs.cc"

//Hubble parameter at scale factor a
double H(double a)  {
	return hh*100./mpc_to_km*pow(Omega_m/pow(a,3) + Omega_l,0.5);
}

//Critical density of the universe
double rho_crit(double a)  {
	return 3*pow(H(a),2)/8/pi/G_grav;
}

double cosmic_time(double a)  {
	int n = 100;
	int i;
	double t; 
	double sf[n];
	double integrand[n];
	
	for (i = 0; i < n; i++)  {
		sf[i]     = i / (n - 1) * a;
		integrand[i] = 1/sf[i]/H(sf[i]);
	}
	
	t = trapz_int(integrand, sf, n);
	
	return t;
	
}