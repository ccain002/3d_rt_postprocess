//healpix index functions
#include <stdio.h>
#include <omp.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "general_funcs.cc"
using namespace std;

#define Ntheta 3
#define Nphi   4

struct pairint  {
	int i1, i2;
};
	
struct paird  {
	double d1, d2;
};

struct tripled  {
	double d1, d2, d3;
};

int f_row(int f)  {
	int frow = floor(f / Nphi);
	return frow;
}

int F1(int f)  {
	return f_row(f) + 2;
}

int F2(int f)  {
	return 2*mod(f, Nphi) - mod(f_row(f), 2) + 1;
}

int get_bit(int n, int i)  {
	return floor(log2(mod(n, pow(2, i+1))));
}

int bits_to_num(int bits[], int i)  {
	int num = 0;
	int j;
	
	for (j = 0; j < i; j++)  {
		num += pow(2*bits[j], j);
	}
	return num;
}

auto get_ij(int pn, int Nside)  {
	int x = 0, y = 0, vv = 0, h2 = 0;
	int m;
	int pnp = mod(pn, pow(Nside, 2));
	int n = max(0, floor(log2(pnp)));
	int i, j;
	int imax = (Ntheta+1)*Nside - 1;

	for (m = 0; m <= n; m++)  {
		int bm = floor(mod(pnp, pow(2, m+1))/pow(2, m));
		if (mod(m, 2) == 0)  {
			x += pow(2,(m/2))*bm;
		}
		else  {
			y += pow(2, (m - 1)/2)*bm;
		}
	}
	
	int f = floor(pn/pow(Nside, 2));
	
	vv = x + y;
	h2 = x - y;

	i = F1(f)*Nside - vv - 1;
	if (i < Nside)  {
		j = (F2(f)*i + h2 + 1)/2;
	}
	else if (i > (Ntheta + 1)*Nside - Nside - 1)  {
		j = (F2(f)*(imax + 1 - i) + h2 + 1)/2;
	}
	else  {
		j = (F2(f)*Nside + h2 + mod(i - Nside + 1, 2))/2;
	}
	
	//testing
	if (j == 0)  {
		j = 4*Nside;
	}  
	
	return pairint {i, j};
}

auto get_zphi(int i, int j, int pn, int Nside)  {
	
	double z, phi;
	int imax = (Ntheta+1)*Nside - 1;
	int f = floor(pn/pow(Nside, 2));
	
	//northern polar cap
	if (i < Nside)  {
		z   = 1. - (1./3.)*pow(i/((double) Nside), 2);
		phi = pi/2./i*(j - 1./2.);
		return paird {z, phi};
	}
	//northern equatorial belt
	else if ( (i >= Nside) && (i <= 2*Nside) )  {
		z   = 4./3. - 2.*i/3./Nside;
		phi = pi/2./Nside*(j - mod(i - Nside + 1, 2)/2.);
		return paird {z, phi};
	}
	//southern equatorial belt
	else if ( (i > 2*Nside) && (i <= 3*Nside) )  {
		z   = -4./3. + 2.*(imax + 1 - i)/3./Nside;
		phi = pi/2./Nside*(j - mod(imax + 1 - i - Nside + 1, 2)/2.);
		return paird {z, phi};
	}
	//southern hemisphere
	else {
		z   = -1. + (1./3.)*pow((imax + 1. - i)/((double) Nside), 2);
		phi = pi/2./(imax + 1 - i)*(j - 1./2.);
		return paird {z, phi};
	}
}

auto get_unit_vector(double z, double phi, double alpha0, double beta0, double gamma0, int direction)  {
	
	double xhat, yhat, zhat;
	double xhatp, yhatp, zhatp;
	double theta = acos(z);
	
	xhat = sin(theta)*cos(phi);
	yhat = sin(theta)*sin(phi);
	zhat = cos(theta);
	
	if (direction == 0)  {
		xhatp = (cos(beta0)*cos(gamma0))*xhat 
			  + (-cos(beta0)*sin(gamma0))*yhat
			  + (sin(beta0))*zhat;
		yhatp = (sin(alpha0)*sin(beta0)*cos(gamma0) + cos(alpha0)*sin(gamma0))*xhat 
			  + (-sin(alpha0)*sin(beta0)*sin(gamma0) + cos(alpha0)*cos(gamma0))*yhat
			  + (-sin(alpha0)*cos(beta0))*zhat;
		zhatp = (-cos(alpha0)*sin(beta0)*cos(gamma0) + sin(alpha0)*sin(gamma0))*xhat
			  + (cos(alpha0)*sin(beta0)*sin(gamma0) + sin(alpha0)*cos(gamma0))*yhat
			  + (cos(alpha0)*cos(beta0))*zhat;
	}
	else  {
		xhatp = (cos(gamma0)*cos(beta0))*xhat
		      + (sin(gamma0)*cos(alpha0) + cos(gamma0)*sin(beta0)*sin(alpha0))*yhat
			  + (sin(gamma0)*sin(alpha0) - cos(gamma0)*sin(beta0)*cos(alpha0))*zhat;
		yhatp = (-sin(gamma0)*cos(beta0))*xhat
			  + (cos(gamma0)*cos(alpha0) - sin(gamma0)*sin(beta0)*sin(alpha0))*yhat
			  + (cos(gamma0)*sin(alpha0) + sin(gamma0)*sin(beta0)*cos(alpha0))*zhat;
		zhatp = (sin(beta0))*xhat
			  + (-cos(beta0)*sin(alpha0))*yhat
			  + (cos(beta0)*cos(alpha0))*zhat;
	}
	return tripled {xhatp, yhatp, zhatp};
}

