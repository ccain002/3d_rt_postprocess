//General functions for coding various physical processes not related to the rt algorithm itself
#include <cstdlib>
#include <math.h>
#include <algorithm>
#include <random>
#include <time.h> 
#include "global_constants.h"
using namespace std; 

int min(int a, int b)  {
	if (a >= b)  {
		return b;
	}
	else  {
		return a;
	}
}

int max(int a, int b)  {
	if (b >= a)  {
		return b;
	}
	else  {
		return a;
	}
}

//Absolute value with double argument
double absd(double x)  {

	if (x >= 0)  {
		return x;
	}
	else  {
		return -x;
	}
}

double mind(double a, double b)  {
	if (a >= b)  {
		return b;
	}
	else  {
		return a;
	}
}

double maxd(double a, double b)  {
	if (b >= a)  {
		return b;
	}
	else  {
		return a;
	}
}

double init_rand()  {
	srand(time(0));
}

double fRand(double a, double b)  {
	
   int v1 = rand() % 1001 + 0;
   double f1 = ((double) v1)/1000.;
   return a + f1*(b - a);
}

int mod(int a, int b) {
	
	if (a < 0) {
		a += b;
	}
	
	int m = a % b;
	
	return m;
}

double remainder(double a, double b)  {
	return a - floor(a/b)*b;
}

//Blackbody radiation function.  
double b_nu(double nu, double T)
{
	return 2*h*pow(nu,3.)/pow(cl,2.)*(1./(exp(h*nu/k_B/T) - 1.));
} 

//power law density function to test the rt algorithm.  
//Use alpha = 0 for constant density plane parallel case.  
double power_law(double r, double r0, double A, double alpha)
{
	if (r != 0)  {
		return A*pow(r/r0, alpha);
	}
	else  {
		return A;
	}
}

double trapz_int(double y[], double x[], int n)  {

	int i;
	double F  = 0;
	double dF = 0;

	for (i = 0; i < n - 1; i++)  {
		dF  = y[i]*(x[i + 1] - x[i]);
		dF += y[i+1]*(x[i + 1] - x[i]);
		F  += dF/2.;
	}
	return F;
}

double cum_trapz_int(double y[], double x[], int n)  {
	
	int i,j;
	double cumint[n];
	cumint[0] = 0.;
	for (i = 1; i < n; i++)  {
		double temp_x[i + 1];
		double temp_y[i + 1];
		for (j = 0; j < i + 1; j++)  {
			temp_x[j] = x[j];
			temp_y[j] = y[j];
		}
		cumint[i] = trapz_int(temp_y, temp_x, i + 1);
	}
	return *cumint;
}

double interpolate(double x[], double y[], double x0, int n)  {
    
	int i;
	double y0 = 0.;
    
	for (i = 0; i < n - 1; i++)  {
		if ( ( (x[i] <= x0) && (x[i+1] >= x0) ) || ( (x[i] >= x0) && (x[i+1] <= x0) ) )  {
			double slope = (y[i+1] - y[i])/(x[i+1] - x[i]);
			y0 = y[i] + slope*(x0 - x[i]);
			break;
		}
		else  {
			double slope = (y[i+1] - y[i])/(x[i+1] - x[i]);
			y0 = y[i] + slope*(x0 - x[i]);
		}
	}
	return y0;
}
