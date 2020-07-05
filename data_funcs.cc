#include <math.h>
#include <vector>
#include "global_variables.cc"
using namespace std;

//Volume-weighted average
double calc_vol_avg(double f[Nx][Ny][Nz])  {
    
	int i,j,k;
	double tot = 0;

	for (i = 0; i < Nx; i++)  {
		for (j = 0; j < Ny; j++)  {
			for (k = 0; k < Nx; k++)  {
				tot += f[i][j][k];
			}
		}
	}
	return tot / Nx / Ny / Nz;
}

//Mass-weighted average
double calc_mass_avg(double f[Nx][Ny][Nz], double den[Nx][Ny][Nz])  {

	int i,j,k;
	double tot    = 0;
	double m_tot  = 0;

	for (i = 0; i < Nx; i++)  {
		for (j = 0; j < Ny; j++)  {
			for (k = 0; k < Nx; k++)  {
				m_tot += den[i][j][k];
				tot   += f[i][j][k]*den[i][j][k];
			}
		}
    }
    return tot/m_tot;
}

int update_treion(void)  {
	for (int i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nx; k++)  {
				if ( (nH1[i][j][k]/nH[i][j][k] < 0.5) && (treion[i][j][k] == 0) )  {
					treion[i][j][k] = t;
				}
			}
		}
	}
}

