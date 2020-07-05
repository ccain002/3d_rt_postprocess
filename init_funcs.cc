//Functions for initializing and updating arrays
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "io_funcs.cc"

int init(void)  {
	int i;
	
	if (input_grid == TRUE)  {
		read_grid_binary();
	}
	else  {
		get_ugas(); //get gas data
		
		//spatial grid   
		for (i = 0; i < Nx; i++)  {
			x[i]  = (i + 0.5)*Lx/Nx*kpc_to_cm;
			dx[i] = Lx/Nx*kpc_to_cm;
		}
		for (i = 0; i < Ny; i++)  {
			y[i]  = (i + 0.5)*Ly/Ny*kpc_to_cm;
			dy[i] = Ly/Ny*kpc_to_cm;
		}
		for (i = 0; i < Nz; i++)  {
			z[i]  = (i + 0.5)*Lz/Nz*kpc_to_cm;
			dz[i] = Lz/Nz*kpc_to_cm;
		}
		#pragma omp parallel
		{
		#pragma omp for 
		for (i = 0; i < Nx; i++)  {
			for (int j = 0; j < Ny; j++)  {
				for (int k = 0; k < Nz; k++)  {
					f_H1[i][j][k]  = fH1_0;
					f_H2[i][j][k]  = 1. - f_H1[i][j][k];
					f_He1[i][j][k] = fHe1_0;
					f_He2[i][j][k] = fHe2_0;
					f_He3[i][j][k] = 1. - f_He1[i][j][k] - f_He2[i][j][k];
					
					temp[i][j][k] = temp_0;
				}
			}
		}
		}
	}
	
	//compute other important quantities
	#pragma omp parallel
	{
	#pragma omp for 
	for (i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				//number densities
				nH[i][j][k]    = (1 - Y)*rho[i][j][k]/m_H;
				nHe[i][j][k]   = Y*rho[i][j][k]/(4*m_H);
				nH1[i][j][k]   = f_H1[i][j][k]*nH[i][j][k];
				nHe1[i][j][k]  = f_He1[i][j][k]*nHe[i][j][k];
				nH2[i][j][k]   = f_H2[i][j][k]*nH[i][j][k];
				nHe2[i][j][k]  = f_He2[i][j][k]*nHe[i][j][k];
				nHe3[i][j][k]  = f_He3[i][j][k]*nHe[i][j][k];
				ne[i][j][k]    = nH2[i][j][k] + nHe2[i][j][k] + 2.*nHe3[i][j][k];
				n_tot[i][j][k] = nH[i][j][k] + nHe[i][j][k] + ne[i][j][k];
				
				nH1_prev[i][j][k]  = nH1[i][j][k];
				nH2_prev[i][j][k]  = nH2[i][j][k];
				nHe1_prev[i][j][k] = nHe1[i][j][k];
				nHe2_prev[i][j][k] = nHe2[i][j][k];
				nHe3_prev[i][j][k] = nHe3[i][j][k];
				
				//on the fly stepping
				f_H1_step[i][j][k]  = f_H1[i][j][k];
				f_H2_step[i][j][k]  = f_H2[i][j][k];
				f_He1_step[i][j][k] = f_He1[i][j][k];
				f_He2_step[i][j][k] = f_He2[i][j][k];
				f_He3_step[i][j][k] = f_He3[i][j][k];
				
				//recombination coefficients
				if (strcmp(rec_case, "A") == 0)  {
					recomb_H2[i][j][k]  = alphaA_H2(temp[i][j][k]);
					recomb_He2[i][j][k] = alphaA_He2(temp[i][j][k]);
					recomb_He3[i][j][k] = alphaA_He3(temp[i][j][k]);
				}
				else  {
					recomb_H2[i][j][k]  = alphaB_H2(temp[i][j][k]);
					recomb_He2[i][j][k] = alphaB_He2(temp[i][j][k]);
					recomb_He3[i][j][k] = alphaB_He3(temp[i][j][k]);
				}
				
				//compute first derivatives of the number densities
				double dnH2_dt  = gamma_H1_tot[i][j][k]*nH1[i][j][k] - recomb_H2[i][j][k]*ne[i][j][k]*nH2[i][j][k];
				double dnHe2_dt = nHe1[i][j][k]*gamma_He1_tot[i][j][k] - recomb_He2[i][j][k]*ne[i][j][k]*nHe2[i][j][k]
								  + recomb_He3[i][j][k]*ne[i][j][k]*nHe3[i][j][k] - nHe2[i][j][k]*gamma_He2_tot[i][j][k];
				double dnHe3_dt = nHe2[i][j][k]*gamma_He2_tot[i][j][k] - recomb_He3[i][j][k]*ne[i][j][k]*nHe3[i][j][k];
				dne_dt[i][j][k]   = dnH2_dt + dnHe2_dt + 2*dnHe3_dt;
			}
		}
	}
	}
}

