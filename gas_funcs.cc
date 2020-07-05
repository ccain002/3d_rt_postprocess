#include <math.h>
#include <omp.h>
#include <stdio.h>
#include "init_funcs.cc"

int update_gamma(void)  {

	int i;

	#pragma omp parallel
	{
	#pragma omp for
	for (i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				gamma_H1_tot[i][j][k]  = sigmapi_H1(nu_phot)*clight[i][j][k]*u_nu[i][j][k]/h/nu_phot;
				gamma_He1_tot[i][j][k] = sigmapi_He1(nu_phot)*clight[i][j][k]*u_nu[i][j][k]/h/nu_phot;
				gamma_He2_tot[i][j][k] = sigmapi_He2(nu_phot)*clight[i][j][k]*u_nu[i][j][k]/h/nu_phot;
				//collisional ionization counts as an effective photoionization rate with coeff
				//cic x ne
				if ( (coll_ion == TRUE) && (temp[i][j][k] >= 1e4) && (temp[i][j][k] <= 1e9) )  {
					gamma_H1_tot[i][j][k]  += cic_H1(temp[i][j][k])*ne[i][j][k];
					gamma_He1_tot[i][j][k] += cic_He1(temp[i][j][k])*ne[i][j][k];
					gamma_He2_tot[i][j][k] += cic_He2(temp[i][j][k])*ne[i][j][k];
				}
			}
		}
	}
	}
}

int update_heat_cool(void)  {

	int i;
	double dheat, dcool;

	#pragma omp parallel
	{
	#pragma omp for
	for (i = 0; i < Nx; i++)  {
		for (int j = 0; j < Nx; j++)  {
			for (int k = 0; k < Nx; k++)  {
				heat_rate[i][j][k] = 0;
				cool_rate[i][j][k] = 0;
				
				heat_rate[i][j][k] += nH1[i][j][k]*sigmapi_H1(nu_phot)*clight[i][j][k]*u_nu[i][j][k]*(h_eV*nu_phot - Eth_H1)/h_eV/nu_phot;
				heat_rate[i][j][k] += nHe1[i][j][k]*sigmapi_He1(nu_phot)*clight[i][j][k]*u_nu[i][j][k]*(h_eV*nu_phot - Eth_He1)/h_eV/nu_phot;
				heat_rate[i][j][k] += nHe2[i][j][k]*sigmapi_He2(nu_phot)*clight[i][j][k]*u_nu[i][j][k]*(h_eV*nu_phot - Eth_He2)/h_eV/nu_phot;
					
				if ( (coll_ion == TRUE) && (temp[i][j][k] >= 1e4) && (temp[i][j][k] <= 1e9) )  {
					//collisional ionization cooling contributions from each ion
					cool_rate[i][j][k] += cicr_H1(temp[i][j][k])*ne[i][j][k]*nH1[i][j][k];
					cool_rate[i][j][k] += cicr_He1(temp[i][j][k])*ne[i][j][k]*nHe1[i][j][k];
					cool_rate[i][j][k] += cicr_He2(temp[i][j][k])*ne[i][j][k]*nHe2[i][j][k];
				}

				//case A or B recombination cooling
				if (recomb_cool == TRUE)  {
					if (strcmp(rec_case, "A") == 0)  {
						cool_rate[i][j][k] += rcrA_H2(temp[i][j][k])*ne[i][j][k]*nH2[i][j][k];
						cool_rate[i][j][k] += rcrA_He2(temp[i][j][k])*ne[i][j][k]*nHe2[i][j][k];
						cool_rate[i][j][k] += rcrA_He3(temp[i][j][k])*ne[i][j][k]*nHe3[i][j][k];
					}
					else if (strcmp(rec_case, "B") == 0) {
						cool_rate[i][j][k] += rcrB_H2(temp[i][j][k])*ne[i][j][k]*nH2[i][j][k];
						cool_rate[i][j][k] += rcrB_He2(temp[i][j][k])*ne[i][j][k]*nHe2[i][j][k];
						cool_rate[i][j][k] += rcrB_He3(temp[i][j][k])*ne[i][j][k]*nHe3[i][j][k];
					}
					else  {
						printf("Warning: Ignoring recombinations\n");
					}
				}
				//collisional excitation cooling contributions from HI, HeI, and HeII
				if (coll_exc_cool == TRUE)  {
					cool_rate[i][j][k] += coll_ex_rate_H1(temp[i][j][k])*ne[i][j][k]*nH1[i][j][k];
					cool_rate[i][j][k] += coll_ex_rate_He1(temp[i][j][k])*pow(ne[i][j][k],2.)*nHe1[i][j][k];
					cool_rate[i][j][k] += coll_ex_rate_He2(temp[i][j][k])*ne[i][j][k]*nHe2[i][j][k];
				}
				//compton heating/cooling.  Counts as a negative heating rate if T > T_cmb.  
				if (compton == TRUE)  {
					heat_rate[i][j][k] += compton_rate(zz)*ne[i][j][k]/rho[i][j][k]*(2.726*(1 + zz) - temp[i][j][k]);
				}
			}
		}
	}
	}
}

int solve_ion(int i, int j, int k)  {
	
	//solve the ionizing balance equation for H usin a 1st order implicit backwards-difference scheme 
	if (nH1_prev[i][j][k] >= nH2_prev[i][j][k])  {
		//Solve for nH2 when there is less ionized gas
		nH2[i][j][k]	= (nH2_prev[i][j][k] + gamma_H1_tot[i][j][k]*nH[i][j][k]*dt)
						/(1. + gamma_H1_tot[i][j][k]*dt + recomb_H2[i][j][k]*ne[i][j][k]*dt);
		nH1[i][j][k] 	= nH[i][j][k] - nH2[i][j][k];
	}
	else  {
		//solve for nH1 when there is less neutral gas
		nH1[i][j][k]	= (nH1_prev[i][j][k] + recomb_H2[i][j][k]*ne[i][j][k]*nH[i][j][k]*dt)
						/(1. + gamma_H1_tot[i][j][k]*dt + recomb_H2[i][j][k]*ne[i][j][k]*dt);
		nH2[i][j][k]  = nH[i][j][k] - nH1[i][j][k];
	}
	
	//solve the ionizing balance equations for He.  Solve the backwards difference equation for the
	//species that have the smaller abundances, and solve closing equation for the remaining species.  
	if ( (nHe1_prev[i][j][k] <= nHe3_prev[i][j][k]) && (nHe2_prev[i][j][k] <= nHe3_prev[i][j][k]) )  {
		
		nHe1[i][j][k] = (nHe1_prev[i][j][k] + recomb_He2[i][j][k]*ne[i][j][k]*(nHe[i][j][k] - nHe3[i][j][k])*dt)
				      /(1. + gamma_He1_tot[i][j][k]*dt + recomb_He2[i][j][k]*ne[i][j][k]*dt);
		nHe2[i][j][k] = (nHe2_prev[i][j][k] + gamma_He1_tot[i][j][k]*(nHe[i][j][k] - nHe3[i][j][k])*dt
					  + recomb_He3[i][j][k]*ne[i][j][k]*(nHe[i] - nHe1[i])*dt)
				      /(1. + (gamma_He1_tot[i][j][k] + gamma_He2_tot[i][j][k])*dt
					  + (recomb_He2[i][j][k] + recomb_He3[i][j][k])*ne[i][j][k]*dt);
		nHe3[i][j][k] = nHe[i][j][k] - nHe1[i][j][k] - nHe2[i][j][k];
	}
	else if( (nHe1_prev[i][j][k] <= nHe2_prev[i][j][k]) && (nHe3_prev[i][j][k] <= nHe2_prev[i][j][k]) )  {
		nHe1[i][j][k] = (nHe1_prev[i][j][k] + recomb_He2[i][j][k]*ne[i][j][k]*(nHe[i][j][k] - nHe3[i][j][k])*dt)
					  /(1. + gamma_He1_tot[i][j][k]*dt + recomb_He2[i][j][k]*ne[i][j][k]*dt);
		nHe3[i][j][k] = (nHe3_prev[i][j][k]+ gamma_He2_tot[i][j][k]*(nHe[i][j][k] - nHe1[i][j][k])*dt)
					  /(1. + gamma_He2_tot[i][j][k]*dt + recomb_He3[i][j][k]*ne[i][j][k]*dt);
		nHe2[i][j][k] = nHe[i][j][k] - nHe1[i][j][k] - nHe3[i][j][k];
	}
	else {
		nHe2[i][j][k] = (nHe2_prev[i][j][k] + gamma_He1_tot[i][j][k]*(nHe[i][j][k] - nHe3[i][j][k])*dt
					  + recomb_He3[i][j][k]*ne[i][j][k]*(nHe[i][j][k] - nHe1[i][j][k])*dt)
					  /(1. + (gamma_He1_tot[i][j][k] + gamma_He2_tot[i][j][k])*dt + (recomb_He2[i][j][k]
					  + recomb_He3[i][j][k])*ne[i][j][k]*dt);
		nHe3[i][j][k] = (nHe3_prev[i][j][k] + gamma_He2_tot[i][j][k]*(nHe[i][j][k] - nHe1[i][j][k])*dt)
					  /(1.+ gamma_He2_tot[i][j][k]*dt + recomb_He3[i][j][k]*ne[i][j][k]*dt);
		nHe1[i][j][k] = nHe[i][j][k] - nHe2[i][j][k] - nHe3[i][j][k];
	}
	
	//update free electron and total number densities
	ne[i][j][k]    = nH2[i][j][k] + nHe2[i][j][k] + 2.*nHe3[i][j][k];
}

int update_chem(void)  {

	int i;

	#pragma omp parallel
	{
	#pragma omp for
	for (i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				
				dne_dt[i][j][k]   = -ne[i][j][k];
				
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
				
				//diaelectric recombination of He II.  
				if ( (temp[i][j][k] >= 3e4) && (temp[i][j][k] <= 1e6) )  {
					recomb_He2[i][j][k] += Dalpha_He2(temp[i][j][k]);
				} 
				
				for (int m = 0; m < 1; m++)  {
					solve_ion(i, j, k);
				}

				n_tot[i][j][k] = nH[i][j][k] + nHe[i][j][k] + ne[i][j][k];
				
				dne_dt[i][j][k] += ne[i][j][k];
				dne_dt[i][j][k]	/= dt;
				
				nH1_prev[i][j][k]  = nH1[i][j][k];
				nH2_prev[i][j][k]  = nH2[i][j][k];
				nHe1_prev[i][j][k] = nHe1[i][j][k];
				nHe2_prev[i][j][k] = nHe2[i][j][k];
				nHe3_prev[i][j][k] = nHe3[i][j][k];
				
				//update abundance fractions
				f_H1[i][j][k]  = nH1[i][j][k]/nH[i][j][k];
				f_H2[i][j][k]  = nH2[i][j][k]/nH[i][j][k];
				f_He1[i][j][k] = nHe1[i][j][k]/nHe[i][j][k];
				f_He2[i][j][k] = nHe2[i][j][k]/nHe[i][j][k];
				f_He3[i][j][k] = nHe3[i][j][k]/nHe[i][j][k];
			}
		}
	}
	}
}

int update_thermal(void)  {

	int i;
	double a = 1/(1 + zz);
	
	if (temp_ev == TRUE)  {
		#pragma omp parallel
		{
		#pragma omp for
		for (i = 0; i < Nx; i++)  {
			for (int j = 0; j < Ny; j++)  {
				for (int k = 0; k < Nz; k++)  {
					temp[i][j][k] 	= (temp[i][j][k] + 2./3./k_B/n_tot[i][j][k]*(heat_rate[i][j][k] - cool_rate[i][j][k])*dt)
									/(1. + 1./n_tot[i][j][k]*dne_dt[i][j][k]*dt);
				}
			}
		}
		}
	}
	else  {
		#pragma omp parallel
		{
		#pragma omp for
		for (i = 0; i < Nx; i++)  {
			for (int j = 0; j < Ny; j++)  {
				for (int k = 0; k < Nz; k++)  {
					if (nH1[i][j][k]/nH[i][j][k] > 0.5)  {
						temp[i][j][k] = temp_0;
					}
					else  {
						temp[i][j][k] = 1e4;
					} 
				}
			}
		}
		}
	}
}

