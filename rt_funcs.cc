#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "ray_tracing.cc"

int update_dt(void)  {
	dt = tstep_factor*dx[0]/cl_factor/cl;
}

int update_step(void)  {
	int i;
	
	#pragma omp parallel
	{
	#pragma omp for
	for (i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				f_H1_step[i][j][k]  = f_H1[i][j][k];
				f_H2_step[i][j][k]  = f_H2[i][j][k];
				f_He1_step[i][j][k] = f_He1[i][j][k];
				f_He2_step[i][j][k] = f_He2[i][j][k];
				f_He3_step[i][j][k] = f_He3[i][j][k];
			}
		}
	}
	}
	t_step = t;
}

int update_gamma_nu(void)  {

	int i;
	
	#pragma omp parallel
	{
	#pragma omp for
	for (i = 0; i < Nx; i++)   {
		for (int j = 0; j < Ny; j++)   {
			for (int k = 0; k < Nz; k++)   {
				//update absorption coefficient
				gamma_nu_H1[i][j][k]  = nH1[i][j][k]*sigmapi_H1(nu_phot);
				gamma_nu_He1[i][j][k] = nHe1[i][j][k]*sigmapi_He1(nu_phot);
				gamma_nu_He2[i][j][k] = nHe2[i][j][k]*sigmapi_He2(nu_phot);
				gamma_nu_tot[i][j][k] = gamma_nu_H1[i][j][k]
										 + gamma_nu_He1[i][j][k] + gamma_nu_He2[i][j][k];
			}
		}
	}
	}
}

int init_rt(void)  {
	int i;
	
	//source field 
	read_source_catalog();

	if (input_grid == TRUE)  {
		printf("Reading ray data\n");
		read_rays_binary();
	}
	if (input_grid == FALSE)  {
		#pragma omp parallel
		{
		#pragma omp for
		for (i = 0; i < Nx; i++)  {
			for (int j = 0; j < Ny; j++)  {
				for (int k = 0; k < Nz; k++)  {
					if (source_lum[i][j][k] != 0)  {
						init_rays(i, j, k);
					}
				}
			}
		}
		}
	}
	ray_dt = 0.;
}

int set_clight(void)  {
	int i;
	#pragma omp parallel
	{
	#pragma omp for
	for (i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				clight[i][j][k] = cl_factor*cl;
			}
		}
	}
	}
}

int update_rays(void)  {
	int i;
	int ray_counter = 0;
	int lmax = 6; //maximum allowed healpix level
	
	ray_dt += dt;
	ray_counter = 0;
	
	//advance all the rays by one time step
	#pragma omp parallel
	{
	#pragma omp for
	for (i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				int x = 0;
				int size = (int) ray_index[i][j][k].size();
				while (x < size)  {
					ray_counter += 1;
					if (ray[i][j][k][x][11] == 0)  {
						advance_ray(i, j, k, x);
					}
					x += 1;
				}
			}
		}
	}
	}

	//remove rays
	#pragma omp parallel
	{
	#pragma omp for
	for (i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				int x = 0;
				while (x < ray_index[i][j][k].size())  {
					if (ray[i][j][k][x][11] == -1)  {
						remove_ray(i, j, k, x);
						x -= 1;
					}
					else  {
						ray[i][j][k][x][11] = 0;
					}
					x += 1;
				}
			}
		}
	}
	}

	//don't cast rays more often than every half light-crossing time
	double dtime = maxd(0.5*dx[0]/cl_factor/cl, dt);

	if (ray_dt >= dtime)  {
		#pragma omp parallel
		{
		#pragma omp for 
		for (i = 0; i < Nx; i++)  {
			for (int j = 0; j < Ny; j++)  {
				for (int k = 0; k < Nz; k++)  {
					int num = (int) ray_index[i][j][k].size();
					int x = 0;
					if ( (num > 0) && (num < 2))  {
						while ( (x < num) && (ray_index[i][j][k][x][0] < lmax) )  {
							split_ray(i, j, k, x);
							x += 1;
						}
					}
					x = 0;
					if (source_lum[i][j][k] != 0)  {
						init_rays(i, j, k);
					}
				}
			}
		}
		}
		ray_dt = 0.;
	}
	
	printf("total rays: \n");
	printf("%d\n", ray_counter);
	printf("%d\n", total_rays);
}

int merge(void)  {
	int i;

	#pragma omp parallel
	{
	#pragma omp for
	for (i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  {
				int num = (int) ray_index[i][j][k].size();
				if ( (num > 12*pow(2, 2*hpx_lvl)) && (source_lum[i][j][k] == 0) )  {
					merge_cell(i, j, k, hpx_lvl);
				}
			}
		}
	}
	}
}










