//Input/output functions
#include <omp.h>
#include <fstream>
#include <stdio.h>
#include <dirent.h> 
#include <string.h>
#include "data_funcs.cc"

int get_ugas(void)  {
	int i, j, k;
	ifstream file (ugas_file, ios::in | ios::binary);

	float  temp_float;
	double om = 0.305147;
	double ol = 1 - om;
   	double ob = 0.0482266;
    double h0 = 0.68;
	double rho_c_cgs = 1.87890e-29;
	
	double rho_unit = rho_c_cgs*om*pow(hh, 2.) * pow(1. + zz, 3.);

	for (k = 0; k < Nz; k++)  {
		for (j = 0; j < Ny; j++)  {
			for (i = 0; i < Nx; i++)  {
				file.read((char*)&temp_float, sizeof(float));
				rho[i][j][k] = (double) temp_float;
				rho[i][j][k] *= rho_unit;
				file.read((char*)&temp_float, sizeof(float));
				file.read((char*)&temp_float, sizeof(float));
				file.read((char*)&temp_float, sizeof(float));
				file.read((char*)&temp_float, sizeof(float));
			}
		}
	}
}

int read_grid_binary(void)  {
	int i, j, k;
	ifstream file (start_file, ios::in | ios::binary);
	
	for (i = 0; i < Nx; i++)  {
		for (j = 0; j < Ny; j++)  {
			for (k = 0; k < Nz; k++)  {
				file.read((char*)&x[i], sizeof(double));
				dx[i] = Lx/Nx*kpc_to_cm;
				file.read((char*)&y[j], sizeof(double));
				dy[j] = Ly/Ny*kpc_to_cm;
				file.read((char*)&z[k], sizeof(double));
				dz[k] = Lz/Nz*kpc_to_cm;
				file.read((char*)&rho[i][j][k], sizeof(double));
				file.read((char*)&temp[i][j][k], sizeof(double));
				file.read((char*)&nH1[i][j][k], sizeof(double));
				file.read((char*)&f_H1[i][j][k], sizeof(double));
				file.read((char*)&f_He1[i][j][k], sizeof(double));
				file.read((char*)&f_H2[i][j][k], sizeof(double));
				file.read((char*)&f_He2[i][j][k], sizeof(double));
				file.read((char*)&f_He3[i][j][k], sizeof(double));
				file.read((char*)&u_nu[i][j][k], sizeof(double));
			}
		}
	}
}

int read_rays_binary(void)  {
	int ix, jx, kx, m, x;
	ifstream file;
	vector<int> index_new;
	vector<double> ray_new;
	index_new.resize(2);
	ray_new.resize(13);
	
	file.open(ray_file, ios::in | ios::binary);
	file.read((char*)&ix, sizeof(int));
	file.read((char*)&jx, sizeof(int));
	file.read((char*)&kx, sizeof(int));
	while (ix != -1)  {
		file.read((char*)&index_new[0], sizeof(int));
		file.read((char*)&index_new[1], sizeof(int));
		for (m = 0; m < 13; m++)  {
			file.read((char*)&ray_new[m], sizeof(double));
		}
		ray_index[ix][jx][kx].push_back(index_new);
		ray[ix][jx][kx].push_back(ray_new);
		vector<double> path{ray_new[6]};
		ray_dep[ix][jx][kx].push_back(path);
		
		file.read((char*)&ix, sizeof(int));
		file.read((char*)&jx, sizeof(int));
		file.read((char*)&kx, sizeof(int));
	}
}

int read_source_catalog(void)  {
	FILE *file = NULL;
	int si, ix, jx, kx;
	double x, y, z;
	double LumGAL;
	int num_src = 184991;
	int num_inc = 184991;
	int CATDIM  = 4;
	double gal_source_catalog[num_src*CATDIM];
	
	file = fopen(source_field, "r");
	if (file == NULL)  {
                printf("Source field not found\n");
                return -1;
        }
	else  {
		fread(gal_source_catalog, sizeof(double), num_src*CATDIM, file);
		for (si = 0; si < num_inc; si++)  {
			LumGAL = gal_source_catalog[si*CATDIM];
			x   = gal_source_catalog[si*CATDIM + 1];
            y   = gal_source_catalog[si*CATDIM + 2];
            z   = gal_source_catalog[si*CATDIM + 3];
			ix  = (int) floor(x/hh/(1. + zz)/Lx*1e3*Nx);
			jx  = (int) floor(y/hh/(1. + zz)/Ly*1e3*Ny);
			kx  = (int) floor(z/hh/(1. + zz)/Lz*1e3*Nz);
			source_lum[ix][jx][kx] += pow(10., 25.5)*fesc*LumGAL*h*nu;
		}
	}
}

int make_output(void)  {

	FILE *file = NULL;

	file = fopen(gas_output, "w");
	fprintf(file, "Gas data\n");
	fclose(file);

	if (input_grid == FALSE)  {
		file = fopen(otf_output, "w");
		fprintf(file, "On-the-fly output\n");
		fprintf(file, "\n");
		fclose(file);
	}
	
	file = fopen(ray_output, "w");
	fprintf(file, "\n");
	fclose(file);
}

int write_rays_binary(void)  {
	int i, j, k, m, x;
	ofstream file;
	int end = -1;
	printf("Writing ray file\n");
	file.open(ray_output, ios::out | ios::binary);
	for (i = 0; i < Nx; i++)  {
		for (j = 0; j < Ny; j++)  {
			for (k = 0; k < Nz; k++)  {
				for (x = 0; x < ray_index[i][j][k].size(); x++)  {
					file.write((char*)&i, sizeof(int));
					file.write((char*)&j, sizeof(int));
					file.write((char*)&k, sizeof(int));
					file.write((char*)&ray_index[i][j][k][x][0], sizeof(int));
					file.write((char*)&ray_index[i][j][k][x][1], sizeof(int));
					for (m = 0; m < 13; m++)  {
						file.write((char*)&ray[i][j][k][x][m], sizeof(double));
					}
				}
			}
		}
	}
	file.write((char*)&end, sizeof(int));
	file.write((char*)&end, sizeof(int));
	file.write((char*)&end, sizeof(int));
}

int write_treion(void)  {
	int i, j, k;
	ofstream file;
	
	printf("Writing treion file\n");
	file.open(tre_output, ios::out | ios::binary);
	
	for (i = 0; i < Nx; i++)  {
		for (j = 0; j < Ny; j++)  {
			for (k = 0; k < Nz; k++)  {
				file.write((char*)&treion[i][j][k], sizeof(double));
			}
		}
	}
}

int write_gas_binary(void)  {
	int i, j, k;
	ofstream file;
	printf("Writing gas file\n");
	file.open(gas_output, ios::out | ios::binary);
	for (i = 0; i < Nx; i++)  {
		for (j = 0; j < Ny; j++)  {
			for (k = 0; k < Nz; k++)  {
				file.write((char*)&x[i], sizeof(double));
				file.write((char*)&y[j], sizeof(double));
				file.write((char*)&z[k], sizeof(double));
				file.write((char*)&rho  [i][j][k], sizeof(double));
				file.write((char*)&temp [i][j][k], sizeof(double));
				file.write((char*)&nH1  [i][j][k], sizeof(double));
				file.write((char*)&f_H1 [i][j][k], sizeof(double));
				file.write((char*)&f_He1[i][j][k], sizeof(double));
				file.write((char*)&f_H2 [i][j][k], sizeof(double));
				file.write((char*)&f_He2[i][j][k], sizeof(double));
				file.write((char*)&f_He3[i][j][k], sizeof(double));
				file.write((char*)&u_nu [i][j][k], sizeof(double));
			}
		}
	}
}

int write_otf(void)  {

	int i;
	FILE *file = NULL;

	file = fopen(otf_output, "a");

	fprintf(file, "Step Number    : %d\n", step);
	fprintf(file, "Time [Myr]     : %le\n", t/yr_to_s/1e6);
	fprintf(file, "Time Step [Myr]: %le\n", dt/yr_to_s/1e6);
	fprintf(file, "\n");
    
	double fH1_avg_vol   = calc_vol_avg(f_H1);
	double fH2_avg_vol   = calc_vol_avg(f_H2);
	double fHe1_avg_vol  = calc_vol_avg(f_He1);
	double fHe2_avg_vol  = calc_vol_avg(f_He2);
	double fHe3_avg_vol  = calc_vol_avg(f_He3);
	
	double fH1_avg_mass  = calc_mass_avg(f_H1, nH);
	double fH2_avg_mass  = calc_mass_avg(f_H2, nH);
	double fHe1_avg_mass = calc_mass_avg(f_He1, nHe);
	double fHe2_avg_mass = calc_mass_avg(f_He2, nHe);
	double fHe3_avg_mass = calc_mass_avg(f_He3, nHe);

	fprintf(file, "Global Averages\n");
	fprintf(file, "fH1  : %le    %le\n", fH1_avg_vol, fH1_avg_mass);
	fprintf(file, "fH2  : %le    %le\n", fH2_avg_vol, fH2_avg_mass);
	fprintf(file, "fHe1 : %le    %le\n", fHe1_avg_vol, fHe1_avg_mass);
	fprintf(file, "fHe2 : %le    %le\n", fHe2_avg_vol, fHe2_avg_mass);
	fprintf(file, "fHe3 : %le    %le\n", fHe3_avg_vol, fHe3_avg_mass);
	fprintf(file, "\n");

	/*int center[3] = {Nx/2, Ny/2, Nz/2};

	printf("computing ifronts\n");
	double *ifront_H1  = calc_ifront(f_H1, f_H1_step, t - t_step, center);
	double *ifront_He1 = calc_ifront(f_He1, f_He1_step, t - t_step, center);
	double *ifront_He3 = calc_ifront(f_He3, f_He3_step, t - t_step, center);

	printf("writing I-fronts\n");
	fprintf(file, "I-fronts\n");
	fprintf(file, "H1   : %le    %le\n", *(ifront_H1 + 0), *(ifront_H1 + 1));
	fprintf(file, "He1  : %le    %le\n", *(ifront_He1 + 0), *(ifront_He1 + 1));
	fprintf(file, "He3  : %le    %le\n", *(ifront_He3 + 0), *(ifront_He3 + 1));
	fprintf(file, "\n");*/

	fclose(file);
}
