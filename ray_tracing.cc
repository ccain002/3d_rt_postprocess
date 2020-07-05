#include <omp.h>
#include <stdio.h>
#include <iostream> 
#include <vector>
#include <array>
#include <algorithm>
#include </home/ccain002/Healpix_3.40/src/cxx/Healpix_cxx/healpix_base.h>
#include "gas_funcs.cc"
using namespace std;

//ray index
//0 = healpix level
//1 = healpix pixel number

//ray
//0, 1, 2 = rotation angles
//3, 4, 5 = x, y, and z unit vector components
//6  = photon count
//7  = attenuation factor
//8  = x position within a cell
//9  = y position
//10 = z position
//11 = update tag
//12 = time

//ray deposition
//0 = photons deposited in cell
//4p + 1, 2, 3 = ip, jp, kp along path
//4p + 4 = distance traveled within cell

int init_rays(int i, int j, int k)  {
	
	int m;
	vector<int> index_init;
	vector<double> ray_init;
	index_init.resize(2);
	ray_init.resize(13);
	int f, pnp;
	double z, phi;
	int l     = hpx_lvl;
	int Nside = pow(2, l);
	double alpha0, beta0, gamma0;

	alpha0 = fRand(0., 2.*pi);
	beta0  = fRand(0., 2.*pi);
	gamma0 = fRand(0., 2.*pi);
	
	for (f = 0; f < 12; f++)  {
		for (pnp = 0; pnp < pow(Nside, 2); pnp++)  {
			
			//testing
			total_rays += 1;
			
			index_init[0] = l;
			index_init[1] = f*pow(Nside, 2) + pnp;
			
			ray_init[0] = alpha0;
			ray_init[1] = beta0;
			ray_init[2] = gamma0;
			
			auto [ix, jx] = get_ij(f*pow(Nside, 2) + pnp, Nside);
			auto [z, phi] = get_zphi(ix, jx, f*pow(Nside, 2) + pnp, Nside);
			auto [xhat, yhat, zhat] = get_unit_vector(z, phi, alpha0, beta0, gamma0, 0);
			
			ray_init[3] = xhat;
			ray_init[4] = yhat;
			ray_init[5] = zhat;
			ray_init[6]  = source_lum[i][j][k]*ray_dt/12./pow(Nside,2.)/h/nu_phot;
			ray_init[7]  = 1.;
			ray_init[8]  = dx[i]/2.;
			ray_init[9]  = dy[j]/2.;
			ray_init[10] = dz[k]/2.;
			ray_init[11] = 0.;
			ray_init[12] = 0.;
			
			ray_index[i][j][k].push_back(index_init);
			ray[i][j][k].push_back(ray_init);
			
			vector<double> path = {source_lum[i][j][k]*ray_dt/12./pow(Nside,2.)/h/nu_phot};
			ray_dep[i][j][k].insert(ray_dep[i][j][k].end(), path);
		}
	}
}

int split_ray(int i, int j, int k, int x)  {
	int n, m;
	vector<int> index_new;
	vector<double> ray_new;
	index_new.resize(2);
	ray_new.resize(13);
	int l  = ray_index[i][j][k][x][0];
	int pn = ray_index[i][j][k][x][1];
	int Nside = pow(2, l + 1);
	double alpha0 = ray[i][j][k][x][0];
	double beta0  = ray[i][j][k][x][1];
	double gamma0 = ray[i][j][k][x][2];
	
	auto [ix, jx] = get_ij(4*pn, Nside);
	auto [z, phi] = get_zphi(ix, jx, 4*pn, Nside);
	auto [xhat, yhat, zhat] = get_unit_vector(z, phi, alpha0, beta0, gamma0, 0);
	
	//update parameters for split ray
	ray_index[i][j][k][x][0]   = l + 1;
	ray_index[i][j][k][x][1]   = 4*pn;
	
	ray[i][j][k][x][3] = xhat;
	ray[i][j][k][x][4] = yhat;
	ray[i][j][k][x][5] = zhat;
	
	ray     [i][j][k][x][6] /= 4.; //divide photon count by 4
	ray_dep [i][j][k][x][0] /= 4.;
	
	//add the new rays at the front.  
	for (n = 1; n < 4; n++)   {
		
		index_new[0] = l + 1;
		index_new[1] = 4*pn + n;
		
		auto [ix, jx] = get_ij(4*pn + n, Nside);
		auto [z, phi] = get_zphi(ix, jx, 4*pn + n, Nside);
		auto [xhat, yhat, zhat] = get_unit_vector(z, phi, alpha0, beta0, gamma0, 0);
		
		ray_new[0] = alpha0;
		ray_new[1] = beta0;
		ray_new[2] = gamma0;
		ray_new[3] = xhat;
		ray_new[4] = yhat;
		ray_new[5] = zhat;
		
		for (m = 6; m < 13; m++)  {
			ray_new[m] = ray[i][j][k][x][m];
		}
		
		ray_index[i][j][k].push_back(index_new);
		ray[i][j][k].push_back(ray_new);
		ray_dep[i][j][k].insert(ray_dep[i][j][k].end(), ray_dep[i][j][k][x]);
		
		//testing
		total_rays += 1;
	}
}

int remove_ray(int i, int j, int k, int x)  {
	ray_index[i][j][k].erase(ray_index[i][j][k].begin() + x);
	ray[i][j][k].erase(ray[i][j][k].begin() + x);
	ray_dep[i][j][k].erase(ray_dep[i][j][k].begin() + x);
}

//attenuate the intensity of the ray over time step dt and .  
vector<double> update_path(int ind, double dist, vector<double> path)  {
	
	//update path vector
	path.insert(path.end(), ind);
	path.insert(path.end(), dist);
	
	return path;
}

int find_avg_unu(void)  {
	//iterate to get the average value of gamma, and therefore the averge u_nu to put into the 
	//photoionization rate calculation.  
	int i, iter;
	
	for (iter = 0; iter < itercount; iter++)  { //iterate
	
		#pragma omp parallel
		{
		#pragma omp for
		for (i = 0; i < Nx; i++)  {
			for (int j = 0; j < Ny; j++)  {
				for (int k = 0; k < Nz; k++)  {
					avg_gamma[i][j][k] = 0.;
					if (iter == 0)  {
						avg_xh1[i][j][k] = f_H1[i][j][k]; //get initial guess for HI fraction
					}
				}
			}
		}
		}
		
		#pragma omp parallel
		{
		#pragma omp for
		for (i = 0; i < Nx; i++)  {
			for (int j = 0; j < Ny; j++)  {
				for (int k = 0; k < Nz; k++)  { //loop over grid cells
					//compute mean photoionization rate in each cell using HI fraction from previous
					//time step as the initial guess
					for (int x = 0; x < ray_dep[i][j][k].size(); x++)  {
						double nphot  = ray_dep[i][j][k][x][0];
						double tau    = 0.;
						double dtau   = 0.;
						for (int p = 1; p < ray_dep[i][j][k][x].size(); p+=2)  {
							int ind  = (int) ray_dep[i][j][k][x][p + 0];
							int ip   = (int) floor(ind/pow(Nx, 2));
							int jp   = (int) floor((ind - ip*pow(Nx, 2))/Ny);
							int kp   = ind - ip*pow(Nx, 2) - jp*Ny;
							dtau = nH[ip][jp][kp]*avg_xh1[ip][jp][kp]*sigmapi_H1(nu_phot)*ray_dep[i][j][k][x][p + 1];
							avg_gamma[ip][jp][kp] += nphot*exp(-tau)/dt/(dx[ip]*dy[jp]*dz[kp])/nH[ip][jp][kp]
												   /avg_xh1[ip][jp][kp] * (1. - exp(-dtau));
							tau += dtau;
						}
						if (iter == itercount - 1)  {
							ray[i][j][k][x][6] *= exp(-tau);
							ray[i][j][k][x][7] *= exp(-tau);
						}
					}
				}
			}
		}
		}
	
		#pragma omp parallel
		{
		#pragma omp for
		for (i = 0; i < Nx; i++)  {
			for (int j = 0; j < Ny; j++)  {
				for (int k = 0; k < Nz; k++)  { //loop over grid cells
					double nelec = nH[i][j][k]*(1. - avg_xh1[i][j][k]);
					double col = 0.;
					if ( (coll_ion == TRUE) && (temp[i][j][k] >= 1e4) && (temp[i][j][k] <= 1e9) )  {
						col = cic_H1(temp[i][j][k])*nelec;
					}
					//solve for average HI fraction assuming constant mean photoionization rate
					double x0    = 1. - nH1[i][j][k]/nH[i][j][k];
					double xeq   = (avg_gamma[i][j][k] + col)
						  /(avg_gamma[i][j][k] + col + recomb_H2[i][j][k]*nelec);
					double ti    = 1./(avg_gamma[i][j][k] + col + recomb_H2[i][j][k]*nelec);
					avg_xh1[i][j][k] = 1. - (xeq + ti/dt*(x0 - xeq)*(1. - exp(-dt/ti)));
				}
			}
		}
		}
	}
	
	//update energy and reset deposition paths
	#pragma omp parallel
	{
	#pragma omp for 
	for (i = 0; i < Nx; i++)  {
		for (int j = 0; j < Ny; j++)  {
			for (int k = 0; k < Nz; k++)  { //loop over grid cells
				u_nu[i][j][k] = avg_gamma[i][j][k]*h*nu_phot/clight[i][j][k]/sigmapi_H1(nu_phot);
				for (int x = 0; x < ray_dep[i][j][k].size(); x++)  {
					//reset paths
					vector<double> path = {ray[i][j][k][x][6]};
					ray_dep[i][j][k][x] = path;
				}
			}
		}
	}
	}
}

int advance_ray(int i, int j, int k, int x)  {
	int m, i_new, j_new, k_new, tag;
	vector<int> index_new;
	vector<double> ray_new;
	index_new.resize(2);
	ray_new.resize(13);
	double tx, ty, tz;
	double tmin = 0.;
	int nx, ny, nz;
	double deltat, tlast, dtau;
	int l         = ray_index[i][j][k][x][0];
	int pn        = ray_index[i][j][k][x][1];
	int Nside     = pow(2, l);
	double alpha0 = ray[i][j][k][x][0];
	double beta0  = ray[i][j][k][x][1];
	double gamma0 = ray[i][j][k][x][2];
	double xhat   = ray[i][j][k][x][3];
	double yhat   = ray[i][j][k][x][4];
	double zhat   = ray[i][j][k][x][5];
	double nphot  = ray[i][j][k][x][6];
	double f      = ray[i][j][k][x][7];
	
	double x0   = ray[i][j][k][x][8];
	double y0   = ray[i][j][k][x][9];
	double z0   = ray[i][j][k][x][10];
	
	//create path vector
	vector<double> path = ray_dep[i][j][k][x];
	
	//ray time
	ray[i][j][k][x][12] += dt;
	deltat = (double) ray[i][j][k][x][12];
	
	if (f <= 1e-4)  { 
		ray[i][j][k][x][11] = -1; //mark for removal
	}
	else  {
		if ( (x0 + clight[i][j][k]*xhat*deltat > dx[i]) || (y0 + clight[i][j][k]*yhat*deltat > dy[j])
		  || (z0 + clight[i][j][k]*zhat*deltat > dz[k]) || (x0 + clight[i][j][k]*xhat*deltat < 0.) 
		  || (y0 + clight[i][j][k]*yhat*deltat < 0.)    || (z0 + clight[i][j][k]*zhat*deltat < 0.) )  {
		
			ray[i][j][k][x][12] = 0.; 
			//alternative algorithm
			i_new = i;
			j_new = j;
			k_new = k;
			tlast = 0;
			while (tlast < deltat)  {
				if (xhat > 0)  {
					tx = (dx[i_new] - ray[i][j][k][x][8])/xhat/clight[i_new][j_new][k_new];
				}
				else if (xhat < 0) {
					tx = -ray[i][j][k][x][8]/xhat/clight[i_new][j_new][k_new];
				}
				else {
					tx = 9e99;
				}
				if (yhat > 0)  {
					ty = (dy[j_new] - ray[i][j][k][x][9])/yhat/clight[i_new][j_new][k_new];
				}
				else if (yhat < 0) {
					ty = -ray[i][j][k][x][9]/yhat/clight[i_new][j_new][k_new];
				}
				else  {
					ty = 9e99;
				}
				if (zhat > 0)  {
					tz = (dz[k_new] - ray[i][j][k][x][10])/zhat/clight[i_new][j_new][k_new];
				}
				else if (zhat < 0) {
					tz = -ray[i][j][k][x][10]/zhat/clight[i_new][j_new][k_new];
				}
				else  {
					tz = 9e99;
				}

				if ( (tx <= ty) && (tx <= tz) )  {
					tmin = tx;
					tag = 0;
				}
				else if ( (ty <= tx) && (ty <= tz) )  {
					tmin = ty;	
					tag = 1;
				}
				else  {
					tmin = tz;
					tag = 2;
				}

				dtau = clight[i_new][j_new][k_new]*tmin*nH1[i_new][j_new][k_new]*sigmapi_H1(nu_phot);
			
				if (dtau >= -1.5*log(1e-4/f))  {
					tlast = 9e99;
				}
				
				if ((tlast + tmin) < deltat)  {
					path = update_path(i_new*pow(Nx, 2) + j_new*Ny + k_new, tmin*clight[i_new][j_new][k_new], path);
					ray[i][j][k][x][8]  += tmin*clight[i_new][j_new][k_new]*xhat;
					ray[i][j][k][x][9]  += tmin*clight[i_new][j_new][k_new]*yhat;
					ray[i][j][k][x][10] += tmin*clight[i_new][j_new][k_new]*zhat;
					if (tag == 0)  {
						if (xhat > 0)  {
							ray[i][j][k][x][8] -= dx[i_new];
							i_new = mod(i_new + 1, Nx);
						}
						else  {
							ray[i][j][k][x][8] += dx[i_new];
							i_new = mod(i_new - 1, Nx);
						}
					}
					else if (tag == 1)  {
						if (yhat > 0)  {
							ray[i][j][k][x][9] -= dy[j_new];
							j_new = mod(j_new + 1, Ny);
						}
						else  {
							ray[i][j][k][x][9] += dy[j_new];
							j_new = mod(j_new - 1, Ny);
						}
						
					}
					else if (tag == 2)  {
						if (zhat > 0)  {
							ray[i][j][k][x][10] -= dz[k_new];
							k_new = mod(k_new + 1, Nz);
						}
						else  {
							ray[i][j][k][x][10] += dz[k_new];
							k_new = mod(k_new - 1, Nz);
						}
					}
					else  {
						printf("Invalid tag\n");
					}
				}
				tlast += tmin;
			}
			
			if (tlast - tmin == 9e99)  {
				path = update_path(i_new*pow(Nx, 2) + j_new*Ny + k_new, 5./nH1[i_new][j_new][k_new]/sigmapi_H1(nu_phot), path);
			}
			else  {
				path = update_path(i_new*pow(Nx, 2) + j_new*Ny + k_new, (deltat - tlast + tmin)*clight[i_new][j_new][k_new], path);
				ray[i][j][k][x][8]  += (deltat - tlast + tmin)*clight[i_new][j_new][k_new]*xhat;
				ray[i][j][k][x][9]  += (deltat - tlast + tmin)*clight[i_new][j_new][k_new]*yhat;
				ray[i][j][k][x][10] += (deltat - tlast + tmin)*clight[i_new][j_new][k_new]*zhat;
			}
			
			index_new[0] = ray_index[i][j][k][x][0];
			index_new[1] = ray_index[i][j][k][x][1];
			
			ray_new[0] = ray[i][j][k][x][0];
			ray_new[1] = ray[i][j][k][x][1];
			ray_new[2] = ray[i][j][k][x][2];
			ray_new[3] = xhat;
			ray_new[4] = yhat;
			ray_new[5] = zhat;
			
			for (m = 6; m < 11; m++)   {
				ray_new[m] = ray[i][j][k][x][m];
			}
			
			ray_new[11] = 1.;
			ray_new[12] = ray[i][j][k][x][12];
	
			//update path
			ray_dep[i_new][j_new][k_new].push_back(path);
				
			ray_index[i_new][j_new][k_new].push_back(index_new);
			ray[i_new][j_new][k_new].push_back(ray_new);
			ray[i][j][k][x][11] = -1.;
		}
	}
}

int merge_cell(int i, int j, int k, int l)  {
	int x, xx, xxx, pix;
	int Nside = pow(2, l);
	int numpix = 12*pow(Nside, 2);
	vector<int> pixels[numpix];
	double xhat, yhat, zhat;
	double xhatp, yhatp, zhatp;
	double alpha0, beta0, gamma0;
	double z, phi;

	alpha0 = fRand(0., 2.*pi);
	beta0  = fRand(0., 2.*pi);
	gamma0 = fRand(0., 2.*pi);
	
	for (x = 0; x < ray_index[i][j][k].size(); x++)  {
		xhat = ray[i][j][k][x][3];
		yhat = ray[i][j][k][x][4];
		zhat = ray[i][j][k][x][5];
		z    = zhat;
		
		ray  [i][j][k][x][0] = alpha0;
		ray  [i][j][k][x][1] = beta0;
		ray  [i][j][k][x][2] = gamma0;
		
		if (yhat >= 0.)  {
			phi  = acos(xhat/pow(1. - pow(zhat, 2.), 0.5));
		}
		else  {
			phi = 2.*pi - acos(xhat/pow(1. - pow(zhat, 2.), 0.5));
		}
		
		auto [xhatp, yhatp, zhatp] = get_unit_vector(z, phi, alpha0, beta0, gamma0, 1);
	
		if (xhatp >= 1.)  {
			xhatp = pow(1. - pow(yhatp, 2.) - pow(zhatp, 2.), 0.5);
		}
		if (yhatp >= 1.)  {
			yhatp = pow(1. - pow(xhatp, 2.) - pow(zhatp, 2.), 0.5);
		}
		if (zhatp >= 1.)  {
			zhatp = pow(1. - pow(xhatp, 2.) - pow(yhatp, 2.), 0.5);
		}

		if (xhatp < 0.)  {
			xhatp = 1e-4;
		}
		if (yhatp < 0.)  {
			yhatp = 1e-4;
		}
		if (zhatp < 0.)  {
			zhatp = 1e-4;
		}
		
		T_Healpix_Base<int> b(l, NEST);
		vec3 v(xhat, yhat, zhat);
		pix = (int) b.vec2pix(v);
		pixels[pix].insert(pixels[pix].end(), x);
	}
	for (pix = 0; pix < numpix; pix++)  {
		if (pixels[pix].size() > 0)  {
			x = pixels[pix][0];
			ray_index[i][j][k][x][0]  = l;
			ray_index[i][j][k][x][1]  = pix;
			
			alpha0 = ray[i][j][k][x][0];
			beta0  = ray[i][j][k][x][1];
			gamma0 = ray[i][j][k][x][2];
			
			auto [ix, jx] = get_ij(pix, Nside);
			auto [z, phi] = get_zphi(ix, jx, pix, Nside);
			auto [xhat, yhat, zhat] = get_unit_vector(z, phi, alpha0, beta0, gamma0, 0);
			
			ray[i][j][k][x][3] = xhat;
			ray[i][j][k][x][4] = yhat;
			ray[i][j][k][x][5] = zhat;
			
			ray[i][j][k][x][7]  *= ray[i][j][k][x][6]; 
			ray[i][j][k][x][8]  *= ray[i][j][k][x][6]; 
			ray[i][j][k][x][9]  *= ray[i][j][k][x][6]; 
			ray[i][j][k][x][10]  *= ray[i][j][k][x][6]; 
			//ray time
			ray[i][j][k][x][12] *= ray[i][j][k][x][6]; 
			
			for (xx = 1; xx < pixels[pix].size(); xx++)  {
				xxx = pixels[pix][xx];
				ray[i][j][k][x][6]     += ray[i][j][k][xxx][6];
				ray_dep[i][j][k][x][0] += ray_dep [i][j][k][xxx][0];
				ray[i][j][k][x][7]     += ray[i][j][k][xxx][7]*ray[i][j][k][xxx][6];
				ray[i][j][k][x][8]     += ray[i][j][k][xxx][8]*ray[i][j][k][xxx][6];
				ray[i][j][k][x][9]     += ray[i][j][k][xxx][9]*ray[i][j][k][xxx][6];
				ray[i][j][k][x][10]    += ray[i][j][k][xxx][10]*ray[i][j][k][xxx][6];
				//ray time
				ray[i][j][k][x][12]    += ray[i][j][k][xxx][12]*ray[i][j][k][xxx][6];
				ray_index[i][j][k][xxx][0] = -1;
			}
			
			ray[i][j][k][x][7] /= ray[i][j][k][x][6]; 
			ray[i][j][k][x][8] /= ray[i][j][k][x][6]; 
			ray[i][j][k][x][9] /= ray[i][j][k][x][6]; 
			ray[i][j][k][x][10] /= ray[i][j][k][x][6];
			ray[i][j][k][x][12] /= ray[i][j][k][x][6]; 
		}
	}
	
	xxx = 0;
	while (xxx < ray_index[i][j][k].size())  {
		if (ray_index[i][j][k][xxx][0] == -1)  {
			remove_ray(i, j, k, xxx);
		}
		else  {
			xxx += 1;
		}
	}
}

