#include <math.h>
#include <vector>
#include "user_inputs.h"

//folders
const char x_int_tables[] = "x_int_tables/";

//default differentials
double dt      = tstep_factor*Lx*kpc_to_cm/Nx/cl_factor/cl; //half a cell crossing time
double ray_dt  = tstep_factor*Lx*kpc_to_cm/Nx/cl_factor/cl;

//time
double t = 0;
double t_step	= 0.;
int    step		= 0;

//speed of light
double clight	[Nx][Ny][Nz];

//grid variables
double x		[Nx]; //x coordinate
double y		[Ny]; //y coordinate
double z		[Nz]; // z coordinate
double dx		[Nx]; //x differential
double dy		[Ny]; //y differential
double dz		[Nz]; //z differential

//number densities
double 	rho		[Nx][Ny][Nz];
double  nH		[Nx][Ny][Nz]; //number density of hydrogen
double  nHe		[Nx][Ny][Nz]; //number density of helium
double  nH1		[Nx][Ny][Nz]; //number density of neutral hydrogen
double  nHe1	[Nx][Ny][Nz]; //number denstiy of neutral helium
double  nH2		[Nx][Ny][Nz]; //number density of ionized hydrogen
double  nHe2	[Nx][Ny][Nz]; //number density of singly ionized helium
double  nHe3	[Nx][Ny][Nz]; //number density of doubly ionized helium
double  ne		[Nx][Ny][Nz]; //number density of free electrons
double  n_tot	[Nx][Ny][Nz]; //total number density

double nH1_prev	[Nx][Ny][Nz]; //number density of neutral hydrogen
double nHe1_prev[Nx][Ny][Nz]; //number denstiy of neutral helium
double nH2_prev	[Nx][Ny][Nz]; //number density of ionized hydrogen
double nHe2_prev[Nx][Ny][Nz]; //number density of singly ionized helium
double nHe3_prev[Nx][Ny][Nz]; //number density of doubly ionized helium

//time derivatives
double dne_dt	[Nx][Ny][Nz];

//abundances
double f_H1		 [Nx][Ny][Nz]; //HI fraction
double f_H1_step [Nx][Ny][Nz];
double f_He1	 [Nx][Ny][Nz]; //HeI fraction
double f_He1_step[Nx][Ny][Nz];
double f_H2		 [Nx][Ny][Nz]; //HII fraction
double f_H2_step [Nx][Ny][Nz];
double f_He2	 [Nx][Ny][Nz]; //HeII fraction
double f_He2_step[Nx][Ny][Nz];
double f_He3	 [Nx][Ny][Nz]; //HeIII fraction
double f_He3_step[Nx][Ny][Nz];

//thermal parameters
double temp	[Nx][Ny][Nz]; //temperature

//photoionization rates
double gamma_H1_tot	[Nx][Ny][Nz]; //photoionization rate of HI
double gamma_He1_tot[Nx][Ny][Nz]; //photoionization rate of HeI
double gamma_He2_tot[Nx][Ny][Nz]; //photionization rate of He2

//recombination coefficients
double recomb_H2 	[Nx][Ny][Nz];
double recomb_He2	[Nx][Ny][Nz];
double recomb_He3	[Nx][Ny][Nz];

//collisional ionization coefficients
double ci_H1		[Nx][Ny][Nz];
double ci_He1		[Nx][Ny][Nz];
double ci_He2		[Nx][Ny][Nz];

//heating/cooling rates
double heat_rate	[Nx][Ny][Nz]; //total heating rate
double cool_rate	[Nx][Ny][Nz]; //total cooling rate

//energy
double u_nu			[Nx][Ny][Nz];  //specific energy density of radiation

//rt equation
double gamma_nu_H1 [Nx][Ny][Nz]; //spectral absorption coefficient of HI
double gamma_nu_He1[Nx][Ny][Nz]; //spectral absorption coefficient of HeI
double gamma_nu_He2[Nx][Ny][Nz]; //spectral absorption coefficient of HeII
double gamma_nu_tot[Nx][Ny][Nz]; //total spectral absorption coefficient

/////////////////////////////////////////////////////////////////////////////////
//RT variables for ray tracing program

int total_rays = 0;

double source_lum[Nx][Ny][Nz]; //luminosity of source cells.  0 for cells that are not a source  

vector<vector<int>> ray_index[Nx][Ny][Nz];
vector<vector<double>> ray[Nx][Ny][Nz];

// Keep track of the propagation information for each ray
vector<vector<double>> ray_dep       [Nx][Ny][Nz];

//testing
double avg_xh1[Nx][Ny][Nz], avg_gamma[Nx][Ny][Nz];
double treion[Nx][Ny][Nz];









