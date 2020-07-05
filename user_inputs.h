#include "rates.cc"

#define FALSE 0
#define TRUE 1

char const rec_case[]   = "B"; //Use case "A" or "B" recombination coefficients

int const input_grid    = TRUE; //User-defined input grid flag

int const recomb_cool   = TRUE; //Recombination cooling flag
int const coll_ion      = TRUE; //Collisional ionization flag
int const coll_exc_cool	= TRUE; //Collisional excitation cooling flag  
int const compton       = FALSE; //compton heating/cooling flag
int const temp_ev       = FALSE; //temperature evolution flag

//input files (if user-defined)
const char start_file[]      = "input_files/gas_test_256_post.txt";
const char ray_file[]        = "input_files/ray_test_256_post.txt";
const char source_field[]    = "input_files/nbody_galaxies_Mmin2e10_z5.8";
const char ugas_file[]       = "input_files/ugas_256_z=5.6";

//output files
const char gas_output[]      = "output_files/gas_test_256_post.txt";
const char ray_output[]      = "output_files/ray_test_256 post.txt";
const char otf_output[]      = "output_files/otf_test_256_post.txt";
const char tre_output[]      = "output_files/treion_256_post.txt";

//grid sizess
int const Nx      = 256; //Number of grid points in the x direction
int const Ny	  = 256; //Number of grid points in the y direction
int const Nz      = 256; //Number of grid points in the z direction

//time stepping
int const    hpx_lvl      = 0; //healpix level for ray casting/merging
double const tstep_factor = 1.0; //fraction of light crossing time per time step
double const cl_factor    = 1.0; //simulation speed of light / true speed of light
int const itercount       = 5; //number of iterations for the photo-ionization rate solver

double const zz       = 5.6; //redshift

double const fesc     = 0.1; //escape fraction for sources
double const Lx       = 200e3/hh/(1. + zz); //Box length in the x direction (in kpc)
double const Ly       = 200e3/hh/(1. + zz); //Box length in the y direction
double const Lz       = 200e3/hh/(1. + zz); //Box length in the z direction
double const t_max    = 10.0;    //Runtime (in Myr)
double const nu_phot  = 1.3*13.6/h_eV; //Photon frequency (in Hz)
double const temp_0   = 1e2; //Initial gas temperature/temperature of neutral gas

double const fH1_0      = 1. - 1.2e-3; //Initial HI fraction
double const fHe1_0     = 1.; //Initial HeI fraction
double const fHe2_0     = 0; //Initial HeII fraction

