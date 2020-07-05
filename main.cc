#include <stdio.h>
#include <time.h>
#include <math.h>
#include "rt_funcs.cc"

int main(void)
{
	clock_t start, end;
	clock_t rt_start, rt_end;
	clock_t iter_start, iter_end;
	double rt_time_used = 0.; 
	double iter_time_used = 0.;
	double cpu_time_used;

	init_rand();

	start = clock();
		
	init();
	printf("Grid initialized.\n");
	printf("Gas initialized.\n");
	update_gamma_nu();    
	set_clight();
	init_rt();
	printf("RT data initialized.\n");
	
	make_output();

	end   = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("%le\n",cpu_time_used);
	
	double t_tot = t_max*yr_to_s*1e6;
	
	while (t < t_tot)  {

		step += 1;
		
		rt_start = clock();
		
		//advance/split rays
		update_rays();
		
		rt_end = clock();
		rt_time_used += ((double) (rt_end - rt_start)) / CLOCKS_PER_SEC;
		
		iter_start = clock();
		
		//update energy
		find_avg_unu();
		
		iter_end = clock();
		iter_time_used += ((double) (iter_end - iter_start)) / CLOCKS_PER_SEC;
		
		rt_start = clock();
		
		//merge rays
		merge();
		
		rt_end = clock();
		rt_time_used += ((double) (rt_end - rt_start)) / CLOCKS_PER_SEC;
		
		//update absorption coefficients
		update_gamma(); //photoionization rates
		
		//update heating and cooling functions
		if ( (temp_ev == TRUE) )  {
			update_heat_cool();
		}
		
		//update chemistry equations
		update_chem();
		
		//update attenuation coefficient
		update_gamma_nu();

		//update heating/cooling rates and then temperature
		update_thermal();
	
		//update time
		t += dt;
	
		printf("Writing on-the-fly output\n");
		write_otf();
		update_step();
			
		//update time steps
		update_dt();
		//testing
		update_treion();
		
		printf("time: ");
		printf("%le\n",t/yr_to_s/1e6);
		
	}
	//write output
	write_gas_binary();
	write_rays_binary();
	write_treion();
	
	end   = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("%le\n",t/yr_to_s/1e6);
	printf("Total time run: ");
	printf("%le\n",cpu_time_used);
	printf("RT time: ");
	printf("%le\n", rt_time_used);
	printf("Iteration time: ");
	printf("%le\n", iter_time_used);
}

