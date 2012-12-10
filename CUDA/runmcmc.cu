/*  runmcmc.cu

    Tim Behrens, Saad Jbabdi, Stam Sotiropoulos, Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */

#include "xfibresoptions.h"
#include <curand.h>
#include "runmcmc_kernels.cu"
#include "sync_check.h"

#include <host_vector.h>
#include <device_vector.h> 

#include <time.h>
#include <sys/time.h>
#include "init_gpu.h"

using namespace Xfibres;

////////////////////////////////////////////////////// 
//   MCMC IN GPU
////////////////////////////////////////////////////// 

void init_Fibres_Multifibres(	//INPUT
				thrust::device_vector<double> 			datam_gpu,
				thrust::device_vector<double> 			params_gpu,
				thrust::device_vector<float> 			tau_gpu,
				thrust::device_vector<double> 			bvals_gpu,
				thrust::device_vector<double> 			alpha_gpu,
				thrust::device_vector<double> 			beta_gpu,
				string 						output_file, 
				//OUTPUT
				thrust::device_vector<FibreGPU>& 		fibres_gpu,
				thrust::device_vector<MultifibreGPU>& 		multifibres_gpu,
				thrust::device_vector<double>&			signals_gpu,
				thrust::device_vector<double>&			isosignals_gpu)
{
	std::ofstream myfile;
	myfile.open (output_file.data(), ios::out | ios::app );
	myfile << "----------------------------------------------------- " << "\n"; 
   	myfile << "------MCMC ALGORITHM PART INITIALITATION IN GPU ----- " << "\n"; 
   	myfile << "----------------------------------------------------- " << "\n"; 	

   	struct timeval t1,t2;
   	double time;
   	gettimeofday(&t1,NULL);

	int nvox = multifibres_gpu.size();

	xfibresOptions& opts = xfibresOptions::getInstance();
	int nfib= opts.nfibres.value();
	int nparams_fit = 2+3*opts.nfibres.value();
	if(opts.modelnum.value()==2) nparams_fit++;
	if(opts.f0.value()) nparams_fit++;

	int blocks = nvox; 
  	dim3 Dim_Grid_MCMC(blocks, 1);
  	dim3 Dim_Block_MCMC(THREADS_BLOCK ,1);	///dimensions for MCMC

	double *datam_ptr = thrust::raw_pointer_cast(datam_gpu.data());
	double *params_ptr = thrust::raw_pointer_cast(params_gpu.data());	
	float *tau_ptr = thrust::raw_pointer_cast(tau_gpu.data());	
	double *bvals_ptr = thrust::raw_pointer_cast(bvals_gpu.data());
	double *alpha_ptr = thrust::raw_pointer_cast(alpha_gpu.data());
	double *beta_ptr = thrust::raw_pointer_cast(beta_gpu.data());
	FibreGPU *fibres_ptr =  thrust::raw_pointer_cast(fibres_gpu.data());
	MultifibreGPU *multifibres_ptr = thrust::raw_pointer_cast(multifibres_gpu.data());
	double *signals_ptr = thrust::raw_pointer_cast(signals_gpu.data());
	double *isosignals_ptr = thrust::raw_pointer_cast(isosignals_gpu.data());

	init_Fibres_Multifibres_kernel<<< Dim_Grid_MCMC, Dim_Block_MCMC>>>(datam_ptr, params_ptr, tau_ptr, bvals_ptr, alpha_ptr, beta_ptr, nvox, nfib, nparams_fit, opts.modelnum.value(), opts.fudge.value(), opts.f0.value(), opts.rician.value(), opts.ardf0.value(), opts.all_ard.value(), opts.no_ard.value(), fibres_ptr, multifibres_ptr, signals_ptr, isosignals_ptr);
	sync_check("init_Fibres_Multifibres_kernel");

	gettimeofday(&t2,NULL);
    	time=timeval_diff(&t2,&t1);
   	myfile << "TIME TOTAL: " << time << " seconds\n"; 
	myfile << "--------------------------------------------" << "\n\n" ; 
}

void runmcmc_burnin(	//INPUT
			thrust::device_vector<double> 			datam_gpu,
			thrust::device_vector<double> 			bvals_gpu,
			thrust::device_vector<double> 			alpha_gpu,
			thrust::device_vector<double> 			beta_gpu,
			double 						seed,
			string 						output_file, 
			//INPUT-OUTPUT
			thrust::device_vector<FibreGPU>& 		fibres_gpu,
			thrust::device_vector<MultifibreGPU>& 		multifibres_gpu,
			thrust::device_vector<double>&			signals_gpu,
			thrust::device_vector<double>&			isosignals_gpu)
{
	xfibresOptions& opts = xfibresOptions::getInstance();
	
	std::ofstream myfile;
	myfile.open (output_file.data(), ios::out | ios::app );
	myfile << "----------------------------------------------------- " << "\n"; 
   	myfile << "--------- MCMC ALGORITHM PART BURNIN IN GPU --------- " << "\n"; 
   	myfile << "----------------------------------------------------- " << "\n"; 	

   	struct timeval t1,t2,t_tot1,t_tot2;
   	double time,timecurand,timemcmc;
   	time=0;
   	timecurand=0;
   	timemcmc=0;

   	gettimeofday(&t_tot1,NULL);

   	size_t free,total;
	
	int nvox = multifibres_gpu.size();
   	int nfib= opts.nfibres.value();
	int nparams;

	if(opts.f0.value()) nparams=3+nfib*3;
	else nparams=2+nfib*3;	
	if(opts.modelnum.value()==2) nparams++;
	if(opts.rician.value()) nparams++;
   
   	unsigned int totalrandoms=(opts.nburn.value() * nvox * nparams);

   	cuMemGetInfo(&free,&total);
   	myfile << "Free memory Before Randoms: "<< free <<  " ---- Total memory: " << total << "\n";
   	//4 bytes each float, 2 randoms arrays, and 80% of total memory at this moment 
   	unsigned int maxrandoms=(free/(4*2))*0.8; 

   	myfile << "Total randoms: " << totalrandoms << "\n"; 
   	myfile << "Max randoms: " << maxrandoms << "\n"; 
   
   	int steps; //num iter if not enough memory
   	int minrandoms; //min num of randoms ensamble
   	minrandoms= nvox * nparams;

   	int iters_step=0;
	int nrandoms=0;	

   	if(totalrandoms>maxrandoms){ 
		iters_step = maxrandoms / minrandoms; 
		nrandoms = iters_step*minrandoms;		//nrandoms for each step
		steps =  (opts.nburn.value()/iters_step);  	//repeat process iterations times, no enough memory for all randoms 			
   	}else{ 
		nrandoms = totalrandoms;
		iters_step= opts.nburn.value();
		steps = 0;
  	 }   

   	myfile << "Process " << opts.nburn.value() << " iterations divided in "<< steps << " steps with "<< iters_step << " iterations in each one" << "\n";    

   	int last_step = opts.nburn.value() - (iters_step*steps);
   	int last_randoms= (last_step*minrandoms); 

   	myfile << "Last step with " << last_step << " iterations" << "\n"; 
	
	thrust::device_vector<float> randomsN_gpu;
	thrust::device_vector<float> randomsU_gpu;	
	randomsN_gpu.resize(nrandoms);
	randomsU_gpu.resize(nrandoms);

   	cuMemGetInfo(&free,&total);
   	myfile << "Free memory after Malloc Randoms: "<< free <<  " ---- Total memory: " << total << "\n";
   
  	int blocks = nvox;        
  	dim3 Dim_Grid(blocks, 1);
  	dim3 Dim_Block(THREADS_BLOCK,1);	//dimensions for MCMC   

   	myfile << "\n" << "NUM BLOCKS: " << blocks << "\n"; 
   	myfile << "THREADS PER BLOCK : " << THREADS_BLOCK << "\n\n"; 	

   	curandGenerator_t gen;
   	curandCreateGenerator(&gen,CURAND_RNG_PSEUDO_DEFAULT);
   	curandSetPseudoRandomGeneratorSeed(gen,seed);

	//get pointers
	double *datam_ptr = thrust::raw_pointer_cast(datam_gpu.data());
	double *bvals_ptr = thrust::raw_pointer_cast(bvals_gpu.data());
	double *alpha_ptr = thrust::raw_pointer_cast(alpha_gpu.data());
	double *beta_ptr = thrust::raw_pointer_cast(beta_gpu.data());
	float *randomsN_ptr = thrust::raw_pointer_cast(randomsN_gpu.data());
	float *randomsU_ptr = thrust::raw_pointer_cast(randomsU_gpu.data());
	FibreGPU *fibres_ptr =  thrust::raw_pointer_cast(fibres_gpu.data());
	MultifibreGPU *multifibres_ptr = thrust::raw_pointer_cast(multifibres_gpu.data());
	double *signals_ptr = thrust::raw_pointer_cast(signals_gpu.data());
	double *isosignals_ptr = thrust::raw_pointer_cast(isosignals_gpu.data());
	
   	for(int i=0;i<steps;i++){

   		gettimeofday(&t1,NULL);

	   	curandGenerateNormal(gen,randomsN_ptr,nrandoms,0,1);
	   	curandGenerateUniform(gen,randomsU_ptr,nrandoms);	//generate randoms

 	   	gettimeofday(&t2,NULL);
    	   	timecurand+=timeval_diff(&t2,&t1);

	   	gettimeofday(&t1,NULL);

	   	runmcmc_burnin_kernel<<< Dim_Grid, Dim_Block >>>(datam_ptr, bvals_ptr, alpha_ptr, beta_ptr, randomsN_ptr, randomsU_ptr, nvox, opts.modelnum.value(), opts.fudge.value(), opts.f0.value(), opts.ardf0.value(), !opts.no_ard.value(), opts.rician.value(), opts.updateproposalevery.value(), iters_step, (i*iters_step), fibres_ptr, multifibres_ptr, signals_ptr, isosignals_ptr);   
	   	sync_check("runmcmc_burnin_kernel");

 	   	gettimeofday(&t2,NULL);
    	   	timemcmc+=timeval_diff(&t2,&t1);
   	}

   	gettimeofday(&t1,NULL);

   	if(nvox!=0){
   		curandGenerateNormal(gen,randomsN_ptr,last_randoms,0,1);
   		curandGenerateUniform(gen,randomsU_ptr,last_randoms); 	//generate randoms
   	}
	
   	gettimeofday(&t2,NULL);
   	timecurand+=timeval_diff(&t2,&t1);

   	gettimeofday(&t1,NULL);

   	if(nvox!=0){
		runmcmc_burnin_kernel<<< Dim_Grid, Dim_Block >>>(datam_ptr, bvals_ptr, alpha_ptr, beta_ptr, randomsN_ptr, randomsU_ptr, nvox, opts.modelnum.value(), opts.fudge.value(), opts.f0.value(), opts.ardf0.value(), !opts.no_ard.value(), opts.rician.value(), opts.updateproposalevery.value(), last_step, (steps*iters_step), fibres_ptr, multifibres_ptr, signals_ptr, isosignals_ptr);   
   		sync_check("runmcmc_burnin_kernel");
   	}

   	gettimeofday(&t2,NULL);
   	timemcmc+=timeval_diff(&t2,&t1);


    	myfile << "TIME CURAND: " << timecurand << " seconds\n"; 
    	myfile << "TIME RUNMCMC: " << timemcmc << " seconds\n"; 
   
   	curandDestroyGenerator(gen);

	gettimeofday(&t_tot2,NULL);
    	time=timeval_diff(&t_tot2,&t_tot1);
   	myfile << "TIME TOTAL: " << time << " seconds\n"; 
	myfile << "--------------------------------------------" << "\n\n" ; 
	myfile.close();

   	sync_check("runmcmc_burnin");
}


void runmcmc_record(	//INPUT
			thrust::device_vector<double> 			datam_gpu,
			thrust::device_vector<double> 			bvals_gpu,
			thrust::device_vector<double> 			alpha_gpu,
			thrust::device_vector<double> 			beta_gpu,
			thrust::device_vector<FibreGPU> 		fibres_gpu,
			thrust::device_vector<MultifibreGPU> 		multifibres_gpu,
			thrust::device_vector<double>			signals_gpu,
			thrust::device_vector<double>			isosignals_gpu,
			double 						seed,
			string 						output_file, 
			//OUTPUT
			thrust::device_vector<int>&			multirecords_gpu,
			thrust::device_vector<float>&			rf0_gpu,
			thrust::device_vector<float>&			rtau_gpu,
			thrust::device_vector<float>&			rs0_gpu,
			thrust::device_vector<float>&			rd_gpu,
			thrust::device_vector<float>&			rdstd_gpu,
			thrust::device_vector<float>&			rth_gpu,
			thrust::device_vector<float>&			rph_gpu,
			thrust::device_vector<float>&			rf_gpu,
			thrust::device_vector<float>&			rlikelihood_energy_gpu)
{
	xfibresOptions& opts = xfibresOptions::getInstance();
	
	std::ofstream myfile;
	myfile.open (output_file.data(), ios::out | ios::app );
	myfile << "----------------------------------------------------- " << "\n"; 
   	myfile << "--------- MCMC ALGORITHM PART RECORD IN GPU --------- " << "\n"; 
   	myfile << "----------------------------------------------------- " << "\n"; 	

   	struct timeval t1,t2,t_tot1,t_tot2;
   	double time,timecurand,timemcmc;
   	time=0;
   	timecurand=0;
   	timemcmc=0;

   	gettimeofday(&t_tot1,NULL);

   	size_t free,total;

	int totalrecords = (opts.njumps.value()/opts.sampleevery.value()); 
	
	int nvox = multifibres_gpu.size();
   	int nfib= opts.nfibres.value();
	int nparams;

	if(opts.f0.value()) nparams=3+nfib*3;
	else nparams=2+nfib*3;	
	if(opts.modelnum.value()==2) nparams++;
	if(opts.rician.value()) nparams++;
   
   	unsigned int totalrandoms=(opts.njumps.value() * nvox * nparams);

   	cuMemGetInfo(&free,&total);
   	myfile << "Free memory Before Randoms: "<< free <<  " ---- Total memory: " << total << "\n";
   	//4 bytes each float, 2 randoms arrays, and 80% of total memory at this moment 
   	unsigned int maxrandoms=(free/(4*2))*0.8; 

   	myfile << "Total randoms: " << totalrandoms << "\n"; 
   	myfile << "Max randoms: " << maxrandoms << "\n"; 
   
   	int steps; //num iter if not enough memory
   	int minrandoms; //min num of randoms ensamble
   	minrandoms= nvox * nparams;

   	int iters_step=0;
	int nrandoms=0;	

   	if(totalrandoms>maxrandoms){ 
		iters_step = maxrandoms / minrandoms; 
		nrandoms = iters_step*minrandoms;		//nrandoms for each step
		steps =  (opts.njumps.value()/iters_step);  	//repeat process iterations times, no enough memory for all randoms 			
   	}else{ 
		nrandoms = totalrandoms;
		iters_step= opts.njumps.value();
		steps = 0;
  	 }   

   	myfile << "Process " << opts.njumps.value() << " iterations divided in "<< steps << " steps with "<< iters_step << " iterations in each one" << "\n";    

   	int last_step = opts.njumps.value() - (iters_step*steps);
   	int last_randoms= (last_step*minrandoms); 

   	myfile << "Last step with " << last_step << " iterations" << "\n"; 
	
	thrust::device_vector<float> randomsN_gpu;
	thrust::device_vector<float> randomsU_gpu;	
	randomsN_gpu.resize(nrandoms);
	randomsU_gpu.resize(nrandoms);

   	cuMemGetInfo(&free,&total);
   	myfile << "Free memory after Malloc Randoms: "<< free <<  " ---- Total memory: " << total << "\n";
   
  	int blocks = nvox;        
  	dim3 Dim_Grid(blocks, 1);
  	dim3 Dim_Block(THREADS_BLOCK,1);	//dimensions for MCMC   

   	myfile << "\n" << "NUM BLOCKS: " << blocks << "\n"; 
   	myfile << "THREADS PER BLOCK : " << THREADS_BLOCK << "\n\n"; 	

   	curandGenerator_t gen;
   	curandCreateGenerator(&gen,CURAND_RNG_PSEUDO_DEFAULT);
   	curandSetPseudoRandomGeneratorSeed(gen,seed);

	//get pointers
	double *datam_ptr = thrust::raw_pointer_cast(datam_gpu.data());
	double *bvals_ptr = thrust::raw_pointer_cast(bvals_gpu.data());
	double *alpha_ptr = thrust::raw_pointer_cast(alpha_gpu.data());
	double *beta_ptr = thrust::raw_pointer_cast(beta_gpu.data());
	float *randomsN_ptr = thrust::raw_pointer_cast(randomsN_gpu.data());
	float *randomsU_ptr = thrust::raw_pointer_cast(randomsU_gpu.data());
	FibreGPU *fibres_ptr =  thrust::raw_pointer_cast(fibres_gpu.data());
	MultifibreGPU *multifibres_ptr = thrust::raw_pointer_cast(multifibres_gpu.data());
	double *signals_ptr = thrust::raw_pointer_cast(signals_gpu.data());
	double *isosignals_ptr = thrust::raw_pointer_cast(isosignals_gpu.data());
	
	int *multirecords_ptr = thrust::raw_pointer_cast(multirecords_gpu.data());
	float *rf0_ptr = thrust::raw_pointer_cast(rf0_gpu.data());
	float *rtau_ptr = thrust::raw_pointer_cast(rtau_gpu.data());
	float *rs0_ptr = thrust::raw_pointer_cast(rs0_gpu.data());
	float *rd_ptr = thrust::raw_pointer_cast(rd_gpu.data());
	float *rdstd_ptr = thrust::raw_pointer_cast(rdstd_gpu.data());
	float *rth_ptr = thrust::raw_pointer_cast(rth_gpu.data());
	float *rph_ptr = thrust::raw_pointer_cast(rph_gpu.data());
	float *rf_ptr = thrust::raw_pointer_cast(rf_gpu.data());
	float *rlikelihood_energy_ptr = thrust::raw_pointer_cast(rlikelihood_energy_gpu.data());

   	for(int i=0;i<steps;i++){

   		gettimeofday(&t1,NULL);

	   	curandGenerateNormal(gen,randomsN_ptr,nrandoms,0,1);
	   	curandGenerateUniform(gen,randomsU_ptr,nrandoms);	//generate randoms

 	   	gettimeofday(&t2,NULL);
    	   	timecurand+=timeval_diff(&t2,&t1);

	   	gettimeofday(&t1,NULL);

	   	runmcmc_record_kernel<<< Dim_Grid, Dim_Block >>>(datam_ptr, bvals_ptr, alpha_ptr, beta_ptr, fibres_ptr, multifibres_ptr, signals_ptr, isosignals_ptr, randomsN_ptr, randomsU_ptr, nvox, opts.modelnum.value(), opts.fudge.value(), opts.f0.value(), opts.ardf0.value(), !opts.no_ard.value(), opts.rician.value(), opts.updateproposalevery.value(), iters_step, (i*iters_step), opts.nburn.value(), opts.sampleevery.value(), totalrecords, multirecords_ptr, rf0_ptr, rtau_ptr, rs0_ptr, rd_ptr, rdstd_ptr, rth_ptr, rph_ptr, rf_ptr, rlikelihood_energy_ptr);   
	   	sync_check("runmcmc_record_kernel");

 	   	gettimeofday(&t2,NULL);
    	   	timemcmc+=timeval_diff(&t2,&t1);
   	}

   	gettimeofday(&t1,NULL);

   	if(nvox!=0){
   		curandGenerateNormal(gen,randomsN_ptr,last_randoms,0,1);
   		curandGenerateUniform(gen,randomsU_ptr,last_randoms); 	//generate randoms
   	}
	
   	gettimeofday(&t2,NULL);
   	timecurand+=timeval_diff(&t2,&t1);

   	gettimeofday(&t1,NULL);

   	if(nvox!=0){
		runmcmc_record_kernel<<< Dim_Grid, Dim_Block >>>(datam_ptr, bvals_ptr, alpha_ptr, beta_ptr, fibres_ptr, multifibres_ptr, signals_ptr, isosignals_ptr, randomsN_ptr, randomsU_ptr, nvox, opts.modelnum.value(), opts.fudge.value(), opts.f0.value(), opts.ardf0.value(), !opts.no_ard.value(), opts.rician.value(), opts.updateproposalevery.value(), last_step, (steps*iters_step), opts.nburn.value(), opts.sampleevery.value(), totalrecords, multirecords_ptr, rf0_ptr, rtau_ptr, rs0_ptr, rd_ptr, rdstd_ptr, rth_ptr, rph_ptr, rf_ptr, rlikelihood_energy_ptr);   
   		sync_check("runmcmc_record_kernel");
   	}

   	gettimeofday(&t2,NULL);
   	timemcmc+=timeval_diff(&t2,&t1);


    	myfile << "TIME CURAND: " << timecurand << " seconds\n"; 
    	myfile << "TIME RUNMCMC: " << timemcmc << " seconds\n"; 
   
   	curandDestroyGenerator(gen);

	gettimeofday(&t_tot2,NULL);
    	time=timeval_diff(&t_tot2,&t_tot1);
   	myfile << "TIME TOTAL: " << time << " seconds\n"; 
	myfile << "--------------------------------------------" << "\n\n" ; 
	myfile.close();

   	sync_check("runmcmc_record");
}
