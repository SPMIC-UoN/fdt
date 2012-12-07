
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
				thrust::device_vector<double>&			isosignals_gpu);

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
			thrust::device_vector<double>&			isosignals_gpu);


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
			thrust::device_vector<float>&			rlikelihood_energy_gpu);
