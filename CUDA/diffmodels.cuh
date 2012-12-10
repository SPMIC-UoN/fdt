/*  diffmodels.cuh

    Tim Behrens, Saad Jbabdi, Stam Sotiropoulos, Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */

#include <device_vector.h>

void fit_PVM_single(	//INPUT
			const vector<ColumnVector> 	datam_vec, 
			const vector<Matrix> 		bvecs_vec,
			const vector<Matrix> 		bvals_vec,
			thrust::device_vector<double> 	datam_gpu, 
			thrust::device_vector<double>	bvecs_gpu, 
			thrust::device_vector<double>	bvals_gpu,	
			bool 				m_include_f0,		
			//OUTPUT
			thrust::device_vector<double>&	params_gpu);

void fit_PVM_single_c(	//INPUT
			const vector<ColumnVector> 	datam_vec, 
			const vector<Matrix> 		bvecs_vec,
			const vector<Matrix> 		bvals_vec,
			thrust::device_vector<double> 	datam_gpu, 
			thrust::device_vector<double>	bvecs_gpu, 
			thrust::device_vector<double>	bvals_gpu,	
			bool 				m_include_f0,		
			//OUTPUT
			thrust::device_vector<double>&	params_gpu);

void fit_PVM_multi(	//INPUT
			thrust::device_vector<double> 	datam_gpu, 
			thrust::device_vector<double>	bvecs_gpu, 
			thrust::device_vector<double>	bvals_gpu,	
			int 				nvox,		
			bool 				m_include_f0,
			//OUTPUT
			thrust::device_vector<double>&	params_gpu);

void calculate_tau(	//INPUT
			thrust::device_vector<double> 	datam_gpu, 
			thrust::device_vector<double>&	params_gpu,
			thrust::device_vector<double>	bvecs_gpu, 
			thrust::device_vector<double>	bvals_gpu,
			thrust::host_vector<int>&	vox_repeat,
			int				nrepeat,		
			string 				output_file,	
			//OUTPUT
			thrust::host_vector<float>&	tau);


__device__ double cf_PVM_single(	//INPUT
					const double*			params,
					const double*			data,
					const double*			bvecs, 
					const double*			bvals,
					const int 			nparams, 
					const bool 			m_include_f0);

__device__ void grad_PVM_single(	//INPUT
					const double*			params,
					const double*			data,
					const double*			bvecs, 
					const double*			bvals,
					const int 			nparams,
					const bool 			m_include_f0,
					// OUTPUT
					double*				grad);

__device__ void hess_PVM_single(	//INPUT
					const double*			params,
					const double*			bvecs, 
					const double*			bvals,
					const int 			nparams,
					const bool 			m_include_f0,
					//OUTPUT
					double*				hess);

__device__ double cf_PVM_single_c(	//INPUT
					const double*			params,
					const double*			data,
					const double*			bvecs, 
					const double*			bvals,
					const int 			nparams, 
					const bool 			m_include_f0);

__device__ void grad_PVM_single_c(	//INPUT
					const double*			params,
					const double*			data,
					const double*			bvecs, 
					const double*			bvals,
					const int 			nparams,
					const bool 			m_include_f0,
					//OUTPUT
					double*				grad);

__device__ void hess_PVM_single_c(	//INPUT
					const double*			params,
					const double*			bvecs, 
					const double*			bvals,
					const int 			nparams,
					const bool 			m_include_f0,
					//OUTPUT
					double*				hess);

__device__ double cf_PVM_multi(		//INPUT
					const double*			params,
					const double*			data,
					const double*			bvecs, 
					const double*			bvals,
					const int 			nparams, 
					const bool 			m_include_f0);

__device__ void grad_PVM_multi(		//INPUT
					const double*			params,
					const double*			data,
					const double*			bvecs, 
					const double*			bvals,
					const int 			nparams,
					const bool 			m_include_f0,
					//OUTPUT
					double*				grad);

__device__ void hess_PVM_multi(		//INPUT
					const double*			params,
					const double*			bvecs, 
					const double*			bvals,
					const int 			nparams,
					const bool 			m_include_f0,
					//OUTPUT
					double*				hess);

