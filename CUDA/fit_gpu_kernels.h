extern "C" __global__ void fit_PVM_single_kernel(	//INPUT
							const double* 		data, 
							const double* 		bvecs, 
							const double* 		bvals, 
							const int 		nvox, 
							const int 		nfib, 
							const bool 		m_include_f0, 
							//INPUT-OUTPUT
							double* 		params);

extern "C" __global__ void fit_PVM_single_c_kernel(	//INPUT
							const double* 		data, 
							const double* 		bvecs, 
							const double* 		bvals, 
							const int 		nvox, 
							const int 		nfib, 
							const bool		m_eval_BIC,
							const bool 		m_include_f0,
							const bool	 	m_return_fanning,
							//INPUT-OUTPUT
							double* 		params);

extern "C" __global__ void fit_PVM_multi_kernel(	//INPUT
							const double* 		data, 
							const double* 		params_PVM_simple_c,
							const double* 		bvecs, 
							const double* 		bvals, 
							const int 		nvox, 
							const int 		nfib, 		
							const bool 		m_include_f0,
							//OUTPUT
							double* 		params);

extern "C" __global__ void get_residuals_PVM_single_kernel(	//INPUT
								const double* 		data, 
								const double* 		params,
								const double* 		bvecs, 
								const double* 		bvals, 
								const int 		nvox, 
								const int 		nfib, 
								const bool 		m_include_f0,
								const bool* 		includes_f0,
								//OUTPUT
								double*			residuals);

extern "C" __global__ void get_residuals_PVM_single_c_kernel(	//INPUT
								const double* 		data, 
								const double* 		params,
								const double* 		bvecs, 
								const double* 		bvals, 
								const int 		nvox, 
								const int 		nfib, 
								const bool 		m_include_f0,
								const bool* 		includes_f0,
								//OUTPUT
								double*			residuals);


extern "C" __global__ void get_residuals_PVM_multi_kernel(	//INPUT
								const double* 		data, 
								const double* 		params,
								const double* 		bvecs, 
								const double* 		bvals, 
								const int 		nvox, 
								const int 		nfib, 
								const bool 		m_include_f0,
								const bool* 		includes_f0,
								//OUTPUT
								double*			residuals);
