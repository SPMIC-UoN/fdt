/*  fit_gpu_kernels.h

    Tim Behrens, Saad Jbabdi, Stam Sotiropoulos, Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */

extern "C" __global__ void fit_PVM_single_kernel(	//INPUT
							const double* 		data, 
							const double* 		bvecs, 
							const double* 		bvals, 
							const int 		nvox, 
							const int		ndirections,
							const int 		nfib,
							const int		nparams, 
							const bool 		m_include_f0, 
							const bool		gradnonlin,
							//INPUT-OUTPUT
							double* 		params);

extern "C" __global__ void fit_PVM_single_c_kernel(	//INPUT
							const double* 		data, 
							const double* 		bvecs, 
							const double* 		bvals, 
							const int 		nvox, 
							const int		ndirections,
							const int 		nfib, 
							const int		nparams,
							const bool		m_eval_BIC,
							const bool 		m_include_f0,
							const bool	 	m_return_fanning,
							const bool		gradnonlin,
							//INPUT-OUTPUT
							double* 		params);

extern "C" __global__ void fit_PVM_multi_kernel(	//INPUT
							const double* 		data, 
							const double* 		params_PVM_simple_c,
							const double* 		bvecs, 
							const double* 		bvals, 
							const int 		nvox, 
							const int		ndirections,
							const int 		nfib, 	
							const int		nparams,	
							const bool 		m_include_f0,
							const bool		gradnonlin,
							//OUTPUT
							double* 		params);

extern "C" __global__ void get_residuals_PVM_single_kernel(	//INPUT
								const double* 		data, 
								const double* 		params,
								const double* 		bvecs, 
								const double* 		bvals, 
								const int 		nvox, 
								const int		ndirections,
								const int 		nfib, 
								const int		nparams,
								const bool 		m_include_f0,
								const bool		gradnonlin,
								const bool* 		includes_f0,								
								//OUTPUT
								double*			residuals);

extern "C" __global__ void get_residuals_PVM_single_c_kernel(	//INPUT
								const double* 		data, 
								const double* 		params,
								const double* 		bvecs, 
								const double* 		bvals, 
								const int 		nvox, 
								const int		ndirections,
								const int 		nfib, 
								const int		nparams,
								const bool 		m_include_f0,
								const bool		gradnonlin,
								const bool* 		includes_f0,								
								//OUTPUT
								double*			residuals);


extern "C" __global__ void get_residuals_PVM_multi_kernel(	//INPUT
								const double* 		data, 
								const double* 		params,
								const double* 		bvecs, 
								const double* 		bvals, 
								const int 		nvox, 
								const int		ndirections,
								const int 		nfib, 
								const int		nparams,
								const bool 		m_include_f0,
								const bool		gradnonlin,
								const bool* 		includes_f0,								
								//OUTPUT
								double*			residuals);

