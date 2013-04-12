/*  diffmodels.cu

    Tim Behrens, Saad Jbabdi, Stam Sotiropoulos, Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */

#include "newmat.h"
#include "newimage/newimageall.h"
#include "xfibresoptions.h"
#include "diffmodels.h"

#include "PVM_single.cu"
#include "PVM_single_c.cu"
#include "PVM_multi.cu"
#include "fit_gpu_kernels.h"
#include "sync_check.h"

#include <time.h>
#include <sys/time.h>
#include "init_gpu.h"

using namespace Xfibres;

////////////////////////////////////////////////////// 
//   FIT IN GPU
////////////////////////////////////////////////////// 
void fit_PVM_single(	//INPUT
			const vector<ColumnVector> 	datam_vec, 
			const vector<Matrix> 		bvecs_vec,
			const vector<Matrix> 		bvals_vec,
			thrust::device_vector<double> 	datam_gpu, 
			thrust::device_vector<double>	bvecs_gpu, 
			thrust::device_vector<double>	bvals_gpu,
			int				ndirections,	
			bool 				m_include_f0,
			string 				output_file,		
			//OUTPUT
			thrust::device_vector<double>&	params_gpu)
{
	std::ofstream myfile;
	myfile.open (output_file.data(), ios::out | ios::app );

	xfibresOptions& opts = xfibresOptions::getInstance();
	int nvox = datam_vec.size();
	int nfib = opts.nfibres.value();
	int nparams;
	if (m_include_f0)
      		nparams = nfib*3 + 3; 
    	else
      		nparams = nfib*3 + 2;

	thrust::host_vector<double> params_host;
	params_host.resize(nvox*nparams);
	
	for(int vox=0;vox<nvox;vox++){
		// initialise with a tensor
  		DTI dti(datam_vec[vox],bvecs_vec[vox],bvals_vec[vox]);
  		dti.linfit();

  		// set starting parameters for nonlinear fitting
  		float _th,_ph;
  		cart2sph(dti.get_v1(),_th,_ph);

  		params_host[vox*nparams] = dti.get_s0();
  		//start(2) = dti.get_md()>0?dti.get_md()*2:0.001; // empirically found that d~2*MD
  		params_host[vox*nparams+1] = dti.get_l1()>0?dti.get_l1():0.002; // empirically found that d~L1
	  	params_host[vox*nparams+2] = dti.get_fa()<1?f2x(dti.get_fa()):f2x(0.95); // first pvf = FA 
	  	params_host[vox*nparams+3] = _th;
	  	params_host[vox*nparams+4] = _ph;
	  	float sumf=x2f(params_host[vox*nparams+2]);
	  	float tmpsumf=sumf;
	  	for(int ii=2,i=5;ii<=nfib;ii++,i+=3){
		    	float denom=2;
		    	do{
		      		params_host[vox*nparams+i] = f2x(x2f(params_host[vox*nparams+i-3])/denom);
		      		denom *= 2;
		      		tmpsumf = sumf + x2f(params_host[vox*nparams+i]);
		    	}while(tmpsumf>=1);
		    	sumf += x2f(params_host[vox*nparams+i]);
		    	cart2sph(dti.get_v(ii),_th,_ph);
		    	params_host[vox*nparams+i+1] = _th;
		    	params_host[vox*nparams+i+2] = _ph;
	  	}
	  	if (m_include_f0)
	    		params_host[vox*nparams+nparams-1]=f2x(FSMALL);
	}

	thrust::copy(params_host.begin(), params_host.end(), params_gpu.begin());	

	int blocks = nvox;
   	dim3 Dim_Grid(blocks,1);
  	dim3 Dim_Block(THREADS_BLOCK_FIT,1);

	int amount_shared = (THREADS_BLOCK_FIT+5*nparams+2*nparams*nparams+4*nfib+8)*sizeof(double)+(2+nparams)*sizeof(int);

	myfile << "Shared Memory Used in fit_PVM_single: " << amount_shared << "\n"; 

	fit_PVM_single_kernel<<<Dim_Grid, Dim_Block, amount_shared>>>(thrust::raw_pointer_cast(datam_gpu.data()), thrust::raw_pointer_cast(bvecs_gpu.data()), thrust::raw_pointer_cast(bvals_gpu.data()) ,nvox, ndirections, nfib, nparams, m_include_f0, thrust::raw_pointer_cast(params_gpu.data()));
	sync_check("fit_PVM_single_kernel");

	myfile.close();
}

void fit_PVM_single_c(	//INPUT
			const vector<ColumnVector> 	datam_vec, 
			const vector<Matrix> 		bvecs_vec,
			const vector<Matrix> 		bvals_vec,
			thrust::device_vector<double> 	datam_gpu, 
			thrust::device_vector<double>	bvecs_gpu, 
			thrust::device_vector<double>	bvals_gpu,
			int				ndirections,		
			bool 				m_include_f0,	
			string 				output_file,	
			//OUTPUT
			thrust::device_vector<double>&	params_gpu)
{
	std::ofstream myfile;
	myfile.open (output_file.data(), ios::out | ios::app );

	xfibresOptions& opts = xfibresOptions::getInstance();
	int nvox = datam_vec.size(); 
	int nfib = opts.nfibres.value();
	int nparams;
	if (m_include_f0)
      		nparams = nfib*3 + 3; 
    	else
      		nparams = nfib*3 + 2;

	thrust::host_vector<double> params_host;
	params_host.resize(nvox*nparams);

	for(int vox=0;vox<nvox;vox++){
		// initialise with a tensor
  		DTI dti(datam_vec[vox],bvecs_vec[vox],bvals_vec[vox]);
  		dti.linfit();

  		// set starting parameters for nonlinear fitting
  		float _th,_ph;
  		cart2sph(dti.get_v1(),_th,_ph);

  		ColumnVector start(nparams);
  		//Initialize the non-linear fitter. Use the DTI estimates for most parameters, apart from the volume fractions
  		start(1) = dti.get_s0();
  		//start(2) = d2lambda(dti.get_md()>0?dti.get_md()*2:0.001); // empirically found that d~2*MD
  		start(2) = d2lambda(dti.get_l1()>0?dti.get_l1():0.002); // empirically found that d~L1
  		start(4) = _th;
  		start(5) = _ph;
  		for(int ii=2,i=6;ii<=nfib;ii++,i+=3){
    			cart2sph(dti.get_v(ii),_th,_ph);
    			start(i+1) = _th;
    			start(i+2) = _ph;
  		}
  
  		// do a better job for initializing the volume fractions
		PVM_single_c pvm(datam_vec[vox],bvecs_vec[vox],bvals_vec[vox],opts.nfibres.value(),false,m_include_f0,false);
  		pvm.fit_pvf(start);

		for(int i=0;i<nparams;i++){ 
			params_host[vox*nparams+i]=start(i+1);
		}
	}

	thrust::copy(params_host.begin(), params_host.end(), params_gpu.begin());	

	int blocks = nvox;
   	dim3 Dim_Grid(blocks,1);
  	dim3 Dim_Block(THREADS_BLOCK_FIT,1);

	int amount_shared = (THREADS_BLOCK_FIT+5*nparams+2*nparams*nparams+4*nfib+nfib*nfib+8)*sizeof(double)+(2+nparams)*sizeof(int);

	myfile << "Shared Memory Used in fit_PVM_single_c: " << amount_shared << "\n"; 

	fit_PVM_single_c_kernel<<<Dim_Grid, Dim_Block, amount_shared>>>(thrust::raw_pointer_cast(datam_gpu.data()), thrust::raw_pointer_cast(bvecs_gpu.data()), thrust::raw_pointer_cast(bvals_gpu.data()) ,nvox, ndirections, nfib, nparams, false, m_include_f0, false, thrust::raw_pointer_cast(params_gpu.data()));
	sync_check("fit_PVM_single_c_kernel");

	myfile.close();
}

void fit_PVM_multi(	//INPUT
			thrust::device_vector<double> 	datam_gpu, 
			thrust::device_vector<double>	bvecs_gpu, 
			thrust::device_vector<double>	bvals_gpu,	
			int 				nvox,	
			int				ndirections,		
			bool 				m_include_f0,
			string 				output_file,
			//OUTPUT
			thrust::device_vector<double>&	params_gpu)
{
	std::ofstream myfile;
	myfile.open (output_file.data(), ios::out | ios::app );

	xfibresOptions& opts = xfibresOptions::getInstance();
	int nfib = opts.nfibres.value();

	int blocks = nvox;
   	dim3 Dim_Grid(blocks,1);
  	dim3 Dim_Block(THREADS_BLOCK_FIT,1);

	int nparams;
	if (m_include_f0)
      		nparams = nfib*3 + 4; 
    	else
      		nparams = nfib*3 + 3;

	thrust::device_vector<double> params_PVM_single_c_gpu; 	//copy params to an auxiliar structure because there are different number of nparams
	params_PVM_single_c_gpu.resize(nvox*nparams);		//between single_c and multi. We must read and write in different structures, 
	thrust::copy(params_gpu.begin(), params_gpu.end(), params_PVM_single_c_gpu.begin());	
								//maybe 1 block finish before other one read their params.

	int amount_shared = (THREADS_BLOCK_FIT+5*nparams+2*nparams*nparams+4*nfib+9)*sizeof(double)+(2+nparams)*sizeof(int);

	myfile << "Shared Memory Used in fit_PVM_multi: " << amount_shared << "\n"; 

	fit_PVM_multi_kernel<<<Dim_Grid, Dim_Block, amount_shared>>>(thrust::raw_pointer_cast(datam_gpu.data()), thrust::raw_pointer_cast(params_PVM_single_c_gpu.data()), thrust::raw_pointer_cast(bvecs_gpu.data()), thrust::raw_pointer_cast(bvals_gpu.data()) ,nvox, ndirections, nfib, nparams, m_include_f0, thrust::raw_pointer_cast(params_gpu.data()));
	sync_check("fit_PVM_multi_kernel");

	myfile.close();
}

void calculate_tau(	//INPUT
			thrust::device_vector<double> 	datam_gpu, 
			thrust::device_vector<double>	params_gpu,
			thrust::device_vector<double>	bvecs_gpu, 
			thrust::device_vector<double>	bvals_gpu,
			thrust::host_vector<int>	vox_repeat,
			int				nrepeat,
			int				ndirections,	
			string 				output_file,	
			//OUTPUT
			thrust::host_vector<float>&	tau)
{
	std::ofstream myfile;
	myfile.open (output_file.data(), ios::out | ios::app );
   	myfile << "--------- CALCULATE TAU/RESIDULAS IN GPU ------------ " << "\n"; 

	struct timeval t1,t2;
   	double time;
   	gettimeofday(&t1,NULL);

	xfibresOptions& opts = xfibresOptions::getInstance();
	int nvox = vox_repeat.size(); 
	int nfib = opts.nfibres.value();

	int nparams;
	if (opts.f0.value())
      		nparams = nfib*3 + 3; 
    	else
      		nparams = nfib*3 + 2;
	if(opts.modelnum.value()==2) nparams++;
	
	thrust::device_vector<bool> includes_f0_gpu;
	includes_f0_gpu.resize(nvox);
	thrust::fill(includes_f0_gpu.begin(), includes_f0_gpu.end(), opts.f0.value());

	if(opts.f0.value()){
		for(int i=0;i<nrepeat;i++){
			includes_f0_gpu[vox_repeat[i]]=false;	//if has been reprocessed f0 will be 0.
		}
	}

	int blocks = nvox;
   	dim3 Dim_Grid(blocks,1);
  	dim3 Dim_Block(THREADS_BLOCK_FIT,1);

	int amount_shared = (nparams+4*nfib+3)*sizeof(double) + sizeof(int);

	thrust::device_vector<double> residuals_gpu;
	residuals_gpu.resize(nvox*ndirections);

	if(opts.modelnum.value()==1){
		if(opts.nonlin.value()){ 
			get_residuals_PVM_single_kernel<<<Dim_Grid, Dim_Block,amount_shared>>>(thrust::raw_pointer_cast(datam_gpu.data()), thrust::raw_pointer_cast(params_gpu.data()), thrust::raw_pointer_cast(bvecs_gpu.data()), thrust::raw_pointer_cast(bvals_gpu.data()), nvox, ndirections, nfib, nparams, opts.f0.value(), thrust::raw_pointer_cast(includes_f0_gpu.data()), thrust::raw_pointer_cast(residuals_gpu.data()));
			sync_check("get_residuals_PVM_single_kernel");

		}else{
			get_residuals_PVM_single_c_kernel<<<Dim_Grid, Dim_Block,amount_shared>>>(thrust::raw_pointer_cast(datam_gpu.data()), thrust::raw_pointer_cast(params_gpu.data()), thrust::raw_pointer_cast(bvecs_gpu.data()), thrust::raw_pointer_cast(bvals_gpu.data()), nvox, ndirections, nfib, nparams, opts.f0.value(), thrust::raw_pointer_cast(includes_f0_gpu.data()), thrust::raw_pointer_cast(residuals_gpu.data()));
			sync_check("get_residuals_PVM_single_c_kernel");
		}
	}else{
      		//model 2 : non-mono-exponential
		get_residuals_PVM_multi_kernel<<<Dim_Grid, Dim_Block,amount_shared>>>(thrust::raw_pointer_cast(datam_gpu.data()), thrust::raw_pointer_cast(params_gpu.data()), thrust::raw_pointer_cast(bvecs_gpu.data()), thrust::raw_pointer_cast(bvals_gpu.data()), nvox, ndirections, nfib, nparams, opts.f0.value(), thrust::raw_pointer_cast(includes_f0_gpu.data()), thrust::raw_pointer_cast(residuals_gpu.data()));
		sync_check("get_residuals_PVM_multi_kernel");
	}

	thrust::host_vector<double> residuals_host;
	residuals_host.resize(nvox*ndirections);
	thrust::copy(residuals_gpu.begin(), residuals_gpu.end(), residuals_host.begin());	

	ColumnVector res(ndirections);
	for(int vox=0;vox<nvox;vox++){	
		for(int i=0;i<ndirections;i++) res(i+1)= residuals_host[vox*ndirections+i];  

		float variance=var(res).AsScalar();
		tau[vox]=1.0/variance;
	}

	gettimeofday(&t2,NULL);
    	time=timeval_diff(&t2,&t1);
   	myfile << "TIME TOTAL: " << time << " seconds\n"; 
	myfile << "-----------------------------------------------------" << "\n\n" ; 	
	myfile.close();			
}

