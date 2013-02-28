/*  PVM_single_c.cu

    Tim Behrens, Saad Jbabdi, Stam Sotiropoulos, Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */

#include "diffmodels_utils.h"
#include "levenberg_marquardt.cu"
#include "options.h"

#include <fstream>

/////////////////////////////////////
/////////////////////////////////////
/// 	    PVM_single_c	  /// 
/////////////////////////////////////
/////////////////////////////////////

__device__ 
inline double isoterm_PVM_single_c(const int pt,const double _d,const double *bvals){
  	return exp(-bvals[pt]*_d);
}

__device__ 
inline double isoterm_lambda_PVM_single_c(const int pt,const double lambda,const double *bvals){
  	return(-2*bvals[pt]*lambda*exp(-bvals[pt]*lambda*lambda));
}

__device__ 
inline double anisoterm_PVM_single_c(const int pt,const double _d,const double3 x, const double *bvecs, const double *bvals){
	double dp = bvecs[pt]*x.x+bvecs[NDIRECTIONS+pt]*x.y+bvecs[(2*NDIRECTIONS)+pt]*x.z;
	return exp(-bvals[pt]*_d*dp*dp);
}

__device__ 
inline double anisoterm_lambda_PVM_single_c(const int pt,const double lambda,const double3 x, const double *bvecs, const double *bvals){
	double dp = bvecs[pt]*x.x+bvecs[NDIRECTIONS+pt]*x.y+bvecs[(2*NDIRECTIONS)+pt]*x.z;
	return(-2*bvals[pt]*lambda*dp*dp*exp(-bvals[pt]*lambda*lambda*dp*dp));
}

__device__ 
inline double anisoterm_th_PVM_single_c(const int pt,const double _d,const double3 x, const double _th,const double _ph,const double *bvecs, const double *bvals){
	double sinth,costh,sinph,cosph;
	sincos(_th,&sinth,&costh);
	sincos(_ph,&sinph,&cosph);
	double dp = bvecs[pt]*x.x+bvecs[NDIRECTIONS+pt]*x.y+bvecs[(2*NDIRECTIONS)+pt]*x.z;
	double dp1 = costh*(bvecs[pt]*cosph+bvecs[NDIRECTIONS+pt]*sinph)-bvecs[(2*NDIRECTIONS)+pt]*sinth;
  	return(-2*bvals[pt]*_d*dp*dp1*exp(-bvals[pt]*_d*dp*dp));
}

__device__ 
inline double anisoterm_ph_PVM_single_c(const int pt,const double _d,const double3 x, const double _th,const double _ph,const double *bvecs, const double *bvals){
	double sinth,sinph,cosph;
	sinth=sin(_th);
	sincos(_ph,&sinph,&cosph);
  	double dp = bvecs[pt]*x.x+bvecs[NDIRECTIONS+pt]*x.y+bvecs[(2*NDIRECTIONS)+pt]*x.z;
	double dp1 = sinth*(-bvecs[pt]*sinph+bvecs[NDIRECTIONS+pt]*cosph);
  	return(-2*bvals[pt]*_d*dp*dp1*exp(-bvals[pt]*_d*dp*dp));
}

//If the sum of the fractions is >1, then zero as many fractions
//as necessary, so that the sum becomes smaller than 1.
//in diffmodel.cc
__device__ void fix_fsum_PVM_single_c(		//INPUT 
						int nfib,
						//INPUT - OUTPUT){
						double *fs)
{
  	double sumf=0.0;
  	for(int i=0;i<nfib;i++){
    		sumf+=fs[i];
    		if(sumf>=1){
      			for(int j=i;j<nfib;j++) 
				fs[j]=FSMALL_gpu;  //make the fraction almost zero
      			break;
    		}
  	}
}

//Returns 1-Sum(f_j), 1<=j<=ii. (ii<=nfib)
//Used for transforming beta to f and vice versa
//in diffmodel.cc
__device__ double partial_fsum_PVM_single_c(double* fs, int ii){
  	double fsum=1.0;
  	for(int j=0;j<ii;j++){
   		fsum-=fs[j];
	}
  	return fsum;
}

//in diffmodel.cc
__device__ void sort_PVM_single_c(int nfib,double* params)
{
	double temp_f, temp_th, temp_ph;
	// Order vector descending using f parameters as index
  	for(int i=1; i<(nfib); i++){ 
    		for(int j=0; j<(nfib-i); j++){ 
      			if (params[2+j*3] < params[2+(j+1)*3]){ 
        			temp_f = params[2+j*3];
				temp_th = params[2+j*3+1];
				temp_ph = params[2+j*3+2];
        			params[2+j*3] = params[2+(j+1)*3]; 
				params[2+j*3+1] = params[2+(j+1)*3+1]; 
				params[2+j*3+2] = params[2+(j+1)*3+2]; 
        			params[2+(j+1)*3] = temp_f; 
				params[2+(j+1)*3+1] = temp_th; 
				params[2+(j+1)*3+2] = temp_ph; 
      			} 
    		} 
  	} 
}

__device__  void fractions_deriv_PVM_single_c(	//INPUT
						const double*	params,
						const double* 	fs, 
						const int	nfib,
						const int	idB,
						//OUTPUT
						double* 	Deriv) 
{
	int nparams_per_fibre=3;
  	double fsum;
	int k=idB%nfib;
	for (int j=0; j<nfib; j++){
		Deriv[j*nfib+k]=0;
    	}

  	int kk = 2+(k*nparams_per_fibre);
	double sinparamkk = sin(2*params[kk]);

	for (int j=0; j<nfib; j++){
		int jj = 2+(j*nparams_per_fibre);
      		if (j==k){
			fsum=1; 
			for (int n=0; n<=(j-1); n++){
	  			fsum-=fs[n];
			}
			Deriv[j*nfib+k]=sinparamkk*fsum;
      		}else if (j>k){
			double sinparam = sin(params[jj]);
			fsum=0;
			for (int n=0; n<=(j-1); n++){
	  			fsum+=Deriv[n*nfib+k];
			}
			Deriv[j*nfib+k]=  -(sinparam*sinparam)*fsum;
      		}
    	}
}

//cost function PVM_single_c
__device__ void cf_PVM_single_c(	//INPUT
					const double*		params,
					const double*		mdata,
					const double*		bvecs, 
					const double*		bvals,
					const int 		nparams,
					const bool 		m_include_f0,
					const int		idB,
					double*			shared,		//shared memory
					double* 		fs,		//shared memory
					double*			x,		//shared memory	
					double 			&_d,		//shared memory
					double 			&sumf,		//shared memory
					//OUTPUT
					double			&cfv)
{
	if(idB<NFIBRES){
		int kk = 2+3*(idB);
		double sinth,costh,sinph,cosph;
		sincos(params[kk+1],&sinth,&costh);
		sincos(params[kk+2],&sinph,&cosph);
		x[idB*3] = sinth*cosph;
    		x[idB*3+1] = sinth*sinph;
    		x[idB*3+2] = costh;
  	}
	if(idB==0){
		_d = lambda2d_gpu(params[1]);
		cfv = 0.0;
		sumf=0;
		double partial_fsum;
		for(int k=0;k<NFIBRES;k++){
			int kk = 2+3*(k);
    			//partial_fsum ///////////
			partial_fsum=1.0;
			for(int j=0;j<k;j++)
				partial_fsum-=fs[j];
    			//////////////////////////
			fs[k] = beta2f_gpu(params[kk])*partial_fsum;
    			sumf += fs[k];
		}
	}

	int ndir = NDIRECTIONS/THREADS_X_BLOCK_FIT;
	if(idB<(NDIRECTIONS%THREADS_X_BLOCK_FIT)) ndir++;
	
	double err;
	double3 x2;
	int dir_iter=idB;

	__syncthreads();
	
	shared[idB]=0;
	for(int dir=0;dir<ndir;dir++){
		err = 0.0;
    		for(int k=0;k<NFIBRES;k++){
			x2.x=x[k*3];
			x2.y=x[k*3+1];
			x2.z=x[k*3+2];	
			err += fs[k]*anisoterm_PVM_single_c(dir_iter,_d,x2,bvecs,bvals); 
    		}
		if(m_include_f0){
			//partial_fsum ///////////
			double partial_fsum=1.0;
			for(int j=0;j<NFIBRES;j++)
				partial_fsum-=fs[j];
	     		//////////////////////////
			double temp_f0=beta2f_gpu(params[nparams-1])*partial_fsum;
			err= (params[0]*((temp_f0+(1-sumf-temp_f0)*isoterm_PVM_single_c(dir_iter,_d,bvals))+err))-mdata[dir_iter];
		}else{
			err = params[0]*((1-sumf)*isoterm_PVM_single_c(dir_iter,_d,bvals)+err)-mdata[dir_iter];
		}
		shared[idB]+= err*err;  
		dir_iter+=THREADS_X_BLOCK_FIT;
  	}  

	__syncthreads();

	if(idB==0){
		for(int i=0;i<THREADS_X_BLOCK_FIT;i++){
			cfv+=shared[i];
		}	
	}	
}


//gradient function PVM_single_c
__device__ void grad_PVM_single_c(	//INPUT
					const double*		params,
					const double*		mdata,
					const double*		bvecs, 
					const double*		bvals,
					const int 		nparams,
					const bool 		m_include_f0,
					const int		idB,		
					double*			shared,		//shared memory
					double* 		fs,		//shared memory
					double*			f_deriv,	//shared memory
					double*			x,		//shared memory
					double 			&_d,		//shared memory
					double 			&sumf,		//shared memory
					//OUTPUT
					double*			grad)
{
	if(idB<NFIBRES){
		int kk = 2+3*(idB);
		double sinth,costh,sinph,cosph;
		sincos(params[kk+1],&sinth,&costh);
		sincos(params[kk+2],&sinph,&cosph);
    		x[idB*3] = sinth*cosph;
    		x[idB*3+1] = sinth*sinph;
    		x[idB*3+2] = costh;
  	}
	if(idB==0){
		_d = lambda2d_gpu(params[1]);
		sumf=0;
		double partial_fsum;
		for(int k=0;k<NFIBRES;k++){
			int kk = 2+3*(k);
    			//partial_fsum ///////////
			partial_fsum=1.0;
			for(int j=0;j<k;j++)
				partial_fsum-=fs[j];
    			//////////////////////////
			fs[k] = beta2f_gpu(params[kk])*partial_fsum;
    			sumf += fs[k];
		}
		for (int p=0;p<nparams;p++) grad[p]=0;
	}

	__syncthreads();

  	if(idB<NFIBRES){ 
		fractions_deriv_PVM_single_c(params,fs,NFIBRES,idB,f_deriv); 
	} 

	int ndir = NDIRECTIONS/THREADS_X_BLOCK_FIT;
	if(idB<(NDIRECTIONS%THREADS_X_BLOCK_FIT)) ndir++;
	int max_dir = NDIRECTIONS/THREADS_X_BLOCK_FIT;
	if(NDIRECTIONS%THREADS_X_BLOCK_FIT) max_dir++;

	double J[NPARAMS];
	double diff;
  	double sig;
	double Iso_term;
	double3 xx;
	int dir_iter=idB;
  	double Aniso_terms[NFIBRES];

	__syncthreads();

  	for(int dir=0;dir<max_dir;dir++){
		for (int p=0; p<nparams; p++) J[p]=0;
		if(dir<ndir){
    			Iso_term=isoterm_PVM_single_c(dir_iter,_d,bvals);  //Precompute some terms for this datapoint
    			for(int k=0;k<NFIBRES;k++){
      				xx.x=x[k*3];
      				xx.y=x[k*3+1];
     				xx.z=x[k*3+2];	
      				Aniso_terms[k]=anisoterm_PVM_single_c(dir_iter,_d,xx,bvecs,bvals);
    			}
    			sig = 0;
    			for(int k=0;k<NFIBRES;k++){
     				int kk = 2+3*(k);
      				xx.x=x[k*3];
      				xx.y=x[k*3+1];
      				xx.z=x[k*3+2];		
      				sig += fs[k]*Aniso_terms[k];
     				J[1] += params[0]*fs[k]*anisoterm_lambda_PVM_single_c(dir_iter,params[1],xx,bvecs,bvals);
     				J[kk] = 0;
      				for (int j=0;j<NFIBRES;j++){
					if(f_deriv[j*NFIBRES+k]!=0){
	  					J[kk] += params[0]*(Aniso_terms[j]-Iso_term)*f_deriv[j*NFIBRES+k]; 
					}
      				}
      				J[kk+1] = params[0]*fs[k]*anisoterm_th_PVM_single_c(dir_iter,_d,xx,params[kk+1],params[kk+2],bvecs,bvals);
      				J[kk+2] = params[0]*fs[k]*anisoterm_ph_PVM_single_c(dir_iter,_d,xx,params[kk+1],params[kk+2],bvecs,bvals);
    			}
    			if(m_include_f0){
				//partial_fsum ///////////
    				double partial_fsum=1.0;
    				for(int j=0;j<(NFIBRES);j++)
					partial_fsum-=fs[j];
				//////////////////////////
				double temp_f0=beta2f_gpu(params[nparams-1])*partial_fsum;

    				//derivative with respect to f0
    				J[nparams-1]= params[0]*(1-Iso_term)*sin(double(2*params[nparams-1]))*partial_fsum; 
				sig=params[0]*((temp_f0+(1-sumf-temp_f0)*Iso_term)+sig);
    				J[1] += params[0]*(1-sumf-temp_f0)*isoterm_lambda_PVM_single_c(dir_iter,params[1],bvals);
    			}else{
				sig = params[0]*((1-sumf)*Iso_term+sig);
	    			J[1] += params[0]*(1-sumf)*isoterm_lambda_PVM_single_c(dir_iter,params[1],bvals);
    			}
    			diff = sig - mdata[dir_iter];
    			J[0] = sig/params[0]; 
		}

		for (int p=0;p<nparams;p++){ 
			shared[idB]=2*J[p]*diff;

			__syncthreads();
			if(idB==0){
				for(int i=0;i<THREADS_X_BLOCK_FIT;i++){
					grad[p] += shared[i];
				}
			}
			__syncthreads(); 
		} 
		dir_iter+=THREADS_X_BLOCK_FIT;
  	}
}


//hessian function PVM_single_c
__device__ void hess_PVM_single_c(	//INPUT
					const double*		params,
					const double*		bvecs, 
					const double*		bvals,
					const int 		nparams,
					const bool 		m_include_f0,
					const int		idB,		
					double*			shared,		//shared memory
					double* 		fs,		//shared memory
					double*			f_deriv,	//shared memory
					double*			x,		//shared memory
					double 			&_d,		//shared memory
					double 			&sumf,		//shared memory
					//OUTPUT
					double*			hess)
{
	if(idB<NFIBRES){
		int kk = 2+3*(idB);
		double sinth,costh,sinph,cosph;
		sincos(params[kk+1],&sinth,&costh);
		sincos(params[kk+2],&sinph,&cosph);
    		x[idB*3] = sinth*cosph;
    		x[idB*3+1] = sinth*sinph;
    		x[idB*3+2] = costh;
  	}
	if(idB==0){
		_d = lambda2d_gpu(params[1]);
		sumf=0;
		double partial_fsum;
		for(int k=0;k<NFIBRES;k++){
			int kk = 2+3*(k);
    			//partial_fsum ///////////
			partial_fsum=1.0;
			for(int j=0;j<k;j++)
				partial_fsum-=fs[j];
    			//////////////////////////
			fs[k] = beta2f_gpu(params[kk])*partial_fsum;
    			sumf += fs[k];
		}
		for (int p=0;p<nparams;p++){
			for (int p2=0;p2<nparams;p2++){ 
				hess[p*nparams+p2] = 0;
			}
		}
	}

	__syncthreads();

  	if(idB<NFIBRES){ 
		fractions_deriv_PVM_single_c(params,fs,NFIBRES,idB,f_deriv); 
	} 

  	int ndir = NDIRECTIONS/THREADS_X_BLOCK_FIT;
	if(idB<(NDIRECTIONS%THREADS_X_BLOCK_FIT)) ndir++;
	int max_dir = NDIRECTIONS/THREADS_X_BLOCK_FIT;
	if(NDIRECTIONS%THREADS_X_BLOCK_FIT) max_dir++;

	double J[NPARAMS];
  	double sig;
	double Iso_term;
	double3 xx;
	int dir_iter=idB;
  	double Aniso_terms[NFIBRES];

	__syncthreads();

  	for(int dir=0;dir<max_dir;dir++){
		for (int p=0; p<nparams; p++) J[p]=0;
		if(dir<ndir){
    			Iso_term=isoterm_PVM_single_c(dir_iter,_d,bvals);  //Precompute some terms for this datapoint
    			for(int k=0;k<NFIBRES;k++){
      				xx.x=x[k*3];
      				xx.y=x[k*3+1];
      				xx.z=x[k*3+2];	
      				Aniso_terms[k]=anisoterm_PVM_single_c(dir_iter,_d,xx,bvecs,bvals);
    			}
    			sig = 0;
    			for(int k=0;k<NFIBRES;k++){
      				int kk = 2+3*(k);
      				xx.x=x[k*3];
      				xx.y=x[k*3+1];
      				xx.z=x[k*3+2];		 
      				sig += fs[k]*Aniso_terms[k];
      				J[1] += params[0]*fs[k]*anisoterm_lambda_PVM_single_c(dir_iter,params[1],xx,bvecs,bvals);	 
      				J[kk] = 0;
      				for (int j=0; j<NFIBRES; j++){
					if (f_deriv[j*NFIBRES+k]!=0)
	  				J[kk] += params[0]*(Aniso_terms[j]-Iso_term)*f_deriv[j*NFIBRES+k]; 
      				}
      				J[kk+1] = params[0]*fs[k]*anisoterm_th_PVM_single_c(dir_iter,_d,xx,params[kk+1],params[kk+2],bvecs,bvals);
      				J[kk+2] = params[0]*fs[k]*anisoterm_ph_PVM_single_c(dir_iter,_d,xx,params[kk+1],params[kk+2],bvecs,bvals);
    			}
    			if(m_include_f0){
				//partial_fsum ///////////
	    			double partial_fsum=1.0;
	    			for(int j=0;j<(NFIBRES);j++)
					partial_fsum-=fs[j];
	    			//////////////////////////
    				double temp_f0=beta2f_gpu(params[nparams-1])*partial_fsum;
    				//derivative with respect to f0
    				J[nparams-1]= params[0]*(1-Iso_term)*sin(double(2*params[nparams-1]))*partial_fsum; 
				sig= params[0]*((temp_f0+(1-sumf-temp_f0)*Iso_term)+sig);
    				J[1] += params[0]*(1-sumf-temp_f0)*isoterm_lambda_PVM_single_c(dir_iter,params[1],bvals);
    			}else{
	    			sig = params[0]*((1-sumf)*Iso_term+sig);
	    			J[1] += params[0]*(1-sumf)*isoterm_lambda_PVM_single_c(dir_iter,params[1],bvals);
    			}
    			J[0] = sig/params[0]; 
		}

		for (int p=0;p<nparams;p++){
			for (int p2=p;p2<nparams;p2++){ 

				shared[idB]=2*(J[p]*J[p2]);
				__syncthreads();
				if(idB==0){
					for(int i=0;i<THREADS_X_BLOCK_FIT;i++){
						hess[p*nparams+p2] += shared[i];
					}
				}
				__syncthreads(); 
			}
		}
		dir_iter+=THREADS_X_BLOCK_FIT;
  	}

	if(idB==0){
  		for (int j=0; j<nparams; j++) {
    			for (int i=j+1; i<nparams; i++) {
     				hess[i*nparams+j]=hess[j*nparams+i];	
    			}
  		}
	}
}
//in diffmodel.cc
extern "C" __global__ void fit_PVM_single_c_kernel(	//INPUT
							const double* 		data, 
							const double* 		bvecs, 
							const double* 		bvals, 
							const int 		nvox, 
							const int 		nfib, 
							const bool		m_eval_BIC,
							const bool 		m_include_f0,
							const bool	 	m_return_fanning,
							//INPUT - OUTPUT
							double* 		params)
{
	int idB = threadIdx.x;
	int idVOX = blockIdx.x;

	__shared__ double myparams[NPARAMS];
	__shared__ int nparams;

	__shared__ double shared[THREADS_X_BLOCK_FIT]; 
	__shared__ double step[NPARAMS];
	__shared__ double grad[NPARAMS];                          
   	__shared__ double hess[NPARAMS*NPARAMS]; 	
	__shared__ double inverse[NPARAMS];
	__shared__ double pcf;
	__shared__ double ncf;
	__shared__ double lambda;
	__shared__ double cftol;
	__shared__ double ltol;
	__shared__ double olambda;
	__shared__ bool success;    
	__shared__ bool end;    

	__shared__ double fs[NFIBRES];
	__shared__ double f_deriv[NFIBRES*NFIBRES];	
  	__shared__ double x[NFIBRES*3];	
	__shared__ double _d;
  	__shared__ double sumf;

	if(idB==0){
		if (m_include_f0)
      			nparams = nfib*3 + 3; 
    		else
      			nparams = nfib*3 + 2;

		for(int i=0;i<nparams;i++){
			myparams[i]=params[(idVOX*nparams)+i];
   		}
	}

	__syncthreads();

	//do the fit
	levenberg_marquardt_PVM_single_c_gpu(&data[idVOX*NDIRECTIONS],&bvecs[idVOX*3*NDIRECTIONS],&bvals[idVOX*NDIRECTIONS],nparams,m_include_f0,idB,step,grad,hess,inverse, pcf,ncf,lambda,cftol,ltol,olambda,success,end,shared,fs,f_deriv,x,_d,sumf,myparams);

	__syncthreads();

	// finalise parameters
	// m_s0-myparams[0] 	m_d-myparams[1] 	m_f-m_th-m_ph-myparams[2,3,4,5, etc..]   	m_f0-myparams[nparams-1]
	
	if(idB==0){
		double m_f[NFIBRES]; 					// for partial_fsum

  		myparams[1] = lambda2d_gpu(myparams[1]); 
  		for(int k=0;k<nfib;k++){
    			int kk = 2 + 3*(k);
    			myparams[kk]  = beta2f_gpu(myparams[kk])*partial_fsum_PVM_single_c(m_f,k);
			m_f[k]=myparams[kk];
  		}
  
  		if (m_include_f0)
    			myparams[nparams-1]= beta2f_gpu(myparams[nparams-1])*partial_fsum_PVM_single_c(m_f,nfib);

		sort_PVM_single_c(nfib,myparams);

		for(int i=0;i<nparams;i++){
			params[(idVOX*nparams)+i] = myparams[i];
   		}
	}
}

//in diffmodel.cc
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
								double*			residuals)
{
	int idB = threadIdx.x;
	int idVOX = blockIdx.x;

	__shared__ double myparams[NPARAMS];
	__shared__ int nparams;
	__shared__ bool my_include_f0;
	__shared__ double val;
  	__shared__ double _d;
  	__shared__ double fs[NFIBRES];
  	__shared__ double x[NFIBRES*3];	
  	__shared__ double sumf;

	double predicted_signal;
	double mydata;

	if(idB==0){
		if (m_include_f0)
      			nparams = nfib*3 + 3; 
    		else
      			nparams = nfib*3 + 2;

		my_include_f0 = includes_f0[idVOX];

		//m_s0-myparams[0]  m_d-myparams[1]  m_d_std-myparams[2]  m_f-m_th-m_ph-myparams[3,4,5,6 etc..]  m_f0-myparams[nparams-1]
  		
  		myparams[0]=params[(idVOX*nparams)+0];
		if(myparams[1]<0)  myparams[1] = 0;	//This can be due to numerical errors..sqrt
  		else myparams[1] = d2lambda_gpu(params[(idVOX*nparams)+1]);

		double partial_fsum;	
  		for(int k=0;k<nfib;k++){
    			int kk = 2+3*k;
			//partial_fsum ///////////
			partial_fsum=1.0;
			for(int j=0;j<k;j++)
				partial_fsum-=fs[j];
	     		//////////////////////////
			fs[k] = params[(idVOX*nparams)+kk];
			double tmpr=fs[k]/partial_fsum;
    			if (tmpr>1.0) tmpr=1; //This can be due to numerical errors
			if (tmpr<0.0) tmpr=0; //This can be due to numerical errors..sqrt
    			myparams[kk]   = f2beta_gpu(tmpr);
    			myparams[kk+1] = params[(idVOX*nparams)+kk+1];
    			myparams[kk+2] = params[(idVOX*nparams)+kk+2];
  		}
  		if (my_include_f0){
			//partial_fsum ///////////
			partial_fsum=1.0;
			for(int j=0;j<NFIBRES;j++)
				partial_fsum-=fs[j];
	     		//////////////////////////	
			double tmpr=params[(idVOX*nparams)+nparams-1]/partial_fsum;
    			if (tmpr>1.0) tmpr=1; //This can be due to numerical errors..asin
			if (tmpr<0.0) tmpr=0; //This can be due to numerical errors..sqrt
    			myparams[nparams-1]= f2beta_gpu(tmpr);	
		}
	}

	__syncthreads();

	if(idB<nfib){
		int kk = 2+3*idB;
		double sinth,costh,sinph,cosph;
		sincos(myparams[kk+1],&sinth,&costh);
		sincos(myparams[kk+2],&sinph,&cosph);
    		x[idB*3] = sinth*cosph;
    		x[idB*3+1] = sinth*sinph;
    		x[idB*3+2] = costh;
  	}

	if(idB==0){
		double partial_fsum;	
		sumf=0;
		for(int k=0;k<nfib;k++){
    			int kk = 2+3*k;
			////// partial_fsum //////
			partial_fsum=1.0;
			for(int j=0;j<k;j++)
				partial_fsum-=fs[j];
    			//////////////////////////
	    		fs[k] = beta2f_gpu(myparams[kk])*partial_fsum;
	    		sumf += fs[k];
		}
		_d = lambda2d_gpu(myparams[1]);
	}

	int ndir = NDIRECTIONS/THREADS_X_BLOCK_FIT;
	if(idB<(NDIRECTIONS%THREADS_X_BLOCK_FIT)) ndir++;
	
	double3 x2;
	int dir_iter=idB; 

	__syncthreads();

	for(int dir=0;dir<ndir;dir++){
		mydata = data[(idVOX*NDIRECTIONS)+dir_iter];
		predicted_signal=0;	//pred = 0;
		val = 0.0;
    		for(int k=0;k<nfib;k++){
			x2.x=x[k*3];
			x2.y=x[k*3+1];
			x2.z=x[k*3+2];	 
      			val += fs[k]*anisoterm_PVM_single_c(dir_iter,_d,x2,&bvecs[idVOX*3*NDIRECTIONS],&bvals[idVOX*NDIRECTIONS]);
    		}	
    		if (my_include_f0){
			//partial_fsum ///////////
			double partial_fsum=1.0;
			for(int j=0;j<NFIBRES;j++)
				partial_fsum-=fs[j];
	     		//////////////////////////
      			double temp_f0= beta2f_gpu(myparams[nparams-1])*partial_fsum;
      			predicted_signal = myparams[0]*(temp_f0+(1-sumf-temp_f0)*isoterm_PVM_single_c(dir_iter,_d,&bvals[idVOX*NDIRECTIONS])+val);
    		} 
    		else
      			predicted_signal = myparams[0]*((1-sumf)*isoterm_PVM_single_c(dir_iter,_d,&bvals[idVOX*NDIRECTIONS])+val); 

		//residuals=m_data-predicted_signal;
		residuals[idVOX*NDIRECTIONS+dir_iter]= mydata - predicted_signal;

		dir_iter+=THREADS_X_BLOCK_FIT;
	}
}
