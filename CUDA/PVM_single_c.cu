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
inline double isoterm_PVM_single_c(const int pt,const double* _d,const double *bvals){
  	return exp(-bvals[pt]**_d);
}

__device__ 
inline double isoterm_lambda_PVM_single_c(const int pt,const double lambda,const double *bvals){
  	return(-2*bvals[pt]*lambda*exp(-bvals[pt]*lambda*lambda));
}

__device__ 
inline double anisoterm_PVM_single_c(const int pt,const double* _d,const double3 x, const double *bvecs, const double *bvals, const int ndirections){
	double dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
	return exp(-bvals[pt]**_d*dp*dp);
}

__device__ 
inline double anisoterm_lambda_PVM_single_c(const int pt,const double lambda,const double3 x, const double *bvecs, const double *bvals, const int ndirections){
	double dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
	return(-2*bvals[pt]*lambda*dp*dp*exp(-bvals[pt]*lambda*lambda*dp*dp));
}

__device__ 
inline double anisoterm_th_PVM_single_c(const int pt,const double* _d,const double3 x, const double _th,const double _ph,const double *bvecs, const double *bvals, const int ndirections){
	double sinth,costh,sinph,cosph;
	sincos(_th,&sinth,&costh);
	sincos(_ph,&sinph,&cosph);
	double dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
	double dp1 = costh*(bvecs[pt]*cosph+bvecs[ndirections+pt]*sinph)-bvecs[(2*ndirections)+pt]*sinth;
  	return(-2*bvals[pt]**_d*dp*dp1*exp(-bvals[pt]**_d*dp*dp));
}

__device__ 
inline double anisoterm_ph_PVM_single_c(const int pt,const double* _d,const double3 x, const double _th,const double _ph,const double *bvecs, const double *bvals, const int ndirections){
	double sinth,sinph,cosph;
	sinth=sin(_th);
	sincos(_ph,&sinph,&cosph);
  	double dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
	double dp1 = sinth*(-bvecs[pt]*sinph+bvecs[ndirections+pt]*cosph);
  	return(-2*bvals[pt]**_d*dp*dp1*exp(-bvals[pt]**_d*dp*dp));
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
						const int	idSubVOX,
						//OUTPUT
						double* 	Deriv) 
{
	int nparams_per_fibre=3;
  	double fsum;
	int k=idSubVOX%nfib;
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
					const int 		ndirections,
					const int		nfib,
					const int 		nparams,
					const bool 		m_include_f0,
					const int		idSubVOX,
					double*			reduction,	//shared memory
					double* 		fs,		//shared memory
					double*			x,		//shared memory	
					double* 		_d,		//shared memory
					double* 		sumf,		//shared memory
					//OUTPUT
					double*			cfv)
{
	if(idSubVOX<nfib){
		int kk = 2+3*(idSubVOX);
		double sinth,costh,sinph,cosph;
		sincos(params[kk+1],&sinth,&costh);
		sincos(params[kk+2],&sinph,&cosph);
		x[idSubVOX*3] = sinth*cosph;
    		x[idSubVOX*3+1] = sinth*sinph;
    		x[idSubVOX*3+2] = costh;
  	}
	if(idSubVOX==0){
		*_d = lambda2d_gpu(params[1]);
		*cfv = 0.0;
		*sumf=0;
		double partial_fsum;
		for(int k=0;k<nfib;k++){
			int kk = 2+3*(k);
    			//partial_fsum ///////////
			partial_fsum=1.0;
			for(int j=0;j<k;j++)
				partial_fsum-=fs[j];
    			//////////////////////////
			fs[k] = beta2f_gpu(params[kk])*partial_fsum;
    			*sumf += fs[k];
		}
	}

	int ndir = ndirections/THREADS_BLOCK_FIT;
	if(idSubVOX<(ndirections%THREADS_BLOCK_FIT)) ndir++;
	
	double err;
	double3 x2;
	int dir_iter=idSubVOX;

	__syncthreads();
	
	reduction[idSubVOX]=0;
	for(int dir=0;dir<ndir;dir++){
		err = 0.0;
    		for(int k=0;k<nfib;k++){
			x2.x=x[k*3];
			x2.y=x[k*3+1];
			x2.z=x[k*3+2];	
			err += fs[k]*anisoterm_PVM_single_c(dir_iter,_d,x2,bvecs,bvals,ndirections); 
    		}
		if(m_include_f0){
			//partial_fsum ///////////
			double partial_fsum=1.0;
			for(int j=0;j<nfib;j++)
				partial_fsum-=fs[j];
	     		//////////////////////////
			double temp_f0=beta2f_gpu(params[nparams-1])*partial_fsum;
			err= (params[0]*((temp_f0+(1-*sumf-temp_f0)*isoterm_PVM_single_c(dir_iter,_d,bvals))+err))-mdata[dir_iter];
		}else{
			err = params[0]*((1-*sumf)*isoterm_PVM_single_c(dir_iter,_d,bvals)+err)-mdata[dir_iter];
		}
		reduction[idSubVOX]+= err*err;  
		dir_iter+=THREADS_BLOCK_FIT;
  	}  

	__syncthreads();

	if(idSubVOX==0){
		for(int i=0;i<THREADS_BLOCK_FIT;i++){
			*cfv+=reduction[i];
		}	
	}	
}


//gradient function PVM_single_c
__device__ void grad_PVM_single_c(	//INPUT
					const double*		params,
					const double*		mdata,
					const double*		bvecs, 
					const double*		bvals,
					const int		ndirections,
					const int		nfib,
					const int 		nparams,
					const bool 		m_include_f0,
					const int		idSubVOX,		
					double*			reduction,	//shared memory
					double* 		fs,		//shared memory
					double*			f_deriv,	//shared memory
					double*			x,		//shared memory
					double* 		_d,		//shared memory
					double* 		sumf,		//shared memory
					//OUTPUT
					double*			grad)
{
	if(idSubVOX<nfib){
		int kk = 2+3*(idSubVOX);
		double sinth,costh,sinph,cosph;
		sincos(params[kk+1],&sinth,&costh);
		sincos(params[kk+2],&sinph,&cosph);
    		x[idSubVOX*3] = sinth*cosph;
    		x[idSubVOX*3+1] = sinth*sinph;
    		x[idSubVOX*3+2] = costh;
  	}
	if(idSubVOX==0){
		*_d = lambda2d_gpu(params[1]);
		*sumf=0;
		double partial_fsum;
		for(int k=0;k<nfib;k++){
			int kk = 2+3*(k);
    			//partial_fsum ///////////
			partial_fsum=1.0;
			for(int j=0;j<k;j++)
				partial_fsum-=fs[j];
    			//////////////////////////
			fs[k] = beta2f_gpu(params[kk])*partial_fsum;
    			*sumf += fs[k];
		}
		for (int p=0;p<nparams;p++) grad[p]=0;
	}

	__syncthreads();

  	if(idSubVOX<nfib){ 
		fractions_deriv_PVM_single_c(params,fs,nfib,idSubVOX,f_deriv); 
	} 

	int ndir = ndirections/THREADS_BLOCK_FIT;
	if(idSubVOX<(ndirections%THREADS_BLOCK_FIT)) ndir++;
	int max_dir = ndirections/THREADS_BLOCK_FIT;
	if(ndirections%THREADS_BLOCK_FIT) max_dir++;

	double J[MAXNPARAMS];
	double diff;
  	double sig;
	double Iso_term;
	double3 xx;
	int dir_iter=idSubVOX;
  	double Aniso_terms[MAXNFIBRES];

	__syncthreads();

  	for(int dir=0;dir<max_dir;dir++){
		for (int p=0; p<nparams; p++) J[p]=0;
		if(dir<ndir){
    			Iso_term=isoterm_PVM_single_c(dir_iter,_d,bvals);  //Precompute some terms for this datapoint
    			for(int k=0;k<nfib;k++){
      				xx.x=x[k*3];
      				xx.y=x[k*3+1];
     				xx.z=x[k*3+2];	
      				Aniso_terms[k]=anisoterm_PVM_single_c(dir_iter,_d,xx,bvecs,bvals,ndirections);
    			}
    			sig = 0;
    			for(int k=0;k<nfib;k++){
     				int kk = 2+3*(k);
      				xx.x=x[k*3];
      				xx.y=x[k*3+1];
      				xx.z=x[k*3+2];		
      				sig += fs[k]*Aniso_terms[k];
     				J[1] += params[0]*fs[k]*anisoterm_lambda_PVM_single_c(dir_iter,params[1],xx,bvecs,bvals,ndirections);
     				J[kk] = 0;
      				for (int j=0;j<nfib;j++){
					if(f_deriv[j*nfib+k]!=0){
	  					J[kk] += params[0]*(Aniso_terms[j]-Iso_term)*f_deriv[j*nfib+k]; 
					}
      				}
      				J[kk+1] = params[0]*fs[k]*anisoterm_th_PVM_single_c(dir_iter,_d,xx,params[kk+1],params[kk+2],bvecs,bvals,ndirections);
      				J[kk+2] = params[0]*fs[k]*anisoterm_ph_PVM_single_c(dir_iter,_d,xx,params[kk+1],params[kk+2],bvecs,bvals,ndirections);
    			}
    			if(m_include_f0){
				//partial_fsum ///////////
    				double partial_fsum=1.0;
    				for(int j=0;j<(nfib);j++)
					partial_fsum-=fs[j];
				//////////////////////////
				double temp_f0=beta2f_gpu(params[nparams-1])*partial_fsum;

    				//derivative with respect to f0
    				J[nparams-1]= params[0]*(1-Iso_term)*sin(double(2*params[nparams-1]))*partial_fsum; 
				sig=params[0]*((temp_f0+(1-*sumf-temp_f0)*Iso_term)+sig);
    				J[1] += params[0]*(1-*sumf-temp_f0)*isoterm_lambda_PVM_single_c(dir_iter,params[1],bvals);
    			}else{
				sig = params[0]*((1-*sumf)*Iso_term+sig);
	    			J[1] += params[0]*(1-*sumf)*isoterm_lambda_PVM_single_c(dir_iter,params[1],bvals);
    			}
    			diff = sig - mdata[dir_iter];
    			J[0] = sig/params[0]; 
		}

		for (int p=0;p<nparams;p++){ 
			reduction[idSubVOX]=2*J[p]*diff;

			__syncthreads();
			if(idSubVOX==0){
				for(int i=0;i<THREADS_BLOCK_FIT;i++){
					grad[p] += reduction[i];
				}
			}
			__syncthreads(); 
		} 
		dir_iter+=THREADS_BLOCK_FIT;
  	}
}


//hessian function PVM_single_c
__device__ void hess_PVM_single_c(	//INPUT
					const double*		params,
					const double*		bvecs, 
					const double*		bvals,
					const int 		ndirections,
					const int		nfib,
					const int 		nparams,
					const bool 		m_include_f0,
					const int		idSubVOX,		
					double*			reduction,	//shared memory
					double* 		fs,		//shared memory
					double*			f_deriv,	//shared memory
					double*			x,		//shared memory
					double* 		_d,		//shared memory
					double* 		sumf,		//shared memory
					//OUTPUT
					double*			hess)
{
	if(idSubVOX<nfib){
		int kk = 2+3*(idSubVOX);
		double sinth,costh,sinph,cosph;
		sincos(params[kk+1],&sinth,&costh);
		sincos(params[kk+2],&sinph,&cosph);
    		x[idSubVOX*3] = sinth*cosph;
    		x[idSubVOX*3+1] = sinth*sinph;
    		x[idSubVOX*3+2] = costh;
  	}
	if(idSubVOX==0){
		*_d = lambda2d_gpu(params[1]);
		*sumf=0;
		double partial_fsum;
		for(int k=0;k<nfib;k++){
			int kk = 2+3*(k);
    			//partial_fsum ///////////
			partial_fsum=1.0;
			for(int j=0;j<k;j++)
				partial_fsum-=fs[j];
    			//////////////////////////
			fs[k] = beta2f_gpu(params[kk])*partial_fsum;
    			*sumf += fs[k];
		}
		for (int p=0;p<nparams;p++){
			for (int p2=0;p2<nparams;p2++){ 
				hess[p*nparams+p2] = 0;
			}
		}
	}

	__syncthreads();

  	if(idSubVOX<nfib){ 
		fractions_deriv_PVM_single_c(params,fs,nfib,idSubVOX,f_deriv); 
	} 

  	int ndir = ndirections/THREADS_BLOCK_FIT;
	if(idSubVOX<(ndirections%THREADS_BLOCK_FIT)) ndir++;
	int max_dir = ndirections/THREADS_BLOCK_FIT;
	if(ndirections%THREADS_BLOCK_FIT) max_dir++;

	double J[MAXNPARAMS];
  	double sig;
	double Iso_term;
	double3 xx;
	int dir_iter=idSubVOX;
  	double Aniso_terms[MAXNFIBRES];

	__syncthreads();

  	for(int dir=0;dir<max_dir;dir++){
		for (int p=0; p<nparams; p++) J[p]=0;
		if(dir<ndir){
    			Iso_term=isoterm_PVM_single_c(dir_iter,_d,bvals);  //Precompute some terms for this datapoint
    			for(int k=0;k<nfib;k++){
      				xx.x=x[k*3];
      				xx.y=x[k*3+1];
      				xx.z=x[k*3+2];	
      				Aniso_terms[k]=anisoterm_PVM_single_c(dir_iter,_d,xx,bvecs,bvals,ndirections);
    			}
    			sig = 0;
    			for(int k=0;k<nfib;k++){
      				int kk = 2+3*(k);
      				xx.x=x[k*3];
      				xx.y=x[k*3+1];
      				xx.z=x[k*3+2];		 
      				sig += fs[k]*Aniso_terms[k];
      				J[1] += params[0]*fs[k]*anisoterm_lambda_PVM_single_c(dir_iter,params[1],xx,bvecs,bvals,ndirections);	 
      				J[kk] = 0;
      				for (int j=0; j<nfib; j++){
					if (f_deriv[j*nfib+k]!=0)
	  				J[kk] += params[0]*(Aniso_terms[j]-Iso_term)*f_deriv[j*nfib+k]; 
      				}
      				J[kk+1] = params[0]*fs[k]*anisoterm_th_PVM_single_c(dir_iter,_d,xx,params[kk+1],params[kk+2],bvecs,bvals,ndirections);
      				J[kk+2] = params[0]*fs[k]*anisoterm_ph_PVM_single_c(dir_iter,_d,xx,params[kk+1],params[kk+2],bvecs,bvals,ndirections);
    			}
    			if(m_include_f0){
				//partial_fsum ///////////
	    			double partial_fsum=1.0;
	    			for(int j=0;j<(nfib);j++)
					partial_fsum-=fs[j];
	    			//////////////////////////
    				double temp_f0=beta2f_gpu(params[nparams-1])*partial_fsum;
    				//derivative with respect to f0
    				J[nparams-1]= params[0]*(1-Iso_term)*sin(double(2*params[nparams-1]))*partial_fsum; 
				sig= params[0]*((temp_f0+(1-*sumf-temp_f0)*Iso_term)+sig);
    				J[1] += params[0]*(1-*sumf-temp_f0)*isoterm_lambda_PVM_single_c(dir_iter,params[1],bvals);
    			}else{
	    			sig = params[0]*((1-*sumf)*Iso_term+sig);
	    			J[1] += params[0]*(1-*sumf)*isoterm_lambda_PVM_single_c(dir_iter,params[1],bvals);
    			}
    			J[0] = sig/params[0]; 
		}

		for (int p=0;p<nparams;p++){
			for (int p2=p;p2<nparams;p2++){ 

				reduction[idSubVOX]=2*(J[p]*J[p2]);
				__syncthreads();
				if(idSubVOX==0){
					for(int i=0;i<THREADS_BLOCK_FIT;i++){
						hess[p*nparams+p2] += reduction[i];
					}
				}
				__syncthreads(); 
			}
		}
		dir_iter+=THREADS_BLOCK_FIT;
  	}

	if(idSubVOX==0){
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
							const int		ndirections, 
							const int 		nfib,
							const int		nparams,
							const bool		m_eval_BIC,
							const bool 		m_include_f0,
							const bool	 	m_return_fanning,
							//INPUT - OUTPUT
							double* 		params)
{
	int idSubVOX = threadIdx.x;
	int idVOX = blockIdx.x;
	int threadsBlock = blockDim.x;

	////////// DYNAMIC SHARED MEMORY ///////////
	extern __shared__ double shared[];
	double* reduction = (double*)shared;				//threadsBlock
	double* myparams = (double*) &reduction[threadsBlock];		//nparams
	double* grad = (double*) &myparams[nparams];			//nparams      
   	double* hess = (double*) &grad[nparams];			//nparams*nparams   
	double* step = (double*) &hess[nparams*nparams];		//nparams      
 	double* inverse = (double*) &step[nparams];			//nparams   
	double* pcf = (double*) &inverse[nparams];			//1   
	double* ncf = (double*) &pcf[1];				//1   
	double* lambda = (double*) &ncf[1];				//1  
	double* cftol = (double*) &lambda[1];				//1  
	double* ltol = (double*) &cftol[1];				//1  
	double* olambda = (double*) &ltol[1];				//1  

	double* fs = (double*) &olambda[1];				//nfib
	double* f_deriv = (double*) &fs[nfib];				//nfib*nfib
  	double* x = (double*) &f_deriv[nfib*nfib];			//nfib*3
	double* _d = (double*) &x[nfib*3];				//1
  	double* sumf = (double*) &_d[1];				//1

	double* C = (double*)&sumf[1];					//nparams*nparams;
	double* el =  (double*)&C[nparams*nparams];			//nparams

	int* indx = (int*)&el[nparams];					//nparams
	int* success = (int*) &indx[nparams];				//1
	int* end = (int*) &success[1];					//1   
	////////// DYNAMIC SHARED MEMORY ///////////

	if(idSubVOX<nparams){
		myparams[idSubVOX]=params[(idVOX*nparams)+idSubVOX];
	}

	__syncthreads();

	//do the fit
	levenberg_marquardt_PVM_single_c_gpu(&data[idVOX*ndirections],&bvecs[idVOX*3*ndirections],&bvals[idVOX*ndirections],ndirections,nfib,nparams,m_include_f0,idSubVOX,step,grad,hess,inverse, pcf,ncf,lambda,cftol,ltol,olambda,success,end,reduction,fs,f_deriv,x,_d,sumf,C,el,indx,myparams);

	__syncthreads();

	// finalise parameters
	// m_s0-myparams[0] 	m_d-myparams[1] 	m_f-m_th-m_ph-myparams[2,3,4,5, etc..]   	m_f0-myparams[nparams-1]
	
	if(idSubVOX==0){
  		myparams[1] = lambda2d_gpu(myparams[1]); 
  		for(int k=0;k<nfib;k++){
    			int kk = 2 + 3*(k);
    			//partial_fsum ///////////
			double partial_fsum=1.0;
			for(int j=0;j<k;j++)
				partial_fsum-=myparams[2 + 3*j];
    			//////////////////////////
    			myparams[kk]  = beta2f_gpu(myparams[kk])*partial_fsum;
  		}
  
  		if (m_include_f0){
			//partial_fsum ///////////
	    		double partial_fsum=1.0;
	    		for(int j=0;j<(nfib);j++){
				partial_fsum-=myparams[2 + 3*j];
			}
	    		//////////////////////////
    			myparams[nparams-1]= beta2f_gpu(myparams[nparams-1])*partial_fsum;
		}
		sort_PVM_single_c(nfib,myparams);
	}
	__syncthreads();

	if(idSubVOX<nparams){
		params[(idVOX*nparams)+idSubVOX] = myparams[idSubVOX];
	}
}

//in diffmodel.cc
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
								const bool* 		includes_f0,
								//OUTPUT
								double*			residuals)
{
	int idSubVOX = threadIdx.x;
	int idVOX = blockIdx.x;
	int threadsBlock = blockDim.x;

	////////// DYNAMIC SHARED MEMORY ///////////
	extern __shared__ double shared[];
	double* myparams = (double*) shared;			//nparams
	double* fs = (double*) &myparams[nparams];		//nfib
  	double* x = (double*) &fs[nfib];			//nfib*3
	double* _d = (double*) &x[nfib*3];			//1
  	double* sumf = (double*) &_d[1];			//1
	int* my_include_f0 = (int*) &sumf[1];			//1		
	////////// DYNAMIC SHARED MEMORY ///////////
	
	double val; 
	double predicted_signal;
	double mydata;
	

	if(idSubVOX==0){
		*my_include_f0 = includes_f0[idVOX];

		//m_s0-myparams[0]  m_d-myparams[1] m_f-m_th-m_ph-myparams[2,3,4,5 etc..]  m_f0-myparams[nparams-1]
  		
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
  		if (*my_include_f0){
			//partial_fsum ///////////
			partial_fsum=1.0;
			for(int j=0;j<nfib;j++)
				partial_fsum-=fs[j];
	     		//////////////////////////	
			double tmpr=params[(idVOX*nparams)+nparams-1]/partial_fsum;
    			if (tmpr>1.0) tmpr=1; //This can be due to numerical errors..asin
			if (tmpr<0.0) tmpr=0; //This can be due to numerical errors..sqrt
    			myparams[nparams-1]= f2beta_gpu(tmpr);	
		}
	}

	__syncthreads();

	if(idSubVOX<nfib){
		int kk = 2+3*idSubVOX;
		double sinth,costh,sinph,cosph;
		sincos(myparams[kk+1],&sinth,&costh);
		sincos(myparams[kk+2],&sinph,&cosph);
    		x[idSubVOX*3] = sinth*cosph;
    		x[idSubVOX*3+1] = sinth*sinph;
    		x[idSubVOX*3+2] = costh;
  	}

	if(idSubVOX==0){
		double partial_fsum;	
		*sumf=0;
		for(int k=0;k<nfib;k++){
    			int kk = 2+3*k;
			////// partial_fsum //////
			partial_fsum=1.0;
			for(int j=0;j<k;j++)
				partial_fsum-=fs[j];
    			//////////////////////////
	    		fs[k] = beta2f_gpu(myparams[kk])*partial_fsum;
	    		*sumf += fs[k];
		}
		*_d = lambda2d_gpu(myparams[1]);
	}

	int ndir = ndirections/threadsBlock;
	if(idSubVOX<(ndirections%threadsBlock)) ndir++;
	
	double3 x2;
	int dir_iter=idSubVOX; 

	__syncthreads();

	for(int dir=0;dir<ndir;dir++){
		mydata = data[(idVOX*ndirections)+dir_iter];
		predicted_signal=0;	//pred = 0;
		val = 0.0;
    		for(int k=0;k<nfib;k++){
			x2.x=x[k*3];
			x2.y=x[k*3+1];
			x2.z=x[k*3+2];	 
      			val += fs[k]*anisoterm_PVM_single_c(dir_iter,_d,x2,&bvecs[idVOX*3*ndirections],&bvals[idVOX*ndirections],ndirections);
    		}	
    		if (*my_include_f0){
			//partial_fsum ///////////
			double partial_fsum=1.0;
			for(int j=0;j<nfib;j++)
				partial_fsum-=fs[j];
	     		//////////////////////////
      			double temp_f0= beta2f_gpu(myparams[nparams-1])*partial_fsum;
      			predicted_signal = myparams[0]*(temp_f0+(1-*sumf-temp_f0)*isoterm_PVM_single_c(dir_iter,_d,&bvals[idVOX*ndirections])+val);
    		}else{
      			predicted_signal = myparams[0]*((1-*sumf)*isoterm_PVM_single_c(dir_iter,_d,&bvals[idVOX*ndirections])+val); 
		}

		//residuals=m_data-predicted_signal;
		residuals[idVOX*ndirections+dir_iter]= mydata - predicted_signal;

		dir_iter+=threadsBlock;
	}
}
