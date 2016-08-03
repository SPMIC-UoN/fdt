/*  PVM_single.cu

    Tim Behrens, Saad Jbabdi, Stam Sotiropoulos, Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */

#include "diffmodels_utils.h"
#include "levenberg_marquardt.cu"
#include "options.h"

//#include <fstream>

/////////////////////////////////////
/////////////////////////////////////
/// 	    PVM_single		  /// 
/////////////////////////////////////
/////////////////////////////////////

__device__ 
inline float isoterm_PVM_single(const int pt,const float* _d,const float *bvals){
  	return expf(-bvals[pt]**_d);
}

__device__ 
inline float isoterm_d_PVM_single(const int pt,const float* _d,const float *bvals){
  	return (-bvals[pt]*expf(-bvals[pt]**_d));
}

__device__ 
inline float anisoterm_PVM_single(const int pt,const float* _d,const float3 x, const float *bvecs, const float *bvals, const int ndirections){
	float dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
	return expf(-bvals[pt]**_d*dp*dp);
}

__device__ 
inline float anisoterm_d_PVM_single(const int pt,const float* _d,const float3 x,const float *bvecs, const float *bvals, const int ndirections){
	float dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
  	return(-bvals[pt]*dp*dp*expf(-bvals[pt]**_d*dp*dp));
}

__device__ 
inline float anisoterm_th_PVM_single(const int pt,const float* _d,const float3 x, const float _th,const float _ph,const float *bvecs, const float *bvals, const int ndirections){
	float sinth,costh,sinph,cosph;
	sincosf(_th,&sinth,&costh);
	sincosf(_ph,&sinph,&cosph);
	float dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
	float dp1 = (costh*(bvecs[pt]*cosph+bvecs[ndirections+pt]*sinph)-bvecs[(2*ndirections)+pt]*sinth);
  	return(-2*bvals[pt]**_d*dp*dp1*expf(-bvals[pt]**_d*dp*dp));
}

__device__ 
inline float anisoterm_ph_PVM_single(const int pt,const float* _d,const float3 x, const float _th,const float _ph,const float *bvecs, const float *bvals, const int ndirections){
	float sinth,sinph,cosph;
	sinth=sinf(_th);
	sincosf(_ph,&sinph,&cosph);
  	float dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
	float dp1 = sinth*(-bvecs[pt]*sinph+bvecs[ndirections+pt]*cosph);
  	return(-2.0f*bvals[pt]**_d*dp*dp1*expf(-bvals[pt]**_d*dp*dp));
}

//in diffmodel.cc
__device__ void fix_fsum_PVM_single(	//INPUT 
					bool m_include_f0, 
					int nfib,
					int nparams,
					//INPUT - OUTPUT){
					float *params)
{
  	float sum=0.0f;
  	if (m_include_f0) 
    		sum=params[nparams-1];
  	for(int i=0;i<nfib;i++){
    		sum += params[2+(i*3)];
    		if(sum>=1.0f){
			for(int j=i;j<nfib;j++)
				params[2+(j*3)]=FSMALL_gpu; 
			break;
		}
  	}
}

//in diffmodel.cc
__device__  void sort_PVM_single(int nfib,float* params)
{
	float temp_f, temp_th, temp_ph;
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

//cost function PVM_single
__device__ void cf_PVM_single(	//INPUT
				const float*		params,
				const float*		mdata,
				const float*		bvecs, 
				const float*		bvals,
				const int		ndirections,
				const int		nfib,
				const int 		nparams,
				const bool 		m_include_f0,
				const int		idSubVOX,
				float*			reduction,	//shared memory
				float* 			fs,		//shared memory
				float*			x,		//shared memory	
				float* 			_d,		//shared memory
				float* 			sumf,		//shared memory
				//OUTPUT
				double*			cfv)
{
	if(idSubVOX<nfib){
		int kk = 2+3*(idSubVOX);
		float sinth,costh,sinph,cosph;
		sincosf(params[kk+1],&sinth,&costh);
		sincosf(params[kk+2],&sinph,&cosph);
    		fs[idSubVOX] = x2f_gpu(params[kk]);
		x[idSubVOX*3] = sinth*cosph;
    		x[idSubVOX*3+1] = sinth*sinph;
    		x[idSubVOX*3+2] = costh;
  	}

	__syncthreads(); 	

	if(idSubVOX==0){
		*_d = fabsf(params[1]);
		*cfv = 0.0f;
		*sumf=0.0f;
		for(int k=0;k<nfib;k++) *sumf+= fs[k];
	}
	
	int ndir = ndirections/THREADS_BLOCK_FIT;
	if(idSubVOX<(ndirections%THREADS_BLOCK_FIT)) ndir++;
	
	float err;
	float3 x2;
	int dir_iter=idSubVOX;

	__syncthreads();

	reduction[idSubVOX]=0.0f;
	for(int dir=0;dir<ndir;dir++){
		err = 0.0f;
    		for(int k=0;k<nfib;k++){
			x2.x=x[k*3];
			x2.y=x[k*3+1];
			x2.z=x[k*3+2];	
			err += fs[k]*anisoterm_PVM_single(dir_iter,_d,x2,bvecs,bvals,ndirections); 
    		}
		if(m_include_f0){
			float temp_f0=x2f_gpu(params[nparams-1]);
			err= (params[0]*((temp_f0+(1.0f-*sumf-temp_f0)*isoterm_PVM_single(dir_iter,_d,bvals))+err))-mdata[dir_iter];
		}else{
			err =  (params[0]*((1.0f-*sumf)*isoterm_PVM_single(dir_iter,_d,bvals)+err))-mdata[dir_iter];
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

//gradient function PVM_single
__device__ void grad_PVM_single(	//INPUT
					const float*		params,
					const float*		mdata,
					const float*		bvecs, 
					const float*		bvals,
					const int 		ndirections,
					const int		nfib,
					const int 		nparams,
					const bool 		m_include_f0,
					const int		idSubVOX,	
					float*			J,		//shared memory	
					float*			reduction,	//shared memory
					float* 			fs,		//shared memory
					float*			x,		//shared memory
					float* 			_d,		//shared memory
					float* 			sumf,		//shared memory
					//OUTPUT
					float*			grad)
{
	if(idSubVOX<nfib){
		int kk = 2+3*(idSubVOX);
		float sinth,costh,sinph,cosph;
		sincosf(params[kk+1],&sinth,&costh);
		sincosf(params[kk+2],&sinph,&cosph);
    		fs[idSubVOX] = x2f_gpu(params[kk]);
    		x[idSubVOX*3] = sinth*cosph;
    		x[idSubVOX*3+1] = sinth*sinph;
    		x[idSubVOX*3+2] = costh;
  	}

	__syncthreads(); 

	if(idSubVOX==0){
		*sumf=0.0f;
		for(int i=0;i<nfib;i++) *sumf+=fs[i];
		*_d = fabsf(params[1]);
		for (int p=0;p<nparams;p++) grad[p]=0.0f;
	}

	int ndir = ndirections/THREADS_BLOCK_FIT;
	if(idSubVOX<(ndirections%THREADS_BLOCK_FIT)) ndir++;
	int max_dir = ndirections/THREADS_BLOCK_FIT;
	if(ndirections%THREADS_BLOCK_FIT) max_dir++;

	float* myJ = &J[idSubVOX*nparams];
	float diff;
  	float sig;
	float3 xx;
	int dir_iter=idSubVOX;

	__syncthreads();

  	for(int dir=0;dir<max_dir;dir++){
		for (int p=0; p<nparams; p++) myJ[p]=0.0f;
		if(dir<ndir){
    			sig = 0.0f;
    			for(int k=0;k<nfib;k++){
      				int kk = 2+3*(k);
      				xx.x=x[k*3];
      				xx.y=x[k*3+1];
      				xx.z=x[k*3+2];			
				sig +=  fs[k]*anisoterm_PVM_single(dir_iter,_d,xx,bvecs,bvals,ndirections);
				myJ[1] +=  (params[1]>0.0f?1.0f:-1.0f)*params[0]*fs[k]*anisoterm_d_PVM_single(dir_iter,_d,xx,bvecs,bvals,ndirections);
      				myJ[kk] = params[0]*(anisoterm_PVM_single(dir_iter,_d,xx,bvecs,bvals,ndirections)-isoterm_PVM_single(dir_iter,_d,bvals)) * two_pi_gpu*sign_gpu(params[kk])*1.0f/(1.0f+params[kk]*params[kk]);
      				myJ[kk+1] = params[0]*fs[k]*anisoterm_th_PVM_single(dir_iter,_d,xx,params[kk+1],params[kk+2],bvecs,bvals,ndirections);
      				myJ[kk+2] = params[0]*fs[k]*anisoterm_ph_PVM_single(dir_iter,_d,xx,params[kk+1],params[kk+2],bvecs,bvals,ndirections);
    			}

    			if(m_include_f0){
				float temp_f0=x2f_gpu(params[nparams-1]);
				myJ[nparams-1]= params[0]*(1.0f-isoterm_PVM_single(dir_iter,_d,bvals))* two_pi_gpu*sign_gpu(params[nparams-1])*1.0f/(1.0f+params[nparams-1]*params[nparams-1]);
				sig= params[0]*((temp_f0+(1.0f-*sumf-temp_f0)*isoterm_PVM_single(dir_iter,_d,bvals))+sig);
    				myJ[1] += (params[1]>0.0f?1.0f:-1.0f)*params[0]*(1.0f-*sumf-temp_f0)*isoterm_d_PVM_single(dir_iter,_d,bvals);
    			}else{
				sig = params[0]*((1.0f-*sumf)*isoterm_PVM_single(dir_iter,_d,bvals)+sig);
				myJ[1] += (params[1]>0.0f?1.0f:-1.0f)*params[0]*(1.0f-*sumf)*isoterm_d_PVM_single(dir_iter,_d,bvals);
    			}
    			diff = sig - mdata[dir_iter];
    			myJ[0] = sig/params[0];
		}

		for (int p=0;p<nparams;p++){ 
			reduction[idSubVOX]=2.0f*myJ[p]*diff;

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

//hessian function PVM_single
__device__ void hess_PVM_single(	//INPUT
					const float*		params,
					const float*		bvecs, 
					const float*		bvals,
					const int		ndirections,
					const int 		nfib,
					const int 		nparams,
					const bool 		m_include_f0,
					const int		idSubVOX,
					float*			J,		//shared memory
					float*			reduction,	//shared memory					
					float* 			fs,		//shared memory
					float*			x,		//shared memory
					float* 			_d,		//shared memory
					float* 			sumf,		//shared memory
					//OUTPUT
					float*			hess)
{
	if(idSubVOX<nfib){
		int kk = 2+3*(idSubVOX);
		float sinth,costh,sinph,cosph;
		sincosf(params[kk+1],&sinth,&costh);
		sincosf(params[kk+2],&sinph,&cosph);
    		fs[idSubVOX] = x2f_gpu(params[kk]);
    		x[idSubVOX*3] = sinth*cosph;
    		x[idSubVOX*3+1] = sinth*sinph;
    		x[idSubVOX*3+2] = costh;
  	}

	__syncthreads(); 

	if(idSubVOX==0){
		*sumf=0.0f;
		for(int i=0;i<nfib;i++) *sumf+=fs[i];
		*_d = fabsf(params[1]);
		for (int p=0;p<nparams;p++){
			for (int p2=0;p2<nparams;p2++){ 
				hess[p*nparams+p2] = 0.0f;
			}
		}
	}

	int ndir = ndirections/THREADS_BLOCK_FIT;
	if(idSubVOX<(ndirections%THREADS_BLOCK_FIT)) ndir++;
	int max_dir = ndirections/THREADS_BLOCK_FIT;
	if(ndirections%THREADS_BLOCK_FIT) max_dir++;

	float* myJ = &J[idSubVOX*nparams];
  	float sig;
	float3 xx;
	int dir_iter=idSubVOX; 

	__syncthreads(); 
	
  	for(int dir=0;dir<max_dir;dir++){
		for (int p=0; p<nparams; p++) myJ[p]=0.0f;
		if(dir<ndir){
    			sig = 0.0f;
    			for(int k=0;k<nfib;k++){
      				int kk = 2+3*(k);
      				xx.x=x[k*3];
      				xx.y=x[k*3+1];
      				xx.z=x[k*3+2];		
				sig += fs[k]*anisoterm_PVM_single(dir_iter,_d,xx,bvecs,bvals,ndirections);
      				myJ[1] += (params[1]>0.0f?1.0f:-1.0f)*params[0]*fs[k]*anisoterm_d_PVM_single(dir_iter,_d,xx,bvecs,bvals,ndirections);
      				myJ[kk] = params[0]*(anisoterm_PVM_single(dir_iter,_d,xx,bvecs,bvals,ndirections)-isoterm_PVM_single(dir_iter,_d,bvals)) * two_pi_gpu*sign_gpu(params[kk])*1.0f/(1.0f+params[kk]*params[kk]);
		      		myJ[kk+1] = params[0]*fs[k]*anisoterm_th_PVM_single(dir_iter,_d,xx,params[kk+1],params[kk+2],bvecs,bvals,ndirections);
		      		myJ[kk+2] = params[0]*fs[k]*anisoterm_ph_PVM_single(dir_iter,_d,xx,params[kk+1],params[kk+2],bvecs,bvals,ndirections);
    			}	

    			if(m_include_f0){
				float temp_f0=x2f_gpu(params[nparams-1]);
				myJ[nparams-1]= params[0]*(1.0f-isoterm_PVM_single(dir_iter,_d,bvals))* two_pi_gpu*sign_gpu(params[nparams-1])*1.0f/(1.0f+params[nparams-1]*params[nparams-1]);
				sig=params[0]*((temp_f0+(1.0f-*sumf-temp_f0)*isoterm_PVM_single(dir_iter,_d,bvals))+sig);
    				myJ[1] += (params[1]>0.0f?1.0f:-1.0f)*params[0]*(1.0f-*sumf-temp_f0)*isoterm_d_PVM_single(dir_iter,_d,bvals);	
    			}else{
				sig = params[0]*((1.0f-*sumf)*isoterm_PVM_single(dir_iter,_d,bvals)+sig);
	    			myJ[1] +=  (params[1]>0.0f?1.0f:-1.0f)*params[0]*(1.0f-*sumf)*isoterm_d_PVM_single(dir_iter,_d,bvals);
    			}   
    			myJ[0] = sig/params[0];
		}

		for (int p=0;p<nparams;p++){
			for (int p2=p;p2<nparams;p2++){ 

				reduction[idSubVOX]=2.0f*(myJ[p]*myJ[p2]);
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
extern "C" __global__ void fit_PVM_single_kernel(	//INPUT
							const float* 		data, 
							const float* 		bvecs,
							const float* 		bvals, 
							const int 		nvox, 
							const int		ndirections,
							const int 		nfib, 
							const int		nparams,
							const bool 		m_include_f0,
							const bool		gradnonlin,
							//INPUT-OUTPUT
							float* 		params)
{
	int idSubVOX = threadIdx.x;
	int idVOX = blockIdx.x;
	int threadsBlock = blockDim.x;

	////////// DYNAMIC SHARED MEMORY ///////////
	extern __shared__ double shared[];
	double* pcf = (double*) shared;					//1   
	double* ncf = (double*) &pcf[1];				//1   
	double* lambda = (double*) &ncf[1];				//1  
	double* cftol = (double*) &lambda[1];				//1  
	double* ltol = (double*) &cftol[1];				//1  
	double* olambda = (double*) &ltol[1];				//1  

	float* J = (float*)&olambda[1];					//threadsBlock*nparams
	float* reduction = (float*)&J[threadsBlock*nparams];		//threadsBlock
	float* myparams = (float*) &reduction[threadsBlock];		//nparams
	float* grad = (float*) &myparams[nparams];			//nparams      
   	float* hess = (float*) &grad[nparams];				//nparams*nparams   
	float* step = (float*) &hess[nparams*nparams];			//nparams      
 	float* inverse = (float*) &step[nparams];			//nparams   

	float* fs = (float*) &inverse[nparams];				//nfib
  	float* x = (float*) &fs[nfib];					//nfib*3
	float* _d = (float*) &x[nfib*3];				//1
  	float* sumf = (float*) &_d[1];					//1

	float* C = (float*)&sumf[1];					//nparams*nparams;
	float* el =  (float*)&C[nparams*nparams];			//nparams

	int* indx = (int*)&el[nparams];					//nparams
	int* success = (int*) &indx[nparams];				//1
	int* end = (int*) &success[1];					//1    
	////////// DYNAMIC SHARED MEMORY ///////////

	if(idSubVOX<nparams){
		myparams[idSubVOX]=params[(idVOX*nparams)+idSubVOX];
	}

	__syncthreads();

	int pos_bvals, pos_bvecs;
	if(gradnonlin){ 
		pos_bvals=idVOX*ndirections;
		pos_bvecs=idVOX*3*ndirections;
	}else{ 
		pos_bvals=0;
		pos_bvecs=0;
	}
	// do the fit
	levenberg_marquardt_PVM_single_gpu(&data[idVOX*ndirections],&bvecs[pos_bvecs],&bvals[pos_bvals],ndirections,nfib,nparams,m_include_f0,idSubVOX,step,grad,hess,inverse, pcf,ncf,lambda,cftol,ltol,olambda,success,end,J,reduction,fs,x,_d,sumf,C,el,indx,myparams);

	__syncthreads();
	
  	// finalise parameters
	//m_s0 in myparams[0] 	m_d in myparams[1] 	m_f-m_th-m_ph in myparams[2,3,4,5, etc..]   	m_f0 in myparams[nparams-1]

	if(idSubVOX==0){  	
  		myparams[1] = fabsf(myparams[1]); 
  		for(int k=1;k<=nfib;k++){
    			int kk = 2 + 3*(k-1);
    			myparams[kk] = x2f_gpu(myparams[kk]);
  		}
  		if(m_include_f0)
    			myparams[nparams-1]=x2f_gpu(myparams[nparams-1]);

  		sort_PVM_single(nfib,myparams);
  		fix_fsum_PVM_single(m_include_f0,nfib,nparams,myparams);
	}
	__syncthreads();

	if(idSubVOX<nparams){
		params[idVOX*nparams+idSubVOX]=myparams[idSubVOX];	
	}
}

//in diffmodel.cc
extern "C" __global__ void get_residuals_PVM_single_kernel(	//INPUT
								const float* 		data, 
								const float* 		params,
								const float* 		bvecs, 
								const float* 		bvals, 
								const int 		nvox, 
								const int		ndirections,
								const int 		nfib, 
								const int		nparams,
								const bool 		m_include_f0,
								const bool		gradnonlin,
								const bool* 		includes_f0,
								//OUTPUT
								float*			residuals)
{
	int idSubVOX = threadIdx.x;
	int idVOX = blockIdx.x;
	int threadsBlock = blockDim.x;

	////////// DYNAMIC SHARED MEMORY ///////////
	extern __shared__ double shared[];
	float* myparams = (float*) shared;			//nparams
	float* fs = (float*) &myparams[nparams];		//nfib
  	float* x = (float*) &fs[nfib];				//nfib*3
	float* _d = (float*) &x[nfib*3];			//1
  	float* sumf = (float*) &_d[1];				//1
	int* my_include_f0 = (int*) &sumf[1];			//1	
	////////// DYNAMIC SHARED MEMORY ///////////

	float val;
	float predicted_signal;
	float mydata;

	if(idSubVOX==0){
		*my_include_f0 = includes_f0[idVOX];

		//m_s0-myparams[0]  m_d-myparams[1]  m_f-m_th-m_ph-myparams[2,3,4,5 etc..]  m_f0-myparams[nparams-1]

		myparams[0]=params[(idVOX*nparams)+0];
		myparams[1]=params[(idVOX*nparams)+1];

  		if (*my_include_f0)
    			myparams[nparams-1]=f2x_gpu(params[(idVOX*nparams)+nparams-1]);
	}

	if(idSubVOX<nfib){
		int kk = 2+3*idSubVOX;
		float sinth,costh,sinph,cosph;

		myparams[kk]   = f2x_gpu(params[(idVOX*nparams)+kk]);
    		myparams[kk+1] = params[(idVOX*nparams)+kk+1];
    		myparams[kk+2] = params[(idVOX*nparams)+kk+2];

		sincosf(myparams[kk+1],&sinth,&costh);
		sincosf(myparams[kk+2],&sinph,&cosph);

    		fs[idSubVOX] = x2f_gpu(myparams[kk]);
    		x[idSubVOX*3] = sinth*cosph;
    		x[idSubVOX*3+1] = sinth*sinph;
    		x[idSubVOX*3+2] = costh;
  	}

	__syncthreads(); 

	if(idSubVOX==0){
		*sumf=0.0f;
		for(int i=0;i<nfib;i++) *sumf+=fs[i];
		*_d = fabsf(myparams[1]);
	}
  	
  	int ndir = ndirections/threadsBlock;
	if(idSubVOX<(ndirections%threadsBlock)) ndir++;
	
	float3 x2;
	int dir_iter=idSubVOX; 

	__syncthreads();

	int pos_bvals, pos_bvecs;
	if(gradnonlin){ 
		pos_bvals=idVOX*ndirections;
		pos_bvecs=idVOX*3*ndirections;
	}else{ 
		pos_bvals=0;
		pos_bvecs=0;
	}

	for(int dir=0;dir<ndir;dir++){
		mydata = data[(idVOX*ndirections)+dir_iter];
  		predicted_signal=0.0f;	//pred = 0;
    		val = 0.0f;
    		for(int k=0;k<nfib;k++){
			x2.x=x[k*3];
			x2.y=x[k*3+1];
			x2.z=x[k*3+2];	 
      			val += fs[k]*anisoterm_PVM_single(dir_iter,_d,x2,&bvecs[pos_bvecs],&bvals[pos_bvals],ndirections);
    		}	
    		if (*my_include_f0){
      			float temp_f0=x2f_gpu(myparams[nparams-1]);
      			predicted_signal = myparams[0]*(temp_f0+(1.0f-*sumf-temp_f0)*isoterm_PVM_single(dir_iter,_d,&bvals[pos_bvals])+val);
    		}else{
      			predicted_signal = myparams[0]*((1.0f-*sumf)*isoterm_PVM_single(dir_iter,_d,&bvals[pos_bvals])+val); 
		}
	
		//residuals=m_data-predicted_signal;
		residuals[idVOX*ndirections+dir_iter]= mydata - predicted_signal;

		dir_iter+=threadsBlock;
  	}
}

