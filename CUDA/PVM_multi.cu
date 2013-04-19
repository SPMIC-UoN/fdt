/*  PVM_multi.cu

    Tim Behrens, Saad Jbabdi, Stam Sotiropoulos, Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */

#include "diffmodels_utils.h"
#include "levenberg_marquardt.cu"
#include "options.h"

/////////////////////////////////////
/////////////////////////////////////
/// 	    PVM_multi	 	  /// 
/////////////////////////////////////
/////////////////////////////////////

__device__ inline double isoterm_PVM_multi(const int pt,const double* _a,const double* _b, const double *bvals){
	return exp(-*_a*log(1+bvals[pt]**_b));
}

__device__ inline double isoterm_a_PVM_multi(const int pt,const double* _a,const double* _b, const double *bvals){
    	return  -log(1+bvals[pt]**_b)*exp(-*_a*log(1+bvals[pt]**_b));
}

__device__ inline double isoterm_b_PVM_multi(const int pt,const double* _a,const double* _b, const double *bvals){
      	return -*_a*bvals[pt]/(1+bvals[pt]**_b)*exp(-*_a*log(1+bvals[pt]**_b));
}

__device__ inline double anisoterm_PVM_multi(const int pt,const double* _a,const double* _b,const double3 x,const double *bvecs, const double *bvals, const int ndirections){
	double dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
  	return exp(-*_a*log(1+bvals[pt]**_b*(dp*dp)));
}
 
__device__ inline double anisoterm_a_PVM_multi(const int pt,const double* _a,const double* _b,const double3 x,const double *bvecs, const double *bvals, const int ndirections){
	double dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
  	return -log(1+bvals[pt]*(dp*dp)**_b)* exp(-*_a*log(1+bvals[pt]*(dp*dp)**_b));
}

__device__ inline double anisoterm_b_PVM_multi(const int pt,const double* _a,const double* _b,const double3 x,const double *bvecs, const double *bvals, const int ndirections){
  	double dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
  	return (-*_a*bvals[pt]*(dp*dp)/ (1+bvals[pt]*(dp*dp)**_b)*exp(-*_a*log(1+bvals[pt]*(dp*dp)**_b)));
}

__device__ inline double anisoterm_th_PVM_multi(const int pt,const double* _a,const double* _b,const double3 x,const double _th,const double _ph,const double *bvecs, const double *bvals, const int ndirections){
	double sinth,costh,sinph,cosph;
	sincos(_th,&sinth,&costh);
	sincos(_ph,&sinph,&cosph);
  	double dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
  	double dp1 = costh* (bvecs[pt]*cosph + bvecs[ndirections+pt]*sinph) - bvecs[(2*ndirections)+pt]*sinth;
	return  (-*_a**_b*bvals[pt]/(1+bvals[pt]*(dp*dp)**_b)*exp(-*_a*log(1+bvals[pt]*(dp*dp)**_b))*2*dp*dp1);	
}

__device__ inline double anisoterm_ph_PVM_multi(const int pt,const double* _a,const double* _b,const double3 x,const double _th,const double _ph,const double *bvecs, const double *bvals, const int ndirections){
	double sinth,sinph,cosph;
	sinth=sin(_th);
	sincos(_ph,&sinph,&cosph);
  	double dp = bvecs[pt]*x.x+bvecs[ndirections+pt]*x.y+bvecs[(2*ndirections)+pt]*x.z;
  	double dp1 = sinth* (-bvecs[pt]*sinph + bvecs[ndirections+pt]*cosph);
  	return  (-*_a**_b*bvals[pt]/(1+bvals[pt]*(dp*dp)**_b)*exp(-*_a*log(1+bvals[pt]*(dp*dp)**_b))*2*dp*dp1);
}

//in diffmodel.cc
__device__ void fix_fsum_PVM_multi(	//INPUT 
					bool m_include_f0, 
					int nfib,
					int nparams,
					//INPUT - OUTPUT){
					double *params)
{
  	double sumf=0;
  	if (m_include_f0) 
    		sumf=params[nparams-1];
  	for(int i=0;i<nfib;i++){
    		if (params[3+(i*3)]==0) 
			params[3+(i*3)]=FSMALL_gpu;
    		sumf+=params[3+(i*3)];
    		if(sumf>=1){
			for(int j=i;j<nfib;j++)
				params[3+(j*3)]=FSMALL_gpu;
			break;
		}
	}
}

//in diffmodel.cc
__device__ void sort_PVM_multi(int nfib,double* params)
{
	double temp_f, temp_th, temp_ph;
	// Order vector descending using f parameters as index
	for(int i=1; i<(nfib); i++){ 
    		for(int j=0; j<(nfib-i); j++){ 
      			if (params[3+j*3] < params[3+(j+1)*3]){ 
        			temp_f = params[3+j*3];
				temp_th = params[3+j*3+1];
				temp_ph = params[3+j*3+2];
        			params[3+j*3] = params[3+(j+1)*3]; 
				params[3+j*3+1] = params[3+(j+1)*3+1]; 
				params[3+j*3+2] = params[3+(j+1)*3+2]; 
        			params[3+(j+1)*3] = temp_f; 
				params[3+(j+1)*3+1] = temp_th; 
				params[3+(j+1)*3+2] = temp_ph; 
      			} 
    		} 
  	} 
}

//cost function PVM_multi
__device__ void cf_PVM_multi(		//INPUT
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
					double*			x,		//shared memory	
					double* 		_a,		//shared memory
					double* 		_b,		//shared memory
					double* 		sumf,		//shared memory
					//OUTPUT
					double*			cfv)
{
	if(idSubVOX<nfib){
		int kk = 3+3*(idSubVOX);
		double sinth,costh,sinph,cosph;
		sincos(params[kk+1],&sinth,&costh);
		sincos(params[kk+2],&sinph,&cosph);
		fs[idSubVOX] = x2f_gpu(params[kk]);
		x[idSubVOX*3] = sinth*cosph;
    		x[idSubVOX*3+1] = sinth*sinph;
    		x[idSubVOX*3+2] = costh;
  	}
	if(idSubVOX==0){
		*_a= abs(params[1]);
		*_b= abs(params[2]); 
		*cfv = 0.0;
		*sumf=0;
		for(int k=0;k<nfib;k++) *sumf+= fs[k];
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
			err += fs[k]*anisoterm_PVM_multi(dir_iter,_a,_b,x2,bvecs,bvals,ndirections); 
    		}
		if(m_include_f0){
			double temp_f0=x2f_gpu(params[nparams-1]);
			err = (abs(params[0])*(temp_f0+((1-*sumf-temp_f0)*isoterm_PVM_multi(dir_iter,_a,_b,bvals)+err)))-mdata[dir_iter];
		}else{
			err = abs(params[0])*((1-*sumf)*isoterm_PVM_multi(dir_iter,_a,_b,bvals)+err)-mdata[dir_iter];
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

//gradient function PVM_multi
__device__ void grad_PVM_multi(		//INPUT
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
					double*			x,		//shared memory	
					double* 		_a,		//shared memory
					double* 		_b,		//shared memory
					double* 		sumf,		//shared memory
					//OUTPUT
					double*			grad)
{	
	if(idSubVOX<nfib){
		int kk = 3+3*(idSubVOX);
		double sinth,costh,sinph,cosph;
		sincos(params[kk+1],&sinth,&costh);
		sincos(params[kk+2],&sinph,&cosph);
		fs[idSubVOX] = x2f_gpu(params[kk]);
		x[idSubVOX*3] = sinth*cosph;
    		x[idSubVOX*3+1] = sinth*sinph;
    		x[idSubVOX*3+2] = costh;
  	}
	if(idSubVOX==0){
		*_a= abs(params[1]);
		*_b= abs(params[2]); 
		*sumf=0;
		for(int k=0;k<nfib;k++) *sumf+= fs[k];
		for (int p=0;p<nparams;p++) grad[p]=0;
	}

  	int ndir = ndirections/THREADS_BLOCK_FIT;
	if(idSubVOX<(ndirections%THREADS_BLOCK_FIT)) ndir++;
	int max_dir = ndirections/THREADS_BLOCK_FIT;
	if(ndirections%THREADS_BLOCK_FIT) max_dir++;

	double J[MAXNPARAMS];
	double diff;
  	double sig;
	double3 xx;
	int dir_iter=idSubVOX;

	__syncthreads();

  	for(int dir=0;dir<max_dir;dir++){
		for (int p=0; p<nparams; p++) J[p]=0;
		if(dir<ndir){
    			sig = 0;
    			for(int k=0;k<nfib;k++){
      				int kk = 3+3*(k);
      				xx.x=x[k*3];
      				xx.y=x[k*3+1];
      				xx.z=x[k*3+2];		
      				sig += fs[k]*anisoterm_PVM_multi(dir_iter,_a,_b,xx,bvecs,bvals,ndirections);
      				J[1] += (params[1]>0?1.0:-1.0)*abs(params[0])*fs[k]*anisoterm_a_PVM_multi(dir_iter,_a,_b,xx,bvecs,bvals,ndirections); 
				J[2] += (params[2]>0?1.0:-1.0)*abs(params[0])*fs[k]*anisoterm_b_PVM_multi(dir_iter,_a,_b,xx,bvecs,bvals,ndirections);
				J[kk] = abs(params[0])*(anisoterm_PVM_multi(dir_iter,_a,_b,xx,bvecs,bvals,ndirections)-isoterm_PVM_multi(dir_iter,_a,_b,bvals))*two_pi_gpu*sign_gpu(params[kk])*1/(1+params[kk]*params[kk]); 
      				J[kk+1] = abs(params[0])*fs[k]*anisoterm_th_PVM_multi(dir_iter,_a,_b,xx,params[kk+1],params[kk+2],bvecs,bvals,ndirections);  
      				J[kk+2] = abs(params[0])*fs[k]*anisoterm_ph_PVM_multi(dir_iter,_a,_b,xx,params[kk+1],params[kk+2],bvecs,bvals,ndirections);
    			}
    			if(m_include_f0){
				double temp_f0=x2f_gpu(params[nparams-1]);
				J[nparams-1]= abs(params[0])*(1-isoterm_PVM_multi(dir_iter,_a,_b,bvals))*two_pi_gpu*sign_gpu(params[nparams-1])*1/(1+params[nparams-1]*params[nparams-1]);
				sig=abs(params[0])*((temp_f0+(1-*sumf-temp_f0)*isoterm_PVM_multi(dir_iter,_a,_b,bvals))+sig);
    				J[1] += (params[1]>0?1.0:-1.0)*abs(params[0])*(1-*sumf-temp_f0)*isoterm_a_PVM_multi(dir_iter,_a,_b,bvals);
				J[2] += (params[2]>0?1.0:-1.0)*abs(params[0])*(1-*sumf-temp_f0)*isoterm_b_PVM_multi(dir_iter,_a,_b,bvals);
    			}else{
	    			sig = abs(params[0]) * ((1-*sumf)*isoterm_PVM_multi(dir_iter,_a,_b,bvals)+sig);
	    			J[1] += (params[1]>0?1.0:-1.0)*abs(params[0])*(1-*sumf)*isoterm_a_PVM_multi(dir_iter,_a,_b,bvals);
	    			J[2] += (params[2]>0?1.0:-1.0)*abs(params[0])*(1-*sumf)*isoterm_b_PVM_multi(dir_iter,_a,_b,bvals);	
    			}
    
    			diff = sig - mdata[dir_iter];
    			J[0] = (params[0]>0?1.0:-1.0)*sig/params[0]; 
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

//hessian function PVM_multi 
__device__ void hess_PVM_multi(		//INPUT
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
					double*			x,		//shared memory	
					double* 		_a,		//shared memory
					double* 		_b,		//shared memory
					double* 		sumf,		//shared memory
					//OUTPUT
					double*			hess)
{
	if(idSubVOX<nfib){
		int kk = 3+3*(idSubVOX);
		double sinth,costh,sinph,cosph;
		sincos(params[kk+1],&sinth,&costh);
		sincos(params[kk+2],&sinph,&cosph);
		fs[idSubVOX] = x2f_gpu(params[kk]);
		x[idSubVOX*3] = sinth*cosph;
    		x[idSubVOX*3+1] = sinth*sinph;
    		x[idSubVOX*3+2] = costh;
  	}
	if(idSubVOX==0){
		*_a= abs(params[1]);
		*_b= abs(params[2]); 
		*sumf=0;
		for(int k=0;k<nfib;k++) *sumf+= fs[k];
		for (int p=0;p<nparams;p++){
			for (int p2=0;p2<nparams;p2++){ 
				hess[p*nparams+p2] = 0;
			}
		}
	}

  	int ndir = ndirections/THREADS_BLOCK_FIT;
	if(idSubVOX<(ndirections%THREADS_BLOCK_FIT)) ndir++;
	int max_dir = ndirections/THREADS_BLOCK_FIT;
	if(ndirections%THREADS_BLOCK_FIT) max_dir++;

	double J[MAXNPARAMS];
  	double sig;
	double3 xx;
	int dir_iter=idSubVOX; 

	__syncthreads(); 
	
  	for(int dir=0;dir<max_dir;dir++){
		for (int p=0; p<nparams; p++) J[p]=0;
		if(dir<ndir){
    			sig = 0;
    			for(int k=0;k<nfib;k++){
      				int kk = 3+3*(k);
      				xx.x=x[k*3];
      				xx.y=x[k*3+1];
      				xx.z=x[k*3+2];		
      				sig += fs[k]*anisoterm_PVM_multi(dir_iter,_a,_b,xx,bvecs,bvals,ndirections);
      				double cov = two_pi_gpu*sign_gpu(params[kk])*1/(1+params[kk]*params[kk]);	
      				J[1] += (params[1]>0?1.0:-1.0)*abs(params[0])*fs[k]*anisoterm_a_PVM_multi(dir_iter,_a,_b,xx,bvecs,bvals,ndirections);
				J[2] += (params[2]>0?1.0:-1.0)*abs(params[0])*fs[k]*anisoterm_b_PVM_multi(dir_iter,_a,_b,xx,bvecs,bvals,ndirections);
				J[kk] = abs(params[0])*(anisoterm_PVM_multi(dir_iter,_a,_b,xx,bvecs,bvals,ndirections)-isoterm_PVM_multi(dir_iter,_a,_b,bvals))*cov;
      				J[kk+1] = abs(params[0])*fs[k]*anisoterm_th_PVM_multi(dir_iter,_a,_b,xx,params[kk+1],params[kk+2],bvecs,bvals,ndirections);
      				J[kk+2] = abs(params[0])*fs[k]*anisoterm_ph_PVM_multi(dir_iter,_a,_b,xx,params[kk+1],params[kk+2],bvecs,bvals,ndirections);
    			}
    			if(m_include_f0){
				double temp_f0=x2f_gpu(params[nparams-1]);
				J[nparams-1]= abs(params[0])*(1-isoterm_PVM_multi(dir_iter,_a,_b,bvals))*two_pi_gpu*sign_gpu(params[nparams-1])*1/(1+params[nparams-1]*params[nparams-1]);
	    			sig = abs(params[0])* (temp_f0+(1-*sumf-temp_f0)*isoterm_PVM_multi(dir_iter,_a,_b,bvals)+sig);
    				J[1] += (params[1]>0?1.0:-1.0)*abs(params[0])*(1-*sumf-temp_f0)*isoterm_a_PVM_multi(dir_iter,_a,_b,bvals);
				J[2] += (params[2]>0?1.0:-1.0)*abs(params[0])*(1-*sumf-temp_f0)*isoterm_b_PVM_multi(dir_iter,_a,_b,bvals);
    			}else{
				sig = abs(params[0])*((1-*sumf)*isoterm_PVM_multi(dir_iter,_a,_b,bvals)+sig);
	    			J[1] += (params[1]>0?1.0:-1.0)*abs(params[0])*(1-*sumf)*isoterm_a_PVM_multi(dir_iter,_a,_b,bvals);
	    			J[2] += (params[2]>0?1.0:-1.0)*abs(params[0])*(1-*sumf)*isoterm_b_PVM_multi(dir_iter,_a,_b,bvals);	
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
extern "C" __global__ void fit_PVM_multi_kernel(	//INPUT
							const double* 		data, 
							const double* 		params_PVM_single_c,
							const double* 		bvecs, 
							const double* 		bvals, 
							const int 		nvox, 
							const int		ndirections,
							const int 		nfib, 	
							const int		nparams,			
							const bool 		m_include_f0,
							const bool		gradnonlin,
							//OUTPUT
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
  	double* x = (double*) &fs[nfib];				//nfib*3
	double* _a = (double*) &x[nfib*3];				//1
	double* _b = (double*) &_a[1];					//1
  	double* sumf = (double*) &_b[1];				//1

	double* C = (double*)&sumf[1];					//nparams*nparams;
	double* el =  (double*)&C[nparams*nparams];			//nparams

	int* indx = (int*)&el[nparams];					//nparams
	int* success = (int*) &indx[nparams];				//1
	int* end = (int*) &success[1];					//1   
	////////// DYNAMIC SHARED MEMORY ///////////

	if(idSubVOX==0){
		
		int nparams_single_c = nparams-1;

		myparams[0] = params_PVM_single_c[(idVOX*nparams_single_c)+0];			//pvm1.get_s0();
  		myparams[1] = 1.0;								//start with d=d_std
  		for(int i=0,ii=3;i<nfib;i++,ii+=3){
    			myparams[ii] = f2x_gpu(params_PVM_single_c[(idVOX*nparams_single_c)+ii-1]);
    			myparams[ii+1] = params_PVM_single_c[(idVOX*nparams_single_c)+ii];
    			myparams[ii+2] = params_PVM_single_c[(idVOX*nparams_single_c)+ii+1];
  		}
		myparams[2] = params_PVM_single_c[(idVOX*nparams_single_c)+1] ; 		//pvm1.get_d();
  		if (m_include_f0)
			myparams[nparams-1]=f2x_gpu(params_PVM_single_c[(idVOX*nparams_single_c)+nparams_single_c-1]);
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
  	//do the fit
	levenberg_marquardt_PVM_multi_gpu(&data[idVOX*ndirections],&bvecs[pos_bvecs],&bvals[pos_bvals],ndirections,nfib,nparams,m_include_f0,idSubVOX,step,grad,hess,inverse, pcf,ncf,lambda,cftol,ltol,olambda,success,end,reduction,fs,x,_a,_b,sumf,C,el,indx,myparams);

	__syncthreads();

  	// finalise parameters
	//m_s0-myparams[0] 	m_d-myparams[1] 	m_d_std-myparams[2]		m_f-m_th-m_ph-myparams[3,4,5,6 etc..]   	m_f0-myparams[nparams-1]

	if(idSubVOX==0){  	
		double aux = myparams[1];

  		myparams[1] = abs(aux*myparams[2]);
		myparams[2] = sqrt(double(abs(aux*myparams[2]*myparams[2])));
  		for(int i=3,k=0;k<nfib;i+=3,k++){
    			myparams[i]  = x2f_gpu(myparams[i]);
  		}
  		if (m_include_f0)
    			myparams[nparams-1]=x2f_gpu(myparams[nparams-1]);

		sort_PVM_multi(nfib,myparams);
  		fix_fsum_PVM_multi(m_include_f0,nfib,nparams,myparams);
	}
	__syncthreads();

	if(idSubVOX<nparams){
		params[(idVOX*nparams)+idSubVOX] = myparams[idSubVOX];
	}
}

//in diffmodel.cc
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
								double*			residuals)
{
	int idSubVOX = threadIdx.x;
	int idVOX = blockIdx.x;

	////////// DYNAMIC SHARED MEMORY ///////////
	extern __shared__ double shared[];
	double* myparams = (double*) shared;			//nparams
	double* fs = (double*) &myparams[nparams];		//nfib
  	double* x = (double*) &fs[nfib];			//nfib*3
	double* _a = (double*) &x[nfib*3];			//1
	double* _b = (double*) &_a[1];				//1
  	double* sumf = (double*) &_b[1];			//1
	int* my_include_f0 = (int*) &sumf[1];			//1	
	////////// DYNAMIC SHARED MEMORY ///////////

	double val;
	double predicted_signal;
	double mydata;

	if(idSubVOX==0){
		*my_include_f0 = includes_f0[idVOX];

  		//m_s0-myparams[0]  m_d-myparams[1]  m_d_std-myparams[2]  m_f-m_th-m_ph-myparams[3,4,5,6 etc..]  m_f0-myparams[nparams-1]

  		myparams[0] = params[(idVOX*nparams)+0];
		double aux1 = params[(idVOX*nparams)+1];
		double aux2 = params[(idVOX*nparams)+2];
		
  		myparams[1] = aux1*aux1/aux2/aux2;		//m_d*m_d/m_d_std/m_d_std;
  		myparams[2] = aux2*aux2/aux1;			//m_d_std*m_d_std/m_d; // =1/beta
  		
  		if (*my_include_f0)
    			myparams[nparams-1]=f2x_gpu(params[(idVOX*nparams)+nparams-1]);
	}

	if(idSubVOX<nfib){
		int kk = 3+3*idSubVOX;
		double sinth,costh,sinph,cosph;
	
		myparams[kk]   = f2x_gpu(params[(idVOX*nparams)+kk]);
    		myparams[kk+1] = params[(idVOX*nparams)+kk+1];
    		myparams[kk+2] = params[(idVOX*nparams)+kk+2];

		sincos(myparams[kk+1],&sinth,&costh);
		sincos(myparams[kk+2],&sinph,&cosph);		
    		fs[idSubVOX] = x2f_gpu(myparams[kk]);
    		x[idSubVOX*3] = sinth*cosph;
    		x[idSubVOX*3+1] = sinth*sinph;
    		x[idSubVOX*3+2] = costh;
  	}

	__syncthreads(); 

	if(idSubVOX==0){
  		*_a = abs(myparams[1]);
  		*_b = abs(myparams[2]);
  		*sumf=0;
  		for(int k=0;k<nfib;k++){
	    		*sumf += fs[k];
		}
  	}
  	
	int ndir = ndirections/THREADS_BLOCK_FIT;
	if(idSubVOX<(ndirections%THREADS_BLOCK_FIT)) ndir++;
	
	double3 x2;
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
  		predicted_signal=0;	//pred = 0;
    		val = 0.0;
    		for(int k=0;k<nfib;k++){
			x2.x=x[k*3];
			x2.y=x[k*3+1];
			x2.z=x[k*3+2];	 
      			val += fs[k]*anisoterm_PVM_multi(dir_iter,_a,_b,x2,&bvecs[pos_bvecs],&bvals[pos_bvals],ndirections);
    		}	
    		if (*my_include_f0){
      			double temp_f0=x2f_gpu(myparams[nparams-1]);
      			predicted_signal = abs(myparams[0])*(temp_f0+(1-*sumf-temp_f0)*isoterm_PVM_multi(dir_iter,_a,_b,&bvals[pos_bvals])+val);
    		}else{
      			predicted_signal = abs(myparams[0])*((1-*sumf)*isoterm_PVM_multi(dir_iter,_a,_b,&bvals[pos_bvals])+val); 
  		}   

		//residuals=m_data-predicted_signal;
		residuals[idVOX*ndirections+dir_iter]= mydata - predicted_signal;

		dir_iter+=THREADS_BLOCK_FIT;
	}
}
