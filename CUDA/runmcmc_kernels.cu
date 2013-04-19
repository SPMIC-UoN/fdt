/*  runmcmc_kernels.cu

    Tim Behrens, Saad Jbabdi, Stam Sotiropoulos, Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */

#include <iostream>
#include <fstream>
#include <stdio.h>
#include "fibre_gpu.h"
#include <math.h>
#include <string.h>
#include <string>
#include <cuda_runtime.h>
#include <cuda.h>
#include <curand.h>

#include <options.h>

#define maxfloat 1e10

__device__
inline void compute_signal(double *signals,double *oldsignals,double mbvals,float* m_d, float* m_dstd,float angtmp, int model){

	oldsignals[0]=signals[0];

	if(model==1 || m_dstd[0]<1e-5){	
		signals[0]=exp(double(-m_d[0]*mbvals*angtmp));
	}else{
		float dbeta= m_d[0]/(m_dstd[0]*m_dstd[0]);
	 	float dalpha=m_d[0]*dbeta;         
           	signals[0]=exp(double(log(double(dbeta/(dbeta+mbvals*angtmp)))*dalpha)); 
	}
}

__device__ 
inline void compute_iso_signal(double *isosignals,double *oldisosignals, double mbvals,float* m_d, float* m_dstd, int model){

	*oldisosignals=*isosignals;

	if(model==1 || m_dstd[0]<1e-5){
	 	*isosignals=exp(-m_d[0]*mbvals);	
	}else{
		float dbeta= m_d[0]/(m_dstd[0]*m_dstd[0]);
	  	float dalpha=m_d[0]*dbeta;
		*isosignals=exp(double(log(double(dbeta/(dbeta+mbvals)))*dalpha));
	}
}

__device__ 
inline void compute_prior(float *m_prior_en, float *m_prior_en_old,float* m_d_prior,float* m_S0_prior,float *m_prior_enf, float* m_f0_prior, float* m_tau_prior, float* m_dstd_prior, int nfib){			
        *m_prior_en_old=*m_prior_en;
	*m_prior_en=*m_d_prior+*m_S0_prior+*m_dstd_prior+*m_tau_prior+*m_f0_prior;
	for(int f=0;f<nfib;f++){
		*m_prior_en=*m_prior_en+m_prior_enf[f];
	}	
}

__device__ 
inline float logIo(const float& x){
    	float y,b;
    	b= fabs(x);
    	if (b<3.75){
      		float a=x/3.75;
      		a*=a;
      		//Bessel function evaluation
		y=1.0+a*(3.5156229+a*(3.0899424+a*(1.2067492+a*(0.2659732+a*(0.0360768+a*0.0045813)))));
      		y=log(double(y));
    	}else{
      		float a=3.75/b; 
      		//Bessel function evaluation
      		//y=(exp(b)/sqrt(b))*(0.39894228+a*(0.01328592+a*(0.00225319+a*(-0.00157565+a*(0.00916281+a*(-0.02057706+a*(0.02635537+a*(-0.01647633+a*0.00392377))))))));
      		//Logarithm of Bessel function

		y=b+log(double((0.39894228+a*(0.01328592+a*(0.00225319+a*(-0.00157565+a*(0.00916281+a*(-0.02057706+a*(0.02635537+a*(-0.01647633+a*0.00392377))))))))/sqrt(b)));
    	}

    	return y;
}

__device__ 
inline void compute_likelihood(int idSubVOX,float* m_S0,float *m_likelihood_en,float *m_f,double *signals,double *isosignals,double *mdata,float* fsum,double *reduction, float* m_f0, const bool rician, float* m_tau,int mydirs, int threadsBlock, int ndirections, int nfib){
	
	double pred[MAXNDIRS_PER_THREAD];

	for(int i=0; i<mydirs; i++){
		pred[i]=0;
	      	for(int f=0;f<nfib;f++){
			pred[i]= pred[i]+m_f[f]*signals[i*nfib+f];
	     	}
	}

	for(int i=0; i<mydirs; i++){
		pred[i]= m_S0[0]*(pred[i]+(1-fsum[0])*isosignals[i]+m_f0[0]); //F0
	}

	reduction[idSubVOX]=0;
	for(int i=0; i<mydirs; i++){	
		if(!rician){
			double diff = mdata[i]-pred[i];
			reduction[idSubVOX] = reduction[idSubVOX]+(diff*diff);
		}else{
			pred[i]= log(mdata[i])+(-0.5*m_tau[0]*(mdata[i]*mdata[i]+pred[i]*pred[i])+logIo(m_tau[0]*pred[i]*mdata[i]));  
			reduction[idSubVOX] = reduction[idSubVOX]+pred[i];
		}
	}

	__syncthreads();

	unsigned int s2=threadsBlock;
	for(unsigned int s=threadsBlock>>1; s>0; s>>=1) {
		if((s2%2)&&(idSubVOX==(s-1))) reduction[idSubVOX]= reduction[idSubVOX] + reduction[idSubVOX + s +1]; 
        	if (idSubVOX < s){
            		reduction[idSubVOX] = reduction[idSubVOX] + reduction[idSubVOX + s];
       	 	}
		s2=s;
        	__syncthreads();
    	}
	if(idSubVOX==0){
		double sumsquares=0;
		sumsquares+=reduction[0];
		if(!rician){ 
		 	*m_likelihood_en=(ndirections/2.0)*log(sumsquares/2.0); 
		}else{
			*m_likelihood_en= -ndirections*log(m_tau[0])-sumsquares;
		}
	}
}
			  
extern "C" __global__ void init_Fibres_Multifibres_kernel(	//INPUT
								const double*			datam,
								const double*			params,
								const float*			tau,
								const double*			bvals,
								const double*			alpha,
								const double*			beta,
								const int			ndirections,
								const int 			nfib,
								const int 			nparams_fit,
								const int 			model,
								const float 			fudgevalue,
								const bool			m_includef0,
								const bool			rician,
								const bool 			m_ardf0,	// opts.ardf0.value()
								const bool 			ard_value,	// opts.all_ard.value()
								const bool 			no_ard_value,	// opts.no_ard.value()
								const bool			gradnonlin,
								//OUTPUT
								FibreGPU*			fibres,
								MultifibreGPU*			multifibres,
								double*				signals,
								double*				isosignals)
{
	int idSubVOX= threadIdx.x;
	int threadsBlock = blockDim.x;
	int idVOX= blockIdx.x;
	bool leader = (idSubVOX==0);

	////////// DYNAMIC SHARED MEMORY ///////////
	extern __shared__ double shared[];
	double* reduction = (double*)shared;				//threadsBlock
	double* tmp = (double*) &reduction[threadsBlock];		//1

	float* m_S0 = (float*) &tmp[1];					//1
	float* m_d = (float*) &m_S0[1];					//1
	float* m_dstd =(float*) &m_d[1];				//1
	float* m_f0 = (float*) &m_dstd[1];				//1
	float* m_tau = (float*) &m_f0[1];				//1
	float* m_th = (float*) &m_tau[1];				//nfib
	float* m_ph = (float*) &m_th[nfib];				//nfib
	float* m_f = (float*) &m_ph[nfib];				//nfib

	float* fsum = (float*) &m_f[nfib];				//1
	float* m_likelihood_en = (float*) &fsum[1];			//1
	float* m_prior_en = (float*) &m_likelihood_en[1];		//1
	////////// DYNAMIC SHARED MEMORY ///////////
	

	double mysig[MAXNDIRS_PER_THREAD*MAXNFIBRES];				
	double myisosig[MAXNDIRS_PER_THREAD];				
	double mydata[MAXNDIRS_PER_THREAD];
	double mybvals[MAXNDIRS_PER_THREAD];	
	double myalpha[MAXNDIRS_PER_THREAD];	
	double mybeta[MAXNDIRS_PER_THREAD];	
 		
	float angtmp[MAXNDIRS_PER_THREAD*MAXNFIBRES];			//this is to pre compute signal
	float cos_alpha_minus_theta;
	float cos_alpha_plus_theta;
	
	// m_s0-params[0] 	m_d-params[1] 	m_f-m_th-m_ph-params[2,3,4,5, etc..]   	m_f0-params[nparams-1]
	if(leader){
		*m_S0 = params[idVOX*nparams_fit];
		multifibres[idVOX].m_S0 = *m_S0;
		multifibres[idVOX].m_S0_prior = 0;
		multifibres[idVOX].m_S0_acc = 0;
		multifibres[idVOX].m_S0_rej = 0;
	
		*m_d=params[idVOX*nparams_fit+1];
		if(*m_d<0 || *m_d>0.008) *m_d=2e-3;			//this is in xfibres...after fit
		multifibres[idVOX].m_d = *m_d;
		multifibres[idVOX].m_d_prior = 0;
		multifibres[idVOX].m_d_acc = 0;
		multifibres[idVOX].m_d_rej = 0;

		if(model==2){ 
			*m_dstd=params[idVOX*nparams_fit+2];
			if(*m_dstd<0 || *m_dstd>0.01) *m_dstd=*m_d/10;	//this is in xfibres...after fit
		}
		else *m_dstd = 0;
		multifibres[idVOX].m_dstd = *m_dstd;
		multifibres[idVOX].m_dstd_prior = 0;
		multifibres[idVOX].m_dstd_acc = 0;
		multifibres[idVOX].m_dstd_rej = 0;

		if (m_includef0) *m_f0=params[idVOX*nparams_fit+nparams_fit-1];
		else *m_f0=0;
		multifibres[idVOX].m_f0 = *m_f0;
		multifibres[idVOX].m_f0_prior = 0;
		multifibres[idVOX].m_f0_acc = 0;
		multifibres[idVOX].m_f0_rej = 0;

		*m_tau = tau[idVOX];
		multifibres[idVOX].m_tau = *m_tau;
		multifibres[idVOX].m_tau_prior = 0;
		multifibres[idVOX].m_tau_acc = 0;
		multifibres[idVOX].m_tau_rej = 0;
	}
	__syncthreads();

	int mydirs = ndirections/threadsBlock;
	int mod = ndirections%threadsBlock;
	if(mod&&(idSubVOX<mod)) mydirs++;

	if(idSubVOX<nfib){
		int add=0;
		if(model==2) add=1;		// if model 2 we have d_std and then 1 more parameter in position 2
		int pos = (idVOX*nfib)+idSubVOX;

		m_th[idSubVOX]=params[idVOX*nparams_fit+2+3*idSubVOX+1+add];
		fibres[pos].m_th = m_th[idSubVOX];
		fibres[pos].m_th_prop = 0.2;
		float m_th_prior = 0;
		fibres[pos].m_th_acc = 0;
		fibres[pos].m_th_rej = 0;
		
		//compute_th_prior();
	      	if(m_th==0){
			m_th_prior=0;
		}else{
			m_th_prior=-log(double(fabs(sin(double(m_th[idSubVOX]))/2)));
	      	}
		fibres[pos].m_th_prior = m_th_prior;
		
		float m_ph_prior=0;	//compute_ph_prior();
		m_ph[idSubVOX]=params[idVOX*nparams_fit+2+3*idSubVOX+2+add];
		fibres[pos].m_ph = m_ph[idSubVOX];
		fibres[pos].m_ph_prop = 0.2;
		fibres[pos].m_ph_prior = 0;	//compute_ph_prior();
		fibres[pos].m_ph_acc = 0;
		fibres[pos].m_ph_rej = 0;

		m_f[idSubVOX] = params[idVOX*nparams_fit+2+3*idSubVOX+add]; 
		fibres[pos].m_f=m_f[idSubVOX];
		fibres[pos].m_f_prop = 0.2;
		float m_f_prior = 0;
		fibres[pos].m_f_acc = 0;
		fibres[pos].m_f_rej = 0;
			
		if(idSubVOX==0){
			fibres[pos].m_lam_jump = ard_value;
		}else{
			fibres[pos].m_lam_jump = !no_ard_value;
		}

		//compute_f_prior();
      		if (m_f[idSubVOX]<=0 | m_f[idSubVOX]>=1 ){
      		}else{
	  		if(fibres[pos].m_lam_jump){              
	    			m_f_prior=log(double(m_f[idSubVOX]));
	  		}else{
	    			m_f_prior=0;
			}
			m_f_prior= fudgevalue* m_f_prior;
      		}
		fibres[pos].m_f_prior = m_f_prior;

		//fibres[vox].m_lam = m_lam; ??
		//fibres[vox].m_lam_prop = 1;
		//fibres[vox].m_lam_prior = 0;
		//compute_lam_prior();

		//compute_prior();
		fibres[pos].m_prior_en= m_th_prior + m_ph_prior + m_f_prior;
	}

	for(int i=0; i<mydirs; i++){	
		int pos = (idVOX*ndirections)+idSubVOX+i*threadsBlock;
		mydata[i] = datam[pos];
	}
	if(gradnonlin){
		for(int i=0; i<mydirs; i++){	
			int pos = (idVOX*ndirections)+idSubVOX+i*threadsBlock;
			mybvals[i] = bvals[pos];
			myalpha[i] = alpha[pos];
			mybeta[i] = beta[pos];
		}
	}else{
		for(int i=0; i<mydirs; i++){	
			int pos = idSubVOX+i*threadsBlock;
			mybvals[i] = bvals[pos];
			myalpha[i] = alpha[pos];
			mybeta[i] = beta[pos];
		}

	}

	__syncthreads();

	//compute_signal_pre
	for(int i=0; i<mydirs; i++){
		for(int f=0;f<nfib;f++){	
			cos_alpha_minus_theta=cos(double(myalpha[i]-m_th[f]));   
			cos_alpha_plus_theta=cos(double(myalpha[i]+m_th[f]));
		     	angtmp[i*nfib+f]= (cos(double(m_ph[f]-mybeta[i]))*(cos_alpha_minus_theta-cos_alpha_plus_theta)/2)+(cos_alpha_minus_theta+cos_alpha_plus_theta)/2;
		 	angtmp[i*nfib+f]=angtmp[i*nfib+f]*angtmp[i*nfib+f];
		}
	}

	//compute_signal()
	double old;
	for(int i=0; i<mydirs; i++){
		for(int f=0;f<nfib;f++){
			compute_signal(&mysig[i*nfib+f],&old,mybvals[i],m_d,m_dstd,angtmp[i*nfib+f],model);
		}
	}

	if(leader){
		//for likehood 
		*fsum=*m_f0;
		for(int g=0;g<nfib;g++){  
			*fsum+=m_f[g];
		}
		//initialise_energies();
	      	//compute_d_prior(); m_d_prior=0; so i don't do nothing, it is already 0
	      	if(model==2){
			//compute_d_std_prior();
			if(*m_dstd<=0 || *m_dstd>0.01){
	      		}else{
				multifibres[idVOX].m_dstd_prior=log(double(*m_dstd));
			}
		}
	      	//compute_tau_prior(); m_tau_prior=0; so it doesn't do nothing, it is already 0
	      	if (m_includef0){
			//compute_f0_prior();
			if (*m_f0<=0 || *m_f0>=1){
	      		}else{
				if(!m_ardf0){}     	//Without ARD
				else              	//With ARD
		  			multifibres[idVOX].m_f0_prior= log(double(*m_f0));
	      		}
		}
	      	//compute_S0_prior(); m_S0_prior=0; so i don't do nothing, it is already 0
		*m_prior_en = 0;
	      	//compute_prior();
	      	//m_prior_en=m_d_prior+m_S0_prior; is 0
	      	if(model==2)
			*m_prior_en= *m_prior_en+multifibres[idVOX].m_dstd_prior;
	      	//if(m_rician) m_prior_en=m_prior_en+m_tau_prior; is 0
	      	if (m_includef0)
			*m_prior_en=*m_prior_en+multifibres[idVOX].m_f0_prior;
	      	for(int fib=0;fib<nfib;fib++){
			*m_prior_en=*m_prior_en+ fibres[idVOX*nfib+fib].m_prior_en;
	      	} 
		multifibres[idVOX].m_prior_en = *m_prior_en;
	}

	//compute_iso_signal()
	for(int i=0; i<mydirs; i++){
		compute_iso_signal(&myisosig[i],&old, mybvals[i],m_d,m_dstd,model);
	}		
 
	__syncthreads();

	//compute_likelihood()
	compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,mysig,myisosig,mydata,fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);

	__syncthreads();

	if(leader){
		multifibres[idVOX].m_likelihood_en = *m_likelihood_en;
	      	
		//compute_energy();	
		multifibres[idVOX].m_energy = *m_prior_en+*m_likelihood_en;

	    	//initialise_props();
	      	multifibres[idVOX].m_S0_prop=multifibres[idVOX].m_S0/10.0; 
	      	multifibres[idVOX].m_d_prop=*m_d/10.0;
	      	multifibres[idVOX].m_dstd_prop=*m_dstd/10.0;
	      	multifibres[idVOX].m_tau_prop=*m_tau/2.0;
	      	multifibres[idVOX].m_f0_prop=0.2;
	}

	for(int i=0; i<mydirs; i++){	
		isosignals[(idVOX*ndirections)+idSubVOX+i*threadsBlock] = myisosig[i];
		for(int f=0;f<nfib;f++){
			signals[(idVOX*ndirections*nfib)+(f*ndirections)+idSubVOX+i*threadsBlock]=mysig[i*nfib+f];
		}
	}
}


extern "C" __global__ void runmcmc_burnin_kernel(	//INPUT 
							const double*			datam,
							const double*			bvals,
							const double*			alpha,
							const double*			beta,
							float*				randomsN,
							float*				randomsU,
							const int			ndirections,
							const int			nfib,
							const int			nparams,
							const int 			model,
							const float 			fudgevalue,
							const bool 			m_include_f0,
							const bool 			m_ardf0,
							const bool 			can_use_ard, 
							const bool 			rician,
							const bool			gradnonlin,
							const int 			updateproposalevery, 	//update every this number of iterations	
							const int 			iterations,		//num of iterations to do this time (maybe is a part of the total)
							const int 			current_iter,		//the number of the current iteration over the total iterations
							//INPUT-OUTPUT
							FibreGPU*			fibres,
							MultifibreGPU*			multifibres,
							double*				signals,
							double*				isosignals)
{	
	int idSubVOX= threadIdx.x;
	int threadsBlock = blockDim.x;
	bool leader = (idSubVOX==0);

	////////// DYNAMIC SHARED MEMORY ///////////
	extern __shared__ double shared[];
	double* reduction = (double*)shared;				//threadsBlock
	double* tmp = (double*) &reduction[threadsBlock];		//1

	float* m_S0 = (float*) &tmp[1];					//1
	float* m_d = (float*) &m_S0[1];					//1
	float* m_dstd =(float*) &m_d[1];				//1
	float* m_f0 = (float*) &m_dstd[1];				//1
	float* m_tau = (float*) &m_f0[1];				//1
	float* m_th = (float*) &m_tau[1];				//nfib
	float* m_ph = (float*) &m_th[nfib];				//nfib
	float* m_f = (float*) &m_ph[nfib];				//nfib

	float* m_S0_prior = (float*) &m_f[nfib];			//1
	float* m_d_prior = (float*) &m_S0_prior[1];			//1
	float* m_dstd_prior = (float*) &m_d_prior[1];			//1
	float* m_f0_prior = (float*) &m_dstd_prior[1];			//1
	float* m_tau_prior = (float*) &m_f0_prior[1];			//1
	float* m_th_prior = (float*) &m_tau_prior[1];			//nfib
	float* m_ph_prior = (float*) &m_th_prior[nfib];			//nfib
	float* m_f_prior = (float*) &m_ph_prior[nfib];			//nfib

	float* m_S0_prop = (float*) &m_f_prior[nfib];			//1
	float* m_d_prop = (float*) &m_S0_prop[1];			//1
	float* m_dstd_prop = (float*) &m_d_prop[1];			//1
	float* m_f0_prop = (float*) &m_dstd_prop[1];			//1
	float* m_tau_prop = (float*) &m_f0_prop[1];			//1
	float* m_th_prop = (float*) &m_tau_prop[1];			//nfib
	float* m_ph_prop = (float*) &m_th_prop[nfib];			//nfib
	float* m_f_prop = (float*) &m_ph_prop[nfib];			//nfib

	float* fsum = (float*) &m_f_prop[nfib];				//1
	float* m_likelihood_en = (float*) &fsum[1];			//1
	float* m_prior_en = (float*) &m_likelihood_en[1];		//1
	float* m_old_prior_en = (float*) &m_prior_en[1];		//1
	float* fm_prior_en = (float*) &m_old_prior_en[1];		//nfib		
	float* fm_old_prior_en = (float*) &fm_prior_en[nfib];		//1		
	float* m_energy = (float*) &fm_old_prior_en[1];			//1
	float* m_old_energy = (float*) &m_energy[1];			//1
	float* old = (float*) &m_old_energy[1];				//2

	float* prerandN = (float*) &old[2];				//nparams
	float* prerandU = (float*) &prerandN[nparams];			//nparams

	int* m_S0_acc = (int*) &prerandU[nparams];			//1
	int* m_d_acc = (int*) &m_S0_acc[1];				//1
	int* m_dstd_acc = (int*) &m_d_acc[1];				//1
	int* m_f0_acc = (int*) &m_dstd_acc[1];				//1
	int* m_tau_acc = (int*) &m_f0_acc[1];				//1
	int* m_th_acc = (int*) &m_tau_acc[1];				//nfib
	int* m_ph_acc = (int*) &m_th_acc[nfib];				//nfib
	int* m_f_acc = (int*) &m_ph_acc[nfib];				//nfib

	int* m_S0_rej = (int*) &m_f_acc[nfib];				//1
	int* m_d_rej = (int*) &m_S0_rej[1];				//1
	int* m_dstd_rej = (int*) &m_d_rej[1];				//1	
	int* m_f0_rej = (int*) &m_dstd_rej[1];				//1
	int* m_tau_rej = (int*) &m_f0_rej[1];				//1	
	int* m_th_rej = (int*) &m_tau_rej[1];				//nfib
	int* m_ph_rej = (int*) &m_th_rej[nfib];				//nfib
	int* m_f_rej = (int*) &m_ph_rej[nfib];				//nfib
	
	int* rejflag = (int*) &m_f_rej[nfib];				//3					
	int* m_lam_jump = (int*) &rejflag[3];				//nfib	
	int* idVOX = (int*) &m_lam_jump[nfib];				//1		
	int* count_update = (int*) &idVOX[1];				//1		
	int* localrand = (int*) &count_update[1];			//1
	////////// DYNAMIC SHARED MEMORY ///////////

	double mysig[MAXNDIRS_PER_THREAD*MAXNFIBRES];			
	double mysigold[MAXNDIRS_PER_THREAD*MAXNFIBRES];			
	double myisosig[MAXNDIRS_PER_THREAD];			
	double myisosigold[MAXNDIRS_PER_THREAD];				
	double mydata[MAXNDIRS_PER_THREAD];
	double mybvals[MAXNDIRS_PER_THREAD];	
	double myalpha[MAXNDIRS_PER_THREAD];	
	double mybeta[MAXNDIRS_PER_THREAD];	
 		
	float angtmp[MAXNDIRS_PER_THREAD*MAXNFIBRES];			//this is to pre compute signal
	float old_angtmp[MAXNDIRS_PER_THREAD*MAXNFIBRES];
	float cos_alpha_minus_theta;
	float cos_alpha_plus_theta;
	
	if (leader){
		*idVOX= blockIdx.x;
		*count_update = current_iter;

		*m_prior_en=multifibres[*idVOX].m_prior_en;

		if(model==2){
			*m_dstd_acc=multifibres[*idVOX].m_dstd_acc;
			*m_dstd_rej=multifibres[*idVOX].m_dstd_rej;
			*m_dstd_prior=multifibres[*idVOX].m_dstd_prior;
			*m_dstd_prop=multifibres[*idVOX].m_dstd_prop;
			*m_dstd=multifibres[*idVOX].m_dstd;
		}else{
			*m_dstd_acc=0;
			*m_dstd_rej=0;
			*m_dstd_prior=0;
			*m_dstd_prop=0;
			*m_dstd=0;
		}
	
		*m_d=multifibres[*idVOX].m_d;
		*m_energy=multifibres[*idVOX].m_energy;
		*m_d_prop=multifibres[*idVOX].m_d_prop;
		*m_d_prior=multifibres[*idVOX].m_d_prior;
		*m_S0_prior=multifibres[*idVOX].m_S0_prior;
		*m_S0=multifibres[*idVOX].m_S0;
		*m_likelihood_en=multifibres[*idVOX].m_likelihood_en;
		*m_d_acc=multifibres[*idVOX].m_d_acc;
		*m_d_rej=multifibres[*idVOX].m_d_rej;
		*m_S0_acc=multifibres[*idVOX].m_S0_acc;
		*m_S0_rej=multifibres[*idVOX].m_S0_rej;
		*m_S0_prop=multifibres[*idVOX].m_S0_prop;

		if(m_include_f0){
			*m_f0_acc=multifibres[*idVOX].m_f0_acc;
			*m_f0_rej=multifibres[*idVOX].m_f0_rej;
			*m_f0_prop=multifibres[*idVOX].m_f0_prop;
			*m_f0_prior=multifibres[*idVOX].m_f0_prior;
			*m_f0=multifibres[*idVOX].m_f0;
		}else{ 
			*m_f0_acc=0;
			*m_f0_rej=0;
			*m_f0_prop=0;
			*m_f0_prior=0;
			*m_f0=0;
		}
				
		if(rician){
			*m_tau_acc=multifibres[*idVOX].m_tau_acc;
			*m_tau_rej=multifibres[*idVOX].m_tau_rej;
			*m_tau_prop=multifibres[*idVOX].m_tau_prop;
			*m_tau_prior=multifibres[*idVOX].m_tau_prior;
			*m_tau=multifibres[*idVOX].m_tau;	
		}else{ 
			*m_tau_acc=0;
			*m_tau_rej=0;
			*m_tau_prop=0;
			*m_tau_prior=0;
			*m_tau=0;
		}
	}
	__syncthreads();

	int mydirs = ndirections/threadsBlock;
	int mod = ndirections%threadsBlock;
	if(mod&&(idSubVOX<mod)) mydirs++;

	if(idSubVOX<nfib){
		int pos = (*idVOX*nfib)+idSubVOX;
		m_th[idSubVOX]=fibres[pos].m_th;
		m_ph[idSubVOX]=fibres[pos].m_ph;
		m_f[idSubVOX]=fibres[pos].m_f;
	
		m_th_acc[idSubVOX]=fibres[pos].m_th_acc;
		m_th_rej[idSubVOX]=fibres[pos].m_th_rej;
		m_ph_acc[idSubVOX]=fibres[pos].m_ph_acc;
		m_ph_rej[idSubVOX]=fibres[pos].m_ph_rej;
		m_f_acc[idSubVOX]=fibres[pos].m_f_acc;
		m_f_rej[idSubVOX]=fibres[pos].m_f_rej;

		fm_prior_en[idSubVOX]=fibres[pos].m_prior_en;
		m_th_prior[idSubVOX]=fibres[pos].m_th_prior;
		m_ph_prior[idSubVOX]=fibres[pos].m_ph_prior;
		m_f_prior[idSubVOX]=fibres[pos].m_f_prior;

		m_th_prop[idSubVOX]=fibres[pos].m_th_prop;
		m_ph_prop[idSubVOX]=fibres[pos].m_ph_prop;
		m_f_prop[idSubVOX]=fibres[pos].m_f_prop;

		m_lam_jump[idSubVOX]=fibres[pos].m_lam_jump;		
	}

	__syncthreads();

	for(int i=0; i<mydirs; i++){	
		int pos = (*idVOX*ndirections)+idSubVOX+i*threadsBlock;
		myisosig[i] = isosignals[pos];
		mydata[i] = datam[pos];
	}
	if(gradnonlin){
		for(int i=0; i<mydirs; i++){	
			int pos = (*idVOX*ndirections)+idSubVOX+i*threadsBlock;
			mybvals[i] = bvals[pos];
			myalpha[i] = alpha[pos];
			mybeta[i] = beta[pos];
		}
	}else{
		for(int i=0; i<mydirs; i++){	
			int pos = idSubVOX+i*threadsBlock;
			mybvals[i] = bvals[pos];
			myalpha[i] = alpha[pos];
			mybeta[i] = beta[pos];
		}

	}
	
	for(int f=0;f<nfib;f++){
		for(int i=0; i<mydirs; i++){	
			mysig[i*nfib+f]= signals[(*idVOX*ndirections*nfib)+(f*ndirections)+idSubVOX+i*threadsBlock];	
		}
	}

	__syncthreads();
		
	//compute_signal_pre
	for(int i=0; i<mydirs; i++){
		for(int f=0;f<nfib;f++){	
			cos_alpha_minus_theta=cos(double(myalpha[i]-m_th[f]));   
			cos_alpha_plus_theta=cos(double(myalpha[i]+m_th[f]));
		     	angtmp[i*nfib+f]= (cos(double(m_ph[f]-mybeta[i]))*(cos_alpha_minus_theta-cos_alpha_plus_theta)/2)+(cos_alpha_minus_theta+cos_alpha_plus_theta)/2;
		 	angtmp[i*nfib+f]=angtmp[i*nfib+f]*angtmp[i*nfib+f];
		}
	}

	for (int niter=0; niter<iterations; niter++){
		//code jump()

		//prefetching randoms
		if (leader){
			*count_update=*count_update+1;
			int idrand = *idVOX*iterations*nparams+(niter*nparams);
			*localrand = 0;
			if(m_include_f0){
				prerandN[*localrand]=randomsN[idrand];
				prerandU[*localrand]=randomsU[idrand];
				idrand++;
				*localrand=*localrand+1;
			}
			if(rician){
				prerandN[*localrand]=randomsN[idrand];
				prerandU[*localrand]=randomsU[idrand];
				idrand++;
				*localrand=*localrand+1;
			}
			prerandN[*localrand]=randomsN[idrand];
			prerandU[*localrand]=randomsU[idrand];
			idrand++;
			*localrand=*localrand+1;

			if(model==2){
				prerandN[*localrand]=randomsN[idrand];
				prerandU[*localrand]=randomsU[idrand];
				idrand++;
				*localrand=*localrand+1;
			}

			prerandN[*localrand]=randomsN[idrand];
			prerandU[*localrand]=randomsU[idrand];
			idrand++;
			*localrand=*localrand+1;

			for(int f=0;f<nfib;f++){
				prerandN[*localrand]=randomsN[idrand];
				prerandU[*localrand]=randomsU[idrand];
				idrand++;
				*localrand=*localrand+1;

				prerandN[*localrand]=randomsN[idrand];
				prerandU[*localrand]=randomsU[idrand];
				idrand++;
				*localrand=*localrand+1;

				prerandN[*localrand]=randomsN[idrand];					
				prerandU[*localrand]=randomsU[idrand];	
				idrand++;	
				*localrand=*localrand+1;	
			}
			*localrand = 0;
		}

////////////////////////////////////////////////////////////////// F0

		if(m_include_f0){
			if (leader){
				//propose_f0
				old[0]=*m_f0;
				*m_f0+=prerandN[*localrand]**m_f0_prop;
				
				//compute_f0_prior()     
				old[1]=*m_f0_prior;
	      			if(*m_f0<=0 || *m_f0 >=1){ 
					rejflag[0]=true;
				}else{ 	
					rejflag[0]=false;
					if(!m_ardf0){
						*m_f0_prior=0;
	      				}else{
						*m_f0_prior=log(double(*m_f0));
					}
				}
				
				//for likehood and reject_f_sum
				*fsum=*m_f0;

				for(int g=0;g<nfib;g++){  
					*fsum+=m_f[g];
				}

				//compute_prior()
				compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,nfib);
				
				//reject_f_sum()
				rejflag[1]=(*fsum>1);
			}

			__syncthreads();

			//compute_likelihood()
			compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,mysig,myisosig,mydata,fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);

			__syncthreads();

			if (leader){
				//compute_energy()
				*m_old_energy=*m_energy;
	      			*m_energy=*m_prior_en+*m_likelihood_en;

				//test_energy
				*tmp=exp(double(*m_old_energy-*m_energy));
				rejflag[2]=(*tmp>prerandU[*localrand]);
				*localrand=*localrand+1;	

				if(!rejflag[0]){
					if(!rejflag[1]){
						if(rejflag[2]){
							//accept_f0()
		  					*m_f0_acc=*m_f0_acc+1;   
						}else{
							//reject_F0()
							*m_f0=old[0];
							*m_f0_prior=old[1];
		      					*m_prior_en=*m_old_prior_en;
		      					*m_f0_rej=*m_f0_rej+1;

							//restore_energy()
      							*m_energy=*m_old_energy;
						}
					}else{ 
						//reject_F0()
						*m_f0=old[0];
						*m_f0_prior=old[1];
		      				*m_prior_en=*m_old_prior_en;
		      				*m_f0_rej=*m_f0_rej+1;

						//restore_energy()
		      				*m_energy=*m_old_energy;
					}
				}else{ 
					//reject_F0()
					*m_f0=old[0];
					*m_f0_prior=old[1];
		      			*m_prior_en=*m_old_prior_en;
		      			*m_f0_rej=*m_f0_rej+1;

					//restore_energy()
		      			*m_energy=*m_old_energy;
				}
			}
		}
			
////////////////////////////////////////////////////////////////// TAU
		if(rician){
			if (leader){
				//propose_tau
				old[0]=*m_tau;
				*m_tau+=prerandN[*localrand]**m_tau_prop;
			
				//compute_tau_prior()     
				old[1]=*m_tau_prior;
	      			if(*m_tau<=0){ 
					rejflag[0]=true;
				}else{ 	
					rejflag[0]=false;
					*m_tau_prior=0;
				}

				//for likehood and reject_f_sum
				*fsum=*m_f0;
				for(int g=0;g<nfib;g++){  
					*fsum+=m_f[g];
				}

				//compute_prior()
				compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,nfib);
			}

			__syncthreads();

			//compute_likelihood()
			compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,mysig,myisosig,mydata,fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);

			__syncthreads();

			if (leader){
				//compute_energy()
				*m_old_energy=*m_energy;
	      			*m_energy=*m_prior_en+*m_likelihood_en;

				//test_energy
				*tmp=exp(double(*m_old_energy-*m_energy));
				rejflag[1]=(*tmp>prerandU[*localrand]);	
				*localrand=*localrand+1;
				
				if(!rejflag[0]){
					if(rejflag[1]){
						//accept_tau()
		  				*m_tau_acc=*m_tau_acc+1;   
					}else{ 
						//reject_tau()
						*m_tau=old[0];
						*m_tau_prior=old[1];
		      				*m_prior_en=*m_old_prior_en;
		      				*m_tau_rej=*m_tau_rej+1;
						//restore_energy()
		      				*m_energy=*m_old_energy;
					}
				}else{ 
					//reject_tau()
					*m_tau=old[0];
					*m_tau_prior=old[1];
		      			*m_prior_en=*m_old_prior_en;
		      			*m_tau_rej=*m_tau_rej+1;

					//restore_energy()
		      			*m_energy=*m_old_energy;
				}
			}	
		}

////////////////////////////////////////////////////////////////// D

		if (leader){
			//propose_d()
			old[0]=*m_d;	
			*m_d+=prerandN[*localrand]**m_d_prop;
				
			//compute_d_prior_f0()      
			old[1]=*m_d_prior;
					
      			if(*m_d<0 || *m_d > 0.008){
				rejflag[0]=true;
			}else{ 
				*m_d_prior=0;
				rejflag[0]=false;
      			}
		}

		__syncthreads();
				
		//compute_signal()
		for(int i=0; i<mydirs; i++){
			for(int f=0;f<nfib;f++){
				compute_signal(&mysig[i*nfib+f],&mysigold[i*nfib+f],mybvals[i],m_d,m_dstd,angtmp[i*nfib+f],model);
			}
		}
		//compute_iso_signal()
		for(int i=0; i<mydirs; i++){
			compute_iso_signal(&myisosig[i],&myisosigold[i], mybvals[i],m_d,m_dstd,model);
		}				

		if (leader){
			//for likehood and reject_f_sum
			*fsum=*m_f0;

			for(int g=0;g<nfib;g++){  
				*fsum+=m_f[g];
			}

			//compute_prior()
			compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,nfib);
		}

		__syncthreads();
				
		//compute_likelihood()
		compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,mysig,myisosig,mydata,fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);	
				
		__syncthreads();
				
		if (leader){
			//compute_energy()
			*m_old_energy=*m_energy;
      			*m_energy=*m_prior_en+*m_likelihood_en;
			//test_energy
			*tmp=exp(double(*m_old_energy-*m_energy));
			rejflag[1]=(*tmp>prerandU[*localrand]);
			*localrand=*localrand+1;
		}
			
		__syncthreads();

       		if(!rejflag[0]){
			if(rejflag[1]){
				//accept_d()
	  			if (leader) *m_d_acc=*m_d_acc+1;   
			}else{
				if (leader){
					//reject_d()
					*m_d=old[0];
					*m_d_prior=old[1];
      					*m_prior_en=*m_old_prior_en;
					*m_d_rej=*m_d_rej+1;
					//restore_energy()
      					*m_energy=*m_old_energy;
				}
      						
				for(int i=0; i<mydirs; i++){
					for(int f=0;f<nfib;f++){
						mysig[i*nfib+f] = mysigold[i*nfib+f];
					}	
				
					myisosig[i]=myisosigold[i];
				}
			}
        	}else{ 
			if (leader){
				//reject_d()
				*m_d=old[0];
				*m_d_prior=old[1];
      				*m_prior_en=*m_old_prior_en;
				*m_d_rej=*m_d_rej+1;
				//restore_energy()
      				*m_energy=*m_old_energy;
			}

      			for(int i=0; i<mydirs; i++){
				for(int f=0;f<nfib;f++){
					mysig[i*nfib+f] = mysigold[i*nfib+f];
				}	
				
				myisosig[i]=myisosigold[i];
			}	
        	}

////////////////////////////////////////////////////////////////// D_STD

		if(model==2){
			if (leader){	
				//propose_d_std
				old[0]=*m_dstd;
				*m_dstd+=prerandN[*localrand]**m_dstd_prop;

				//compute_d_std_prior()     
				old[1]=*m_dstd_prior;
	      			if(*m_dstd<=0 || *m_dstd > 0.01){ 
					rejflag[0]=true;
				}else{ 	
					rejflag[0]=false;
					*m_dstd_prior=log(double(*m_dstd));
				}
			}

			__syncthreads();

			//compute_signal()
			for(int i=0; i<mydirs; i++){
				for(int f=0;f<nfib;f++){
					compute_signal(&mysig[i*nfib+f],&mysigold[i*nfib+f],mybvals[i],m_d,m_dstd,angtmp[i*nfib+f],model);
				}
			}

			//compute_iso_signal()
			for(int i=0; i<mydirs; i++){
				compute_iso_signal(&myisosig[i],&myisosigold[i], mybvals[i],m_d,m_dstd,model);
			}

			if (leader){	
				//for likehood and reject_f_sum
				*fsum=*m_f0;
				for(int g=0;g<nfib;g++){  
					*fsum+=m_f[g];
				}

				//compute_prior()
				compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,nfib);
			}

			__syncthreads();

			//compute_likelihood()
			compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,mysig,myisosig,mydata,fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);

			__syncthreads();

			if (leader){
				//compute_energy()
				*m_old_energy=*m_energy;
	      			*m_energy=*m_prior_en+*m_likelihood_en;

				//test_energy
				*tmp=exp(double(*m_old_energy-*m_energy));
				rejflag[1]=(*tmp>prerandU[*localrand]);
				*localrand=*localrand+1;					
			}

			__syncthreads();
				
			if(!rejflag[0]){
				if(rejflag[1]){
					//accept_dstd()
		  			if (leader) *m_dstd_acc=*m_dstd_acc+1;   
				}else{ 
					if (leader){
						//reject_dstd()
						*m_dstd=old[0];
						*m_dstd_prior=old[1];
		      				*m_prior_en=*m_old_prior_en;
		      				*m_dstd_rej=*m_dstd_rej+1;

						//restore_energy()
		      				*m_energy=*m_old_energy;
					}

					for(int i=0; i<mydirs; i++){
						for(int f=0;f<nfib;f++){
							mysig[i*nfib+f] = mysigold[i*nfib+f];
						}
				
						myisosig[i]=myisosigold[i];
					}
				}
			}else{ 
				if (leader){
					//reject_dstd()
					*m_dstd=old[0];
					*m_dstd_prior=old[1];
		      			*m_prior_en=*m_old_prior_en;
		      			*m_dstd_rej=*m_dstd_rej+1;

					//restore_energy()
		      			*m_energy=*m_old_energy;
				}

				for(int i=0; i<mydirs; i++){
					for(int f=0;f<nfib;f++){
						mysig[i*nfib+f] = mysigold[i*nfib+f];
					}	
				
					myisosig[i]=myisosigold[i];
				}
			}
		}

////////////////////////////////////////////////////////////////// S0

		if (leader){
			//propose_S0()
			old[0]=*m_S0;
			*m_S0+=prerandN[*localrand]**m_S0_prop;
				
			//compute_S0_prior()
			old[1]=*m_S0_prior;
        		if(*m_S0<0) rejflag[0]=true;
        		else{    
				*m_S0_prior=0;
	  			rejflag[0]=false;
        		}
				
			//for likehood and reject_f_sum
			*fsum=*m_f0;
			for(int g=0;g<nfib;g++){  
				*fsum+=m_f[g];
			}

			//compute_prior()
			compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,nfib);
		}

		__syncthreads();

		//compute_likelihood()
		compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,mysig,myisosig,mydata,fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);

		__syncthreads();

		if (leader){
			//compute_energy()
			*m_old_energy=*m_energy;
      			*m_energy=*m_prior_en+*m_likelihood_en;

			//test_energy
			*tmp=exp(double(*m_old_energy-*m_energy));
			rejflag[1]=(*tmp>prerandU[*localrand]);
			*localrand=*localrand+1;

        		if(!rejflag[0]){
				if(rejflag[1]){
					//accept_S0()
	  				*m_S0_acc=*m_S0_acc+1;   
				}else{
					//reject_S0()
					*m_S0=old[0];
					*m_S0_prior=old[1];
      					*m_prior_en=*m_old_prior_en;
      					*m_S0_rej=*m_S0_rej+1;

					//restore_energy()
      					*m_energy=*m_old_energy;
				}
        		}else{ 
				//reject_S0()
				*m_S0=old[0];
				*m_S0_prior=old[1];
	      			*m_prior_en=*m_old_prior_en;
	      			*m_S0_rej=*m_S0_rej+1;

				//restore_energy()
      				*m_energy=*m_old_energy;
			}
        	}

////////////////////////////////////////////////////////////////////////////     TH

     		for(int fibre=0;fibre<nfib;fibre++){  
			if (leader){ 
				//propose_th()
				old[0]=m_th[fibre];
				m_th[fibre]+=prerandN[*localrand]*m_th_prop[fibre];
					
				//compute_th_prior()
				old[1]=m_th_prior[fibre];
      	   			if(m_th[fibre]==0){
					m_th_prior[fibre]=0;
		   		}else{
					m_th_prior[fibre]=-log(double(fabs(sin(double(m_th[fibre]))/2)));
	      	   		}
		  		//rejflag[0]=false; /////////////////always false
		
				//compute_prior()
				*fm_old_prior_en=fm_prior_en[fibre];
	      	   		fm_prior_en[fibre]=m_th_prior[fibre]+m_ph_prior[fibre]+m_f_prior[fibre];	
			}

			__syncthreads();

			//compute_signal()
			//compute_signal_pre	
			for(int i=0; i<mydirs; i++){
				cos_alpha_minus_theta=cos(double(myalpha[i]-m_th[fibre]));   
				cos_alpha_plus_theta=cos(double(myalpha[i]+m_th[fibre]));
				old_angtmp[i*nfib+fibre]=angtmp[i*nfib+fibre];
				angtmp[i*nfib+fibre]= (cos(double(m_ph[fibre]-mybeta[i]))*(cos_alpha_minus_theta-cos_alpha_plus_theta)/2)+(cos_alpha_minus_theta+cos_alpha_plus_theta)/2;
		 	   	angtmp[i*nfib+fibre]=angtmp[i*nfib+fibre]*angtmp[i*nfib+fibre];

				compute_signal(&mysig[i*nfib+fibre],&mysigold[i*nfib+fibre],mybvals[i],m_d,m_dstd,angtmp[i*nfib+fibre],model);
			}

			if (leader){	
				//for likehood and reject_f_sum
				*fsum=*m_f0;
				for(int g=0;g<nfib;g++){  
					*fsum+=m_f[g];
				}

				//compute_prior()
				compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,nfib);	
			}

			__syncthreads();

			//compute_likelihood()
			compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,mysig,myisosig,mydata,fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);

			__syncthreads();

			if (leader){ 
				//compute_energy()
				*m_old_energy=*m_energy;
	      			*m_energy=*m_prior_en+*m_likelihood_en;

				//test_energy
				*tmp=exp(double(*m_old_energy-*m_energy));
				rejflag[1]=(*tmp>prerandU[*localrand]);
				*localrand=*localrand+1;					
			}

			__syncthreads();
			
			if(rejflag[1]){
				//accept_th
		  		if (leader) m_th_acc[fibre]++;   
				}else{

				if (leader){
					//reject_th()
					m_th[fibre]=old[0];
					m_th_prior[fibre]=old[1];
      					fm_prior_en[fibre]=*fm_old_prior_en;
					m_th_rej[fibre]++;						
						
					//restore_energy()
					*m_prior_en=*m_old_prior_en;
	      				*m_energy=*m_old_energy;
				}

				//compute_signal_pre undo
				for(int i=0; i<mydirs; i++){
					angtmp[i*nfib+fibre]=old_angtmp[i*nfib+fibre];
					mysig[i*nfib+fibre] = mysigold[i*nfib+fibre];						
      				}
			}
			__syncthreads();
		
///////////////////////////////////////     PH

			if (leader){
				//propose_ph()
				old[0]=m_ph[fibre];
				m_ph[fibre]+=prerandN[*localrand]*m_ph_prop[fibre];

				//compute_ph_prior()
				old[1]=m_ph_prior[fibre];
      				m_ph_prior[fibre]=0;
      				//rejflag[0]=false;

				//compute_prior()
				*fm_old_prior_en=fm_prior_en[fibre];
      	   			fm_prior_en[fibre]=m_th_prior[fibre]+m_ph_prior[fibre]+m_f_prior[fibre];
			}

			__syncthreads();
				
			//compute_signal()
			//compute_signal_pre
			for(int i=0; i<mydirs; i++){
				cos_alpha_minus_theta=cos(double(myalpha[i]-m_th[fibre]));   
			  	cos_alpha_plus_theta=cos(double(myalpha[i]+m_th[fibre]));
				old_angtmp[i*nfib+fibre]=angtmp[i*nfib+fibre];
		     	   	angtmp[i*nfib+fibre]= (cos(double(m_ph[fibre]-mybeta[i]))*(cos_alpha_minus_theta-cos_alpha_plus_theta)/2)+(cos_alpha_minus_theta+cos_alpha_plus_theta)/2;
		 	   	angtmp[i*nfib+fibre]=angtmp[i*nfib+fibre]*angtmp[i*nfib+fibre];
	
				compute_signal(&mysig[i*nfib+fibre],&mysigold[i*nfib+fibre],mybvals[i],m_d,m_dstd,angtmp[i*nfib+fibre],model);
			}

			if (leader){
				//for likehood and reject_f_sum
				*fsum=*m_f0;
				for(int g=0;g<nfib;g++){  
					*fsum+=m_f[g];
				}

				//compute_prior()
				compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,nfib);
			}

			__syncthreads();

			//compute_likelihood()
			compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,mysig,myisosig,mydata,fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);


			__syncthreads();

			if (leader){
				//compute_energy()
				*m_old_energy=*m_energy;
	      			*m_energy=*m_prior_en+*m_likelihood_en;

				//test_energy
				*tmp=exp(double(*m_old_energy-*m_energy));
				rejflag[1]=(*tmp>prerandU[*localrand]);
				*localrand=*localrand+1;
			}

			__syncthreads();

			//if(!rejflag[0]){
			if(rejflag[1]){
				//accept_ph()
		  		if (leader) m_ph_acc[fibre]++;   
			}else{
				if (leader){
					//reject_ph()
					m_ph[fibre]=old[0];
					m_ph_prior[fibre]=old[1];
      					fm_prior_en[fibre]=*fm_old_prior_en;
					m_ph_rej[fibre]++;						
						
					//restore_energy()
					*m_prior_en=*m_old_prior_en;
	      				*m_energy=*m_old_energy;
				}
				//compute_signal_pre undo
				for(int i=0; i<mydirs; i++){
					angtmp[i*nfib+fibre]=old_angtmp[i*nfib+fibre];

					mysig[i*nfib+fibre] = mysigold[i*nfib+fibre];	
				}
			}

			__syncthreads();

////////////////////////////////////////////             F

			if (leader){
				//propose_f()
				old[0]=m_f[fibre];
				m_f[fibre]+=prerandN[*localrand]*m_f_prop[fibre];

	     			//compute_f_prior()
	        		old[1]=m_f_prior[fibre];
				if (m_f[fibre]<=0 || m_f[fibre]>=1) rejflag[0]=true;
	        		else{
		      			if(!can_use_ard ){
		  				m_f_prior[fibre]=0;
					}else{
		  				if(m_lam_jump[fibre]){
							m_f_prior[fibre]=log(double(m_f[fibre]));
						}else{
		    					m_f_prior[fibre]=0;
		  				}
					}
					m_f_prior[fibre]=fudgevalue*m_f_prior[fibre];
					rejflag[0]=false;
	      			}

				//compute_prior()
				*fm_old_prior_en=fm_prior_en[fibre];
      	   			fm_prior_en[fibre]=m_th_prior[fibre]+m_ph_prior[fibre]+m_f_prior[fibre];
						
				//for likehood and reject_f_sum
				*fsum=*m_f0;
				for(int g=0;g<nfib;g++){  
					*fsum+=m_f[g];
				}
				//reject_f_sum()
				rejflag[1]=(*fsum>1);

				//compute_prior()
				compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,nfib);	
			}

			__syncthreads();

			//compute_likelihood()
			compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,mysig,myisosig,mydata,fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);	

			__syncthreads();

			if (leader){
				//compute_energy()
				*m_old_energy=*m_energy;
		      		*m_energy=*m_prior_en+*m_likelihood_en;

				//test_energy
				*tmp=exp(double(*m_old_energy-*m_energy));
				rejflag[2]=(*tmp>prerandU[*localrand]);
				*localrand=*localrand+1;			

		      		if(!rejflag[0]){
					if(!rejflag[1]){
						if(rejflag[2]){
							//accept_f()
			  				m_f_acc[fibre]++;   
						}else{
							//reject_f()
							m_f[fibre]=old[0];
							m_f_prior[fibre]=old[1];
	      						fm_prior_en[fibre]=*fm_old_prior_en;
	      						m_f_rej[fibre]++;				
					
							//restore_energy()
							*m_prior_en=*m_old_prior_en;
		      					*m_energy=*m_old_energy;
						}
					}else{ 
						//reject_f()
						m_f[fibre]=old[0];
						m_f_prior[fibre]=old[1];
	      					fm_prior_en[fibre]=*fm_old_prior_en;
	      					m_f_rej[fibre]++;

						//restore_energy()
						*m_prior_en=*m_old_prior_en;
		      				*m_energy=*m_old_energy;	
					}
				}else{
					//reject_f()
					m_f[fibre]=old[0];
					m_f_prior[fibre]=old[1];
		      			fm_prior_en[fibre]=*fm_old_prior_en;
		      			m_f_rej[fibre]++;

					//restore_energy()
					*m_prior_en=*m_old_prior_en;
		      			*m_energy=*m_old_energy;
				}
			}
			__syncthreads();	

        	}//end while nfib

        	if(((*count_update%updateproposalevery)==0)&&leader){
			//m_multifibre.update_proposals();
			*m_d_prop*=sqrt(float(*m_d_acc+1)/float(*m_d_rej+1));
			*m_d_prop=min(*m_d_prop,maxfloat);

			if(rician){
				*m_tau_prop*=sqrt(float(*m_tau_acc+1)/float(*m_tau_rej+1));
				*m_tau_prop=min(*m_tau_prop,maxfloat);
				*m_tau_acc=0; 
				*m_tau_rej=0;	
			}

			if(m_include_f0){
				*m_f0_prop*=sqrt(float(*m_f0_acc+1)/float(*m_f0_rej+1));
				*m_f0_prop=min(*m_f0_prop,maxfloat);
				*m_f0_acc=0; 
				*m_f0_rej=0;	
			}	

			if(model==2){
				*m_dstd_prop*=sqrt(float(*m_dstd_acc+1)/float(*m_dstd_rej+1));
				*m_dstd_prop=min(*m_dstd_prop,maxfloat);
				*m_dstd_acc=0; 
				*m_dstd_rej=0;	
			}

			*m_S0_prop*=sqrt(float(*m_S0_acc+1)/float(*m_S0_rej+1));
			*m_S0_prop=min(*m_S0_prop,maxfloat);
			*m_d_acc=0; 
			*m_d_rej=0;
			*m_S0_acc=0; 
			*m_S0_rej=0;
			for(int f=0; f<nfib;f++){
				//m_fibres[f].update_proposals();
				m_th_prop[f]*=sqrt(float(m_th_acc[f]+1)/float(m_th_rej[f]+1));
				m_th_prop[f]=min(m_th_prop[f],maxfloat);
		      		m_ph_prop[f]*=sqrt(float(m_ph_acc[f]+1)/float(m_ph_rej[f]+1));
		      		m_ph_prop[f]=min(m_ph_prop[f],maxfloat);
		      		m_f_prop[f]*=sqrt(float(m_f_acc[f]+1)/float(m_f_rej[f]+1));
		      		m_f_prop[f]=min(m_f_prop[f],maxfloat);
			      
		      		m_th_acc[f]=0; 
		      		m_th_rej[f]=0;
		      		m_ph_acc[f]=0; 
		      		m_ph_rej[f]=0;
		      		m_f_acc[f]=0; 
		      		m_f_rej[f]=0;
			}
		}

		__syncthreads();	

        } //end while iterations

	if(leader){
		multifibres[*idVOX].m_S0=*m_S0;
		multifibres[*idVOX].m_S0_prior=*m_S0_prior;
		multifibres[*idVOX].m_S0_prop=*m_S0_prop;
		multifibres[*idVOX].m_S0_acc=*m_S0_acc;
		multifibres[*idVOX].m_S0_rej=*m_S0_rej;

		multifibres[*idVOX].m_d=*m_d;
		multifibres[*idVOX].m_d_prior=*m_d_prior;
		multifibres[*idVOX].m_d_prop=*m_d_prop;
		multifibres[*idVOX].m_d_acc=*m_d_acc;
		multifibres[*idVOX].m_d_rej=*m_d_rej;
	
		multifibres[*idVOX].m_prior_en=*m_prior_en;
		multifibres[*idVOX].m_energy=*m_energy;
		multifibres[*idVOX].m_likelihood_en=*m_likelihood_en;

		if(m_include_f0){
			multifibres[*idVOX].m_f0_prior=*m_f0_prior;
			multifibres[*idVOX].m_f0=*m_f0;
			multifibres[*idVOX].m_f0_acc=*m_f0_acc;
			multifibres[*idVOX].m_f0_rej=*m_f0_rej;
			multifibres[*idVOX].m_f0_prop=*m_f0_prop;
		}
		if(rician){
			multifibres[*idVOX].m_tau_prior=*m_tau_prior;
			multifibres[*idVOX].m_tau=*m_tau;
			multifibres[*idVOX].m_tau_acc=*m_tau_acc;
			multifibres[*idVOX].m_tau_rej=*m_tau_rej;
			multifibres[*idVOX].m_tau_prop=*m_tau_prop;
		}
		if(model==2){
			multifibres[*idVOX].m_dstd_prior=*m_dstd_prior;
			multifibres[*idVOX].m_dstd=*m_dstd;
			multifibres[*idVOX].m_dstd_acc=*m_dstd_acc;
			multifibres[*idVOX].m_dstd_rej=*m_dstd_rej;
			multifibres[*idVOX].m_dstd_prop=*m_dstd_prop;
		}
	}
	
	if(idSubVOX<nfib){
		int pos = (*idVOX*nfib)+idSubVOX;
	
		fibres[pos].m_th=m_th[idSubVOX];
		fibres[pos].m_ph=m_ph[idSubVOX];
		fibres[pos].m_f=m_f[idSubVOX];

		fibres[pos].m_th_acc=m_th_acc[idSubVOX];
		fibres[pos].m_th_rej=m_th_rej[idSubVOX];
		fibres[pos].m_ph_acc=m_ph_acc[idSubVOX];
		fibres[pos].m_ph_rej=m_ph_rej[idSubVOX];
		fibres[pos].m_f_acc=m_f_acc[idSubVOX];
		fibres[pos].m_f_rej=m_f_rej[idSubVOX];

		fibres[pos].m_prior_en=fm_prior_en[idSubVOX];
		fibres[pos].m_th_prior=m_th_prior[idSubVOX];
		fibres[pos].m_ph_prior=m_ph_prior[idSubVOX];
		fibres[pos].m_f_prior=m_f_prior[idSubVOX];

		fibres[pos].m_th_prop=m_th_prop[idSubVOX];
		fibres[pos].m_ph_prop=m_ph_prop[idSubVOX];
		fibres[pos].m_f_prop=m_f_prop[idSubVOX];

		fibres[pos].m_lam_jump=m_lam_jump[idSubVOX];		
	}

	for(int i=0; i<mydirs; i++){	
		isosignals[(*idVOX*ndirections)+idSubVOX+i*threadsBlock] = myisosig[i];

		for(int f=0;f<nfib;f++){
			signals[(*idVOX*ndirections*nfib)+(f*ndirections)+idSubVOX+i*threadsBlock]=mysig[i*nfib+f];
		}
	}
}

extern "C" __global__ void runmcmc_record_kernel(	//INPUT 
							const double*			datam,
							const double*			bvals,
							const double*			alpha,
							const double*			beta,
							//INPUT-OUTPUT
							FibreGPU*			fibres,
							MultifibreGPU*			multifibres,
							double*				signals,
							double*				isosignals,
							//INPUT 
							float*				randomsN,
							float*				randomsU,
							const int			ndirections,
							const int			nfib,
							const int			nparams,
							const int 			model,
							const float 			fudgevalue,
							const bool 			m_include_f0,
							const bool 			m_ardf0,
							const bool 			can_use_ard, 
							const bool 			rician,
							const bool			gradnonlin,
							const int 			updateproposalevery, 	//update every this number of iterations	
							const int 			iterations,		//num of iterations to do this time (maybe is a part of the total)	
							const int 			current_iter,		//the number of the current iteration over the total iterations
							const int 			iters_burnin,		//iters in burin, we need it to continue the updates at the correct time. 
							const int 			record_every, 		//record every this number
							const int 			totalrecords,		//total number of records to do
							
							//OUTPUT
							float*				rf0,			//records of parameters
							float*				rtau,
							float*				rs0,
							float*				rd,
							float*				rdstd,
							float*				rth,
							float*				rph, 
							float*				rf)

{
	int idSubVOX= threadIdx.x;
	int threadsBlock = blockDim.x;
	bool leader = (idSubVOX==0);

	////////// DYNAMIC SHARED MEMORY ///////////
	extern __shared__ double shared[];
	double* reduction = (double*)shared;				//threadsBlock
	double* tmp = (double*) &reduction[threadsBlock];		//1

	float* m_S0 = (float*) &tmp[1];					//1
	float* m_d = (float*) &m_S0[1];					//1
	float* m_dstd =(float*) &m_d[1];				//1
	float* m_f0 = (float*) &m_dstd[1];				//1
	float* m_tau = (float*) &m_f0[1];				//1
	float* m_th = (float*) &m_tau[1];				//nfib
	float* m_ph = (float*) &m_th[nfib];				//nfib
	float* m_f = (float*) &m_ph[nfib];				//nfib

	float* m_S0_prior = (float*) &m_f[nfib];			//1
	float* m_d_prior = (float*) &m_S0_prior[1];			//1
	float* m_dstd_prior = (float*) &m_d_prior[1];			//1
	float* m_f0_prior = (float*) &m_dstd_prior[1];			//1
	float* m_tau_prior = (float*) &m_f0_prior[1];			//1
	float* m_th_prior = (float*) &m_tau_prior[1];			//nfib
	float* m_ph_prior = (float*) &m_th_prior[nfib];			//nfib
	float* m_f_prior = (float*) &m_ph_prior[nfib];			//nfib

	float* m_S0_prop = (float*) &m_f_prior[nfib];			//1
	float* m_d_prop = (float*) &m_S0_prop[1];			//1
	float* m_dstd_prop = (float*) &m_d_prop[1];			//1
	float* m_f0_prop = (float*) &m_dstd_prop[1];			//1
	float* m_tau_prop = (float*) &m_f0_prop[1];			//1
	float* m_th_prop = (float*) &m_tau_prop[1];			//nfib
	float* m_ph_prop = (float*) &m_th_prop[nfib];			//nfib
	float* m_f_prop = (float*) &m_ph_prop[nfib];			//nfib

	float* fsum = (float*) &m_f_prop[nfib];				//1
	float* m_likelihood_en = (float*) &fsum[1];			//1
	float* m_prior_en = (float*) &m_likelihood_en[1];		//1
	float* m_old_prior_en = (float*) &m_prior_en[1];		//1
	float* fm_prior_en = (float*) &m_old_prior_en[1];		//nfib		
	float* fm_old_prior_en = (float*) &fm_prior_en[nfib];		//1		
	float* m_energy = (float*) &fm_old_prior_en[1];			//1
	float* m_old_energy = (float*) &m_energy[1];			//1
	float* old = (float*) &m_old_energy[1];				//2

	float* prerandN = (float*) &old[2];				//nparams
	float* prerandU = (float*) &prerandN[nparams];			//nparams

	int* m_S0_acc = (int*) &prerandU[nparams];			//1
	int* m_d_acc = (int*) &m_S0_acc[1];				//1
	int* m_dstd_acc = (int*) &m_d_acc[1];				//1
	int* m_f0_acc = (int*) &m_dstd_acc[1];				//1
	int* m_tau_acc = (int*) &m_f0_acc[1];				//1
	int* m_th_acc = (int*) &m_tau_acc[1];				//nfib
	int* m_ph_acc = (int*) &m_th_acc[nfib];				//nfib
	int* m_f_acc = (int*) &m_ph_acc[nfib];				//nfib

	int* m_S0_rej = (int*) &m_f_acc[nfib];				//1
	int* m_d_rej = (int*) &m_S0_rej[1];				//1
	int* m_dstd_rej = (int*) &m_d_rej[1];				//1	
	int* m_f0_rej = (int*) &m_dstd_rej[1];				//1
	int* m_tau_rej = (int*) &m_f0_rej[1];				//1	
	int* m_th_rej = (int*) &m_tau_rej[1];				//nfib
	int* m_ph_rej = (int*) &m_th_rej[nfib];				//nfib
	int* m_f_rej = (int*) &m_ph_rej[nfib];				//nfib
	
	int* rejflag = (int*) &m_f_rej[nfib];				//3					
	int* m_lam_jump = (int*) &rejflag[3];				//nfib	
	int* idVOX = (int*) &m_lam_jump[nfib];				//1		
	int* count_update = (int*) &idVOX[1];				//1	
	int* recordcount = (int*) &count_update[1];			//1
	int* sample = (int*) &recordcount[1];				//1
	int* localrand = (int*) &sample[1];				//1
	////////// DYNAMIC SHARED MEMORY ///////////

	double mysig[MAXNDIRS_PER_THREAD*MAXNFIBRES];			
	double mysigold[MAXNDIRS_PER_THREAD*MAXNFIBRES];			
	double myisosig[MAXNDIRS_PER_THREAD];			
	double myisosigold[MAXNDIRS_PER_THREAD];				
	double mydata[MAXNDIRS_PER_THREAD];
	double mybvals[MAXNDIRS_PER_THREAD];	
	double myalpha[MAXNDIRS_PER_THREAD];	
	double mybeta[MAXNDIRS_PER_THREAD];	
 		
	float angtmp[MAXNDIRS_PER_THREAD*MAXNFIBRES];			//this is to pre compute signal
	float old_angtmp[MAXNDIRS_PER_THREAD*MAXNFIBRES];
	float cos_alpha_minus_theta;
	float cos_alpha_plus_theta;

	if (leader){
		*idVOX= blockIdx.x;
		*count_update = current_iter + iters_burnin;	//count for updates
		*recordcount= current_iter;	
		*sample=1+(current_iter/record_every);		//the next number of sample.....the index start in 0

		*m_prior_en=multifibres[*idVOX].m_prior_en;

		if(model==2){
			*m_dstd_acc=multifibres[*idVOX].m_dstd_acc;
			*m_dstd_rej=multifibres[*idVOX].m_dstd_rej;
			*m_dstd_prior=multifibres[*idVOX].m_dstd_prior;
			*m_dstd_prop=multifibres[*idVOX].m_dstd_prop;
			*m_dstd=multifibres[*idVOX].m_dstd;
		}else{
			*m_dstd_acc=0;
			*m_dstd_rej=0;
			*m_dstd_prior=0;
			*m_dstd_prop=0;
			*m_dstd=0;
		}
	
		*m_d=multifibres[*idVOX].m_d;
		*m_energy=multifibres[*idVOX].m_energy;
		*m_d_prop=multifibres[*idVOX].m_d_prop;
		*m_d_prior=multifibres[*idVOX].m_d_prior;
		*m_S0_prior=multifibres[*idVOX].m_S0_prior;
		*m_S0=multifibres[*idVOX].m_S0;
		*m_likelihood_en=multifibres[*idVOX].m_likelihood_en;
		*m_d_acc=multifibres[*idVOX].m_d_acc;
		*m_d_rej=multifibres[*idVOX].m_d_rej;
		*m_S0_acc=multifibres[*idVOX].m_S0_acc;
		*m_S0_rej=multifibres[*idVOX].m_S0_rej;
		*m_S0_prop=multifibres[*idVOX].m_S0_prop;

		if(m_include_f0){
			*m_f0_acc=multifibres[*idVOX].m_f0_acc;
			*m_f0_rej=multifibres[*idVOX].m_f0_rej;
			*m_f0_prop=multifibres[*idVOX].m_f0_prop;
			*m_f0_prior=multifibres[*idVOX].m_f0_prior;
			*m_f0=multifibres[*idVOX].m_f0;
		}else{ 
			*m_f0_acc=0;
			*m_f0_rej=0;
			*m_f0_prop=0;
			*m_f0_prior=0;
			*m_f0=0;
		}
				
		if(rician){
			*m_tau_acc=multifibres[*idVOX].m_tau_acc;
			*m_tau_rej=multifibres[*idVOX].m_tau_rej;
			*m_tau_prop=multifibres[*idVOX].m_tau_prop;
			*m_tau_prior=multifibres[*idVOX].m_tau_prior;
			*m_tau=multifibres[*idVOX].m_tau;	
		}else{ 
			*m_tau_acc=0;
			*m_tau_rej=0;
			*m_tau_prop=0;
			*m_tau_prior=0;
			*m_tau=0;
		}
	}
	__syncthreads();

	int mydirs = ndirections/threadsBlock;
	int mod = ndirections%threadsBlock;
	if(mod&&(idSubVOX<mod)) mydirs++;

	if(idSubVOX<nfib){
		int pos = (*idVOX*nfib)+idSubVOX;
		m_th[idSubVOX]=fibres[pos].m_th;
		m_ph[idSubVOX]=fibres[pos].m_ph;
		m_f[idSubVOX]=fibres[pos].m_f;
	
		m_th_acc[idSubVOX]=fibres[pos].m_th_acc;
		m_th_rej[idSubVOX]=fibres[pos].m_th_rej;
		m_ph_acc[idSubVOX]=fibres[pos].m_ph_acc;
		m_ph_rej[idSubVOX]=fibres[pos].m_ph_rej;
		m_f_acc[idSubVOX]=fibres[pos].m_f_acc;
		m_f_rej[idSubVOX]=fibres[pos].m_f_rej;

		fm_prior_en[idSubVOX]=fibres[pos].m_prior_en;
		m_th_prior[idSubVOX]=fibres[pos].m_th_prior;
		m_ph_prior[idSubVOX]=fibres[pos].m_ph_prior;
		m_f_prior[idSubVOX]=fibres[pos].m_f_prior;

		m_th_prop[idSubVOX]=fibres[pos].m_th_prop;
		m_ph_prop[idSubVOX]=fibres[pos].m_ph_prop;
		m_f_prop[idSubVOX]=fibres[pos].m_f_prop;

		m_lam_jump[idSubVOX]=fibres[pos].m_lam_jump;		
	}

	__syncthreads();

	for(int i=0; i<mydirs; i++){	
		int pos = (*idVOX*ndirections)+idSubVOX+i*threadsBlock;
		myisosig[i] = isosignals[pos];
		mydata[i] = datam[pos];
	}
	if(gradnonlin){
		for(int i=0; i<mydirs; i++){	
			int pos = (*idVOX*ndirections)+idSubVOX+i*threadsBlock;
			mybvals[i] = bvals[pos];
			myalpha[i] = alpha[pos];
			mybeta[i] = beta[pos];
		}
	}else{
		for(int i=0; i<mydirs; i++){	
			int pos = idSubVOX+i*threadsBlock;
			mybvals[i] = bvals[pos];
			myalpha[i] = alpha[pos];
			mybeta[i] = beta[pos];
		}

	}

	for(int f=0;f<nfib;f++){
		for(int i=0; i<mydirs; i++){	
			mysig[i*nfib+f]= signals[(*idVOX*ndirections*nfib)+(f*ndirections)+idSubVOX+i*threadsBlock];	
		}
	}

	__syncthreads();
	
	//compute_signal_pre
	for(int i=0; i<mydirs; i++){
		for(int f=0;f<nfib;f++){	
			cos_alpha_minus_theta=cos(double(myalpha[i]-m_th[f]));   
			cos_alpha_plus_theta=cos(double(myalpha[i]+m_th[f]));
		     	angtmp[i*nfib+f]= (cos(double(m_ph[f]-mybeta[i]))*(cos_alpha_minus_theta-cos_alpha_plus_theta)/2)+(cos_alpha_minus_theta+cos_alpha_plus_theta)/2;
		 	angtmp[i*nfib+f]=angtmp[i*nfib+f]*angtmp[i*nfib+f];
		}
	}

	for (int niter=0; niter<iterations; niter++){
		//code jump()

		//prefetching randoms
		if (leader){
			*count_update=*count_update+1;
			*recordcount=*recordcount+1;
			int idrand = *idVOX*iterations*nparams+(niter*nparams);
			*localrand = 0;
			if(m_include_f0){
				prerandN[*localrand]=randomsN[idrand];
				prerandU[*localrand]=randomsU[idrand];
				idrand++;
				*localrand=*localrand+1;
			}
			if(rician){
				prerandN[*localrand]=randomsN[idrand];
				prerandU[*localrand]=randomsU[idrand];
				idrand++;
				*localrand=*localrand+1;
			}
			prerandN[*localrand]=randomsN[idrand];
			prerandU[*localrand]=randomsU[idrand];
			idrand++;
			*localrand=*localrand+1;

			if(model==2){
				prerandN[*localrand]=randomsN[idrand];
				prerandU[*localrand]=randomsU[idrand];
				idrand++;
				*localrand=*localrand+1;
			}

			prerandN[*localrand]=randomsN[idrand];
			prerandU[*localrand]=randomsU[idrand];
			idrand++;
			*localrand=*localrand+1;

			for(int f=0;f<nfib;f++){
				prerandN[*localrand]=randomsN[idrand];
				prerandU[*localrand]=randomsU[idrand];
				idrand++;
				*localrand=*localrand+1;

				prerandN[*localrand]=randomsN[idrand];
				prerandU[*localrand]=randomsU[idrand];
				idrand++;
				*localrand=*localrand+1;

				prerandN[*localrand]=randomsN[idrand];					
				prerandU[*localrand]=randomsU[idrand];	
				idrand++;	
				*localrand=*localrand+1;	
			}
			*localrand = 0;
		}

////////////////////////////////////////////////////////////////// F0

		if(m_include_f0){
			if (leader){
				//propose_f0
				old[0]=*m_f0;
				*m_f0+=prerandN[*localrand]**m_f0_prop;
				
				//compute_f0_prior()     
				old[1]=*m_f0_prior;
	      			if(*m_f0<=0 || *m_f0 >=1){ 
					rejflag[0]=true;
				}else{ 	
					rejflag[0]=false;
					if(!m_ardf0){
						*m_f0_prior=0;
	      				}else{
						*m_f0_prior=log(double(*m_f0));
					}
				}
				
				//for likehood and reject_f_sum
				*fsum=*m_f0;

				for(int g=0;g<nfib;g++){  
					*fsum+=m_f[g];
				}

				//compute_prior()
				compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,nfib);
				
				//reject_f_sum()
				rejflag[1]=(*fsum>1);
			}

			__syncthreads();

			//compute_likelihood()
			compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,mysig,myisosig,mydata,fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);

			__syncthreads();

			if (leader){
				//compute_energy()
				*m_old_energy=*m_energy;
	      			*m_energy=*m_prior_en+*m_likelihood_en;

				//test_energy
				*tmp=exp(double(*m_old_energy-*m_energy));
				rejflag[2]=(*tmp>prerandU[*localrand]);
				*localrand=*localrand+1;	

				if(!rejflag[0]){
					if(!rejflag[1]){
						if(rejflag[2]){
							//accept_f0()
		  					*m_f0_acc=*m_f0_acc+1;   
						}else{
							//reject_F0()
							*m_f0=old[0];
							*m_f0_prior=old[1];
		      					*m_prior_en=*m_old_prior_en;
		      					*m_f0_rej=*m_f0_rej+1;

							//restore_energy()
      							*m_energy=*m_old_energy;
						}
					}else{ 
						//reject_F0()
						*m_f0=old[0];
						*m_f0_prior=old[1];
		      				*m_prior_en=*m_old_prior_en;
		      				*m_f0_rej=*m_f0_rej+1;

						//restore_energy()
		      				*m_energy=*m_old_energy;
					}
				}else{ 
					//reject_F0()
					*m_f0=old[0];
					*m_f0_prior=old[1];
		      			*m_prior_en=*m_old_prior_en;
		      			*m_f0_rej=*m_f0_rej+1;

					//restore_energy()
		      			*m_energy=*m_old_energy;
				}
			}
		}
			
////////////////////////////////////////////////////////////////// TAU
		if(rician){
			if (leader){
				//propose_tau
				old[0]=*m_tau; 
				*m_tau+=prerandN[*localrand]**m_tau_prop;
			
				//compute_tau_prior()     
				old[1]=*m_tau_prior;
	      			if(*m_tau<=0){ 
					rejflag[0]=true;
				}else{ 	
					rejflag[0]=false;
					*m_tau_prior=0;
				}

				//for likehood and reject_f_sum
				*fsum=*m_f0;
				for(int g=0;g<nfib;g++){  
					*fsum+=m_f[g];
				}

				//compute_prior()
				compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,nfib);
			}

			__syncthreads();

			//compute_likelihood()
			compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,mysig,myisosig,mydata,fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);

			__syncthreads();

			if (leader){ 
				//compute_energy()
				*m_old_energy=*m_energy;
	      			*m_energy=*m_prior_en+*m_likelihood_en;

				//test_energy
				*tmp=exp(double(*m_old_energy-*m_energy));
				rejflag[1]=(*tmp>prerandU[*localrand]);	
				*localrand=*localrand+1;
				
				if(!rejflag[0]){
					if(rejflag[1]){
						//accept_tau()
		  				*m_tau_acc=*m_tau_acc+1;   
					}else{ 
						//reject_tau()
						*m_tau=old[0];
						*m_tau_prior=old[1];
		      				*m_prior_en=*m_old_prior_en;
		      				*m_tau_rej=*m_tau_rej+1;
						//restore_energy()
		      				*m_energy=*m_old_energy;
					}
				}else{ 
					//reject_tau()
					*m_tau=old[0];
					*m_tau_prior=old[1];
		      			*m_prior_en=*m_old_prior_en;
		      			*m_tau_rej=*m_tau_rej+1;

					//restore_energy()
		      			*m_energy=*m_old_energy;
				}
			}	
		}

////////////////////////////////////////////////////////////////// D

		if (leader){
			//propose_d()
			old[0]=*m_d;	
			*m_d+=prerandN[*localrand]**m_d_prop;
				
			//compute_d_prior_f0()      
			old[1]=*m_d_prior;
					
      			if(*m_d<0 || *m_d > 0.008){
				rejflag[0]=true;
			}else{ 
				*m_d_prior=0;
				rejflag[0]=false;
      			}
		}

		__syncthreads();
				
		//compute_signal()
		for(int i=0; i<mydirs; i++){
			for(int f=0;f<nfib;f++){
				compute_signal(&mysig[i*nfib+f],&mysigold[i*nfib+f],mybvals[i],m_d,m_dstd,angtmp[i*nfib+f],model);
			}
		}
		//compute_iso_signal()
		for(int i=0; i<mydirs; i++){
			compute_iso_signal(&myisosig[i],&myisosigold[i], mybvals[i],m_d,m_dstd,model);
		}				

		if (leader){
			//for likehood and reject_f_sum
			*fsum=*m_f0;

			for(int g=0;g<nfib;g++){  
				*fsum+=m_f[g];
			}

			//compute_prior()
			compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,nfib);
		}

		__syncthreads();
				
		//compute_likelihood()
		compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,mysig,myisosig,mydata,fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);	
				
		__syncthreads();
				
		if (leader){
			//compute_energy()
			*m_old_energy=*m_energy;
      			*m_energy=*m_prior_en+*m_likelihood_en;
			//test_energy
			*tmp=exp(double(*m_old_energy-*m_energy));
			rejflag[1]=(*tmp>prerandU[*localrand]);
			*localrand=*localrand+1;
		}
			
		__syncthreads();

       		if(!rejflag[0]){
			if(rejflag[1]){
				//accept_d()
	  			if (leader) *m_d_acc=*m_d_acc+1;   
			}else{
				if (leader){
					//reject_d()
					*m_d=old[0];
					*m_d_prior=old[1];
      					*m_prior_en=*m_old_prior_en;
					*m_d_rej=*m_d_rej+1;
					//restore_energy()
      					*m_energy=*m_old_energy;
				}
      						
				for(int i=0; i<mydirs; i++){
					for(int f=0;f<nfib;f++){
						mysig[i*nfib+f] = mysigold[i*nfib+f];
					}	
				
					myisosig[i]=myisosigold[i];
				}
			}
        	}else{ 
			if (leader){
				//reject_d()
				*m_d=old[0];
				*m_d_prior=old[1];
      				*m_prior_en=*m_old_prior_en;
				*m_d_rej=*m_d_rej+1;
				//restore_energy()
      				*m_energy=*m_old_energy;
			}

      			for(int i=0; i<mydirs; i++){
				for(int f=0;f<nfib;f++){
					mysig[i*nfib+f] = mysigold[i*nfib+f];
				}	
				
				myisosig[i]=myisosigold[i];
			}	
        	}

////////////////////////////////////////////////////////////////// D_STD

		if(model==2){
			if (leader){	
				//propose_d_std
				old[0]=*m_dstd;
				*m_dstd+=prerandN[*localrand]**m_dstd_prop;

				//compute_d_std_prior()     
				old[1]=*m_dstd_prior;
	      			if(*m_dstd<=0 || *m_dstd > 0.01){ 
					rejflag[0]=true;
				}else{ 	
					rejflag[0]=false;
					*m_dstd_prior=log(double(*m_dstd));
				}
			}

			__syncthreads();

			//compute_signal()
			for(int i=0; i<mydirs; i++){
				for(int f=0;f<nfib;f++){
					compute_signal(&mysig[i*nfib+f],&mysigold[i*nfib+f],mybvals[i],m_d,m_dstd,angtmp[i*nfib+f],model);
				}
			}

			//compute_iso_signal()
			for(int i=0; i<mydirs; i++){
				compute_iso_signal(&myisosig[i],&myisosigold[i], mybvals[i],m_d,m_dstd,model);
			}

			if (leader){	
				//for likehood and reject_f_sum
				*fsum=*m_f0;
				for(int g=0;g<nfib;g++){  
					*fsum+=m_f[g];
				}

				//compute_prior()
				compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,nfib);
			}

			__syncthreads();

			//compute_likelihood()
			compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,mysig,myisosig,mydata,fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);

			__syncthreads();

			if (leader){
				//compute_energy()
				*m_old_energy=*m_energy;
	      			*m_energy=*m_prior_en+*m_likelihood_en;

				//test_energy
				*tmp=exp(double(*m_old_energy-*m_energy));
				rejflag[1]=(*tmp>prerandU[*localrand]);
				*localrand=*localrand+1;					
			}

			__syncthreads();
				
			if(!rejflag[0]){
				if(rejflag[1]){
					//accept_dstd()
		  			if (leader) *m_dstd_acc=*m_dstd_acc+1;   
				}else{ 
					if (leader){
						//reject_dstd()
						*m_dstd=old[0];
						*m_dstd_prior=old[1];
		      				*m_prior_en=*m_old_prior_en;
		      				*m_dstd_rej=*m_dstd_rej+1;

						//restore_energy()
		      				*m_energy=*m_old_energy;
					}

					for(int i=0; i<mydirs; i++){
						for(int f=0;f<nfib;f++){
							mysig[i*nfib+f] = mysigold[i*nfib+f];
						}
				
						myisosig[i]=myisosigold[i];
					}
				}
			}else{ 
				if (leader){
					//reject_dstd()
					*m_dstd=old[0];
					*m_dstd_prior=old[1];
		      			*m_prior_en=*m_old_prior_en;
		      			*m_dstd_rej=*m_dstd_rej+1;

					//restore_energy()
		      			*m_energy=*m_old_energy;
				}

				for(int i=0; i<mydirs; i++){
					for(int f=0;f<nfib;f++){
						mysig[i*nfib+f] = mysigold[i*nfib+f];
					}	
				
					myisosig[i]=myisosigold[i];
				}
			}
		}

////////////////////////////////////////////////////////////////// S0

		if (leader){
			//propose_S0()
			old[0]=*m_S0;
			*m_S0+=prerandN[*localrand]**m_S0_prop;
				
			//compute_S0_prior()
			old[1]=*m_S0_prior;
        		if(*m_S0<0) rejflag[0]=true;
        		else{    
				*m_S0_prior=0;
	  			rejflag[0]=false;
        		}
				
			//for likehood and reject_f_sum
			*fsum=*m_f0;
			for(int g=0;g<nfib;g++){  
				*fsum+=m_f[g];
			}

			//compute_prior()
			compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,nfib);
		}

		__syncthreads();

		//compute_likelihood()
		compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,mysig,myisosig,mydata,fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);

		__syncthreads();

		if (leader){
			//compute_energy()
			*m_old_energy=*m_energy;
      			*m_energy=*m_prior_en+*m_likelihood_en;

			//test_energy
			*tmp=exp(double(*m_old_energy-*m_energy));
			rejflag[1]=(*tmp>prerandU[*localrand]);
			*localrand=*localrand+1;

        		if(!rejflag[0]){
				if(rejflag[1]){
					//accept_S0()
	  				*m_S0_acc=*m_S0_acc+1;   
				}else{
					//reject_S0()
					*m_S0=old[0];
					*m_S0_prior=old[1];
      					*m_prior_en=*m_old_prior_en;
      					*m_S0_rej=*m_S0_rej+1;

					//restore_energy()
      					*m_energy=*m_old_energy;
				}
        		}else{ 
				//reject_S0()
				*m_S0=old[0];
				*m_S0_prior=old[1];
	      			*m_prior_en=*m_old_prior_en;
	      			*m_S0_rej=*m_S0_rej+1;

				//restore_energy()
      				*m_energy=*m_old_energy;
			}
        	}

////////////////////////////////////////////////////////////////////////////     TH

     		for(int fibre=0;fibre<nfib;fibre++){  
			if (leader){ 
				//propose_th()
				old[0]=m_th[fibre];
				m_th[fibre]+=prerandN[*localrand]*m_th_prop[fibre];
					
				//compute_th_prior()
				old[1]=m_th_prior[fibre];
      	   			if(m_th[fibre]==0){
					m_th_prior[fibre]=0;
		   		}else{
					m_th_prior[fibre]=-log(double(fabs(sin(double(m_th[fibre]))/2)));
	      	   		}
		  		//rejflag[0]=false; /////////////////always false
		
				//compute_prior()
				*fm_old_prior_en=fm_prior_en[fibre];
	      	   		fm_prior_en[fibre]=m_th_prior[fibre]+m_ph_prior[fibre]+m_f_prior[fibre];	
			}

			__syncthreads();

			//compute_signal()
			//compute_signal_pre	
			for(int i=0; i<mydirs; i++){
				cos_alpha_minus_theta=cos(double(myalpha[i]-m_th[fibre]));   
				cos_alpha_plus_theta=cos(double(myalpha[i]+m_th[fibre]));
				old_angtmp[i*nfib+fibre]=angtmp[i*nfib+fibre];
				angtmp[i*nfib+fibre]= (cos(double(m_ph[fibre]-mybeta[i]))*(cos_alpha_minus_theta-cos_alpha_plus_theta)/2)+(cos_alpha_minus_theta+cos_alpha_plus_theta)/2;
		 	   	angtmp[i*nfib+fibre]=angtmp[i*nfib+fibre]*angtmp[i*nfib+fibre];

				compute_signal(&mysig[i*nfib+fibre],&mysigold[i*nfib+fibre],mybvals[i],m_d,m_dstd,angtmp[i*nfib+fibre],model);
			}

			if (leader){	
				//for likehood and reject_f_sum
				*fsum=*m_f0;
				for(int g=0;g<nfib;g++){  
					*fsum+=m_f[g];
				}

				//compute_prior()
				compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,nfib);	
			}

			__syncthreads();

			//compute_likelihood()
			compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,mysig,myisosig,mydata,fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);

			__syncthreads();

			if (leader){ 
				//compute_energy()
				*m_old_energy=*m_energy;
	      			*m_energy=*m_prior_en+*m_likelihood_en;

				//test_energy
				*tmp=exp(double(*m_old_energy-*m_energy));
				rejflag[1]=(*tmp>prerandU[*localrand]);
				*localrand=*localrand+1;					
			}

			__syncthreads();
			
			if(rejflag[1]){
				//accept_th
		  		if (leader) m_th_acc[fibre]++;   
				}else{

				if (leader){
					//reject_th()
					m_th[fibre]=old[0];
					m_th_prior[fibre]=old[1];
      					fm_prior_en[fibre]=*fm_old_prior_en;
					m_th_rej[fibre]++;						
						
					//restore_energy()
					*m_prior_en=*m_old_prior_en;
	      				*m_energy=*m_old_energy;
				}

				//compute_signal_pre undo
				for(int i=0; i<mydirs; i++){
					angtmp[i*nfib+fibre]=old_angtmp[i*nfib+fibre];
					mysig[i*nfib+fibre] = mysigold[i*nfib+fibre];						
      				}
			}
			__syncthreads();
		
///////////////////////////////////////     PH

			if (leader){
				//propose_ph()
				old[0]=m_ph[fibre];
				m_ph[fibre]+=prerandN[*localrand]*m_ph_prop[fibre];

				//compute_ph_prior()
				old[1]=m_ph_prior[fibre];
      				m_ph_prior[fibre]=0;
      				//rejflag[0]=false;

				//compute_prior()
				*fm_old_prior_en=fm_prior_en[fibre];
      	   			fm_prior_en[fibre]=m_th_prior[fibre]+m_ph_prior[fibre]+m_f_prior[fibre];
			}

			__syncthreads();
				
			//compute_signal()
			//compute_signal_pre
			for(int i=0; i<mydirs; i++){
				cos_alpha_minus_theta=cos(double(myalpha[i]-m_th[fibre]));   
			  	cos_alpha_plus_theta=cos(double(myalpha[i]+m_th[fibre]));
				old_angtmp[i*nfib+fibre]=angtmp[i*nfib+fibre];
		     	   	angtmp[i*nfib+fibre]= (cos(double(m_ph[fibre]-mybeta[i]))*(cos_alpha_minus_theta-cos_alpha_plus_theta)/2)+(cos_alpha_minus_theta+cos_alpha_plus_theta)/2;
		 	   	angtmp[i*nfib+fibre]=angtmp[i*nfib+fibre]*angtmp[i*nfib+fibre];
	
				compute_signal(&mysig[i*nfib+fibre],&mysigold[i*nfib+fibre],mybvals[i],m_d,m_dstd,angtmp[i*nfib+fibre],model);
			}

			if (leader){
				//for likehood and reject_f_sum
				*fsum=*m_f0;
				for(int g=0;g<nfib;g++){  
					*fsum+=m_f[g];
				}

				//compute_prior()
				compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,nfib);
			}

			__syncthreads();

			//compute_likelihood()
			compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,mysig,myisosig,mydata,fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);


			__syncthreads();

			if (leader){
				//compute_energy()
				*m_old_energy=*m_energy;
	      			*m_energy=*m_prior_en+*m_likelihood_en;

				//test_energy
				*tmp=exp(double(*m_old_energy-*m_energy));
				rejflag[1]=(*tmp>prerandU[*localrand]);
				*localrand=*localrand+1;
			}

			__syncthreads();

			//if(!rejflag[0]){
			if(rejflag[1]){
				//accept_ph()
		  		if (leader) m_ph_acc[fibre]++;   
			}else{
				if (leader){
					//reject_ph()
					m_ph[fibre]=old[0];
					m_ph_prior[fibre]=old[1];
      					fm_prior_en[fibre]=*fm_old_prior_en;
					m_ph_rej[fibre]++;						
						
					//restore_energy()
					*m_prior_en=*m_old_prior_en;
	      				*m_energy=*m_old_energy;
				}
				//compute_signal_pre undo
				for(int i=0; i<mydirs; i++){
					angtmp[i*nfib+fibre]=old_angtmp[i*nfib+fibre];

					mysig[i*nfib+fibre] = mysigold[i*nfib+fibre];	
				}
			}

			__syncthreads();

////////////////////////////////////////////             F

			if (leader){
				//propose_f()
				old[0]=m_f[fibre];
				m_f[fibre]+=prerandN[*localrand]*m_f_prop[fibre];

	     			//compute_f_prior()
	        		old[1]=m_f_prior[fibre];
				if (m_f[fibre]<=0 || m_f[fibre]>=1) rejflag[0]=true;
	        		else{
		      			if(!can_use_ard ){
		  				m_f_prior[fibre]=0;
					}else{
		  				if(m_lam_jump[fibre]){
							m_f_prior[fibre]=log(double(m_f[fibre]));
						}else{
		    					m_f_prior[fibre]=0;
		  				}
					}
					m_f_prior[fibre]=fudgevalue*m_f_prior[fibre];
					rejflag[0]=false;
	      			}

				//compute_prior()
				*fm_old_prior_en=fm_prior_en[fibre];
      	   			fm_prior_en[fibre]=m_th_prior[fibre]+m_ph_prior[fibre]+m_f_prior[fibre];
						
				//for likehood and reject_f_sum
				*fsum=*m_f0;
				for(int g=0;g<nfib;g++){  
					*fsum+=m_f[g];
				}
				//reject_f_sum()
				rejflag[1]=(*fsum>1);

				//compute_prior()
				compute_prior(m_prior_en,m_old_prior_en,m_d_prior,m_S0_prior,fm_prior_en,m_f0_prior,m_tau_prior,m_dstd_prior,nfib);	
			}

			__syncthreads();

			//compute_likelihood()
			compute_likelihood(idSubVOX,m_S0,m_likelihood_en,m_f,mysig,myisosig,mydata,fsum,reduction,m_f0,rician,m_tau,mydirs,threadsBlock,ndirections,nfib);	

			__syncthreads();

			if (leader){
				//compute_energy()
				*m_old_energy=*m_energy;
		      		*m_energy=*m_prior_en+*m_likelihood_en;

				//test_energy
				*tmp=exp(double(*m_old_energy-*m_energy));
				rejflag[2]=(*tmp>prerandU[*localrand]);
				*localrand=*localrand+1;			

		      		if(!rejflag[0]){
					if(!rejflag[1]){
						if(rejflag[2]){
							//accept_f()
			  				m_f_acc[fibre]++;   
						}else{
							//reject_f()
							m_f[fibre]=old[0];
							m_f_prior[fibre]=old[1];
	      						fm_prior_en[fibre]=*fm_old_prior_en;
	      						m_f_rej[fibre]++;				
					
							//restore_energy()
							*m_prior_en=*m_old_prior_en;
		      					*m_energy=*m_old_energy;
						}
					}else{ 
						//reject_f()
						m_f[fibre]=old[0];
						m_f_prior[fibre]=old[1];
	      					fm_prior_en[fibre]=*fm_old_prior_en;
	      					m_f_rej[fibre]++;

						//restore_energy()
						*m_prior_en=*m_old_prior_en;
		      				*m_energy=*m_old_energy;	
					}
				}else{
					//reject_f()
					m_f[fibre]=old[0];
					m_f_prior[fibre]=old[1];
		      			fm_prior_en[fibre]=*fm_old_prior_en;
		      			m_f_rej[fibre]++;

					//restore_energy()
					*m_prior_en=*m_old_prior_en;
		      			*m_energy=*m_old_energy;
				}
			}
			__syncthreads();	

        	}//end while nfib

		if(((*recordcount%record_every)==0)&&leader){
			rd[(*idVOX*totalrecords)+*sample-1]= *m_d;
			if(m_include_f0) rf0[(*idVOX*totalrecords)+*sample-1]= *m_f0;
			if(rician) rtau[(*idVOX*totalrecords)+*sample-1]= *m_tau;
			if(model==2) rdstd[(*idVOX*totalrecords)+*sample-1]= *m_dstd;
			rs0[(*idVOX*totalrecords)+*sample-1]= *m_S0;
			for(int j=0;j<nfib;j++){
				rth[(*idVOX*totalrecords*nfib)+(j*totalrecords)+*sample-1]=m_th[j];
				rph[(*idVOX*totalrecords*nfib)+(j*totalrecords)+*sample-1]=m_ph[j];
				rf[(*idVOX*totalrecords*nfib)+(j*totalrecords)+*sample-1]=m_f[j];
			}
			*sample=*sample+1;
        	}

		if(((*count_update%updateproposalevery)==0)&&leader){
			//m_multifibre.update_proposals();
			*m_d_prop*=sqrt(float(*m_d_acc+1)/float(*m_d_rej+1));
			*m_d_prop=min(*m_d_prop,maxfloat);

			if(rician){
				*m_tau_prop*=sqrt(float(*m_tau_acc+1)/float(*m_tau_rej+1));
				*m_tau_prop=min(*m_tau_prop,maxfloat);
				*m_tau_acc=0; 
				*m_tau_rej=0;	
			}

			if(m_include_f0){
				*m_f0_prop*=sqrt(float(*m_f0_acc+1)/float(*m_f0_rej+1));
				*m_f0_prop=min(*m_f0_prop,maxfloat);
				*m_f0_acc=0; 
				*m_f0_rej=0;	
			}	

			if(model==2){
				*m_dstd_prop*=sqrt(float(*m_dstd_acc+1)/float(*m_dstd_rej+1));
				*m_dstd_prop=min(*m_dstd_prop,maxfloat);
				*m_dstd_acc=0; 
				*m_dstd_rej=0;	
			}

			*m_S0_prop*=sqrt(float(*m_S0_acc+1)/float(*m_S0_rej+1));
			*m_S0_prop=min(*m_S0_prop,maxfloat);
			*m_d_acc=0; 
			*m_d_rej=0;
			*m_S0_acc=0; 
			*m_S0_rej=0;
			for(int f=0; f<nfib;f++){
				//m_fibres[f].update_proposals();
				m_th_prop[f]*=sqrt(float(m_th_acc[f]+1)/float(m_th_rej[f]+1));
				m_th_prop[f]=min(m_th_prop[f],maxfloat);
		      		m_ph_prop[f]*=sqrt(float(m_ph_acc[f]+1)/float(m_ph_rej[f]+1));
		      		m_ph_prop[f]=min(m_ph_prop[f],maxfloat);
		      		m_f_prop[f]*=sqrt(float(m_f_acc[f]+1)/float(m_f_rej[f]+1));
		      		m_f_prop[f]=min(m_f_prop[f],maxfloat);
			      
		      		m_th_acc[f]=0; 
		      		m_th_rej[f]=0;
		      		m_ph_acc[f]=0; 
		      		m_ph_rej[f]=0;
		      		m_f_acc[f]=0; 
		      		m_f_rej[f]=0;
			}
		}

		__syncthreads();	

        } //end while iterations

	if(leader){
		multifibres[*idVOX].m_S0=*m_S0;
		multifibres[*idVOX].m_S0_prior=*m_S0_prior;
		multifibres[*idVOX].m_S0_prop=*m_S0_prop;
		multifibres[*idVOX].m_S0_acc=*m_S0_acc;
		multifibres[*idVOX].m_S0_rej=*m_S0_rej;

		multifibres[*idVOX].m_d=*m_d;
		multifibres[*idVOX].m_d_prior=*m_d_prior;
		multifibres[*idVOX].m_d_prop=*m_d_prop;
		multifibres[*idVOX].m_d_acc=*m_d_acc;
		multifibres[*idVOX].m_d_rej=*m_d_rej;
	
		multifibres[*idVOX].m_prior_en=*m_prior_en;
		multifibres[*idVOX].m_energy=*m_energy;
		multifibres[*idVOX].m_likelihood_en=*m_likelihood_en;

		if(m_include_f0){
			multifibres[*idVOX].m_f0_prior=*m_f0_prior;
			multifibres[*idVOX].m_f0=*m_f0;
			multifibres[*idVOX].m_f0_acc=*m_f0_acc;
			multifibres[*idVOX].m_f0_rej=*m_f0_rej;
			multifibres[*idVOX].m_f0_prop=*m_f0_prop;
		}
		if(rician){
			multifibres[*idVOX].m_tau_prior=*m_tau_prior;
			multifibres[*idVOX].m_tau=*m_tau;
			multifibres[*idVOX].m_tau_acc=*m_tau_acc;
			multifibres[*idVOX].m_tau_rej=*m_tau_rej;
			multifibres[*idVOX].m_tau_prop=*m_tau_prop;
		}
		if(model==2){
			multifibres[*idVOX].m_dstd_prior=*m_dstd_prior;
			multifibres[*idVOX].m_dstd=*m_dstd;
			multifibres[*idVOX].m_dstd_acc=*m_dstd_acc;
			multifibres[*idVOX].m_dstd_rej=*m_dstd_rej;
			multifibres[*idVOX].m_dstd_prop=*m_dstd_prop;
		}
	}
	
	if(idSubVOX<nfib){
		int pos = (*idVOX*nfib)+idSubVOX;
	
		fibres[pos].m_th=m_th[idSubVOX];
		fibres[pos].m_ph=m_ph[idSubVOX];
		fibres[pos].m_f=m_f[idSubVOX];

		fibres[pos].m_th_acc=m_th_acc[idSubVOX];
		fibres[pos].m_th_rej=m_th_rej[idSubVOX];
		fibres[pos].m_ph_acc=m_ph_acc[idSubVOX];
		fibres[pos].m_ph_rej=m_ph_rej[idSubVOX];
		fibres[pos].m_f_acc=m_f_acc[idSubVOX];
		fibres[pos].m_f_rej=m_f_rej[idSubVOX];

		fibres[pos].m_prior_en=fm_prior_en[idSubVOX];
		fibres[pos].m_th_prior=m_th_prior[idSubVOX];
		fibres[pos].m_ph_prior=m_ph_prior[idSubVOX];
		fibres[pos].m_f_prior=m_f_prior[idSubVOX];

		fibres[pos].m_th_prop=m_th_prop[idSubVOX];
		fibres[pos].m_ph_prop=m_ph_prop[idSubVOX];
		fibres[pos].m_f_prop=m_f_prop[idSubVOX];

		fibres[pos].m_lam_jump=m_lam_jump[idSubVOX];		
	}

	for(int i=0; i<mydirs; i++){	
		isosignals[(*idVOX*ndirections)+idSubVOX+i*threadsBlock] = myisosig[i];

		for(int f=0;f<nfib;f++){
			signals[(*idVOX*ndirections*nfib)+(f*ndirections)+idSubVOX+i*threadsBlock]=mysig[i*nfib+f];
		}
	}
}

