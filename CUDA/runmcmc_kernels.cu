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

__device__ 
static const float maxfloat=1e10;

__device__
inline void compute_signal(double *signals,double *oldsignals,double mbvals,float m_d,int model, float m_dstd,float angtmp){

	oldsignals[0]=signals[0];

	if(model==1 || m_dstd<1e-5){	
		signals[0]=exp(double(-m_d*mbvals*angtmp));
	}else{
		float dbeta=__ddiv_rn(m_d,(m_dstd*m_dstd));
	 	float dalpha=m_d*dbeta;         
           	signals[0]=exp(double(log(double(__ddiv_rn(dbeta,__dadd_rn(dbeta,mbvals*angtmp))))*dalpha)); 
	}
}

__device__ 
inline void compute_iso_signal(double *isosignals,double *oldisosignals, double mbvals,float m_d, float m_dstd, int model){

	*oldisosignals=*isosignals;

	 if(model==1 || m_dstd<1e-5){
	 	*isosignals=exp(-m_d*mbvals);	
	}else{
		float dbeta=__ddiv_rn(m_d,(m_dstd*m_dstd));
	  	float dalpha=m_d*dbeta;
		*isosignals=exp(double(log(double(__ddiv_rn(dbeta,__dadd_rn(dbeta,mbvals))))*dalpha));
	}
}

__device__ 
inline void compute_prior(float *m_prior_en, float *m_prior_en_old,float m_d_prior,float m_S0_prior,float *m_prior_enf, float m_f0_prior, float m_tau_prior, float m_dstd_prior){			
        *m_prior_en_old=*m_prior_en;
	*m_prior_en=m_d_prior+m_S0_prior;
	*m_prior_en=*m_prior_en+m_dstd_prior;
	*m_prior_en=*m_prior_en+m_tau_prior;
	*m_prior_en=*m_prior_en+m_f0_prior;	
	for(int f=0;f<NFIBRES;f++){
		*m_prior_en=*m_prior_en+m_prior_enf[f];
	}	
}

__device__ 
inline float logIo(const float& x){
    	float y,b;
    	b=std::fabs(x);
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
inline void compute_likelihood(bool leader,int idrest,float m_d,float m_S0,float *m_likelihood_en,float *m_f,double *signals,double *isosignals,double *mdata,float fsum,double *reduction, double m_f0, const bool rician, double m_tau,int ndirs){
	

	double sumsquares=0;
	double pred[MAXNDIRS_PER_THREAD];

	for(int i=0; i<ndirs; i++){
		pred[i]=0;
	      	for(int f=0;f<NFIBRES;f++){
			pred[i]=__dadd_rn (pred[i],__dmul_rn (m_f[f],signals[i*NFIBRES+f]));
	     	}
	}

	for(int i=0; i<ndirs; i++){
		pred[i]=m_S0*__dadd_rn(__dadd_rn (pred[i],__dmul_rn ((1-fsum),isosignals[i])),m_f0); //F0
	}

	reduction[idrest]=0;
	for(int i=0; i<ndirs; i++){	
		if(!rician) 
			reduction[idrest] = __dadd_rn(reduction[idrest],__dmul_rn(__dadd_rn(mdata[i],-pred[i]),__dadd_rn(mdata[i],-pred[i])));
		else{
			pred[i]=__dadd_rn(log(mdata[i]),__dadd_rn(-0.5*m_tau*__dadd_rn(mdata[i]*mdata[i],pred[i]*pred[i]),logIo(m_tau*pred[i]*mdata[i])));  
			reduction[idrest] = __dadd_rn(reduction[idrest],pred[i]);
		}
	}

	__syncthreads();

	unsigned int s2=THREADS_BLOCK;
	for(unsigned int s=THREADS_BLOCK/2; s>0; s>>=1) {
		if((s2%2)&&(idrest==(s-1))) reduction[idrest]= reduction[idrest] + reduction[idrest + s +1]; 
        	if (idrest < s){
            		reduction[idrest] = reduction[idrest] + reduction[idrest + s];
       	 	}
		s2=s;
        	__syncthreads();
    	}
	if(leader){
		sumsquares+=reduction[0];
		if(!rician){ 
		 	*m_likelihood_en=(NDIRECTIONS/2.0)*log(sumsquares/2.0); 
		}else{
			*m_likelihood_en=__dadd_rn(-NDIRECTIONS*log(m_tau),-sumsquares);
		}
	}

	//REDUCTION WITH ONLY ONE THREAD. THE LEADER
	/*if(leader){
		if(!rician){ 
			for(int i=0;i<THREADS_BLOCK;i++){
				sumsquares+=reduction[i]; //normal sum
			}
      			*m_likelihood_en=(NDIRECTIONS/2.0)*log(double(sumsquares/2.0)); 
		}else{
			for(int i=0;i<THREADS_BLOCK;i++){
				sumsquares+=reduction[i]; //normal sum
			}
			*m_likelihood_en=__dadd_rn(-NDIRECTIONS*log(m_tau),-sumsquares);
		}
	}*/
}
						  
extern "C" __global__ void init_Fibres_Multifibres_kernel(	//INPUT
								const double*			data,
								const double*			params,
								const float*			tau,
								const double*			bvals,
								const double*			alpha,
								const double*			beta,
								const int 			nvox,	
								const int 			nfib,
								const int 			nparams_fit,
								const int 			model,
								const float 			fudgevalue,
								const bool			m_includef0,
								const bool			m_rician,
								const bool 			m_ardf0,	// opts.ardf0.value()
								const bool 			ard_value,	// opts.all_ard.value()
								const bool 			no_ard_value,	// opts.no_ard.value()
								//OUTPUT
								FibreGPU*			fibres,
								MultifibreGPU*			multifibres,
								double*				signals,
								double*				isosignals)
{
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	if (id >= nvox*blockDim.x) { return; }

	int add=0;
	if(model==2) add=1;		// if model 2 we have d_std and then 1 more parameter in position 2

	int idmult= blockIdx.x;		// this is the voxel number
	int idrest= threadIdx.x;	// and the number of thread inside each voxel
	bool leader = (idrest==0);	// is the first one in the voxel?

	int ndirs = 1;
	if(idrest<(NDIRECTIONS%THREADS_BLOCK)) ndirs=MAXNDIRS_PER_THREAD;
	else if((NDIRECTIONS%THREADS_BLOCK)==0) ndirs=MAXNDIRS_PER_THREAD;
	else ndirs=(MAXNDIRS_PER_THREAD-1);

	if(idrest>=NDIRECTIONS)  { return; }

	__shared__ float fsum;					
	__shared__ double reduction[THREADS_BLOCK];		
	__shared__ float m_d;					
	__shared__ float m_S0;					
	__shared__ float m_f0;					
	__shared__ float m_dstd;
	__shared__ float m_tau;
	__shared__ float m_th[NFIBRES];
	__shared__ float m_ph[NFIBRES];
	__shared__ float m_f[NFIBRES];
	__shared__ float m_likelihood_en;
	__shared__ double mydata[MAXNDIRS_PER_THREAD*THREADS_BLOCK];

	__shared__ float m_prior_en;

	double mysig[MAXNDIRS_PER_THREAD*NFIBRES];			
	double myisosig[MAXNDIRS_PER_THREAD];	
	
	for(int i=0; i<ndirs; i++){	
		mydata[idrest*MAXNDIRS_PER_THREAD+i] = data[(idmult*NDIRECTIONS)+idrest+i*THREADS_BLOCK];
	}	

	// m_s0-params[0] 	m_d-params[1] 	m_f-m_th-m_ph-params[2,3,4,5, etc..]   	m_f0-params[nparams-1]

	if(leader){
		m_S0 = params[idmult*nparams_fit];
		multifibres[idmult].m_S0 = m_S0;
		multifibres[idmult].m_S0_prior = 0;
		multifibres[idmult].m_S0_acc = 0;
		multifibres[idmult].m_S0_rej = 0;
	
		m_d=params[idmult*nparams_fit+1];
		if(m_d<0 || m_d>0.008) m_d=2e-3;			//this is in xfibres...after fit
		multifibres[idmult].m_d = m_d;
		multifibres[idmult].m_d_prior = 0;
		multifibres[idmult].m_d_acc = 0;
		multifibres[idmult].m_d_rej = 0;

		if(model==2){ 
			m_dstd=params[idmult*nparams_fit+2];
			if(m_dstd<0 || m_dstd>0.01) m_dstd=m_d/10;	//this is in xfibres...after fit
		}
		else m_dstd = 0;
		multifibres[idmult].m_dstd = m_dstd;
		multifibres[idmult].m_dstd_prior = 0;
		multifibres[idmult].m_dstd_acc = 0;
		multifibres[idmult].m_dstd_rej = 0;

		if (m_includef0) m_f0=params[idmult*nparams_fit+nparams_fit-1];
		else m_f0=0;
		multifibres[idmult].m_f0 = m_f0;
		multifibres[idmult].m_f0_prior = 0;
		multifibres[idmult].m_f0_acc = 0;
		multifibres[idmult].m_f0_rej = 0;

		m_tau = tau[idmult];
		multifibres[idmult].m_tau = m_tau;
		multifibres[idmult].m_tau_prior = 0;
		multifibres[idmult].m_tau_acc = 0;
		multifibres[idmult].m_tau_rej = 0;
           
		for(int fib=0;fib<nfib;fib++){
			m_th[fib]=params[idmult*nparams_fit+2+3*fib+1+add];
			fibres[idmult*nfib+fib].m_th = m_th[fib];
			fibres[idmult*nfib+fib].m_th_prop = 0.2;
			float m_th_prior = 0;
			fibres[idmult*nfib+fib].m_th_acc = 0;
			fibres[idmult*nfib+fib].m_th_rej = 0;
			//compute_th_prior();
	      		if(m_th==0){
				m_th_prior=0;
			}else{
				m_th_prior=-log(double(fabs(sin(double(m_th[fib]))/2)));
	      		}
			fibres[idmult*nfib+fib].m_th_prior = m_th_prior;
		
			float m_ph_prior=0;	//compute_ph_prior();
			m_ph[fib]=params[idmult*nparams_fit+2+3*fib+2+add];
			fibres[idmult*nfib+fib].m_ph = m_ph[fib];
			fibres[idmult*nfib+fib].m_ph_prop = 0.2;
			fibres[idmult*nfib+fib].m_ph_prior = 0;	//compute_ph_prior();
			fibres[idmult*nfib+fib].m_ph_acc = 0;
			fibres[idmult*nfib+fib].m_ph_rej = 0;

			m_f[fib] = params[idmult*nparams_fit+2+3*fib+add]; 
			fibres[idmult*nfib+fib].m_f=m_f[fib];
			fibres[idmult*nfib+fib].m_f_prop = 0.2;
			float m_f_prior = 0;
			fibres[idmult*nfib+fib].m_f_acc = 0;
			fibres[idmult*nfib+fib].m_f_rej = 0;
			
			if(fib==0){
				fibres[idmult*nfib+fib].m_lam_jump = ard_value;
			}else{
				fibres[idmult*nfib+fib].m_lam_jump = !no_ard_value;
			}

			//compute_f_prior();
      			if (m_f[fib]<=0 | m_f[fib]>=1 ){
      			}else{
	  			if(fibres[idmult*nfib+fib].m_lam_jump){              
	    				m_f_prior=log(double(m_f[fib]));
	  			}else{
	    				m_f_prior=0;
				}
				m_f_prior= fudgevalue* m_f_prior;
      			}
			fibres[idmult*nfib+fib].m_f_prior = m_f_prior;

			//fibres[vox].m_lam = m_lam; ??
			//fibres[vox].m_lam_prop = 1;
			//fibres[vox].m_lam_prior = 0;
			//compute_lam_prior();

			//compute_prior();
			fibres[idmult*nfib+fib].m_prior_en= m_th_prior + m_ph_prior + m_f_prior;
		}
		//for likehood 
		fsum=m_f0;
		for(int g=0;g<NFIBRES;g++){  
			fsum+=m_f[g];
		}

	}
	__syncthreads();

	float angtmp[MAXNDIRS_PER_THREAD*NFIBRES];
	float cos_alpha_minus_theta;
	float cos_alpha_plus_theta;

	double old;

	//compute_signal_pre
	for(int i=0; i<ndirs; i++){
		for(int f=0;f<NFIBRES;f++){				
			cos_alpha_minus_theta=cos(double(alpha[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]-m_th[f]));   
			cos_alpha_plus_theta=cos(double(alpha[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]+m_th[f]));
		     	angtmp[i*NFIBRES+f]= __dadd_rn(__ddiv_rn (cos(double(m_ph[f]-beta[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]))*__dadd_rn(cos_alpha_minus_theta,-cos_alpha_plus_theta),2) ,__ddiv_rn (__dadd_rn(cos_alpha_minus_theta,cos_alpha_plus_theta),2));
		 	angtmp[i*NFIBRES+f]=angtmp[i*NFIBRES+f]*angtmp[i*NFIBRES+f];
		}
	}

	//compute_signal()
	for(int i=0; i<ndirs; i++){
		for(int f=0;f<NFIBRES;f++){
			compute_signal(&mysig[i*NFIBRES+f],&old,bvals[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK],m_d,model,m_dstd,angtmp[i*NFIBRES+f]);
		}
	}

	if(leader){
		//initialise_energies();
	      	//compute_d_prior(); m_d_prior=0; so i don't do nothing, it is already 0
	      	if(model==2){
			//compute_d_std_prior();
			if(m_dstd<=0 || m_dstd>0.01){
	      		}else{
				multifibres[idmult].m_dstd_prior=log(double(m_dstd));
			}
		}
	      	//compute_tau_prior(); m_tau_prior=0; so i don't do nothing, it is already 0
	      	if (m_includef0){
			//compute_f0_prior();
			if (m_f0<=0 || m_f0>=1){
	      		}else{
				if(!m_ardf0){}     	//Without ARD
				else              	//With ARD
		  			multifibres[idmult].m_f0_prior= log(double(m_f0));
	      		}
		}
	      	//compute_S0_prior(); m_S0_prior=0; so i don't do nothing, it is already 0
		m_prior_en = 0;
	      	//compute_prior();
	      	//m_prior_en=m_d_prior+m_S0_prior; is 0
	      	if(model==2)
			m_prior_en= m_prior_en+multifibres[idmult].m_dstd_prior;
	      	//if(m_rician) m_prior_en=m_prior_en+m_tau_prior; is 0
	      	if (m_includef0)
			m_prior_en=m_prior_en+multifibres[idmult].m_f0_prior;
	      	for(int fib=0;fib<nfib;fib++){
			m_prior_en=m_prior_en+ fibres[idmult*nfib+fib].m_prior_en;
	      	} 
		multifibres[idmult].m_prior_en = m_prior_en;
	}
      	//compute_iso_signal()
	for(int i=0; i<ndirs; i++){
		compute_iso_signal(&myisosig[i],&old, bvals[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK],m_d,m_dstd,model);
	}		
	__syncthreads();
      	//compute_likelihood()
	compute_likelihood(leader,idrest,m_d,m_S0,&m_likelihood_en,&m_f[0],mysig,myisosig,&mydata[idrest*MAXNDIRS_PER_THREAD],fsum,&reduction[0],m_f0,m_rician,m_tau,ndirs);
	__syncthreads();

	if(leader){
		multifibres[idmult].m_likelihood_en = m_likelihood_en;
	      	
		//compute_energy();	
		multifibres[idmult].m_energy = m_prior_en+m_likelihood_en;

	    	//initialise_props();
	      	multifibres[idmult].m_S0_prop=multifibres[idmult].m_S0/10.0; 
	      	multifibres[idmult].m_d_prop=m_d/10.0;
	      	multifibres[idmult].m_dstd_prop=m_dstd/10.0;
	      	multifibres[idmult].m_tau_prop=m_tau/2.0;
	      	multifibres[idmult].m_f0_prop=0.2;
	}

	for(int i=0; i<ndirs; i++){	
		isosignals[(idmult*NDIRECTIONS)+idrest+i*THREADS_BLOCK] = myisosig[i];
		for(int f=0;f<NFIBRES;f++){
			signals[(idmult*NDIRECTIONS*NFIBRES)+(f*NDIRECTIONS)+idrest+i*THREADS_BLOCK]=mysig[i*NFIBRES+f];
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
							const int 			nvox,
							const int 			model,
							const float 			fudgevalue,
							const bool 			m_include_f0,
							const bool 			m_ardf0,
							const bool 			can_use_ard, 
							const bool 			rician,
							const int 			updateproposalevery, 	//update every this number of iterations	
							const int 			iterations,		//num of iterations to do this time (maybe are a part of the total)
							const int 			current_iter,		//the number of the current iteration over the total iterations
							//INPUT-OUTPUT
							FibreGPU*			fibres,
							MultifibreGPU*		multifibres,
							double*				signals,
							double*				isosignals)
{	
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	if (id >= nvox*blockDim.x) { return; }

	int add=0;
	int add2=0;
	int add3=0;

	if(m_include_f0) add=1;
	if(rician) add2=1;
	if(model==2) add3=1;

	int idmult= blockIdx.x;
	int idrest= threadIdx.x;
	bool leader = (idrest==0);

	int ndirs = 1;
	if(idrest<(NDIRECTIONS%THREADS_BLOCK)) ndirs=MAXNDIRS_PER_THREAD;
	else if((NDIRECTIONS%THREADS_BLOCK)==0) ndirs=MAXNDIRS_PER_THREAD;
	else ndirs=(MAXNDIRS_PER_THREAD-1);

	if(idrest>=NDIRECTIONS)  { return; }

	__shared__ bool rejflag;				//max fibres 4
	__shared__ bool rejflag2;				//max 1024 directions / 64 = 16 MAXDIR
	__shared__ bool rejflag3;				//max shared used as base: 8KB of mydata + 
	__shared__ float fsum;					//64 doubles = 0.5 KB
	__shared__ double reduction[THREADS_BLOCK];		//23 floats + NFIBRES(4)*10 floats + 
	__shared__ float m_d;					//NPARAMS(16)*2 floats = 95 floats, 0.37 KB
	__shared__ float m_S0;					//3 bool + NFIBRES bool 
	__shared__ float m_f0;					//10 int + NFIBRES*6 int = 0.16KB 
	__shared__ float m_tau;					// MAX TOTAL 9.03 KB
	__shared__ float m_dstd;
	__shared__ float m_th[NFIBRES];
	__shared__ float m_ph[NFIBRES];
	__shared__ float m_f[NFIBRES];
	__shared__ float m_likelihood_en;

	__shared__ double tmp;

	double mysig[MAXNDIRS_PER_THREAD*NFIBRES];			
	double mysigold[MAXNDIRS_PER_THREAD*NFIBRES];			
	double myisosig[MAXNDIRS_PER_THREAD];			
	double myisosigold[MAXNDIRS_PER_THREAD];				
	__shared__ double mydata[MAXNDIRS_PER_THREAD*THREADS_BLOCK];	//save local memory, array max 1024 * 8 bytes = 8KB 

	__shared__ float old1;
	__shared__ float old2;
 		
	int idrand;
	int count_update;

	__shared__ float m_prior_en;
	__shared__ float m_old_prior_en;

	__shared__ float m_d_prior;
	__shared__ float m_S0_prior;
	__shared__ float m_f0_prior;
	__shared__ float m_tau_prior;
	__shared__ float m_dstd_prior;
	__shared__ float fm_prior_en[NFIBRES];		//m_prior_en for each fibre
	__shared__ float fm_old_prior_en;		//and old
	__shared__ float m_th_prior[NFIBRES];
	__shared__ float m_ph_prior[NFIBRES];
	__shared__ float m_f_prior[NFIBRES];

	__shared__ float m_energy;
	__shared__ float m_old_energy;

	__shared__ int m_d_acc;
	__shared__ int m_d_rej;
	__shared__ int m_S0_acc;
	__shared__ int m_S0_rej;
	__shared__ int m_f0_acc;
	__shared__ int m_f0_rej;
	__shared__ int m_tau_acc;
	__shared__ int m_tau_rej;	
	__shared__ int m_dstd_acc;
	__shared__ int m_dstd_rej;	

	__shared__ int m_th_acc[NFIBRES]; 
	__shared__ int m_th_rej[NFIBRES]; 
	__shared__ int m_ph_acc[NFIBRES]; 
	__shared__ int m_ph_rej[NFIBRES]; 
	__shared__ int m_f_acc[NFIBRES]; 
	__shared__ int m_f_rej[NFIBRES]; 
			
	__shared__ float m_d_prop;
	__shared__ float m_S0_prop;
	__shared__ float m_f0_prop;
	__shared__ float m_tau_prop;
	__shared__ float m_dstd_prop;
	__shared__ float m_th_prop[NFIBRES];
	__shared__ float m_ph_prop[NFIBRES];	
	__shared__ float m_f_prop[NFIBRES];	

	__shared__ bool m_lam_jump[NFIBRES];

	__shared__ float prerandN[NPARAMS];
	__shared__ float prerandU[NPARAMS];

	float angtmp[MAXNDIRS_PER_THREAD*NFIBRES];			//this is for pre compute signal
	float old_angtmp[MAXNDIRS_PER_THREAD*NFIBRES];
	float cos_alpha_minus_theta;
	float cos_alpha_plus_theta;
	
	if (leader){
		idrand = idmult*iterations*NPARAMS;
		count_update = current_iter;
		m_prior_en=multifibres[idmult].m_prior_en;

		if(model==2){
			m_dstd_acc=multifibres[idmult].m_dstd_acc;
			m_dstd_rej=multifibres[idmult].m_dstd_rej;
			m_dstd_prior=multifibres[idmult].m_dstd_prior;
			m_dstd_prop=multifibres[idmult].m_dstd_prop;
			m_dstd=multifibres[idmult].m_dstd;
		}else{
			m_dstd_acc=0;
			m_dstd_rej=0;
			m_dstd_prior=0;
			m_dstd_prop=0;
			m_dstd=0;
		}
	
		m_d=multifibres[idmult].m_d;
		m_energy=multifibres[idmult].m_energy;
		m_d_prop=multifibres[idmult].m_d_prop;
		m_d_prior=multifibres[idmult].m_d_prior;
		m_S0_prior=multifibres[idmult].m_S0_prior;
		m_S0=multifibres[idmult].m_S0;
		m_likelihood_en=multifibres[idmult].m_likelihood_en;
		m_d_acc=multifibres[idmult].m_d_acc;
		m_d_rej=multifibres[idmult].m_d_rej;
		m_S0_acc=multifibres[idmult].m_S0_acc;
		m_S0_rej=multifibres[idmult].m_S0_rej;
		m_S0_prop=multifibres[idmult].m_S0_prop;

		if(m_include_f0){
			m_f0_acc=multifibres[idmult].m_f0_acc;
			m_f0_rej=multifibres[idmult].m_f0_rej;
			m_f0_prop=multifibres[idmult].m_f0_prop;
			m_f0_prior=multifibres[idmult].m_f0_prior;
			m_f0=multifibres[idmult].m_f0;
		}else{ 
			m_f0_acc=0;
			m_f0_rej=0;
			m_f0_prop=0;
			m_f0_prior=0;
			m_f0=0;
		}
				
		
		if(rician){
			m_tau_acc=multifibres[idmult].m_tau_acc;
			m_tau_rej=multifibres[idmult].m_tau_rej;
			m_tau_prop=multifibres[idmult].m_tau_prop;
			m_tau_prior=multifibres[idmult].m_tau_prior;
			m_tau=multifibres[idmult].m_tau;	
		}else{ 
			m_tau_acc=0;
			m_tau_rej=0;
			m_tau_prop=0;
			m_tau_prior=0;
			m_tau=0;
		}


		for(int f=0;f<NFIBRES;f++){
			m_th[f]=fibres[(idmult*NFIBRES)+f].m_th;
			m_ph[f]=fibres[(idmult*NFIBRES)+f].m_ph;
			m_f[f]=fibres[(idmult*NFIBRES)+f].m_f;
	
			m_th_acc[f]=fibres[(idmult*NFIBRES)+f].m_th_acc;
			m_th_rej[f]=fibres[(idmult*NFIBRES)+f].m_th_rej;
			m_ph_acc[f]=fibres[(idmult*NFIBRES)+f].m_ph_acc;
			m_ph_rej[f]=fibres[(idmult*NFIBRES)+f].m_ph_rej;
			m_f_acc[f]=fibres[(idmult*NFIBRES)+f].m_f_acc;
			m_f_rej[f]=fibres[(idmult*NFIBRES)+f].m_f_rej;

			fm_prior_en[f]=fibres[(idmult*NFIBRES)+f].m_prior_en;
			m_th_prior[f]=fibres[(idmult*NFIBRES)+f].m_th_prior;
			m_ph_prior[f]=fibres[(idmult*NFIBRES)+f].m_ph_prior;
			m_f_prior[f]=fibres[(idmult*NFIBRES)+f].m_f_prior;

			m_th_prop[f]=fibres[(idmult*NFIBRES)+f].m_th_prop;
			m_ph_prop[f]=fibres[(idmult*NFIBRES)+f].m_ph_prop;
			m_f_prop[f]=fibres[(idmult*NFIBRES)+f].m_f_prop;

			m_lam_jump[f]=fibres[(idmult*NFIBRES)+f].m_lam_jump;		
		}

	}

	for(int i=0; i<ndirs; i++){	
		myisosig[i] = isosignals[(idmult*NDIRECTIONS)+idrest+i*THREADS_BLOCK];
		mydata[idrest*MAXNDIRS_PER_THREAD+i] = datam[(idmult*NDIRECTIONS)+idrest+i*THREADS_BLOCK];
	}

	for(int f=0;f<NFIBRES;f++){
		for(int i=0; i<ndirs; i++){	
			mysig[i*NFIBRES+f]= signals[(idmult*NDIRECTIONS*NFIBRES)+(f*NDIRECTIONS)+idrest+i*THREADS_BLOCK];	
		}
	}

	__syncthreads();
		
		//compute_signal_pre
		for(int i=0; i<ndirs; i++){
			for(int f=0;f<NFIBRES;f++){				
				cos_alpha_minus_theta=cos(double(alpha[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]-m_th[f]));   
			  	cos_alpha_plus_theta=cos(double(alpha[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]+m_th[f]));
		     	   	angtmp[i*NFIBRES+f]= __dadd_rn(__ddiv_rn (cos(double(m_ph[f]-beta[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]))*__dadd_rn(cos_alpha_minus_theta,-cos_alpha_plus_theta),2) ,__ddiv_rn (__dadd_rn(cos_alpha_minus_theta,cos_alpha_plus_theta),2));
		 	   	angtmp[i*NFIBRES+f]=angtmp[i*NFIBRES+f]*angtmp[i*NFIBRES+f];
			}
		}

		for (int niter=0; niter<iterations; niter++){
		//code jump()

			//prefetching randoms
			if (leader){
				if(m_include_f0){
					prerandN[0]=randomsN[idrand+(niter*NPARAMS)];
					prerandU[0]=randomsU[idrand+(niter*NPARAMS)];
				}
				if(rician){
					prerandN[add]=randomsN[idrand+(niter*NPARAMS)+add];
					prerandU[add]=randomsU[idrand+(niter*NPARAMS)+add];
				}
				prerandN[add+add2]=randomsN[idrand+(niter*NPARAMS)+add+add2];
				prerandU[add+add2]=randomsU[idrand+(niter*NPARAMS)+add+add2];

				if(model==2){
					prerandN[1+add+add2]=randomsN[idrand+(niter*NPARAMS)+1+add+add2];
					prerandU[1+add+add2]=randomsU[idrand+(niter*NPARAMS)+1+add+add2];
				}

				prerandN[1+add+add2+add3]=randomsN[idrand+(niter*NPARAMS)+1+add+add2+add3];
				prerandU[1+add+add2+add3]=randomsU[idrand+(niter*NPARAMS)+1+add+add2+add3];

				for(int f=0;f<NFIBRES;f++){
					prerandN[2+(f*3)+add+add2+add3]=randomsN[idrand+(niter*NPARAMS)+2+(f*3)+add+add2+add3];
					prerandN[3+(f*3)+add+add2+add3]=randomsN[idrand+(niter*NPARAMS)+3+(f*3)+add+add2+add3];
					prerandN[4+(f*3)+add+add2+add3]=randomsN[idrand+(niter*NPARAMS)+4+(f*3)+add+add2+add3];

					prerandU[2+(f*3)+add+add2+add3]=randomsU[idrand+(niter*NPARAMS)+2+(f*3)+add+add2+add3];
					prerandU[3+(f*3)+add+add2+add3]=randomsU[idrand+(niter*NPARAMS)+3+(f*3)+add+add2+add3];
					prerandU[4+(f*3)+add+add2+add3]=randomsU[idrand+(niter*NPARAMS)+4+(f*3)+add+add2+add3];		
				}
			}

////////////////////////////////////////////////////////////////// F0

			if(m_include_f0){
				if (leader){
					//propose_f0
					old1=m_f0;
					m_f0+=prerandN[0]*m_f0_prop;
				
					//compute_f0_prior()     
					old2=m_f0_prior;
	      				if(m_f0<=0 || m_f0 >=1){ 
						rejflag=true;
					}else{ 	
						rejflag=false;
						if(!m_ardf0){
							m_f0_prior=0;
	      					}else{
							m_f0_prior=log(double(m_f0));
						}
					}
				
					//for likehood and reject_f_sum
					fsum=m_f0;

					for(int g=0;g<NFIBRES;g++){  
						fsum+=m_f[g];
					}

					//compute_prior()
					compute_prior(&m_prior_en,&m_old_prior_en,m_d_prior,m_S0_prior,&fm_prior_en[0],m_f0_prior,m_tau_prior,m_dstd_prior);
				
					//reject_f_sum()
					rejflag2=(fsum>1);
				}

				__syncthreads();

				//compute_likelihood()
				compute_likelihood(leader,idrest,m_d,m_S0,&m_likelihood_en,&m_f[0],mysig,myisosig,&mydata[idrest*MAXNDIRS_PER_THREAD],fsum,&reduction[0],m_f0,rician,m_tau,ndirs);

				__syncthreads();

				if (leader){
					//compute_energy()
					m_old_energy=m_energy;
	      				m_energy=m_prior_en+m_likelihood_en;

					//test_energy
					tmp=exp(double(m_old_energy-m_energy));
					rejflag3=(tmp>prerandU[0]);	

					if(!rejflag){
						if(!rejflag2){
							if(rejflag3){
								//accept_f0()
		  						m_f0_acc++;   
							}else{
								//reject_F0()
								m_f0=old1;
								m_f0_prior=old2;
		      						m_prior_en=m_old_prior_en;
		      						m_f0_rej++;


								//restore_energy()
	      							m_energy=m_old_energy;
							}
						}else{ 
							//reject_F0()
							m_f0=old1;
							m_f0_prior=old2;
		      					m_prior_en=m_old_prior_en;
		      					m_f0_rej++;

							//restore_energy()
		      					m_energy=m_old_energy;
						}
					}else{ 
						//reject_F0()
						m_f0=old1;
						m_f0_prior=old2;
		      				m_prior_en=m_old_prior_en;
		      				m_f0_rej++;

						//restore_energy()
		      				m_energy=m_old_energy;
					}
				}
			}

////////////////////////////////////////////////////////////////// TAU
			if(rician){
				if (leader){
					//propose_tau
					old1=m_tau;
					m_tau+=prerandN[add]*m_tau_prop;
			
					//compute_tau_prior()     
					old2=m_tau_prior;
	      				if(m_tau<=0){ rejflag=true;
					}else{ 	
						rejflag=false;
						m_tau_prior=0;
					}

					//for likehood and reject_f_sum
					fsum=m_f0;
					for(int g=0;g<NFIBRES;g++){  
						fsum+=m_f[g];
					}

					//compute_prior()
					compute_prior(&m_prior_en,&m_old_prior_en,m_d_prior,m_S0_prior,&fm_prior_en[0],m_f0_prior,m_tau_prior,m_dstd_prior);
				}

				__syncthreads();

				//compute_likelihood()
				compute_likelihood(leader,idrest,m_d,m_S0,&m_likelihood_en,&m_f[0],mysig,myisosig,&mydata[idrest*MAXNDIRS_PER_THREAD],fsum,&reduction[0],m_f0,rician,m_tau,ndirs);

				__syncthreads();

				if (leader){
					//compute_energy()
					m_old_energy=m_energy;
	      				m_energy=m_prior_en+m_likelihood_en;

					//test_energy
					tmp=exp(double(m_old_energy-m_energy));
					rejflag2=(tmp>prerandU[add]);	
				
					if(!rejflag){
						if(rejflag2){
							//accept_tau()
		  					m_tau_acc++;   
						}else{ 
							//reject_tau()
							m_tau=old1;
							m_tau_prior=old2;
		      					m_prior_en=m_old_prior_en;
		      					m_tau_rej++;
							//restore_energy()
		      					m_energy=m_old_energy;
						}
					}else{ 
							//reject_tau()
							m_tau=old1;
							m_tau_prior=old2;
		      					m_prior_en=m_old_prior_en;
		      					m_tau_rej++;

							//restore_energy()
		      					m_energy=m_old_energy;
					}
				}
			}

////////////////////////////////////////////////////////////////// D

			if (leader){
				//propose_d()
				old1=m_d;	
				m_d+=prerandN[add+add2]*m_d_prop;
				
				//compute_d_prior_f0()      
				old2=m_d_prior;
					
      				if(m_d<0 || m_d > 0.008){
					rejflag=true;
				}else{ 
					m_d_prior=0;
					rejflag=false;
      				}
			}

			__syncthreads();
				
				//compute_signal()
				for(int i=0; i<ndirs; i++){
					for(int f=0;f<NFIBRES;f++){
						compute_signal(&mysig[i*NFIBRES+f],&mysigold[i*NFIBRES+f],bvals[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK],m_d,model,m_dstd,angtmp[i*NFIBRES+f]);
					
					}
				}
				//compute_iso_signal()
				for(int i=0; i<ndirs; i++){
					compute_iso_signal(&myisosig[i],&myisosigold[i], bvals[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK],m_d,m_dstd,model);
				}				

				if (leader){

					//for likehood and reject_f_sum
					fsum=m_f0;

					for(int g=0;g<NFIBRES;g++){  
						fsum+=m_f[g];
					}

					//compute_prior()
					compute_prior(&m_prior_en,&m_old_prior_en,m_d_prior,m_S0_prior,&fm_prior_en[0],m_f0_prior,m_tau_prior,m_dstd_prior);
				}

				__syncthreads();
				
				//compute_likelihood()
				compute_likelihood(leader,idrest,m_d,m_S0,&m_likelihood_en,&m_f[0],mysig,myisosig,&mydata[idrest*MAXNDIRS_PER_THREAD],fsum,&reduction[0],m_f0,rician,m_tau,ndirs);	
				
				__syncthreads();
				
				if (leader){
					//compute_energy()
					m_old_energy=m_energy;
      					m_energy=m_prior_en+m_likelihood_en;
					//test_energy
					tmp=exp(double(m_old_energy-m_energy));
					rejflag2=(tmp>prerandU[add+add2]);
				
				}
			
				__syncthreads();

       			if(!rejflag){
				if(rejflag2){
					//accept_d()
	  				if (leader) m_d_acc++;   
				}else{
					if (leader){
						//reject_d()
						m_d=old1;
						m_d_prior=old2;
      						m_prior_en=m_old_prior_en;
						m_d_rej++;
						//restore_energy()
      						m_energy=m_old_energy;
					}
      						
					for(int i=0; i<ndirs; i++){
						for(int f=0;f<NFIBRES;f++){
							mysig[i*NFIBRES+f] = mysigold[i*NFIBRES+f];
						}	
				
						myisosig[i]=myisosigold[i];
					}
				}
        		}else{ 
				if (leader){
				//reject_d()
					m_d=old1;
					m_d_prior=old2;
      					m_prior_en=m_old_prior_en;
					m_d_rej++;
					//restore_energy()
      					m_energy=m_old_energy;
				}

      				for(int i=0; i<ndirs; i++){
					for(int f=0;f<NFIBRES;f++){
						mysig[i*NFIBRES+f] = mysigold[i*NFIBRES+f];
					}	
				
					myisosig[i]=myisosigold[i];
				}
      				
        		}

////////////////////////////////////////////////////////////////// D_STD

			if(model==2){
				if (leader){	
					//propose_d_std
					old1=m_dstd;
					m_dstd+=prerandN[1+add+add2]*m_dstd_prop;

					//compute_d_std_prior()     
					old2=m_dstd_prior;
	      				if(m_dstd<=0 || m_dstd > 0.01){ 
						rejflag=true;
					}else{ 	
						rejflag=false;
						m_dstd_prior=log(double(m_dstd));
					}

				}

				__syncthreads();

				//compute_signal()
				for(int i=0; i<ndirs; i++){
					for(int f=0;f<NFIBRES;f++){
						compute_signal(&mysig[i*NFIBRES+f],&mysigold[i*NFIBRES+f],bvals[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK],m_d,model,m_dstd,angtmp[i*NFIBRES+f]);
					}
				}

				//compute_iso_signal()
				for(int i=0; i<ndirs; i++){
					compute_iso_signal(&myisosig[i],&myisosigold[i], bvals[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK],m_d,m_dstd,model);
				}

				if (leader){
				
					//for likehood and reject_f_sum
					fsum=m_f0;
					for(int g=0;g<NFIBRES;g++){  
						fsum+=m_f[g];
					}

					//compute_prior()
					compute_prior(&m_prior_en,&m_old_prior_en,m_d_prior,m_S0_prior,&fm_prior_en[0],m_f0_prior,m_tau_prior,m_dstd_prior);

				}

				__syncthreads();

				//compute_likelihood()
				compute_likelihood(leader,idrest,m_d,m_S0,&m_likelihood_en,&m_f[0],mysig,myisosig,&mydata[idrest*MAXNDIRS_PER_THREAD],fsum,&reduction[0],m_f0,rician,m_tau,ndirs);

				__syncthreads();

				if (leader){
					//compute_energy()
					m_old_energy=m_energy;
	      				m_energy=m_prior_en+m_likelihood_en;

					//test_energy
					tmp=exp(double(m_old_energy-m_energy));
					rejflag2=(tmp>prerandU[1+add+add2]);					
				}

				__syncthreads();
				
				if(!rejflag){
					if(rejflag2){
						//accept_dstd()
		  				if (leader) m_dstd_acc++;   
					}else{ 
						if (leader){
							//reject_dstd()
							m_dstd=old1;
							m_dstd_prior=old2;
		      					m_prior_en=m_old_prior_en;
		      					m_dstd_rej++;

							//restore_energy()
		      					m_energy=m_old_energy;
						}

						for(int i=0; i<ndirs; i++){
							for(int f=0;f<NFIBRES;f++){
								mysig[i*NFIBRES+f] = mysigold[i*NFIBRES+f];
							}	
				
							myisosig[i]=myisosigold[i];
						}
					}
				}else{ 
					if (leader){
						//reject_dstd()
						m_dstd=old1;
						m_dstd_prior=old2;
		      				m_prior_en=m_old_prior_en;
		      				m_dstd_rej++;

						//restore_energy()
		      				m_energy=m_old_energy;
					}

					for(int i=0; i<ndirs; i++){
						for(int f=0;f<NFIBRES;f++){
							mysig[i*NFIBRES+f] = mysigold[i*NFIBRES+f];
						}	
				
						myisosig[i]=myisosigold[i];
					}
				}
			}

////////////////////////////////////////////////////////////////// S0

			if (leader){
			
				//propose_S0()
				old1=m_S0;
				m_S0+=prerandN[1+add+add2+add3]*m_S0_prop;
				
				//compute_S0_prior()
				old2=m_S0_prior;
        			if(m_S0<0) rejflag=true;
        			else{    
					m_S0_prior=0;
	  				rejflag=false;
        			}
				
				//for likehood and reject_f_sum
				fsum=m_f0;
				for(int g=0;g<NFIBRES;g++){  
					fsum+=m_f[g];
				}

				//compute_prior()
				compute_prior(&m_prior_en,&m_old_prior_en,m_d_prior,m_S0_prior,&fm_prior_en[0],m_f0_prior,m_tau_prior,m_dstd_prior);
			}

			__syncthreads();

			//compute_likelihood()
			compute_likelihood(leader,idrest,m_d,m_S0,&m_likelihood_en,&m_f[0],mysig,myisosig,&mydata[idrest*MAXNDIRS_PER_THREAD],fsum,&reduction[0],m_f0,rician,m_tau,ndirs);

			__syncthreads();

			if (leader){
				//compute_energy()
				m_old_energy=m_energy;
      				m_energy=m_prior_en+m_likelihood_en;

				//test_energy
				tmp=exp(double(m_old_energy-m_energy));
				rejflag2=(tmp>prerandU[1+add+add2+add3]);

        			if(!rejflag){
					if(rejflag2){
						//accept_S0()
	  					m_S0_acc++;   
					}else{
						//reject_S0()
						m_S0=old1;
						m_S0_prior=old2;
      						m_prior_en=m_old_prior_en;
      						m_S0_rej++;

						//restore_energy()
      						m_energy=m_old_energy;
					}
        			}else{ 
					//reject_S0()
					m_S0=old1;
					m_S0_prior=old2;
	      				m_prior_en=m_old_prior_en;
	      				m_S0_rej++;

					//restore_energy()
      					m_energy=m_old_energy;
				}
        		}

////////////////////////////////////////////////////////////////////////////     TH

     			for(int fibra=0;fibra<NFIBRES;fibra++){  
				if (leader){ 

					//propose_th()
					old1=m_th[fibra];
					m_th[fibra]+=prerandN[2+(fibra*3)+add+add2+add3]*m_th_prop[fibra];
					
					//compute_th_prior()
					old2=m_th_prior[fibra];
      	   				if(m_th[fibra]==0){
						m_th_prior[fibra]=0;
		   			}else{
						m_th_prior[fibra]=-log(double(fabs(sin(double(m_th[fibra]))/2)));
	      	   			}
		  			//rejflag=false; /////////////////always false
		
					//compute_prior()
					fm_old_prior_en=fm_prior_en[fibra];
	      	   			fm_prior_en[fibra]=m_th_prior[fibra]+m_ph_prior[fibra]+m_f_prior[fibra];
					
				}

				__syncthreads();

				//compute_signal()
				//compute_signal_pre	
				for(int i=0; i<ndirs; i++){
					cos_alpha_minus_theta=cos(double(alpha[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]-m_th[fibra]));   
			  		cos_alpha_plus_theta=cos(double(alpha[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]+m_th[fibra]));
					old_angtmp[i*NFIBRES+fibra]=angtmp[i*NFIBRES+fibra];
		     	   		angtmp[i*NFIBRES+fibra]= __dadd_rn(__ddiv_rn (cos(double(m_ph[fibra]-beta[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]))*__dadd_rn(cos_alpha_minus_theta,-cos_alpha_plus_theta),2) ,__ddiv_rn (__dadd_rn(cos_alpha_minus_theta,cos_alpha_plus_theta),2));
		 	   		angtmp[i*NFIBRES+fibra]=angtmp[i*NFIBRES+fibra]*angtmp[i*NFIBRES+fibra];

					compute_signal(&mysig[i*NFIBRES+fibra],&mysigold[i*NFIBRES+fibra],bvals[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK],m_d,model,m_dstd,angtmp[i*NFIBRES+fibra]);
				}

				if (leader){
					
					//for likehood and reject_f_sum
					fsum=m_f0;
					for(int g=0;g<NFIBRES;g++){  
						fsum+=m_f[g];
					}

					//compute_prior()
					compute_prior(&m_prior_en,&m_old_prior_en,m_d_prior,m_S0_prior,&fm_prior_en[0],m_f0_prior,m_tau_prior,m_dstd_prior);
				
				}

				__syncthreads();

				//compute_likelihood()
				compute_likelihood(leader,idrest,m_d,m_S0,&m_likelihood_en,&m_f[0],mysig,myisosig,&mydata[idrest*MAXNDIRS_PER_THREAD],fsum,&reduction[0],m_f0,rician,m_tau,ndirs);

				__syncthreads();

				if (leader){ 
					//compute_energy()
					m_old_energy=m_energy;
	      				m_energy=m_prior_en+m_likelihood_en;

					//test_energy
					tmp=exp(double(m_old_energy-m_energy));
					rejflag2=(tmp>prerandU[2+(fibra*3)+add+add2+add3]);
					
				}

				__syncthreads();
			
					if(rejflag2){
						//accept_th
		  				if (leader) m_th_acc[fibra]++;   
					}else{

						if (leader){
						//reject_th()
							m_th[fibra]=old1;
							m_th_prior[fibra]=old2;
      							fm_prior_en[fibra]=fm_old_prior_en;
							m_th_rej[fibra]++;						
						
							//restore_energy()
							m_prior_en=m_old_prior_en;
	      						m_energy=m_old_energy;

						}

						//compute_signal_pre undo
						for(int i=0; i<ndirs; i++){
							angtmp[i*NFIBRES+fibra]=old_angtmp[i*NFIBRES+fibra];
							mysig[i*NFIBRES+fibra] = mysigold[i*NFIBRES+fibra];						
      						}
					}
				__syncthreads();
		
///////////////////////////////////////     PH

				if (leader){

					//propose_ph()
					old1=m_ph[fibra];
					m_ph[fibra]+=prerandN[3+(fibra*3)+add+add2+add3]*m_ph_prop[fibra];

					//compute_ph_prior()
					old2=m_ph_prior[fibra];
      					m_ph_prior[fibra]=0;
      					//rejflag=false;

					//compute_prior()
					fm_old_prior_en=fm_prior_en[fibra];
      	   				fm_prior_en[fibra]=m_th_prior[fibra]+m_ph_prior[fibra]+m_f_prior[fibra];

				}

				__syncthreads();
				
				//compute_signal()
				//compute_signal_pre
				for(int i=0; i<ndirs; i++){
					cos_alpha_minus_theta=cos(double(alpha[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]-m_th[fibra]));   
			  		cos_alpha_plus_theta=cos(double(alpha[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]+m_th[fibra]));
					old_angtmp[i*NFIBRES+fibra]=angtmp[i*NFIBRES+fibra];
		     	   		angtmp[i*NFIBRES+fibra]= __dadd_rn(__ddiv_rn (cos(double(m_ph[fibra]-beta[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]))*__dadd_rn(cos_alpha_minus_theta,-cos_alpha_plus_theta),2) ,__ddiv_rn (__dadd_rn(cos_alpha_minus_theta,cos_alpha_plus_theta),2));
		 	   		angtmp[i*NFIBRES+fibra]=angtmp[i*NFIBRES+fibra]*angtmp[i*NFIBRES+fibra];
	
					compute_signal(&mysig[i*NFIBRES+fibra],&mysigold[i*NFIBRES+fibra],bvals[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK],m_d,model,m_dstd,angtmp[i*NFIBRES+fibra]);
				}

				if (leader){
					//for likehood and reject_f_sum
					fsum=m_f0;
					for(int g=0;g<NFIBRES;g++){  
						fsum+=m_f[g];
					}

					//compute_prior()
					compute_prior(&m_prior_en,&m_old_prior_en,m_d_prior,m_S0_prior,&fm_prior_en[0],m_f0_prior,m_tau_prior,m_dstd_prior);
				}

				__syncthreads();

				//compute_likelihood()
				compute_likelihood(leader,idrest,m_d,m_S0,&m_likelihood_en,&m_f[0],mysig,myisosig,&mydata[idrest*MAXNDIRS_PER_THREAD],fsum,&reduction[0],m_f0,rician,m_tau,ndirs);


				__syncthreads();

				if (leader){
					//compute_energy()
					m_old_energy=m_energy;
	      				m_energy=m_prior_en+m_likelihood_en;

					//test_energy
					tmp=exp(double(m_old_energy-m_energy));
					rejflag2=(tmp>prerandU[3+(fibra*3)+add+add2+add3]);
				}

					__syncthreads();

				//if(!rejflag){
		        
					if(rejflag2){
						//accept_ph()
		  				if (leader) m_ph_acc[fibra]++;   
					}else{
						if (leader){
						//reject_ph()
							m_ph[fibra]=old1;
							m_ph_prior[fibra]=old2;
      							fm_prior_en[fibra]=fm_old_prior_en;
							m_ph_rej[fibra]++;						
						
							//restore_energy()
							m_prior_en=m_old_prior_en;
	      						m_energy=m_old_energy;
						}
						//compute_signal_pre undo
						for(int i=0; i<ndirs; i++){
							angtmp[i*NFIBRES+fibra]=old_angtmp[i*NFIBRES+fibra];

							mysig[i*NFIBRES+fibra] = mysigold[i*NFIBRES+fibra];	
						}
					}

				__syncthreads();
			

////////////////////////////////////////////             F

				if (leader){
					//propose_f()
					old1=m_f[fibra];
					m_f[fibra]+=prerandN[4+(fibra*3)+add+add2+add3]*m_f_prop[fibra];

	     				//compute_f_prior()
	        			old2=m_f_prior[fibra];
					if (m_f[fibra]<=0 || m_f[fibra]>=1) rejflag=true;
	        			else{
		      				if(!can_use_ard ){
		  					m_f_prior[fibra]=0;
						}else{
		  					if(m_lam_jump[fibra]){
								m_f_prior[fibra]=log(double(m_f[fibra]));
							}else{
		    						m_f_prior[fibra]=0;
		  					}
						}
						m_f_prior[fibra]=fudgevalue*m_f_prior[fibra];
						rejflag=false;
	      				}

					//compute_prior()
					fm_old_prior_en=fm_prior_en[fibra];
      	   				fm_prior_en[fibra]=m_th_prior[fibra]+m_ph_prior[fibra]+m_f_prior[fibra];
					
					
					//for likehood and reject_f_sum
					fsum=m_f0;
					for(int g=0;g<NFIBRES;g++){  
						fsum+=m_f[g];
					}
					//reject_f_sum()
					rejflag2=(fsum>1);

					//compute_prior()
					compute_prior(&m_prior_en,&m_old_prior_en,m_d_prior,m_S0_prior,&fm_prior_en[0],m_f0_prior,m_tau_prior,m_dstd_prior);	
				}

				__syncthreads();

				//compute_likelihood()
				compute_likelihood(leader,idrest,m_d,m_S0,&m_likelihood_en,&m_f[0],mysig,myisosig,&mydata[idrest*MAXNDIRS_PER_THREAD],fsum,&reduction[0],m_f0,rician,m_tau,ndirs);	

				__syncthreads();


				if (leader){
					//compute_energy()
					m_old_energy=m_energy;
		      			m_energy=m_prior_en+m_likelihood_en;

					//test_energy
					tmp=exp(double(m_old_energy-m_energy));
					rejflag3=(tmp>prerandU[4+(fibra*3)+add+add2+add3]);	

		      			if(!rejflag){
						if(!rejflag2){
							if(rejflag3){
								//accept_f()
			  					m_f_acc[fibra]++;   
							}else{
								//reject_f()
								m_f[fibra]=old1;
								m_f_prior[fibra]=old2;
	      							fm_prior_en[fibra]=fm_old_prior_en;
	      							m_f_rej[fibra]++;				
						
								//restore_energy()
								m_prior_en=m_old_prior_en;
		      						m_energy=m_old_energy;
							}
						}else{ 
							//reject_f()
							m_f[fibra]=old1;
							m_f_prior[fibra]=old2;
	      						fm_prior_en[fibra]=fm_old_prior_en;
	      						m_f_rej[fibra]++;

							//restore_energy()
							m_prior_en=m_old_prior_en;
		      					m_energy=m_old_energy;	
						}
					}else{
						//reject_f()
						m_f[fibra]=old1;
						m_f_prior[fibra]=old2;
		      				fm_prior_en[fibra]=fm_old_prior_en;
		      				m_f_rej[fibra]++;

						//restore_energy()
						m_prior_en=m_old_prior_en;
		      				m_energy=m_old_energy;
					}
				}
				__syncthreads();
			

        		}//end while NFIBRES

        		count_update++;

        		if(((count_update%updateproposalevery)==0)&&leader){
				//m_multifibre.update_proposals();
				m_d_prop*=sqrt(float(m_d_acc+1)/float(m_d_rej+1));
				m_d_prop=min(m_d_prop,maxfloat);

				if(rician){
					m_tau_prop*=sqrt(float(m_tau_acc+1)/float(m_tau_rej+1));
					m_tau_prop=min(m_tau_prop,maxfloat);
					m_tau_acc=0; 
					m_tau_rej=0;	
				}

				if(m_include_f0){
					m_f0_prop*=sqrt(float(m_f0_acc+1)/float(m_f0_rej+1));
					m_f0_prop=min(m_f0_prop,maxfloat);
					m_f0_acc=0; 
					m_f0_rej=0;	
				}	

				if(model==2){
					m_dstd_prop*=sqrt(float(m_dstd_acc+1)/float(m_dstd_rej+1));
					m_dstd_prop=min(m_dstd_prop,maxfloat);
					m_dstd_acc=0; 
					m_dstd_rej=0;	
				}

				m_S0_prop*=sqrt(float(m_S0_acc+1)/float(m_S0_rej+1));
				m_S0_prop=min(m_S0_prop,maxfloat);
				m_d_acc=0; 
				m_d_rej=0;
				m_S0_acc=0; 
				m_S0_rej=0;
				for(int f=0; f<NFIBRES;f++){
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
			multifibres[idmult].m_prior_en=m_prior_en;
			multifibres[idmult].m_d=m_d;
			multifibres[idmult].m_energy=m_energy;		
			multifibres[idmult].m_d_prop=m_d_prop;
			multifibres[idmult].m_d_prior=m_d_prior;

			for(int f=0;f<NFIBRES;f++){
				fibres[(idmult*NFIBRES)+f].m_th=m_th[f];
				fibres[(idmult*NFIBRES)+f].m_ph=m_ph[f];
				fibres[(idmult*NFIBRES)+f].m_f=m_f[f];
			}

			multifibres[idmult].m_S0_prior=m_S0_prior;
			multifibres[idmult].m_S0=m_S0;
			multifibres[idmult].m_likelihood_en=m_likelihood_en;
			multifibres[idmult].m_d_acc=m_d_acc;
			multifibres[idmult].m_d_rej=m_d_rej;
			multifibres[idmult].m_S0_acc=m_S0_acc;
			multifibres[idmult].m_S0_rej=m_S0_rej;
			multifibres[idmult].m_S0_prop=m_S0_prop;

			if(m_include_f0){
				multifibres[idmult].m_f0_prior=m_f0_prior;
				multifibres[idmult].m_f0=m_f0;
				multifibres[idmult].m_f0_acc=m_f0_acc;
				multifibres[idmult].m_f0_rej=m_f0_rej;
				multifibres[idmult].m_f0_prop=m_f0_prop;
			}
			if(rician){
				multifibres[idmult].m_tau_prior=m_tau_prior;
				multifibres[idmult].m_tau=m_tau;
				multifibres[idmult].m_tau_acc=m_tau_acc;
				multifibres[idmult].m_tau_rej=m_tau_rej;
				multifibres[idmult].m_tau_prop=m_tau_prop;
			}
			if(model==2){
				multifibres[idmult].m_dstd_prior=m_dstd_prior;
				multifibres[idmult].m_dstd=m_dstd;
				multifibres[idmult].m_dstd_acc=m_dstd_acc;
				multifibres[idmult].m_dstd_rej=m_dstd_rej;
				multifibres[idmult].m_dstd_prop=m_dstd_prop;
			}
	
			for(int f=0;f<NFIBRES;f++){
				fibres[(idmult*NFIBRES)+f].m_th_acc=m_th_acc[f];
				fibres[(idmult*NFIBRES)+f].m_th_rej=m_th_rej[f];
				fibres[(idmult*NFIBRES)+f].m_ph_acc=m_ph_acc[f];
				fibres[(idmult*NFIBRES)+f].m_ph_rej=m_ph_rej[f];
				fibres[(idmult*NFIBRES)+f].m_f_acc=m_f_acc[f];
				fibres[(idmult*NFIBRES)+f].m_f_rej=m_f_rej[f];

				fibres[(idmult*NFIBRES)+f].m_prior_en=fm_prior_en[f];
				fibres[(idmult*NFIBRES)+f].m_th_prior=m_th_prior[f];
				fibres[(idmult*NFIBRES)+f].m_ph_prior=m_ph_prior[f];
				fibres[(idmult*NFIBRES)+f].m_f_prior=m_f_prior[f];

				fibres[(idmult*NFIBRES)+f].m_th_prop=m_th_prop[f];
				fibres[(idmult*NFIBRES)+f].m_ph_prop=m_ph_prop[f];
				fibres[(idmult*NFIBRES)+f].m_f_prop=m_f_prop[f];

				fibres[(idmult*NFIBRES)+f].m_lam_jump=m_lam_jump[f];		
			}
		}

		for(int i=0; i<ndirs; i++){	
			isosignals[(idmult*NDIRECTIONS)+idrest+i*THREADS_BLOCK] = myisosig[i];

			for(int f=0;f<NFIBRES;f++){
				signals[(idmult*NDIRECTIONS*NFIBRES)+(f*NDIRECTIONS)+idrest+i*THREADS_BLOCK]=mysig[i*NFIBRES+f];
			}
		}
	
}

extern "C" __global__ void runmcmc_record_kernel(	//INPUT 
							const double*			datam,
							const double*			bvals,
							const double*			alpha,
							const double*			beta,
							const FibreGPU*			fibres,
							const MultifibreGPU*		multifibres,
							double*				signals,
							double*				isosignals,
							float*				randomsN,
							float*				randomsU,
							const int 			nvox,
							const int 			model,
							const float 			fudgevalue,
							const bool 			m_include_f0,
							const bool 			m_ardf0,
							const bool 			can_use_ard, 
							const bool 			rician,
							const int 			updateproposalevery, 	//update every this number of iterations	
							const int 			iterations,		//num of iterations to do this time (maybe are a part of the total)	
							const int 			current_iter,		//the number of the current iteration over the total iterations
							const int 			iters_burnin,		//iters in burin, we need it for continue the updates at the correct time. 
							const int 			record_every, 		//record every this number
							const int 			totalrecords,		//total number of records to do
							
							//OUTPUT
							int*			multirecords,		//id for each sample in each voxel
							float*				rf0,			//records od parameters
							float*				rtau,
							float*				rs0,
							float*				rd,
							float*				rdstd,
							float*				rth,
							float*				rph, 
							float*				rf,
							float*				rlikelihood_energy)

{
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	if (id >= nvox*blockDim.x) { return; }	

	int add=0;
	int add2=0;
	int add3=0;

	if(m_include_f0) add=1;
	if(rician) add2=1;
	if(model==2) add3=1;

	int idmult= blockIdx.x;
	int idrest= threadIdx.x;
	bool leader = ((id%THREADS_BLOCK)==0);	

	int ndirs = 1;
	if(idrest<(NDIRECTIONS%THREADS_BLOCK)) ndirs=MAXNDIRS_PER_THREAD;
	else if((NDIRECTIONS%THREADS_BLOCK)==0) ndirs=MAXNDIRS_PER_THREAD;
	else ndirs=(MAXNDIRS_PER_THREAD-1);

	if(idrest>=NDIRECTIONS) { return; }

	__shared__ bool rejflag;
	__shared__ bool rejflag2;
	__shared__ bool rejflag3;
	__shared__ float fsum;
	__shared__ double reduction[THREADS_BLOCK];
	__shared__ float m_d;
	__shared__ float m_S0;	
	__shared__ float m_f0;
	__shared__ float m_tau;
	__shared__ float m_dstd;
	__shared__ float m_th[NFIBRES];
	__shared__ float m_ph[NFIBRES];
	__shared__ float m_f[NFIBRES];
	__shared__ float m_likelihood_en;

	__shared__ double tmp;

	double mysig[MAXNDIRS_PER_THREAD*NFIBRES];			
	double mysigold[MAXNDIRS_PER_THREAD*NFIBRES];			
	double myisosig[MAXNDIRS_PER_THREAD];			
	double myisosigold[MAXNDIRS_PER_THREAD];		
	__shared__ double mydata[MAXNDIRS_PER_THREAD*THREADS_BLOCK];				

	__shared__ float old1;
	__shared__ float old2;
 		
	int idrand;
	int count_update;
	int recordcount;
	int sample;

	__shared__ float m_prior_en;
	__shared__ float m_old_prior_en;

	__shared__ float m_d_prior;	
	__shared__ float m_S0_prior;
	__shared__ float m_f0_prior;
	__shared__ float m_tau_prior;
	__shared__ float m_dstd_prior;

	__shared__ float fm_prior_en[NFIBRES];
	__shared__ float fm_old_prior_en;

	__shared__ float m_th_prior[NFIBRES];
	__shared__ float m_ph_prior[NFIBRES];
	__shared__ float m_f_prior[NFIBRES];

	__shared__ float m_energy;
	__shared__ float m_old_energy;

	__shared__ int m_d_acc;
	__shared__ int m_d_rej;
	
	__shared__ int m_S0_acc;
	__shared__ int m_S0_rej;

	__shared__ int m_f0_acc;
	__shared__ int m_f0_rej;
	__shared__ int m_tau_acc;
	__shared__ int m_tau_rej;
	__shared__ int m_dstd_acc;
	__shared__ int m_dstd_rej;

	__shared__ int m_th_acc[NFIBRES]; 
	__shared__ int m_th_rej[NFIBRES]; 
	__shared__ int m_ph_acc[NFIBRES]; 
	__shared__ int m_ph_rej[NFIBRES]; 
	__shared__ int m_f_acc[NFIBRES]; 
	__shared__ int m_f_rej[NFIBRES]; 
		
	__shared__ float m_d_prop;

	__shared__ float m_S0_prop;
	__shared__ float m_f0_prop;
	__shared__ float m_tau_prop;
	__shared__ float m_dstd_prop;

	__shared__ float m_th_prop[NFIBRES];
	__shared__ float m_ph_prop[NFIBRES];	
	__shared__ float m_f_prop[NFIBRES];	

	__shared__ bool m_lam_jump[NFIBRES];

	__shared__ float prerandN[NPARAMS];
	__shared__ float prerandU[NPARAMS];

	float angtmp[MAXNDIRS_PER_THREAD*NFIBRES];	//for pre compute sigmal
	float old_angtmp[MAXNDIRS_PER_THREAD*NFIBRES];
	float cos_alpha_minus_theta;
	float cos_alpha_plus_theta;

	if (leader){
		idrand = idmult*iterations*NPARAMS;

		count_update = current_iter + iters_burnin;	//count for updates
		recordcount= current_iter;	
		sample=1+(current_iter/record_every);		//the next number of sample.....the index start in 0

		m_prior_en=multifibres[idmult].m_prior_en;
		m_d=multifibres[idmult].m_d;
		m_energy=multifibres[idmult].m_energy;
		m_d_prop=multifibres[idmult].m_d_prop;
		m_d_prior=multifibres[idmult].m_d_prior;
		m_S0_prior=multifibres[idmult].m_S0_prior;
		m_S0=multifibres[idmult].m_S0;
		m_likelihood_en=multifibres[idmult].m_likelihood_en;
		m_d_acc=multifibres[idmult].m_d_acc;
		m_d_rej=multifibres[idmult].m_d_rej;
		m_S0_acc=multifibres[idmult].m_S0_acc;
		m_S0_rej=multifibres[idmult].m_S0_rej;
		m_S0_prop=multifibres[idmult].m_S0_prop;

		if(model==2){
			m_dstd_acc=multifibres[idmult].m_dstd_acc;
			m_dstd_rej=multifibres[idmult].m_dstd_rej;
			m_dstd_prior=multifibres[idmult].m_dstd_prior;
			m_dstd_prop=multifibres[idmult].m_dstd_prop;
			m_dstd=multifibres[idmult].m_dstd;
		}else{
			m_dstd_acc=0;
			m_dstd_rej=0;
			m_dstd_prior=0;
			m_dstd_prop=0;
			m_dstd=0;
		}
		
		if(m_include_f0){
			m_f0_acc=multifibres[idmult].m_f0_acc;
			m_f0_rej=multifibres[idmult].m_f0_rej;
			m_f0_prop=multifibres[idmult].m_f0_prop;
			m_f0_prior=multifibres[idmult].m_f0_prior;
			m_f0=multifibres[idmult].m_f0;
		}else{
			m_f0_acc=0;
			m_f0_rej=0;
			m_f0_prop=0;
			m_f0_prior=0;	
			m_f0=0;
		}		
		if(rician){
			m_tau_acc=multifibres[idmult].m_tau_acc;
			m_tau_rej=multifibres[idmult].m_tau_rej;
			m_tau_prop=multifibres[idmult].m_tau_prop;
			m_tau_prior=multifibres[idmult].m_tau_prior;
			m_tau=multifibres[idmult].m_tau;
		}else{ 
			m_tau_acc=0;
			m_tau_rej=0;
			m_tau_prop=0;
			m_tau_prior=0;
			m_tau=0;
		}

		for(int f=0;f<NFIBRES;f++){
			m_th[f]=fibres[(idmult*NFIBRES)+f].m_th;
			m_ph[f]=fibres[(idmult*NFIBRES)+f].m_ph;
			m_f[f]=fibres[(idmult*NFIBRES)+f].m_f;
			m_th_acc[f]=fibres[(idmult*NFIBRES)+f].m_th_acc;
			m_th_rej[f]=fibres[(idmult*NFIBRES)+f].m_th_rej;
			m_ph_acc[f]=fibres[(idmult*NFIBRES)+f].m_ph_acc;
			m_ph_rej[f]=fibres[(idmult*NFIBRES)+f].m_ph_rej;
			m_f_acc[f]=fibres[(idmult*NFIBRES)+f].m_f_acc;
			m_f_rej[f]=fibres[(idmult*NFIBRES)+f].m_f_rej;

			fm_prior_en[f]=fibres[(idmult*NFIBRES)+f].m_prior_en;
			m_th_prior[f]=fibres[(idmult*NFIBRES)+f].m_th_prior;
			m_ph_prior[f]=fibres[(idmult*NFIBRES)+f].m_ph_prior;
			m_f_prior[f]=fibres[(idmult*NFIBRES)+f].m_f_prior;

			m_th_prop[f]=fibres[(idmult*NFIBRES)+f].m_th_prop;
			m_ph_prop[f]=fibres[(idmult*NFIBRES)+f].m_ph_prop;
			m_f_prop[f]=fibres[(idmult*NFIBRES)+f].m_f_prop;

			m_lam_jump[f]=fibres[(idmult*NFIBRES)+f].m_lam_jump;		
		}
	}
			
	for(int i=0; i<ndirs; i++){	
		myisosig[i] = isosignals[(idmult*NDIRECTIONS)+idrest+i*THREADS_BLOCK];
		mydata[idrest*MAXNDIRS_PER_THREAD+i] = datam[(idmult*NDIRECTIONS)+idrest+i*THREADS_BLOCK];
	}

	for(int f=0;f<NFIBRES;f++){
		for(int i=0; i<ndirs; i++){	
			mysig[i*NFIBRES+f]= signals[(idmult*NDIRECTIONS*NFIBRES)+(f*NDIRECTIONS)+idrest+i*THREADS_BLOCK];
		}
	}

	__syncthreads();

	//compute_signal_pre
	for(int i=0; i<ndirs; i++){
		for(int f=0;f<NFIBRES;f++){				
			cos_alpha_minus_theta=cos(double(alpha[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]-m_th[f]));   
		  	cos_alpha_plus_theta=cos(double(alpha[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]+m_th[f]));
		   	angtmp[i*NFIBRES+f]= __dadd_rn(__ddiv_rn (cos(double(m_ph[f]-beta[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]))*__dadd_rn(cos_alpha_minus_theta,-cos_alpha_plus_theta),2) ,__ddiv_rn (__dadd_rn(cos_alpha_minus_theta,cos_alpha_plus_theta),2));
		 	angtmp[i*NFIBRES+f]=angtmp[i*NFIBRES+f]*angtmp[i*NFIBRES+f];
		}
	}

	for (int niter=0; niter<iterations; niter++){

	//code jump()

	//prefetching randoms
	if (leader){
		if(m_include_f0){
			prerandN[0]=randomsN[idrand+(niter*NPARAMS)];
			prerandU[0]=randomsU[idrand+(niter*NPARAMS)];
		}
		if(rician){
			prerandN[add]=randomsN[idrand+(niter*NPARAMS)+add];
			prerandU[add]=randomsU[idrand+(niter*NPARAMS)+add];
		}
		prerandN[add+add2]=randomsN[idrand+(niter*NPARAMS)+add+add2];
		prerandU[add+add2]=randomsU[idrand+(niter*NPARAMS)+add+add2];

		if(model==2){
			prerandN[1+add+add2]=randomsN[idrand+(niter*NPARAMS)+1+add+add2];
			prerandU[1+add+add2]=randomsU[idrand+(niter*NPARAMS)+1+add+add2];
		}

		prerandN[1+add+add2+add3]=randomsN[idrand+(niter*NPARAMS)+1+add+add2+add3];
		prerandU[1+add+add2+add3]=randomsU[idrand+(niter*NPARAMS)+1+add+add2+add3];

		for(int f=0;f<NFIBRES;f++){
			prerandN[2+(f*3)+add+add2+add3]=randomsN[idrand+(niter*NPARAMS)+2+(f*3)+add+add2+add3];
			prerandN[3+(f*3)+add+add2+add3]=randomsN[idrand+(niter*NPARAMS)+3+(f*3)+add+add2+add3];
			prerandN[4+(f*3)+add+add2+add3]=randomsN[idrand+(niter*NPARAMS)+4+(f*3)+add+add2+add3];

			prerandU[2+(f*3)+add+add2+add3]=randomsU[idrand+(niter*NPARAMS)+2+(f*3)+add+add2+add3];
			prerandU[3+(f*3)+add+add2+add3]=randomsU[idrand+(niter*NPARAMS)+3+(f*3)+add+add2+add3];
			prerandU[4+(f*3)+add+add2+add3]=randomsU[idrand+(niter*NPARAMS)+4+(f*3)+add+add2+add3];		
		}
	}

////////////////////////////////////////////////////////////////// F0

		if(m_include_f0){
			if (leader){
				//propose_f0
				old1=m_f0;
				m_f0+=prerandN[0]*m_f0_prop;
					
				//compute_f0_prior()     
				old2=m_f0_prior;
	      			if(m_f0<=0 || m_f0 >=1){ rejflag=true;
				}else{ 	
					rejflag=false;
					if(!m_ardf0){
						m_f0_prior=0;
	      				}else{
						m_f0_prior=log(double(m_f0));
					}
				}

				//for likehood and reject_f_sum
				fsum=m_f0;
				for(int g=0;g<NFIBRES;g++){  
					fsum+=m_f[g];
				}

				//compute_prior()
				compute_prior(&m_prior_en,&m_old_prior_en,m_d_prior,m_S0_prior,&fm_prior_en[0],m_f0_prior,m_tau_prior,m_dstd_prior);

				//reject_f_sum()
				rejflag2=(fsum>1);	
			
			}

			__syncthreads();

			//compute_likelihood()
			compute_likelihood(leader,idrest,m_d,m_S0,&m_likelihood_en,&m_f[0],mysig,myisosig, &mydata[idrest*MAXNDIRS_PER_THREAD],fsum,&reduction[0],m_f0,rician,m_tau,ndirs);

			__syncthreads();

			if (leader){
				//compute_energy()
				m_old_energy=m_energy;
	      			m_energy=m_prior_en+m_likelihood_en;

				//test_energy
				tmp=exp(double(m_old_energy-m_energy));
				rejflag3=(tmp>prerandU[0]);
				
				if(!rejflag){
					if(!rejflag2){
						if(rejflag3){
							//accept_f0()
		  					m_f0_acc++;   
						}else{
							//reject_F0()
							m_f0=old1;
							m_f0_prior=old2;
		      					m_prior_en=m_old_prior_en;
		      					m_f0_rej++;

							//restore_energy()
	      						m_energy=m_old_energy;
						}
					}else{ 
						//reject_F0()
						m_f0=old1;
						m_f0_prior=old2;
		      				m_prior_en=m_old_prior_en;
		      				m_f0_rej++;

						//restore_energy()
		      				m_energy=m_old_energy;
					}
				}else{ 
					//reject_F0()
					m_f0=old1;
					m_f0_prior=old2;
		      			m_prior_en=m_old_prior_en;
		      			m_f0_rej++;

					//restore_energy()
		      			m_energy=m_old_energy;
				}
			}
		}

////////////////////////////////////////////////////////////////// TAU

		if(rician){
			if (leader){
				//propose_tau
				old1=m_tau;
				m_tau+=prerandN[add]*m_tau_prop;
			
				//compute_tau_prior()     
				old2=m_tau_prior;
	      			if(m_tau<=0){ rejflag=true;
				}else{ 	
					rejflag=false;
					m_tau_prior=0;
				}
				
				//for likehood and reject_f_sum
				fsum=m_f0;
				for(int g=0;g<NFIBRES;g++){  
					fsum+=m_f[g];
				}

				//compute_prior()
				compute_prior(&m_prior_en,&m_old_prior_en,m_d_prior,m_S0_prior,&fm_prior_en[0],m_f0_prior,m_tau_prior,m_dstd_prior);
			}
			__syncthreads();

			//compute_likelihood()
			compute_likelihood(leader,idrest,m_d,m_S0,&m_likelihood_en,&m_f[0],mysig,myisosig, &mydata[idrest*MAXNDIRS_PER_THREAD],fsum,&reduction[0],m_f0,rician,m_tau,ndirs);

			__syncthreads();

			if (leader){
				//compute_energy()
				m_old_energy=m_energy;
	      			m_energy=m_prior_en+m_likelihood_en;

				//test_energy
				tmp=exp(double(m_old_energy-m_energy));
				rejflag2=(tmp>prerandU[add]);	
				
				if(!rejflag){
					if(rejflag2){
						//accept_tau()
		  				m_tau_acc++;   
					}else{ 
						//reject_tau()
						m_tau=old1;
						m_tau_prior=old2;
		      				m_prior_en=m_old_prior_en;
		      				m_tau_rej++;

						//restore_energy()
		      				m_energy=m_old_energy;
					}
				}else{ 
					//reject_tau()
					m_tau=old1;
					m_tau_prior=old2;
		      			m_prior_en=m_old_prior_en;
		      			m_tau_rej++;

					//restore_energy()
		      			m_energy=m_old_energy;
				}
			}
		}

////////////////////////////////////////////////////////////////// D

		if (leader){
			//propose_d()	
			old1=m_d;	
			m_d+=prerandN[add+add2]*m_d_prop;

			//compute_d_prior_f0()      
			old2=m_d_prior;
      			if(m_d<0 || m_d > 0.008 ){ 
				rejflag=true;
			}else{ 
				m_d_prior=0;
				rejflag=false;
      			}
		}

		__syncthreads();

		//compute_signal()
		for(int i=0; i<ndirs; i++){
			for(int f=0;f<NFIBRES;f++){
				compute_signal(&mysig[i*NFIBRES+f],&mysigold[i*NFIBRES+f],bvals[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK],m_d,model,m_dstd,angtmp[i*NFIBRES+f]);
			}
		}
						
		//compute_iso_signal()
		for(int i=0; i<ndirs; i++){
			compute_iso_signal(&myisosig[i],&myisosigold[i], bvals[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK],m_d,m_dstd,model);
		}

		if (leader){
			//for likehood and reject_f_sum
			fsum=m_f0;
			for(int g=0;g<NFIBRES;g++){  
				fsum+=m_f[g];
			}

			//compute_prior()
			compute_prior(&m_prior_en,&m_old_prior_en,m_d_prior,m_S0_prior,&fm_prior_en[0],m_f0_prior,m_tau_prior,m_dstd_prior);
		}

		__syncthreads();

		//compute_likelihood()
		compute_likelihood(leader,idrest,m_d,m_S0,&m_likelihood_en,&m_f[0],mysig,myisosig,&mydata[idrest*MAXNDIRS_PER_THREAD],fsum,&reduction[0],m_f0,rician,m_tau,ndirs);
				
		__syncthreads();

		if (leader){
			//compute_energy()
			m_old_energy=m_energy;
      			m_energy=m_prior_en+m_likelihood_en;

			//test_energy
			tmp=exp(double(m_old_energy-m_energy));
			rejflag2=(tmp>prerandU[add+add2]);
		}
		__syncthreads();

       		if(!rejflag){
			if(rejflag2){
				//accept_d()
	  			if (leader) m_d_acc++;   
			}else{
				if (leader){
					//reject_d()
					m_d=old1;
					m_d_prior=old2;
      					m_prior_en=m_old_prior_en;
					m_d_rej++;
					//restore_energy()
      					m_energy=m_old_energy;
				}
      						
				for(int i=0; i<ndirs; i++){
					for(int f=0;f<NFIBRES;f++){
						mysig[i*NFIBRES+f] = mysigold[i*NFIBRES+f];
					}	
					myisosig[i]=myisosigold[i];
				}
			}
        	}else{ 
			if (leader){
				//reject_d()
				m_d=old1;
				m_d_prior=old2;
      				m_prior_en=m_old_prior_en;
				m_d_rej++;
				//restore_energy()
      				m_energy=m_old_energy;
			}
      			for(int i=0; i<ndirs; i++){
				for(int f=0;f<NFIBRES;f++){
					mysig[i*NFIBRES+f] = mysigold[i*NFIBRES+f];
				}	
				myisosig[i]=myisosigold[i];
			}	
        	}

////////////////////////////////////////////////////////////////// D_STD

		if(model==2){
			if (leader){
				//propose_d_std
				old1=m_dstd;
				m_dstd+=prerandN[1+add+add2]*m_dstd_prop;
				
				//compute_d_std_prior()     
				old2=m_dstd_prior;
	      			if(m_dstd<=0 || m_dstd > 0.01){ 
					rejflag=true;
				}else{ 	
					rejflag=false;
					m_dstd_prior=log(double(m_dstd));
				}
			}

			__syncthreads();

			//compute_signal()
			for(int i=0; i<ndirs; i++){
				for(int f=0;f<NFIBRES;f++){
					compute_signal(&mysig[i*NFIBRES+f],&mysigold[i*NFIBRES+f],bvals[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK],m_d,model,m_dstd,angtmp[i*NFIBRES+f]);
				}
			}

			//compute_iso_signal()
			for(int i=0; i<ndirs; i++){
				compute_iso_signal(&myisosig[i],&myisosigold[i], bvals[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK],m_d,m_dstd,model);
			}

			if (leader){
				//for likehood and reject_f_sum
				fsum=m_f0;
				for(int g=0;g<NFIBRES;g++){  
					fsum+=m_f[g];
				}

				//compute_prior()
				compute_prior(&m_prior_en,&m_old_prior_en,m_d_prior,m_S0_prior,&fm_prior_en[0],m_f0_prior,m_tau_prior,m_dstd_prior);
			}

			__syncthreads();

			//compute_likelihood()
			compute_likelihood(leader,idrest,m_d,m_S0,&m_likelihood_en,&m_f[0],mysig,myisosig, &mydata[idrest*MAXNDIRS_PER_THREAD],fsum,&reduction[0],m_f0,rician,m_tau,ndirs);

			__syncthreads();

			if (leader){
				//compute_energy()
				m_old_energy=m_energy;
	      			m_energy=m_prior_en+m_likelihood_en;

				//test_energy
				tmp=exp(double(m_old_energy-m_energy));
				rejflag2=(tmp>prerandU[1+add+add2]);
			}
				
			if(!rejflag){
				if(rejflag2){
					//accept_dstd()
		  			if (leader) m_dstd_acc++;   
				}else{ 
					if (leader){
						//reject_dstd()
						m_dstd=old1;
						m_dstd_prior=old2;
		      				m_prior_en=m_old_prior_en;
		      				m_dstd_rej++;

						//restore_energy()
		      				m_energy=m_old_energy;
					}
					for(int i=0; i<ndirs; i++){
						for(int f=0;f<NFIBRES;f++){
							mysig[i*NFIBRES+f] = mysigold[i*NFIBRES+f];
						}	
						myisosig[i]=myisosigold[i];
					}
				}
			}else{ 
				if (leader){
					//reject_dstd()
					m_dstd=old1;
					m_dstd_prior=old2;
		      			m_prior_en=m_old_prior_en;
		      			m_dstd_rej++;

					//restore_energy()
		      			m_energy=m_old_energy;
				}
				for(int i=0; i<ndirs; i++){
					for(int f=0;f<NFIBRES;f++){
						mysig[i*NFIBRES+f] = mysigold[i*NFIBRES+f];
					}	
					myisosig[i]=myisosigold[i];
				}
			}
		}

////////////////////////////////////////////////////////////////// S0

		if (leader){
			//propose_S0()
			old1=m_S0;
			m_S0+=prerandN[1+add+add2+add3]*m_S0_prop;
				
			//compute_S0_prior()
			old2=m_S0_prior;
        		if(m_S0<0) rejflag=true;
        		else{    
				m_S0_prior=0;
	  			rejflag=false;
        		}
			//for likehood and reject_f_sum
			fsum=m_f0;
			for(int g=0;g<NFIBRES;g++){  
				fsum+=m_f[g];
			}

			//compute_prior()
			compute_prior(&m_prior_en,&m_old_prior_en,m_d_prior,m_S0_prior,&fm_prior_en[0],m_f0_prior,m_tau_prior,m_dstd_prior);
		}

		__syncthreads();

		//compute_likelihood()
		compute_likelihood(leader,idrest,m_d,m_S0,&m_likelihood_en,&m_f[0],mysig,myisosig,&mydata[idrest*MAXNDIRS_PER_THREAD],fsum,&reduction[0],m_f0,rician,m_tau,ndirs);

		__syncthreads();

		if (leader){
			//compute_energy()
			m_old_energy=m_energy;
      			m_energy=m_prior_en+m_likelihood_en;

			//test_energy
			tmp=exp(double(m_old_energy-m_energy));
			rejflag2=(tmp>prerandU[1+add+add2+add3]);

        		if(!rejflag){
				if(rejflag2){
					//accept_S0()
	  				m_S0_acc++;   
				}else{
					//reject_S0()
					m_S0=old1;
					m_S0_prior=old2;
      					m_prior_en=m_old_prior_en;
      					m_S0_rej++;

					//restore_energy()
      					m_energy=m_old_energy;
				}
        		}else{ 
				//reject_S0()
				m_S0=old1;
				m_S0_prior=old2;
	      			m_prior_en=m_old_prior_en;
	      			m_S0_rej++;

				//restore_energy()
      				m_energy=m_old_energy;
			}
        	}

////////////////////////////////////////////////////////////////////////////     TH

     		for(int fibra=0;fibra<NFIBRES;fibra++){  
			if (leader){ 
				//propose_th()
				old1=m_th[fibra];	
				m_th[fibra]+=prerandN[2+(fibra*3)+add+add2+add3]*m_th_prop[fibra];
		
				//compute_th_prior()
				old2=m_th_prior[fibra];
      	   			if(m_th[fibra]==0){
					m_th_prior[fibra]=0;
		   		}else{
					m_th_prior[fibra]=-log(double(fabs(sin(double(m_th[fibra]))/2)));
	      	   		}
		  		//rejflag=false; /////////////////always false
		
				//compute_prior()
				fm_old_prior_en=fm_prior_en[fibra];
	      	   		fm_prior_en[fibra]=m_th_prior[fibra]+m_ph_prior[fibra]+m_f_prior[fibra];
			}

			__syncthreads();

			//compute_signal()
			//compute_signal_pre	
			for(int i=0; i<ndirs; i++){
				cos_alpha_minus_theta=cos(double(alpha[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]-m_th[fibra]));   
			  	cos_alpha_plus_theta=cos(double(alpha[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]+m_th[fibra]));
				old_angtmp[i*NFIBRES+fibra]=angtmp[i*NFIBRES+fibra];
		     	   	angtmp[i*NFIBRES+fibra]= __dadd_rn(__ddiv_rn (cos(double(m_ph[fibra]-beta[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]))*__dadd_rn(cos_alpha_minus_theta,-cos_alpha_plus_theta),2) ,__ddiv_rn (__dadd_rn(cos_alpha_minus_theta,cos_alpha_plus_theta),2));
		 	   	angtmp[i*NFIBRES+fibra]=angtmp[i*NFIBRES+fibra]*angtmp[i*NFIBRES+fibra];

				compute_signal(&mysig[i*NFIBRES+fibra],&mysigold[i*NFIBRES+fibra],bvals[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK],m_d,model,m_dstd,angtmp[i*NFIBRES+fibra]);
			}

			if (leader){
				//for likehood and reject_f_sum
				fsum=m_f0;
				for(int g=0;g<NFIBRES;g++){  
					fsum+=m_f[g];
				}
	
				//compute_prior()
				compute_prior(&m_prior_en,&m_old_prior_en,m_d_prior,m_S0_prior,&fm_prior_en[0],m_f0_prior,m_tau_prior,m_dstd_prior);
			}

			__syncthreads();

			//compute_likelihood()
			compute_likelihood(leader,idrest,m_d,m_S0,&m_likelihood_en,&m_f[0],mysig,myisosig,&mydata[idrest*MAXNDIRS_PER_THREAD],fsum,&reduction[0],m_f0,rician,m_tau,ndirs);

			__syncthreads();
					
			if (leader){ 
				//compute_energy()
				m_old_energy=m_energy;
	      			m_energy=m_prior_en+m_likelihood_en;

				//test_energy
				tmp=exp(double(m_old_energy-m_energy));
				rejflag2=(tmp>prerandU[2+(fibra*3)+add+add2+add3]);
			}

			__syncthreads();
			
          		//if(!rejflag){
				if(rejflag2){
					//accept_th
		  			if (leader) m_th_acc[fibra]++;   
				}else{
					if (leader){
						//reject_th()
						m_th[fibra]=old1;
						m_th_prior[fibra]=old2;
      						fm_prior_en[fibra]=fm_old_prior_en;
						m_th_rej[fibra]++;						
						
						//restore_energy()
						m_prior_en=m_old_prior_en;
	      					m_energy=m_old_energy;
					}
					//compute_signal_pre undo
					for(int i=0; i<ndirs; i++){
						angtmp[i*NFIBRES+fibra]=old_angtmp[i*NFIBRES+fibra];
						mysig[i*NFIBRES+fibra] = mysigold[i*NFIBRES+fibra];						
      					}
				}
		
			__syncthreads();

///////////////////////////////////////     PH

			if (leader){
				//propose_ph()
				old1=m_ph[fibra];
				m_ph[fibra]+=prerandN[3+(fibra*3)+add+add2+add3]*m_ph_prop[fibra];
					
				//compute_ph_prior()
				old2=m_ph_prior[fibra];
      				m_ph_prior[fibra]=0;
      				//rejflag=false;

				//compute_prior()
				fm_old_prior_en=fm_prior_en[fibra];
      	   			fm_prior_en[fibra]=m_th_prior[fibra]+m_ph_prior[fibra]+m_f_prior[fibra];
			}

			__syncthreads();
				
			//compute_signal()
			//compute_signal_pre
			for(int i=0; i<ndirs; i++){
				cos_alpha_minus_theta=cos(double(alpha[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]-m_th[fibra]));   
			  	cos_alpha_plus_theta=cos(double(alpha[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]+m_th[fibra]));
				old_angtmp[i*NFIBRES+fibra]=angtmp[i*NFIBRES+fibra];
		     	   	angtmp[i*NFIBRES+fibra]= __dadd_rn(__ddiv_rn (cos(double(m_ph[fibra]-beta[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK]))*__dadd_rn(cos_alpha_minus_theta,-cos_alpha_plus_theta),2) ,__ddiv_rn (__dadd_rn(cos_alpha_minus_theta,cos_alpha_plus_theta),2));
		 	   	angtmp[i*NFIBRES+fibra]=angtmp[i*NFIBRES+fibra]*angtmp[i*NFIBRES+fibra];

				compute_signal(&mysig[i*NFIBRES+fibra],&mysigold[i*NFIBRES+fibra],bvals[idmult*NDIRECTIONS+idrest+i*THREADS_BLOCK],m_d,model,m_dstd,angtmp[i*NFIBRES+fibra]);
			}

			if (leader){
				//for likehood and reject_f_sum
				fsum=m_f0;
				for(int g=0;g<NFIBRES;g++){  
					fsum+=m_f[g];
				}

				//compute_prior()
				compute_prior(&m_prior_en,&m_old_prior_en,m_d_prior,m_S0_prior,&fm_prior_en[0],m_f0_prior,m_tau_prior,m_dstd_prior);
			}

			__syncthreads();

			//compute_likelihood()
			compute_likelihood(leader,idrest,m_d,m_S0,&m_likelihood_en,&m_f[0],mysig,myisosig,&mydata[idrest*MAXNDIRS_PER_THREAD],fsum,&reduction[0],m_f0,rician,m_tau,ndirs);

			__syncthreads();

			if (leader){
				//compute_energy()
				m_old_energy=m_energy;
	      			m_energy=m_prior_en+m_likelihood_en;

				//test_energy
				tmp=exp(double(m_old_energy-m_energy));
				rejflag2=(tmp>prerandU[3+(fibra*3)+add+add2+add3]);	
			}

			__syncthreads();
		        
			//if(!rejflag){
				if(rejflag2){
					//accept_ph()
		  			if (leader) m_ph_acc[fibra]++;   
				}else{
					if (leader){
						//reject_ph()
						m_ph[fibra]=old1;
						m_ph_prior[fibra]=old2;
      						fm_prior_en[fibra]=fm_old_prior_en;
						m_ph_rej[fibra]++;						
						
						//restore_energy()
						m_prior_en=m_old_prior_en;
	      					m_energy=m_old_energy;
					}
					//compute_signal_pre undo
					for(int i=0; i<ndirs; i++){
						angtmp[i*NFIBRES+fibra]=old_angtmp[i*NFIBRES+fibra];
						mysig[i*NFIBRES+fibra] = mysigold[i*NFIBRES+fibra];	
					}
				}
				
			__syncthreads();

////////////////////////////////////////////             F

			if (leader){
				//propose_f()
				old1=m_f[fibra];
				m_f[fibra]+=prerandN[4+(fibra*3)+add+add2+add3]*m_f_prop[fibra];	
					
	     			//compute_f_prior()
	        		old2=m_f_prior[fibra];
				if (m_f[fibra]<=0 | m_f[fibra]>=1) rejflag=true;
	        		else{
		      			if(!can_use_ard ){
		  				m_f_prior[fibra]=0;
					}else{
		  				if(m_lam_jump[fibra]){
							m_f_prior[fibra]=log(double(m_f[fibra]));
						}else{
		    					m_f_prior[fibra]=0;
		  				}
					}
					m_f_prior[fibra]=fudgevalue*m_f_prior[fibra];
					rejflag=false;
	      			}

				//compute_prior()
				fm_old_prior_en=fm_prior_en[fibra];
      	   			fm_prior_en[fibra]=m_th_prior[fibra]+m_ph_prior[fibra]+m_f_prior[fibra];
				
				//for likehood and reject_f_sum
				fsum=m_f0;
				for(int g=0;g<NFIBRES;g++){  
					fsum+=m_f[g];
				}
				//reject_f_sum()
				rejflag2=(fsum>1);

				//compute_prior()
				compute_prior(&m_prior_en,&m_old_prior_en,m_d_prior,m_S0_prior,&fm_prior_en[0],m_f0_prior,m_tau_prior,m_dstd_prior);
			}

			__syncthreads();

			//compute_likelihood()
			compute_likelihood(leader,idrest,m_d,m_S0,&m_likelihood_en,&m_f[0],mysig,myisosig,&mydata[idrest*MAXNDIRS_PER_THREAD],fsum,&reduction[0],m_f0,rician,m_tau,ndirs);

			__syncthreads();
						
			if (leader){
			//compute_energy()
			m_old_energy=m_energy;
		      	m_energy=m_prior_en+m_likelihood_en;

			//test_energy
			tmp=exp(double(m_old_energy-m_energy));
			rejflag3=(tmp>prerandU[4+(fibra*3)+add+add2+add3]);

		      		if(!rejflag){
					if(!rejflag2){
						if(rejflag3){
							//accept_f()
			  				m_f_acc[fibra]++;   
						}else{
							//reject_f()
							m_f[fibra]=old1;
							m_f_prior[fibra]=old2;
	      						fm_prior_en[fibra]=fm_old_prior_en;
	      						m_f_rej[fibra]++;	
						
							//restore_energy()
							m_prior_en=m_old_prior_en;
		      					m_energy=m_old_energy;	
						}
					}else{ 
						//reject_f()	
						m_f[fibra]=old1;
						m_f_prior[fibra]=old2;
	      					fm_prior_en[fibra]=fm_old_prior_en;
	      					m_f_rej[fibra]++;		

						//restore_energy()
						m_prior_en=m_old_prior_en;
		      				m_energy=m_old_energy;	
					}
				}else{
					//reject_f()
					m_f[fibra]=old1;
					m_f_prior[fibra]=old2;
	      				fm_prior_en[fibra]=fm_old_prior_en;
	      				m_f_rej[fibra]++;

					//restore_energy()
					m_prior_en=m_old_prior_en;
		      			m_energy=m_old_energy;		
				}
			}
				
			__syncthreads();

        	}//end while NFIBRES

        	count_update++;
		recordcount++;

		if(((recordcount%record_every)==0)&&leader){
			rd[(idmult*totalrecords)+sample-1]= m_d;
			if(m_include_f0) rf0[(idmult*totalrecords)+sample-1]= m_f0;
			if(rician) rtau[(idmult*totalrecords)+sample-1]= m_tau;
			if(model==2) rdstd[(idmult*totalrecords)+sample-1]= m_dstd;
			rs0[(idmult*totalrecords)+sample-1]=m_S0;
			rlikelihood_energy[(idmult*totalrecords)+sample-1]=m_likelihood_en;
			for(int j=0;j<NFIBRES;j++){
				rth[(idmult*totalrecords*NFIBRES)+(j*totalrecords)+sample-1]=m_th[j];
				rph[(idmult*totalrecords*NFIBRES)+(j*totalrecords)+sample-1]=m_ph[j];
				rf[(idmult*totalrecords*NFIBRES)+(j*totalrecords)+sample-1]=m_f[j];
			}
			multirecords[(idmult*totalrecords)+sample-1]=sample;
			sample++;
        	}

        	if(((count_update%updateproposalevery)==0)&&leader){
			//m_multifibre.update_proposals();
			m_d_prop*=sqrt(float(m_d_acc+1)/float(m_d_rej+1));
			m_d_prop=min(m_d_prop,maxfloat);

			if(rician){
				m_tau_prop*=sqrt(float(m_tau_acc+1)/float(m_tau_rej+1));
				m_tau_prop=min(m_tau_prop,maxfloat);
				m_tau_acc=0; 
				m_tau_rej=0;	
			}
			if(m_include_f0){
				m_f0_prop*=sqrt(float(m_f0_acc+1)/float(m_f0_rej+1));
				m_f0_prop=min(m_f0_prop,maxfloat);
				m_f0_acc=0; 
				m_f0_rej=0;
			}
			if(model==2){
				m_dstd_prop*=sqrt(float(m_dstd_acc+1)/float(m_dstd_rej+1));
				m_dstd_prop=min(m_dstd_prop,maxfloat);
				m_dstd_acc=0; 
				m_dstd_rej=0;	
			}
			m_S0_prop*=sqrt(float(m_S0_acc+1)/float(m_S0_rej+1));
			m_S0_prop=min(m_S0_prop,maxfloat);
			m_d_acc=0; 
			m_d_rej=0;
			m_S0_acc=0; 
			m_S0_rej=0;
			for(int f=0; f<NFIBRES;f++){
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

	/*if(leader){
		multifibres[idmult].m_prior_en=m_prior_en;
		multifibres[idmult].m_d=m_d;

		multifibres[idmult].m_energy=m_energy;		

		multifibres[idmult].m_d_prop=m_d_prop;
		multifibres[idmult].m_d_prior=m_d_prior;

		for(int f=0;f<NFIBRES;f++){
			fibres[(idmult*NFIBRES)+f].m_th=m_th[f];
			fibres[(idmult*NFIBRES)+f].m_ph=m_ph[f];
			fibres[(idmult*NFIBRES)+f].m_f=m_f[f];
		}

		multifibres[idmult].m_S0_prior=m_S0_prior;
	
		multifibres[idmult].m_S0=m_S0;

		multifibres[idmult].m_likelihood_en=m_likelihood_en;

		multifibres[idmult].m_d_acc=m_d_acc;
		multifibres[idmult].m_d_rej=m_d_rej;
			
		multifibres[idmult].m_S0_acc=m_S0_acc;
		multifibres[idmult].m_S0_rej=m_S0_rej;

		multifibres[idmult].m_S0_prop=m_S0_prop;

		if(m_include_f0){
			multifibres[idmult].m_f0_prior=m_f0_prior;
			multifibres[idmult].m_f0=m_f0;
			multifibres[idmult].m_f0_acc=m_f0_acc;
			multifibres[idmult].m_f0_rej=m_f0_rej;
			multifibres[idmult].m_f0_prop=m_f0_prop;	
		}

		if(rician){
			multifibres[idmult].m_tau_prior=m_tau_prior;
			multifibres[idmult].m_tau=m_tau;
			multifibres[idmult].m_tau_acc=m_tau_acc;
			multifibres[idmult].m_tau_rej=m_tau_rej;
			multifibres[idmult].m_tau_prop=m_tau_prop;	
		}

		if(model==2){
			multifibres[idmult].m_dstd_prior=m_dstd_prior;
			multifibres[idmult].m_dstd=m_dstd;
			multifibres[idmult].m_dstd_acc=m_dstd_acc;
			multifibres[idmult].m_dstd_rej=m_dstd_rej;
			multifibres[idmult].m_dstd_prop=m_dstd_prop;
		}
	
		for(int f=0;f<NFIBRES;f++){
			fibres[(idmult*NFIBRES)+f].m_th_acc=m_th_acc[f];
			fibres[(idmult*NFIBRES)+f].m_th_rej=m_th_rej[f];
			fibres[(idmult*NFIBRES)+f].m_ph_acc=m_ph_acc[f];
			fibres[(idmult*NFIBRES)+f].m_ph_rej=m_ph_rej[f];
			fibres[(idmult*NFIBRES)+f].m_f_acc=m_f_acc[f];
			fibres[(idmult*NFIBRES)+f].m_f_rej=m_f_rej[f];

			fibres[(idmult*NFIBRES)+f].m_prior_en=fm_prior_en[f];
			fibres[(idmult*NFIBRES)+f].m_th_prior=m_th_prior[f];
			fibres[(idmult*NFIBRES)+f].m_ph_prior=m_ph_prior[f];
			fibres[(idmult*NFIBRES)+f].m_f_prior=m_f_prior[f];

			fibres[(idmult*NFIBRES)+f].m_th_prop=m_th_prop[f];
			fibres[(idmult*NFIBRES)+f].m_ph_prop=m_ph_prop[f];
			fibres[(idmult*NFIBRES)+f].m_f_prop=m_f_prop[f];

			fibres[(idmult*NFIBRES)+f].m_lam_jump=m_lam_jump[f];			
		}
	}

	for(int i=0; i<ndirs; i++){	
		isosignals[(idmult*NDIRECTIONS)+idrest+i*THREADS_BLOCK] = myisosig[i];

		for(int f=0;f<NFIBRES;f++){
			signals[(idmult*NDIRECTIONS*NFIBRES)+(f*NDIRECTIONS)+idrest+i*THREADS_BLOCK]=mysig[i*NFIBRES+f];
		}
	}*/
}

