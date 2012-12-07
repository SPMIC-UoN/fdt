#ifndef __LEVENBERG
#define __LEVENBERG

#include "solver_mult_inverse.cu"
#include "diffmodels.cuh"
#include "options.h"

//CPU version in nonlin.h
__device__ const double EPS_gpu = 2.0e-16;       	//Losely based on NRinC 20.1

//CPU version in nonlin.cpp
__device__ inline bool zero_cf_diff_conv(double cfo,double cfn,double cftol){
  	return(2.0*abs(cfo-cfn) <= cftol*(abs(cfo)+abs(cfn)+EPS_gpu));
}

__device__ void levenberg_marquardt_PVM_single_gpu(	//INPUT
							const double*		mydata, 
							const double*		bvecs, 
							const double*		bvals, 
							const int 		nparams,
							const bool 		m_include_f0,
							//INPUT-OUTPUT
							double*			myparams)
{
   	double pcf;
   	double lambda=0.1;
   	double cftol=1.0e-8;
   	double ltol=1.0e20;                  

   	bool success = true;             
   	double olambda = 0.0;              
   	double grad[NPARAMS];                          
   	double hess[NPARAMS*NPARAMS];   
   	double inverse[NPARAMS];
   	double step[NPARAMS];	

   	double ncf=0;

   	int maxiter=200;
   	int niter=0;     

   	pcf=cf_PVM_single(myparams,mydata,bvecs,bvals,nparams,m_include_f0); 
	
   	while (!(success&&niter++ >= maxiter)){ 	//if success we not increase niter (first condition is true)
							//function cost has decreise, we have advanced.
   		if(success){
    			grad_PVM_single(myparams,mydata,bvecs,bvals,nparams,m_include_f0,grad);   
    			hess_PVM_single(myparams,bvecs,bvals,nparams,m_include_f0,hess);  
    		}

    		for (int i=0; i<nparams; i++) {                         
			hess[(i*nparams)+i]+=lambda-olambda;	
    		}

    		solver(hess,grad,nparams,inverse);

    		for (int i=0;i<nparams;i++){
			step[i]=-inverse[i];		
    		}

   		for(int i=0;i<nparams;i++){
			step[i]=myparams[i]+step[i];
   		}

   		ncf = cf_PVM_single(step,mydata,bvecs,bvals,nparams,m_include_f0);

   		if (success = (ncf < pcf)) {
			olambda = 0.0;
        		for(int i=0;i<nparams;i++){
				myparams[i]=step[i];
   			}
        		lambda=lambda/10.0;

			if (zero_cf_diff_conv(pcf,ncf,cftol)){
				return;
			}
			pcf=ncf;
    		}else{
			olambda=lambda;
			lambda=lambda*10.0;
			if(lambda> ltol){ 
				return;
			}
    		}	
   	}
	
	return;
}

__device__ void levenberg_marquardt_PVM_single_c_gpu(	//INPUT
							const double*		mydata, 
							const double*		bvecs, 
							const double*		bvals, 
							const int 		nparams,
							const bool 		m_include_f0,
							//INPUT-OUTPUT
							double*			myparams)
{
   	double pcf;
   	double lambda=0.1;
   	double cftol=1.0e-8;
   	double ltol=1.0e20;                  

   	bool success = true;             
   	double olambda = 0.0;              
   	double grad[NPARAMS];                          
   	double hess[NPARAMS*NPARAMS];   
   	double inverse[NPARAMS];
   	double step[NPARAMS];	

   	double ncf=0;

   	int maxiter=200;
   	int niter=0;     

   	pcf=cf_PVM_single_c(myparams,mydata,bvecs,bvals,nparams,m_include_f0);
	
   	while (!(success&&niter++ >= maxiter)){ 	//if success we not increase niter (first condition is true)
							//function cost has decreise, we have advanced.
   		if(success){
    			grad_PVM_single_c(myparams,mydata,bvecs,bvals,nparams,m_include_f0,grad);   
    			hess_PVM_single_c(myparams,bvecs,bvals,nparams,m_include_f0,hess);  
    		}

    		for (int i=0; i<nparams; i++) {                         
			hess[(i*nparams)+i]+=lambda-olambda;	
    		}

    		solver(hess,grad,nparams,inverse);

    		for (int i=0;i<nparams;i++){
			step[i]=-inverse[i];		
    		}

   		for(int i=0;i<nparams;i++){
			step[i]=myparams[i]+step[i];
   		}

   		ncf = cf_PVM_single_c(step,mydata,bvecs,bvals,nparams,m_include_f0);
		
   		if (success = (ncf < pcf)) {
			olambda = 0.0;
        		for(int i=0;i<nparams;i++){
				myparams[i]=step[i];
   			}
        		lambda=lambda/10.0;

			if (zero_cf_diff_conv(pcf,ncf,cftol)){
				return;
			}
			pcf=ncf;
    		}else{
			olambda=lambda;
			lambda=lambda*10.0;
			if(lambda> ltol){ 
				return;
			}
    		}	
   	}
	
	return;
}


__device__ void levenberg_marquardt_PVM_multi_gpu(	//INPUT
							const double*		mydata, 
							const double*		bvecs, 
							const double*		bvals, 
							const int 		nparams,
							const bool 		m_include_f0,
							//INPUT-OUTPUT
							double*			myparams)
{
   	double pcf;
   	double lambda=0.1;
   	double cftol=1.0e-8;
   	double ltol=1.0e20;                  

   	bool success = true;             
   	double olambda = 0.0;              
   	double grad[NPARAMS];                          
   	double hess[NPARAMS*NPARAMS];   
   	double inverse[NPARAMS];
   	double step[NPARAMS];	

   	double ncf=0;

   	int maxiter=200;
   	int niter=0;     

   	pcf=cf_PVM_multi(myparams,mydata,bvecs,bvals,nparams,m_include_f0);
	
   	while (!(success&&niter++ >= maxiter)){ 	//if success we not increase niter (first condition is true)
							//function cost has decreise, we have advanced.
   		if(success){
    			grad_PVM_multi(myparams,mydata,bvecs,bvals,nparams,m_include_f0,grad);   
    			hess_PVM_multi(myparams,bvecs,bvals,nparams,m_include_f0,hess);  
    		}

    		for (int i=0; i<nparams; i++) {                         
			hess[(i*nparams)+i]+=lambda-olambda;	
    		}

    		solver(hess,grad,nparams,inverse);

    		for (int i=0;i<nparams;i++){
			step[i]=-inverse[i];		
    		}

   		for(int i=0;i<nparams;i++){
			step[i]=myparams[i]+step[i];
   		}

   		ncf = cf_PVM_multi(step,mydata,bvecs,bvals,nparams,m_include_f0);

   		if (success = (ncf < pcf)) {
			olambda = 0.0;
        		for(int i=0;i<nparams;i++){
				myparams[i]=step[i];
   			}
        		lambda=lambda/10.0;

			if (zero_cf_diff_conv(pcf,ncf,cftol)){
				return;
			}
			pcf=ncf;
    		}else{
			olambda=lambda;
			lambda=lambda*10.0;
			if(lambda> ltol){ 
				return;
			}
    		}	
   	}
	
	return;
}
#endif
