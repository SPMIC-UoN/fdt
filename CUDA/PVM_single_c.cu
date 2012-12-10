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
  	return exp(double(-bvals[pt]*_d));
}

__device__ 
inline double isoterm_lambda_PVM_single_c(const int pt,const double lambda,const double *bvals){
  	return(-2*bvals[pt]*lambda*exp(double(-bvals[pt]*lambda*lambda)));
}

__device__ 
inline double anisoterm_PVM_single_c(const int pt,const double _d,const double3 x, const double *bvecs, const double *bvals){
	double dp = bvecs[pt]*x.x+bvecs[NDIRECTIONS+pt]*x.y+bvecs[(2*NDIRECTIONS)+pt]*x.z;
	return exp(double(-bvals[pt]*_d*dp*dp));
}

__device__ 
inline double anisoterm_lambda_PVM_single_c(const int pt,const double lambda,const double3 x, const double *bvecs, const double *bvals){
	double dp = bvecs[pt]*x.x+bvecs[NDIRECTIONS+pt]*x.y+bvecs[(2*NDIRECTIONS)+pt]*x.z;
	return(-2*bvals[pt]*lambda*dp*dp*exp(double(-bvals[pt]*lambda*lambda*dp*dp)));
}

__device__ 
inline double anisoterm_th_PVM_single_c(const int pt,const double _d,const double3 x, const double _th,const double _ph,const double *bvecs, const double *bvals){

	double dp = bvecs[pt]*x.x+bvecs[NDIRECTIONS+pt]*x.y+bvecs[(2*NDIRECTIONS)+pt]*x.z;
	double dp1 = cos(double(_th))*(bvecs[pt]*cos(double(_ph))+bvecs[NDIRECTIONS+pt]*sin(double(_ph)))-bvecs[(2*NDIRECTIONS)+pt]*sin(double(_th));
  	return(-2*bvals[pt]*_d*dp*dp1*exp(double(-bvals[pt]*_d*dp*dp)));
}

__device__ 
inline double anisoterm_ph_PVM_single_c(const int pt,const double _d,const double3 x, const double _th,const double _ph,const double *bvecs, const double *bvals){
  	double dp = bvecs[pt]*x.x+bvecs[NDIRECTIONS+pt]*x.y+bvecs[(2*NDIRECTIONS)+pt]*x.z;
	double dp1 = sin(double(_th))*(-bvecs[pt]*sin(double(_ph))+bvecs[NDIRECTIONS+pt]*cos(double(_ph)));
  	return(-2*bvals[pt]*_d*dp*dp1*exp(double(-bvals[pt]*_d*dp*dp)));
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
__device__ void sort_PVM_single_c(int nfib,int nparams,double* params)
{
	double temp_f, temp_th, temp_ph;
	// Order vector descending using f parameters as index
  	for(int i=1; i<(nfib); i++){ 
    		for(int j=0; j<(nfib-i); j++){ 
      			if (params[2+j*3] < params[2+i*3]){ 
        			temp_f = params[2+j*3];
				temp_th = params[2+j*3+1];
				temp_ph = params[2+j*3+2];
        			params[2+j*3] = params[2+i*3]; 
				params[2+j*3+1] = params[2+i*3+1]; 
				params[2+j*3+2] = params[2+i*3+2]; 
        			params[2+i*3] = temp_f; 
				params[2+i*3+1] = temp_th; 
				params[2+i*3+2] = temp_ph; 
      			} 
    		} 
  	} 

	//if (m_return_fanning){
     	 //	fantmp=m_fanning_angles;
      	//	Hess_vec_tmp=m_invprHes_e1;
      	//	Hess=m_Hessian;
  	//}

	//if (m_return_fanning){
     	 	//m_fanning_angles(i)=fantmp(fvals[ii].second);
      		//m_invprHes_e1[i-1]=Hess_vec_tmp[fvals[ii].second-1];
      		//m_Hessian[i-1]=Hess[fvals[ii].second-1];
    	//}
}

//in diffmodels.cc -- for calculate residuals
__device__ void  forwardModel_PVM_single_c(	//INPUT
						const double* 		p,
						const double*		bvecs, 
						const double*		bvals,
						const int		nfib,
						const int 		nparams,
						const bool 		m_include_f0,
						//OUTPUT
						double*		 	predicted_signal)
{
  	for(int i=0;i<NDIRECTIONS;i++){
		predicted_signal[i]=0;		//pred = 0;
	}
  	double val;
  	double _d = lambda2d_gpu(p[1]);
  	////////////////////////////////////
  	double fs[NFIBRES];
  	double x[NFIBRES*3];	
  	double sumf=0;
	double3 x2;
	double partial_fsum;
  	for(int k=0;k<nfib;k++){
    		int kk = 2+3*k;
		////// partial_fsum //////
		partial_fsum=1.0;
		for(int j=0;j<k;j++)
			partial_fsum-=fs[j];
    		//////////////////////////
	    	fs[k] = beta2f_gpu(p[kk])*partial_fsum;
	    	sumf += fs[k];
		x[k*3] = sin(p[kk+1])*cos(p[kk+2]);
    		x[k*3+1] = sin(p[kk+1])*sin(p[kk+2]);
    		x[k*3+2] = cos(p[kk+1]);
  	}
  	////////////////////////////////////
  	for(int i=0;i<NDIRECTIONS;i++){
    		val = 0.0;
    		for(int k=0;k<nfib;k++){
			x2.x=x[k*3];
			x2.y=x[k*3+1];
			x2.z=x[k*3+2];	 
      			val += fs[k]*anisoterm_PVM_single_c(i,_d,x2,bvecs,bvals);
    		}	
    		if (m_include_f0){
			//partial_fsum ///////////
			partial_fsum=1.0;
			for(int j=0;j<NFIBRES;j++)
				partial_fsum-=fs[j];
	     		//////////////////////////
      			double temp_f0= beta2f_gpu(p[nparams-1])*partial_fsum;
      			predicted_signal[i] = p[0]*(temp_f0+(1-sumf-temp_f0)*isoterm_PVM_single_c(i,_d,bvals)+val);
    		} 
    		else
      			predicted_signal[i] = p[0]*((1-sumf)*isoterm_PVM_single_c(i,_d,bvals)+val); 
  	}
}

//in diffmodels.cc -- for calculate residuals
__device__ void get_prediction_PVM_single_c(	//INPUT
						const double*	params,
						const double*	bvecs, 
						const double*	bvals,
						const int 	nfib,
						const int 	nparams,
						const bool 	m_include_f0,
						//OUTPUT
						double* 	predicted_signal)
{
	//m_s0-myparams[0] 	m_d-myparams[1] 	m_d_std-myparams[2]		m_f-m_th-m_ph-myparams[3,4,5,6 etc..]   	m_f0-myparams[nparams-1]
  	double p[NPARAMS];
  	p[0] = params[0];
	if(params[1]<0)  p[1] = 0;	//This can be due to numerical errors..sqrt
  	else p[1] = d2lambda_gpu(params[1]);	
	double partial_fsum;	
	double fs[NFIBRES];
  	for(int k=0;k<nfib;k++){
    		int kk = 2+3*k;
		//partial_fsum ///////////
		partial_fsum=1.0;
		for(int j=0;j<k;j++)
			partial_fsum-=fs[j];
	     	//////////////////////////
		fs[k] = params[kk];
		double tmpr=params[kk]/partial_fsum;
    		if (tmpr>1.0) tmpr=1; //This can be due to numerical errors
		if (tmpr<0.0) tmpr=0; //This can be due to numerical errors..sqrt
    		p[kk]   = f2beta_gpu(tmpr);
    		p[kk+1] = params[kk+1];
    		p[kk+2] = params[kk+2];
  	}
  	if (m_include_f0){
		//partial_fsum ///////////
		partial_fsum=1.0;
		for(int j=0;j<NFIBRES;j++)
			partial_fsum-=fs[j];
	     	//////////////////////////	
		double tmpr=params[nparams-1]/partial_fsum;
    		if (tmpr>1.0) tmpr=1; //This can be due to numerical errors..asin
		if (tmpr<0.0) tmpr=0; //This can be due to numerical errors..sqrt
    		p[nparams-1]= f2beta_gpu(tmpr);	
	}
  	forwardModel_PVM_single_c(p,bvecs,bvals,nfib,nparams,m_include_f0,predicted_signal);
}

//cost function PVM_single_c
__device__ double cf_PVM_single_c(	//INPUT
					const double*		params,
					const double*		mdata,
					const double*		bvecs, 
					const double*		bvals,
					const int 		nparams,
					const bool 		m_include_f0)
{
	double cfv = 0.0;
  	double err;
	double _d = lambda2d_gpu(params[1]);
	double fs[NFIBRES];    
	double x[NFIBRES*3];	
	double sumf=0;
	double3 x2;
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
    		
    		x[k*3] = sin(params[kk+1])*cos(params[kk+2]);
    		x[k*3+1] = sin(params[kk+1])*sin(params[kk+2]);
    		x[k*3+2] = cos(params[kk+1]);
  	}
	
	for(int i=0;i<NDIRECTIONS;i++){
    		err = 0.0;
    		for(int k=0;k<NFIBRES;k++){
			x2.x=x[k*3];
			x2.y=x[k*3+1];
			x2.z=x[k*3+2];	
			err += fs[k]*anisoterm_PVM_single_c(i,_d,x2,bvecs,bvals); 
    		}
		if(m_include_f0){
			//partial_fsum ///////////
			partial_fsum=1.0;
			for(int j=0;j<NFIBRES;j++)
				partial_fsum-=fs[j];
	     		//////////////////////////
	      		double temp_f0=beta2f_gpu(params[nparams-1])*partial_fsum;
			err=(params[0]*((temp_f0+(1-sumf-temp_f0)*isoterm_PVM_single_c(i,_d,bvals))+err))-mdata[i];
		}else{

			err = params[0]*((1-sumf)*isoterm_PVM_single_c(i,_d,bvals)+err)-mdata[i];
		
		}
		cfv += err*err;  	
  	}  

	return(cfv);
}


//gradient function PVM_single_c
__device__ void grad_PVM_single_c(	//INPUT
					const double*		params,
					const double*		mdata,
					const double*		bvecs, 
					const double*		bvals,
					const int 		nparams,
					const bool 		m_include_f0,
					//OUTPUT
					double*			grad)
{
  	double _d = lambda2d_gpu(params[1]);
  	double fs[NFIBRES];
  	double bs[NFIBRES];
  	double x[NFIBRES*3];	
  	double3 xx;		
  	double sumf=0;
  	double partial_fsum;

  	for(int k=0;k<NFIBRES;k++){
    		int kk = 2+3*(k);
    		bs[k]=params[kk];
    		//partial_fsum ///////////
		partial_fsum=1.0;
		for(int j=0;j<k;j++){
			partial_fsum-=fs[j];
		}
   		 //////////////////////////

    		fs[k] = beta2f_gpu(params[kk])*partial_fsum;
    		sumf += fs[k];
    		x[k*3] = sin(params[kk+1])*cos(params[kk+2]);
    		x[k*3+1] = sin(params[kk+1])*sin(params[kk+2]);
    		x[k*3+2] = cos(params[kk+1]);
  	}

  	////////// fraction deriv //////////////
  	// f_deriv=fractions_deriv(nfib, fs, bs);  
  	//////////////////////////////////////
  	double f_deriv[NFIBRES*NFIBRES];
  	double fsum;
  	for (int j=0; j<NFIBRES; j++){
   		 for (int k=0; k<NFIBRES; k++){
			f_deriv[j*NFIBRES+k]=0;
    		}
  	}  

  	for (int j=0; j<NFIBRES; j++){
    		for (int k=0; k<NFIBRES; k++){
      			if (j==k){
				fsum=1; 
				for (int n=0; n<=(j-1); n++)
	  				fsum-=fs[n];
	  			f_deriv[j*NFIBRES+k] = sin(double(2*bs[k]))*fsum;
      			}else if (j>k){
				fsum=0;
				for (int n=0; n<=(j-1); n++)
	  				fsum += f_deriv[n*NFIBRES+k];
				f_deriv[j*NFIBRES+k] += -sin(bs[j])*sin(bs[j])*fsum;

      			}
    		} 	   
  	}
  	///////////////////////////////
  	double J[NPARAMS];
  	double diff;
  	double sig,Iso_term;
  	double Aniso_terms[NFIBRES];

	for (int p=0;p<nparams;p++) grad[p]=0;
  
  	for(int i=0;i<NDIRECTIONS;i++){
    		Iso_term=isoterm_PVM_single_c(i,_d,bvals);  //Precompute some terms for this datapoint
    		for(int k=0;k<NFIBRES;k++){
      			xx.x=x[k*3];
      			xx.y=x[k*3+1];
     			xx.z=x[k*3+2];	
      			Aniso_terms[k]=anisoterm_PVM_single_c(i,_d,xx,bvecs,bvals);
    		}
    		sig = 0;
    		for(int a=0;a<nparams;a++) J[a]=0;
    		for(int k=0;k<NFIBRES;k++){
     			int kk = 2+3*(k);
      			xx.x=x[k*3];
      			xx.y=x[k*3+1];
      			xx.z=x[k*3+2];		
      			sig += fs[k]*Aniso_terms[k];
     			J[1] += params[0]*fs[k]*anisoterm_lambda_PVM_single_c(i,params[1],xx,bvecs,bvals);
     			J[kk] = 0;
      			for (int j=0;j<NFIBRES;j++){
				if(f_deriv[j*NFIBRES+k]!=0){
	  				J[kk] += params[0]*(Aniso_terms[j]-Iso_term)*f_deriv[j*NFIBRES+k]; 
				}
      			}
      			J[kk+1] = params[0]*fs[k]*anisoterm_th_PVM_single_c(i,_d,xx,params[kk+1],params[kk+2],bvecs,bvals);
      			J[kk+2] = params[0]*fs[k]*anisoterm_ph_PVM_single_c(i,_d,xx,params[kk+1],params[kk+2],bvecs,bvals);
    		}
    		if(m_include_f0){
			//partial_fsum ///////////
    			partial_fsum=1.0;
    			for(int j=0;j<(NFIBRES);j++)
				partial_fsum-=fs[j];
			//////////////////////////
			double temp_f0=beta2f_gpu(params[nparams-1])*partial_fsum;

    			//derivative with respect to f0
    			J[nparams-1]= params[0]*(1-Iso_term)*sin(double(2*params[nparams-1]))*partial_fsum; 
			sig=params[0]*((temp_f0+(1-sumf-temp_f0)*Iso_term)+sig);
    			J[1] += params[0]*(1-sumf-temp_f0)*isoterm_lambda_PVM_single_c(i,params[1],bvals);
    		}else{
			sig = params[0]*((1-sumf)*Iso_term+sig);
	    		J[1] += params[0]*(1-sumf)*isoterm_lambda_PVM_single_c(i,params[1],bvals);
    		}
    		diff = sig - mdata[i];
    		J[0] = sig/params[0]; 

		for (int p=0;p<nparams;p++) grad[p] += 2*J[p]*diff;
  	}
}


//hessian function PVM_single_c
__device__ void hess_PVM_single_c(	//INPUT
					const double*		params,
					const double*		bvecs, 
					const double*		bvals,
					const int 		nparams,
					const bool 		m_include_f0,
					//OUTPUT
					double*			hess)
{
  	double _d = lambda2d_gpu(params[1]);
  	double fs[NFIBRES];
  	double bs[NFIBRES];
  	double x[NFIBRES*3];	
  	double3 xx;
  	double sumf=0;
  	double partial_fsum;

  	for(int k=0;k<NFIBRES;k++){
    		int kk = 2+3*(k);
    		bs[k]=params[kk];
    		//partial_fsum ///////////
		partial_fsum=1.0;
		for(int j=0;j<k;j++)
			partial_fsum-=fs[j];
    		//////////////////////////
    		fs[k] = beta2f_gpu(params[kk])*partial_fsum;
    		sumf += fs[k];
    		x[k*3] = sin(params[kk+1])*cos(params[kk+2]);
    		x[k*3+1] = sin(params[kk+1])*sin(params[kk+2]);
    		x[k*3+2] = cos(params[kk+1]);
  	}

  	////////// fraction deriv //////////////
  	// f_deriv=fractions_deriv(nfib, fs, bs);  
  	//////////////////////////////////////
  	double f_deriv[NFIBRES*NFIBRES];
  	double fsum;
  	for (int j=0; j<NFIBRES; j++){
    		for (int k=0; k<NFIBRES; k++){
			f_deriv[j*NFIBRES+k]=0;
    		}
  	}

  	for (int j=0; j<NFIBRES; j++){
    		for (int k=0; k<NFIBRES; k++){
      			if (j==k){
				fsum=1; 
				for (int n=0; n<=(j-1); n++)
	  			fsum-=fs[n];
	  			f_deriv[j*NFIBRES+k] = sin(double(2*bs[k]))*fsum;
      			}
      			else if (j>k){
				fsum=0;
				for (int n=0; n<=(j-1); n++)
	  				fsum += f_deriv[n*NFIBRES+k];
				f_deriv[j*NFIBRES+k] += -sin(bs[j])*sin(bs[j])*fsum;
      			}
    		} 	   
  	}
  	///////////////////////////////
  	double J[NPARAMS];
  	double sig,Iso_term;
  	double Aniso_terms[NFIBRES];

	for (int p=0;p<nparams;p++){
		for (int p2=0;p2<nparams;p2++){ 
			hess[p*nparams+p2] = 0;
		}
	}

  	for(int i=0;i<NDIRECTIONS;i++){
    		Iso_term=isoterm_PVM_single_c(i,_d,bvals);  //Precompute some terms for this datapoint
    		for(int k=0;k<NFIBRES;k++){
      			xx.x=x[k*3];
      			xx.y=x[k*3+1];
      			xx.z=x[k*3+2];	
      			Aniso_terms[k]=anisoterm_PVM_single_c(i,_d,xx,bvecs,bvals);
    		}
    		sig = 0;
    		for(int a=0;a<nparams;a++) J[a]=0;
    		for(int k=0;k<NFIBRES;k++){
      			int kk = 2+3*(k);
      			xx.x=x[k*3];
      			xx.y=x[k*3+1];
      			xx.z=x[k*3+2];		 
      			sig += fs[k]*Aniso_terms[k];
      			J[1] += params[0]*fs[k]*anisoterm_lambda_PVM_single_c(i,params[1],xx,bvecs,bvals);	 
      			J[kk] = 0;
      			for (int j=0; j<NFIBRES; j++){
				if (f_deriv[j*NFIBRES+k]!=0)
	  			J[kk] += params[0]*(Aniso_terms[j]-Iso_term)*f_deriv[j*NFIBRES+k]; 
      			}
      			J[kk+1] = params[0]*fs[k]*anisoterm_th_PVM_single_c(i,_d,xx,params[kk+1],params[kk+2],bvecs,bvals);
      			J[kk+2] = params[0]*fs[k]*anisoterm_ph_PVM_single_c(i,_d,xx,params[kk+1],params[kk+2],bvecs,bvals);
    		}
    		if(m_include_f0){
			//partial_fsum ///////////
	    		partial_fsum=1.0;
	    		for(int j=0;j<(NFIBRES);j++)
				partial_fsum-=fs[j];
	    		//////////////////////////
    			double temp_f0=beta2f_gpu(params[nparams-1])*partial_fsum;
    			//derivative with respect to f0
    			J[nparams-1]= params[0]*(1-Iso_term)*sin(double(2*params[nparams-1]))*partial_fsum; 
			sig= params[0]*((temp_f0+(1-sumf-temp_f0)*Iso_term)+sig);
    			J[1] += params[0]*(1-sumf-temp_f0)*isoterm_lambda_PVM_single_c(i,params[1],bvals);
    		}else{
	    		sig = params[0]*((1-sumf)*Iso_term+sig);
	    		J[1] += params[0]*(1-sumf)*isoterm_lambda_PVM_single_c(i,params[1],bvals);
    		}
    		J[0] = sig/params[0]; 

		for (int p=0;p<nparams;p++){
			for (int p2=p;p2<nparams;p2++){ 
				hess[p*nparams+p2] += 2*(J[p]*J[p2]);
			}
		}
  	}

  	for (int j=0; j<nparams; j++) {
    		for (int i=j+1; i<nparams; i++) {
     			hess[i*nparams+j]=hess[j*nparams+i];	
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
	int id = blockIdx.x * blockDim.x + threadIdx.x;	
   	if (id >=nvox) { return; }	

	int nparams;
	if (m_include_f0)
      		nparams = nfib*3 + 3; 
    	else
      		nparams = nfib*3 + 2;

   	double myparams[NPARAMS];
   	double mydata[NDIRECTIONS];

	for(int i=0;i<nparams;i++){
		myparams[i]=params[(id*nparams)+i];
   	}

   	for(int i=0;i<NDIRECTIONS;i++){
		mydata[i]=data[(id*NDIRECTIONS)+i];
   	}

	//if(id==1) for(int i=0;i<nparams;i++)printf("START[%i]: %.20f\n",i,myparams[i]);
	//do the fit
	levenberg_marquardt_PVM_single_c_gpu(mydata, &bvecs[id*3*NDIRECTIONS], &bvals[id*NDIRECTIONS], nparams, m_include_f0, myparams);

	//double m_BIC;
	//if (m_eval_BIC){  
    	//	double RSS= cf_PVM_single_c(myparams,mydata,&bvecs[id*3*NDIRECTIONS], &bvals[id*NDIRECTIONS], nparams,m_include_f0); // get the sum of squared residuals
    	//	m_BIC=NDIRECTIONS*log(double(RSS/NDIRECTIONS))+log(double(NDIRECTIONS))*nparams; // evaluate BIC
  	//}   NOT USED at the moment

	// finalise parameters
	// m_s0-myparams[0] 	m_d-myparams[1] 	m_f-m_th-m_ph-myparams[2,3,4,5, etc..]   	m_f0-myparams[nparams-1]
	
	double m_f[NFIBRES]; 					// for partial_fsum

  	myparams[1] = lambda2d_gpu(myparams[1]); 
  	for(int k=0;k<nfib;k++){
    		int kk = 2 + 3*(k);
    		myparams[kk]  = beta2f_gpu(myparams[kk])*partial_fsum_PVM_single_c(m_f,k);
		m_f[k]=myparams[kk];
  	}
  
  	//if (m_return_fanning)
    		//Fanning_angles_from_Hessian(); NOT USED at the moment
  
  	if (m_include_f0)
    		myparams[nparams-1]= beta2f_gpu(myparams[nparams-1])*partial_fsum_PVM_single_c(m_f,nfib);

	sort_PVM_single_c(nfib,nparams,myparams);

	for(int i=0;i<nparams;i++){
		params[(id*nparams)+i] = myparams[i];
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
	int id = blockIdx.x * blockDim.x + threadIdx.x;	
   	if (id >=nvox) { return; }	

	int nparams;
	if (m_include_f0)
      		nparams = nfib*3 + 3; 
    	else
      		nparams = nfib*3 + 2;

	bool my_include_f0 = includes_f0[id];

   	double myparams[NPARAMS];
   	double mydata[NDIRECTIONS];

	for(int i=0;i<nparams;i++){
		myparams[i]=params[(id*nparams)+i];
   	}

   	for(int i=0;i<NDIRECTIONS;i++){
		mydata[i]=data[(id*NDIRECTIONS)+i];
   	}

	double predicted_signal[NDIRECTIONS];

	get_prediction_PVM_single_c(myparams, &bvecs[id*3*NDIRECTIONS], &bvals[id*NDIRECTIONS], nfib, nparams, my_include_f0, predicted_signal);

	for(int i=0;i<NDIRECTIONS;i++){		//residuals=m_data-predicted_signal;
		residuals[id*NDIRECTIONS+i]= mydata[i] - predicted_signal[i];
	}
}

