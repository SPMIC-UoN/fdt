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

__device__ inline double isoterm_PVM_multi(const int pt,const double _a,const double _b, const double *bvals){
	return exp(double(-_a*log(double(1+bvals[pt]*_b))));
}

__device__ inline double isoterm_a_PVM_multi(const int pt,const double _a,const double _b, const double *bvals){
    	return  -log(double(1+bvals[pt]*_b))*exp(double(-_a*log(double(1+bvals[pt]*_b))));
}

__device__ inline double isoterm_b_PVM_multi(const int pt,const double _a,const double _b, const double *bvals){
      	return -_a*bvals[pt]/(1+bvals[pt]*_b)*exp(double(-_a*log(double(1+bvals[pt]*_b))));
}

__device__ inline double anisoterm_PVM_multi(const int pt,const double _a,const double _b,const double3 x,const double *bvecs, const double *bvals){
	double dp = bvecs[pt]*x.x+bvecs[NDIRECTIONS+pt]*x.y+bvecs[(2*NDIRECTIONS)+pt]*x.z;
  	return exp(double(-_a*log(double(1+bvals[pt]*_b*(dp*dp)))));
}
 
__device__ inline double anisoterm_a_PVM_multi(const int pt,const double _a,const double _b,const double3 x,const double *bvecs, const double *bvals){
	double dp = bvecs[pt]*x.x+bvecs[NDIRECTIONS+pt]*x.y+bvecs[(2*NDIRECTIONS)+pt]*x.z;
  	return -log(double(1+bvals[pt]*(dp*dp)*_b))* exp(double(-_a*log(double(1+bvals[pt]*(dp*dp)*_b))));
}

__device__ inline double anisoterm_b_PVM_multi(const int pt,const double _a,const double _b,const double3 x,const double *bvecs, const double *bvals){
  	double dp = bvecs[pt]*x.x+bvecs[NDIRECTIONS+pt]*x.y+bvecs[(2*NDIRECTIONS)+pt]*x.z;
  	return (-_a*bvals[pt]*(dp*dp)/ (1+bvals[pt]*(dp*dp)*_b)*exp(double(-_a*log(double(1+bvals[pt]*(dp*dp)*_b)))));
}

__device__ inline double anisoterm_th_PVM_multi(const int pt,const double _a,const double _b,const double3 x,const double _th,const double _ph,const double *bvecs, const double *bvals){
  	double dp = bvecs[pt]*x.x+bvecs[NDIRECTIONS+pt]*x.y+bvecs[(2*NDIRECTIONS)+pt]*x.z;
  	double dp1 = cos(double(_th))* (bvecs[pt]*cos(double(_ph)) + bvecs[NDIRECTIONS+pt]*sin(double(_ph))) - bvecs[(2*NDIRECTIONS)+pt]*sin(double(_th));
	return  (-_a*_b*bvals[pt]/(1+bvals[pt]*(dp*dp)*_b)*exp(double(-_a*log(double(1+bvals[pt]*(dp*dp)*_b))))*2*dp*dp1);	
}

__device__ inline double anisoterm_ph_PVM_multi(const int pt,const double _a,const double _b,const double3 x,const double _th,const double _ph,const double *bvecs, const double *bvals){
  	double dp = bvecs[pt]*x.x+bvecs[NDIRECTIONS+pt]*x.y+bvecs[(2*NDIRECTIONS)+pt]*x.z;
  	double dp1 = sin(double(_th))* (-bvecs[pt]*sin(double(_ph)) + bvecs[NDIRECTIONS+pt]*cos(double(_ph)));
  	return  (-_a*_b*bvals[pt]/(1+bvals[pt]*(dp*dp)*_b)*exp(double(-_a*log(double(1+bvals[pt]*(dp*dp)*_b))))*2*dp*dp1);
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
__device__ void sort_PVM_multi(int nfib,int nparams,double* params)
{
	double temp_f, temp_th, temp_ph;
	// Order vector descending using f parameters as index
	for(int i=1; i<(nfib); i++){ 
    		for(int j=0; j<(nfib-i); j++){ 
      			if (params[3+j*3] < params[3+i*3]){ 
        			temp_f = params[3+j*3];
				temp_th = params[3+j*3+1];
				temp_ph = params[3+j*3+2];
        			params[3+j*3] = params[3+i*3]; 
				params[3+j*3+1] = params[3+i*3+1]; 
				params[3+j*3+2] = params[3+i*3+2]; 
        			params[3+i*3] = temp_f; 
				params[3+i*3+1] = temp_th; 
				params[3+i*3+2] = temp_ph; 
      			} 
    		} 
  	} 
}

//in diffmodels.cc -- for calculate residuals
__device__ void  forwardModel_PVM_multi(	//INPUT
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
  	double _a = abs(p[1]);
  	double _b = abs(p[2]);
  	////////////////////////////////////
  	double fs[NFIBRES];
  	double x[NFIBRES*3];	
  	double sumf=0;
	double3 x2;
  	for(int k=0;k<nfib;k++){
    		int kk = 3+3*k;
	    	fs[k] = x2f_gpu(p[kk]);
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
      			val += fs[k]*anisoterm_PVM_multi(i,_a,_b,x2,bvecs,bvals);
    		}	
    		if (m_include_f0){
      			double temp_f0=x2f_gpu(p[nparams-1]);
      			predicted_signal[i] = abs(p[0])*(temp_f0+(1-sumf-temp_f0)*isoterm_PVM_multi(i,_a,_b,bvals)+val);
    		} 
    		else
      			predicted_signal[i] = abs(p[0])*((1-sumf)*isoterm_PVM_multi(i,_a,_b,bvals)+val); 
  		}   
}

//in diffmodels.cc -- for calculate residuals
__device__ void get_prediction_PVM_multi(	//INPUT
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
  	p[1] = params[1]*params[1]/params[2]/params[2];		//m_d*m_d/m_d_std/m_d_std;
  	p[2] = params[2]*params[2]/params[1];			//m_d_std*m_d_std/m_d; // =1/beta
  	for(int k=0;k<nfib;k++){
    		int kk = 3+3*k;
    		p[kk]   = f2x_gpu(params[kk]);
    		p[kk+1] = params[kk+1];
    		p[kk+2] = params[kk+2];
  	}
  	if (m_include_f0)
    		p[nparams-1]=f2x_gpu(params[nparams-1]);
  	forwardModel_PVM_multi(p,bvecs,bvals,nfib,nparams,m_include_f0,predicted_signal);
}

//cost function PVM_multi
__device__ double cf_PVM_multi(		//INPUT
					const double*		params,
					const double*		mdata,
					const double*		bvecs, 
					const double*		bvals,
					const int 		nparams,
					const bool 		m_include_f0)
{
	double cfv = 0.0;
  	double err;	
	double _a= abs(params[1]);
	double _b= abs(params[2]); 
	double fs[NFIBRES];    
	double x[NFIBRES*3];	
	double sumf=0;
	double3 x2;

	for(int k=0;k<NFIBRES;k++){
		int kk = 3+3*(k);
		fs[k] = x2f_gpu(params[kk]);
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
			err += fs[k]*anisoterm_PVM_multi(i,_a,_b,x2,bvecs,bvals); 
    		}
		if(m_include_f0){
			double temp_f0=x2f_gpu(params[nparams-1]);
			err = (abs(params[0])*(temp_f0+((1-sumf-temp_f0)*isoterm_PVM_multi(i,_a,_b,bvals)+err)))-mdata[i];
		}else{
			err = abs(params[0])*((1-sumf)*isoterm_PVM_multi(i,_a,_b,bvals)+err)-mdata[i];
		}
		cfv += err*err;
  	}  
	return(cfv);
}


//gradient function PVM_multi
__device__ void grad_PVM_multi(		//INPUT
					const double*		params,
					const double*		mdata,
					const double*		bvecs, 
					const double*		bvals,
					const int 		nparams,
					const bool 		m_include_f0,
					//OUTPUT
					double*			grad)
{
  	double _a= abs(params[1]);
  	double _b= abs(params[2]);
  	double fs[NFIBRES];
  	double x[NFIBRES*3];	
  	double3 xx;		
  	double sumf=0;

  	for(int k=0;k<NFIBRES;k++){
    		int kk = 3+3*(k);	
    		fs[k] = x2f_gpu(params[kk]);
    		sumf += fs[k];
    		x[k*3] = sin(params[kk+1])*cos(params[kk+2]);
    		x[k*3+1] = sin(params[kk+1])*sin(params[kk+2]);
    		x[k*3+2] = cos(params[kk+1]);
  	}

  	double J[NPARAMS];
  	double diff;
  	double sig;
  
	for (int p=0;p<nparams;p++) grad[p]=0;

  	for(int i=0;i<NDIRECTIONS;i++){
    		sig = 0;
    		for(int a=0;a<nparams;a++) J[a]=0;
    		for(int k=0;k<NFIBRES;k++){
      			int kk = 3+3*(k);
      			xx.x=x[k*3];
      			xx.y=x[k*3+1];
      			xx.z=x[k*3+2];		
      			sig += fs[k]*anisoterm_PVM_multi(i,_a,_b,xx,bvecs,bvals);
      			J[1] += (params[1]>0?1.0:-1.0)*abs(params[0])*fs[k]*anisoterm_a_PVM_multi(i,_a,_b,xx,bvecs,bvals); 
			J[2] += (params[2]>0?1.0:-1.0)*abs(params[0])*fs[k]*anisoterm_b_PVM_multi(i,_a,_b,xx,bvecs,bvals);
			J[kk] = abs(params[0])*(anisoterm_PVM_multi(i,_a,_b,xx,bvecs,bvals)-isoterm_PVM_multi(i,_a,_b,bvals))*two_pi_gpu*sign_gpu(params[kk])*1/(1+params[kk]*params[kk]); 
      			J[kk+1] = abs(params[0])*fs[k]*anisoterm_th_PVM_multi(i,_a,_b,xx,params[kk+1],params[kk+2],bvecs,bvals);  
      			J[kk+2] = abs(params[0])*fs[k]*anisoterm_ph_PVM_multi(i,_a,_b,xx,params[kk+1],params[kk+2],bvecs,bvals);
    		}
    		if(m_include_f0){
			double temp_f0=x2f_gpu(params[nparams-1]);
			J[nparams-1]= abs(params[0])*(1-isoterm_PVM_multi(i,_a,_b,bvals))*two_pi_gpu*sign_gpu(params[nparams-1])*1/(1+params[nparams-1]*params[nparams-1]);
			sig=abs(params[0])*((temp_f0+(1-sumf-temp_f0)*isoterm_PVM_multi(i,_a,_b,bvals))+sig);
    			J[1] += (params[1]>0?1.0:-1.0)*abs(params[0])*(1-sumf-temp_f0)*isoterm_a_PVM_multi(i,_a,_b,bvals);
			J[2] += (params[2]>0?1.0:-1.0)*abs(params[0])*(1-sumf-temp_f0)*isoterm_b_PVM_multi(i,_a,_b,bvals);
    		}else{
	    		sig = abs(params[0]) * ((1-sumf)*isoterm_PVM_multi(i,_a,_b,bvals)+sig);
	    		J[1] += (params[1]>0?1.0:-1.0)*abs(params[0])*(1-sumf)*isoterm_a_PVM_multi(i,_a,_b,bvals);
	    		J[2] += (params[2]>0?1.0:-1.0)*abs(params[0])*(1-sumf)*isoterm_b_PVM_multi(i,_a,_b,bvals);	
    		}
    
    		diff = sig - mdata[i];
    		J[0] = (params[0]>0?1.0:-1.0)*sig/params[0]; 

		for (int p=0;p<nparams;p++) grad[p] += 2*J[p]*diff;
  	}
}


//hessian function PVM_multi 
__device__ void hess_PVM_multi(		//INPUT
					const double*		params,
					const double*		bvecs, 
					const double*		bvals,
					const int 		nparams,
					const bool 		m_include_f0,
					double*			hess)
{
  	double _a= abs(params[1]);
  	double _b= abs(params[2]);
  	double fs[NFIBRES];
  	double x[NFIBRES*3];	
  	double3 xx;
  	double sumf=0;
 
  	for(int k=0;k<NFIBRES;k++){
    		int kk = 3+3*(k);	
    		fs[k] = x2f_gpu(params[kk]);
    		sumf += fs[k];
    		x[k*3] = sin(params[kk+1])*cos(params[kk+2]);
    		x[k*3+1] = sin(params[kk+1])*sin(params[kk+2]);
    		x[k*3+2] = cos(params[kk+1]);
  	}

  	double J[NPARAMS];
  	double sig;

	for (int p=0;p<nparams;p++){
		for (int p2=0;p2<nparams;p2++){ 
			hess[p*nparams+p2] = 0;
		}
	}
  
  	for(int i=0;i<NDIRECTIONS;i++){
    		sig = 0;
    		for(int a=0;a<nparams;a++) J[a]=0;
    		for(int k=0;k<NFIBRES;k++){
      			int kk = 3+3*(k);
      			xx.x=x[k*3];
      			xx.y=x[k*3+1];
      			xx.z=x[k*3+2];		
      			sig += fs[k]*anisoterm_PVM_multi(i,_a,_b,xx,bvecs,bvals);
      			double cov = two_pi_gpu*sign_gpu(params[kk])*1/(1+params[kk]*params[kk]);	
      			J[1] += (params[1]>0?1.0:-1.0)*abs(params[0])*fs[k]*anisoterm_a_PVM_multi(i,_a,_b,xx,bvecs,bvals);
			J[2] += (params[2]>0?1.0:-1.0)*abs(params[0])*fs[k]*anisoterm_b_PVM_multi(i,_a,_b,xx,bvecs,bvals);
			J[kk] = abs(params[0])*(anisoterm_PVM_multi(i,_a,_b,xx,bvecs,bvals)-isoterm_PVM_multi(i,_a,_b,bvals))*cov;
      			J[kk+1] = abs(params[0])*fs[k]*anisoterm_th_PVM_multi(i,_a,_b,xx,params[kk+1],params[kk+2],bvecs,bvals);
      			J[kk+2] = abs(params[0])*fs[k]*anisoterm_ph_PVM_multi(i,_a,_b,xx,params[kk+1],params[kk+2],bvecs,bvals);
    		}
    		if(m_include_f0){
			double temp_f0=x2f_gpu(params[nparams-1]);
			J[nparams-1]= abs(params[0])*(1-isoterm_PVM_multi(i,_a,_b,bvals))*two_pi_gpu*sign_gpu(params[nparams-1])*1/(1+params[nparams-1]*params[nparams-1]);
	    		sig = abs(params[0])* (temp_f0+(1-sumf-temp_f0)*isoterm_PVM_multi(i,_a,_b,bvals)+sig);
    			J[1] += (params[1]>0?1.0:-1.0)*abs(params[0])*(1-sumf-temp_f0)*isoterm_a_PVM_multi(i,_a,_b,bvals);
			J[2] += (params[2]>0?1.0:-1.0)*abs(params[0])*(1-sumf-temp_f0)*isoterm_b_PVM_multi(i,_a,_b,bvals);
    		}else{
			sig = abs(params[0])*((1-sumf)*isoterm_PVM_multi(i,_a,_b,bvals)+sig);
	    		J[1] += (params[1]>0?1.0:-1.0)*abs(params[0])*(1-sumf)*isoterm_a_PVM_multi(i,_a,_b,bvals);
	    		J[2] += (params[2]>0?1.0:-1.0)*abs(params[0])*(1-sumf)*isoterm_b_PVM_multi(i,_a,_b,bvals);	
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
extern "C" __global__ void fit_PVM_multi_kernel(	//INPUT
							const double* 		data, 
							const double* 		params_PVM_single_c,
							const double* 		bvecs, 
							const double* 		bvals, 
							const int 		nvox, 
							const int 		nfib, 				
							const bool 		m_include_f0,
							//OUTPUT
							double* 		params)
{
	int id = blockIdx.x * blockDim.x + threadIdx.x;	
  	if (id >=nvox) { return; }	
	
	int nparams;
	if (m_include_f0)
      		nparams = nfib*3 + 4; 
    	else
      		nparams = nfib*3 + 3;

	int nparams_single_c = nparams-1;

   	double myparams[NPARAMS];
	double myparams_aux[NPARAMS];
   	double mydata[NDIRECTIONS];

   	for(int i=0;i<NDIRECTIONS;i++){
		mydata[i]=data[(id*NDIRECTIONS)+i];
   	}

	for(int i=0;i<nparams_single_c;i++){
		myparams_aux[i]=params_PVM_single_c[(id*nparams_single_c)+i];
   	}

  	myparams[0] = myparams_aux[0];			//pvm1.get_s0();
  	myparams[1] = 1.0;				//start with d=d_std
  	for(int i=0,ii=3;i<nfib;i++,ii+=3){
    		myparams[ii] = f2x_gpu(myparams_aux[ii-1]);
    		myparams[ii+1] = myparams_aux[ii];
    		myparams[ii+2] = myparams_aux[ii+1];
  	}
	myparams[2] = myparams_aux[1];   		//pvm1.get_d();
  	if (m_include_f0)
		myparams[nparams-1]=f2x_gpu(myparams_aux[nparams_single_c-1]);

  	//do the fit
	levenberg_marquardt_PVM_multi_gpu(mydata, &bvecs[id*3*NDIRECTIONS], &bvals[id*NDIRECTIONS], nparams, m_include_f0, myparams);

  	// finalise parameters
	//m_s0-myparams[0] 	m_d-myparams[1] 	m_d_std-myparams[2]		m_f-m_th-m_ph-myparams[3,4,5,6 etc..]   	m_f0-myparams[nparams-1]

	myparams_aux[1] = myparams[1];

  	myparams[1] = abs(myparams_aux[1]*myparams[2]);
  	myparams[2] = sqrt(double(abs(myparams_aux[1]*myparams[2]*myparams[2])));
  	for(int i=3,k=0;k<nfib;i+=3,k++){
    		myparams[i]  = x2f_gpu(myparams[i]);
  	}
  	if (m_include_f0)
    		myparams[nparams-1]=x2f_gpu(myparams[nparams-1]);

	sort_PVM_multi(nfib,nparams,myparams);
  	fix_fsum_PVM_multi(m_include_f0,nfib,nparams,myparams);

	for(int i=0;i<nparams;i++){
		params[(id*nparams)+i] = myparams[i];
		//printf("PARAM[%i]: %.20f\n",i,myparams[i]);
   	}
}

//in diffmodel.cc
extern "C" __global__ void get_residuals_PVM_multi_kernel(	//INPUT
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
      		nparams = nfib*3 + 4; 
    	else
      		nparams = nfib*3 + 3;

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

	get_prediction_PVM_multi(myparams, &bvecs[id*3*NDIRECTIONS], &bvals[id*NDIRECTIONS], nfib, nparams, my_include_f0, predicted_signal);

	for(int i=0;i<NDIRECTIONS;i++){		//residuals=m_data-predicted_signal;
		residuals[id*NDIRECTIONS+i]= mydata[i] - predicted_signal[i];
	}
}

