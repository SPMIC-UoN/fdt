/*  diffmodels_utils.h

    Tim Behrens, Saad Jbabdi, Stam Sotiropoulos, Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */

#ifndef __DIFFMODELS_UTILS
#define __DIFFMODELS_UTILS

//defined in diffmodels.h
#define two_pi_gpu 0.636619772f

//defined in diffmodels.h
#define FSMALL_gpu 0.001f

//defined in diffmodels.h
#define f2beta_gpu(f) (asin(sqrt(f)))

//defined in diffmodels.h
#define d2lambda_gpu(d) (sqrt(d)) 

//defined in diffmodels.h
//#define beta2f_gpu(beta) (pow(sinf(beta),2.0))
__device__ inline  float beta2f_gpu(float beta){ 
	float sinbeta= sinf(beta);
	return sinbeta*sinbeta;
}

//defined in diffmodels.h
#define lambda2d_gpu(lambda) (lambda*lambda) 

//defined in diffmodels.h
__device__ inline float f2x_gpu(float x){
	return tan(x/two_pi); 	
}

//defined in diffmodels.h
__device__ inline float x2f_gpu(float x){         
	return fabsf(two_pi*atan(x));
}

//defined in miscmaths.h 
__device__ inline  int sign_gpu(int x){ if (x>0) return 1; else { if (x<0) return -1; else return 0; } }
__device__ inline int sign_gpu(float x){ if (x>0) return 1; else { if (x<0) return -1; else return 0; } }
__device__ inline int sign_gpu(double x){ if (x>0) return 1; else { if (x<0) return -1; else return 0; } }


#endif
