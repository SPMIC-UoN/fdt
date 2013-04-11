/*  solver_mult_inverse.cu

    Tim Behrens, Saad Jbabdi, Stam Sotiropoulos, Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */

#include "options.h"

//X = A.i() * B . Used in Levenberg-Marquardt
//MATRIX INVERSE AS NEWMAT LU SOLVER
//implemented in NEWMAT:newmat7.cpp GeneralSolvI.
__device__ void solver(	//INPUT
			double *A, 
			double *P,
			int length,
			//TO USE
			double *C,
			double *el,
			int *indx,	
			//OUTPUT
			double *B)  
{  
	//double C[NPARAMS*NPARAMS];

	for(int i=0;i<length;i++){
		for(int j=0;j<length;j++){
			C[i*length+j]=A[i*length+j];
		}
	}
	
 	bool d=true; 
  	//int indx[NPARAMS];

   	double* akk = C;   
	double big = fabs(*akk); 
	int mu = 0; 
	double* ai = akk; 
	int k;

	for (k = 1; k<length; k++){
      		ai += length; 
		const double trybig = fabs(*ai);
      		if (big < trybig){ 
			big = trybig; 
			mu = k; 
		}
   	}

   	if(length) for (k = 0;;){

		indx[k] = mu;
		if (mu != k){
         		double* a1 = C + length*k; 
			double* a2 = C + length*mu; 
			d = !d;
         		int j = length;
         		while (j--){ 
				const double temp = *a1; 
				*a1++ = *a2; 
				*a2++ = temp; 
			}
      		}

      		double diag = *akk; 
		big = 0; 
		mu = k + 1;
      		if (diag != 0){
         		ai = akk; 
			int i = length - k - 1;
         		while (i--){
            			ai += length; 
				double* al = ai; 
				double mult = *al / diag; 
				*al = mult;
            			int l = length - k - 1; 
				double* aj = akk;
				if (l-- != 0){
				
					double aux=al[1]-(mult* *(++aj));
					*(++al) = aux;
					//*(++al) = __dadd_rn (*al,-mult* *(++aj)); //FAIL in cuda 4.2 compiler
					
               				const double trybig = fabs(*al);
               				if (big < trybig){ 
						big = trybig; 
						mu = length - i - 1; 
					}
               				while (l--){ 
						double aux= al[1]-(mult* *(++aj));
						*(++al) = aux;
						//*(++al) = __dadd_rn (*al,-mult* *(++aj)); //FAIL in cuda 4.2 compiler
					}
           			 }
         		}
      		}
      		if (++k == length) break;      
      		akk += length + 1;
   	}


//////////////////////////////

	//double el[NPARAMS];

	for(int e=0;e<length;e++){
		el[e]=P[e];		
    	}
		
   	int j;
	int ii = length; 
	int ip;    
	double temp;
	int i;
     
	for (i=0; i<length; i++){
 		ip = indx[i]; 
		temp = el[ip]; 
		el[ip] = el[i];
		el[i] = temp;
      		if (temp != 0.0) { ii = i; break; }
   	}
	
  	double* bi; 
	double* ai2;
   	i = ii + 1;

  	if (i < length){
      		bi = el + ii; 
		ai2 = C + ii + i * length;
      		for (;;){
         		int ip = indx[i]; 
			double sum = el[ip]; 
			el[ip] = el[i];
         		double* aij = ai2; 
			double* bj = bi; 
			j = i - ii;
         		while (j--){ 
				sum -=  *aij++* *bj++; 
			}
         		el[i] = sum;
         		if (++i == length) break;
         		ai2 += length;
      		}
   	}

   	ai2 = C + length*length;

   	for (i = length - 1; i >= 0; i--){
      		double* bj = el+i; 
		ai2 -= length; 
		double* ajx = ai2+i;
      		double sum = *bj; 
		double diag = *ajx;
      		j = length - i; 
		while(--j){ 
			sum -= *(++ajx)* *(++bj);  
		}
      		el[i] = sum / diag;
			
   	}
	for(int e=0;e<length;e++){
		B[e]=el[e];
    	}
}

