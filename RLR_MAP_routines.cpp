// Nonlinear estimation of diffusion tensor
// using a Gaussian or Rician noise model
//
// RLR_MAP_routines.cpp 
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2011 University of Oxford 
//
// These routines were originally written so as to
// be compiled into a mex-file callable from within
// Matlab. Hence it was written using Matlab/FORTRAN
// style matrices, was written in "pure" C (i.e. did
// not use any C++ features), used LAPACK for matrix
// inversion etc etc. Since the project was considered
// a "dead end" but the software still useful I
// decided on a minimal re-write. Hence all the
// original C quirkiness has been retained and a minimal
// interface has been added.
//

#include <cstdlib>
#include <iostream>
#include <string>
#include <ctime>
#include <cmath>
#include "newmat.h"
#include "newmatio.h"
#include "RLR_MAP_routines.h"

using namespace std;
using namespace NEWMAT;
namespace NLDTIFIT {

/*
   This file contains routines used for ML/MAP estimation
   of DTI data. It works directly with a spectral decomposition
   D=R*L*R' of the diffusion tensor and estimates three angles
   (R being a function of these) and three eigenvalues (the 
   diagonal elements of L).
   There is a quite strict procedural division between code
   for the Rician noise model and that for the gaussian noise
   model.
*/
/*
   Given a set of gradients, b-values, data-points and initial
   guesses it will evaluate the MAP estimate of the parameters
   pertaining to our RLR data model with a Rician noise-model
   as described in the paper.
*/
int fit_RLR(/* Input */
            double   *g,      /* Gradients, nx3 Matlab matrix. */
            double   *b,      /* b-values. */
            double   *y,      /* Data. */
            int      n,       /* # of data points. */
            int      nvox,    /* # of voxels. */
            double   *pmu,    /* Means of priors. */
            double   *pvar,   /* Variances of priors. */
            int      *pi,     /* Indicies (into theta) of priors. */
            int      np,      /* # of parameters with priors. */
            int      nm,      /* Noise model 1->Gauss, 2->Rician. */
            /* Input/Output */
            double   *theta,  /* Parameters. */
            double   *ei,     /* Expected intensities. */
            double   *hist)   /* History of log-likelihoods. */
{
   double   gg[3];       /* Used to repack gradients. A drag! */
   double   L[3];        /* Diagonal 3x3 matrix of eigenvalues. */
   double   R[9];        /* 3x3 rotation matrix. */
   double   f[MLEN];     /* f_l as defined in paper. */
   double   z[MLEN*MVOX];/* z_l as defined in paper. */
   double   ef[MLEN];    /* exp(f_l). */
   double   e2f[MLEN];   /* exp(2*f_l). */
   double   ll=0.0;      /* Current log-likelihood. */
   double   bll=0.0;     /* Best log-likelihood.    */
   /*
   These are a set of variables intended to hold various
   instances of the gradient and the information matrix.
   For the single voxel case we attempt to estimate 8
   parameters and and the sizes of corresponding variables 
   are 8 and 8*8 respectively. For the more general case
   where we attempt to estimate parameters for multiple
   voxels, pooling the variance estimate, we attempt to
   estimate 1+7*nvox parameters. 
   */
   double   grad[1+7*MVOX];               /* Gradient of cost-function. */
   double   finf[(1+7*MVOX)*(1+7*MVOX)];  /* Iformation matrix. */
   double   btheta[1+7*MVOX];             /* Best theta. */
   double   bgrad[1+7*MVOX];              /* Best gradient. */
   double   bfinf[(1+7*MVOX)*(1+7*MVOX)]; /* Best information matrix. */
   double   ttheta[8];                    /* Tmp storage for single voxel parameters. */
   double   dw = 0.1;                     /* Extra weight to diagonal. */
   double   es2 = 0.0;
   int      i=0, l=0;
   int      vox=0;
   int      iter=0;
   int      npar = 1+7*nvox;
   int      indx[8];

   // Newmat Matrix and vectors for solving equation system
   // Replaces the original LAPACK based solution.
   Matrix       nmH(npar,npar);
   ColumnVector nmGrad(npar);
   ColumnVector nmDTheta(npar);         

   if (nm==GAUS && nvox!=1)
   {
      printf("\nMulti-voxel fitting not supported for Gaussian noise model (pointless)");
      return(-1);
   }
   if (nm==RICE)
   {
      for (vox=0; vox<nvox; vox++)
      {
	 for (l=0; l<n; l++)
	 {
	    if (y[vox*n+l] < 0.0)
	    {
	       printf("\nNegative data, voxel #%d, measurement #%d",vox+1,l+1);
               return(-1);
            }
	 }
      }
   }
  
   memset(hist,0,MAXITER*sizeof(double));
   fix_phi2(theta,nvox);  /* Make sure phi_2 within valid range. */

   for (iter=0; iter<MAXITER; iter++)
   {
      /* 
         Get likelihood. In the process we calculate also
         values that may be useful for later calculating 
         the gradient and the information matrix.
      */
      ll = 0;
      for (vox=0; vox<nvox; vox++)
      {
	 indx[0] = 0; ttheta[0] = theta[0];
         for (i=1; i<8; i++) 
         {
            indx[i] = vox*7 + i; 
            ttheta[i] = theta[indx[i]];
         }
         es2 = exp(ttheta[S2_I]);
         for (i=0; i<3; i++) {L[i] = exp(ttheta[L_I+i]);} /* L */
         make_R(&(ttheta[PHI_I]),R);                      /* R */
         for (l=0; l<n; l++)
         {
            for (i=0; i<3; i++) {gg[i] = g[gindx(l,i,n)];}
            f[l] = make_f(gg,L,R,b[l]);      /* Calculate f_l. */
         }
         if (nm==RICE)
         {
            make_ei(ttheta,f,n,&(ei[vox*n]));                                    /* Expected intensities. */
            make_z(&(ei[vox*n]),&(y[vox*n]),es2,n,&(z[vox*MLEN]));               /* Calculate z_l. */
            ll += log_like_ric(ttheta,&(y[vox*n]),&(z[vox*MLEN]),&(ei[vox*n]),n);/* neg-log-likelihood. */
         }
         else if (nm==GAUS) 
         {
	    make_ei_ef_e2f(theta,f,n,ei,ef,e2f);
	    ll += log_like_gauss(theta,y,ei,n);        /* Gaussian neg-log-likelihood. */
         }
	 ll += priors_ll(ttheta,pmu,pvar,pi,np);       /* Add effect of priors. */
      }

      /*
	 See if the last step was successful. 
      */
      if (!iter || (ll < bll)) /* If in right direction or first iteration. */
      {
	 bll = ll;
         dw *= 0.1;

         memset(grad,0,npar*sizeof(double));
         memset(finf,0,SQR(npar)*sizeof(double));
         for (vox=0; vox<nvox; vox++)
         {
	    indx[0] = 0; ttheta[0] = theta[0];
            for (i=1; i<8; i++) 
            {
               indx[i] = vox*7 + i; 
               ttheta[i] = theta[indx[i]];
            }
            for (i=0; i<3; i++) {L[i] = exp(ttheta[L_I+i]);} /* L */
            make_R(&(ttheta[PHI_I]),R);                      /* R */

            /* New gradient and Hessian. */

            if (nm==RICE) 
            {
               make_ric_grad(ttheta,g,L,R,b,&(y[vox*n]),&(z[vox*MLEN]),&(ei[vox*n]),n,npar,indx,grad,finf);
            }
            else if (nm==GAUS) 
            {
               make_gauss_grad(theta,g,L,R,b,y,ei,ef,e2f,n,grad,finf);
            }

            /* Add contribution from priors. */

            priors_grad(theta,pmu,pvar,pi,indx,np,npar,grad,finf);
         }

         memcpy(btheta,theta,npar*sizeof(double));
         memcpy(bgrad,grad,npar*sizeof(double));
         memcpy(bfinf,finf,SQR(npar)*sizeof(double));
      }
      else
      {
	 ll = bll;
	 dw *= 10;
         memcpy(theta,btheta,npar*sizeof(double));
         memcpy(grad,bgrad,npar*sizeof(double));
         memcpy(finf,bfinf,SQR(npar)*sizeof(double));
         if (dw > 1e20) {break;}
      }
      hist[iter] = ll;

      /*
         Make information matrix more diagonal dominant.
      */
      for (i=0; i<npar; i++) {finf[Hindx(i,i,npar)] *= (1.0 + dw);}
      //
      // Repack into NEWMAT matrix and vector
      //
      for (int ii=0; ii<npar; ii++) { nmGrad(ii+1) = grad[ii]; for (int jj=0; jj<npar; jj++) nmH(ii+1,jj+1) = finf[Hindx(ii,jj,npar)]; }
      /*
         Solve for new values of theta.
      */
      try {
        nmDTheta = nmH.i() * nmGrad;
        for (int ii=0; ii<npar; ii++) theta[ii] -= nmDTheta(ii+1);
        fix_phi2(theta,nvox); /* Fix phi_2 range error. */
      }
      catch(...) {
        // Deliberately empty. Doing nothing means theta is unchanged
        // and the information matrix will be made more diagonal dominant.
      }
   }

   /*
      See if we want to use latest theta or not.
   */
   ll = 0;
   for (vox=0; vox<nvox; vox++)
   {
      indx[0] = 0; ttheta[0] = theta[0];
      for (i=1; i<8; i++) 
      {
         indx[i] = vox*7 + i; 
         ttheta[i] = theta[indx[i]];
      }
      es2 = exp(ttheta[S2_I]);
      for (i=0; i<3; i++) {L[i] = exp(ttheta[L_I+i]);} /* L */
      make_R(&(ttheta[PHI_I]),R);                      /* R */
      for (l=0; l<n; l++)
      {
         for (i=0; i<3; i++) {gg[i] = g[gindx(l,i,n)];}
         f[l] = make_f(gg,L,R,b[l]);      /* Calculate f_l. */
      }

      if (nm==RICE)
      {
         make_ei(ttheta,f,n,&(ei[vox*n]));                       /* Calculate expected intensities. */
         make_z(&(ei[vox*n]),&(y[vox*n]),es2,n,&(z[vox*MLEN]));               /* Calculate z_l. */
         ll += log_like_ric(ttheta,&(y[vox*n]),&(z[vox*MLEN]),&(ei[vox*n]),n); /* Rician neg-log-likelihood. */
      }
      else if (nm==GAUS) 
      {
         make_ei_ef_e2f(theta,f,n,ei,ef,e2f);
         ll += log_like_gauss(theta,y,ei,n);         /* Gaussian neg-log-likelihood. */
      }
      ll += priors_ll(theta,pmu,pvar,pi,np);        /* Add effect of priors. */
   }

   if (bll < ll) 
   { 
      memcpy(theta,btheta,npar*sizeof(double)); 
      /*
         Use best theta to calculate expected intensities.
         Mainly for debugging purposes.
      */
      for (vox=0; vox<nvox; vox++)
      {
         indx[0] = 0; ttheta[0] = theta[0];
         for (i=1; i<8; i++) 
         {
            indx[i] = vox*7 + i; 
            ttheta[i] = theta[indx[i]];
         }
         es2 = exp(ttheta[S2_I]);
         for (i=0; i<3; i++) {L[i] = exp(ttheta[L_I+i]);} /* L */
         make_R(&(ttheta[PHI_I]),R);                      /* R */
         for (l=0; l<n; l++)
         {
            for (i=0; i<3; i++) {gg[i] = g[gindx(l,i,n)];}
            f[l] = make_f(gg,L,R,b[l]);             /* Calculate f_l. */
         }
         make_ei(ttheta,f,n,&(ei[vox*n]));                     /* Calculate expected intensities. */
      }
   }

   return(1);         
}

/*
   Given a set of gradients, b-values, data-points and parameters
   it will return the log-likelihood its gradient, information
   matrix and expected data values.
*/
int get_ric_ll_info(/* Input */
                    double  *theta,
                    double  *y,
                    double  *g,
                    double  *b,
                    int     n,
                    /* Input/output */
                    double  *ll,
                    double  *grad,
                    double  *info,
                    double  *ey)
{
   int      i= 0;
   int      l = 0;
   int      indx[8];
   double   es2 = 0.0;
   double   gg[3];       /* Used to repack gradients. A drag! */
   double   L[3];        /* Diagonal 3x3 matrix of eigenvalues. */
   double   R[9];        /* 3x3 rotation matrix. */
   double   f[MLEN];     /* f_l as defined in paper. */
   double   z[MLEN*MVOX];/* z_l as defined in paper. */

   memset(grad,0,8*sizeof(double));
   memset(info,0,8*8*sizeof(double));
   memset(ey,0,n*sizeof(double));
   es2 = exp(theta[S2_I]);
   for (i=0; i<8; i++) 
   {
      indx[i] = i; 
   }

   for (i=0; i<3; i++) {L[i] = exp(theta[L_I+i]);} /* L */
   make_R(&(theta[PHI_I]),R);                      /* R */
   for (l=0; l<n; l++)
   {
      for (i=0; i<3; i++) {gg[i] = g[gindx(l,i,n)];}
      f[l] = make_f(gg,L,R,b[l]);         /* Calculate f_l. */
   }
   make_ei(theta,f,n,ey);                 /* Expected intensities. */
   make_z(ey,y,es2,n,z);                  /* Calculate z_l. */
   *ll = log_like_ric(theta,y,z,ey,n);    /* Rician neg-log-likelihood. */
   make_ric_grad(theta,g,L,R,b,y,z,ey,n,8,indx,grad,info);


   return(1);
}

int get_gaus_ll_info(/* Input */
                     double  *theta,
                     double  *y,
                     double  *g,
                     double  *b,
                     int     n,
                     /* Input/output */
                     double  *ll,
                     double  *grad,
                     double  *info,
                     double  *ey)
{
   int      i= 0;
   int      l = 0;
   double   es2 = 0.0;
   double   gg[3];       /* Used to repack gradients. A drag! */
   double   L[3];        /* Diagonal 3x3 matrix of eigenvalues. */
   double   R[9];        /* 3x3 rotation matrix. */
   double   f[MLEN];     /* f_l as defined in paper. */
   double   ef[MLEN];    /* exp(f_l). */
   double   e2f[MLEN];   /* exp(2*f_l). */

   memset(grad,0,8*sizeof(double));
   memset(info,0,8*8*sizeof(double));
   memset(ey,0,n*sizeof(double));
                      
   es2 = exp(theta[S2_I]);
   for (i=0; i<3; i++) {L[i] = exp(theta[L_I+i]);} /* L */
   make_R(&(theta[PHI_I]),R);                      /* R */
   for (l=0; l<n; l++)
   {
      for (i=0; i<3; i++) {gg[i] = g[gindx(l,i,n)];}
      f[l] = make_f(gg,L,R,b[l]);         /* Calculate f_l. */
   }
   make_ei_ef_e2f(theta,f,n,ey,ef,e2f);
   *ll = log_like_gauss(theta,y,ey,n);     /* Gaussian neg-log-likelihood. */
   make_gauss_grad(theta,g,L,R,b,y,ey,ef,e2f,n,grad,info);

   return(1);
}

/*
   Calculate the gradient and information matrix for a Rician noise-model.
   There is a quite ugly redundancy in the input parameters. In principle
   we could have calculated L, R, zl and ei given the other parameters.
   However, there are good efficiency reasons for separating the calculation
   of the likelihood and its derivatives. But that means that in the name of
   efficiency we also need to calculate outside these routines and pass to
   them any calculations that are common to them.
*/
void make_ric_grad(/* Input */
                   double   *theta, /* Parameters. */
                   double   *g,     /* Gradient vectors. */
                   double   *L,     /* L matrix, diagonal represented as 3x1. */
                   double   *R,     /* Rotation matrix. */
                   double   *b,     /* b-values. */
                   double   *y,     /* data. */
                   double   *z,     /* z_l as defined in paper. */
                   double   *ei,    /* Expected intensities. */
                   int      n,      /* # of data points. */
                   int      npar,   /* Total # of parameters. */
                   int      *indx,  /* Indicies into grad and finf. */
                   /* Output */
                   double   *grad,  /* Gradient, 8x1. */
                   double   *finf)  /* Information matrix, 8x8. */
{
   double  dfdphi[3];
   double  dfdl[3];
   double  *dRdphi[3];
   double  dRdx[9];
   double  dRdy[9];
   double  dRdz[9];
   double  gg[3];
   double  es2 = 0.0;
   double  ei2_es2 = 0.0;
   double  ei2_es2_2 = 0.0;
   double  zl_I1I0 = 0.0;
   double  te = 0.0;
   int     l=0,i=0,j=0;

   grad[indx[S2_I]] += ((double) n);
   finf[Hindx(indx[S2_I],indx[S2_I],npar)] += ((double) n);
   
   es2 = exp(theta[S2_I]);

   make_dR(&(theta[PHI_I]),dRdx,dRdy,dRdz);
   dRdphi[0]=dRdx; dRdphi[1]=dRdy; dRdphi[2]=dRdz;
   
   for (l=0; l<n; l++)
   {
      /*  
         Repack g, assume g is nx3 Matlab matrix. 
      */
      for (i=0; i<3; i++) {gg[i] = g[gindx(l,i,n)];}
      /*
         Calculate df_l/dphi_i and df_l/dl_i
      */
      for (i=0; i<3; i++)
      {
	 dfdphi[i] = make_dfdphi(gg,L,R,dRdphi[i],b[l]);
         dfdl[i] = make_dfdl(gg,L[i],R,b[l],i);
      }
      ei2_es2 = (SQR(ei[l])) / es2;     /* Square of expected intensity / sigma^2. */
      ei2_es2_2 = SQR(ei2_es2);          /* And squared. */
      zl_I1I0 = z[l] * make_I1I0(z[l]);
      te = make_tricky(ei[l],es2);       /* Tricky expectation. */
      /*
         Start building the gradient.
      */
      grad[indx[S2_I]] += (zl_I1I0 - 0.5*((SQR(y[l])/es2) + ei2_es2));
      grad[indx[S0_I]] += (ei2_es2 - zl_I1I0);
      for (i=0; i<3; i++)
      {
	grad[indx[PHI_I+i]] += (dfdphi[i] * (zl_I1I0 - ei2_es2));
	grad[indx[L_I+i]] += (dfdl[i] * (zl_I1I0 - ei2_es2));
      }
      /*
         And the lower left part of information matrix.
         This (sadly) imposes some assumptions on the order
         of the parameters in theta. I assume 
         theta = [s2 s0 phi lambda].
      */
      /* Partials containing s2. */
      finf[Hindx(indx[S2_I],indx[S2_I],npar)] += (te - ei2_es2_2 - ei2_es2);
      finf[Hindx(indx[S0_I],indx[S2_I],npar)] += (ei2_es2_2 + ei2_es2 - te);
      for (i=0; i<3; i++)
      {
	 finf[Hindx(indx[PHI_I+i],indx[S2_I],npar)] += dfdphi[i] * (te - ei2_es2_2 - ei2_es2);
	 finf[Hindx(indx[L_I+i],indx[S2_I],npar)] += dfdl[i] * (te - ei2_es2_2 - ei2_es2);
      }
      /* Partials containing s0. */
      finf[Hindx(indx[S0_I],indx[S0_I],npar)] += (te - ei2_es2_2);
      for (i=0; i<3; i++)
      {
	 finf[Hindx(indx[PHI_I+i],indx[S0_I],npar)] += dfdphi[i] * (ei2_es2_2 - te);
	 finf[Hindx(indx[L_I+i],indx[S0_I],npar)] += dfdl[i] * (ei2_es2_2 - te);
      }
      /* Partials containing phi. */
      for (i=0; i<3; i++)
      {
	 for (j=i; j<3; j++)
	 {
	    finf[Hindx(indx[PHI_I+j],indx[PHI_I+i],npar)] += dfdphi[i]*dfdphi[j]*(te - ei2_es2_2);
         }
         for (j=0; j<3; j++)
	 {
	    finf[Hindx(indx[L_I+j],indx[PHI_I+i],npar)] += dfdphi[i]*dfdl[j]*(te - ei2_es2_2);
         }
      }
      /* Partials containing lambda. */
      for (i=0; i<3; i++)
      {
	 for (j=i; j<3; j++)
	 {
	    finf[Hindx(indx[L_I+j],indx[L_I+i],npar)] += dfdl[i]*dfdl[j]*(te - ei2_es2_2);
         }
      }
   }
   /*
      Mirror information matrix.
   */
   for (i=0; i<npar; i++)
   {
      for (j=i+1; j<npar; j++)
      {
	 finf[Hindx(i,j,npar)] = finf[Hindx(j,i,npar)];
      }
   }

   /* Hurrah! */

   return; 
}


/*
   Calculate the gradient and information matrix for a Gaussian noise-model.
   The same comments as made for the rician case are pretty much valid here too.
*/
void make_gauss_grad(/* Input */
                     double   *theta, /* Parameters. */
                     double   *g,     /* Gradient vectors. */
                     double   *L,     /* L matrix, diagonal represented as 3x1. */
                     double   *R,     /* Rotation matrix. */
                     double   *b,     /* b-values. */
                     double   *y,     /* data. */
                     double   *ei,    /* Expected intensities. */
                     double   *ef,    /* exp(-fl). */
                     double   *e2f,   /* exp(-2*fl). */
                     int      n,      /* # of data points. */
                     /* Output */
                     double   *grad,  /* Gradient, 8x1. */
                     double   *finf)  /* Information matrix, 8x8. */
{
   double  dfdphi[3];
   double  dfdl[3];
   double  *dRdphi[3];
   double  dRdx[9];
   double  dRdy[9];
   double  dRdz[9];
   double  gg[3];
   double  es2 = 0.0;
   double  es0_es2 = 0.0;
   double  e2s0_es2 = 0.0;
   double  diff = 0.0;
   double  ef_diff = 0.0;
   double  e2f_l = 0.0;
   int     l=0,i=0,j=0;

   memset(grad,0,8*sizeof(double));
   memset(finf,0,8*8*sizeof(double));

   finf[Hindx(S2_I,S2_I,8)] = (((double) n) / 2.0);
   
   es2 = exp(theta[S2_I]);
   es0_es2 = exp(theta[S0_I]) / es2;
   e2s0_es2 = exp(2.0*theta[S0_I]) / es2;

   make_dR(&(theta[PHI_I]),dRdx,dRdy,dRdz);
   dRdphi[0]=dRdx; dRdphi[1]=dRdy; dRdphi[2]=dRdz; 
   
   for (l=0; l<n; l++)
   {
      /*  
         Repack g, assume g is nx3 Matlab matrix. 
      */
      for (i=0; i<3; i++) {gg[i] = g[gindx(l,i,n)];}
      /*
         Calculate df_l/dphi_i and df_l/dl_i
      */
      for (i=0; i<3; i++)
      {
	 dfdphi[i] = make_dfdphi(gg,L,R,dRdphi[i],b[l]);
         dfdl[i] = make_dfdl(gg,L[i],R,b[l],i);
      }
      /*
         Start building the gradient.
      */
      diff = y[l]-ei[l];
      ef_diff = ef[l]*diff;
      e2f_l = e2f[l];
      grad[S2_I] -= SQR(diff);
      grad[S0_I] -= ef_diff;
      for (i=0; i<3; i++)
      {
	 grad[PHI_I+i] += dfdphi[i] * ef_diff;
	 grad[L_I+i] += dfdl[i] * ef_diff;
      }
      /*
         And the lower left part of information matrix.
         This (sadly) imposes some assumptions on the order
         of the parameters in theta. I assume 
         theta = [s2 s0 phi lambda].
      */
      /* The only non-zero partial containing s2 is alreday done. */
      /* Partials containing s0. */
      finf[Hindx(S0_I,S0_I,8)] += e2f_l;
      for (i=0; i<3; i++)
      {
	 finf[Hindx(PHI_I+i,S0_I,8)] -= dfdphi[i] * e2f_l;
	 finf[Hindx(L_I+i,S0_I,8)] -= dfdl[i] * e2f_l;
      }
      /* Partials containing phi. */
      for (i=0; i<3; i++)
      {
	 for (j=i; j<3; j++)
	 {
	    finf[Hindx(PHI_I+j,PHI_I+i,8)] += dfdphi[i]*dfdphi[j]*e2f_l;
         }
         for (j=0; j<3; j++)
	 {
	    finf[Hindx(L_I+j,PHI_I+i,8)] += dfdphi[i]*dfdl[j]*e2f_l;
         }
      }
      /* Partials containing lambda. */
      for (i=0; i<3; i++)
      {
	 for (j=i; j<3; j++)
	 {
	    finf[Hindx(L_I+j,L_I+i,8)] += dfdl[i]*dfdl[j]*e2f_l;
         }
      }
   }

   /* Fix some scalings and stuff. */

   grad[S2_I] /= (2.0 * es2);   
   grad[S2_I] += (((double) n) / 2.0);
   for (i=S0_I; i<8; i++)
   {
      grad[i] *= es0_es2;
      for (j=i; j<8; j++)
      {
	 finf[Hindx(j,i,8)] *= e2s0_es2;
      }
   }

   /*
      Mirror information matrix.
   */
   for (i=0; i<8; i++)
   {
      for (j=i+1; j<8; j++)
      {
	 finf[Hindx(i,j,8)] = finf[Hindx(j,i,8)];
      }
   }

   return; 
}


/* 
   Calculate the efects of priors on the log-likelihood.
*/
double priors_ll(/* Input. */
                 double   *theta,  /* Parameters. */
                 double   *mu,     /* Means of priors. */
                 double   *var,    /* Variance of priors. */
                 int      *indx,   /* Indicies (into theta) of priors. */
                 int      n)       /* # of parameters that have priors. */
{
   double  ll=0.0; 
   int     i=0;

   /* Always add prior pertaining to phi_2. */

   ll -= log(cos(theta[PHI_I+1]));

   /* Then add priors pertaining to "optional priors". */

   for (i=0; i<n; i++)
   {
      ll += ((1.0/(2.0*var[i])) * SQR(mu[i] - theta[indx[i]]));
   }

   return(ll);
}

/* 
   Add the effects of priors onto gradient
   and information matrix.
*/
void priors_grad(/* Input. */
                 double   *theta, /* Parameters. */
                 double   *mu,    /* Means of priors. */
                 double   *var,   /* Variance of priors. */
                 int      *tindx, /* Indicies (into theta) of priors. */
                 int      *indx,  /* Indicies (into grad) of priors. */
                 int      n,      /* # of parameters that have priors. */
                 int      npar,   /* Total # of parameters. */
                 /* Output. */
                 double   *grad,  /* Gradient vector npar*1. */
                 double   *finf)  /* Information matrix npar*npar. */
{
   double  s2 = sin(theta[indx[PHI_I+1]]);
   double  c2 = cos(theta[indx[PHI_I+1]]);
   int     i = 0;
   
   /* Always add prior pertaining to phi_2. */

   grad[indx[PHI_I+1]] += (s2/c2);
   finf[Hindx(indx[PHI_I+1],indx[PHI_I+1],npar)] += (1.0 + SQR(s2/c2));

   /* Then add "optional" priors. */

   for (i=0; i<n; i++)
   {
      grad[indx[tindx[i]]] -= ((1.0/var[i]) * (mu[i] - theta[indx[tindx[i]]]));    
      finf[Hindx(indx[tindx[i]],indx[tindx[i]],npar)] += (1.0/var[i]);
   }

   return;
}

/*
   Create rotation matrix R from three angles phi.
*/
double *make_R(/* Input */
               double  *phi,  /* Angles in radians. */
               /* Output */
               double  *R)    /* Rotation matrix. */
{
   double s1 = sin(phi[0]);
   double c1 = cos(phi[0]);
   double s2 = sin(phi[1]);
   double c2 = cos(phi[1]);
   double s3 = sin(phi[2]);
   double c3 = cos(phi[2]);
   
   R[Rindx(0,0)] = c3*c2;
   R[Rindx(1,0)] = s3*c2;
   R[Rindx(2,0)] = s2;
   R[Rindx(0,1)] = -s3*c1-c3*s2*s1;
   R[Rindx(1,1)] = c3*c1-s3*s2*s1;
   R[Rindx(2,1)] = c2*s1;
   R[Rindx(0,2)] = s3*s1-c3*s2*c1;
   R[Rindx(1,2)] = -c3*s1-s3*s2*c1;
   R[Rindx(2,2)] = c2*c1;

   return(R);
}

/*
   Create partial derivatives of rotation 
   matrix R from three angles phi.
*/
void make_dR(/* Input */
             double  *phi,  /* Angles in radians. */
             /* Output */
             double  *dRdx,
             double  *dRdy,
             double  *dRdz)
{
   double s1 = sin(phi[0]);
   double c1 = cos(phi[0]);
   double s2 = sin(phi[1]);
   double c2 = cos(phi[1]);
   double s3 = sin(phi[2]);
   double c3 = cos(phi[2]);

   dRdx[Rindx(0,0)] = 0.0;
   dRdx[Rindx(1,0)] = 0.0;
   dRdx[Rindx(2,0)] = 0.0;
   dRdx[Rindx(0,1)] = s3*s1-c3*s2*c1;
   dRdx[Rindx(1,1)] = -c3*s1-s3*s2*c1;
   dRdx[Rindx(2,1)] = c2*c1;
   dRdx[Rindx(0,2)] = s3*c1+c3*s2*s1;
   dRdx[Rindx(1,2)] = -c3*c1+s3*s2*s1;
   dRdx[Rindx(2,2)] = -c2*s1;

   dRdy[Rindx(0,0)] = -c3*s2;
   dRdy[Rindx(1,0)] = -s3*s2;
   dRdy[Rindx(2,0)] = c2;
   dRdy[Rindx(0,1)] = -c3*c2*s1;
   dRdy[Rindx(1,1)] = -s3*c2*s1;
   dRdy[Rindx(2,1)] = -s2*s1;
   dRdy[Rindx(0,2)] = -c3*c2*c1;
   dRdy[Rindx(1,2)] = -s3*c2*c1;
   dRdy[Rindx(2,2)] = -s2*c1;

   dRdz[Rindx(0,0)] = -s3*c2;
   dRdz[Rindx(1,0)] = c3*c2;
   dRdz[Rindx(2,0)] = 0.0;
   dRdz[Rindx(0,1)] = -c3*c1+s3*s2*s1;
   dRdz[Rindx(1,1)] = -s3*c1-c3*s2*s1;
   dRdz[Rindx(2,1)] = 0.0;
   dRdz[Rindx(0,2)] = c3*s1+s3*s2*c1;
   dRdz[Rindx(1,2)] = s3*s1-c3*s2*c1;
   dRdz[Rindx(2,2)] = 0.0;

   return; 
}   

/* 
   Efficient calculation of the vector g'*R. Note
   that this is the same vector as R'*g (transposed)
*/
double *make_gR(/* Input */
                double  *g,   /* 3x1 vector. */
                double  *R,   /* 3x3 matrix. */
                /* Output. */
                double  *gR)  /* g'*R or R'*g */
{
   gR[0] = g[0]*R[Rindx(0,0)] + g[1]*R[Rindx(1,0)] + g[2]*R[Rindx(2,0)];
   gR[1] = g[0]*R[Rindx(0,1)] + g[1]*R[Rindx(1,1)] + g[2]*R[Rindx(2,1)];
   gR[2] = g[0]*R[Rindx(0,2)] + g[1]*R[Rindx(1,2)] + g[2]*R[Rindx(2,2)];

   return(gR);
}

/*
   Efficient calculation of the scalar g'*R*L*R'*g
   where g is a 3x1 vector, R a full 3x3 matrix and
   L a 3x3 diagonal matrix (represented by a 3x1 vector).
*/
double make_gRLRg(/* Input */
                  double   *g,
                  double   *L,
                  double   *R)
{
   double  gR[3];
   
   make_gR(g,R,gR);

   return(SQR(gR[0])*L[0] + SQR(gR[1])*L[1] + SQR(gR[2])*L[2]);   
}

/*
   Efficient calculation of the scalar g'*(R1*L*R2'+R2*L*R1')*g
   where g is a 3x1 vector, R1 and R2 are full 3x3 matrices and
   L a 3x3 diagonal matrix (represented by a 3x1 vector).
*/
double make_gR1LR2g(/* Input */
                    double   *g,
                    double   *L,
                    double   *R1,
                    double   *R2)
{
   double  gR1[3];
   double  gR2[3];
   
   make_gR(g,R1,gR1);
   make_gR(g,R2,gR2);

   return(2.0*(gR1[0]*gR2[0]*L[0] + gR1[1]*gR2[1]*L[1] + gR1[2]*gR2[2]*L[2]));   
}

/*
   Efficient (could actually be more efficient) calculation of 
   the scalar g'*R*dL*R'*g where g is a 3x1 vector, R is a full 
   3x3 matrix and L a 3x3 diagonal matrix with a single non-zero
   element. L is represented by that elemnt (l) and the index (li)
   of the non-zero element.
*/
double make_gRdLRg(/* Input */
                   double   *g,
                   double   *R,
                   double   l,
                   int      li)
{
   double  gR[3];

   make_gR(g,R,gR);

   return(l * SQR(gR[li]));
}

/*
   Calculates f_l as defined in paper.
*/
double make_f(/* Input */
              double   *g,
	      double   *L,
	      double   *R,
	      double   b) 
{
  return(b * make_gRLRg(g,L,R));
}

/*
   Calculates z_l as defined in paper.
*/
double *make_z(/* Input */
               double   *ei, /* Expected intensities */
               double   *y,  /* Data (observed intensities). */
               double   es2, /* The estimate of the variance.*/ 
               int      n,   /* Length of the data. */
               /* Output. */
               double   *z)              
{
   int   i=0;

   for (i=0; i<n; i++)
   {
      z[i] = (y[i]*ei[i])/es2;
   }

   return(z);
}

/*
   Calculates \frac{\partial f_l}{\partial \phi}
   as defined in paper. 
*/
double make_dfdphi(/* Input */
                   double   *g,
                   double   *L,
                   double   *R,
                   double   *dRdphi,
                   double   b)
{
   return(b * make_gR1LR2g(g,L,dRdphi,R));
}

/*
   Calculates \frac{\partial f_l}{\partial \lambda}
   as defined in paper. 
*/
double make_dfdl(/* Input */
                 double   *g,
                 double   dRdl,
                 double   *R,
                 double   b,
                 int      dRdl_i)
{
   return(b * make_gRdLRg(g,R,dRdl,dRdl_i));
}

/*
   Calculates the neg-log-likelihood of the data
   assuming a Rician distribution.
*/
double log_like_ric(/* Input */
                    double  *theta, /* Parameters. */
                    double  *y,     /* Data. */
                    double  *z,    /* z_l as defined in paper. */
                    double  *ei,    /* Expected intensities. */
                    int     n)      /* Number of data points. */
{
   double   ll = 0.0;
   int      i = 0;

   for (i=0; i<n; i++) {ll += SQR(y[i]) + SQR(ei[i]);}
   ll /= (2.0 * exp(theta[S2_I]));
   ll += ((double) n) * theta[S2_I];
   for (i=0; i<n; i++)
   {
      ll -= (log(y[i]) + logI0(z[i]));
   }

   return(ll);
}

/*
   Calculates the neg-log-likelihood of the data
   assuming a gaussian distribution.
*/
double log_like_gauss(/* Input */
                      double  *theta, /* Parameters. */
                      double  *y,     /* Data. */
                      double  *ei,    /* Expected intensities. */
                      int     n)      /* Number of data points. */
{
   double   ll = 0.0;
   int      i = 0;

   for (i=0; i<n; i++) {ll += SQR(y[i] - ei[i]);}
   ll /= (2.0 * exp(theta[S2_I]));
   ll += (((double) n)/2.0) * (log(2.0*M_PI)+theta[S2_I]);

   return(ll);
}

/*
   Expected intensities.
*/
double *make_ei(/* Input */
                double   *theta,
                double   *fl,
                int      n,
                /* Output */
                double   *ei)
{
   double   tg[3];
   double   es0 = 0.0;
   int      i = 0;

   es0 = exp(theta[S0_I]);

   for (i=0; i<n; i++)
   {
      ei[i] = es0 * exp(-fl[i]);
   }
   return(ei);
}

/*
   Expected intensities, exp(-fl) and exp(-2*fl).
   These are needed to calculate the gradient and
   information matrix for the Gaussian case.
*/
void make_ei_ef_e2f(/* Input */
                    double   *theta,
                    double   *fl,
                    int      n,
                    /* Output */
                    double   *ei,
                    double   *ef,
                    double   *e2f)
{
   double   tg[3];
   double   es0 = 0.0;
   int      i = 0;

   es0 = exp(theta[S0_I]);

   for (i=0; i<n; i++)
   {
      ef[i] = exp(-fl[i]);
      e2f[i] = SQR(ef[i]);
      ei[i] = es0 * ef[i];
   }
   return;
}

/*
   Returns the log of the incomplete zeroth 
   order Bessel function I_0(x). Constants
   from Abramowitz & Stegun.
*/
double logI0(double  x)
{
   double   t = 0.0;
   double   t2 = 0.0;
   double   tpow = 0.0;
   double   li0 = 0.0;

   if (fabs(x) < 3.75)
   {
      t2 = SQR(x / 3.75);
      li0 = 1;
      li0 += 3.5156229 * (tpow = t2);
      li0 += 3.0899424 * (tpow *= t2);
      li0 += 1.2067492 * (tpow *= t2);
      li0 += 0.2659732 * (tpow *= t2);
      li0 += 0.0360768 * (tpow *= t2);
      li0 += 0.0045813 * (tpow *= t2);
      li0 = log(li0);
   }
   else
   {
      t = 3.75 / x;   /* Inverse of t. */      
      li0 = 0.39894228;
      li0 += 0.01328592 * (tpow = t);
      li0 += 0.00225319 * (tpow *= t);
      li0 -= 0.00157565 * (tpow *= t);
      li0 += 0.00916281 * (tpow *= t);
      li0 -= 0.02057706 * (tpow *= t);
      li0 += 0.02635537 * (tpow *= t);
      li0 -= 0.01647633 * (tpow *= t);
      li0 += 0.00392377 * (tpow *= t);

      li0 = log(li0) - 0.5*log(x) + x;
   }

   return(li0);
}

/*
   Returns the log of the incomplete first 
   order Bessel function I_1(x). Constants
   from Abramowitz & Stegun.
*/
double logI1(double  x)
{
   double   t = 0.0;
   double   t2 = 0.0;
   double   tpow = 0.0;
   double   li1 = 0.0;

   if (fabs(x) < 3.75)
   {
      t2 = SQR(x / 3.75);
      li1 = 0.5;
      li1 += 0.87890594 * (tpow = t2);
      li1 += 0.51498869 * (tpow *= t2);
      li1 += 0.15084934 * (tpow *= t2);
      li1 += 0.02658733 * (tpow *= t2);
      li1 += 0.00301532 * (tpow *= t2);
      li1 += 0.00032411 * (tpow *= t2);
      li1 = log(li1) + log(x);
   }
   else
   {
      t = 3.75 / x;   /* Inverse of t. */      
      li1 = 0.39894228;
      li1 -= 0.03988024 * (tpow = t);
      li1 -= 0.00362018 * (tpow *= t);
      li1 += 0.00163801 * (tpow *= t);
      li1 -= 0.01031555 * (tpow *= t);
      li1 += 0.02282967 * (tpow *= t);
      li1 -= 0.02895312 * (tpow *= t);
      li1 += 0.01787654 * (tpow *= t);
      li1 -= 0.00420059 * (tpow *= t);

      li1 = log(li1) - 0.5*log(x) + x;
   }

   return(li1);
}

/*
   Calculates I_1(1)/I_0(x)
*/
double make_I1I0(double   x)
{
   return(exp(logI1(x) - logI0(x)));
}

/*
   Calculates the expectation we could not
   work out analytically.
*/
double make_tricky(double  A,
                   double  s2)
{
   static double b1[] = {0.0, 0.02071246146216, -0.27072813226432, 0.95032573184007, 0.82139386588719};
   static double b2[] = {-0.38563380475534, -0.18815651963905, 1.06556710921996, -0.00835947755286, 1.00035792923264};
   static double b3[] = {-0.50446345284581, 0.00031463569030, 0.99999186450441, 0.00000008937218, 0.99999999964833};
   double x = A/sqrt(s2);

   if (x <= 1.0) {return(b1[0] + x*(b1[1] + x*(b1[2] + x*(b1[3] + x*b1[4]))));}
   else if (x <= 10) {return(b2[0] + x*(b2[1] + x*(b2[2] + x*(b2[3] + x*b2[4]))));}
   else {return(b3[0] + x*(b3[1] + x*(b3[2] + x*(b3[3] + x*b3[4]))));}
}

/*
   This routine will make sure phi_2 (rotation around y-axis)
   stays within the range -pi/2 < phi_2 < pi/2. It does so in
   a very simplistic way, by resetting it to zero. In reality
   there exist a phi_n and L_n such that
 
   g'*R(phi_n)*L_n*R(phi_n)'*g = g'*R(phi)*L*R(phi)'*g

   for all possible vectors g, where phi_n and L_n is the new
   set of angles and eigenvalues with phi_2 within bounds. But
   I haven't been able to find that set. Anyone who wants to try,
   please do.
   It will also set phi_1 and phi_3 to be in the range 0--2*pi,
   but that's more of cosmetics really.
*/

void fix_phi2(double   *theta,
              int      nvox)
{
   int  i=0;
   int  j=0;

   for (i=0; i<nvox; i++)
   {
      if (fabs(theta[i*7+3]) >= M_PI_2)  /* phi_2. */
      {
	 theta[i*7+3] = 0.0;
      }
      for (j=2; j<5; j+=2)
      {
         if (fabs(theta[i*7+j]) > 2.0*M_PI)
         {
	    theta[i*7+j] = fmod(theta[i*7+j],2.0*M_PI);
            if (theta[i*7+j] > 0.0)
	    {
	       theta[i*7+j] += 2.0*M_PI;
            }
	 }
      }
   }
   return;
}

} // end namespace NLDTIFIT
