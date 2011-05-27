// Nonlinear estimation of diffusion tensor
// using a Gaussian or Rician noise model
//
// RLR_MAP_routines.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2011 University of Oxford 
//
// These routines were originally written so as to
// be compiled into a mex-file callable from within
// Matlab. Hence it was written using Matlab/FORTRAN
// style "matrices", was written in "pure" C (i.e. did
// not use any C++ features), used LAPACK for matrix
// inversion etc etc. Since the project was considered
// a "dead end" but the software still useful I
// decided on a minimal re-write. Hence all the
// original C quirkiness has been retained and a minmal
// interface has been added.
//

namespace NLDTIFIT {

#ifndef M_PI
#define M_PI 	   3.14159265358979323846   /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923   /* pi/2 */
#endif

/* Need macros for different platforms here. */

#define DGESV  dgesv_ 

/* Used to index Hessian/Information matrix. */
#ifndef Hindx
#define Hindx(R,C,SZ) ((C)*(SZ) + (R)) /* Column first for Fortran (and Matlab). */
#endif

/* Used to index rotation matrix R. */
#ifndef Rindx
#define Rindx(R,C) ((C)*3 + (R)) /* Column first for Fortran (and Matlab). */
#endif

/* Used to index gradient matrix g. */
#ifndef gindx
#define gindx(R,C,M) (((C)*(M)) + (R)) /* Column first for Fortran (and Matlab). */
#endif

/* Used to index eigenvalue matrix L. */
#ifndef Lindx
#define Lindx(R,C) ((C)*3 + (R)) /* Column first for Fortran (and Matlab). */
#endif

#ifndef SQR
#define SQR(A) ((A)*(A))
#endif

#ifndef S2_I
#define S2_I   0
#define S0_I   1
#define PHI_I  2
#define L_I    5
#endif

#ifndef MLEN   /* Longest "time-series" we ever expect. */
#define MLEN   1024
#endif

#ifndef MVOX   /* Largest # of voxels estimated together. */
#define MVOX   27
#endif

#ifndef MAXITER
#define MAXITER  100
#endif

#ifndef GAUSS
#define GAUS  1
#define RICE  2
#endif


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
                    double  *ey);

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
                     double  *ey);
                     
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
                   double   *finf); /* Information matrix, 8x8. */

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
                     double   *finf); /* Information matrix, 8x8. */

double priors_ll(/* Input. */
                 double   *theta,  /* Parameters. */
                 double   *mu,     /* Means of priors. */
                 double   *var,    /* Variance of priors. */
                 int      *indx,   /* Indicies (into theta) of priors. */
                 int      n);      /* # of parameters that have priors. */

void priors_grad(/* Input. */
                 double   *theta, /* Parameters. */
                 double   *mu,    /* Means of priors. */
                 double   *var,   /* Variance of priors. */
                 int      *tindx,  /* Indicies (into theta) of priors. */
                 int      *indx, /* Indicies (into grad) of priors. */
                 int      n,      /* # of parameters that have priors. */
                 int      npar,   /* Total # of parameters. */
                 /* Output. */
                 double   *grad,  /* Gradient vector npar*1. */
                 double   *finf); /* Information matrix npar*npar. */

double *make_R(/* Input */
               double  *phi,  /* Angles in radians. */
               /* Output */
               double  *R);   /* Rotation matrix. */

void make_dR(/* Input */
             double  *phi,  /* Angles in radians. */
             /* Output */
             double  *dRdx,
             double  *dRdy,
             double  *dRdz);

double *make_gR(/* Input */
                double  *g,   /* 3x1 vector. */
                double  *R,   /* 3x3 matrix. */
                /* Output. */
                double  *gR); /* g'*R or R'*g */

double make_gRLRg(/* Input */
                  double   *g,
                  double   *L,
                  double   *R);

double make_gR1LR2g(/* Input */
                    double   *g,
                    double   *L,
                    double   *R1,
                    double   *R2);

double make_gRdLRg(/* Input */
                   double   *g,
                   double   *R,
                   double   l,
                   int      li);

double make_f(/* Input */
              double   *g,
	      double   *L,
	      double   *R,
	      double   b);

double *make_z(/* Input */
               double   *ei, /* Expected intensities */
               double   *y,  /* Data (observed intensities). */
               double   es2, /* The estimate of the variance.*/ 
               int      n,   /* Length of the data. */
               /* Output. */
               double   *z);             

double make_dfdphi(/* Input */
                   double   *g,
                   double   *L,
                   double   *R,
                   double   *dRdphi,
                   double   b);

double make_dfdl(/* Input */
                 double   *g,
                 double   dRdl,
                 double   *R,
                 double   b,
                 int      dRdl_i);

double log_like_ric(/* Input */
                    double  *theta, /* Parameters. */
                    double  *y,     /* Data. */
                    double  *z,    /* z_l as defined in paper. */
                    double  *ei,    /* Expected intensities. */
                    int     n);     /* Number of data points. */

double log_like_gauss(/* Input */
                      double  *theta, /* Parameters. */
                      double  *y,     /* Data. */
                      double  *ei,    /* Expected intensities. */
                      int     n);     /* Number of data points. */

double *make_ei(/* Input */
                double   *theta,
                double   *fl,
                int      n,
                /* Output */
                double   *ei);

void make_ei_ef_e2f(/* Input */
                    double   *theta,
                    double   *fl,
                    int      n,
                    /* Output */
                    double   *ei,
                    double   *ef,
                    double   *e2f);

double logI0(double  x);

double logI1(double  x);

double make_I1I0(double   x);

double make_tricky(double  A,
                   double  s2);

void fix_phi2(double   *theta,
              int      nvox);

} // end namespace NLDTIFIT
