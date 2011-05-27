// Nonlinear estimation of diffusion tensor
// using a Gaussian or Rician noise model
//
// nldtifit 
//
// Jesper Andersson & Stam Rastapopoulos, FMRIB Image Analysis Group
//
// Copyright (C) 2011 University of Oxford 
//

#include <cstdlib>
#include <iostream>
#include <string>
#include <ctime>
#include <cmath>
#include "newmat.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "RLR_MAP_routines.h"
#include "nldtifit.h"

using namespace std;
using namespace NEWMAT;
using namespace NEWIMAGE;

int main(int   argc,
         char  *argv[])
{
  // To be written by Stam
}

namespace NLDTIFIT {

NEWIMAGE::volume4D<float> nldtifit(// Data, b=0 and diffusion weighted
                                   const NEWIMAGE::volume4D<float>&   scans,
                                   // Initial guesses of the parameters in the order s2, s0, R1, R2, R3, L1, L2, L3.
                                   // N.B. that s2, s0 and L* must all be positive
                                   const NEWIMAGE::volume4D<float>&   inpar,
                                   // Mask indicating where we want calculations to be performed
                                   const NEWIMAGE::volume<char>&      mask,
                                   // nx3 matrix with gradient directions in same order as in scans
                                   const NEWMAT::Matrix&              bvec,
                                   // nx1 vector with b-values in same order as in scans. N.B. that wa can have many different b-values
                                   const NEWMAT::ColumnVector&        bval,
                                   // Indicates what neighbourhood of voxels we want to base s2 estimate on.
                                   NbHood                             nbhood,
                                   // Gaussian or Rician
                                   NoiseModel                         nm)
{
  verify_input(scans,inpar,mask,bvec,bval,nbhood,nm);

  NEWIMAGE::volume4D<float>  ovol = inpar;
  ovol = 0.0;
  int n = scans.tsize();
  double *y = new double[nbhood*n];
  double *theta = new double[nbhood*7+1];
  double *g = new double[3*n];
  double *b = new double[n];
  repack_g_and_b(bvec,bval,g,b);
  int    np = 3;
  double pmu[3] = {-7.29, -7.30, -7.31};
  double pvar[3] = {4.0, 4.0, 4.0};
  int    pi[3] = {5, 6, 7};
  double *ei = new double[nbhood*n];
  double hist[MAXITER];

  for (int k=0; k<scans.zsize(); k++) {
    for (int j=0; j<scans.ysize(); j++) {
      for (int i=0; i<scans.xsize(); i++) {
        if (mask(i,j,k)) {
          int nvox = load_data_and_initial_guesses(scans,inpar,mask,i,j,k,nbhood,y,theta);
          if (fit_RLR(g,b,y,n,nvox,pmu,pvar,pi,np,nm,theta,ei,hist) > 0) {
            deal_results(theta,i,j,k,ovol);
	  }
	}
      }
    }
  }

  delete[] y; 
  delete[] theta;
  delete[] g; 
  delete[] b; 
  delete[] ei; 

  return(ovol); 
}

void verify_input(const NEWIMAGE::volume4D<float>&   scans,
                  const NEWIMAGE::volume4D<float>&   inpar,
                  const NEWIMAGE::volume<char>&      mask,
                  const NEWMAT::Matrix&              bvec,
                  const NEWMAT::ColumnVector&        bval,
                  NbHood                             nbhood,
                  NoiseModel                         nm)
{
  int n = scans.tsize();
  if (!samesize(scans[0],inpar[0])) throw nldtifitException("verify_input: Data and intial guesses must have same dimensions");
  if (!samesize(scans[0],mask)) throw nldtifitException("verify_input: Data and mask must have same dimensions");
  if (inpar.tsize() != 8) throw nldtifitException("verify_input: There must be initial guesses for eight parameters per voxel");
  if (bvec.Nrows() != n) throw nldtifitException("verify_input: The number of gradient vectors must match the number of scans");
  if (bvec.Ncols() != 3) throw nldtifitException("verify_input: Each gradient vector must have three elements");
  if (bval.Nrows() != n) throw nldtifitException("verify_input: The number of b-values must match the number of scans");
  if (nbhood > NONE && nm == GAUSSIAN) throw nldtifitException("verify_input: You cannot combine a neighbourhood for s2 with Gaussian noise model");
}

void deal_results(// Input
                  const double                      *theta,
                  int                               ii,
                  int                               jj,
                  int                               kk,
                  // Output
                  NEWIMAGE::volume4D<float>&        ovol)
{
  for (int i=0; i<8; i++) ovol(ii,jj,kk,i) = theta[i];
}

void repack_g_and_b(// Input 
                    const NEWMAT::Matrix&         bvec,
                    const NEWMAT::ColumnVector&   bval,
		    // Output
		    double                        *g,
                    double                        *b)
{
  for (int r=1; r<=bvec.Nrows(); r++) {
    *b++ = bval(r);
    for (int c=1; c<=bvec.Ncols(); c++) {
      *g++ = bvec(r,c);
    }
  }
  return;
}

int load_data_and_initial_guesses(// Input
                                  const NEWIMAGE::volume4D<float>&  scans,
                                  const NEWIMAGE::volume4D<float>&  inpar,
                                  const NEWIMAGE::volume<char>&     mask,
                                  int                               ii,
                                  int                               jj,
                                  int                               kk,
                                  NbHood                            nbhood,
                                  // Output
                                  double                            *y,
                                  double                            *ig)
{
  int nvox = 0;
  if (mask(ii,jj,kk)) { // This SHOULD be true
    y = load_data_single_voxel(scans,ii,jj,kk,y);
    *ig++ = inpar(ii,jj,kk,0);
    ig = load_guess_single_voxel(inpar,ii,jj,kk,ig);
    nvox++;
    if (nbhood > NONE) { // Add face voxels
      if (ii && mask(ii-1,jj,kk)) {y = load_data_single_voxel(scans,ii-1,jj,kk,y); ig = load_guess_single_voxel(inpar,ii-1,jj,kk,ig); nvox++;}
      if (jj && mask(ii,jj-1,kk)) {y = load_data_single_voxel(scans,ii,jj-1,kk,y); ig = load_guess_single_voxel(inpar,ii,jj-1,kk,ig); nvox++;}
      if (kk && mask(ii,jj,kk-1)) {y = load_data_single_voxel(scans,ii,jj,kk-1,y); ig = load_guess_single_voxel(inpar,ii,jj,kk-1,ig); nvox++;}
      if (ii<scans.xsize()-1 && mask(ii+1,jj,kk)) {y = load_data_single_voxel(scans,ii+1,jj,kk,y); ig = load_guess_single_voxel(inpar,ii+1,jj,kk,ig); nvox++;}
      if (jj<scans.ysize()-1 && mask(ii,jj+1,kk)) {y = load_data_single_voxel(scans,ii,jj+1,kk,y); ig = load_guess_single_voxel(inpar,ii,jj+1,kk,ig); nvox++;}
      if (kk<scans.zsize()-1 && mask(ii,jj,kk+1)) {y = load_data_single_voxel(scans,ii,jj,kk+1,y); ig = load_guess_single_voxel(inpar,ii,jj,kk+1,ig); nvox++;}
    }
    if (nbhood > FACE) { // Add edge voxels
      if (ii && jj && mask(ii-1,jj-1,kk)) {y = load_data_single_voxel(scans,ii-1,jj-1,kk,y); ig = load_guess_single_voxel(inpar,ii-1,jj-1,kk,ig); nvox++;}
      if (ii && kk && mask(ii-1,jj,kk-1)) {y = load_data_single_voxel(scans,ii-1,jj,kk-1,y); ig = load_guess_single_voxel(inpar,ii-1,jj,kk-1,ig); nvox++;}
      if (ii && jj<scans.ysize()-1 && mask(ii-1,jj+1,kk)) {y = load_data_single_voxel(scans,ii-1,jj+1,kk,y); ig = load_guess_single_voxel(inpar,ii-1,jj+1,kk,ig); nvox++;}
      if (ii && kk<scans.zsize()-1 && mask(ii-1,jj,kk+1)) {y = load_data_single_voxel(scans,ii-1,jj,kk+1,y); ig = load_guess_single_voxel(inpar,ii-1,jj,kk+1,ig); nvox++;}
      if (ii<scans.xsize()-1 && jj && mask(ii+1,jj-1,kk)) {y = load_data_single_voxel(scans,ii+1,jj-1,kk,y); ig = load_guess_single_voxel(inpar,ii+1,jj-1,kk,ig); nvox++;}
      if (ii<scans.xsize()-1 && kk && mask(ii+1,jj,kk-1)) {y = load_data_single_voxel(scans,ii+1,jj,kk-1,y); ig = load_guess_single_voxel(inpar,ii+1,jj,kk-1,ig); nvox++;}
      if (ii<scans.xsize()-1 && jj<scans.ysize()-1 && mask(ii+1,jj+1,kk)) {y = load_data_single_voxel(scans,ii+1,jj+1,kk,y); ig = load_guess_single_voxel(inpar,ii+1,jj+1,kk,ig); nvox++;}
      if (ii<scans.xsize()-1 && kk<scans.zsize()-1 && mask(ii+1,jj,kk+1)) {y = load_data_single_voxel(scans,ii+1,jj,kk+1,y); ig = load_guess_single_voxel(inpar,ii+1,jj,kk+1,ig); nvox++;}
      if (jj && kk && mask(ii,jj-1,kk-1)) {y = load_data_single_voxel(scans,ii,jj-1,kk-1,y); ig = load_guess_single_voxel(inpar,ii,jj-1,kk-1,ig); nvox++;}
      if (jj && kk<scans.zsize()-1 && mask(ii,jj-1,kk+1)) {y = load_data_single_voxel(scans,ii,jj-1,kk+1,y); ig = load_guess_single_voxel(inpar,ii,jj-1,kk+1,ig); nvox++;}
      if (jj<scans.ysize()-1 && kk && mask(ii,jj+1,kk-1)) {y = load_data_single_voxel(scans,ii,jj+1,kk-1,y); ig = load_guess_single_voxel(inpar,ii,jj+1,kk-1,ig); nvox++;}
      if (jj<scans.ysize()-1 && kk<scans.zsize()-1 && mask(ii,jj+1,kk+1)) {y = load_data_single_voxel(scans,ii,jj+1,kk+1,y); ig = load_guess_single_voxel(inpar,ii,jj+1,kk+1,ig); nvox++;}
    }
    if (nbhood > EDGE) {
      if (ii && jj && kk && mask(ii-1,jj-1,kk-1)) {y = load_data_single_voxel(scans,ii-1,jj-1,kk-1,y); ig = load_guess_single_voxel(inpar,ii-1,jj-1,kk-1,ig); nvox++;}
      if (ii<scans.xsize()-1 && jj && kk && mask(ii+1,jj-1,kk-1)) {y = load_data_single_voxel(scans,ii+1,jj-1,kk-1,y); ig = load_guess_single_voxel(inpar,ii+1,jj-1,kk-1,ig); nvox++;}
      if (ii && jj<scans.ysize()-1 && kk && mask(ii-1,jj+1,kk-1)) {y = load_data_single_voxel(scans,ii-1,jj+1,kk-1,y); ig = load_guess_single_voxel(inpar,ii-1,jj+1,kk-1,ig); nvox++;}
      if (ii && jj && kk<scans.zsize()-1 && mask(ii-1,jj-1,kk+1)) {y = load_data_single_voxel(scans,ii-1,jj-1,kk+1,y); ig = load_guess_single_voxel(inpar,ii-1,jj-1,kk+1,ig); nvox++;}
      if (ii<scans.xsize()-1 && jj<scans.ysize()-1 && kk && mask(ii+1,jj+1,kk-1)) {y = load_data_single_voxel(scans,ii+1,jj+1,kk-1,y); ig = load_guess_single_voxel(inpar,ii+1,jj+1,kk-1,ig); nvox++;}
      if (ii<scans.xsize()-1 && jj && kk<scans.zsize()-1 && mask(ii+1,jj-1,kk+1)) {y = load_data_single_voxel(scans,ii+1,jj-1,kk+1,y); ig = load_guess_single_voxel(inpar,ii+1,jj-1,kk+1,ig); nvox++;}
      if (ii && jj<scans.ysize()-1 && kk<scans.zsize()-1 && mask(ii-1,jj+1,kk+1)) {y = load_data_single_voxel(scans,ii-1,jj+1,kk+1,y); ig = load_guess_single_voxel(inpar,ii-1,jj+1,kk+1,ig); nvox++;}
      if (ii<scans.xsize()-1 && jj<scans.ysize()-1 && kk<scans.zsize()-1 && mask(ii+1,jj+1,kk+1)) {y = load_data_single_voxel(scans,ii+1,jj+1,kk+1,y); ig = load_guess_single_voxel(inpar,ii+1,jj+1,kk+1,ig); nvox++;}
    }
  }
  return(nvox);
}

double *load_data_single_voxel(const NEWIMAGE::volume4D<float>&  scans,
                               int                               ii,
                               int                               jj,
                               int                               kk,
                               // Output
                               double                            *y)
{
  for (int i = 0; i<scans.tsize(); i++) *y++ = scans(ii,jj,kk,i);
  return(y);
}

double *load_guess_single_voxel(const NEWIMAGE::volume4D<float>&  inpar,
                                int                               ii,
                                int                               jj,
                                int                               kk,
                                // Output
                                double                            *ig)
{
  for (int i = 1; i<inpar.tsize(); i++) *ig++ = inpar(ii,jj,kk,i); // N.B. starting at 1 (after s2)
  return(ig);
}

} // end namespace NLDTIFIT
