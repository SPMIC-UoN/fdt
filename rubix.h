/*  rubix.h: Classes utilized in RubiX MCMC storage and handling  */
/*  Stamatios Sotiropoulos, FMRIB Analysis Group */
/*  Copyright (C) 2012 University of Oxford  */
/*  CCOPYRIGHT  */


#if !defined(rubix_h)
#define rubix_h

#include <iostream>
#include <fstream>
#include <iomanip>
#define WANT_STREAM
#define WANT_MATH
#include <string>
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "newimage/newimageall.h"
#include "stdlib.h"

using namespace NEWMAT;
using namespace MISCMATHS;
using namespace NEWIMAGE;

namespace RUBIX{

  ////////////////////////////////////////////
  //       MCMC SAMPLE STORAGE
  ////////////////////////////////////////////
  //Storing Samples for parameters of all HR voxels 
  class HRSamples{
    Matrix m_dsamples;              //Variables for storing MCMC samples of all voxels in the HRgrid
    Matrix m_d_stdsamples;
    Matrix m_S0samples;
    Matrix m_tausamples;
    vector<Matrix> m_thsamples;
    vector<Matrix> m_phsamples;
    vector<Matrix> m_fsamples;
  
    //for storing means
    RowVector m_mean_dsamples;    //Storing mean_samples for all voxels in the HRgrid
    RowVector m_mean_d_stdsamples;
    RowVector m_mean_S0samples;
    RowVector m_mean_tausamples;
    vector<Matrix> m_dyadic_vectors;
    vector<RowVector> m_mean_fsamples;
  
    int m_nsamps;
    const int m_njumps;
    const int m_sample_every;
    const int m_numfibres;
    const bool m_rician;
    const int m_modelnum;
    //const string m_logdir;
    
  public:
  HRSamples(int nvoxels, const int njumps, const int sample_every, const int numfibres, const bool rician=false, const int modelnum=1):
    m_njumps(njumps),m_sample_every(sample_every), m_numfibres(numfibres), m_rician(rician), m_modelnum(modelnum){
      int count=0;
      int nsamples=0;
    
      for(int i=0;i<m_njumps; i++){
	count++;
	if(count==m_sample_every){
	  count=0; 
	  nsamples++;
	}
      }
      m_nsamps=nsamples;

      m_dsamples.ReSize(nsamples,nvoxels);  m_dsamples=0;
      m_S0samples.ReSize(nsamples,nvoxels); m_S0samples=0;

      m_mean_dsamples.ReSize(nvoxels);      m_mean_dsamples=0;
      m_mean_S0samples.ReSize(nvoxels);     m_mean_S0samples=0;
      if (m_rician){
	m_tausamples.ReSize(nsamples,nvoxels);  m_tausamples=0;
	m_mean_tausamples.ReSize(nvoxels);      m_mean_tausamples=0;
      }
      if (m_modelnum==2){
	m_d_stdsamples.ReSize(nsamples,nvoxels);  m_d_stdsamples=0;
	m_mean_d_stdsamples.ReSize(nvoxels);      m_mean_d_stdsamples=0;
      }
      Matrix tmpvecs(3,nvoxels);  tmpvecs=0;  
      for(int f=0; f<m_numfibres; f++){
	m_thsamples.push_back(m_S0samples);   
	m_phsamples.push_back(m_S0samples);
	m_fsamples.push_back(m_S0samples);  
	m_dyadic_vectors.push_back(tmpvecs); 
	m_mean_fsamples.push_back(m_mean_S0samples);
      }
    }

    ~HRSamples(){}

    void record(const HRvoxel& HRv, int vox, int samp); //Store parameters for a certain sample at a certain HR voxel
    void finish_voxel(int vox);                   //Get the mean samples for a voxel once jumping has finished
    void save(const volume<float>& mask);         //Save samples for all voxels
  };




  ////////////////////////////////////////////
  //       MCMC SAMPLE STORAGE
  ////////////////////////////////////////////
  //Storing Samples for parameters at the Low-Res level (e.g. priors, Low-res parameters)
  class LRSamples{
    vector<Matrix> m_thsamples;           //Variables for storing MCMC samples of all voxels in the LRgrid
    vector<Matrix> m_phsamples;
    vector<Matrix> m_ksamples;
    Matrix m_S0samples;
    Matrix m_tauLRsamples;
    Matrix m_lik_energy;
    Matrix m_prior_energy;

    //for storing means
    vector<Matrix> m_dyadic_vectors;
    vector<RowVector> m_mean_ksamples;
    RowVector m_mean_S0samples;
    RowVector m_mean_tausamples;

    int m_nsamps;
    const int m_njumps;
    const int m_sample_every;
    const int m_Nmodes;
    const bool m_rician;
    //const string m_logdir;
    
  public:
  LRSamples(int nvoxels, const int njumps, const int sample_every, const int Nmodes, const bool rician=false):
    m_njumps(njumps),m_sample_every(sample_every), m_Nmodes(Nmodes), m_rician(rician){
      int count=0;
      int nsamples=0;
    
      for(int i=0;i<m_njumps; i++){
	count++;
	if(count==m_sample_every){
	  count=0; 
	  nsamples++;
	}
      }
      m_nsamps=nsamples;

      m_S0samples.ReSize(nsamples,nvoxels); m_S0samples=0;
      m_lik_energy.ReSize(nsamples,nvoxels); m_lik_energy=0;
      m_prior_energy.ReSize(nsamples,nvoxels); m_prior_energy=0;
      m_mean_S0samples.ReSize(nvoxels);     m_mean_S0samples=0;
      if (m_rician){
	m_tauLRsamples.ReSize(nsamples,nvoxels); m_tauLRsamples=0;
	m_mean_tausamples.ReSize(nvoxels);     m_mean_tausamples=0;
      }
      Matrix tmpvecs(3,nvoxels);  tmpvecs=0;  
    
      for(int f=0; f<m_Nmodes; f++){
	m_thsamples.push_back(m_S0samples);   
	m_phsamples.push_back(m_S0samples);
	m_ksamples.push_back(m_S0samples);  
	m_dyadic_vectors.push_back(tmpvecs); 
	m_mean_ksamples.push_back(m_mean_S0samples);
      }
    }

    ~LRSamples(){}

    void record(const LRvoxel& LRv, int vox, int samp); //Store parameters for a certain sample at a certain LR voxel
    void finish_voxel(int vox);                   //Get the mean samples for a voxel once jumping has finished
    void save(const volume<float>& mask);         //Save samples for all voxels
  };



  //////////////////////////////////////////////
  //       MCMC HANDLING for a single LR voxel
  //////////////////////////////////////////////
  class LRVoxelManager{
    rubixOptions& opts;
    HRSamples& m_HRsamples;         //keep MCMC samples of the parameters of all voxels inferred at High-res grid 
    LRSamples& m_LRsamples;         //keep MCMC samples of the parameters of all voxels inferred at Low-res grid 
    int m_LRvoxnumber;
    ColumnVector m_HRvoxnumber;
    LRvoxel m_LRv;
    const ColumnVector& m_dataLR;    //Low-Res Data for the specific LR voxel 
    const vector<ColumnVector>& m_dataHR; //High-Res Data for all contained HR voxels
    const Matrix& m_bvecsLR;         //bvecs at Low-Res    (3 x LR_NumPoints)
    const Matrix& m_bvalsLR;         //bvalues at Low-Res  (1 x HR_NumPoints)
    const Matrix& m_bvecsHR;         //bvecs at High-Res   (3 x HR_NumPoints)
    const Matrix& m_bvalsHR;         //bvalues at High-Res (1 x HR_NumPoints)
    const ColumnVector& m_HRweights; //Holds the volume fraction each HR voxel occupies out of the LR one
  public:
    //Constructor
  LRVoxelManager(HRSamples& Hsamples, LRSamples& Lsamples, int LRvoxnum, ColumnVector& HRvoxnum, 
		   const ColumnVector& dataLR,const vector<ColumnVector>& dataHR, 
		 const Matrix& bvecsLR, const Matrix& bvalsLR, const Matrix& bvecsHR, const Matrix& bvalsHR, const ColumnVector& HRweights):
    opts(rubixOptions::getInstance()), m_HRsamples(Hsamples), m_LRsamples(Lsamples), m_LRvoxnumber(LRvoxnum),m_HRvoxnumber(HRvoxnum), 
      m_LRv(bvecsHR, bvalsHR, bvecsLR, bvalsLR, dataLR, dataHR, opts.nfibres.value(), opts.nmodes.value(), HRweights, opts.modelnum.value(), opts.fudge.value(),opts.all_ard.value(), opts.no_ard.value(),opts.kappa_ard.value(), opts.fsumPrior.value(), opts. dPrior.value(), opts.rician.value()),
      m_dataLR(dataLR), m_dataHR(dataHR),m_bvecsLR(bvecsLR), m_bvalsLR(bvalsLR), m_bvecsHR(bvecsHR), m_bvalsHR(bvalsHR), m_HRweights(HRweights) { } 
    
    ~LRVoxelManager() { }

    void initialise(); //Initialise all parameters for an LR voxel
    void runmcmc(); //Run MCMC for an LR voxel 
};


}


#endif
