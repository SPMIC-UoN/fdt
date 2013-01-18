/* 
   RubiX
   Stamatios Sotiropoulos, Saad Jbabdi and Tim Behrens  - FMRIB Image Analysis Group
   Copyright (C) 2012 University of Oxford  
   CCOPYRIGHT  
*/


#include <iostream>
#include <fstream>
#include <iomanip>
#define WANT_STREAM
#define WANT_MATH
#include <string>
#include <math.h>
#include "utils/log.h"
#include "utils/tracer_plus.h"
#include "miscmaths/miscprob.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/nonlin.h"
#include "newimage/newimageall.h"
#include "stdlib.h"
#include "rubixoptions.h"
#include "rubixvox.h"
#include "rubix.h"
#include "diffmodels.h"

using namespace RUBIX;
using namespace Utilities;
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace MISCMATHS;


////////////////////////////////////////////
//     HRSamples: MCMC SAMPLE STORAGE
////////////////////////////////////////////

//Store parameters for a certain sample at a certain HR voxel
void HRSamples::record(const HRvoxel& HRv, int vox, int samp){
  m_dsamples(samp,vox)=HRv.get_d();
  m_S0samples(samp,vox)=HRv.get_S0();
  if (m_rician)
    m_tausamples(samp,vox)=HRv.get_tau();
  if (m_modelnum==2)
    m_d_stdsamples(samp,vox)=HRv.get_d_std();

  for(int f=0; f<m_numfibres; f++){
    m_thsamples[f](samp,vox)=HRv.fibres()[f].get_th();
    m_phsamples[f](samp,vox)=HRv.fibres()[f].get_ph();
    m_fsamples[f](samp,vox)=HRv.fibres()[f].get_f();
  }
}


//Get the mean samples for a voxel once jumping has finished
void HRSamples::finish_voxel(int vox){

  m_mean_dsamples(vox)=m_dsamples.Column(vox).Sum()/m_nsamps;
  m_mean_S0samples(vox)=m_S0samples.Column(vox).Sum()/m_nsamps;

  if (m_rician)
    m_mean_tausamples(vox)=m_tausamples.Column(vox).Sum()/m_nsamps;
  if (m_modelnum==2)
    m_mean_d_stdsamples(vox)=m_d_stdsamples.Column(vox).Sum()/m_nsamps;
  
  for(int f=0; f<m_numfibres; f++){  //for each fibre of the voxel
    m_mean_fsamples[f](vox)=m_fsamples[f].Column(vox).Sum()/m_nsamps;

    SymmetricMatrix m_dyad(3); m_dyad=0;
    ColumnVector m_vec(3);
    DiagonalMatrix dyad_D; //eigenvalues
    Matrix dyad_V; //eigenvectors

    //Get the sum of dyadic tensors across samples
    for (int n=1; n<=m_nsamps; n++){
      float th=m_thsamples[f](n,vox);
      float ph=m_phsamples[f](n,vox);
      m_vec << sin(th)*cos(ph) << sin(th)*sin(ph)<<cos(th) ;
      m_dyad << m_dyad+m_vec*m_vec.t(); 
    }
    m_dyad=m_dyad/m_nsamps;

    //Eigendecompose the mean dyadic tensor to get the mean orientation across samples
    EigenValues(m_dyad,dyad_D,dyad_V);
    int maxeig;
    if(dyad_D(1)>dyad_D(2)){
      if(dyad_D(1)>dyad_D(3)) maxeig=1;
      else maxeig=3;
    }
    else{
      if(dyad_D(2)>dyad_D(3)) maxeig=2;
      else maxeig=3;
    }
    m_dyadic_vectors[f](1,vox)=dyad_V(1,maxeig);
    m_dyadic_vectors[f](2,vox)=dyad_V(2,maxeig);
    m_dyadic_vectors[f](3,vox)=dyad_V(3,maxeig);
  }
}


//save samples for all voxels
void HRSamples::save(const volume<float>& mask){
  volume4D<float> tmp;
  //So that I can sort the output fibres into
  // files ordered by fibre fractional volume..
  vector<Matrix> thsamples_out=m_thsamples;
  vector<Matrix> phsamples_out=m_phsamples;
  vector<Matrix> fsamples_out=m_fsamples;
    
  vector<Matrix> dyadic_vectors_out=m_dyadic_vectors;
  vector<Matrix> mean_fsamples_out;
  for(unsigned int f=0;f<m_mean_fsamples.size();f++)
    mean_fsamples_out.push_back(m_mean_fsamples[f]);

  Log& logger = LogSingleton::getInstance();
  tmp.setmatrix(m_mean_dsamples,mask);
  tmp.setDisplayMaximumMinimum(tmp.max(),0);
  save_volume(tmp[0],logger.appendDir("mean_dsamples"));
  
  tmp.setmatrix(m_mean_S0samples,mask);
  tmp.setDisplayMaximumMinimum(tmp.max(),0);
  save_volume(tmp[0],logger.appendDir("mean_S0samples"));
  
  if (m_rician){
    tmp.setmatrix(m_mean_tausamples,mask);
    tmp.setDisplayMaximumMinimum(tmp.max(),0);
    save_volume(tmp[0],logger.appendDir("mean_tausamples"));
  }

  if (m_modelnum==2){
    tmp.setmatrix(m_mean_d_stdsamples,mask);
    tmp.setDisplayMaximumMinimum(tmp.max(),0);
    save_volume(tmp[0],logger.appendDir("mean_dstd_samples"));
  }

  //Sort the output based on mean_fsamples
  vector<Matrix> sumf;
  for(int f=0; f<m_numfibres; f++){
    Matrix tmp=sum(m_fsamples[f],1);
    sumf.push_back(tmp);
  }  
  for(int vox=1;vox<=m_dsamples.Ncols();vox++){
    vector<pair<float,int> > sfs;
    pair<float,int> ftmp;
      
    for(int f=0; f<m_numfibres; f++){
      ftmp.first=sumf[f](1,vox);
      ftmp.second=f;
      sfs.push_back(ftmp);
    }
    sort(sfs.begin(),sfs.end());
      
    for(int samp=1;samp<=m_dsamples.Nrows();samp++){
      for(int f=0; f<m_numfibres; f++){;
	thsamples_out[f](samp,vox)=m_thsamples[sfs[(sfs.size()-1)-f].second](samp,vox);
	phsamples_out[f](samp,vox)=m_phsamples[sfs[(sfs.size()-1)-f].second](samp,vox);
	fsamples_out[f](samp,vox)=m_fsamples[sfs[(sfs.size()-1)-f].second](samp,vox);
      }
    }
      
    for(int f=0;f<m_numfibres;f++){
      mean_fsamples_out[f](1,vox)=m_mean_fsamples[sfs[(sfs.size()-1)-f].second](vox);
      dyadic_vectors_out[f](1,vox)=m_dyadic_vectors[sfs[(sfs.size()-1)-f].second](1,vox);
      dyadic_vectors_out[f](2,vox)=m_dyadic_vectors[sfs[(sfs.size()-1)-f].second](2,vox);
      dyadic_vectors_out[f](3,vox)=m_dyadic_vectors[sfs[(sfs.size()-1)-f].second](3,vox);
    }
  }
  
  tmp.setmatrix(m_dsamples,mask);
  tmp.setDisplayMaximumMinimum(tmp.max(),tmp.min());
  save_volume4D(tmp,logger.appendDir("d_samples"));

  // save the sorted fibres
  for(int f=0; f<m_numfibres; f++){
    tmp.setmatrix(thsamples_out[f],mask);
    tmp.setDisplayMaximumMinimum(tmp.max(),tmp.min());
    string oname="th"+num2str(f+1)+"samples";
    save_volume4D(tmp,logger.appendDir(oname));
      
    tmp.setmatrix(phsamples_out[f],mask);
    tmp.setDisplayMaximumMinimum(tmp.max(),tmp.min());
    oname="ph"+num2str(f+1)+"samples";
    save_volume4D(tmp,logger.appendDir(oname));
   
    tmp.setmatrix(fsamples_out[f],mask);
    tmp.setDisplayMaximumMinimum(1,0);
    oname="f"+num2str(f+1)+"samples";
    save_volume4D(tmp,logger.appendDir(oname));

    tmp.setmatrix(mean_fsamples_out[f],mask);
    tmp.setDisplayMaximumMinimum(1,0);
    oname="mean_f"+num2str(f+1)+"samples";
    save_volume(tmp[0],logger.appendDir(oname));
      
    tmp.setmatrix(dyadic_vectors_out[f],mask);
    tmp.setDisplayMaximumMinimum(1,-1);
    oname="dyads"+num2str(f+1);
    save_volume4D(tmp,logger.appendDir(oname));
  }
}




////////////////////////////////////////////
//     LRSamples: MCMC SAMPLE STORAGE
////////////////////////////////////////////
//Store parameters for a certain sample at a certain LR voxel
void LRSamples::record(const LRvoxel& LRv, int vox, int samp){
  m_S0samples(samp,vox)=LRv.get_S0LR();

  if (m_rician)
    m_tauLRsamples(samp,vox)=LRv.get_tauLR();
  if (m_fsumPrior)
    m_sumfsamples(samp,vox)=LRv.get_mean_fsum();
  if (m_dPrior)
    m_meandsamples(samp,vox)=LRv.get_meand();

  m_lik_energy(samp,vox)=LRv.get_likelihood_energy();
  m_prior_energy(samp,vox)=LRv.get_prior_energy();

  for(int m=0; m<m_Nmodes; m++){
    m_thsamples[m](samp,vox)=LRv.PModes()[m].get_th();
    m_phsamples[m](samp,vox)=LRv.PModes()[m].get_ph();
    m_ksamples[m](samp,vox)=1.0/LRv.PModes()[m].get_invkappa();
  }
}


//Get the mean samples for a voxel once jumping has finished
void LRSamples::finish_voxel(int vox){
  m_mean_S0samples(vox)=m_S0samples.Column(vox).Sum()/m_nsamps;
  if (m_rician)
    m_mean_tausamples(vox)=m_tauLRsamples.Column(vox).Sum()/m_nsamps;
  if (m_fsumPrior)
    m_mean_sumfsamples(vox)=m_sumfsamples.Column(vox).Sum()/m_nsamps;
  if (m_dPrior)
    m_mean_meandsamples(vox)=m_meandsamples.Column(vox).Sum()/m_nsamps;

  for(int m=0; m<m_Nmodes; m++){
    m_mean_ksamples[m](vox)=m_ksamples[m].Column(vox).Sum()/m_nsamps;

    DiagonalMatrix dyad_D; //eigenvalues
    Matrix dyad_V; //eigenvectors
    ColumnVector m_vec(3); 
    SymmetricMatrix m_dyad(3); m_dyad=0;
  
    for (int n=1; n<=m_nsamps; n++){
      float th=m_thsamples[m](n,vox);
      float ph=m_phsamples[m](n,vox);
      m_vec << sin(th)*cos(ph) << sin(th)*sin(ph)<<cos(th) ;
      m_dyad << m_dyad+m_vec*m_vec.t(); 
    }
    m_dyad=m_dyad/m_nsamps;

    //Eigendecompose the mean dyadic tensor to get the mean orientation across samples
    EigenValues(m_dyad,dyad_D,dyad_V);
    int maxeig;
    if(dyad_D(1)>dyad_D(2)){
      if(dyad_D(1)>dyad_D(3)) maxeig=1;
      else maxeig=3;
    }
    else{
      if(dyad_D(2)>dyad_D(3)) maxeig=2;
      else maxeig=3;
    }
    m_dyadic_vectors[m](1,vox)=dyad_V(1,maxeig);
    m_dyadic_vectors[m](2,vox)=dyad_V(2,maxeig);
    m_dyadic_vectors[m](3,vox)=dyad_V(3,maxeig);
  }
}


//save samples for all voxels
void LRSamples::save(const volume<float>& mask){
  volume4D<float> tmp;
  //So that I can sort the output fibres into
  // files ordered by dispersion..
  vector<Matrix> thsamples_out=m_thsamples;
  vector<Matrix> phsamples_out=m_phsamples;
  vector<Matrix> ksamples_out=m_ksamples;
    
  vector<Matrix> dyadic_vectors_out=m_dyadic_vectors;
  vector<Matrix> mean_ksamples_out;
  for(unsigned int f=0;f<m_mean_ksamples.size();f++)
    mean_ksamples_out.push_back(m_mean_ksamples[f]);

  Log& logger = LogSingleton::getInstance();
  tmp.setmatrix(m_mean_S0samples,mask);
  tmp.setDisplayMaximumMinimum(tmp.max(),0);
  save_volume(tmp[0],logger.appendDir("mean_S0_LRsamples"));
  
  //tmp.setmatrix(m_lik_energy,mask);
  //tmp.setDisplayMaximumMinimum(tmp.max(),0);
  //save_volume4D(tmp,logger.appendDir("En_Lik_LRsamples"));

  //tmp.setmatrix(m_prior_energy,mask);
  //tmp.setDisplayMaximumMinimum(tmp.max(),0);
  //save_volume4D(tmp,logger.appendDir("En_Prior_LRsamples"));

  if (m_rician){
    tmp.setmatrix(m_mean_tausamples,mask);
    tmp.setDisplayMaximumMinimum(tmp.max(),0);
    save_volume(tmp[0],logger.appendDir("mean_tau_LRsamples"));
  }
  if (m_fsumPrior){
    tmp.setmatrix(m_mean_sumfsamples,mask);
    tmp.setDisplayMaximumMinimum(1,0);
    save_volume(tmp[0],logger.appendDir("mean_fsumPriorMode_LRsamples"));
  }
  if (m_dPrior){
    tmp.setmatrix(m_mean_meandsamples,mask);
    tmp.setDisplayMaximumMinimum(tmp.max(),0);
    save_volume(tmp[0],logger.appendDir("mean_dPriorMode_LRsamples"));
  }


  //Sort the output based on mean_invksamples
  vector<Matrix> sumk;
  for(int f=0; f<m_Nmodes; f++){
    Matrix tmp=sum(m_ksamples[f],1);
    sumk.push_back(tmp);
  }  
  for(int vox=1;vox<=m_S0samples.Ncols();vox++){
    vector<pair<float,int> > sfs;
    pair<float,int> ftmp;
      
    for(int f=0; f<m_Nmodes; f++){
      ftmp.first=sumk[f](1,vox);
      ftmp.second=f;
      sfs.push_back(ftmp);
    }
    sort(sfs.begin(),sfs.end());
      
    for(int samp=1; samp<=m_nsamps; samp++){
      for(int f=0; f<m_Nmodes; f++){;
	thsamples_out[f](samp,vox)=m_thsamples[sfs[(sfs.size()-1)-f].second](samp,vox);
	phsamples_out[f](samp,vox)=m_phsamples[sfs[(sfs.size()-1)-f].second](samp,vox);
	ksamples_out[f](samp,vox)=m_ksamples[sfs[(sfs.size()-1)-f].second](samp,vox);
      }
    }
      
    for(int f=0;f<m_Nmodes;f++){
      mean_ksamples_out[f](1,vox)=m_mean_ksamples[sfs[(sfs.size()-1)-f].second](vox);
      dyadic_vectors_out[f](1,vox)=m_dyadic_vectors[sfs[(sfs.size()-1)-f].second](1,vox);
      dyadic_vectors_out[f](2,vox)=m_dyadic_vectors[sfs[(sfs.size()-1)-f].second](2,vox);
      dyadic_vectors_out[f](3,vox)=m_dyadic_vectors[sfs[(sfs.size()-1)-f].second](3,vox);
    }
  }
  
  // save the sorted fibres
  for(int f=0; f<m_Nmodes; f++){
    tmp.setmatrix(thsamples_out[f],mask);
    tmp.setDisplayMaximumMinimum(tmp.max(),tmp.min());
    string oname="Pth"+num2str(f+1)+"samples";
    save_volume4D(tmp,logger.appendDir(oname));
      
    tmp.setmatrix(phsamples_out[f],mask);
    tmp.setDisplayMaximumMinimum(tmp.max(),tmp.min());
    oname="Pph"+num2str(f+1)+"samples";
    save_volume4D(tmp,logger.appendDir(oname));
   
    tmp.setmatrix(ksamples_out[f],mask);
    tmp.setDisplayMaximumMinimum(1,0);
    oname="Pk"+num2str(f+1)+"samples";
    save_volume4D(tmp,logger.appendDir(oname));

    tmp.setmatrix(mean_ksamples_out[f],mask);
    tmp.setDisplayMaximumMinimum(1,0);
    oname="mean_Pk"+num2str(f+1)+"samples";
    save_volume(tmp[0],logger.appendDir(oname));
      
    tmp.setmatrix(dyadic_vectors_out[f],mask);
    tmp.setDisplayMaximumMinimum(1,-1);
    oname="Pdyads"+num2str(f+1);
    save_volume4D(tmp,logger.appendDir(oname));
  }
}



//////////////////////////////////////////////////////////////
//    LRVoxelManager:MCMC HANDLING for a single LR voxel
/////////////////////////////////////////////////////////////


void LRVoxelManager::initialise(){
  float pvmS0, pvmd,pvmd_std;
  ColumnVector pvmf,pvmth,pvmph,pvm2invk,pvm2th,pvm2ph,predicted_signal;
  
  //Initialise each HR voxel using the HR data
  float sumd=0, sumd2=0, sumf=0;
  for (int n=0; n<m_HRvoxnumber.Nrows(); n++){
    if (opts.modelnum.value()==1){ //Model 1

      PVM_single_c pvm(m_dataHR[n],m_bvecsHR,m_bvalsHR,opts.nfibres.value());
      pvm.fit(); // this will give th,ph,f in the correct order
      
      pvmf  = pvm.get_f();  pvmth = pvm.get_th(); pvmph = pvm.get_ph();
      pvmS0 = fabs(pvm.get_s0()); pvmd  = pvm.get_d();  predicted_signal=pvm.get_prediction();
      if(pvmd<0 || pvmd>0.01) pvmd=2e-3;

       // DTI dti1(m_dataHR[n],m_bvecsHR,m_bvalsHR);
      //dti1.linfit();
      //pvmS0=fabs(dti1.get_s0());

      m_LRv.set_HRparams(n,pvmd,pvmS0,pvmth,pvmph,pvmf);
    }
    else{  //Model 2
      PVM_multi pvm(m_dataHR[n],m_bvecsHR,m_bvalsHR,opts.nfibres.value());
      pvm.fit();
      
      pvmf  = pvm.get_f();  pvmth = pvm.get_th(); pvmph = pvm.get_ph(); pvmd_std=pvm.get_d_std();
      pvmS0 =fabs(pvm.get_s0()); pvmd  = pvm.get_d();  predicted_signal=pvm.get_prediction();
      if(pvmd<0 || pvmd>0.01) pvmd=2e-3;
      if(pvmd_std<0 || pvmd_std>0.01) pvmd_std=pvmd/10;

      //   DTI dti1(m_dataHR[n],m_bvecsHR,m_bvalsHR);
      //dti1.linfit();
      //pvmS0=fabs(dti1.get_s0());

      m_LRv.set_HRparams(n,pvmd,pvmd_std,pvmS0,pvmth,pvmph,pvmf); 
    } 
   
    if (opts.rician.value()){  //If using Rician Energy, initialize tau, using the variance of the initial fit residuals
      ColumnVector residuals;
      residuals=m_dataHR[n]-predicted_signal;
      float tau=1.0/var(residuals).AsScalar();
      m_LRv.set_tauHR(n,tau);
    }
    sumd+=pvmd;
    sumd2+=pvmd*pvmd;
    sumf+=pvmf.Sum();
  } 

  sumd/=m_HRvoxnumber.Nrows();
  m_LRv.set_meand(sumd);
  m_LRv.set_stdevd(sumd/100);//sqrt((sumd2-m_HRvoxnumber.Nrows()*sumd*sumd)/(m_HRvoxnumber.Nrows()-1.0)));
  
  m_LRv.set_mean_fsum(0.6);
  sumf/=m_HRvoxnumber.Nrows();
  //m_LRv.set_mean_fsum(sumf); //does not make a big difference compared to initialising with a constant sumf=0.6
  m_LRv.set_stdev_fsum(0.01);

  //Initialise the orientation prior parameters using the LR data
  if (opts.nmodes.value()>0){
    PVM_single_c pvm2(m_dataLR,m_bvecsLR,m_bvalsLR,opts.nmodes.value());
    pvm2.fit();
    
    pvm2invk  = pvm2.get_f();   
    pvm2th = pvm2.get_th();
    pvm2ph = pvm2.get_ph();
    for (int m=1; m<=pvm2th.Nrows(); m++)
      pvm2invk(m)=0.02;

    m_LRv.set_Priorparams(pvm2th,pvm2ph,pvm2invk);
  }

  //Initialise the rest LR params
  //DTI dti2(m_dataLR,m_bvecsLR,m_bvalsLR);
  //dti2.linfit();
  PVM_single_c pvm3(m_dataLR,m_bvecsLR,m_bvalsLR,1);
  pvm3.fit();  
  m_LRv.set_S0LR(fabs(pvm3.get_s0()));
    
  if (opts.rician.value()){  //If using Rician Energy, initialize tau, using the variance of the initial fit residuals
    ColumnVector residuals;
    predicted_signal=pvm3.get_prediction();
    residuals=m_dataLR-predicted_signal;
    float tau=1.0/var(residuals).AsScalar();
    m_LRv.set_tauLR(tau);
  } 
  
  m_LRv.update_Orient_hyp_prior();
  m_LRv.initialise_energies(); 
}
 


//Run MCMC for an LRvoxel
void LRVoxelManager::runmcmc(){
  int count=0, recordcount=0,sample=1;
 
  for(int i=0; i<opts.nburn.value(); i++){   //Burn-In
    m_LRv.jump();
    count++;
    if(count==opts.updateproposalevery.value()){
      m_LRv.update_proposals();
      count=0;
    }
  }
    
  for( int i =0;i<opts.njumps.value();i++){  //Sampling
    m_LRv.jump();
    count++;
    recordcount++;
    if(recordcount==opts.sampleevery.value()){
      for (int n=1; n<=m_HRvoxnumber.Nrows(); n++) //for each HR voxel 
	m_HRsamples.record(m_LRv.HRvoxels()[n-1],(int)m_HRvoxnumber(n),sample);
      m_LRsamples.record(m_LRv,(int)m_LRvoxnumber,sample);
      
      sample++;
      recordcount=0;
    }
    if(count==opts.updateproposalevery.value()){
      m_LRv.update_proposals();
      count=0;
    }
  }

  for (int n=1; n<=m_HRvoxnumber.Nrows(); n++)   //for each HR voxel 
    m_HRsamples.finish_voxel((int)m_HRvoxnumber(n)); 
  m_LRsamples.finish_voxel((int)m_LRvoxnumber);
      
}


//For a voxel (x,y,z) at Low_res, returns that overlapping HR voxels coordinates
ReturnMatrix get_HRindices(const int x,const int y,const int z, const int xratio, const int yratio, const int zratio){
  Matrix HRindices(xratio*yratio*zratio,3);  //num_HR_voxels x 3 coordinates per voxel

  int count=1;
  for (int Hx=1; Hx<=xratio; Hx++)
    for (int Hy=1; Hy<=yratio; Hy++)
      for (int Hz=1; Hz<=zratio; Hz++){
	HRindices(count,1)=(x+1)*xratio-Hx;
	HRindices(count,2)=(y+1)*yratio-Hy;
	HRindices(count,3)=(z+1)*zratio-Hz;
	count++;
      }

  HRindices.Release();
  return HRindices;
}


//Using a Low-Res mask, create an HR mask, after finding the corresponding HR voxels
volume<float> createHR_mask(const volume<float>& maskLR, const volume4D<float>& dataHR){
  volume<float> maskHR(dataHR.xsize(),dataHR.ysize(),dataHR.zsize());
  copybasicproperties(dataHR,maskHR);

  Matrix temp_indices;
  maskHR=0;
  float xratio=maskLR.xdim()/dataHR.xdim();
  float yratio=maskLR.ydim()/dataHR.ydim();
  float zratio=maskLR.zdim()/dataHR.zdim();

  for (int i=0; i<maskLR.zsize(); i++)
    for (int j=0; j<maskLR.ysize(); j++)
      for (int k=0; k<maskLR.xsize(); k++)
	if (maskLR(k,j,i)!=0){
	  temp_indices=get_HRindices(k,j,i,(int)xratio,(int)yratio,(int)zratio);
          for (int n=1; n<=temp_indices.Nrows(); n++)
	    maskHR.value((int)temp_indices(n,1),(int)temp_indices(n,2),(int)temp_indices(n,3))=1;
	}
  return maskHR;
}



////////////////////////////////////////////
//       MAIN
////////////////////////////////////////////
  
int main(int argc, char *argv[])
{
  try{  
    // Setup logging:
    Log& logger = LogSingleton::getInstance();
    rubixOptions& opts = rubixOptions::getInstance();
    opts.parse_command_line(argc,argv,logger);
    srand(rubixOptions::getInstance().seed.value());
    Matrix datamLR, bvalsLR,bvecsLR,matrix2volkeyLR;
    Matrix datamHR, bvalsHR,bvecsHR,matrix2volkeyHR;
    volume<float> maskLR, maskHR;
    volume<int> vol2matrixkeyLR, vol2matrixkeyHR;

    bvalsLR=read_ascii_matrix(opts.LRbvalsfile.value());
    bvecsLR=read_ascii_matrix(opts.LRbvecsfile.value());
    if(bvecsLR.Nrows()>3) bvecsLR=bvecsLR.t();
    if(bvalsLR.Nrows()>1) bvalsLR=bvalsLR.t();
    for(int i=1;i<=bvecsLR.Ncols();i++){
      float tmpsum=sqrt(bvecsLR(1,i)*bvecsLR(1,i)+bvecsLR(2,i)*bvecsLR(2,i)+bvecsLR(3,i)*bvecsLR(3,i));
      if(tmpsum!=0){
	bvecsLR(1,i)=bvecsLR(1,i)/tmpsum;
	bvecsLR(2,i)=bvecsLR(2,i)/tmpsum;
	bvecsLR(3,i)=bvecsLR(3,i)/tmpsum;
      }  
    }

    bvalsHR=read_ascii_matrix(opts.HRbvalsfile.value());
    bvecsHR=read_ascii_matrix(opts.HRbvecsfile.value());
    if(bvecsHR.Nrows()>3) bvecsHR=bvecsHR.t();
    if(bvalsHR.Nrows()>1) bvalsHR=bvalsHR.t();
    for(int i=1;i<=bvecsHR.Ncols();i++){
      float tmpsum=sqrt(bvecsHR(1,i)*bvecsHR(1,i)+bvecsHR(2,i)*bvecsHR(2,i)+bvecsHR(3,i)*bvecsHR(3,i));
      if(tmpsum!=0){
	bvecsHR(1,i)=bvecsHR(1,i)/tmpsum;
	bvecsHR(2,i)=bvecsHR(2,i)/tmpsum;
	bvecsHR(3,i)=bvecsHR(3,i)/tmpsum;
      }  
    }
    
    volume4D<float> dataLR,dataHR;
    read_volume4D(dataLR,opts.LRdatafile.value());
    read_volume(maskLR,opts.LRmaskfile.value());
    datamLR=dataLR.matrix(maskLR); 
    matrix2volkeyLR=dataLR.matrix2volkey(maskLR);
    vol2matrixkeyLR=dataLR.vol2matrixkey(maskLR);
 
    read_volume4D(dataHR,opts.HRdatafile.value());
    maskHR=createHR_mask(maskLR, dataHR);
    save_volume(maskHR,logger.appendDir("HRbrain_mask"));  //Generate and save an HR mask using the LR mask
    
    datamHR=dataHR.matrix(maskHR);
    matrix2volkeyHR=dataHR.matrix2volkey(maskHR);
    vol2matrixkeyHR=dataHR.vol2matrixkey(maskHR);
 
    float xratio=maskLR.xdim()/maskHR.xdim();
    float yratio=maskLR.ydim()/maskHR.ydim();
    float zratio=maskLR.zdim()/maskHR.zdim();

    HRSamples HRsampl(datamHR.Ncols(), opts.njumps.value(), opts.sampleevery.value(), opts.nfibres.value(), opts.rician.value(), opts.modelnum.value());
    LRSamples LRsampl(datamLR.Ncols(), opts.njumps.value(), opts.sampleevery.value(), opts.nmodes.value(), opts.rician.value(),opts.fsumPrior.value(),opts.dPrior.value());

    //dHR.push_back(datamHR.Column(1)); dHR.push_back(datamHR.Column(2));
    //dHR.push_back(datamHR.Column(3)); dHR.push_back(datamHR.Column(4));
    //HRvoxnum<<1<<2<<3<<4; HRweights<<0.25<<0.25<<0.25<<0.25;
    
    for(int vox=1;vox<=datamLR.Ncols();vox++){  //For each LR voxel
      //  if (vox==19){  //vox==6
      ColumnVector dLR; Matrix HRindices;
      dLR=datamLR.Column(vox);                  //Find the corresponding HR ones
      HRindices=get_HRindices((int)matrix2volkeyLR(vox,1),(int)matrix2volkeyLR(vox,2),(int)matrix2volkeyLR(vox,3),(int)xratio,(int)yratio,(int)zratio);
      //cout<<endl<<"S0LR: "<<dLR(1)<<endl<<endl; //Debugging code
     
      ColumnVector HRvoxnum(HRindices.Nrows()), HRweights(HRindices.Nrows()); 
      vector<ColumnVector> dHR;
      
     for (int n=1; n<=HRindices.Nrows(); n++){
	HRweights(n)=1.0/HRindices.Nrows();
	HRvoxnum(n)=vol2matrixkeyHR((int)HRindices(n,1),(int)HRindices(n,2),(int)HRindices(n,3));
	ColumnVector tmp_vec=datamHR.Column((int)HRvoxnum(n));
 	//cout<<(int)HRindices(n,1)<<" "<<(int)HRindices(n,2)<<" "<<(int)HRindices(n,3)<<"  S0HR:"<<tmp_vec(1)<<endl;
	dHR.push_back(datamHR.Column((int)HRvoxnum(n)));
      }
      cout <<vox<<"/"<<datamLR.Ncols()<<endl;
      LRVoxelManager  vm(HRsampl,LRsampl,vox,HRvoxnum, dLR,dHR, bvecsLR, bvalsLR, bvecsHR, bvalsHR, HRweights);
      vm.initialise();
      vm.runmcmc();
      //      } 
    }
    HRsampl.save(maskHR);
    LRsampl.save(maskLR);
  }

  catch(Exception& e){
    cerr << endl << e.what() << endl;
  }
  catch(X_OptionError& e){
      cerr << endl << e.what() << endl;
  }

  return 0;
}
