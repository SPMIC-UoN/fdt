/* Xfibres Diffusion Partial Volume Model  

    Tim Behrens - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <strstream>
#define WANT_STREAM
#define WANT_MATH
//  #include "newmatap.h"
//  #include "newmatio.h"
#include <string>
#include <math.h>
#include "utils/log.h"
#include "diff_pvmoptions.h"
#include "utils/tracer_plus.h"
#include "miscmaths/miscprob.h"
#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"
#include "stdlib.h"
#include "fibre.h"
#include "xfibresoptions.h"

//#include "bint/model.h"
//#include "bint/lsmcmcmanager.h"
//#include "bint/lslaplacemanager.h"

using namespace FIBRE;
using namespace  Xfibres;
using namespace Utilities;
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace MISCMATHS;



const float maxfloat=1e10;
const float minfloat=1e-10;
const float maxlogfloat=23;
const float minlogfloat=-23;


inline float min(float a,float b){
  return a<b ? a:b;}
inline float max(float a,float b){
  return a>b ? a:b;}
inline Matrix Anis()
{ 
  Matrix A(3,3);
  A << 1 << 0 << 0
    << 0 << 0 << 0
    << 0 << 0 << 0;
  return A;
}

inline Matrix Is()
{ 
  Matrix I(3,3);
  I << 1 << 0 << 0
    << 0 << 1 << 0
    << 0 << 0 << 1;
  return I;
}

inline ColumnVector Cross(const ColumnVector& A,const ColumnVector& B)
{
  ColumnVector res(3);
  res << A(2)*B(3)-A(3)*B(2)
      << A(3)*B(1)-A(1)*B(3)
      << A(1)*B(2)-B(1)*A(2);
  return res;
}

inline Matrix Cross(const Matrix& A,const Matrix& B)
{
  Matrix res(3,1);
  res << A(2,1)*B(3,1)-A(3,1)*B(2,1)
      << A(3,1)*B(1,1)-A(1,1)*B(3,1)
      << A(1,1)*B(2,1)-B(1,1)*A(2,1);
  return res;
}

float mod(float a, float b){
  while(a>b){a=a-b;}
  while(a<0){a=a+b;} 
  return a;
}


Matrix form_Amat(const Matrix& r,const Matrix& b)
{
  Matrix A(r.Ncols(),7);
  Matrix tmpvec(3,1), tmpmat;
  
  for( int i = 1; i <= r.Ncols(); i++){
    tmpvec << r(1,i) << r(2,i) << r(3,i);
    tmpmat = tmpvec*tmpvec.t()*b(1,i);
    A(i,1) = tmpmat(1,1);
    A(i,2) = 2*tmpmat(1,2);
    A(i,3) = 2*tmpmat(1,3);
    A(i,4) = tmpmat(2,2);
    A(i,5) = 2*tmpmat(2,3);
    A(i,6) = tmpmat(3,3);
    A(i,7) = 1;
  }
  return A;
}

inline SymmetricMatrix vec2tens(ColumnVector& Vec){
  SymmetricMatrix tens(3);
  tens(1,1)=Vec(1);
  tens(2,1)=Vec(2);
  tens(3,1)=Vec(3);
  tens(2,2)=Vec(4);
  tens(3,2)=Vec(5);
  tens(3,3)=Vec(6);
  return tens;
}



class Samples{
  xfibresOptions& opts;
  Matrix m_dsamples;
  Matrix m_S0samples;
  vector<Matrix> m_thsamples;
  vector<Matrix> m_phsamples;
  vector<Matrix> m_fsamples;
  vector<Matrix> m_lamsamples;
public:

  Samples(int nvoxels):opts(xfibresOptions::getInstance()){
    int count=0;
    int nsamples=0;
  
    for(int i=0;i<opts.njumps.value();i++){
      count++;
      if(count==opts.sampleevery.value()){
	count=0;nsamples++;
      }
    }
 

    m_dsamples.ReSize(nsamples,nvoxels);
    m_dsamples=0;
    m_S0samples.ReSize(nsamples,nvoxels);
    m_S0samples=0;
   
    for(int f=0;f<opts.nfibres.value();f++){
      m_thsamples.push_back(m_S0samples);
      m_phsamples.push_back(m_S0samples);
      m_fsamples.push_back(m_S0samples);
      m_lamsamples.push_back(m_S0samples);

    }
 
  }
  
  
  void record(Multifibre& mfib, int vox, int samp){
    m_dsamples(samp,vox)=mfib.get_d();
    m_S0samples(samp,vox)=mfib.get_S0();
    for(int f=0;f<opts.nfibres.value();f++){
      m_thsamples[f](samp,vox)=mfib.fibres()[f].get_th();
      m_phsamples[f](samp,vox)=mfib.fibres()[f].get_ph();
      m_fsamples[f](samp,vox)=mfib.fibres()[f].get_f();
      m_lamsamples[f](samp,vox)=mfib.fibres()[f].get_lam();
      
    }
  }
  
  void save(const volume<float>& mask){
    volume4D<float> tmp;
    //So that I can sort the output fibres into
    // files ordered by fibre fractional volume..
    vector<Matrix> thsamples_out=m_thsamples;
    vector<Matrix> phsamples_out=m_phsamples;
    vector<Matrix> fsamples_out=m_fsamples;
    vector<Matrix> lamsamples_out=m_lamsamples;

    Log& logger = LogSingleton::getInstance();
    tmp.setmatrix(m_dsamples,mask);
    save_volume4D(tmp,logger.appendDir("dsamples"));
    tmp.setmatrix(m_S0samples,mask);
    save_volume4D(tmp,logger.appendDir("S0samples"));
    
    //Sort the output based on mean_fsamples
    // 
    vector<Matrix> sumf;
    cerr<<"cock"<<endl;
    for(int f=0;f<opts.nfibres.value();f++){
      Matrix tmp=sum(m_fsamples[f],1);
      sumf.push_back(tmp);
    }  
    cerr<<"cock2"<<endl;
    for(int vox=1;vox<=m_dsamples.Ncols();vox++){
      vector<pair<float,int> > sfs;
      pair<float,int> ftmp;
      
      for(int f=0;f<opts.nfibres.value();f++){
	ftmp.first=sumf[f](1,vox);
	ftmp.second=f;
	sfs.push_back(ftmp);
      }
      sort(sfs.begin(),sfs.end());
      
      for(int samp=1;samp<=m_dsamples.Nrows();samp++){
	for(int f=0;f<opts.nfibres.value();f++){
 	  thsamples_out[f](samp,vox)=m_thsamples[sfs[(sfs.size()-1)-f].second](samp,vox);
	  phsamples_out[f](samp,vox)=m_phsamples[sfs[(sfs.size()-1)-f].second](samp,vox);
	  fsamples_out[f](samp,vox)=m_fsamples[sfs[(sfs.size()-1)-f].second](samp,vox);
	  lamsamples_out[f](samp,vox)=m_lamsamples[sfs[(sfs.size()-1)-f].second](samp,vox);
	}
      
      }
      
      
    }
    cerr<<"cock3"<<endl;
    // save the sorted fibres
    for(int f=0;f<opts.nfibres.value();f++){
      //      element_mod_n(thsamples_out[f],M_PI);
      //      element_mod_n(phsamples_out[f],2*M_PI);
      tmp.setmatrix(thsamples_out[f],mask);
      string oname="th"+num2str(f+1)+"samples";
      save_volume4D(tmp,logger.appendDir(oname));
      tmp.setmatrix(phsamples_out[f],mask);
      oname="ph"+num2str(f+1)+"samples";
      save_volume4D(tmp,logger.appendDir(oname));
      tmp.setmatrix(fsamples_out[f],mask);
      oname="f"+num2str(f+1)+"samples";
      save_volume4D(tmp,logger.appendDir(oname));
      tmp.setmatrix(lamsamples_out[f],mask);
      oname="lam"+num2str(f+1)+"samples";
      save_volume4D(tmp,logger.appendDir(oname));
    }
  }
  
};











class xfibresVoxelManager{
 
  xfibresOptions& opts;
  
  Samples& m_samples;
  int m_voxelnumber;
  const ColumnVector m_data;
  const ColumnVector& m_alpha;
  const ColumnVector& m_beta;
  const Matrix& m_bvals; 
  Multifibre m_multifibre;
 public:
  xfibresVoxelManager(const ColumnVector& data,const ColumnVector& alpha, 
		      const ColumnVector& beta, const Matrix& b,
		      Samples& samples,int voxelnumber):
    opts(xfibresOptions::getInstance()), 
    m_samples(samples),m_voxelnumber(voxelnumber),m_data(data), 
    m_alpha(alpha), m_beta(beta), m_bvals(b), 
    m_multifibre(m_data,m_alpha,m_beta,m_bvals,opts.nfibres.value()){ }
  
   
  
  void initialise(const Matrix& Amat){
    //initialising 
    ColumnVector logS(m_data.Nrows()),tmp(m_data.Nrows()),Dvec(7),dir(3);
    SymmetricMatrix tens;   
    DiagonalMatrix Dd;  
    Matrix Vd;  
    float mDd,fsquared;
    float th,ph,f,D,S0;
    for ( int i = 1; i <= logS.Nrows(); i++)
      {
	if(m_data(i)>0){
	  logS(i)=log(m_data(i));
	}
	else{
	  logS(i)=0;
	}
      }
    Dvec = -pinv(Amat)*logS;
    if(  Dvec(7) >  -maxlogfloat  ){ 
      S0=exp(-Dvec(7));
    }
    else{
      S0=m_data.MaximumAbsoluteValue();
    }

    for ( int i = 1; i <= logS.Nrows(); i++)
      {
	if(S0<m_data.Sum()/m_data.Nrows()){ S0=m_data.MaximumAbsoluteValue();  }
	logS(i)=(m_data(i)/S0)>0.01 ? log(m_data(i)):log(0.01*S0);
      }

    Dvec = -pinv(Amat)*logS;
    S0=exp(-Dvec(7));
    if(S0<m_data.Sum()/m_data.Nrows()){ S0=m_data.Sum()/m_data.Nrows();  }
    tens = vec2tens(Dvec);
    EigenValues(tens,Dd,Vd);
    mDd = Dd.Sum()/Dd.Nrows();
    int maxind = Dd(1) > Dd(2) ? 1:2;   //finding maximum eigenvalue
    maxind = Dd(maxind) > Dd(3) ? maxind:3;
    dir << Vd(1,maxind) << Vd(2,maxind) << Vd(3,maxind);
    cart2sph(dir,th,ph);
    th= mod(th,M_PI);
    ph= mod(ph,2*M_PI);
    D = Dd(maxind);

    float numer=1.5*((Dd(1)-mDd)*(Dd(1)-mDd)+(Dd(2)-mDd)*(Dd(2)-mDd)+(Dd(3)-mDd)*(Dd(3)-mDd));
    float denom=(Dd(1)*Dd(1)+Dd(2)*Dd(2)+Dd(3)*Dd(3));
    if(denom>0) fsquared=numer/denom;
    else fsquared=0;
    if(fsquared>0){f=sqrt(fsquared);}
    else{f=0;}
    if(f>=0.95) f=0.95;
    if(f<=0.001) f=0.001;
    if(D<=0) D=2e-3;
    m_multifibre.set_d(D);
    m_multifibre.set_S0(S0);
    if(opts.nfibres.value()>0){
      m_multifibre.addfibre(th,ph,f,1);
      for(int i=2; i<=opts.nfibres.value(); i++){
	 m_multifibre.addfibre();
      }
    
    }
    m_multifibre.initialise_energies();
    m_multifibre.initialise_props();
  }
 

  void runmcmc(){
    int count=0, recordcount=0,sample=1;//sample will index a newmat matrix 
    for( int i =0;i<opts.nburn.value();i++){
      m_multifibre.jump();
      count++;
      if(count==opts.updateproposalevery.value()){
	m_multifibre.update_proposals();
	count=0;
      }
    }
    
    for( int i =0;i<opts.njumps.value();i++){
      m_multifibre.jump();
      count++;
      recordcount++;
      if(recordcount==opts.sampleevery.value()){
	m_samples.record(m_multifibre,m_voxelnumber,sample);
	sample++;
	recordcount=0;
      }
      if(count==opts.updateproposalevery.value()){
	m_multifibre.update_proposals();
	count=0;
	
      }
    }
    
    
  }
    
};


  
int main(int argc, char *argv[])
{
  try{  

    // Setup logging:
    Log& logger = LogSingleton::getInstance();
    cerr<<"here"<<endl;
    xfibresOptions& opts = xfibresOptions::getInstance();
    opts.parse_command_line(argc,argv,logger);
    srand(xfibresOptions::getInstance().seed.value());
    
    
    Matrix datam, bvals,bvecs;
    volume<float> mask;
    bvals=read_ascii_matrix(opts.bvalsfile.value());
    bvecs=read_ascii_matrix(opts.bvecsfile.value());
    
    {//scope in which the data exists in 4D format;
      volume4D<float> data;
      read_volume4D(data,opts.datafile.value());
      read_volume(mask,opts.maskfile.value());
      datam=data.matrix(mask);  
    }
    cerr<<"ok"<<endl;
    Matrix Amat;
    ColumnVector alpha, beta;
    Amat=form_Amat(bvecs,bvals);
    cart2sph(bvecs,alpha,beta);
    Samples samples(datam.Ncols());
    cerr<<"ok2"<<endl;
    cerr<<datam.Ncols()<<endl;
    for(int vox=1;vox<=datam.Ncols();vox++){
      cerr <<vox<<"/"<<datam.Ncols()<<endl;
      xfibresVoxelManager  vm(datam.Column(vox),alpha,beta,bvals,samples,vox);
      vm.initialise(Amat);
      vm.runmcmc();
    }
    
    samples.save(mask);

  }
  catch(Exception& e) 
    {
      cerr << endl << e.what() << endl;
    }
  catch(X_OptionError& e) 
    {
      cerr << endl << e.what() << endl;
    }

  return 0;
}
