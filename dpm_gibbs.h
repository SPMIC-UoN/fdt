#if !defined(_DPM_GIBBS_H)
#define _DPM_GIBBS_H

#include "gibbs.h"
#include "dpmOptions.h"
#include "newran/newran.h"
#include "miscmaths/miscprob.h"
#include <stdlib.h>
#include <stdio.h>
#include <cmath>


using namespace NEWMAT;
using namespace NEWRAN;
using namespace MISCMATHS;
using namespace DPM;
using namespace std;


// Gaussian-InverWishart distribution
// p(mu,sigma)=det(sigma)^(-(nu+d)/2-1)exp(-trace(Nu*inv(sigma))/2 -kappa/2*(mu-m_mu)'inv(sigma)(mu-m_mu))
class GaussianWishart{
  private:
  friend std::ostream& operator << (ostream& o,GaussianWishart& g);
 protected:
  ColumnVector       m_mu;
  SymmetricMatrix    m_Nu;
  float              m_kappa;
  int                m_dof;
  int                m_dim;
  
  ColumnVector       m_smu;     // sample mean
  SymmetricMatrix    m_ssigma;  // sample covariance

  
 public:
  GaussianWishart(){}
  GaussianWishart(const int dim):m_dim(dim){
    m_mu.ReSize(m_dim);
    m_Nu.ReSize(m_dim);
  }
  GaussianWishart(const ColumnVector& mu,const SymmetricMatrix& Nu,const int dof,const float& kappa):
    m_mu(mu),m_Nu(Nu),m_kappa(kappa),m_dof(dof){
    m_dim=m_mu.Nrows();
    sample();
  }
  ~GaussianWishart(){}
  inline ColumnVector get_mu()const{return m_mu;}
  void set_mu(const ColumnVector& mu){m_mu=mu;}
  inline SymmetricMatrix get_Nu()const{return m_Nu;}
  void set_Nu(const SymmetricMatrix& Nu){m_Nu=Nu;}
  void set_kappa(const float& kappa){m_kappa=kappa;}
  inline float get_kappa()const{return m_kappa;}
  inline int get_dof()const{return m_dof;}
  
  void postupdate(const vector<ColumnVector>& data,const GaussianWishart& gw0){
    ColumnVector mdat(m_dim);
    SymmetricMatrix S(m_dim),SS(m_dim);
    
    float n = (float)data.size();
    m_dof   = gw0.get_dof()   + int(n);
    m_kappa = gw0.get_kappa() + n;
    mdat=0;S=0,SS=0;
    for(int i=0;i<int(n);i++){
      SS << data[i]*data[i].t();
      S  += SS;
      mdat += data[i];
    }
    mdat /= n;
    
    SS << S -n*mdat*mdat.t();
    SS << SS + gw0.get_kappa()*n/m_kappa * (mdat-gw0.get_mu())*(mdat-gw0.get_mu()).t();
    
    m_mu    = ( gw0.get_kappa()*gw0.get_mu() + n*mdat )/m_kappa;
    m_Nu   << gw0.get_Nu() + SS;
    
    sample();
  }
  void postupdate(const Matrix& data,const GaussianWishart& gw0){
    ColumnVector mdat(m_dim);
    SymmetricMatrix S(m_dim),SS(m_dim);
    
    float n = (float)data.Nrows();
    m_dof   = gw0.get_dof()   + int(n);
    m_kappa = gw0.get_kappa() + n;
    mdat=0;S=0,SS=0;
    for(int i=1;i<=int(n);i++){
      SS << data.Row(i).t()*data.Row(i);
      S  += SS;
      mdat += data.Row(i).t();
    }
    mdat /= n;
    
    SS << S -n*mdat*mdat.t();
    SS << SS + gw0.get_kappa()*n/m_kappa * (mdat-gw0.get_mu())*(mdat-gw0.get_mu()).t();
    
    m_mu    = ( gw0.get_kappa()*gw0.get_mu() + n*mdat )/m_kappa;
    m_Nu   << gw0.get_Nu() + SS;
    
    sample();
  }
  void sample(ColumnVector& mu,SymmetricMatrix& sigma){
    sigma = iwishrnd(m_Nu.i(),m_dof);
    mu    = mvnrnd(m_mu.t(),sigma/m_kappa).t();
  }
  void sample(){
    m_ssigma = iwishrnd(m_Nu.i(),m_dof);
    m_smu    = mvnrnd(m_mu.t(),m_ssigma/m_kappa).t();
  }
  void print(ostream& os)const{ 
    os << "Gaussian-InverseWishart distribution" << endl;
    os << "mean       : " << m_mu.t();
    os << "variance   : " << m_Nu.Row(1);
    for(int i=2;i<=m_dim;i++)
      os << "             "<<m_Nu.Row(i);
    os << "dof        : "<<m_dof<<endl;
    os << "kappa      : "<<m_kappa<<endl;
    os << "sample mu  : "<<m_smu.t();
    os << "sample var : "<<m_ssigma.Row(1);
    for(int i=2;i<=m_dim;i++)
      os << "             "<<m_ssigma.Row(i);
    os << "-----------------------------------"<<endl;
  }
  ColumnVector get_smu()const{return m_smu;}
  SymmetricMatrix get_ssigma()const{return m_ssigma;}
  GaussianWishart& operator=(const GaussianWishart& rhs){
    m_mu     = rhs.m_mu;
    m_Nu     = rhs.m_Nu;
    m_kappa  = rhs.m_kappa;
    m_dof    = rhs.m_dof;
    m_dim    = rhs.m_dim;
    
    m_smu    = rhs.m_smu;
    m_ssigma = rhs.m_ssigma;

    return *this;
  }

};

//bool compare(const pair<int,float> &p1,const pair<int,float> &p2){
//return (p1.second < p2.second) ? true : false;
//}


class DPM_GibbsSampler : public GibbsSampler
{
 private:
  friend std::ostream& operator << (ostream& o,DPM_GibbsSampler& g);
  
 protected:
  DPM::dpmOptions& opts;

  // parameters            ------> estimated via gibb's sampling
  float                    m_alpha;
  vector<GaussianWishart>  m_gw;
  vector<int>              m_z;
  // hyperparameters       ------> estimated via gibb's sampling
  GaussianWishart          m_gw0;
  // hyper-hyperparameters ------> these are the only fixed parameters
  float                    m_a0;         // a0 = 1
  float                    m_b0;         // b0 = 1
  ColumnVector             m_m0;         // m0 = mean(data)
  SymmetricMatrix          m_S0;         // S0 = cov(data)
  SymmetricMatrix          m_N0;         // inv(cov(data))/(nu0-d-1)^2
  int                      m_n0;
  // data-related quantities
  vector<int>              m_classnd;
  int                      m_k;
  double                   m_margintbase;
  
  vector< pair<float,int> > randindex;

  // samples
  vector<float>            m_sample_alpha;
  vector<int>              m_sample_k;
  vector<double>           m_sample_likelihood;
  double                   m_likelihood;
  int                      m_nsamples;
  vector<float>            m_mean_z;

  const Matrix&            m_data;

public:
  DPM_GibbsSampler(const Matrix& data,int numiter,int burnin,int sampleevery):
    GibbsSampler(numiter,burnin,sampleevery),
      opts(DPM::dpmOptions::getInstance()),m_data(data){
      m_n = m_data.Nrows();
      m_d = m_data.Ncols();

    m_nsamples = (int)floor( (numiter - burnin) / sampleevery );

    m_sample_alpha.resize(m_nsamples);
    m_sample_k.resize(m_nsamples);
    m_sample_likelihood.resize(m_nsamples);
    m_mean_z.resize(m_n);

    Random::Set(rand() / float(RAND_MAX));
  }
  ~DPM_GibbsSampler(){}

  // parent class function definitions
  void sample_parameters();
  void sample_hyperparameters();

  // initialisation functions
  void init();
  void init_oneperdata();
  void init_onebigclass();
  void init_kmeans(const int k=10);
  void init_random(const int k=10);
  
  // sample model parameters
  void sample_z();
  void sample_gw();
  void sample_gw0();
  void sample_alpha();

  // utils
  double marglik(const ColumnVector&,const int);
  double margint(const ColumnVector&);
  void do_kmeans();
  ReturnMatrix get_dataindex();
  ReturnMatrix get_mldataindex();

  int get_numclass()const{return m_k;}

  // io
  void print(ostream& os){
    os << "-------fixed parameters-------"<<endl;
    os << "a0     = "<<m_a0<<endl;
    os << "b0     = "<<m_b0<<endl;
    os << "nu0    = "<<m_gw0.get_dof()<<endl;
    os << "m0     = "<<m_m0.t();
    os << "S0     = "<<m_S0.Row(1);
    for(int i=2;i<=m_S0.Ncols();i++)
      os << "         "<<m_S0.Row(i);
    os << "N0     = "<<m_N0.Row(1);
    for(int i=2;i<=m_N0.Ncols();i++)
      os << "         "<<m_N0.Row(i);
    os << "-------hyper-parameters-------"<<endl;
    os << "k      = "<<m_k<<endl;
    os << "alpha  = "<<m_alpha<<endl;
    os << "mu0    = "<<m_gw0.get_mu().t();
    os << "Nu0    = "<<m_gw0.get_Nu().Row(1);
    for(int i=2;i<=m_N0.Ncols();i++)
      os << "         "<<m_gw0.get_Nu().Row(i);
    os << "kappa0 = "<<m_gw0.get_kappa()<<endl;
    os << "-------class-parameters-------"<<endl;
    for(int i=0;i<m_k;i++){
      //os << "cluster "<<i<<endl;
      //os << "n\t=\t"<<m_classnd[i]<<endl;
      os <<m_classnd[i]<<" ";
      //os << m_gw[i];
      //os << endl;
    }
    os << endl;
  }  
  void save(){
    string logsamples   = opts.logfile.value() + ".samples";
    string logmeans     = opts.logfile.value() + ".means";
    string logvariances = opts.logfile.value() + ".variances";
    string zzz = opts.logfile.value() + ".z";
    
    ofstream of_s(logsamples.c_str());
    ofstream of_m(logmeans.c_str());
    ofstream of_v(logvariances.c_str());
    ofstream of_z(zzz.c_str());
    
    double evidence=0;
    double maxlog=0;
    
    of_s << "k\talpha\tlik\n";
    for (unsigned int i=0;i<m_sample_likelihood.size();i++){
      //OUT(i);
      of_s << m_sample_k[i]          << "\t"
	   << m_sample_alpha[i]      << "\t"
	   << m_sample_likelihood[i] << "\n";
      if(m_sample_likelihood[i]>maxlog)
	maxlog=m_sample_likelihood[i];
    }
    // compute evidence
    for(unsigned int i=0;i<m_sample_likelihood.size();i++){
      evidence += std::exp(m_sample_likelihood[i]-maxlog);
    }
    
    // store means and variances
    for(int k=0;k<m_k;k++){
      of_m << m_gw[k].get_smu().t();
      of_v << m_gw[k].get_ssigma().t();
    }
    
    evidence = -log((float)m_sample_likelihood.size()) + maxlog + log(evidence);
    cout<<m_k<<" ";
    cout<<evidence<<endl;
    
    ColumnVector mlz(m_n);
    mlz = get_mldataindex();
    of_z << mlz;
    
    cout<<"final k="<< mlz.MaximumAbsoluteValue()<<endl;
    
  }
  void record(const int samp){
    //cout<<"record sample "<<samp<<endl;
    m_sample_likelihood[samp] = m_likelihood;
    m_sample_k[samp]          = m_k;
    m_sample_alpha[samp]      = m_alpha;
    for(int i=0;i<m_n;i++)
      m_mean_z[i] += m_z[i];
  }
  
  
  
};


#endif
