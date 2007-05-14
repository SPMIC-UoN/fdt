#if !defined(_DPM_GIBBS_H)
#define _DPM_GIBBS_H

#include "gibbs.h"
#include "dpmOptions.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include <stdlib.h>
#include <stdio.h>
#include <newmat.h>
#include <newran.h>
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
std::ostream& operator << (ostream& o,GaussianWishart& g){
  g.print(o);
  return o;
}

bool compare(const pair<int,float> &p1,const pair<int,float> &p2){
  return (p1.second < p2.second) ? true : false;
}


class DPM_GibbsSampler : public GibbsSampler
{
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
  float                    m_margintbase;
  
  vector< pair<int,float> > randindex;

  // samples
  vector<float>            m_sample_alpha;
  vector<int>              m_sample_k;
  vector<double>           m_sample_likelihood;
  double                   m_likelihood;
  int                      m_nsamples;

public:
  DPM_GibbsSampler(const Matrix& data,int numiter,int burnin,int sampleevery):
    GibbsSampler(data,numiter,burnin,sampleevery),
    opts(DPM::dpmOptions::getInstance()){
    m_nsamples = (int)floor( (numiter - burnin) / sampleevery );

    m_sample_alpha.resize(m_nsamples);
    m_sample_k.resize(m_nsamples);
    m_sample_likelihood.resize(m_nsamples);
  }
  ~DPM_GibbsSampler(){}

  void init(){
    // fix fixed parameters
    m_a0       = 1.0;
    m_b0       = 1.0;
    m_S0       << 1000.0*Identity(m_d);//cov(m_data);
    m_N0       << Identity(m_d);///(m_nu0-m_d-1);
    m_m0       = mean(m_data,1).t();
    m_n0       = 1; 

    // initialise all other parameters
    m_alpha    = 1.0;
    m_k        = opts.numclass.value();

    float kappa0   = 1.0;
    int nu0        = m_d;
    SymmetricMatrix Nu0(m_d);
    Nu0 << m_n0*m_N0;//.01*m_d*Identity(m_d);//cov(m_data);//*(m_nu0-m_d-1);
    ColumnVector mu0(m_d); 
    mu0 = m_m0;
    
    m_gw0      = GaussianWishart(mu0,Nu0,nu0,kappa0);

    
    if(m_k < 0){
      if(opts.init_class.value() == "oneperdata"){
	cout << "Initialise with one class per data"<<endl;
	m_k = m_n;
	// set parameters
	for(int i=1;i<=m_n;i++){
	  GaussianWishart gw(m_d);
	  gw.postupdate(m_data.SubMatrix(i,i,1,m_d),m_gw0);
	  m_gw.push_back(gw);
	  m_z.push_back(i-1);         // one data point per class
	  m_classnd.push_back(1);
	}
	cout << *this;
      }
      else if (opts.init_class.value() == "one"){
	cout << "initialise with one big class"<<endl;
	m_k = 1;
	// initialise with one big class
	GaussianWishart gw(m_d);
	gw.postupdate(m_data,m_gw0);
	m_gw.push_back(gw);
	for(int i=0;i<m_data.Nrows();i++)m_z.push_back(0);
	m_classnd.push_back(m_data.Nrows());
	cout << *this;
      }
      else{ // kmeans initialisation
	cout << "Initialise using kmeans" << endl;
	m_z.resize(m_n);
	m_k = 10;
	kmeans(m_data,m_z,m_k);
	cout<<"done"<<endl;
	for(int k=1;k<=m_k;k++){
	  GaussianWishart gw(m_d);
	  vector<ColumnVector> dat;
	  for(int i=1;i<=m_n;i++)
	    if(m_z[i-1] == k){
	      dat.push_back(m_data.Row(i).t());
	      m_z[i-1] -- ;
	    }
	  gw.postupdate(dat,m_gw0);
	  m_gw.push_back(gw);
	  m_classnd.push_back((int)dat.size());
	}
	cout << *this;
      }
    }
    else{
      m_z.resize(m_n);
      kmeans(m_data,m_z,m_k);

      for(int k=1;k<=m_k;k++){
	GaussianWishart gw(m_d);
	vector<ColumnVector> dat;
	for(int i=1;i<=m_n;i++)
	  if(m_z[i-1] == k){
	    dat.push_back(m_data.Row(i).t());
	    m_z[i-1] -- ;
	  }
	//OUT(dat.size());
	gw.postupdate(dat,m_gw0);
	m_gw.push_back(gw);
	m_classnd.push_back((int)dat.size());
      }

    }

    m_margintbase = m_d/2*(log(m_gw0.get_kappa()/(1+m_gw0.get_kappa()))-log(M_PI)) 
      + lgam(float(nu0+1)/2.0) -lgam(float(nu0+1-m_d)/2.0);

    OUT(m_margintbase);
    //print();


    // randomised index for loop over data items
    randindex.resize(m_n);

  }
  void sample_parameters(){
    // sample indicators
    //cout<<"sample z"<<endl;
    sample_z();
    // sample mean and variance of each class
    //cout<<"sample gw"<<endl;
    sample_gw();

    cout << *this;

  }
  void sample_hyperparameters(){
    // sample hyperpriors
    //cout<<"sample gw0"<<endl;
    sample_gw0();
    // sample alpha
    //cout<<"sample alpha"<<endl;
    sample_alpha();
  }
  void randomise(vector< pair<int,float> >& r){
    for(int i=0;i<m_n;i++){
      pair<int,float> p(i,rand());
      r[i]=p;
    }
    sort(r.begin(),r.end(),compare);
  }
  // sample indicator variables
  void sample_z(){
    ColumnVector datapoint(m_d);

    randomise(randindex);

    float extra_finite1 = opts.numclass.value() < 0 ? 0.0 : m_alpha/float(m_k);
    float extra_finite2 = opts.numclass.value() < 0 ? 1.0 : 0.0;

    vector< pair<int,float> >::iterator iter;
    for(iter=randindex.begin(); iter!=randindex.end(); iter++){
      ColumnVector cumsum(m_k+1);

      datapoint=m_data.Row((*iter).first+1).t();
      int oldz=m_z[(*iter).first],newz=oldz;
      m_classnd[oldz] -= 1;

      //cout<<"-----"<<endl;
      // compute class weights
      float sum=0.0;
      for(int k=0;k<m_k;k++){
	sum += (m_classnd[k]+extra_finite1)*marglik(datapoint,k);
	cumsum(k+1) = sum;
      }
      sum += m_alpha*margint(datapoint) * extra_finite2;
      cumsum(m_k+1) = sum;
      // sample z using the weights
      float U=rand()/float(RAND_MAX);
      U *= sum;
      for(int k=1;k<=m_k+1;k++){
	if(U<cumsum(k)){
	  newz=k-1;
	  break;
	}
      }
      m_z[(*iter).first] = newz;
      
      //cout<<"-----"<<endl;

      if( newz >= m_k ){ // add a new class
	//cout<<"ADD A NEW CLASS"<<endl;
	m_k++;
	m_classnd.push_back(1);
	GaussianWishart gw(m_d);
	gw.postupdate(datapoint.t(),m_gw0);
	m_gw.push_back(gw);
      }
      else{
	m_classnd[newz] += 1;
      }
    }// end loop over data points

    // delete empty classes if in infinite mode
    if(opts.numclass.value()<0){
      for(int k=m_k-1;k>=0;k--){
	if(m_classnd[k] == 0){
	  for(int i=0;i<m_n;i++)
	    if(m_z[i]>k)m_z[i]--;
	  for(int kk=k;kk<m_k-1;kk++){
	    m_classnd[kk]=m_classnd[kk+1];
	    m_gw[kk]=m_gw[kk+1];
	  }
	  m_classnd.pop_back();
	  m_gw.pop_back();
	  m_k--;
	}
      }
    }

  }

  void sample_gw(){
    // update classes posteriors
    vector< vector<ColumnVector> > data;
    data.resize(m_k);    
    
    m_likelihood = 0;
    for(int i=0;i<m_n;i++){
      data[ m_z[i] ].push_back(m_data.Row(i+1).t());
      m_likelihood += -std::log(marglik(m_data.Row(i+1).t(),m_z[i]));
    }

    for(int k=0;k<m_k;k++){
      if(data[k].size()>0)
	m_gw[k].postupdate(data[k],m_gw0);
    }
    
  }


  void sample_gw0(){
    SymmetricMatrix Nu0(m_d),A(m_d),S(m_d);
    ColumnVector a(m_d),mu0(m_d);    
    float B=0;
    
    A=0;a=0;
    for(int k=0;k<m_k;k++){
      S = m_gw[k].get_ssigma().i();
      a += S*m_gw[k].get_smu();
      A << A+S;
      B += ((m_gw[k].get_smu()-m_gw0.get_mu()).t()*S*(m_gw[k].get_smu()-m_gw0.get_mu())).AsScalar();
    }
    S << A+m_N0.i();
    A << (A+m_S0.i()).i();
    a = A*(a+m_S0.i()*m_m0);
    

    Nu0 = wishrnd(S.i(),(m_k+1)*m_gw0.get_dof());
    mu0 = mvnrnd(a.t(),A).t();

    m_gw0.set_Nu(Nu0);
    m_gw0.set_mu(mu0);
      
    Gamma G(1+m_k*m_d/2); //G.Set(rand()/float(RAND_MAX));
    m_gw0.set_kappa(G.Next()*2/(1+B));
    //m_gw0.set_kappa(1.0);

    
  }

  // sample from alpha using additional variable eta
  void sample_alpha(){
    float eta,prop;
    float ak=m_a0+m_k-1,bn;

    Gamma G1(ak+1);       //G1.Set(rand()/float(RAND_MAX));
    Gamma G2(ak);         //G2.Set(rand()/float(RAND_MAX));
    Gamma B1(m_alpha+1);  //B1.Set(rand()/float(RAND_MAX));
    Gamma B2(m_n);        //B2.Set(rand()/float(RAND_MAX));

    eta  = B1.Next();
    eta /= (eta+B2.Next());
    bn   = m_b0-std::log(eta);

    prop=ak/(ak+m_n*bn);
    m_alpha=(prop*G1.Next()+(1-prop)*G2.Next())*bn;
  }
  
  float marglik(const ColumnVector& data,const int k){
    float res;
    res=normpdf(data,m_gw[k].get_smu(),m_gw[k].get_ssigma());

    //OUT(res);

    return res;
  }
  float margint(const ColumnVector& data){
    LogAndSign ld;
    float res=m_margintbase;

    ld = m_gw0.get_Nu().LogDeterminant();
    res += ld.LogValue()*m_gw0.get_dof()/2;

    SymmetricMatrix A(m_d);
    A << m_gw0.get_Nu()+m_gw0.get_kappa()/(1+m_gw0.get_kappa())*(data-m_gw0.get_mu())*(data-m_gw0.get_mu()).t();
    ld = A.LogDeterminant();
    res -= ld.LogValue()*(m_gw0.get_dof()+1)/2;

    //OUT(exp(res));

    return std::exp(res);
  }
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
      os << "cluster "<<i<<endl;
      os << "n\t=\t"<<m_classnd[i]<<endl;
      os << m_gw[i];
      os << endl;
    }
  }

  
  void save(){
    string logsamples   = opts.logfile.value() + ".samples";
    string logmeans     = opts.logfile.value() + ".means";
    string logvariances = opts.logfile.value() + ".variances";

    ofstream of_s(logsamples.c_str());
    ofstream of_m(logmeans.c_str());
    ofstream of_v(logvariances.c_str());
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
    cout<<evidence<<endl;;
    
  }
  void record(const int samp){
    //cout<<"record sample "<<samp<<endl;
    m_sample_likelihood[samp] = m_likelihood;
    m_sample_k[samp]          = m_k;
    m_sample_alpha[samp]      = m_alpha;
  }

  void kmeans(const Matrix& data,vector<int>& z,const int k){
    int numiter = 100;
    int n = data.Nrows();
    int d = data.Ncols();

    Matrix means(d,k),newmeans(d,k);
    ColumnVector nmeans(k);
    //z.resize(n);
    
    means=0;
    nmeans=0;

    //    cout<<"inside kmeans"<<endl;
    // initialise random
    vector< pair<int,float> > rindex(n);
    randomise(rindex);
    vector<pair<int,float> >::iterator riter;
    int nn=0,cl=1,nperclass=(int)(float(n)/float(k));
    for(riter=rindex.begin();riter!=rindex.end();riter++){
      means.Column(cl) += data.Row((*riter).first+1).t();
      nmeans(cl) += 1;
      z[(*riter).first]=cl;
      nn++;
      if(nn>=nperclass && cl<k){
	nn=0;
	cl++;
      }
    }
    for(int m=1;m<=k;m++)
      means.Column(m) /= nmeans(m);

    //cout<<"kmeans init"<<endl;
    //for(int i=0;i<n;i++)
    //cout<<z[i]<<" ";
    //cout<<endl;

    // iterate
    for(int iter=0;iter<numiter;iter++){
      // loop over datapoints and attribute z for closest mean
      newmeans=0;
      nmeans=0;
      for(int i=1;i<=n;i++){
	float mindist=1E20,dist=0;
	int mm=1;
	for(int m=1;m<=k;m++){
	  dist = (means.Column(m)-data.Row(i).t()).SumSquare();
	  if( dist<mindist){
	    mindist=dist;
	    mm = m;
	  }
	}
	z[i] = mm;
	newmeans.Column(mm) += data.Row(i).t();
	nmeans(mm) += 1;
      }

      // compute means
      for(int m=1;m<=k;m++){
	if(nmeans(m)==0){
	  kmeans(data,z,k);
	  return;
	}
	newmeans.Column(m) /= nmeans(m);
      }
      means = newmeans;
    }

    //OUT(n);
    //OUT(d);
    //OUT(nmeans.t());
    //OUT(newmeans);

    //cout<<"kmeans end"<<endl;
    //for(int i=0;i<n;i++)
    //cout<<z[i]<<" ";
    //cout<<endl;
  }    

};


#endif
