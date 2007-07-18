#include "dpm_gibbs.h"

bool compare(const pair<float,int> &r1,const pair<float,int> &r2){
  return (r1.first<r2.first);
}

void randomise(vector< pair<float,int> >& r){
  for(unsigned int i=0;i<r.size();i++){
    pair<float,int> p(rand()/float(RAND_MAX),i);
    r[i]=p;
  }
  sort(r.begin(),r.end(),compare);
  
}
std::ostream& operator << (ostream& o,DPM_GibbsSampler& g){
  g.print(o);
  return o;
}
std::ostream& operator << (ostream& o,GaussianWishart& g){
  g.print(o);
  return o;
}


void DPM_GibbsSampler::init(){
  // fix fixed parameters
  m_a0       = 1.0;
  m_b0       = 1.0E8;
  m_S0       << 1000.0*Identity(m_d);//cov(m_data);
  m_N0       << Identity(m_d);///(m_nu0-m_d-1);
  m_m0       = mean(m_data,1).t();
  m_n0       = 1; 

  // initialise all other parameters
  m_alpha    = 1.0;
  m_k        = opts.numclass.value();

  // class hyper parameters
  float kappa0   = 1.0;
  int nu0        = m_d;
  SymmetricMatrix Nu0(m_d);
  Nu0 << m_n0*m_N0;//.01*m_d*Identity(m_d);//cov(m_data);//*(m_nu0-m_d-1);
  ColumnVector mu0(m_d); 
  mu0 = m_m0;
  
  m_gw0      = GaussianWishart(mu0,Nu0,nu0,kappa0);
  
  // class parameters
  if(opts.numclass.value() < 0){ // infinite mixture case
    if(opts.init_class.value() == "oneperdata"){
      if(opts.verbose.value())
	cout << "Initialise with one class per data"<<endl;
      init_oneperdata();
    }
    else if (opts.init_class.value() == "one"){
      if(opts.verbose.value())
	cout << "initialise with one big class"<<endl;
      init_onebigclass();
    }
    else if (opts.init_class.value() == "kmeans"){
      if(opts.verbose.value())
	cout << "Initialise using kmeans" << endl;
      init_kmeans();
    }
    else{ // random
      cout << "Random initialisation using 10 classes" << endl;
      init_random();
    }
  }
  else{ // finite mixture case
    init_kmeans(opts.numclass.value());
  }

  // calculate part of the marginalisation over class mean/variance
  // this part doesn't change through the iterations
  m_margintbase = m_d/2*(log(m_gw0.get_kappa()/(1+m_gw0.get_kappa()))-log(M_PI)) 
    + lgam(float(nu0+1)/2.0) -lgam(float(nu0+1-m_d)/2.0);
    
  // randomised index for loop over data items
  randindex.resize(m_n);
  
  //cout << *this;
}
// different initialisation schemes
void DPM_GibbsSampler::init_oneperdata(){
  m_k = m_n;
  // set parameters
  for(int i=1;i<=m_n;i++){
    GaussianWishart gw(m_d);
    gw.postupdate(m_data.SubMatrix(i,i,1,m_d),m_gw0);
    m_gw.push_back(gw);
    m_z.push_back(i-1);
    m_classnd.push_back(1);
  }
}
void DPM_GibbsSampler::init_onebigclass(){
  m_k = 1;
  GaussianWishart gw(m_d);
  gw.postupdate(m_data,m_gw0);
  m_gw.push_back(gw);
  for(int i=0;i<m_data.Nrows();i++)m_z.push_back(0);
  m_classnd.push_back(m_data.Nrows());
}
void DPM_GibbsSampler::init_kmeans(const int k){
  m_k=k;
  m_z.resize(m_n);
  do_kmeans();
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
}
void DPM_GibbsSampler::init_random(const int k){
  m_k=k;
  m_z.resize(m_n);
  vector< pair<float,int> > rindex(m_n);
  randomise(rindex);
  vector<pair<float,int> >::iterator riter;
  int nn=0,cl=1,nperclass=(int)(float(m_n)/float(m_k));
  for(riter=rindex.begin();riter!=rindex.end();++riter){
    m_z[(*riter).second]=cl;
    nn++;
    if(nn>=nperclass && cl<m_k){
      nn=0;
      cl++;
    }
  }
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
}
  

void DPM_GibbsSampler::sample_parameters(){
  cout << *this;

  // sample indicators
  //cout<<"sample z"<<endl;
  sample_z();
  // sample mean and variance of each class
  //cout<<"sample gw"<<endl;
  sample_gw();
}
void DPM_GibbsSampler::sample_hyperparameters(){
  // sample hyperpriors
  //cout<<"sample gw0"<<endl;
  sample_gw0();
  // sample alpha
  //cout<<"sample alpha"<<endl;
  sample_alpha();
}
// sample indicator variables
void DPM_GibbsSampler::sample_z(){
  ColumnVector datapoint(m_d);
  randomise(randindex);

  // if finite gaussian mixture, do not add new classes
  float extra_finite1 = opts.numclass.value() < 0 ? 0.0 : m_alpha/float(m_k);
  float extra_finite2 = opts.numclass.value() < 0 ? 1.0 : 0.0;
  
  vector< pair<float,int> >::iterator iter;
  for(iter=randindex.begin(); iter!=randindex.end(); ++iter){
    ColumnVector cumsum(m_k+1);
    ColumnVector w(m_k+1);
    
    datapoint=m_data.Row((*iter).second+1).t();
    int oldz=m_z[(*iter).second],newz=oldz;
    m_classnd[oldz] -= 1;
    
    // compute class weights
    double sum=0.0;
    for(int k=0;k<m_k;k++){
      w(k+1) = exp(log(m_classnd[k]+extra_finite1)+marglik(datapoint,k));
      sum += exp(log(m_classnd[k]+extra_finite1)+marglik(datapoint,k));
      cumsum(k+1) = sum;
    }
    w(m_k+1) = m_alpha*exp(margint(datapoint)) * extra_finite2;
    sum += m_alpha*exp(margint(datapoint)) * extra_finite2;
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
    m_z[(*iter).second] = newz;

    if( newz >= m_k ){ // add a new class
      m_k++;
      m_classnd.push_back(1);
      GaussianWishart gw(m_d);
      gw.postupdate(datapoint.t(),m_gw0);
      m_gw.push_back(gw);
    }
    else{
      m_classnd[newz] += 1;
    }
    //cout << " chosen cluster: "<<(*iter).second<<",oldz="<<oldz<<",newz="<<newz;
    //cout << ",w="<<w(newz+1)<<",nold="<<m_classnd[oldz]<<"n="<<m_classnd[newz]<<endl;

  }// end loop over data points
  //cout<<"end data"<<endl<<endl;

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

void DPM_GibbsSampler::sample_gw(){
  // update classes posteriors
  vector< vector<ColumnVector> > data;
  data.resize(m_k);    
  
  // calculate likelihood
  m_likelihood = 0;
  for(int i=0;i<m_n;i++){
    data[ m_z[i] ].push_back(m_data.Row(i+1).t());
    m_likelihood += -marglik(m_data.Row(i+1).t(),m_z[i]);
  }

  for(int k=0;k<m_k;k++){
    if(data[k].size()>0)
      m_gw[k].postupdate(data[k],m_gw0);
  }
}


void DPM_GibbsSampler::sample_gw0(){
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
void DPM_GibbsSampler::sample_alpha(){
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
  m_alpha=(prop*G1.Next()+(1-prop)*G2.Next())/bn;
  //m_alpha=.00000001;
}
  
double DPM_GibbsSampler::marglik(const ColumnVector& data,const int k){
  double res=0.0;
  LogAndSign ld=(2*M_PI*m_gw[k].get_ssigma()).LogDeterminant();

  res -= 0.5*(ld.LogValue()
	      +((data-m_gw[k].get_smu()).t()
	      *m_gw[k].get_ssigma().i()
		*(data-m_gw[k].get_smu())).AsScalar());

  return res;
}
double DPM_GibbsSampler::margint(const ColumnVector& data){
  LogAndSign ld;
  double res=m_margintbase;
  
  ld = m_gw0.get_Nu().LogDeterminant();
  res += ld.LogValue()*m_gw0.get_dof()/2;
  
  SymmetricMatrix A(m_d);
  A << m_gw0.get_Nu()+m_gw0.get_kappa()/(1+m_gw0.get_kappa())*(data-m_gw0.get_mu())*(data-m_gw0.get_mu()).t();
  ld = A.LogDeterminant();
  res -= ld.LogValue()*(m_gw0.get_dof()+1)/2;

  return res;
}


// utils
void DPM_GibbsSampler::do_kmeans(){
  int numiter = 100;
  
  Matrix means(m_d,m_k),newmeans(m_d,m_k);
  ColumnVector nmeans(m_k);
  
  means=0;
  nmeans=0;
  
  //    cout<<"inside kmeans"<<endl;
  // initialise random
  vector< pair<float,int> > rindex(m_n);
  randomise(rindex);
  vector<pair<float,int> >::iterator riter;
  int nn=0,cl=1,nperclass=(int)(float(m_n)/float(m_k));
  for(riter=rindex.begin();riter!=rindex.end();++riter){
    means.Column(cl) += m_data.Row((*riter).second+1).t();
    nmeans(cl) += 1;
    m_z[(*riter).second]=cl;
    nn++;
    if(nn>=nperclass && cl<m_k){
      nn=0;
      cl++;
    }
  }
  for(int m=1;m<=m_k;m++)
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
    for(int i=1;i<=m_n;i++){
      float mindist=1E20,dist=0;
      int mm=1;
      for(int m=1;m<=m_k;m++){
	dist = (means.Column(m)-m_data.Row(i).t()).SumSquare();
	if( dist<mindist){
	  mindist=dist;
	  mm = m;
	}
      }
      m_z[i] = mm;
      newmeans.Column(mm) += m_data.Row(i).t();
      nmeans(mm) += 1;
    }
    
    // compute means
    for(int m=1;m<=m_k;m++){
      if(nmeans(m)==0){
	if(opts.numclass.value()<0) m_k--;
	do_kmeans();
	return;
      }
      newmeans.Column(m) /= nmeans(m);
    }
    means = newmeans;
  }
  
  
  //cout<<"kmeans end"<<endl;
  //for(int i=0;i<n;i++)
  //cout<<z[i]<<" ";
  //cout<<endl;
}    

ReturnMatrix DPM_GibbsSampler::get_dataindex(){
    ColumnVector index(m_n);
    for(unsigned int i=0;i<m_z.size();i++)
      index(i+1) = m_z[i];
    index.Release();
    return index;
}
ReturnMatrix DPM_GibbsSampler::get_mldataindex(){
  ColumnVector index(m_n);
  double lik,tmplik;
  for(int i=1;i<=m_n;i++){
    lik=0.0;tmplik=0;index(i) = 0;
    for(int k=0;k<m_k;k++){
      tmplik = m_classnd[k]*marglik(m_data.Row(i).t(),k);
      if(tmplik>lik && m_classnd[k]>3){
	lik = tmplik;
	index(i) = k+1;
      }
    }
  }
  index.Release();
  return index;
}
