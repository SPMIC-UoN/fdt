/*  Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(_GIBBS_H)
#define _GIBBS_H

#include "newmat.h"
#include "newran.h"
#include "miscmaths/miscmaths.h"

using namespace NEWMAT;
using namespace NEWRAN;

class GibbsSampler
{
 private:
  friend std::ostream& operator << (ostream& o,GibbsSampler& g);
 protected:

  int m_numiter;
  int m_burnin;
  int m_sampleevery;
  int m_n;
  int m_d;

  const Matrix& m_data;

 public:
  GibbsSampler(const Matrix& data,int numiter,int burnin,int sampleevery):
    m_numiter(numiter),m_burnin(burnin),m_sampleevery(sampleevery),m_data(data){
    
    m_n=data.Nrows();
    m_d=data.Ncols();
  }
  virtual ~GibbsSampler(){}

  virtual void init() = 0 ;
  virtual void record(const int) = 0;
  virtual void sample_parameters() = 0;
  virtual void sample_hyperparameters() = 0;
  virtual void print(ostream&) = 0;

  void  run(){
    Random::Set(rand() / float(RAND_MAX));

    int recordcount=0;

    // burnin period (no sampling)
    cout<<"burnin"<<endl;
    for(int i=0;i<m_burnin;i++){
      //cout<<"-----------"<<endl;
      sample_parameters();
      sample_hyperparameters();
      //print();
    }

    // m_numiter=2;
    
    // after burnin, sample ervery "sampleevery"
    cout<<"gibbs"<<endl;
    int samp = 0;
    for(int i=m_burnin;i<m_numiter;i++){
      sample_parameters();
      sample_hyperparameters();
      
      //print();
	
      recordcount++;
      
      if(recordcount==m_sampleevery){
	//cout<<"record"<<endl;
	record(samp);samp++;
	recordcount=0;
      }
    }
    
    
  }

};

std::ostream& operator << (ostream& o,GibbsSampler& g){
  g.print(o);
  return o;
}


#endif
