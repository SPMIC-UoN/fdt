/*  Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(_GIBBS_H)
#define _GIBBS_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>

using namespace std;

class GibbsSampler
{
 protected:

  int m_numiter;
  int m_burnin;
  int m_sampleevery;
  int m_n;
  int m_d;

 public:
  GibbsSampler(int numiter,int burnin,int sampleevery):
    m_numiter(numiter),m_burnin(burnin),m_sampleevery(sampleevery){}
  virtual ~GibbsSampler(){}

  virtual void init() = 0 ;
  virtual void record(const int) = 0;
  virtual void sample_parameters() = 0;
  virtual void sample_hyperparameters() = 0;

  void  run(){

    int recordcount=0;

    // burnin period (no sampling)
    cout<<"burnin"<<endl;
    for(int i=0;i<m_burnin;i++){
      //cout<<"-----------"<<endl;
      sample_parameters();
      sample_hyperparameters();
      
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



#endif
