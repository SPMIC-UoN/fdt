/*  Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */

#ifndef __FIBRE_H_
#define __FIBRE_H_


#include <iostream>
#include "stdlib.h"
#include "libprob.h"
#include <cmath>
#include "miscmaths/miscprob.h"

using namespace std; 
using namespace NEWMAT;
using namespace MISCMATHS;

namespace FIBRE{
  
  class Fibre{
    float m_th;
    float m_ph;
    float m_f;
    float m_lam;
    float m_th_prop;
    float m_ph_prop;
    float m_f_prop;
    float m_lam_prop;
    float m_th_old;
    float m_ph_old;
    float m_f_old;
    float m_lam_old;
    float m_th_prior;
    float m_ph_prior;
    float m_f_prior;
    float m_lam_prior;
    float m_th_old_prior;
    float m_ph_old_prior;
    float m_f_old_prior;
    float m_lam_old_prior;
    float m_prior_en;
    float m_old_prior_en;
    int m_th_acc; 
    int m_th_rej;
    int m_ph_acc;
    int m_ph_rej; 
    int m_f_acc;
    int m_f_rej;
    int m_lam_acc; 
    int m_lam_rej;
    ColumnVector m_Signal; 
    ColumnVector m_Signal_old; 
    const float& m_d;
    const ColumnVector& m_alpha;
    const ColumnVector& m_beta;
    const Matrix& m_bvals;
  public:
    //constructors::
    Fibre(const float& d, const ColumnVector& alpha, 
	  const ColumnVector& beta, const Matrix& bvals):
      m_d(d), m_alpha(alpha), m_beta(beta), m_bvals(bvals){

      m_th=0;
      m_th_old=m_th;
      m_ph=0;
      m_ph_old=m_ph;
      m_f=0;
      m_f_old=m_f;
      m_lam=1;
      m_lam_old=m_lam;
      m_th_prop=0.2;
      m_ph_prop=0.2;
      m_f_prop=0.2;
      m_lam_prop=1;

      m_th_prior=0;
      compute_th_prior();

      m_ph_prior=0;
      compute_ph_prior();
      
      m_f_prior=0;
      compute_f_prior();
      
      m_lam_prior=0;
      compute_lam_prior();
      
      m_th_acc=0; m_th_rej=0;
      m_ph_acc=0; m_ph_rej=0;
      m_f_acc=0; m_f_rej=0;
      m_lam_acc=0; m_lam_rej=0;
      m_Signal.ReSize(alpha.Nrows());
      m_Signal=0;
      m_Signal_old=m_Signal;
    }
    Fibre(const float& d, const ColumnVector& alpha, 
	  const ColumnVector& beta, const Matrix& bvals, 
	  const float& th, const float& ph, const float& f, 
	  const float& lam, const int Npts) : 
      m_th(th), m_ph(ph), m_f(f), m_lam(lam), m_d(d), 
      m_alpha(alpha), m_beta(beta), m_bvals(bvals)
     {

      m_th_old=m_th;
      m_ph_old=m_ph;
      m_f_old=m_f;
      m_lam_old=m_lam;
      m_th_prop=0.2;
      m_ph_prop=0.2;
      m_f_prop=0.2;
      m_lam_prop=1;

      m_th_prior=0;
      compute_th_prior();

      m_ph_prior=0;
      compute_ph_prior();
      
      m_f_prior=0;
      compute_f_prior();
      
      m_lam_prior=0;
      compute_lam_prior();
      
      m_th_acc=0; m_th_rej=0;
      m_ph_acc=0; m_ph_rej=0;
      m_f_acc=0; m_f_rej=0;
      m_lam_acc=0; m_lam_rej=0;
      m_Signal.ReSize(alpha.Nrows());
      m_Signal=0;
      m_Signal_old=m_Signal;
    }
    ~Fibre(){}
    
    inline float get_th() const{ return m_th;}
    inline void set_th(const float th){ m_th=th; }
    
    inline float get_ph() const{ return m_ph;}
    inline void set_ph(const float ph){ m_ph=ph; }
    
    inline float get_f() const{ return m_f;}
    inline void set_f(const float f){ m_f=f; }
    
    inline ColumnVector getSignal() const{  //flitney had "inline ColumnVector&" here
      return m_Signal;                      // What is the difference?
    }
    
    inline void restoreSignal() {
      m_Signal=m_Signal_old;
    }
    inline void setSignal(const ColumnVector& Signal){
      m_Signal=Signal;
    }
    
    inline void setSignal(const int i, const float val){
      m_Signal(i)=val;
    }

    inline float get_prior() const{ return m_prior_en;}
    

    inline void update_proposals(){
      m_th_prop*=sqrt(float(m_th_acc+1)/float(m_th_rej+1));
      m_ph_prop*=sqrt(float(m_ph_acc+1)/float(m_ph_rej+1));
      m_f_prop*=sqrt(float(m_f_acc+1)/float(m_f_rej+1));
      m_lam_prop*=sqrt(float(m_lam_acc+1)/float(m_lam_rej+1));
      m_th_acc=0; 
      m_th_rej=0;
      m_ph_acc=0; 
      m_ph_rej=0;
      m_f_acc=0; 
      m_f_rej=0;
      m_lam_acc=0; 
      m_lam_rej=0;
    }
    
    
    inline bool compute_th_prior(){
      m_th_old_prior=m_th_prior;
      m_th_prior=-log(fabs(sin(m_th)/2));
      return false; //instant rejection flag
    }
    inline bool compute_ph_prior(){
      m_ph_old_prior=m_ph_prior;
      m_ph_prior=0;
      return false;
    }
    inline bool compute_f_prior(){
      //note(gamma(lam+1)/(gamma(1)*gamma(lam)) = lam
      // the following is a beta distribution with alpha=0
      m_f_old_prior=m_f_prior;
      if (m_f<0 | m_f>=1 )
	return true;
      else{
	m_f_prior=-(log(m_lam) + (m_lam-1)*log(1-m_f));
	return false;
      }
    }
    
    inline bool compute_lam_prior(){
      m_lam_old_prior=m_lam_prior;
      if(m_lam <0 | m_lam > 1e16)
	return true;
      else{
	m_lam_prior=0;
	return false;
      }
    }
    
    inline void compute_prior(){
      m_old_prior_en=m_prior_en;
      m_prior_en=m_th_prior+m_ph_prior+m_f_prior+m_lam_prior;
    }

     void compute_signal(){
       m_Signal_old=m_Signal;
       for (int i = 1; i <= m_alpha.Nrows(); i++){
 	float angtmp=cos(m_ph-m_beta(i))*sin(m_alpha(i))*sin(m_th) + cos(m_alpha(i))*cos(m_th);
 	angtmp=angtmp*angtmp;
 	m_Signal(i)=exp(-m_d*m_bvals(1,i)*angtmp);
       }
     }


    inline bool propose_th(){
      m_th_old=m_th;
      m_th+=normrnd().AsScalar()*m_th_prop;
      bool rejflag=compute_th_prior();//inside this it stores the old prior
      compute_prior();
      compute_signal();
      return rejflag;
    };
    
    inline void accept_th(){
      m_th_acc++;      
    }
    
    inline void reject_th(){
      m_th=m_th_old;
      m_th_prior=m_th_old_prior;
      m_prior_en=m_old_prior_en;
      m_Signal=m_Signal_old;//Is there a better way of doing this??
      m_th_rej++;
    }
    
    inline bool propose_ph(){
      m_ph_old=m_ph;
      m_ph+=normrnd().AsScalar()*m_ph_prop;
       bool rejflag=compute_ph_prior();//inside this it stores the old prior
      compute_prior();
      compute_signal();
      return rejflag;
    };
    
    inline void accept_ph(){
      m_ph_acc++;
    }
    
    inline void reject_ph(){
      m_ph=m_ph_old;
      m_ph_prior=m_ph_old_prior;
      m_prior_en=m_old_prior_en;
      m_Signal=m_Signal_old;//Is there a better way of doing this??
      m_ph_rej++;
    }
    
    inline bool propose_f(){
      m_f_old=m_f;
      m_f+=normrnd().AsScalar()*m_f_prop;
      bool rejflag=compute_f_prior();
      compute_prior();
      return rejflag;
    };
    
    inline void accept_f(){
      m_f_acc++;
    }
    
    inline void reject_f(){
      m_f=m_f_old;
      m_f_prior=m_f_old_prior;
      m_prior_en=m_old_prior_en;
      m_f_rej++;
    }
    
    inline bool propose_lam(){
      m_lam_old=m_lam;
      m_lam+=normrnd().AsScalar()*m_lam_prop;
      bool rejflag=compute_lam_prior();
      compute_f_prior();
      compute_prior();
      return rejflag;
    };
    
    inline void accept_lam(){
      m_lam_acc++;
    }
    
    inline void reject_lam(){
      m_lam=m_lam_old;
      m_lam_prior=m_lam_old_prior;
      m_prior_en=m_old_prior_en;
      m_lam_rej++;
    }
    
    
    

    friend  ostream& operator<<(ostream& ostr,const Fibre& p);

  };
//overload <<
  inline ostream& operator<<(ostream& ostr,const Fibre& p){
    ostr<<p.m_th<<" "<<p.m_ph<<" "<<p.m_f<<endl;
    return ostr;
  }


  class Multifibre{
    vector<Fibre> m_fibres;
    float m_d;
    float m_d_old;
    float m_d_prop;
    float m_d_prior; 
    float m_d_old_prior;
    float m_d_acc;
    float m_d_rej;
    float m_S0;
    float m_S0_old;
    float m_S0_prop;
    float m_S0_prior;
    float m_S0_old_prior;
    float m_S0_acc;
    float m_S0_rej;
    float m_prior_en;
    float m_old_prior_en;
    float m_likelihood_en;
    float m_old_likelihood_en;
    float m_energy;
    float m_old_energy;
    
    const ColumnVector& m_data;
    const ColumnVector& m_alpha;
    const ColumnVector& m_beta;
    const Matrix& m_bvals;
    
  public:
    Multifibre(const ColumnVector& data,const ColumnVector& alpha, 
	       const ColumnVector& beta, const Matrix& b, const Matrix& Amat, int N ):
    m_data(data), m_alpha(alpha), m_beta(beta), m_bvals(b){
      
      //      initialise(Amat,N);
    }
    
    ~Multifibre(){}
    
    /*
     void initialise(const Matrix& Amat, int N){
      //initialising 
      ColumnVector logS(m_data.Nrows()),tmp(m_data.Nrows()),Dvec(7),dir(3);
      SymmetricMatrix tens;   //Basser's Diffusion Tensor;
      DiagonalMatrix Dd;   //eigenvalues
      Matrix Vd;   //eigenvectors
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
      if(  Dvec(7) >  -23  ){ //23=maxlogfloat
	S0=exp(-Dvec(7));
      }
      else{
	S0=m_data.MaximumAbsoluteValue();
      }

      for ( int i = 1; i <= logS.Nrows(); i++)
	{
	  if(S0<m_data.Sum()/m_data.Nrows()){ S0=m_data.MaximumAbsoluteValue();  }
	  logS(i)=(m_data(i)/S0)>0.01 ? log((i)):log(0.01*S0);
	}
      Dvec = -pinv(Amat)*logS;
      S0=exp(-Dvec(7));
      if(S0<m_data.Sum()/S.Nrows()){ S0=m_data.Sum()/m_data.Nrows();  }
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
      
      float numer=1.5(*(Dd(1)-mDd)*(Dd(1)-mDd)+(Dd(2)-mDd)*(Dd(2)-mDd)+(Dd(3)-mDd)*(Dd(3)-mDd));
      float denom=(Dd(1)*Dd(1)+Dd(2)*Dd(2)+Dd(3)*Dd(3));
      if(denom>0) fsquared=numer/denom;
      else fsquared=0;
      if(fsquared>0){f=sqrt(fsquared);}
      else{f=0;}
      if(f>=0.95) f=0.95;
      if(f<=0.001) f=0.001;
      
      
      m_d=D; m_d_old=m_d;
      m_S0=S0; m_S0_old=m_S0;
      
      
      
	       
      
    }
    
    */

    inline bool compute_d_prior(){
      m_d_old_prior=m_d_prior;
      if(m_d<0)
	return true;
      else{
	m_d_prior=0;
	return false;
      }
    }
    
    inline bool compute_S0_prior(){
      m_S0_old_prior=m_S0_prior;
      if(m_S0<0) return true;
      else
	{    
	  m_S0_prior=0;
	  return false;
	}
    }

    inline void compute_prior(){
      m_old_prior_en=m_prior_en;
      m_prior_en=m_d_prior+m_S0_prior;
      for(unsigned int f=0;f<m_fibres.size(); f++){
	m_prior_en=m_prior_en+m_fibres[f].get_prior();
      } 
    }
    
    void compute_likelihood(){
      m_old_likelihood_en=m_likelihood_en;
      ColumnVector pred(m_alpha.Nrows());
      pred=0;
      float fsum=0;
      for(unsigned int f=0;f<m_fibres.size();f++){
	pred=pred+m_fibres[f].get_f()*m_fibres[f].getSignal();
	fsum+=m_fibres[f].get_f();
      }
      for(int i=1;i<=pred.Nrows();i++){
	pred(i)=pred(i)+(1-fsum)*exp(-m_d*m_bvals(1,i));
      }
      pred=pred*m_S0;
      float sumsquares=(m_data-pred).SumSquare();
      m_likelihood_en=(m_data.Nrows()/2)*log(sumsquares/2);
    }
    
  
    inline  void  compute_energy(){
      m_old_energy=m_energy;
      m_energy=m_prior_en+m_likelihood_en;
    }
    
    bool test_energy(){
      float tmp=exp(m_old_energy-m_energy);
      return (tmp>unifrnd().AsScalar());
    }
    
    inline void restore_energy(){
      m_energy=m_old_energy;
    }
    

    inline bool propose_d(){
      m_d_old=m_d;
      m_d+=normrnd().AsScalar()*m_d_prop;
      bool rejflag=compute_d_prior();//inside this it stores the old prior
      for(unsigned int f=0;f<m_fibres.size();f++)
	m_fibres[f].compute_signal();
      return rejflag;
    };
    
    inline void accept_d(){
      m_d_acc++;      
    }
    
    inline void reject_d(){
      m_d=m_d_old;
      m_d_prior=m_d_old_prior;
      m_prior_en=m_old_prior_en;
      for(unsigned int f=0;f<m_fibres.size();f++)
	m_fibres[f].restoreSignal();
      m_d_rej++;
    }
    
    inline bool propose_S0(){
      m_S0_old=m_S0;
      m_S0+=normrnd().AsScalar()*m_S0_prop;
      bool rejflag=compute_S0_prior();//inside this it stores the old prior
      return rejflag;
    };
    
    inline void accept_S0(){
      m_S0_acc++;      
    }
    
    inline void reject_S0(){
      m_S0=m_S0_old;
      m_S0_prior=m_S0_old_prior;
      m_prior_en=m_old_prior_en;
      m_S0_rej++;
    }
    
    inline void update_proposals(){
      m_d_prop*=sqrt(float(m_d_acc+1)/float(m_d_rej+1));
      m_S0_prop*=sqrt(float(m_S0_acc+1)/float(m_S0_rej+1));
      m_d_acc=0; 
      m_d_rej=0;
      m_S0_acc=0; 
      m_S0_rej=0;
    }

    
    
    

    void jump(){
      
      if(!propose_d()){
	compute_prior();
	compute_likelihood();
	compute_energy();
	if(test_energy())
	  accept_d();
	else
	  reject_d();
      }
      else{ 
	reject_d();
      }
      
      if(!propose_S0()){
	compute_prior();
	compute_likelihood();
	compute_energy();
	if(test_energy())
	  accept_S0();
	else
	  reject_S0();
      }
      else{
	reject_S0();
      }
      
      for(unsigned int f=0;f<m_fibres.size();f++){
	if(!m_fibres[f].propose_th()){
	  compute_prior();
	  compute_likelihood();
	  compute_energy();
	  if(test_energy())
	    m_fibres[f].accept_th();
	  else
	    m_fibres[f].reject_th();
	}
	else {
	  m_fibres[f].reject_th();
	}
	
	if(!m_fibres[f].propose_ph()){
	compute_prior();
	compute_likelihood();
	compute_energy();
	if(test_energy())
	  m_fibres[f].accept_ph();
	else
	  m_fibres[f].reject_ph();
	}
	else{
	  m_fibres[f].reject_ph();
	}
	
	if(!m_fibres[f].propose_f()){
	  compute_prior();
	  compute_likelihood();
	  compute_energy();
	  if(test_energy())
	    m_fibres[f].accept_f();
	  else
	    m_fibres[f].reject_f();
	}
	else{
	  m_fibres[f].reject_f();
	}
	

	if(!m_fibres[f].propose_lam()){
	compute_prior();
	//	compute_likelihood(); //lam does not affect the likelihood.
	compute_energy();
	if(test_energy())
	  m_fibres[f].accept_lam();
	else
	  m_fibres[f].reject_lam();
	}
	else{
	  m_fibres[f].reject_lam();
	}
      }


    }
    


};
  
  
}

#endif
