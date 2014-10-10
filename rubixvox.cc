/*  rubixvox.cc: Classes utilized in RubiX    */
/*  Stamatios Sotiropoulos, FMRIB Analysis Group */
/*  Copyright (C) 2012 University of Oxford  */
/*  CCOPYRIGHT  */


#include "rubixvox.h"

namespace RUBIX{

//Returns the natural log of the 0th order modified Bessel function of first kind for an argument x
//Follows the exponential implementation of the Bessel function in Numerical Recipes, Ch. 6
float logIo(const float x){
  float y,b;

  b=std::fabs(x);
  if (b<3.75){
    float a=x/3.75;
    a*=a;
    //Bessel function evaluation
    y=1.0+a*(3.5156229+a*(3.0899424+a*(1.2067492+a*(0.2659732+a*(0.0360768+a*0.0045813)))));
    y=std::log(y);
  }
  else{
    float a=3.75/b; 
    //Bessel function evaluation
    //y=(exp(b)/sqrt(b))*(0.39894228+a*(0.01328592+a*(0.00225319+a*(-0.00157565+a*(0.00916281+a*(-0.02057706+a*(0.02635537+a*(-0.01647633+a*0.00392377))))))));
    //Logarithm of Bessel function
    y=b+std::log((0.39894228+a*(0.01328592+a*(0.00225319+a*(-0.00157565+a*(0.00916281+a*(-0.02057706+a*(0.02635537+a*(-0.01647633+a*0.00392377))))))))/std::sqrt(b));
  }
  return y;
}



//////////////////////////////////////////////////////////////////////
//       RFibre: Models one anisotropic compartment in an HR voxel
//////////////////////////////////////////////////////////////////////

//Adapt the standard deviation of the proposal distributions during MCMC execution 
//to avoid over-rejection/over-acceptance of samples
void RFibre::update_proposals(){  
  m_th_prop*=sqrt(float(m_th_acc+1)/float(m_th_rej+1));
  m_th_prop=min(m_th_prop,maxfloat);
  m_ph_prop*=sqrt(float(m_ph_acc+1)/float(m_ph_rej+1));
  m_ph_prop=min(m_ph_prop,maxfloat);
  m_f_prop*=sqrt(float(m_f_acc+1)/float(m_f_rej+1));
  m_f_prop=min(m_f_prop,maxfloat);
  m_th_acc=0; 
  m_th_rej=0;
  m_ph_acc=0; 
  m_ph_rej=0;
  m_f_acc=0; 
  m_f_rej=0;
}


//Call that after initializing the parameters
void RFibre::initialise_energies(){
  compute_th_ph_prior();
  compute_f_prior();
  compute_prior();
  compute_signal();
}  


//conditional prior for th and ph, use hyperparameters
bool RFibre::compute_th_ph_prior(){
  m_th_ph_old_prior=m_th_ph_prior;
  m_th_ph_prior=0;
  if (m_Orient_hyp_prior.Nrows()!=0){// && m_f>0.05){ //If an informative orientation prior is used   
    ColumnVector prior_vec(m_Orient_hyp_prior.Nrows());   
    for (int m=1; m<=m_Orient_hyp_prior.Nrows(); m++){ //for each mode of the prior
      float ct=m_Orient_hyp_prior(m,5);  
      float dot=std::min(fabs(m_Orient_hyp_prior(m,1)*m_vec(1)+m_Orient_hyp_prior(m,2)*m_vec(2)+m_Orient_hyp_prior(m,3)*m_vec(3)),1.0);
      prior_vec(m)=(1.0/m_Orient_hyp_prior(m,4))*dot*dot-ct;                       //This is the argument of the Watson (with the logged normalization constant incorporated) 
      //m_th_ph_prior=m_th_ph_prior+exp((1.0/m_Orient_hyp_prior(m,4))*dot*dot-ct); //compute the respective Watson pdf, this is numerically unstable
    }
    if (m_Orient_hyp_prior.Nrows()==1)  //For one mode the prior energy is -log(exp(prior_vec(1))
      m_th_ph_prior=-prior_vec(1);                 
    else{                               //For many modes, instead of calculating -log(Sum(exp(xi))) which is numerically unstable
      float maxarg=prior_vec.Maximum(); //calculate -m-log(Sum(exp(xi-m))), with m=max(xi)
      for (int m=1; m<=m_Orient_hyp_prior.Nrows(); m++) //For each mode of the prior
	m_th_ph_prior=m_th_ph_prior+exp(prior_vec(m)-maxarg);
      m_th_ph_prior=-maxarg-log(m_th_ph_prior); //Convert to -log Energy
    }
    //m_th_ph_prior=m_th_ph_prior-log(0.5*fabs(sin(m_th)));  //Correct with the Jacobian of the transformation! We jump spherical coordinates instead of Cartesian!!
    m_th_ph_prior=m_th_ph_prior-log(fabs(sin(m_th)));        //Correct with the Jacobian of the transformation! We jump spherical coordinates instead of Cartesian!!
    return false; //set instant rejection flag to false
  }
  else{  //if no orientation prior imposed, use uniform prior
    if (m_th==0) 
      m_th_ph_prior=0;
    else
      m_th_ph_prior=-log(0.5*fabs(sin(m_th)));
    return false; //instant rejection flag
  }
}


//ARD prior on f (if f_ard==true)
bool RFibre::compute_f_prior(){
  m_f_old_prior=m_f_prior;
  if (m_f<=0 || m_f>=1 )
    return true;
  else{
    if(!f_ard)
      m_f_prior=0;
    else
      m_f_prior=std::log(m_f);
    m_f_prior=m_ardfudge*m_f_prior;
    return false;
  }
}


//Compute Joint Prior
void RFibre::compute_prior(){
  m_old_prior_en=m_prior_en;
  m_prior_en=m_th_ph_prior+m_f_prior;
}


//Used when orientation prior params are jumped
void RFibre::restore_th_ph_prior(){
  m_th_ph_prior=m_th_ph_old_prior; 
  m_prior_en=m_old_prior_en; 
}


//Compute model predicted signal at High and Low Res, only due to the anisotropic compartment 
void RFibre::compute_signal(){
  m_SignalHR_old=m_SignalHR;
  m_SignalLR_old=m_SignalLR;
       
  if(m_modelnum==1 || m_d_std<1e-5){
    for (int i=1; i<=m_bvecsHR.Ncols(); i++){
      float angtmp=m_vec(1)*m_bvecsHR(1,i)+m_vec(2)*m_bvecsHR(2,i)+m_vec(3)*m_bvecsHR(3,i);
      angtmp=angtmp*angtmp;	  
      m_SignalHR(i)=exp(-m_d*m_bvalsHR(1,i)*angtmp);
    }

    for (int i=1; i<=m_bvecsLR.Ncols(); i++){
      float angtmp=m_vec(1)*m_bvecsLR(1,i)+m_vec(2)*m_bvecsLR(2,i)+m_vec(3)*m_bvecsLR(3,i);
      angtmp=angtmp*angtmp;	  
      m_SignalLR(i)=exp(-m_d*m_bvalsLR(1,i)*angtmp);
    }
  }
  else if(m_modelnum==2){
    //float dbeta=m_d/(m_d_std*m_d_std);
    //float dalpha=m_d*dbeta;     
    float sig2=m_d_std*m_d_std;
    float dalpha=m_d*m_d/sig2;                        
               
    for (int i=1; i<=m_bvecsHR.Ncols(); i++){
      float angtmp=m_vec(1)*m_bvecsHR(1,i)+m_vec(2)*m_bvecsHR(2,i)+m_vec(3)*m_bvecsHR(3,i);
      angtmp=angtmp*angtmp;	  
      //m_SignalHR(i)=exp(log(dbeta/(dbeta + m_bvalsHR(1,i)*angtmp))*dalpha);
      m_SignalHR(i)=exp(log(m_d/(m_d + m_bvalsHR(1,i)*angtmp*sig2))*dalpha); // more stable

    }

    for (int i=1; i<=m_bvecsLR.Ncols(); i++){
      float angtmp=m_vec(1)*m_bvecsLR(1,i)+m_vec(2)*m_bvecsLR(2,i)+m_vec(3)*m_bvecsLR(3,i);
      angtmp=angtmp*angtmp;	  
      //m_SignalLR(i)=exp(log(dbeta/(dbeta + m_bvalsLR(1,i)*angtmp))*dalpha);
      m_SignalLR(i)=exp(log(m_d/(m_d + m_bvalsLR(1,i)*angtmp*sig2))*dalpha); // more stable
    }
  }
}



bool RFibre::propose_th(){
  m_th_old=m_th;
  m_th+=normrnd().AsScalar()*m_th_prop;
  m_vec_old=m_vec;
  m_vec<<sin(m_th)*cos(m_ph) <<sin(m_th)*sin(m_ph) <<cos(m_th);
  bool rejflag=compute_th_ph_prior();//inside this it stores the old prior
  compute_prior();
  compute_signal();
  return rejflag;
}


void RFibre::reject_th(){
  m_th=m_th_old; 
  m_th_ph_prior=m_th_ph_old_prior; 
  m_vec=m_vec_old;
  m_prior_en=m_old_prior_en; 
  restoreSignals(); 
  m_th_rej++;
}
    

bool RFibre::propose_ph(){
  m_ph_old=m_ph;
  m_ph+=normrnd().AsScalar()*m_ph_prop;
  m_vec_old=m_vec;
  m_vec<<sin(m_th)*cos(m_ph) <<sin(m_th)*sin(m_ph) <<cos(m_th);
  bool rejflag=compute_th_ph_prior();//inside this it stores the old prior
  compute_prior();
  compute_signal();
  return rejflag;
}
  
  
void RFibre::reject_ph(){
  m_ph=m_ph_old; 
  m_th_ph_prior=m_th_ph_old_prior; 
  m_vec=m_vec_old;
  m_prior_en=m_old_prior_en; 
  restoreSignals(); 
  m_ph_rej++;
}
    

bool RFibre::propose_f(){
  m_f_old=m_f;
  m_f+=normrnd().AsScalar()*m_f_prop;
  bool rejflag=compute_f_prior();
  compute_prior();
  return rejflag;
}


void RFibre::reject_f(){
  m_f=m_f_old; 
  m_f_prior=m_f_old_prior;
  m_prior_en=m_old_prior_en; 
  m_f_rej++;
}


void RFibre::report() const {
  OUT(m_th);        OUT(m_th_old);       OUT(m_th_prop);     OUT(m_th_acc);    OUT(m_th_rej);
  OUT(m_th_ph_prior); OUT(m_th_ph_old_prior);
  OUT(m_ph);        OUT(m_ph_old);       OUT(m_ph_prop);     OUT(m_ph_acc);    OUT(m_ph_rej);
  OUT(m_f);         OUT(m_f_old);        OUT(m_f_prop);      OUT(m_f_acc);     OUT(m_f_rej);
  OUT(m_f_prior);   OUT(m_f_old_prior);
  OUT(m_prior_en);  OUT(m_old_prior_en);
}




////////////////////////////////////////////////
//       HRvoxel
////////////////////////////////////////////////

//Initialize energies, signals and standard deviations for the proposal distributions
//Call that after initializing the parameters
void HRvoxel::initialise_energies_props(){
  m_S0_prop=m_S0/10.0;
  m_d_prop=m_d/10.0;
  for (int f=0; f<m_numfibres; f++)
    m_fibres[f].initialise_energies();
  if (m_rician){
    compute_tau_prior();
    m_tau_prop=m_tau/2.0;
  }
  if (m_modelnum==2){
    compute_d_std_prior();
    m_d_std_prop=m_d_std/10.0;
  }
  compute_S0_prior();
  compute_d_prior();
  if (m_fsumPrior_ON)
    compute_fsum_prior();
  compute_prior();
  compute_iso_signal();                  
  compute_signal();
}

 
bool HRvoxel::compute_d_prior(){
  m_d_old_prior=m_d_prior;
  if(m_d<=0)
    return true;
  else{
    //m_d_prior=0;
    
    if (m_dPrior_ON)
      m_d_prior=0.5*(m_d-m_mean_d)*(m_d-m_mean_d)/(m_stdev_d*m_stdev_d)+log(m_stdev_d); //Use a Gaussian centered around the neighbourhood meand
    else{  //Use a Gamma Prior centered around 1e-3
      float alpha=3.0; float beta=2000;
      m_d_prior=(1.0-alpha)*log(m_d)+beta*m_d; //Gamma_prior: pow(m_d,alpha-1.0)*exp(-beta*m_d)
    }
    return false;
  }
}
       

bool HRvoxel::compute_d_std_prior(){
  m_d_std_old_prior=m_d_std_prior;
  if(m_d_std<=0 || m_d_std>0.01)
    return true;
  else{
    //m_d_std_prior=0;
    m_d_std_prior=std::log(m_d_std); //Use ARD
    return false;
  }
}


bool HRvoxel::compute_S0_prior(){
  m_S0_old_prior=m_S0_prior;
  if(m_S0<0) return true;
  else{    
    m_S0_prior=0;
    return false;
  }
}


bool HRvoxel::compute_tau_prior(){
  m_tau_old_prior=m_tau_prior;
  if(m_tau<=0) return true;
  else{    
    m_tau_prior=0;
    return false;
  }
}


//Check if sum of volume fractions is >1
bool HRvoxel::reject_f_sum(){
  float fsum=0;//m_f0;
  for(int f=0; f<m_numfibres; f++){
    fsum+=m_fibres[f].get_f();	
  }
  return fsum>1; //true if sum(f) > 1 and therefore, we should reject f
}


//Compute Joint Prior
void HRvoxel::compute_prior(){ 
  m_old_prior_en=m_prior_en;
  m_prior_en=m_d_prior+m_S0_prior;
  if (m_fsumPrior_ON)
    m_prior_en=m_prior_en+m_fsum_prior;
  if(m_rician)
    m_prior_en=m_prior_en+m_tau_prior;
  if(m_modelnum==2)
    m_prior_en=m_prior_en+m_d_std_prior;
  for(int f=0; f<m_numfibres; f++){
    m_prior_en=m_prior_en+m_fibres[f].get_prior();
  } 
}


//Used to update the conditional orientation prior
//when the prior parameters are jumped
void HRvoxel::update_th_ph_prior(){
  for(int f=0; f<m_numfibres; f++){
    m_fibres[f].compute_th_ph_prior();
    m_fibres[f].compute_prior();
  }
  compute_prior();
}


void HRvoxel::restore_th_ph_prior(){
  for(int f=0; f<m_numfibres; f++)
    m_fibres[f].restore_th_ph_prior(); 
  m_prior_en=m_old_prior_en;
} 


void HRvoxel::update_d_prior(){
  compute_d_prior();
  compute_prior();
}


void HRvoxel::restore_d_prior(){
  m_d_prior=m_d_old_prior;
  m_prior_en=m_old_prior_en;
}


void HRvoxel::compute_fsum_prior(){
  m_fsum_old_prior=m_fsum_prior;
  float fsum=0;
  for(int f=0; f<m_numfibres; f++){
    fsum+=m_fibres[f].get_f();	
  }    //Gaussian centered around LRvox neighbourhood mean
  m_fsum_prior=0.5*(fsum-m_mean_fsum)*(fsum-m_mean_fsum)/(m_stdev_fsum*m_stdev_fsum)+log(m_stdev_fsum);
  //m_fsum_prior=0.5*(fsum-0.6)*(fsum-0.6)/(0.1*0.1);
}


void HRvoxel::update_fsum_prior(){
  compute_fsum_prior();
  compute_prior();
}


void HRvoxel::restore_fsum_prior(){
  m_fsum_prior=m_fsum_old_prior;
  m_prior_en=m_old_prior_en;
}



//Compute the predicted signal from the isotropic compartment only
void HRvoxel::compute_iso_signal(){              
  m_iso_SignalHR_old=m_iso_SignalHR;
  m_iso_SignalLR_old=m_iso_SignalLR;

  if(m_modelnum==1 || m_d_std<1e-5){
    for(int i=1; i<=m_bvecsHR.Ncols(); i++)
      m_iso_SignalHR(i)=exp(-m_d*m_bvalsHR(1,i));
    for(int i=1; i<=m_bvecsLR.Ncols(); i++)
      m_iso_SignalLR(i)=exp(-m_d*m_bvalsLR(1,i));
  }
  else if (m_modelnum==2){
    //float dbeta=m_d/(m_d_std*m_d_std);
    //float dalpha=m_d*dbeta;	  
    float sig2=m_d_std*m_d_std;
    float dalpha=m_d*m_d/sig2;                        
    
    for(int i=1; i<=m_bvecsHR.Ncols(); i++)
      //m_iso_SignalHR(i)=exp(log(dbeta/(dbeta+m_bvalsHR(1,i)))*dalpha);
      m_iso_SignalHR(i)=exp(log(m_d/(m_d+m_bvalsHR(1,i)*sig2))*dalpha); //more stable
    for(int i=1; i<=m_bvecsLR.Ncols(); i++)
      //m_iso_SignalLR(i)=exp(log(dbeta/(dbeta+m_bvalsLR(1,i)))*dalpha); 
      m_iso_SignalLR(i)=exp(log(m_d/(m_d+m_bvalsLR(1,i)*sig2))*dalpha); //more stable
  }
}


//Compute the total predicted signal at High and Low Res measurement points
void HRvoxel::compute_signal(){
  m_SignalHR_old=m_SignalHR;
  m_SignalLR_old=m_SignalLR;
  float fsum=0;//m_f0;
  m_SignalHR=0; m_SignalLR=0;
  for(int f=0; f<m_numfibres; f++){   //Signal from the anisotropic compartments
    m_SignalHR=m_SignalHR+m_fibres[f].get_f()*m_fibres[f].getSignalHR();
    m_SignalLR=m_SignalLR+m_fibres[f].get_f()*m_fibres[f].getSignalLR();
    fsum+=m_fibres[f].get_f();               //Total anisotropic volume fraction
  }

  for(int i=1;i<=m_bvecsHR.Ncols();i++)                                     //Add the signal from the isotropic compartment 
    m_SignalHR(i)=m_S0*(m_SignalHR(i)+(1-fsum)*m_iso_SignalHR(i));//+m_f0); //and multiply by S0 to get the total signal  

  for(int i=1;i<=m_bvecsLR.Ncols();i++)                                     //Add the signal from the isotropic compartment 
    m_SignalLR(i)=m_S0*(m_SignalLR(i)+(1-fsum)*m_iso_SignalLR(i));//+m_f0); //and multiply by S0 to get the total signal 
}
 


bool HRvoxel::propose_d(){
  bool rejflag;
  m_d_old=m_d;
  m_d+=normrnd().AsScalar()*m_d_prop;
  rejflag=compute_d_prior();
  for(int f=0; f<m_numfibres; f++)
    m_fibres[f].compute_signal();
  compute_iso_signal();        
  return rejflag;
}


void HRvoxel::reject_d(){
  m_d=m_d_old; 
  m_d_prior=m_d_old_prior;
  for(int f=0; f<m_numfibres; f++)
    m_fibres[f].restoreSignals();
  m_iso_SignalHR=m_iso_SignalHR_old;  
  m_iso_SignalLR=m_iso_SignalLR_old;
  m_d_rej++;
}


bool HRvoxel::propose_d_std(){
  bool rejflag;
  m_d_std_old=m_d_std;
  m_d_std+=normrnd().AsScalar()*m_d_std_prop;
  rejflag=compute_d_std_prior();
  for(int f=0; f<m_numfibres; f++)
    m_fibres[f].compute_signal();
  compute_iso_signal();        
  return rejflag;
}


void HRvoxel::reject_d_std(){
  m_d_std=m_d_std_old; 
  m_d_std_prior=m_d_std_old_prior;
  for(int f=0; f<m_numfibres; f++)
    m_fibres[f].restoreSignals();
  m_iso_SignalHR=m_iso_SignalHR_old;  
  m_iso_SignalLR=m_iso_SignalLR_old;
  m_d_std_rej++;
}


bool HRvoxel::propose_S0(){
  m_S0_old=m_S0;
  m_S0+=normrnd().AsScalar()*m_S0_prop;
  bool rejflag=compute_S0_prior();//inside this it stores the old prior
  return rejflag;
}

    
void HRvoxel::reject_S0(){
  m_S0=m_S0_old;
  m_S0_prior=m_S0_old_prior;
  m_S0_rej++;
}


bool HRvoxel::propose_tau(){
  m_tau_old=m_tau;
  m_tau+=normrnd().AsScalar()*m_tau_prop;
  bool rejflag=compute_tau_prior();//inside this it stores the old prior
  return rejflag;
}

    
void HRvoxel::reject_tau(){
  m_tau=m_tau_old;
  m_tau_prior=m_tau_old_prior;
  m_tau_rej++;
}


void HRvoxel::restore_prior_totsignal(){
  m_prior_en=m_old_prior_en;
  m_SignalLR=m_SignalLR_old;
  m_SignalHR=m_SignalHR_old;
}

void HRvoxel::restore_prior(){
  m_prior_en=m_old_prior_en;
}


//Adapt standard deviation of proposal distributions during MCMC execution 
//to avoid over-rejection/over-acceptance of MCMC samples
void HRvoxel::update_proposals(){
  m_d_prop*=sqrt(float(m_d_acc+1)/float(m_d_rej+1));
  m_d_prop=min(m_d_prop,maxfloat);
  m_d_acc=0; 
  m_d_rej=0;
  if (!m_noS0jump){
    m_S0_prop*=sqrt(float(m_S0_acc+1)/float(m_S0_rej+1));
    m_S0_prop=min(m_S0_prop,maxfloat);
    m_S0_acc=0; 
    m_S0_rej=0;
  }
  if (m_rician){
    m_tau_prop*=sqrt(float(m_tau_acc+1)/float(m_tau_rej+1));
    m_tau_prop=min(m_tau_prop,maxfloat);
    m_tau_acc=0;  m_tau_rej=0;
  }
  if(m_modelnum==2){
    m_d_std_prop*=sqrt(float(m_d_std_acc+1)/float(m_d_std_rej+1));
    m_d_std_prop=min(m_d_std_prop,maxfloat);
    m_d_std_acc=0;  m_d_std_rej=0;
  }
  for(unsigned int f=0; f<m_fibres.size();f++)
    m_fibres[f].update_proposals();
}



void HRvoxel::report() const{
  OUT(m_d);           OUT(m_d_old);           OUT(m_d_prop);
  OUT(m_d_prior);     OUT(m_d_old_prior);     OUT(m_d_acc);          OUT(m_d_rej);
  OUT(m_d_std);       OUT(m_d_std_old);       OUT(m_d_std_prop);
  OUT(m_d_std_prior); OUT(m_d_std_old_prior); OUT(m_d_std_acc);      OUT(m_d_std_rej);
  OUT(m_S0);          OUT(m_S0_old);          OUT(m_S0_prop);        
  OUT(m_S0_prior);    OUT(m_S0_old_prior);    OUT(m_S0_acc);         OUT(m_S0_rej);
  OUT(m_prior_en);    OUT(m_old_prior_en);    OUT(m_SignalHR.t());   OUT(m_SignalLR.t());
  for (int i=0; i<m_numfibres; i++){
    cout <<"fibre "<<i<<endl;
    m_fibres[i].report();} 
}




//Class that models a single mode of the orientation prior
/////////////////////////////
//   Orient_Prior_Mode     //
/////////////////////////////

//Call that after initialising the parameters
void Orient_Prior_Mode::initialise_energies(){
  compute_th_prior();
  compute_ph_prior(); 
  compute_invkappa_prior(); 
  compute_prior();   
}


//Uniform Prior on theta
bool Orient_Prior_Mode::compute_th_prior(){
  m_th_old_prior=m_th_prior;
  if(m_th==0){ m_th_prior=0; }
  else{
   m_th_prior=-log(fabs(sin(m_th)/2.0));
  }
  return false; //instant rejection flag
}


//Uniform Prior on phi
bool Orient_Prior_Mode::compute_ph_prior(){
  m_ph_old_prior=m_ph_prior;
  m_ph_prior=0;
  return false;
}


//ARD Prior on invkappa
bool Orient_Prior_Mode::compute_invkappa_prior(){
  m_invkappa_old_prior=m_invkappa_prior;
  if(m_invkappa<=0)   
    return true;
  else{
    if (m_kappa_ard)
      m_invkappa_prior=log(m_invkappa);
    else
      m_invkappa_prior=0;
    return false;
  }
}


//Compute Joint Prior energy
void Orient_Prior_Mode::compute_prior(){
  m_old_prior_en=m_prior_en;
  m_prior_en=m_th_prior+m_ph_prior+m_invkappa_prior;
}      


//Filter out extreme values of invkappa to ensure calculations are numerically stable 
void Orient_Prior_Mode::filter_invkappa(){
  const float MINinvk=1e-7;  //1e-6
  if (m_invkappa<MINinvk)
    m_invkappa=MINinvk; 
}

  
  

//Compute the normalization constant of a Watson distribution
//with dispersion index 1/kappa. Use a saddlepoint approximation (see Kume & Wood, 2005)
float Orient_Prior_Mode::compute_Watson_norm() {
  float k,ct,R,t,Bm,K2,K3,K4,T;
  
  k=1.0/m_invkappa;
  R=sqrt(4.0*k*k+9.0-4.0*k);
  t=-0.25*(R+3+2.0*k); 
  Bm=1.0/(R+3+2*k); 
  K2=1.0/(2.0*(k+t)*(k+t))+1.0/(t*t);
  K3=-1.0/((k+t)*(k+t)*(k+t))-2.0/(t*t*t);
  K4=3.0/((k+t)*(k+t)*(k+t)*(k+t))+6.0/(t*t*t*t);
  T=K4/(8.0*K2*K2)-5.0*K3*K3/(24.0*K2*K2*K2);
  //ct=4.0*M_PI*sqrt(Bm/R)*exp(-t+T);
  ct=-t+T+log(4.0*M_PI*sqrt(Bm/R));   //compute the log of the normalization constant, to avoid exponential overflow
  return ct;
}


bool Orient_Prior_Mode::propose_th(){
  m_th_old=m_th;
  m_th+=normrnd().AsScalar()*m_th_prop;
  m_vec_old=m_vec;
  m_vec<<sin(m_th)*cos(m_ph) <<sin(m_th)*sin(m_ph) <<cos(m_th);
  bool rejflag=compute_th_prior();//inside this it stores the old prior
  compute_prior();
  return rejflag;
}


void Orient_Prior_Mode::reject_th(){
  m_th=m_th_old; 
  m_th_prior=m_th_old_prior; 
  m_vec=m_vec_old;
  m_prior_en=m_old_prior_en; 
  m_th_rej++;
}
    

bool Orient_Prior_Mode::propose_ph(){
  m_ph_old=m_ph;
  m_ph+=normrnd().AsScalar()*m_ph_prop;
  m_vec_old=m_vec;
  m_vec<<sin(m_th)*cos(m_ph) <<sin(m_th)*sin(m_ph) <<cos(m_th);
  bool rejflag=compute_ph_prior();//inside this it stores the old prior
  compute_prior();
  return rejflag;
}


void Orient_Prior_Mode::reject_ph(){
  m_ph=m_ph_old; 
  m_ph_prior=m_ph_old_prior; 
  m_vec=m_vec_old;
  m_prior_en=m_old_prior_en; 
  m_ph_rej++;
}
 

bool Orient_Prior_Mode::propose_invkappa(){
  m_invkappa_old=m_invkappa;
  m_invkappa+=normrnd().AsScalar()*m_invkappa_prop;
  filter_invkappa();
  m_Watson_norm_old=m_Watson_norm;
  m_Watson_norm=compute_Watson_norm();
  bool rejflag=compute_invkappa_prior();
  compute_prior();
  return rejflag;
}


void Orient_Prior_Mode::reject_invkappa(){
  m_invkappa=m_invkappa_old; 
  m_invkappa_prior=m_invkappa_old_prior;
  m_Watson_norm=m_Watson_norm_old; 
  m_prior_en=m_old_prior_en; 
  m_invkappa_rej++;
}


void Orient_Prior_Mode::update_proposals(){
  m_th_prop*=sqrt(float(m_th_acc+1)/float(m_th_rej+1));
  m_th_prop=min(m_th_prop,maxfloat);
  m_ph_prop*=sqrt(float(m_ph_acc+1)/float(m_ph_rej+1));
  m_ph_prop=min(m_ph_prop,maxfloat);
  m_invkappa_prop*=sqrt(float(m_invkappa_acc+1)/float(m_invkappa_rej+1));
  m_invkappa_prop=min(m_invkappa_prop,maxfloat);
  m_th_acc=0;  m_th_rej=0;
  m_ph_acc=0;  m_ph_rej=0;
  m_invkappa_acc=0; m_invkappa_rej=0;
}





////////////////////////////////////////////////
//       LRvoxel
////////////////////////////////////////////////

//Call that after initialising all parameters
void LRvoxel::initialise_energies(){
  m_S0LR_prop=m_S0LR/10.0;
  compute_S0LR_prior();

  if (m_dPrior_ON){
    m_mean_d_prop=m_mean_d/10.0;
    m_stdev_d_prop=m_stdev_d/10.0;
    compute_meand_prior();   compute_stdevd_prior();
  }
  if (m_fsumPrior_ON){
    compute_mean_fsum_prior();   
    compute_stdev_fsum_prior();
  }

  if (m_rician){
    m_tauLR_prop=m_tauLR/2.0;
    compute_tauLR_prior();
  }
  for (int m=0; m<m_Nmodes; m++)
    m_PModes[m].initialise_energies();
  for (unsigned int n=0; n<m_dataHR.size(); n++)
    m_HRvoxels[n].initialise_energies_props();
  compute_prior();
  compute_likelihood();
  compute_posterior();
  //  cout<<"Likelihood Energy:"<<m_likelihood_en<<endl;
  //cout<<"Prior Energy:"<<m_prior_en<<endl;
  //cout<<"Posterior Energy:"<<m_posterior_en<<endl<<endl;
}


//Update the matrix that keeps information
//on the orientation prior parameters
void LRvoxel::update_Orient_hyp_prior(){
  for (int m=0; m<m_Nmodes; m++){
    ColumnVector temp;
    temp=m_PModes[m].getVec();
    m_Orient_hyp_prior(m+1,1)=temp(1);     
    m_Orient_hyp_prior(m+1,2)=temp(2);
    m_Orient_hyp_prior(m+1,3)=temp(3);
    m_Orient_hyp_prior(m+1,4)=m_PModes[m].get_invkappa();
    m_Orient_hyp_prior(m+1,5)=m_PModes[m].get_Watson_norm();
  }
}


//Update the matrix that keeps information
//on the orientation prior parameters. Update only
//mode M (0<=M<N_modes)
void LRvoxel::update_Orient_hyp_prior(int M){
    ColumnVector temp;
    temp=m_PModes[M].getVec();
    m_Orient_hyp_prior(M+1,1)=temp(1);     
    m_Orient_hyp_prior(M+1,2)=temp(2);
    m_Orient_hyp_prior(M+1,3)=temp(3);
    m_Orient_hyp_prior(M+1,4)=m_PModes[M].get_invkappa();
    m_Orient_hyp_prior(M+1,5)=m_PModes[M].get_Watson_norm();
}



//Set params for a single HRvoxel with index 0 <= n < N (model 1)
void LRvoxel::set_HRparams(const int n, const float d, const float S0, const ColumnVector& th, const ColumnVector& ph, const ColumnVector& f){
  m_HRvoxels[n].set_d(d);
  m_HRvoxels[n].set_S0(S0);
  if (m_allard && !m_Noard){
    m_HRvoxels[n].addfibre(th(1),ph(1),f(1),true);    //ARD on for the first fibre
  }
  else{
    m_HRvoxels[n].addfibre(th(1),ph(1),f(1),false);   //ARD off for the first fibre
  }
  for (int m=2; m<=m_numfibres; m++){
    if (m_Noard){
      m_HRvoxels[n].addfibre(th(m),ph(m),f(m),false); //ARD off for the other fibres
    }
    else{
      m_HRvoxels[n].addfibre(th(m),ph(m),f(m),true);  //ARD on for the other fibres
    }
  }
}


//Set params for a single HRvoxel with index 0 <= n < N (model 2)
void LRvoxel::set_HRparams(const int n, const float d, const float d_std,const float S0, const ColumnVector& th, const ColumnVector& ph, const ColumnVector& f){
  m_HRvoxels[n].set_d(d);
  m_HRvoxels[n].set_d_std(d_std);
  m_HRvoxels[n].set_S0(S0);
  if (m_allard && !m_Noard)
    m_HRvoxels[n].addfibre(th(1),ph(1),f(1),true);  //ARD on for the first fibre
  else
    m_HRvoxels[n].addfibre(th(1),ph(1),f(1),false); //ARD off for the first fibre
  for (int m=2; m<=m_numfibres; m++){
    if (m_Noard)
      m_HRvoxels[n].addfibre(th(m),ph(m),f(m),false); //ARD off for the other fibres
    else
      m_HRvoxels[n].addfibre(th(m),ph(m),f(m),true);  //ARD on for the other fibres
  }
}



//Set params for all modes of orientation priors
void LRvoxel::set_Priorparams(const ColumnVector& th, const ColumnVector& ph, const ColumnVector& invkappa){
  for (int m=1; m<=m_Nmodes; m++){
    m_PModes[m-1].set_th(th(m));
    m_PModes[m-1].set_ph(ph(m));
    m_PModes[m-1].set_invkappa(invkappa(m));
  }
}



bool LRvoxel::propose_S0LR(){
  m_S0LR_old=m_S0LR;
  m_S0LR+=normrnd().AsScalar()*m_S0LR_prop;
  bool rejflag=compute_S0LR_prior();//inside this it stores the old prior
  return rejflag;
}

    
void LRvoxel::reject_S0LR(){
  m_S0LR=m_S0LR_old;
  m_S0LR_prior=m_S0LR_old_prior;
  m_S0LR_rej++;
}


bool LRvoxel::compute_S0LR_prior(){
  m_S0LR_old_prior=m_S0LR_prior;
  if(m_S0LR<0) return true;
  else{    
    m_S0LR_prior=0;
    return false;
  }
}



bool LRvoxel::propose_meand(){
  m_mean_d_old=m_mean_d;
  m_mean_d+=normrnd().AsScalar()*m_mean_d_prop;
  bool rejflag=compute_meand_prior();//inside this it stores the old prior
  return rejflag;
}

void LRvoxel::reject_meand(){
  m_mean_d=m_mean_d_old;
  m_mean_d_prior=m_mean_d_old_prior;
  m_mean_d_rej++;
}


bool LRvoxel::compute_meand_prior(){
  m_mean_d_old_prior=m_mean_d_prior;
  if(m_mean_d<=0) return true;
  else{    
    //m_mean_d_prior=0;
    float alpha=3.0;
    float beta=2000;  //Gamma_prior centered around 1e-3
    m_mean_d_prior=(1.0-alpha)*log(m_mean_d)+beta*m_mean_d; //Gamma_prior: pow(m_mean_d,alpha-1.0)*exp(-beta*m_mean_d)

    return false;
  }
}


bool LRvoxel::propose_stdevd(){
  m_stdev_d_old=m_stdev_d;
  m_stdev_d+=normrnd().AsScalar()*m_stdev_d_prop;
  bool rejflag=compute_stdevd_prior();//inside this it stores the old prior
  return rejflag;
}


void LRvoxel::reject_stdevd(){
  m_stdev_d=m_stdev_d_old;
  m_stdev_d_prior=m_stdev_d_old_prior;
  m_stdev_d_rej++;
}

bool LRvoxel::compute_stdevd_prior(){
  m_stdev_d_old_prior=m_stdev_d_prior;
  if(m_stdev_d<=0 || m_stdev_d>=0.002) return true;
  else{    
    m_stdev_d_prior=0;
    //m_stdev_d_prior=log(m_stdev_d); //ARD prior
    return false;
  }
}


bool LRvoxel::propose_mean_fsum(){
  m_mean_fsum_old=m_mean_fsum;
  m_mean_fsum+=normrnd().AsScalar()*m_mean_fsum_prop;
  bool rejflag=compute_mean_fsum_prior();//inside this it stores the old prior
  return rejflag;
}


void LRvoxel::reject_mean_fsum(){
  m_mean_fsum=m_mean_fsum_old;
  m_mean_fsum_prior=m_mean_fsum_old_prior;
  m_mean_fsum_rej++;
}


bool LRvoxel::compute_mean_fsum_prior(){
  m_mean_fsum_old_prior=m_mean_fsum_prior;
  if(m_mean_fsum<=0 || m_mean_fsum>=1) return true;
  else{    
    //m_mean_fsum_prior=0;
    m_mean_fsum_prior=-log(m_mean_fsum)-log(1-m_mean_fsum);   //Beta distribution with a=b=2: fsum^(a-1)*(1-fsum)^(b-1)
    return false;
  }
}


bool LRvoxel::propose_stdev_fsum(){
  m_stdev_fsum_old=m_stdev_fsum;
  m_stdev_fsum+=normrnd().AsScalar()*m_stdev_fsum_prop;
  bool rejflag=compute_stdev_fsum_prior();//inside this it stores the old prior
  return rejflag;
}


void LRvoxel::reject_stdev_fsum(){
  m_stdev_fsum=m_stdev_fsum_old;
  m_stdev_fsum_prior=m_stdev_fsum_old_prior;
  m_stdev_fsum_rej++;
}


bool LRvoxel::compute_stdev_fsum_prior(){
  m_stdev_fsum_old_prior=m_stdev_fsum_prior;
  if(m_stdev_fsum<=0 || m_stdev_fsum>=0.1) return true; //Without restriction to (0,0.1), very large values (e.g. 100's) are preferred for stdev, so that the fsum_prior becomes uniform 
  else{    
    m_stdev_fsum_prior=0;      
    //m_stdev_fsum_prior=log(m_stdev_fsum);      //ARD prior  
    //m_stdev_fsum_prior=-log(25*m_stdev_fsum)-log(1-5*m_stdev_fsum);      //Beta distribution with a=2,b=2, defined on [0,0.2].  
    //m_stdev_fsum_prior=-2*log(1-m_stdev_fsum);   //Beta distribution with a=1, b=3. Similar to an ARD, but is naturally restricted to [0,1] 
    return false;
  }
}


bool LRvoxel::propose_tauLR(){
  m_tauLR_old=m_tauLR;
  m_tauLR+=normrnd().AsScalar()*m_tauLR_prop;
  bool rejflag=compute_tauLR_prior();//inside this it stores the old prior
  return rejflag;
}

    
void LRvoxel::reject_tauLR(){
  m_tauLR=m_tauLR_old;
  m_tauLR_prior=m_tauLR_old_prior;
  m_tauLR_rej++;
}


bool LRvoxel::compute_tauLR_prior(){
  m_tauLR_old_prior=m_tauLR_prior;
  if(m_tauLR<=0) return true;
  else{    
    m_tauLR_prior=0;
    return false;
  }
}


void LRvoxel::compute_prior(){
  m_old_prior_en=m_prior_en;  
  m_prior_en=0;

  m_prior_en=m_S0LR_prior;
  if (m_rician)
    m_prior_en+=m_tauLR_prior;
  for (int m=0; m<m_Nmodes; m++)
    m_prior_en+=m_PModes[m].get_prior();
    
  for (unsigned int m=0; m<m_dataHR.size(); m++)
    m_prior_en+=m_HRvoxels[m].get_prior();
  
  if (m_dPrior_ON)
    m_prior_en+=m_mean_d_prior+m_stdev_d_prior;
  if (m_fsumPrior_ON)
    m_prior_en+=m_mean_fsum_prior+m_stdev_fsum_prior;
}


void LRvoxel::compute_likelihood(){
  ColumnVector SLRpred, SHRpred, Diff;
  SLRpred.ReSize(m_bvecsLR.Ncols()); 
  SHRpred.ReSize(m_bvecsHR[0].Ncols()); SHRpred=0;
  float likLR, likHR;

  m_old_likelihood_en=m_likelihood_en;

  if (!m_rician){ //Gaussian Likelihood Energy Calculation
    //likelihood of Low-Res data 
    SLRpred=0;
    for (unsigned int m=0; m<m_dataHR.size(); m++) //add attenuations
      SLRpred+=m_HRweights(m+1)*m_HRvoxels[m].getSignalLR()/m_HRvoxels[m].get_S0();
    SLRpred=m_S0LR*SLRpred;
    SLRpred=m_dataLR-SLRpred;
    likLR=0.5*m_bvecsLR.Ncols()*log(0.5*SLRpred.SumSquare()); 

    //likelihood of High-Res data
    likHR=0; 
    for (unsigned int m=0; m<m_dataHR.size(); m++){
      SHRpred=m_HRvoxels[m].getSignalHR();
      SHRpred=m_dataHR[m]-SHRpred;
      likHR+=log(0.5*SHRpred.SumSquare());
    }
    likHR*=0.5*m_bvecsHR[0].Ncols();      

    m_likelihood_en=likLR+likHR;
  }
  else{  //Rician Likelihood Energy Calculation
    //likelihood of Low-Res data 
    SLRpred=0;
    for (unsigned int m=0; m<m_dataHR.size(); m++) //add attenuations
      SLRpred+=m_HRweights(m+1)*m_HRvoxels[m].getSignalLR()/m_HRvoxels[m].get_S0();
    SLRpred=m_S0LR*SLRpred;
    likLR=-m_bvecsLR.Ncols()*log(m_tauLR);  
    for (int k=1; k<=m_bvecsLR.Ncols(); k++)
      likLR-=m_logdataLR(k)-0.5*m_tauLR*(m_dataLR(k)*m_dataLR(k)+SLRpred(k)*SLRpred(k))+logIo(m_tauLR*m_dataLR(k)*SLRpred(k));

    //likelihood of High-Res data
    likHR=0; 
    for (unsigned int m=0; m<m_dataHR.size(); m++){    
      float m_tauHR=m_HRvoxels[m].get_tau();
      SHRpred=m_HRvoxels[m].getSignalHR();
      likHR+=m_bvecsHR[0].Ncols()*log(m_tauHR);
      for (int k=1; k<=m_bvecsHR[0].Ncols(); k++)
	likHR+=m_logdataHR[m](k)-0.5*m_tauHR*(m_dataHR[m](k)*m_dataHR[m](k)+SHRpred(k)*SHRpred(k))+logIo(m_tauHR*m_dataHR[m](k)*SHRpred(k));
    }
    m_likelihood_en=likLR-likHR;
  }
}


void LRvoxel::compute_posterior(){
  m_old_posterior_en=m_posterior_en;
  m_posterior_en=m_prior_en+m_likelihood_en;
}


//Test whether the candidate state will be accepted 
bool LRvoxel::test_energy() const{
  float tmp=exp(m_old_posterior_en-m_posterior_en);
  return (tmp>unifrnd().AsScalar());
}
 

void LRvoxel::restore_energies(){
  m_prior_en=m_old_prior_en;
  m_likelihood_en=m_old_likelihood_en;
  m_posterior_en=m_old_posterior_en;
}


void LRvoxel::restore_Prior_Posterior(){
  m_prior_en=m_old_prior_en;
  m_posterior_en=m_old_posterior_en;
}


void LRvoxel::update_proposals(){
  for(int m=0; m<m_Nmodes; m++)
    m_PModes[m].update_proposals();
  for(unsigned int n=0; n<m_dataHR.size(); n++)
    m_HRvoxels[n].update_proposals();

  if (!m_noS0jump){
    m_S0LR_prop*=sqrt(float(m_S0LR_acc+1)/float(m_S0LR_rej+1));
    m_S0LR_prop=min(m_S0LR_prop,maxfloat);
    m_S0LR_acc=0; m_S0LR_rej=0;
  }
  if (m_dPrior_ON){
    m_mean_d_prop*=sqrt(float(m_mean_d_acc+1)/float(m_mean_d_rej+1));
    m_mean_d_prop=min(m_mean_d_prop,maxfloat);
    m_mean_d_acc=0; m_mean_d_rej=0;

    m_stdev_d_prop*=sqrt(float(m_stdev_d_acc+1)/float(m_stdev_d_rej+1));
    m_stdev_d_prop=min(m_stdev_d_prop,maxfloat);
    m_stdev_d_acc=0; m_stdev_d_rej=0;
  }

  if (m_fsumPrior_ON){
    m_mean_fsum_prop*=sqrt(float(m_mean_fsum_acc+1)/float(m_mean_fsum_rej+1));
    m_mean_fsum_prop=min(m_mean_fsum_prop,maxfloat);
    m_mean_fsum_acc=0; m_mean_fsum_rej=0;

    m_stdev_fsum_prop*=sqrt(float(m_stdev_fsum_acc+1)/float(m_stdev_fsum_rej+1));
    m_stdev_fsum_prop=min(m_stdev_fsum_prop,maxfloat);
    m_stdev_fsum_acc=0; m_stdev_fsum_rej=0;
  }

  if (m_rician){
    m_tauLR_prop*=sqrt(float(m_tauLR_acc+1)/float(m_tauLR_rej+1));
    m_tauLR_prop=min(m_tauLR_prop,maxfloat);
    m_tauLR_acc=0; m_tauLR_rej=0;
  }
}


//Single MCMC iteration with all parameters jumped
void LRvoxel::jump(){

  //Jump LRvoxel parameters first
 
  if (!m_noS0jump){
    if(!propose_S0LR()){    //Try S0LR 
      compute_prior();
      compute_likelihood();
      compute_posterior();
      if(test_energy()){
	accept_S0LR();
      }
      else{
	restore_energies();
	reject_S0LR();
      }
    }
    else 
      reject_S0LR();
  }

  if (m_dPrior_ON){
    if(!propose_meand()){    //Try mean_d 
      for (unsigned int n=0; n<m_dataHR.size(); n++) //update nuisance_prior for each HR voxel
	m_HRvoxels[n].update_d_prior();
    
      compute_prior();
      compute_posterior();
      if(test_energy()){
	accept_meand();
      }
      else{
	restore_Prior_Posterior();
	for (unsigned int n=0; n<m_dataHR.size(); n++) //restore nuisance_prior_energy for each HR voxel
	  m_HRvoxels[n].restore_d_prior();
	reject_meand();
      }
    }
    else 
      reject_meand();
 
    if(!propose_stdevd()){    //Try stdev_d 
      for (unsigned int n=0; n<m_dataHR.size(); n++) //update nuisance_prior energy for each HR voxel
	m_HRvoxels[n].update_d_prior();
    
      compute_prior();
      compute_posterior();
      if(test_energy()){
	accept_stdevd();
      }
      else{
	restore_Prior_Posterior();
	for (unsigned int n=0; n<m_dataHR.size(); n++) //restore nuisance_prior_energy for each HR voxel
	  m_HRvoxels[n].restore_d_prior();
	reject_stdevd();
      }
    }
    else 
      reject_stdevd();
  }

  if (m_fsumPrior_ON){
    if(!propose_mean_fsum()){    //Try mean_fsum
      for (unsigned int n=0; n<m_dataHR.size(); n++) //update nuisance_prior for each HR voxel
	m_HRvoxels[n].update_fsum_prior();
    
      compute_prior();
      compute_posterior();
      if(test_energy()){
	accept_mean_fsum();
      }
      else{
	restore_Prior_Posterior();
	for (unsigned int n=0; n<m_dataHR.size(); n++) //restore nuisance_prior_energy for each HR voxel
	  m_HRvoxels[n].restore_fsum_prior();
	reject_mean_fsum();
      }
    }
    else 
      reject_mean_fsum();

    if(!propose_stdev_fsum()){    //Try stdev_fsum 
      for (unsigned int n=0; n<m_dataHR.size(); n++) //update nuisance_prior energy for each HR voxel
	m_HRvoxels[n].update_fsum_prior();
    
      compute_prior();
      compute_posterior();
      if(test_energy()){
	accept_stdev_fsum();
      }
      else{
	restore_Prior_Posterior();
	for (unsigned int n=0; n<m_dataHR.size(); n++) //restore nuisance_prior_energy for each HR voxel
	  m_HRvoxels[n].restore_fsum_prior();
	reject_stdev_fsum();
      }
    }
    else 
      reject_stdev_fsum();
  }

  if (m_rician){
    if(!propose_tauLR()){    //Try tauLR 
      compute_prior();
      compute_likelihood();
      compute_posterior();
      if(test_energy()){
	accept_tauLR();
      }
      else{
	restore_energies();
	reject_tauLR();
      }
    }
    else 
      reject_tauLR();
  }

  //For each mode of the Orientation Prior
  for(int m=0; m<m_Nmodes; m++){ 
  
    if(!m_PModes[m].propose_th()){    //Try theta 
      update_Orient_hyp_prior(m);
      for (unsigned int n=0; n<m_dataHR.size(); n++) //update th_ph_prior for each HR voxel
	m_HRvoxels[n].update_th_ph_prior();
	
      compute_prior();
      //compute_likelihood();
      compute_posterior();
      if(test_energy()){
	m_PModes[m].accept_th();
      }
      else{
	//restore_energies();
	restore_Prior_Posterior();
	for (unsigned int n=0; n<m_dataHR.size(); n++) //update th_ph_prior for each HR voxel
	  m_HRvoxels[n].restore_th_ph_prior();
	m_PModes[m].reject_th();
	update_Orient_hyp_prior(m);
      }
    }
    else 
      m_PModes[m].reject_th();
    

    if(!m_PModes[m].propose_ph()){    //Try phi 
      update_Orient_hyp_prior(m);
      for (unsigned int n=0; n<m_dataHR.size(); n++) //update th_ph_prior for each HR voxel
	m_HRvoxels[n].update_th_ph_prior();
	
      compute_prior();
      //compute_likelihood();
      compute_posterior();
      if(test_energy()){
	m_PModes[m].accept_ph();
      }
      else{
	//restore_energies();
	restore_Prior_Posterior();
	for (unsigned int n=0; n<m_dataHR.size(); n++) //update th_ph_prior for each HR voxel
	  m_HRvoxels[n].restore_th_ph_prior();
	m_PModes[m].reject_ph();
	update_Orient_hyp_prior(m);
      }
    }
    else 
      m_PModes[m].reject_ph();


    if(!m_PModes[m].propose_invkappa()){    //Try invkappa 
      update_Orient_hyp_prior(m);
      for (unsigned int n=0; n<m_dataHR.size(); n++) //update th_ph_prior for each HR voxel
	m_HRvoxels[n].update_th_ph_prior();
	
      compute_prior();
      //compute_likelihood();
      compute_posterior();
      if(test_energy()){
	m_PModes[m].accept_invkappa();
      }
      else{
	//restore_energies();
	restore_Prior_Posterior();
	for (unsigned int n=0; n<m_dataHR.size(); n++) //update th_ph_prior for each HR voxel
	  m_HRvoxels[n].restore_th_ph_prior();
	m_PModes[m].reject_invkappa();
	update_Orient_hyp_prior(m);
      }
    }
    else 
      m_PModes[m].reject_invkappa();
  } //end Orientation Prior params


  //For each HRvoxel
  for (unsigned int n=0; n<m_dataHR.size(); n++){

    if(!m_HRvoxels[n].propose_d()){    //Try d 
      m_HRvoxels[n].compute_prior();
      m_HRvoxels[n].compute_signal();
      compute_prior();
      compute_likelihood();
      compute_posterior();
      if(test_energy()){
	m_HRvoxels[n].accept_d();
      }
      else{
	restore_energies();
	m_HRvoxels[n].restore_prior_totsignal();
	m_HRvoxels[n].reject_d();
      }
    }
    else 
      m_HRvoxels[n].reject_d();
      

    if (m_modelnum==2){
      if(!m_HRvoxels[n].propose_d_std()){    //Try d_std
	m_HRvoxels[n].compute_prior();
	m_HRvoxels[n].compute_signal();
	compute_prior();
	compute_likelihood();
	compute_posterior();
	if(test_energy()){
	  m_HRvoxels[n].accept_d_std();
	}
	else{
	  restore_energies();
	  m_HRvoxels[n].restore_prior_totsignal();
	  m_HRvoxels[n].reject_d_std();
	}
      }
      else 
	m_HRvoxels[n].reject_d_std();
    }

    if (!m_noS0jump){
      if(!m_HRvoxels[n].propose_S0()){    //Try S0
	m_HRvoxels[n].compute_prior();
	m_HRvoxels[n].compute_signal();
	compute_prior();
	compute_likelihood();
	compute_posterior();
	if(test_energy())
	  m_HRvoxels[n].accept_S0();
	else{
	  restore_energies();
	  m_HRvoxels[n].restore_prior_totsignal();
	  m_HRvoxels[n].reject_S0();
	}
      }
      else
	m_HRvoxels[n].reject_S0();
    }

    if (m_rician){
      if(!m_HRvoxels[n].propose_tau()){    //Try tau_HR
	m_HRvoxels[n].compute_prior();
	compute_prior();
	compute_likelihood();
	compute_posterior();
	if(test_energy())
	  m_HRvoxels[n].accept_tau();
	else{
	  restore_energies();
	  m_HRvoxels[n].restore_prior();
	  m_HRvoxels[n].reject_tau();
	}
      }
      else
	m_HRvoxels[n].reject_tau();
    }

    for(int f=0; f<m_numfibres; f++){   //For each fibre in the HR voxel
      if(!m_HRvoxels[n].propose_th(f)){ //Try theta
	m_HRvoxels[n].compute_prior();
	m_HRvoxels[n].compute_signal();
      	compute_prior();
	compute_likelihood();
	compute_posterior();
	if(test_energy())
	    m_HRvoxels[n].accept_th(f);
	else{
	  restore_energies();
	  m_HRvoxels[n].restore_prior_totsignal();
	  m_HRvoxels[n].reject_th(f);
	}
      }
      else 
	m_HRvoxels[n].reject_th(f);
	
	
      if(!m_HRvoxels[n].propose_ph(f)){ //Try phi
	m_HRvoxels[n].compute_prior();
	m_HRvoxels[n].compute_signal();
      	compute_prior();
	compute_likelihood();
	compute_posterior();
	if(test_energy())
	    m_HRvoxels[n].accept_ph(f);
	else{
	  restore_energies();
	  m_HRvoxels[n].restore_prior_totsignal();
	  m_HRvoxels[n].reject_ph(f);
	}
      }
      else 
	m_HRvoxels[n].reject_ph(f);
	 

      if(!m_HRvoxels[n].propose_f(f)){  //Try f 
	if(!m_HRvoxels[n].reject_f_sum()){
	  m_HRvoxels[n].compute_fsum_prior();
	  m_HRvoxels[n].compute_prior();
	  m_HRvoxels[n].compute_signal();
	  compute_prior();
	  compute_likelihood();
	  compute_posterior();
	  if(test_energy())
	    m_HRvoxels[n].accept_f(f);
	  else{
	    restore_energies();
	    m_HRvoxels[n].restore_prior_totsignal();
	    m_HRvoxels[n].restore_fsum_prior();
	    m_HRvoxels[n].reject_f(f);
	  }
	}
	else   //else for rejectin fsum>1
	  m_HRvoxels[n].reject_f(f);
      }
      else    //else for rejecting rejflag returned from propose_f()
	m_HRvoxels[n].reject_f(f);
    } //end anisotropic compartment parameters
   
  } //end HR voxel parameters
}

} //end namespace
