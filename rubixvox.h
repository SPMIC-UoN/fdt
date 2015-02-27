/*  rubixvox.h: Classes utilized in RubiX    */

/*  Stam Sotiropoulos, FMRIB Analysis Group */

/*  Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(rubixvox_h)
#define rubixvox_h

#include <iostream>
#include <fstream>
#include <iomanip>
#define WANT_STREAM
#define WANT_MATH
#include <string>
#include <math.h>
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "stdlib.h"
#include "utils/log.h"

using namespace std; 
using namespace NEWMAT;
using namespace MISCMATHS;

#define UPPERDIFF 0.005

namespace RUBIX{
  const float maxfloat=1e10;
  const float minfloat=1e-10;
  const float maxlogfloat=23;
  const float minlogfloat=-23;

 //////////////////////////////////////////////////////////////////////
 //       RFibre: Models one anisotropic compartment in an HR voxel
 //////////////////////////////////////////////////////////////////////
 class RFibre{
    float m_th;                   //Current/candidate MCMC state
    float m_ph;
    float m_f;
    ColumnVector m_vec;           //Holds the current/candidate orientation in cartesian coordinates
    ColumnVector m_vec_old;
    float m_th_prop;              //Standard deviation for Gaussian proposal distributions of parameters
    float m_ph_prop;
    float m_f_prop;
    float m_th_old;               //Last accepted value. If a sample is rejected this value is restored
    float m_ph_old;
    float m_f_old;
    float m_th_ph_prior;          //Priors for the model parameters 
    float m_f_prior;
    float m_th_ph_old_prior;
    float m_f_old_prior;
    float m_prior_en;             //Joint Prior 
    float m_old_prior_en;
    int m_th_acc;     
    int m_th_rej;
    int m_ph_acc;
    int m_ph_rej; 
    int m_f_acc;
    int m_f_rej;
    bool f_ard;                     //By default ARD is on, on the volume fraction f
    ColumnVector m_SignalHR;        //Vector that stores the predicted signal from the anisotropic compartment during the candidate/current MCMC state at High-Res measurement points
    ColumnVector m_SignalHR_old;
    ColumnVector m_SignalLR;        //Vector that stores the predicted signal from the anisotropic compartment during the candidate/current MCMC state at Low-Res measurement points 
    ColumnVector m_SignalLR_old;

    const Matrix& m_Orient_hyp_prior;//Matrix Nmodes x 5 that contains the hyperparameters for the orientation prior 
                                     //columns 1-3 contains the (x,y,z) coordinates for the mode, 4th column contains the invkappa value, 5th the Watson normalization constant
    const float m_ardfudge;
    const float& m_d;
    const float& m_d_std;
    const float& m_R;               //R=Perpendicular/Axial diffusivity in model3 
    const Matrix& m_bvecsHR;        //bvecs at High-Res   (3 x HR_NumPoints)
    const Matrix& m_bvalsHR;        //bvalues at High-Res (1 x HR_NumPoints)
    const Matrix& m_bvecsLR;        //bvecs at Low-Res    (3 x LR_NumPoints)
    const Matrix& m_bvalsLR;        //bvalues at Low-Res  (1 x HR_NumPoints)
    const int m_modelnum;           //1 for deconvolution with sticks, 2 for deconvolution with sticks and a Gamma distribution of diffusivities, 3 for deconvolution with zeppelins

 public:
 //constructor
 RFibre(const float th, const float ph,const float f, const Matrix& bvecsHR, const Matrix& bvalsHR, 
	const Matrix& bvecsLR, const Matrix& bvalsLR, const Matrix& Orient_hyp_prior,
	const float& d, const float& d_std, const float& R, const bool ard=true, const int modelnum=1, const float ardfudge=1):
    m_th(th), m_ph(ph), m_f(f), f_ard(ard), m_Orient_hyp_prior(Orient_hyp_prior),  m_ardfudge(ardfudge), m_d(d),
      m_d_std(d_std), m_R(R), m_bvecsHR(bvecsHR), m_bvalsHR(bvalsHR), m_bvecsLR(bvecsLR), m_bvalsLR(bvalsLR), m_modelnum(modelnum){
      m_th_old=m_th; m_ph_old=m_ph;
      m_vec.ReSize(3);
      m_vec<<sin(m_th)*cos(m_ph) <<sin(m_th)*sin(m_ph) <<cos(m_th);
      m_vec_old=m_vec;
      m_f_old=m_f;
      m_th_prop=0.2; m_ph_prop=0.2; m_f_prop=m_f/5.0;
	
      m_th_ph_prior=0; m_th_ph_old_prior=0;
      m_f_prior=0; m_f_old_prior=0;
	
      m_SignalHR.ReSize(m_bvecsHR.Ncols());  m_SignalHR=0;  m_SignalHR_old=m_SignalHR;
      m_SignalLR.ReSize(m_bvecsLR.Ncols());  m_SignalLR=0;  m_SignalLR_old=m_SignalLR;

      m_th_acc=0; m_th_rej=0;
      m_ph_acc=0; m_ph_rej=0;
      m_f_acc=0; m_f_rej=0;
    }
    
    ~RFibre(){}
    
    inline float get_th() const{ return m_th;}
    inline float get_ph() const{ return m_ph;}
    inline const ColumnVector& getVec() const { return m_vec; }
    inline float get_f() const{ return m_f;}
    //inline void set_th(const float th){ m_th=th;   m_vec<<sin(m_th)*cos(m_ph) <<sin(m_th)*sin(m_ph) <<cos(m_th); }
    //inline void set_ph(const float ph){ m_ph=ph;   m_vec<<sin(m_th)*cos(m_ph) <<sin(m_th)*sin(m_ph) <<cos(m_th); }
    //inline void set_f(const float f){ m_f=f; }

    inline const ColumnVector& getSignalHR() const{ return m_SignalHR; }
    inline const ColumnVector& getSignalLR() const{ return m_SignalLR; }
    inline void restoreSignals() { m_SignalHR=m_SignalHR_old; m_SignalLR=m_SignalLR_old; }
    inline float get_prior() const{ return m_prior_en;}
    
    void initialise_energies(); //Initialize energies, signals 
    void update_proposals(); 
    bool compute_th_ph_prior();
    bool compute_f_prior();
    void compute_prior();
    void compute_signal(); //Compute model predicted signal, only due to the anisotropic compartment 
    void restore_th_ph_prior();

    bool propose_th();
    inline void accept_th(){ m_th_acc++; }   
    void reject_th();
    
    bool propose_ph();
    inline void accept_ph(){ m_ph_acc++; }
    void reject_ph();
    
    bool propose_f();
    inline void accept_f(){ m_f_acc++; }
    void reject_f();

    void report() const;

    RFibre& operator=(const RFibre& rhs){
      m_th=rhs.m_th; m_ph=rhs.m_ph; m_f=rhs.m_f;   
      m_vec=rhs.m_vec;  m_vec_old=rhs.m_vec_old;   
      m_th_prop=rhs.m_th_prop; m_ph_prop=rhs.m_ph_prop; m_f_prop=rhs.m_f_prop;
      m_th_old=rhs.m_th_old; m_ph_old=rhs.m_ph_old; m_f_old=rhs.m_f_old;
      m_th_ph_prior=rhs.m_th_ph_prior; m_f_prior=rhs.m_f_prior;
      m_th_ph_old_prior=rhs.m_th_ph_old_prior; m_f_old_prior=rhs.m_f_old_prior;
      m_prior_en=rhs.m_prior_en; m_old_prior_en=rhs.m_old_prior_en;
      m_th_acc=rhs.m_th_acc; m_th_rej=rhs.m_th_rej; 
      m_ph_acc=rhs.m_ph_acc; m_ph_rej=rhs.m_ph_rej;
      m_f_acc=rhs.m_f_acc;   m_f_rej=rhs.m_f_rej;
      f_ard=rhs.f_ard;
      m_SignalHR=rhs.m_SignalHR; m_SignalHR_old=rhs.m_SignalHR_old;
      m_SignalLR=rhs.m_SignalLR; m_SignalLR_old=rhs.m_SignalLR_old;
      return *this;
    } 
    
  };



 //////////////////////////////
 //       HRvoxel            //
 //////////////////////////////
  class HRvoxel{
    vector<RFibre> m_fibres;
    float m_d;
    float m_d_old;
    float m_d_prop;
    float m_d_prior; 
    float m_d_old_prior;
    float m_d_acc;
    float m_d_rej;
 
    float m_d_std;
    float m_d_std_old;
    float m_d_std_prop;
    float m_d_std_prior; 
    float m_d_std_old_prior;
    float m_d_std_acc;
    float m_d_std_rej;

    float m_R;   //anisotropy parameter in model3 - R=l2/l1 for the zeppelin compartment 
    float m_R_old;
    float m_R_prop;
    float m_R_prior; 
    float m_R_old_prior;
    float m_R_acc;
    float m_R_rej;

    float m_S0;
    float m_S0_old;
    float m_S0_prop;
    float m_S0_prior;
    float m_S0_old_prior;
    float m_S0_acc;
    float m_S0_rej;
    
    float m_tau;                    //Noise at a High-res voxel (precision)
    float m_tau_old;
    float m_tau_prop;
    float m_tau_prior;
    float m_tau_old_prior;
    float m_tau_acc;
    float m_tau_rej;

    const float& m_mean_d;        //hyperparameters for d prior
    const float& m_stdev_d;
    const float& m_mean_fsum;     //hyperparameters for fsum prior
    const float& m_stdev_fsum;
    const float& m_R_priormean;   //Parameters to use for the Gaussian prior on R. Mean
    const float& m_R_priorstd;    //and variance
    const float& m_R_priorfudge;  //and fudge factor

    float m_fsum_prior;
    float m_fsum_old_prior;

    float m_prior_en;               //Joint Prior
    float m_old_prior_en;
    ColumnVector m_iso_SignalHR;    //Vector that stores the predicted signal from the isotropic compartment only
    ColumnVector m_iso_SignalHR_old; 
    ColumnVector m_iso_SignalLR;    //Vector that stores the predicted signal from the isotropic compartment only 
    ColumnVector m_iso_SignalLR_old; 
 
    ColumnVector m_SignalHR;        //Vector that stores the total predicted signal from the specific voxel at High-Res measurement points
    ColumnVector m_SignalHR_old;
    ColumnVector m_SignalLR;        //Vector that stores the total predicted signal from the specific voxel at Low-Res measurement points
    ColumnVector m_SignalLR_old;

    const Matrix& m_Orient_hyp_prior;//Matrix Nmodes x 5 that contains the hyperparameters for the orientation prior 
    
    const Matrix& m_bvecsHR;        //bvecs at High-Res   (3 x HR_NumPoints)
    const Matrix& m_bvalsHR;        //bvalues at High-Res (1 x HR_NumPoints)
    const Matrix& m_bvecsLR;        //bvecs at Low-Res    (3 x LR_NumPoints)
    const Matrix& m_bvalsLR;        //bvalues at Low-Res  (1 x HR_NumPoints)
    const int m_modelnum;           //1 for single-shell, 2 for multi-shell model
    const int m_numfibres;          //number of fibres in this HR voxel
    const float m_ardfudge;
    const bool m_fsumPrior_ON;      //Flag that indicates whether a prior on the fsum will be imposed (based on the LRvox neighbourhood mean)
    const bool m_dPrior_ON;         //Flag that indicates whether a prior on diffusivity will be imposed (based on the LRvox neighbourhood mean)
    const bool m_rician;            //Indicates whether Rician Noise model is used
    const bool m_noS0jump;          //Indicates whether S0 parameters will be kept constant during MCMC 
    //const ColumnVector& LRvox_inter; //array with indices on LR voxels intersected by this HR
 
  public:
    //Constructor
    HRvoxel(const Matrix& bvecsHR, const Matrix& bHR, 
	    const Matrix& bvecsLR, const Matrix& bLR, 
	    const Matrix& Orient_hyp_prior,
	    const float& mean_d, const float& stdev_d, const float& mean_fsum, const float& stdev_fsum, const float& Rmean, const float& Rstd, const float& Rfudge,
	    const int N, const int modelnum=1,const float ardfudge=1, const bool fsumPrior_ON=false, const bool dPrior_ON=false, const bool rician=false, const bool noS0jump=false):
    m_mean_d(mean_d), m_stdev_d(stdev_d), m_mean_fsum(mean_fsum), m_stdev_fsum(stdev_fsum), m_R_priormean(Rmean), m_R_priorstd(Rstd), m_R_priorfudge(Rfudge),
    m_Orient_hyp_prior(Orient_hyp_prior), 
      m_bvecsHR(bvecsHR), m_bvalsHR(bHR), m_bvecsLR(bvecsLR), m_bvalsLR(bLR), m_modelnum(modelnum), m_numfibres(N),  m_ardfudge(ardfudge), m_fsumPrior_ON(fsumPrior_ON), m_dPrior_ON(dPrior_ON), m_rician(rician), m_noS0jump(noS0jump) {
      
      //Initialize vectors that keep the signal from the isotropic compartment
      m_iso_SignalHR.ReSize(m_bvecsHR.Ncols()); m_iso_SignalHR=0; m_iso_SignalHR_old=m_iso_SignalHR; 
      m_SignalHR.ReSize(m_bvecsHR.Ncols()); m_SignalHR=0; m_SignalHR_old=m_SignalHR; 
      m_iso_SignalLR.ReSize(m_bvecsLR.Ncols()); m_iso_SignalLR=0; m_iso_SignalLR_old=m_iso_SignalLR;            
      m_SignalLR.ReSize(m_bvecsLR.Ncols()); m_SignalLR=0; m_SignalLR_old=m_SignalLR; 
     
      m_d_acc=0; m_d_rej=0; m_d_prior=0; m_d_old_prior=0;              
      m_d_std_acc=0; m_d_std_rej=0; m_d_std_prior=0; m_d_std_old_prior=0;  
      m_R_acc=0; m_R_rej=0; m_R_prior=0; m_R_old_prior=0;              
      m_S0_acc=0; m_S0_rej=0; m_S0_prior=0; m_S0_old_prior=0;
      m_d=0; m_d_old=0;  m_d_std=0; m_d_std_old=0; m_R=0; m_R_old=0; m_S0=0; m_S0_old=0; 
      m_tau_acc=0; m_tau_rej=0; m_tau_prior=0; m_tau_old_prior=0; m_tau=0; m_tau_old=0; 
      m_prior_en=0; m_old_prior_en=0;
      m_fsum_prior=0; m_fsum_old_prior=0;
    }

    //Destructor
    ~HRvoxel(){}


    const vector<RFibre>& fibres() const{ return m_fibres; }
    void addfibre(const float th, const float ph, const float f,const bool use_ard=true) {
      RFibre fib(th,ph,f,m_bvecsHR,m_bvalsHR,m_bvecsLR,m_bvalsLR,m_Orient_hyp_prior,m_d,m_d_std,m_R,use_ard,m_modelnum,m_ardfudge);
      m_fibres.push_back(fib);
    }
	
    void initialise_energies_props(); //Initialize energies, signals and standard deviations for the proposal distributions
    inline float get_d() const{ return m_d;}
    inline void set_d(const float d){ m_d=d; }
    inline float get_d_std() const{ return m_d_std; }
    inline void set_d_std(const float d_std){ m_d_std=d_std; }
    inline float get_R() const{ return m_R; }
    inline void set_R(const float R){ m_R=R; }
    inline float get_S0() const{ return m_S0; }
    inline void set_S0(const float S0){ m_S0=S0; }
    inline float get_tau() const{ return m_tau; }
    inline void set_tau(const float tau){ m_tau=tau; }
    inline float get_prior() const { return m_prior_en; }

    inline const ColumnVector& getSignalHR() const{ return m_SignalHR; }
    inline const ColumnVector& getSignalLR() const{ return m_SignalLR; }

    bool compute_d_prior();
    bool compute_d_std_prior(); 
    bool compute_R_prior(); 
    bool compute_S0_prior(); 
    bool compute_tau_prior(); 
    bool reject_f_sum();       //Check if sum of volume fractions is >1

    void update_d_prior();
    void restore_d_prior();

    void compute_fsum_prior();
    void update_fsum_prior();
    void restore_fsum_prior();

    void compute_prior();      //Compute Joint Prior energy
    void compute_iso_signal(); //Compute the predicted signal from the isotropic compartment only
    void compute_signal();     //Compute the total predicted signal

    bool propose_d();
    void accept_d(){ m_d_acc++; }
    void reject_d();

    bool propose_d_std();
    void accept_d_std(){ m_d_std_acc++; }
    void reject_d_std();

    bool propose_R();
    void accept_R(){ m_R_acc++; }
    void reject_R();

    bool propose_S0();
    void accept_S0(){ m_S0_acc++; }
    void reject_S0();

    bool propose_tau();
    void accept_tau(){ m_tau_acc++; }
    void reject_tau();

    bool propose_th(int n){ return m_fibres[n].propose_th(); }  //n is the fibre index from 0 to N-1
    void accept_th(int n) { m_fibres[n].accept_th(); }
    void reject_th(int n) { m_fibres[n].reject_th(); }
    bool propose_ph(int n){ return m_fibres[n].propose_ph(); }  
    void accept_ph(int n) { m_fibres[n].accept_ph(); }
    void reject_ph(int n) { m_fibres[n].reject_ph(); }
    bool propose_f(int n) { return m_fibres[n].propose_f(); }
    void accept_f(int n)  { m_fibres[n].accept_f(); }
    void reject_f(int n)  { m_fibres[n].reject_f(); }

    void restore_prior_totsignal();
    void restore_prior();
    
    void update_th_ph_prior();  //Used to update the conditional orientation prior when the prior parameters are jumped
    void restore_th_ph_prior();

    void update_proposals();   //Adapt standard deviation of proposal distributions during MCMC execution 
    void report() const;

    HRvoxel& operator=(const HRvoxel& rhs){
      m_fibres=rhs.m_fibres;
      m_d=rhs.m_d;   m_d_old=rhs.m_d_old;   m_d_prop=rhs.m_d_prop;
      m_d_prior=rhs.m_d_prior;   m_d_old_prior=rhs.m_d_old_prior; 
      m_d_acc=rhs.m_d_acc;       m_d_rej=rhs.m_d_rej;
      m_d_std=rhs.m_d_std;   m_d_std_old=rhs.m_d_std_old;   m_d_std_prop=rhs.m_d_std_prop;
      m_d_std_prior=rhs.m_d_std_prior;   m_d_std_old_prior=rhs.m_d_std_old_prior; 
      m_d_std_acc=rhs.m_d_std_acc;       m_d_std_rej=rhs.m_d_std_rej;
      m_R=rhs.m_R;   m_R_old=rhs.m_R_old;   m_R_prop=rhs.m_R_prop;
      m_R_prior=rhs.m_R_prior;   m_R_old_prior=rhs.m_R_old_prior; 
      m_R_acc=rhs.m_R_acc;       m_R_rej=rhs.m_R_rej;

      m_fsum_prior=rhs.m_fsum_prior;     m_fsum_old_prior=rhs.m_fsum_old_prior;
      m_S0=rhs.m_S0;   m_S0_old=rhs.m_S0_old;   m_S0_prop=rhs.m_S0_prop;
      m_S0_prior=rhs.m_S0_prior;   m_S0_old_prior=rhs.m_S0_old_prior; 
      m_S0_acc=rhs.m_S0_acc;       m_S0_rej=rhs.m_S0_rej;
      m_tau=rhs.m_tau;   m_tau_old=rhs.m_tau_old;   m_tau_prop=rhs.m_tau_prop;
      m_tau_prior=rhs.m_tau_prior;   m_tau_old_prior=rhs.m_tau_old_prior; 
      m_tau_acc=rhs.m_tau_acc;       m_tau_rej=rhs.m_tau_rej;
      m_prior_en=rhs.m_prior_en;   m_old_prior_en=rhs.m_old_prior_en;
      m_iso_SignalHR=rhs.m_iso_SignalHR;       m_iso_SignalHR_old=rhs.m_iso_SignalHR_old;
      m_iso_SignalLR=rhs.m_iso_SignalLR;       m_iso_SignalLR_old=rhs.m_iso_SignalLR_old;
      m_SignalHR=rhs.m_SignalHR;       m_SignalHR_old=rhs.m_SignalHR_old;
      m_SignalLR=rhs.m_SignalLR;       m_SignalLR_old=rhs.m_SignalLR_old;
      return *this;
    }
  };




  /////////////////////////////////////////////////////////////////////////////
  //     Orient_Prior_Mode: Models a single mode of the orientation prior   ///
  //     The orientation prior is a sum of Watson distributions centred     ///
  //     around different modes                                             ///
  /////////////////////////////////////////////////////////////////////////////
  class Orient_Prior_Mode{
    float m_th;
    float m_th_old;
    float m_th_prop;
    float m_th_prior;
    float m_th_old_prior;
    float m_th_acc;
    float m_th_rej;

    float m_ph;
    float m_ph_old;
    float m_ph_prop;
    float m_ph_prior;
    float m_ph_old_prior;
    float m_ph_acc;
    float m_ph_rej;

    float m_invkappa;               //Dispersion Index of a Watson distribution (1/kappa actually)
    float m_invkappa_old;
    float m_invkappa_prop;
    float m_invkappa_prior;
    float m_invkappa_old_prior;
    float m_invkappa_acc;
    float m_invkappa_rej;

    float m_Watson_norm;           //this is a function of current m_kappa
    float m_Watson_norm_old;
    ColumnVector m_vec;
    ColumnVector m_vec_old;
    float m_prior_en;              //Joint HyperPrior
    float m_old_prior_en;

    const bool m_kappa_ard;        //Flag for setting ARD on the dispersion index
  
  public:
    //Constructor
  Orient_Prior_Mode(bool kappa_ard): m_kappa_ard(kappa_ard) {
     m_th=M_PI/2; m_th_old=m_th;
     m_ph=0; m_ph_old=m_ph;
     m_vec.ReSize(3);
     m_vec<<sin(m_th)*cos(m_ph) <<sin(m_th)*sin(m_ph) <<cos(m_th);
     m_vec_old=m_vec;
     m_invkappa=0.02; m_invkappa_old=m_invkappa;
     m_Watson_norm=compute_Watson_norm(); m_Watson_norm_old=m_Watson_norm;
     m_th_prop=0.2;  m_ph_prop=0.2; m_invkappa_prop=0.2;
	
     m_th_prior=0; m_th_old_prior=0; 
     m_ph_prior=0; m_ph_old_prior=0; 
     m_invkappa_prior=0; m_invkappa_old_prior=0;
     m_prior_en=0; m_old_prior_en=0;
     
     m_th_acc=0; m_th_rej=0;
     m_ph_acc=0; m_ph_rej=0;
     m_invkappa_acc=0; m_invkappa_rej=0;
    }

    //Destructor
    ~Orient_Prior_Mode(){}

    inline float get_th() const{ return m_th;}
    inline float get_ph() const{ return m_ph;}
    inline float get_invkappa() const{ return m_invkappa;}
    inline void set_th(const float th) { m_th=th; m_vec<< sin(m_th)*cos(m_ph) << sin(m_th)*sin(m_ph) << cos(m_th);}
    inline void set_ph(const float ph) { m_ph=ph; m_vec<< sin(m_th)*cos(m_ph) << sin(m_th)*sin(m_ph) << cos(m_th);}
    inline void set_invkappa(const float invkappa) { m_invkappa=invkappa; filter_invkappa(); m_Watson_norm=compute_Watson_norm(); }
    void filter_invkappa();    //Filter out extreme values of invkappa to ensure calculations are numerically stable 

    inline float get_prior() const{ return m_prior_en;}
    inline const ColumnVector& getVec() const { return m_vec; }
    inline float get_Watson_norm() const{ return m_Watson_norm; }
    float compute_Watson_norm();
    
    void initialise_energies();
    bool compute_th_prior();
    bool compute_ph_prior(); 
    bool compute_invkappa_prior(); 
    void compute_prior();      //Compute Joint Prior energy
    void update_proposals();

    bool propose_th();
    inline void accept_th(){ m_th_acc++; }   
    void reject_th();
    
    bool propose_ph();
    inline void accept_ph(){ m_ph_acc++; }
    void reject_ph();
    
    bool propose_invkappa();
    inline void accept_invkappa(){ m_invkappa_acc++; }
    void reject_invkappa();

    Orient_Prior_Mode& operator=(const Orient_Prior_Mode& rhs){
      m_th=rhs.m_th; 
      m_th_old=rhs.m_th_old;  
      m_th_prop=rhs.m_th_prop;
      m_th_prior=rhs.m_th_prior;  m_th_old_prior=rhs.m_th_old_prior;
      m_th_acc=rhs.m_th_acc;      m_th_rej=rhs.m_th_rej;
      m_ph=rhs.m_ph; m_ph_old=rhs.m_ph_old;  m_ph_prop=rhs.m_ph_prop;
      m_ph_prior=rhs.m_ph_prior;  m_ph_old_prior=rhs.m_ph_old_prior;
      m_ph_acc=rhs.m_ph_acc;      m_ph_rej=rhs.m_ph_rej;
      m_invkappa=rhs.m_invkappa;  m_invkappa_old=rhs.m_invkappa_old;  m_invkappa_prop=rhs.m_invkappa_prop;
      m_invkappa_prior=rhs.m_invkappa_prior;  m_invkappa_old_prior=rhs.m_invkappa_old_prior;
      m_invkappa_acc=rhs.m_invkappa_acc;      m_invkappa_rej=rhs.m_invkappa_rej;
      m_Watson_norm=rhs.m_Watson_norm;        m_Watson_norm_old=rhs.m_Watson_norm_old;
      m_vec=rhs.m_vec;           m_vec_old=rhs.m_vec_old; 
      m_prior_en=rhs.m_prior_en; m_old_prior_en=rhs.m_old_prior_en;
      return *this;
    }
  };



  //////////////////////////////
  //       LRvoxel            //
  //////////////////////////////
  class LRvoxel{
    vector<HRvoxel> m_HRvoxels;
    vector<Orient_Prior_Mode> m_PModes;
    
    float m_tauLR;                //models noise at Low-res 
    float m_tauLR_old;
    float m_tauLR_prop;
    float m_tauLR_prior;
    float m_tauLR_old_prior;
    float m_tauLR_acc;
    float m_tauLR_rej;

    float m_S0LR;                 //models S0 intensity at the Low-Res acquisition
    float m_S0LR_old;
    float m_S0LR_prop;
    float m_S0LR_prior;
    float m_S0LR_old_prior;
    float m_S0LR_acc;
    float m_S0LR_rej;
    Matrix m_Orient_hyp_prior;    //Matrix Nmodes x 5 that contains the hyperparameters for the orientation prior 
                                  //columns 1-3 contains the (x,y,z) coordinates for the mode, 4th column contains the invkappa value, 5th the Watson normalization constant

    float m_mean_d;               //models mean of d hyperprior for the High-Res voxels intersected by the LR_voxel
    float m_mean_d_old;
    float m_mean_d_prop;
    float m_mean_d_prior;
    float m_mean_d_old_prior;
    float m_mean_d_acc;
    float m_mean_d_rej;

    float m_stdev_d;               //models std_dev of d hyperprior for the High-Res voxels intersected by the LR_voxel
    float m_stdev_d_old;
    float m_stdev_d_prop;
    float m_stdev_d_prior;
    float m_stdev_d_old_prior;
    float m_stdev_d_acc;
    float m_stdev_d_rej;

    float m_mean_fsum;               //models mean of fsum hyperprior for the High-Res voxels intersected by the LR_voxel
    float m_mean_fsum_old;
    float m_mean_fsum_prop;
    float m_mean_fsum_prior;
    float m_mean_fsum_old_prior;
    float m_mean_fsum_acc;
    float m_mean_fsum_rej;

    float m_stdev_fsum;               //models std_dev of fsum hyperprior for the High-Res voxels intersected by the LR_voxel
    float m_stdev_fsum_old;
    float m_stdev_fsum_prop;
    float m_stdev_fsum_prior;
    float m_stdev_fsum_old_prior;
    float m_stdev_fsum_acc;
    float m_stdev_fsum_rej;

    float m_prior_en;             //Joint Prior Energy
    float m_old_prior_en;
    float m_likelihood_en;        //Likelihood Energy
    float m_old_likelihood_en;
    float m_posterior_en;         //Posterior Energy
    float m_old_posterior_en;
   
    ColumnVector m_logdataLR;       //Log of Low-Res Data for the specific LR voxel (use it in Rician Energy)
    vector<ColumnVector> m_logdataHR; //Log of High-Res Data for all contained HR voxels (use it in Rician Energy)

    const ColumnVector& m_dataLR;   //Low-Res Data for the specific LR voxel 
    const vector<ColumnVector>& m_dataHR; //High-Res Data for all HRvoxels within a LRvoxel
    const vector<Matrix>& m_bvecsHR; //bvecs at High-Res   (HRvoxels within a LRvoxel x 3 x HR_NumPoints)
    const vector<Matrix>& m_bvalsHR; //bvalues at High-Res (HRvoxels within a LRvoxel x 3 x HR_NumPoints)
    const Matrix& m_bvecsLR;        //bvecs at Low-Res    (3 x LR_NumPoints)
    const Matrix& m_bvalsLR;        //bvalues at Low-Res  (1 x HR_NumPoints)
    /////////////////// User-defined parameters follow  ///////////////////////
    const int m_modelnum;           //1 for single-shell, 2 for multi-shell model, 3 for zeppelins
    const int m_PVmodelnum;         //Patial volume model, 1: for sum of attenuations, 2: for attenuation of sums
    const int m_numfibres;          //Number of fibres in each HR voxel
    const int m_Nmodes;             //Number of modes for the Orientation Prior
    const float m_ardfudge;
    const bool m_allard;            //Flag for setting ARD on for all fibres in all HR voxels
    const bool m_Noard;             //Flag for setting ARD off for all fibres in all HR voxels
    const bool m_kappa_ard;         //Flag for setting ARD on the dispersion index
    const bool m_fsumPrior_ON;      //Flag for setting on a prior on fsums across intersected HR voxels
    const bool m_dPrior_ON;         //Flag for setting on a prior on the diffusivity across intersected HR voxels
    const bool m_rician;            //Flag for using a Rician noise model 
    const bool m_noS0jump;          //Indicates whether S0 parameters will be kept constant during MCMC 
    const float m_R_priormean;      //Parameters to use for the Gaussian prior on R. Mean
    const float m_R_priorstd;       //and variance
    const float m_R_priorfudge;     //and fudge
    const ColumnVector& m_HRweights;//Holds the volume fraction each HR voxel occupies out of the LR one
 
  public:
    //Constructor
    LRvoxel(const vector<Matrix>& bvecsHR, const vector<Matrix>& bHR, 
	    const Matrix& bvecsLR, const Matrix& bLR, 
	    const ColumnVector& dataLR, const vector<ColumnVector>& dataHR, const int N, const int Nmodes, const ColumnVector& HRweights, const int modelnum=1, const int PVmodelnum=1, const float ardfudge=1, const bool allard=false, const bool Noard=false, const bool kappa_ard=false, const bool fsumPrior_ON=false, const bool dPrior_ON=false, const bool rician=false, const bool noS0jump=false, const float Rmean=0.13, const float Rstd=0.03,const float Rfudge=0):
      m_dataLR(dataLR), m_dataHR(dataHR), m_bvecsHR(bvecsHR), m_bvalsHR(bHR), m_bvecsLR(bvecsLR), m_bvalsLR(bLR),
	m_modelnum(modelnum), m_PVmodelnum(PVmodelnum),m_numfibres(N), m_Nmodes(Nmodes),m_ardfudge(ardfudge), m_allard(allard), m_Noard(Noard), m_kappa_ard(kappa_ard), m_fsumPrior_ON(fsumPrior_ON), m_dPrior_ON(dPrior_ON), m_rician(rician), m_noS0jump(noS0jump),  m_R_priormean(Rmean), m_R_priorstd(Rstd), m_R_priorfudge(Rfudge), m_HRweights(HRweights) {
    
      m_S0LR=0; m_S0LR_old=0; m_S0LR_prior=0; m_S0LR_old_prior=0; m_S0LR_acc=0; m_S0LR_rej=0; m_S0LR_prop=0.2;
      m_tauLR=0; m_tauLR_old=0; m_tauLR_prior=0; m_tauLR_old_prior=0; m_tauLR_acc=0; m_tauLR_rej=0; m_tauLR_prop=0.2;

      m_mean_d=0; m_mean_d_old=0; m_mean_d_prior=0; m_mean_d_old_prior=0; m_mean_d_acc=0; m_mean_d_rej=0; m_mean_d_prop=0.001;
      m_stdev_d=0; m_stdev_d_old=0; m_stdev_d_prior=0; m_stdev_d_old_prior=0; m_stdev_d_acc=0; m_stdev_d_rej=0; m_stdev_d_prop=0.001;

      m_mean_fsum=0; m_mean_fsum_old=0; m_mean_fsum_prior=0; m_mean_fsum_old_prior=0; m_mean_fsum_acc=0; m_mean_fsum_rej=0; m_mean_fsum_prop=0.1;
      m_stdev_fsum=0; m_stdev_fsum_old=0; m_stdev_fsum_prior=0; m_stdev_fsum_old_prior=0; m_stdev_fsum_acc=0; m_stdev_fsum_rej=0; m_stdev_fsum_prop=0.1;

      m_prior_en=0; m_old_prior_en=0; 
      m_likelihood_en=0; m_old_likelihood_en=0;
      m_posterior_en=0; m_old_posterior_en=0;

      for (int m=1; m<=m_Nmodes; m++){                //Add Modes for the orientation Prior
      	Orient_Prior_Mode pMod(m_kappa_ard); 
      	m_PModes.push_back(pMod);
      }
      m_Orient_hyp_prior.ReSize(m_Nmodes,5);

      for (unsigned int n=0; n<m_dataHR.size(); n++){ //Add HRvoxel Objects
      	HRvoxel HRv(m_bvecsHR[n], m_bvalsHR[n], m_bvecsLR, m_bvalsLR, m_Orient_hyp_prior, m_mean_d, m_stdev_d, m_mean_fsum, m_stdev_fsum, m_R_priormean, m_R_priorstd, m_R_priorfudge, m_numfibres, m_modelnum, m_ardfudge, m_fsumPrior_ON, m_dPrior_ON, m_rician, m_noS0jump);
      	
	m_HRvoxels.push_back(HRv);
      }

      // ofstream myfile;     //Debugging code
      //myfile.open("/Users/stam/Rubix_data/Energies.txt", ios::trunc);
      //myfile.close();
  
      if (m_rician){//Store the log of the data for energy calculations
	m_logdataLR=dataLR; m_logdataHR=dataHR;
	for (int m=1; m<=dataLR.Nrows(); m++){  
	  if (dataLR(m)!=0)
	    m_logdataLR(m)=log(dataLR(m));
	  else
	    m_logdataLR(m)=minlogfloat;
	}
	for (unsigned int n=0; n<dataHR.size(); n++)
	  for (int m=1; m<=dataHR[n].Nrows(); m++){
	    if (dataHR[n](m)!=0)
	      m_logdataHR[n](m)=log(dataHR[n](m));
	    else
	      m_logdataHR[n](m)=minlogfloat;
	  }
      }
    }

    //Destructor
      ~LRvoxel(){ }
    
    const vector<HRvoxel>& HRvoxels() const {return m_HRvoxels;}
    const vector<Orient_Prior_Mode>& PModes() const {return m_PModes;} 

    void update_Orient_hyp_prior();       //Update the matrix that keeps information on the orientation prior parameters
    void update_Orient_hyp_prior(int M);  //Update the entry only for a specific Mode 0<=M<N_modes
    void set_HRparams(const int n, const float d, const float S0, const ColumnVector& th, const ColumnVector& ph, const ColumnVector& f);  //Set params for a single HR voxel
    void set_HRparams(const int n, const float d, const float d_std, const float S0, const ColumnVector& th, const ColumnVector& ph, const ColumnVector& f);  
    void set_HRparams(const int n, const float d, const float d_std, const float S0, const float R, const ColumnVector& th, const ColumnVector& ph, const ColumnVector& f); 
    void set_Priorparams(const ColumnVector& th, const ColumnVector& ph, const ColumnVector& invkappa);  //Set params for all modes of orientation priors
    float get_S0LR() const { return m_S0LR; }
    void set_S0LR(const float S0)  { m_S0LR=S0; }
    float get_tauLR() const { return m_tauLR; }
    void set_tauLR(const float tau)  { m_tauLR=tau; }
    void set_tauHR(const int n, const float tau)  { m_HRvoxels[n].set_tau(tau); }
    float get_likelihood_energy() const { return m_likelihood_en; }
    float get_prior_energy() const { return m_prior_en; }
    void initialise_energies();
    
    bool propose_S0LR();
    void accept_S0LR(){ m_S0LR_acc++; }
    void reject_S0LR();
    bool compute_S0LR_prior();

    bool propose_tauLR();
    void accept_tauLR(){ m_tauLR_acc++; }
    void reject_tauLR();
    bool compute_tauLR_prior();

    void compute_prior();
    void compute_likelihood();
    void compute_posterior();
    bool test_energy() const;
    void restore_energies();
    void restore_Prior_Posterior();
    void update_proposals();
    void jump(); //Single MCMC iteration with all parameters jumped

    
    bool propose_meand();
    void accept_meand(){ m_mean_d_acc++; }
    void reject_meand();
    bool compute_meand_prior();
    bool propose_stdevd();
    void accept_stdevd(){ m_stdev_d_acc++; }
    void reject_stdevd();
    bool compute_stdevd_prior();
    void set_meand(const float meand)  { m_mean_d=meand; }
    void set_stdevd(const float stdevd)  { m_stdev_d=stdevd; }
    float get_meand() const{ return m_mean_d;}
    float get_stdevd() const{ return m_stdev_d;}

    bool propose_mean_fsum();
    void accept_mean_fsum(){ m_mean_fsum_acc++; }
    void reject_mean_fsum();
    bool compute_mean_fsum_prior();
    bool propose_stdev_fsum();
    void accept_stdev_fsum(){ m_stdev_fsum_acc++; }
    void reject_stdev_fsum();
    bool compute_stdev_fsum_prior();
    void set_mean_fsum(const float fsum)  { m_mean_fsum=fsum; }
    void set_stdev_fsum(const float stdev_fsum)  { m_stdev_fsum=stdev_fsum; }
    float get_mean_fsum() const{ return m_mean_fsum;}
    float get_stdev_fsum() const{ return m_stdev_fsum;}

    
    LRvoxel& operator=(const LRvoxel& rhs){    
      m_HRvoxels=rhs.m_HRvoxels;            m_PModes=rhs.m_PModes;
      m_tauLR=rhs.m_tauLR;                  m_tauLR_old=rhs.m_tauLR_old;
      m_tauLR_prop=rhs.m_tauLR_prop;
      m_tauLR_prior=rhs.m_tauLR_prior;      m_tauLR_old_prior=rhs.m_tauLR_old_prior;
      m_tauLR_acc=rhs.m_tauLR_acc;          m_tauLR_rej=rhs.m_tauLR_rej;
      m_S0LR=rhs.m_S0LR;                    m_S0LR_old=rhs.m_S0LR_old;
      m_S0LR_prop=rhs.m_S0LR_prop;
      m_S0LR_prior=rhs.m_S0LR_prior;        m_S0LR_old_prior=rhs.m_S0LR_old_prior;
      m_S0LR_acc=rhs.m_S0LR_acc;            m_S0LR_rej=rhs.m_S0LR_rej;
      m_Orient_hyp_prior=rhs.m_Orient_hyp_prior;
      m_prior_en=rhs.m_prior_en;            m_old_prior_en=rhs.m_old_prior_en;
      m_likelihood_en=rhs.m_likelihood_en;  m_old_likelihood_en=rhs.m_old_likelihood_en;
      m_posterior_en=rhs.m_posterior_en;    m_old_posterior_en=rhs.m_old_posterior_en;
      m_logdataLR=rhs.m_logdataLR;          m_logdataHR=rhs.m_logdataHR;  

      m_mean_d=rhs.m_mean_d;                m_mean_d_old=rhs.m_mean_d_old;
      m_mean_d_prop=rhs.m_mean_d_prop;
      m_mean_d_prior=rhs.m_mean_d_prior;    m_mean_d_old_prior=rhs.m_mean_d_old_prior;
      m_mean_d_acc=rhs.m_mean_d_acc;        m_mean_d_rej=rhs.m_mean_d_rej;
      m_stdev_d=rhs.m_stdev_d;              m_stdev_d_old=rhs.m_stdev_d_old;
      m_stdev_d_prop=rhs.m_stdev_d_prop;
      m_stdev_d_prior=rhs.m_stdev_d_prior;  m_stdev_d_old_prior=rhs.m_stdev_d_old_prior;
      m_stdev_d_acc=rhs.m_stdev_d_acc;      m_stdev_d_rej=rhs.m_stdev_d_rej;

      m_mean_fsum=rhs.m_mean_fsum;                m_mean_fsum_old=rhs.m_mean_fsum_old;
      m_mean_fsum_prop=rhs.m_mean_fsum_prop;
      m_mean_fsum_prior=rhs.m_mean_fsum_prior;    m_mean_fsum_old_prior=rhs.m_mean_fsum_old_prior;
      m_mean_fsum_acc=rhs.m_mean_fsum_acc;        m_mean_fsum_rej=rhs.m_mean_fsum_rej;
      m_stdev_fsum=rhs.m_stdev_fsum;              m_stdev_fsum_old=rhs.m_stdev_fsum_old;
      m_stdev_fsum_prop=rhs.m_stdev_fsum_prop;
      m_stdev_fsum_prior=rhs.m_stdev_fsum_prior;  m_stdev_fsum_old_prior=rhs.m_stdev_fsum_old_prior;
      m_stdev_fsum_acc=rhs.m_stdev_fsum_acc;      m_stdev_fsum_rej=rhs.m_stdev_fsum_rej;
     
      return *this;
    }
  };



}

#endif
