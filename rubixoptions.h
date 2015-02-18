/*  rubixoptions.h

    Stam Sotiropoulos, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(rubixptions_h)
#define rubixptions_h

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"
#include "utils/tracer_plus.h"
//#include "newmatall.h"
using namespace Utilities;

namespace RUBIX{
  
  class rubixOptions {
  public:
    static rubixOptions& getInstance();
    ~rubixOptions() { delete gopt; }
  
    Option<bool> verbose;
    Option<bool> help;
    Option<string> logdir;
    Option<bool> forcedir;
    Option<string> LRdatafile;
    Option<string> LRmaskfile;
    Option<string> LRbvecsfile;
    Option<string> LRbvalsfile;
    Option<string> HRdatafile;
    //    Option<string> HRmaskfile;
    Option<string> HRbvecsfile;
    Option<string> HRbvalsfile;
    Option<int> nfibres;
    Option<int> nmodes;  //number of modes in the orientation prior
    Option<int> modelnum;
    Option<float> fudge;
    Option<int> njumps;
    Option<int> nburn;
    Option<int> sampleevery;
    Option<int> updateproposalevery;
    Option<int> seed;
    Option<bool> no_ard;
    Option<bool> all_ard;
    Option<bool> kappa_ard;
    Option<bool> fsumPrior;
    Option<bool> dPrior;
    Option<bool> rician;
    Option<bool> noS0jump;
    Option<string> LRgrad_file;
    Option<string> HRgrad_file;
    Option<float> R_prior_mean;  //setting the prior for model's 3 ratio of perp. to parallel diffusivity
    Option<float> R_prior_std;
    Option<float> R_prior_fudge; //if used (set to a positive number), an ARD for R is used for the high diffusivity regions with the requested fudge factor

    void parse_command_line(int argc, char** argv,  Log& logger);
  
  private:
    rubixOptions();  
    const rubixOptions& operator=(rubixOptions&);
    rubixOptions(rubixOptions&);

    OptionParser options; 
    static rubixOptions* gopt;
  };

  inline rubixOptions& rubixOptions::getInstance(){
    if(gopt == NULL)
      gopt = new rubixOptions();
   
    return *gopt;
  }

  inline rubixOptions::rubixOptions() :
   verbose(string("-V,--verbose"), false, 
	   string("switch on diagnostic messages"), 
	   false, no_argument),
   help(string("-h,--help"), false,
	string("display this message"),
	false, no_argument),
   logdir(string("--ld,--logdir"), string("logdir"),
	  string("log directory (default is logdir)"),
	  false, requires_argument),
   forcedir(string("--forcedir"),false,string("Use the actual directory name given - i.e. don't add + to make a new directory"),false,no_argument),
   LRdatafile(string("--dLR"), string("data_LR"),
	    string("Low-Res data file"),
	    true, requires_argument),  
   LRmaskfile(string("--mLR"), string("nodif_brain_mask_LR"),
	    string("Low-Res mask file"),
	    true, requires_argument),
   LRbvecsfile(string("--rLR"), string("bvecs_LR"),
	     string("Low-Res b vectors file"),
	     true, requires_argument),  
   LRbvalsfile(string("--bLR"), string("bvals_LR"),
	     string("Low-Res b values file"),
	     true, requires_argument),
   HRdatafile(string("--dHR"), string("data_HR"),
	    string("High-Res data file"),
	    true, requires_argument),  
   //  HRmaskfile(string("--mHR"), string("nodif_brain_mask_HR"),
   //	    string("High-Res mask file"),
   //	    true, requires_argument),
   HRbvecsfile(string("--rHR"), string("bvecs_HR"),
	     string("High-Res b vectors file"),
	     true, requires_argument),  
   HRbvalsfile(string("--bHR"), string("bvals_HR"),
	     string("High-Res b values file"),
	     true, requires_argument),
   nfibres(string("--nf,--nfibres"),1,
	   string("Maximum number of fibres to fit in each HR voxel (default 1)"),
	   false,requires_argument),
   nmodes(string("--nM,--nmodes"),2,
	   string("Number of modes for the orientation prior (default 2)"),
	   false,requires_argument),
   modelnum(string("--model"),1,
	    string("\tWhich deconvolution model to use. 1:With sticks (default), 2:With sticks and a range of diffusivities, 3:With zeppelins"),
	    false,requires_argument),
   fudge(string("--fudge"),1,
	 string("\tARD fudge factor"),
	 false,requires_argument),
   njumps(string("--nj,--njumps"),1500,
	  string("Num of jumps to be made by MCMC (default is 1500)"),
	  false,requires_argument),
   nburn(string("--bi,--burnin"),0,
	 string("Total num of jumps at start of MCMC to be discarded (default is 0)"),
	 false,requires_argument),
   sampleevery(string("--se"),30,
	       string("\tNum of jumps for each sample (MCMC) (default is 30)"),
	       false,requires_argument),
   updateproposalevery(string("--upe"),40,
		       string("\tNum of jumps for each update to the proposal density std (MCMC) (default is 40)"),
		       false,requires_argument),
   seed(string("--seed"),8665904,string("\tseed for pseudo random number generator"),
	false,requires_argument),
   no_ard(string("--noard"),false,string("\tTurn ARD off on all fibres"),
    	  false,no_argument),
   all_ard(string("--allard"),false,string("Turn ARD on for all fibres (default: on for all but first)"),
    	   false,no_argument),
   kappa_ard(string("--ard_kappa"),false,string("Turn ARD on for the dispersion of all orientation prior modes (default: off)"),
    	   false,no_argument),
   fsumPrior(string("--fsumPrior"),false,string("Turn on prior for the fsums across an LR voxel (default: off)"),
    	   false,no_argument),
   dPrior(string("--dPrior"),false,string("Turn on prior for the diffusivity across an LR voxel  (default: off)"),
    	   false,no_argument),
   rician(string("--rician"),false,string("Use Rician Noise model (default is Gaussian)"),
    	   false,no_argument),
   noS0jump(string("--noS0jump"),false,string("Do not jump S0 parameters and keep them to their ML estimates (default: off)"),
    	   false,no_argument),
   LRgrad_file(string("--gLR"), string("grad_devLR"),
	     string("\tLR Gradient Nonlinearity Tensor"),
	     false, requires_argument),  
   HRgrad_file(string("--gHR"), string("grad_devHR"),
	     string("\tHR Gradient Nonlinearity Tensor"),
	     false, requires_argument),  
   R_prior_mean(string("--Rmean"),0.13,string("\tSet the prior mean for R of model3 (default:0.13- Must be<0.5)"),false, requires_argument),
   R_prior_std(string("--Rstd"),0.03,string("\tSet the prior standard deviation for R of model3 (default:0.03)"),false, requires_argument),
   R_prior_fudge(string("--Rfudge"),0,string("If set(>0), an ARD prior is used for R with the requested fudge factor"),false, requires_argument),
   options("RubiX v1.0", "rubix --help (for list of options)\n")
     {
       try {
       options.add(verbose);
       options.add(help);
       options.add(logdir);
       options.add(forcedir);
       options.add(LRdatafile);
       options.add(LRmaskfile);
       options.add(LRbvecsfile);
       options.add(LRbvalsfile);
       options.add(HRdatafile);
       //options.add(HRmaskfile);
       options.add(HRbvecsfile);
       options.add(HRbvalsfile);
       options.add(nfibres);
       options.add(nmodes);
       options.add(modelnum);
       options.add(fudge);
       options.add(njumps);
       options.add(nburn);
       options.add(sampleevery);
       options.add(updateproposalevery);
       options.add(seed);
       options.add(no_ard);
       options.add(all_ard);
       options.add(kappa_ard);
       options.add(fsumPrior);
       options.add(dPrior);
       options.add(rician);
       options.add(noS0jump);
       options.add(LRgrad_file);
       options.add(HRgrad_file);
       options.add(R_prior_mean);
       options.add(R_prior_std);
       options.add(R_prior_fudge);
     }
     catch(X_OptionError& e) {
       options.usage();
       cerr << endl << e.what() << endl;
     } 
     catch(std::exception &e) {
       cerr << e.what() << endl;
     }    
   }
}

#endif
