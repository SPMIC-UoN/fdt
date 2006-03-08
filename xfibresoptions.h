/*  BpmOptions.h

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(xfibresOptions_h)
#define xfibresOptions_h

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

namespace Xfibres {

class xfibresOptions {
 public:
  static xfibresOptions& getInstance();
  ~xfibresOptions() { delete gopt; }
  
  Option<bool> verbose;
  Option<bool> help;
  Option<string> logdir;
  Option<bool> forcedir;
  Option<string> datafile;
  Option<string> maskfile;
  Option<string> bvecsfile;
  Option<string> bvalsfile;
  Option<int> nfibres;
  Option<float> fudge;
  Option<int> njumps;
  Option<int> nburn;
  Option<int> nburn_noard;
  Option<int> sampleevery;
  Option<int> updateproposalevery;
  Option<int> seed;
  Option<bool> no_ard;
  void parse_command_line(int argc, char** argv,  Log& logger);
  
 private:
  xfibresOptions();  
  const xfibresOptions& operator=(xfibresOptions&);
  xfibresOptions(xfibresOptions&);

  OptionParser options; 
      
  static xfibresOptions* gopt;
  
};

 inline xfibresOptions& xfibresOptions::getInstance(){
   if(gopt == NULL)
     gopt = new xfibresOptions();
   
   return *gopt;
 }

 inline xfibresOptions::xfibresOptions() :
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
  datafile(string("-k,--data,--datafile"), string("data"),
	      string("data file"),
	      true, requires_argument),  
  maskfile(string("-m,--mask, --maskfile"), string("nodif_brain_mask"),
	    string("mask file"),
	    true, requires_argument),
  bvecsfile(string("-r,--bvecs"), string("bvecs"),
	     string("b vectors file"),
	     true, requires_argument),  
  bvalsfile(string("-b,--bvals"), string("bvals"),
	     string("b values file"),
	     true, requires_argument), 
  nfibres(string("--nf,--nfibres"),1,
	 string("Maximum nukmber of fibres to fit in each voxel (default 1)"),
	 false,requires_argument),
  fudge(string("--fudge"),1,
	 string("ARD fudge factor"),
	 false,requires_argument),
  njumps(string("--nj,--njumps"),5000,
	 string("Num of jumps to be made by MCMC (default is 5000)"),
	 false,requires_argument),
  nburn(string("--bi,--burnin"),1,
	string("Total num of jumps at start of MCMC to be discarded"),
	false,requires_argument),
  nburn_noard(string("--bn,--burnin_noard"),0,
	string("num of burnin jumps before the ard is imposed"),
	false,requires_argument),
  sampleevery(string("--se,--sampleevery"),1,
	string("Num of jumps for each sample (MCMC) (default is 1)"),
	false,requires_argument),
  updateproposalevery(string("--upe,--updateproposalevery"),40,
	string("Num of jumps for each update to the proposal density std (MCMC) (default is 40)"),
	false,requires_argument),
  seed(string("--seed"),8665904,string("seed for pseudo random number generator"),
       false,requires_argument),
  no_ard(string("--noard"),false,string("Turn ARD off on all fibres"),
       false,no_argument),
   options("xfibres", "xfibres -k <filename>\n xfibres --verbose\n")
   {
     
    
     try {
       options.add(verbose);
       options.add(help);
       options.add(logdir);
       options.add(forcedir);
       options.add(datafile);
       options.add(maskfile);
       options.add(bvecsfile);
       options.add(bvalsfile);
       options.add(nfibres);
       options.add(fudge);
       options.add(njumps);
       options.add(nburn);
       options.add(nburn_noard);
       options.add(sampleevery);
       options.add(updateproposalevery);
       options.add(seed);
       options.add(no_ard);
       
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





