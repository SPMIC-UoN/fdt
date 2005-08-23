/*  probtrackOptions.h

    Tim Behrens, FMRIB Image Analysis Group

    Copyright (C) 2004 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(probtrackOptions_h)
#define probtrackOptions_h

#include <string>
#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"

//#include "newmatall.h"
using namespace Utilities;

namespace TRACT {

class probtrackOptions {
 public:
  static probtrackOptions& getInstance();
  ~probtrackOptions() { delete gopt; }
  
  Option<int> verbose;
  Option<bool> help;
  Option<string> basename;
  Option<string> maskfile;
  Option<string> seedfile; 
  Option<string> mode;
  Option<string> targetfile;
  Option<string> outfile;
  Option<string> rubbishfile;
  FmribOption<string> uberrubbishfile;
  Option<string> seeds_to_dti;
  FmribOption<string> skipmask;
  Option<string> seedref;
  Option<string> mask2;
  Option<string> meshfile;
  FmribOption<string> lrmask;
  Option<string> logdir; 
  Option<bool> forcedir;
  Option<int> nparticles;
  Option<int> nsteps;
  Option<float> c_thr;
  Option<float> steplength;
  Option<bool> loopcheck;
  Option<bool> usef;
  Option<bool> modeuler;
  Option<int> rseed;
  void parse_command_line(int argc, char** argv,Log& logger);
  void modecheck();
  void modehelp();
  void matrixmodehelp();
  void status();
 private:
  probtrackOptions();  
  const probtrackOptions& operator=(probtrackOptions&);
  probtrackOptions(probtrackOptions&);

  OptionParser options; 
      
  static probtrackOptions* gopt;
  
};


 inline probtrackOptions& probtrackOptions::getInstance(){
   if(gopt == NULL)
     gopt = new probtrackOptions();
   
   return *gopt;
 }

 inline probtrackOptions::probtrackOptions() :
  verbose(string("-V,--verbose"), 0, 
	  string("verbose level, [0-2]"), 
	  false, requires_argument),
   help(string("-h,--help"), false,
	string("display this message"),
	false, no_argument),
   basename(string("-s,--samples"), string("DTI"),
	       string("basename for samples files"),
	       true, requires_argument),  
   maskfile(string("-m,--mask"), string("mask"),
	    string("Bet binary mask file in diffusion space"),
	    true, requires_argument),
   seedfile(string("-x,--seed"), string("Seed"),
	    string("Seed volume or voxel depending on mode"),
	    true, requires_argument),
   mode(string("--mode"), string("simple"),
	string("tracking mode -  To see list of options and explanations, use --mode=help"),
	    true, requires_argument),
   targetfile(string("--targetmasks"), string("cmasks"),
	    string("File containing a list of target masks - required for seeds_to_targets mode"),
	    false, requires_argument),
   outfile(string("-o,--out"), string(""),
	    string("Output file"),
	    false, requires_argument),
   rubbishfile(string("--rubbish"), string(""),
	    string("Rubbish file"),
	    false, requires_argument),
  uberrubbishfile(string("--rubbishkill"), string(""),
	    string("Rubbish and kill file - does not include the _whole path_ of a streamline that reaches this maskfile"),
	    false, requires_argument),
   seeds_to_dti(string("--xfm"), string(""),
	      string("Transform Matrix taking seed space to DTI space default is to use the identity"),false, requires_argument),
  skipmask(string("--no_integrity"), string(""),
	   string("no explanation needed"),
	   false, requires_argument),
  seedref(string("--seedref"), string(""),
	 string("Reference vol to define seed space in simple mode - diffusion space assumed if absent"),
	 false, requires_argument),
  mask2(string("--mask2"), string(""),
	 string("second mask in twomask_symm mode. Waypoint mask or ascii list of masks in wapoints mode."),
       false, requires_argument),
  meshfile(string("--mesh"), string(""),
	 string(""),
       false, requires_argument),
  lrmask(string("--lrmask"), string(""),
	 string("low resolution binary brain mask for stroring connectivity distribution in matrix2 mode."),
       false, requires_argument),
  logdir(string("--dir"), string(""),
	    string("Directory to put the final volumes in - code makes this directory"),
	    false, requires_argument),
  forcedir(string("--forcedir"), false,
	 string("Use the actual directory name given - i.e. don't add + to make a new directory"),
	 false, no_argument),
  nparticles(string("-P,--nsamples"), 10000,
	 string("Number of samples"),
	 false, requires_argument),
   nsteps(string("-S,--nsteps"), 2000,
	    string("Number of steps per sample"),
	    false, requires_argument),
   c_thr(string("-c,--cthr"), 0.2, 
	 string("Curvature threshold"), 
	 false, requires_argument),
   steplength(string("--steplength"), 0.5, 
	 string("steplength"), 
	 false, requires_argument),
   loopcheck(string("-l,--loopcheck"), false, 
	 string("perform loopchecks on paths - slower, but allows lower curvature threshold"), 
	 false, no_argument),
   usef(string("-f,--usef"), false, 
	 string("Use anisotropy to constrain tracking"), 
	 false, no_argument),
  modeuler(string("--modeuler"), false, 
	   string("Use modified euler streamlining"), 
	   false, no_argument),
   rseed(string("--rseed"), 12345,
	 string("Random seed"),
	 false, requires_argument), 
   options("probtrack","probtrack -s <basename> -m <maskname> -x <seedfile> -o <output> --targetmasks=<textfile>\n probtrack --help\n")
   {
     
    
     try {
       options.add(verbose);
       options.add(help);
       options.add(basename);
       options.add(maskfile);
       options.add(seedfile); 
       options.add(mode);
       options.add(targetfile);
       options.add(skipmask);
       options.add(mask2);
       options.add(meshfile);
       options.add(lrmask);
       options.add(seedref);
       options.add(logdir); 
       options.add(forcedir); 
       options.add(outfile);
       options.add(rubbishfile);
       options.add(uberrubbishfile);
       options.add(seeds_to_dti);
       options.add(nparticles);
       options.add(nsteps);
       options.add(c_thr);
       options.add(steplength);
       options.add(loopcheck);
       options.add(usef);
       options.add(modeuler);
       options.add(rseed);
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







