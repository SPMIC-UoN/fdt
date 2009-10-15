/*  probtrackxOptions.h

    Tim Behrens, FMRIB Image Analysis Group

    Copyright (C) 2004 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(probtrackxOptions_h)
#define probtrackxOptions_h

#include <string> 
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"
#include "commonopts.h"

//#include "newmatall.h"
using namespace Utilities;

namespace TRACT {

class probtrackxOptions {
 public:
  static probtrackxOptions& getInstance();
  ~probtrackxOptions() { delete gopt; }
  
  Option<int> verbose;
  Option<bool> help;
  Option<string> basename;
  Option<string> maskfile;
  Option<string> seedfile; 
  Option<string> mode;
  Option<string> targetfile;
  Option<bool> simpleout;
  Option<bool> pathdist;
  Option<bool> s2tout;
  FmribOption<bool> matrix1out;
  Option<bool> matrix2out;
  FmribOption<bool> matrix3out;
  FmribOption<string> maskmatrix3;
  FmribOption<bool> maskmatrixout;
  Option<string> outfile;
  Option<string> rubbishfile;
  Option<string> stopfile;
  Option<string> prefdirfile;
  Option<string> seeds_to_dti;
  Option<string> dti_to_seeds;
  FmribOption<string> skipmask;
  Option<string> seedref;
  Option<string> mask2;
  Option<string> waypoints;
  Option<bool> network;
  Option<string> meshfile;
  Option<string> lrmask;
  Option<string> logdir; 
  Option<bool> forcedir;
  Option<int> nparticles;
  Option<int> nsteps;
  Option<float> distthresh;
  Option<float> c_thr;
  FmribOption<float> fibthresh;
  Option<float> steplength;
  Option<bool> loopcheck;
  Option<bool> usef;
  Option<int> randfib;
  Option<int> fibst;
  Option<bool> modeuler;
  Option<int> rseed;
  Option<bool> seedcountastext;
  Option<bool> splitmatrix2;

  void parse_command_line(int argc, char** argv,Log& logger);
  void modecheck();
  void modehelp();
  void matrixmodehelp();
  void status();
 private:
  probtrackxOptions();  
  const probtrackxOptions& operator=(probtrackxOptions&);
  probtrackxOptions(probtrackxOptions&);

  OptionParser options; 
      
  static probtrackxOptions* gopt;
  
};


 inline probtrackxOptions& probtrackxOptions::getInstance(){
   if(gopt == NULL)
     gopt = new probtrackxOptions();
   
   return *gopt;
 }

 inline probtrackxOptions::probtrackxOptions() :
  verbose(string("-V,--verbose"), 0, 
	  string("verbose level, [0-2]"), 
	  false, requires_argument),
   help(string("-h,--help"), false,
	string("display this message"),
	false, no_argument),
   basename(string("-s,--samples"), string("merged"),
	       string("basename for samples files"),
	       true, requires_argument),  
   maskfile(string("-m,--mask"), string("mask"),
	    string("Bet binary mask file in diffusion space"),
	    true, requires_argument),
   seedfile(string("-x,--seed"), string("Seed"),
	    string("Seed volume, or voxel, or ascii file with multiple volumes, or freesurfer label file"),
	    true, requires_argument),
   mode(string("--mode"), string(""),
	string("use --mode=simple for single seed voxel"),
	    false, requires_argument),
  targetfile(string("--targetmasks"), string("cmasks"),
	    string("File containing a list of target masks - required for seeds_to_targets classification"),
	    false, requires_argument),
  simpleout(string("--opd"), false,
	    string("output path distribution"),
	    false, no_argument), 
  pathdist(string("--pd"), false,
	   string("Correct path distribution for the length of the pathways"),
	   false, no_argument), 
  s2tout(string("--os2t"), false,
	 string("output seeds to targets"),
	 false, no_argument),
  matrix1out(string("--omatrix1"), false,
	  string("output matrix1"),
	  false, no_argument), 
  matrix2out(string("--omatrix2"), false,
	  string("output matrix2"),
	  false, no_argument), 
  matrix3out(string("--omatrix3"), false,
	  string("output matrix3 (NxN connectivity matrix)"),
	  false, no_argument), 
  maskmatrix3(string("--mask3"), "",
	  string("mask used for NxN connectivity matrix"),
	  false, requires_argument), 
  maskmatrixout(string("--omaskmatrix"), false,
		string("output maskmatrix"),
		false, no_argument), 
   outfile(string("-o,--out"), string("fdt_paths"),
	   string("Output file (default='fdt_paths')"),
	   false, requires_argument),
   rubbishfile(string("--avoid"), string(""),
	       string("Reject pathways passing through locations given by this mask"),
	       false, requires_argument),
   stopfile(string("--stop"), string(""),
	       string("Stop tracking at locations given by this mask file"),
	       false, requires_argument),
   prefdirfile(string("--prefdir"), string(""),
	       string("prefered orientation preset in a 4D mask"),
	       false, requires_argument),
   seeds_to_dti(string("--xfm"), string(""),
		string("Transform taking seed space to DTI space (either FLIRT matrix or FNIRT warpfield) - default is identity"),
		false, requires_argument),
   dti_to_seeds(string("--invxfm"), string(""),
		string("Transform taking DTI space to seed space (compulsory when using a warpfield for seeds_to_dti)"),
		false, requires_argument),
   skipmask(string("--no_integrity"), string(""),
	   string("no explanation needed"),
	   false, requires_argument),
  seedref(string("--seedref"), string(""),
	 string("Reference vol to define seed space in simple mode - diffusion space assumed if absent"),
	 false, requires_argument),
  mask2(string("--mask2"), string(""),
	 string("second mask in twomask_symm mode."),
       false, requires_argument),
 waypoints(string("--waypoints"), string(""),
	 string("Waypoint mask or ascii list of waypoint masks - only keep paths going through ALL the masks"),
       false, requires_argument),
 network(string("--network"), false,
	 string("Activate network mode - only keep paths going through at least one seed mask (required if multiple seed masks)"),
       false, no_argument),
   meshfile(string("--mesh"), string(""),
	 string("Freesurfer-type surface descriptor (in ascii format)"),
       false, requires_argument),
  lrmask(string("--lrmask"), string(""),
	 string("low resolution binary brain mask for stroring connectivity distribution in matrix2 mode"),
       false, requires_argument),
  logdir(string("--dir"), string("logdir"),
	    string("Directory to put the final volumes in - code makes this directory - default='logdir'"),
	    false, requires_argument),
  forcedir(string("--forcedir"), false,
	 string("Use the actual directory name given - i.e. don't add + to make a new directory"),
	 false, no_argument),
  nparticles(string("-P,--nsamples"), 5000,
	 string("Number of samples - default=5000"),
	 false, requires_argument),
   nsteps(string("-S,--nsteps"), 2000,
	    string("Number of steps per sample - default=2000"),
	    false, requires_argument),
   distthresh(string("--distthresh"), 0,
	    string("Discards samples shorter than this threshold (in mm - default=0)"),
	    false, requires_argument),
   c_thr(string("-c,--cthr"), 0.2, 
	 string("Curvature threshold - default=0.2"), 
	 false, requires_argument),
  fibthresh(string("--fibthresh"), 0.01, 
	    string("volume fraction before subsidary fibre orientations are considered - default=0.01"), 
	 false, requires_argument),
   steplength(string("--steplength"), 0.5, 
	 string("steplength in mm - default=0.5"), 
	 false, requires_argument),
   loopcheck(string("-l,--loopcheck"), false, 
	 string("perform loopchecks on paths - slower, but allows lower curvature threshold"), 
	 false, no_argument),
   usef(string("-f,--usef"), false, 
	 string("Use anisotropy to constrain tracking"), 
	 false, no_argument),
  randfib(string("--randfib"), 0, 
	 string("Set to 1 to randomly sample initial fibres. Set to 2 to sample in proportion to f"), 
	 false, no_argument),
  fibst(string("--fibst"),1, 
	 string("Force a starting fibre for tracking - default=1, i.e. first fibre orientation"), 
	 false, requires_argument),
  modeuler(string("--modeuler"), false, 
	   string("Use modified euler streamlining"), 
	   false, no_argument),
  rseed(string("--rseed"), 12345,
	string("Random seed"),
	false, requires_argument), 
  seedcountastext(string("--seedcountastext"), false,
		  string("Output seed-to-target counts as a text file (useful when seeding from a mesh)"),
		  false, no_argument), 
  splitmatrix2(string("--splitmatrix2"), false,
		  string("split matrix 2 along seed dimension (in case it is too large)"),
		  false, no_argument), 
   options("probtrackx","probtrackx -s <basename> -m <maskname> -x <seedfile> -o <output> --targetmasks=<textfile>\n probtrackx --help\n")
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
       options.add(waypoints);
       options.add(network);
       options.add(meshfile);
       options.add(lrmask);
       options.add(seedref);
       options.add(logdir); 
       options.add(forcedir); 
       options.add(simpleout);
       options.add(pathdist);
       options.add(s2tout);
       options.add(matrix1out);
       options.add(matrix2out);
       options.add(matrix3out);
       options.add(maskmatrix3);
       options.add(maskmatrixout);
       options.add(outfile);
       options.add(rubbishfile);
       options.add(stopfile);
       options.add(prefdirfile);
       options.add(seeds_to_dti);
       options.add(dti_to_seeds);
       options.add(nparticles);
       options.add(nsteps);
       options.add(distthresh);
       options.add(c_thr);
       options.add(fibthresh);
       options.add(steplength);
       options.add(loopcheck);
       options.add(usef);
       options.add(randfib);
       options.add(fibst);
       options.add(modeuler);
       options.add(rseed);
       options.add(seedcountastext);
       options.add(splitmatrix2);
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







