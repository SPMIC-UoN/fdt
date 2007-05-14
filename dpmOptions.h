#if !defined(dpmOptions_h)
#define dpmOptions_h

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"

using namespace Utilities;

namespace DPM {

class dpmOptions {
 public:
  static dpmOptions& getInstance();
  ~dpmOptions() { delete gopt; }

  Option<bool>   help;
  Option<string> datafile;
  Option<string> logfile;
  Option<string> init_class;
  Option<int>    numclass;
 
  Option<int>    numiter;
  Option<int>    burnin;
  Option<int>    sampleevery;
  void parse_command_line(int argc,char** argv,Log& logger);

 private:
  dpmOptions();  
  const dpmOptions& operator=(dpmOptions&);
  dpmOptions(dpmOptions&);

  OptionParser options;
      
  static dpmOptions* gopt;
  
};


 inline dpmOptions& dpmOptions::getInstance(){
   if(gopt == NULL)
     gopt = new dpmOptions();
   
   return *gopt;
 }

 inline dpmOptions::dpmOptions() :
   help(string("-h,--help"), false,
	string("display this message"),
	false,no_argument),
   datafile(string("-d,--data"), string(""),
	    string("data file"),
	    true,requires_argument),
   logfile(string("-o,--out"), string(""),
	    string("output file"),
	    true, requires_argument),
   init_class(string("--ic,--initclass"), "oneperdata",
	    string("data labelling initialisation"),
	    false, requires_argument),
   numclass(string("-k,--numclass"),-1,
	    string("fix number of classes - default=infinite"),
	    false,requires_argument),
   numiter(string("--ni,--numiter"),2000,
	    string("number of iterations - default=2000"),
	   false,requires_argument),
   burnin(string("--bi,--burnin"),1000,
	  string("number of iterations before sampling - default=1000"),
	  false,requires_argument),
   sampleevery(string("--se,--sampleevery"),1,
	       string("sampling frequency - default=1"),
	       false,requires_argument),
   options("dpm","dpm -d data -o logfile")
   {
     
    
     try {
       options.add(help);
       options.add(datafile);
       options.add(logfile);
       options.add(init_class);
       options.add(numclass);
       options.add(numiter);
       options.add(burnin);
       options.add(sampleevery);
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
