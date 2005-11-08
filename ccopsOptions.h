/*  BpmOptions.h

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(ccopsOptions_h)
#define ccopsOptions_h

#include <string>
#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
//#include "newmatall.h"
using namespace Utilities;

namespace CCOPS {

class ccopsOptions {
 public:
  static ccopsOptions& getInstance();
  ~ccopsOptions() { delete gopt; }
  
  Option<bool> help;
  Option<string> inmatrix;
  Option<string> basename;
  Option<string> excl_mask;
  Option<bool> reord1;
  Option<bool> reord2;
  Option<int> bin;
  Option<float> power;
  bool parse_command_line(int argc, char** argv);
  
 private:
  ccopsOptions();  
  const ccopsOptions& operator=(ccopsOptions&);
  ccopsOptions(ccopsOptions&);

  OptionParser options; 
      
  static ccopsOptions* gopt;
  
};

 inline ccopsOptions& ccopsOptions::getInstance(){
   if(gopt == NULL)
     gopt = new ccopsOptions();
   
   return *gopt;
 }

 inline ccopsOptions::ccopsOptions() :
   help(string("-h,--help"), false,
	string("display this message"),
	false, no_argument),
   inmatrix(string("-i,--in"), string(""),
	       string("input matrix"),
	       true, requires_argument),  
   basename(string("-b,--basename"), string(""),
	       string("Output basename"),
	       true, requires_argument),
   excl_mask(string("-x"), string(""),
	     string("exclusion mask"),
	     false, requires_argument),  
   reord1(string("--r1"), bool(false),
	     string("do seedspace reordering (default no)"),
	     false, no_argument), 
   reord2(string("--r2"), bool(false),
	     string("do tractspace reordering (default no)"),
	     false, no_argument), 
   bin(string("--bin"), 0, 
	 string("binarise at (default 0 - no binarisation)"), 
	 false, requires_argument),
   power(string("-p,--power"), 1, 
	 string("power to raise the correlation matrix to (default 1)"), 
	 false, requires_argument),
   options("ccops","")
   {
     
    
     try {
       options.add(help);
       options.add(inmatrix);
       options.add(excl_mask);
       options.add(reord1);
       options.add(reord2);
       options.add(bin);
       options.add(power);
       
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





