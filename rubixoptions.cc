/*  rubixoptions.cc

    Stam Sotiropoulos, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */

#define WANT_STREAM
#define WANT_MATH

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "rubixoptions.h"
#include "utils/log.h"
#include "utils/tracer_plus.h"

using namespace Utilities;

namespace RUBIX {
  
  rubixOptions* rubixOptions::gopt = NULL;
  
  void rubixOptions::parse_command_line(int argc, char** argv, Log& logger){
    Tracer_Plus("RubiXOptions::parse_command_line");

    // do once to establish log directory name
    for(int a = options.parse_command_line(argc, argv); a < argc; a++);
  
    if(help.value() || ! options.check_compulsory_arguments()){
      options.usage();
      //throw Exception("Not all of the compulsory arguments have been provided");
      exit(2);
    }
    else{
      // setup logger directory
      if(forcedir.value())
	logger.setthenmakeDir(logdir.value());
      else
	logger.makeDir(logdir.value());
    
      cout << "Log directory is: " << logger.getDir() << endl;
    
      // do again so that options are logged
      for(int a = 0; a < argc; a++)
	logger.str() << argv[a] << " ";
      logger.str() << endl << "---------------------------------------------" << endl << endl;
    }      
  }

}
