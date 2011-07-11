/*  Copyright (C) 2004 University of Oxford  */

/*  CCOPYRIGHT  */

#include <iostream>
#include <fstream>
#include "newimage/newimageall.h"
#include "utils/log.h"
#include "probtrackx.h"
#include "utils/tracer_plus.h"

using namespace std;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace TRACT;


int main ( int argc, char **argv ){
  //Tracer_Plus::settimingon();

  probtrackxOptions& opts =probtrackxOptions::getInstance();
  Log& logger = LogSingleton::getInstance();
  opts.parse_command_line(argc,argv,logger);
  srand(opts.rseed.value());
  
  if(opts.verbose.value()>0){
    opts.status();
  }
  if(opts.simple.value()){
    if( opts.matrix1out.value() || opts.matrix3out.value()){
      cerr<<"Error: cannot use matrix1 and matrix3 in simple mode"<<endl;
      exit(1);
    }
    cout<<"Running in simple mode"<<endl;
    track();
  }
  else{
    if(!opts.network.value()){
      cout<<"Running in seedmask mode"<<endl;
      seedmask();
    }
    else{
      cout<<"Running in network mode"<<endl;
      nmasks();
    }
  }

  //Tracer_Plus::dump_times(logger.getDir());

  return 0;
}















