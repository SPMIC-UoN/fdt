/*  Copyright (C) 2007 University of Oxford  */

/* S. Jbabdi */

/*  CCOPYRIGHT  */

#include <stdio.h>
#include "dpm_gibbs.h"

using namespace DPM;
using namespace Utilities;


int main (int argc, char *argv[]){

  Log& logger = LogSingleton::getInstance();

  dpmOptions& opts = dpmOptions::getInstance();
  opts.parse_command_line(argc,argv,logger);

  // read input files
  Matrix data; // data is nxd
  data=read_ascii_matrix(opts.datafile.value());

  
  // create Gibb's sampler instance
  //cout <<"instanciate DPM_GibbsSampler"<<endl;
  DPM_GibbsSampler gs(data,
		      opts.numiter.value(),opts.burnin.value(),opts.sampleevery.value());

  //cout<<"initialisation"<<endl;
  gs.init();
  //cout<<"running..."<<endl;
  gs.run();


  // save output files
  gs.save();
}

