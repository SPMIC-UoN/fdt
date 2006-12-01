/*  Copyright (C) 2004 University of Oxford  */

/*  CCOPYRIGHT  */

#include "ptx_simple.h"
#include "streamlines.h"

using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;
using namespace Utilities;
using namespace PARTICLE;
using namespace mesh;


void track(){
  probtrackxOptions& opts =probtrackxOptions::getInstance();
  
  ////////////////////////////
  Log& logger = LogSingleton::getInstance();
  if(opts.verbose.value()>1){
    logger.makeDir("particles","particle0",true,false);
  }
  
  volume<int> seedref;
  if(opts.seedref.value()!=""){
    read_volume(seedref,opts.seedref.value());
  }
  else{
    read_volume(seedref,opts.maskfile.value());
  }
  
  Streamliner stline;
  Counter counter(seedref,stline);
  counter.initialise();
  Seedmanager seedmanager(counter);
  
  
  Matrix Seeds = read_ascii_matrix(opts.seedfile.value());
  for(int SN=1; SN<=Seeds.Nrows();SN++){
    float xst=Seeds(SN,1);
    float yst=Seeds(SN,2);
    float zst=Seeds(SN,3);
    seedmanager.run(xst,yst,zst,false,0);
    string add=num2str(Seeds(SN,1))+(string)"_"+num2str(Seeds(SN,2))+(string)"_"+num2str(Seeds(SN,3));
    
    counter.save_pathdist(add);
    counter.reset_prob();
  } //Close Seed number Loop
}
