/*  Copyright (C) 2004 University of Oxford  */

/*  CCOPYRIGHT  */

#include "ptx_seedmask.h"
#include "streamlines.h"

using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;
using namespace Utilities;
using namespace PARTICLE;
using namespace TRACTVOLS;
using namespace mesh;



void seedmask()
{ 
  probtrackxOptions& opts =probtrackxOptions::getInstance();

  ////////////////////////////////
  //  Log& logger = LogSingleton::getInstance();
  volume<int> seeds;
  read_volume(seeds,opts.seedfile.value());
  Streamliner stline;
  Counter counter(seeds,stline);
  counter.initialise();
  Seedmanager seedmanager(counter);
  
  for(int z=0;z<seeds.zsize();z++){
    for(int y=0;y<seeds.ysize();y++){
      for(int x=0;x<seeds.xsize();x++){
	if(seeds(x,y,z)>0){
	  seedmanager.run(x,y,z); 
	}
      }
    }
  }
  
  counter.save();
  
}


