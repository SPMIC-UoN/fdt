/*  Copyright (C) 2004 University of Oxford  */

/*  CCOPYRIGHT  */

#include "ptx_twomasks.h"
#include "streamlines.h"

using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;
using namespace Utilities;
using namespace PARTICLE;
using namespace TRACTVOLS;
using namespace mesh;



void twomasks()
{ 
  probtrackOptions& opts =probtrackOptions::getInstance();

  ////////////////////////////////
  //  Log& logger = LogSingleton::getInstance();
  volume<int> seeds,seeds2;
  read_volume(seeds,opts.seedfile.value());
  read_volume(seeds2,opts.mask2.value());
  Streamliner stline;
  Counter counter(seeds,stline);
  counter.initialise();
  Seedmanager seedmanager(counter);
  
  stline.add_waymask(seeds2);  
  for(int z=0;z<seeds.zsize();z++){
    for(int y=0;y<seeds.ysize();y++){
      for(int x=0;x<seeds.xsize();x++){
	if(seeds(x,y,z)>0){
	  seedmanager.run(x,y,z); 
	}
      }
    }
  }
  stline.pop_waymasks();
  
  
  stline.add_waymask(seeds);
  for(int z=0;z<seeds2.zsize();z++){
    for(int y=0;y<seeds2.ysize();y++){
      for(int x=0;x<seeds2.xsize();x++){
	if(seeds2(x,y,z)>0){
	  seedmanager.run(x,y,z); 
	}
      }
    }
  }
  

  counter.save();
  
}


