/*  Copyright (C) 2004 University of Oxford  */

/*  CCOPYRIGHT  */

#include "ptx_twomasks.h"
#include "streamlines.h"

using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;
using namespace Utilities;
using namespace PARTICLE;
using namespace mesh;



void twomasks()
{ 
  probtrackxOptions& opts =probtrackxOptions::getInstance();

  ////////////////////////////////
  //  Log& logger = LogSingleton::getInstance();
  volume<int> seeds,seeds2;
  read_volume(seeds,opts.seedfile.value());
  read_volume(seeds2,opts.mask2.value());


  // correct for non-1 values
  // should be seeds.binarise();
  for(int z=seeds.minz();z<=seeds.maxz();z++)
    for(int y=seeds.miny();y<=seeds.maxy();y++)
      for(int x=seeds.minx();x<=seeds.maxx();x++)
	if(seeds(x,y,z)!=0)
	  seeds(x,y,z)=1;
  for(int z=seeds2.minz();z<=seeds2.maxz();z++)
    for(int y=seeds2.miny();y<=seeds2.maxy();y++)
      for(int x=seeds2.minx();x<=seeds2.maxx();x++)
	if(seeds2(x,y,z)!=0)
	  seeds2(x,y,z)=1;
  




  Streamliner stline;
  Counter counter(seeds,stline);
  counter.initialise();
  Seedmanager seedmanager(counter);
  
  stline.add_waymask(seeds2);  
  for(int z=0;z<seeds.zsize();z++){
    for(int y=0;y<seeds.ysize();y++){
      for(int x=0;x<seeds.xsize();x++){
	if(seeds(x,y,z)>0){
	  seedmanager.run(x,y,z,false,-1); 
	}
      }
    }
  }
  stline.pop_waymasks();
  
  stline.add_waymask(seeds);
  cout<<"added"<<endl;
  for(int z=0;z<seeds2.zsize();z++){
    for(int y=0;y<seeds2.ysize();y++){
      for(int x=0;x<seeds2.xsize();x++){
	if(seeds2(x,y,z)>0){
	  seedmanager.run(x,y,z,false,seeds2(x,y,z)-1); 
	}
      }
    }
  }
  
  
  counter.save();
  
}


