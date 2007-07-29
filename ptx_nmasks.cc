/*  Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT  */

#include "ptx_nmasks.h"
#include "streamlines.h"

using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;
using namespace Utilities;
using namespace PARTICLE;
using namespace mesh;



void nmasks()
{ 
  probtrackxOptions& opts =probtrackxOptions::getInstance();

  ////////////////////////////////
  //  Log& logger = LogSingleton::getInstance();
  vector< volume<int> > seeds;
  volume<int> tmpvol;
  vector<string> masknames;

  if(fsl_imageexists(opts.seedfile.value())){
    cerr << "Seed file must be a text file with a list of seeds in multiple masks mode" << endl;
    exit(0);
  }
  else
    read_masks(masknames,opts.seedfile.value());

  for(unsigned int m=0;m<masknames.size();m++){
    read_volume(tmpvol,masknames[m]);
    tmpvol=NEWIMAGE::abs(tmpvol);
    tmpvol.binarise(0,tmpvol.max()+1,exclusive);
    seeds.push_back(tmpvol);
  }

  Streamliner stline;
  Counter counter(seeds[0],stline);
  counter.initialise();
  Seedmanager seedmanager(counter);

  for(unsigned int i=0;i<seeds.size();i++){
    stline.clear_waymasks();
    tmpvol=0;
    for(unsigned int j=0;j<seeds.size();j++){
      if(j!=i)
	tmpvol += seeds[j];
      //	stline.add_waymask(seeds[j]);
    }
    stline.add_waymask(tmpvol);

    for(int z=0;z<seeds[i].zsize();z++){
      for(int y=0;y<seeds[i].ysize();y++){
	for(int x=0;x<seeds[i].xsize();x++){
	  if(seeds[i](x,y,z)!=0){
	    seedmanager.run(x,y,z,false,-1); 
	  }
	}
      }
    } 
  }
  
  counter.save();
  
  cout<<"finished"<<endl;
}


