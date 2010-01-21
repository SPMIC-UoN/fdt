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
  
  volume<float> seedref;
  if(opts.seedref.value()!=""){
    read_volume(seedref,opts.seedref.value());
  }
  else{
    read_volume(seedref,opts.maskfile.value());
  }

  Streamliner stline(seedref);
  Counter counter(seedref,stline);
  counter.initialise();
  Seedmanager seedmanager(counter);
    
  Matrix Seeds = read_ascii_matrix(opts.seedfile.value());
  if(Seeds.Ncols()!=3 && Seeds.Nrows()==3)
	Seeds=Seeds.t();
  
  // convert coordinates from nifti (external) to newimage (internal)
  //   conventions - Note: for radiological files this should do nothing
  for (int n=1; n<=Seeds.Nrows(); n++) {
    ColumnVector v(4);
    v << Seeds(n,1) << Seeds(n,2) << Seeds(n,3) << 1.0;
    v = seedref.niftivox2newimagevox_mat() * v;
    Seeds(n,1) = v(1);  Seeds(n,2) = v(2);  Seeds(n,3) = v(3);
  }

  int keeptot=0;
  for(int SN=1; SN<=Seeds.Nrows();SN++){
    float xst=Seeds(SN,1);
    float yst=Seeds(SN,2);
    float zst=Seeds(SN,3);
    keeptot += seedmanager.run(xst,yst,zst,false,0);
    string add="_"+num2str(Seeds(SN,1))+(string)"_"+num2str(Seeds(SN,2))+(string)"_"+num2str(Seeds(SN,3));
    
    counter.save_pathdist(add);
    counter.count_seed();

    counter.reset_prob();
  } //Close Seed number Loop
  
  counter.save();

  cout<<"finished"<<endl;
}
