/*  Copyright (C) 2004 University of Oxford  */

/*  CCOPYRIGHT  */

#include "ptx_nmasks.h"
#include "streamlines.h"
#include <time.h>

using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;
using namespace Utilities;
using namespace PARTICLE;
using namespace mesh;



void nmasks()
{ 
  probtrackxOptions& opts =probtrackxOptions::getInstance();
  
  // we need a reference volume for CSV
  // (in case seeds are a list of surfaces)
  volume<short int> refvol;
  if(opts.seedref.value()!="")
    read_volume(refvol,opts.seedref.value());
  else
    read_volume(refvol,opts.maskfile.value());

  cout<<"load seeds"<<endl;  
  // check if seeds are a single volume
  if(fsl_imageexists(opts.seedfile.value())){
    cerr<<"Seed is a volume file - please turn off --network option as it is not needed"<<endl;
    exit(1);
  }

  CSV seeds(refvol);
  seeds.set_convention(opts.meshspace.value());
  seeds.load_rois(opts.seedfile.value());
 
  if(seeds.nVols()==0 && opts.seedref.value()==""){
    cerr<<"Warning: need to set a reference volume when defining a surface-based seed"<<endl;
  }
  // get roiind
  vector<int> roiind(seeds.nRois());int iv=0,is=0;
  for(int i=0;i<seeds.nRois();i++){
    if(seeds.isVol(i)){roiind[iv]=i;iv++;}
    else{roiind[seeds.nVols()+is]=i;is++;}
  }
  
  Streamliner  stline      (seeds);
  Counter      counter     (stline);
  counter.initialise();
  Seedmanager  seedmanager (counter);

  vector<int> keeptotal(seeds.nRois(),0);
  int cnt=-1;

  time_t _time;
  _time=time(NULL);
  // seed from volume-like ROIs
  if(seeds.nVols()>0){
    cout << "Volume seeds" << endl;
    for(int roi=1;roi<=seeds.nVols();roi++){
      cout<<"volume "<<roi-1<<endl;

      cnt++;
      counter.volumeSeeding(roi);
      seedmanager.get_stline().load_waymasks(opts.seedfile.value(),roiind[cnt]);
      if(!opts.waycond.set())
	seedmanager.get_stline().set_waycond("OR");

      for(int z=0;z<seeds.zsize();z++){
	if(opts.verbose.value()>=1)
	  cout <<"sl "<<z<<endl;
	for(int y=0;y<seeds.ysize();y++){
	  for(int x=0;x<seeds.xsize();x++){
	    if(seeds.isInRoi(x,y,z,roi)){
	      counter.updateSeedLocation(seeds.get_loc(roi-1,x,y,z));
	      if(opts.verbose.value()>=1){
		cout <<"run"<<endl;
		cout <<x<<" "<<y<<" "<<z<<endl;
	      }
	      keeptotal[roiind[cnt]] += seedmanager.run((float)x,(float)y,(float)z,
					   false,-1,opts.sampvox.value());
	    }
	  }
	}
      }
      
    }

  }

  // seed from surface-like ROIs
  if(seeds.nSurfs()>0){
    cout << "Surface seeds" << endl;
    ColumnVector pos,dir;
    for(int i=0;i<seeds.nSurfs();i++){
      cout<<"surface "<<i<<endl;
      cnt++;
      counter.surfaceSeeding(i);
      
      seedmanager.get_stline().load_waymasks(opts.seedfile.value(),roiind[cnt]);
      if(!opts.waycond.set())
	seedmanager.get_stline().set_waycond("OR");      

      // inform user if whole surface is used or not
      if( seeds.nActVertices(i) != seeds.nVertices(i) ){
	cout << "  Using a subset of the vertices labelled active (i.e. non zero value)" << endl;
	cout << "   set all values to 0 or non-zero to use entire surface" << endl;
      }

      for(int p=0;p<seeds.get_mesh(i).nvertices();p++){
	// check if active point	
	if(seeds.get_mesh(i).get_point(p)->get_value()==0.0)
	  continue;

	counter.updateSeedLocation(seeds.get_loc(i,p));
	pos=seeds.get_vertex_as_vox(i,p);
	dir=seeds.get_normal_as_vox(i,p);

	if(opts.meshspace.value()=="caret")
	  dir*=-1; // normals in caret point away from the brain
	if(opts.verbose.value()>=1){
	  cout <<"run"<<endl;
	  cout <<pos(1)<<" "<<pos(2)<<" "<<pos(3)<<endl;
	}
	keeptotal[roiind[cnt]] += seedmanager.run(pos(1),pos(2),pos(3),
				     false,-1,false);
	//keeptotal += seedmanager.run(pos(1),pos(2),pos(3),
	//		     false,-1,false,dir);

      }
    }
  }

  cout<<endl<<"time spent tracking: "<<(time(NULL)-_time)<<" seconds"<<endl<<endl;

  // save results
  cout << "save results" << endl;
  counter.save_total(keeptotal);  
  counter.save();

  cout<<"finished"<<endl;
}


