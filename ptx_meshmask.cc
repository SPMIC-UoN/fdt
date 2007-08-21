/*  Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT  */

#include "ptx_meshmask.h"
#include "streamlines.h"

using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;
using namespace Utilities;
using namespace PARTICLE;
using namespace mesh;



void meshmask()
{ 
  probtrackxOptions& opts =probtrackxOptions::getInstance();


  // load seed mesh
  Mesh mseeds;
  mseeds.load(opts.meshfile.value());
  mseeds.load_fs_label(opts.seedfile.value());

  // internally create seed mask in voxel space
  volume<int> seeds;
  if(opts.seedref.value()!="")
    read_volume(seeds,opts.seedref.value());
  else
    read_volume(seeds,opts.maskfile.value());
  seeds=0;



  ////////////////////////////////
  //  Log& logger = LogSingleton::getInstance();
  Streamliner stline;
  Counter counter(seeds,stline);
  counter.initialise();
  Seedmanager seedmanager(counter);

  Matrix Seeds_to_DTI;
  if(opts.seeds_to_dti.value()!=""){
    read_ascii_matrix(Seeds_to_DTI,opts.seeds_to_dti.value());
  }
  else{
    Seeds_to_DTI=Identity(4);
  }
  ColumnVector mni_origin(3),fs_coord_mm(3),xyz_dti;
  mni_origin << 92 << 128 << 37;
  for(vector<Mpoint*>::iterator i = mseeds._points.begin();i!=mseeds._points.end();i++){
    if((*i)->get_value() > 0){
      fs_coord_mm<<(*i)->get_coord().X<<(*i)->get_coord().Y<<(*i)->get_coord().Z; 
      xyz_dti=mni_to_imgvox(fs_coord_mm,mni_origin,Seeds_to_DTI,counter.get_streamline().get_vols().dimensions());
      float x=xyz_dti(1);float y=xyz_dti(2);float z=xyz_dti(3);

      seeds(round(x),round(y),round(z)) = 1;
    }
  }



  int keeptotal=0;
  for(int z=0;z<seeds.zsize();z++){
    cout <<"sl "<<z<<endl;
    for(int y=0;y<seeds.ysize();y++){
      for(int x=0;x<seeds.xsize();x++){
	if(seeds(x,y,z)>0){
	  cout <<"run"<<endl;
	  keeptotal += seedmanager.run(x,y,z); 
	}
      }
    }
  }

  counter.save_total(keeptotal);  
  counter.save();

  cout<<"finished"<<endl;
}


