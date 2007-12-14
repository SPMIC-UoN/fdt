/*  Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT  */

#include "ptx_meshmask.h"
#include "streamlines.h"
#include "shapeModel/shapeModel.h"
using namespace shapemodel;
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

  Matrix mm_to_vox(4,4);
  mm_to_vox << -1 << 0 << 0 <<  (seeds.xsize()-1)/2
	    <<  0 << 0 << -1 << (seeds.zsize()-1)/2
	    <<  0 << 1 << 0 <<  (seeds.ysize()-1)/2
	    <<  0 << 0 << 0 << 1;

// shapeModel * model1 = new shapeModel;
// model1->load_vtk(opts.meshfile.value(),1);

//  model1->meshReg(&mseeds, mm_to_vox);
//  model1->setShapeMesh(0,mseeds);
//  model1->setImageParameters(seeds.xsize(),seeds.ysize(),seeds.zsize(),seeds.xdim(),seeds.ydim(), seeds.zdim());

//  Mesh m = model1->getTranslatedMesh(0);

//   for(vector<Mpoint*>::iterator i = m._points.begin();i!=m._points.end();i++){
//     cout<<(*i)->get_coord().X<<" "<<(*i)->get_coord().Y<<" "<<(*i)->get_coord().Z<<endl; 

//   }
  ColumnVector fs_coord_mm(4),xyz_vox,seeddim(3);
  seeddim << seeds.xdim() << seeds.ydim() << seeds.zdim();
  ColumnVector dir(3);
  int keeptotal=0;

  for(vector<Mpoint*>::iterator i = mseeds._points.begin();i!=mseeds._points.end();i++){
    //if((*i)->get_value() > 0){

      fs_coord_mm<<(*i)->get_coord().X<<(*i)->get_coord().Y<<(*i)->get_coord().Z << 1.0; 
      //      xyz_vox = seeds.qform_mat().i()*fs_coord_mm
      xyz_vox = mm_to_vox*fs_coord_mm;

      float x=xyz_vox(1);float y=xyz_vox(2);float z=xyz_vox(3);
      Pt newPt(x,y,z);
                (*i)->_update_coord = newPt;

      seeds(round(x),round(y),round(z)) = 1;
    
      cout <<"run"<<endl;
      dir << (*i)->local_normal().X << (*i)->local_normal().Y << (*i)->local_normal().Z;
      keeptotal += seedmanager.run(round(x),round(y),round(z),true,-1,dir); 
	
      //}
  }
  mseeds.update();
  mseeds.save("test.vtk",3);
  //return;

//   for(int z=0;z<seeds.zsize();z++){
//     cout <<"sl "<<z<<endl;
//     for(int y=0;y<seeds.ysize();y++){
//       for(int x=0;x<seeds.xsize();x++){
// 	if(seeds(x,y,z)>0){
// 	  cout <<"run"<<endl;
// 	  dir << (*i)->local_normal().X << (*i)->local_normal().Y << (*i)->local_normal().Z;
// 	  keeptotal += seedmanager.run(x,y,z,true,-1,dir); 
// 	}
//       }
//     }
//   }

  counter.save_total(keeptotal);  
  counter.save();

  cout<<"finished"<<endl;
}


