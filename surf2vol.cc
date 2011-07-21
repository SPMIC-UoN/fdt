/*  Copyright (C) 2004 University of Oxford  */

/*  CCOPYRIGHT  */


#include "newimage/newimageio.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "miscmaths/SpMat.h"
#include "meshclass/meshclass.h"
#include "csv.h"

using namespace NEWIMAGE;
using namespace MISCMATHS;
using namespace mesh;


int main(int argc, char** argv){

  if(argc<2){
    cout<<"surf2vol <surf> <refvol> <outvol> <convention>"<<endl;
    exit(1);
  }

  volume<short int> refvol;
  read_volume(refvol,argv[2]);

  CSV csv(refvol);
  csv.set_convention(argv[4]);
  csv.load_rois(argv[1]);
  csv.save_surfvol(argv[3],true);

  return 0;
}
