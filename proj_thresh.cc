/*  Copyright (C) 2004 University of Oxford  */

/*  CCOPYRIGHT  */

#include <iostream>
#include <fstream>
#include "newimage/newimageall.h"
#include <vector>
using namespace std;
using namespace NEWIMAGE;

int main ( int argc, char **argv ){
  if(argc<3){
    cerr<<"usage:proj_thresh <lots of volumes> threshold"<<endl;
    exit(0);
  }
  vector<volume<float> > tmpvec;
  //  vector<volume<float> > outvec;
  float thresh=atof(argv[argc-1]);
  tmpvec.reserve(argc-2);
  //  outvec.reserve(argc-2);
  volume<float> tmp;
  cout<<"number of inputs "<<argc-2<<endl;
  for(int i=1;i<=argc-2;i++){
    cout<<i<<" "<<argv[i]<<endl;
    read_volume(tmp,argv[i]);
    tmpvec.push_back(tmp);
  }
  cerr<<"threshold "<<thresh<<endl;

  volume<float> total;
  volume<float> total_thresh;
  total=tmp*0;

  for(unsigned int i=0;i<tmpvec.size();i++){
    total+=tmpvec[i];
  }
  
  total_thresh=binarise(total,thresh);
  save_volume(total,"total");
  for(unsigned int i=0;i<tmpvec.size();i++){
    tmp=divide(tmpvec[i],total,total_thresh);
    string outname =argv[i+1];
    make_basename(outname);
    string thrname="_thr_"+num2str(thresh);
    save_volume(tmp,outname+"_proj_seg"+thrname);
  }
  
}










