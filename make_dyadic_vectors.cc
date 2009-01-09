/*  Copyright (C) 2004 University of Oxford  */

/*  CCOPYRIGHT  */

#include <iostream>
#include <fstream>
#include "newimage/newimageall.h"
#include <vector>
using namespace std;
using namespace NEWIMAGE;

int main ( int argc, char **argv ){
  if(argc<4){
    cerr<<"usage: make_dyadic_vectors <theta_vol> <phi_vol> [mask] <output>"<<endl;
    cerr<<"[mask] is optional"<<endl;
    exit(1);
  }

  volume4D<float> ths,phs;
  read_volume4D(ths,argv[1]);
  read_volume4D(phs,argv[2]);
  volume<float> mask;
  string oname;
  if(argc==5){
    oname=argv[4];
    read_volume(mask,argv[3]);
  }
  else{
    mask=ths[0]*0+1;
    oname=argv[3];
  }
  volume4D<float> dyadic_vecs(ths.xsize(),ths.ysize(),ths.zsize(),3);
  dyadic_vecs=0;
  copybasicproperties(ths[0],dyadic_vecs[0]);
  SymmetricMatrix dyad(3);dyad=0;
  ColumnVector dir(3);
  
  DiagonalMatrix dyad_D; //eigenvalues
  Matrix dyad_V; //eigenvectors

  for(int k=ths.minz();k<=ths.maxz();k++){
    for(int j=ths.miny();j<=ths.maxy();j++){
      for(int i=ths.minx();i<=ths.maxx();i++){
	if(mask(i,j,k)>0){
	  dyad=0;
	  for(int s=ths.mint();s<=ths.maxt();s++){		   
	  
	    float th=ths(i,j,k,s);
	    float ph=phs(i,j,k,s);
	    dir(1)=sin(th)*cos(ph);
	    dir(2)=sin(th)*sin(ph);
	    dir(3)=cos(th);
	    dyad << dyad+dir*dir.t();
	  }
	  
	  EigenValues(dyad,dyad_D,dyad_V);
	  int maxeig;
	  if(dyad_D(1)>dyad_D(2)){
	    if(dyad_D(1)>dyad_D(3)) maxeig=1;
	    else maxeig=3;
	  }
	  else{
	    if(dyad_D(2)>dyad_D(3)) maxeig=2;
	    else maxeig=3;
	  }
	  dyadic_vecs(i,j,k,0)=dyad_V(1,maxeig);
	  dyadic_vecs(i,j,k,1)=dyad_V(2,maxeig);
	  dyadic_vecs(i,j,k,2)=dyad_V(3,maxeig);
	}
	else{
	  dyadic_vecs(i,j,k,0)=0;
	  dyadic_vecs(i,j,k,1)=0;
	  dyadic_vecs(i,j,k,2)=0;
	}
      }
    }
  }
  
  save_volume4D(dyadic_vecs,oname);
}









