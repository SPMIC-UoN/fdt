#include <iostream>
#include <fstream>
#include "newimage/newimageall.h"
#include <vector>
using namespace std;
using namespace NEWIMAGE;

int main ( int argc, char **argv ){
  if(argc<4){
    cerr<<"usage: make_dyadic_vectors <theta_vol> <phi_vol> <output>"<<endl;
    exit(0);
  }
  volume4D<float> ths,phs;
  volumeinfo tempinfo;
  read_volume4D(ths,argv[1],tempinfo);
  read_volume4D(phs,argv[2]);
  volume4D<float> dyadic_vecs(ths.xsize(),ths.ysize(),ths.zsize(),3);
  dyadic_vecs=0;
  copybasicproperties(ths,dyadic_vecs);
  SymmetricMatrix dyad(3);dyad=0;
  ColumnVector dir(3);
  
  DiagonalMatrix dyad_D; //eigenvalues
  Matrix dyad_V; //eigenvectors

  for(int k=ths.minz();k<=ths.maxz();k++){
    for(int j=ths.miny();j<=ths.maxy();j++){
      for(int i=ths.minx();i<=ths.maxx();i++){
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
      
    }
  }
  
  save_volume4D(dyadic_vecs,argv[3],tempinfo);
}









