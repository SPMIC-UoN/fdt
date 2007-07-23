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
    cerr<<"usage: find_the_biggest <lots of volumes> output"<<endl;
    cerr<<"output is index in order of inputs"<<endl;
    exit(1);
  }
  vector<volume<float> > tmpvec;
  tmpvec.reserve(argc-2);
  volume<float> tmp;
  cout<<"number of inputs "<<argc-2<<endl;
  cout<<"Indices"<<endl;
  for(int i=1;i<=argc-2;i++){
    cout<<i<<" "<<argv[i]<<endl;
    read_volume(tmp,argv[i]);
    tmpvec.push_back(tmp);
  }
  volume<int> output(tmp.xsize(),tmp.ysize(),tmp.zsize());
  copybasicproperties(tmp,output);output=0;

  for(int z=tmp.minz();z<=tmp.maxz();z++){
    for(int y=tmp.miny();y<=tmp.maxy();y++){
      for(int x=tmp.minx();x<=tmp.maxx();x++){
	RowVector bum(argc-2);
	Matrix bugger;
	ColumnVector index;
	for( int i=0;i<argc-2;i++ ){
	    bum(i+1)=tmpvec[i](x,y,z);
	} 
	bugger=max(bum,index);
	bool flag=true;
	if(index.AsScalar()==1){
	  // Check to see whether they're all zero.
	  flag=false;
	  for(int i=1;i<=argc-2;i++ ){
	    if(bum(i)!=0){
	      flag=true;break;
	    }
	      
	  }
	}
	if(flag)
	  output(x,y,z)=(int)index.AsScalar();
	else
	  output(x,y,z)=0;
      }  
    }
  }
  save_volume(output,argv[argc-1]);
}









