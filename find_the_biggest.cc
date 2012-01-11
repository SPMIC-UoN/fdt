/*  Copyright (C) 2004 University of Oxford  */

/*  CCOPYRIGHT  */

#include <iostream>
#include <fstream>
#include "meshclass/meshclass.h"
#include "newimage/newimageall.h"
#include "fslvtkio/fslvtkio.h"
#include <vector>
using namespace std;
using namespace NEWIMAGE;
using namespace mesh;
using namespace fslvtkio;

void biggest_from_volumes(vector<string> innames,string oname){
  vector<volume<float> > tmpvec;
  tmpvec.reserve(innames.size());
  volume<float> tmp;
  cout<<"number of inputs "<<innames.size()<<endl;
  cout<<"Indices"<<endl;
  for(unsigned int i=0;i<innames.size();i++){
    cout<<i+1<<" "<<innames[i]<<endl;
    read_volume(tmp,innames[i]);
    tmpvec.push_back(tmp);
  }
  volume<int> output(tmp.xsize(),tmp.ysize(),tmp.zsize());
  copybasicproperties(tmp,output);output=0;

  for(int z=tmp.minz();z<=tmp.maxz();z++){
    for(int y=tmp.miny();y<=tmp.maxy();y++){
      for(int x=tmp.minx();x<=tmp.maxx();x++){
	RowVector bum(innames.size());
	Matrix bugger;
	ColumnVector index;
	for(unsigned int i=0;i<innames.size();i++ ){
	    bum(i+1)=tmpvec[i](x,y,z);
	} 
	bugger=max(bum,index);
	if(bum.MaximumAbsoluteValue()!=0)
	  output(x,y,z)=(int)index.AsScalar();
	else
	  output(x,y,z)=0;
      }  
    }
  }

  output.setDisplayMaximumMinimum(innames.size(),0);
  save_volume(output,oname);

}

void biggest_from_surfaces(vector<string> innames,string oname){
  vector<Mesh> meshes;
  Mesh m;int nv=0,nt=0;

  cout<<"number of inputs "<<innames.size()<<endl;
  cout<<"Indices"<<endl;
  for(unsigned int i=0;i<innames.size();i++){
    cout<<i+1<<" "<<innames[i]<<endl;
    nt++;
    m.load(innames[i]);
    meshes.push_back(m);
    // test size
    if(i==0)nv=m.nvertices();
    else{
      if(m.nvertices()!=nv){
	cerr<<"error: number of vertices incompatible amongst input"<<endl;
	exit(1);
      }
    }
  }
    
  fslvtkIO * fout = new fslvtkIO();
  fout->setMesh(meshes[0]);
  Matrix allvals(nv,nt);
  for(unsigned int i=0;i<meshes.size();i++){
    int cnt=1;
    for(vector<Mpoint*>::iterator it = meshes[i]._points.begin();it!=meshes[i]._points.end();it++){
      allvals(cnt,i+1)=(*it)->get_value();
      cnt++;	
    }
  }
  //OUT(allvals.t());
  // get max
  Matrix values(nv,1);
  values=0;
  for(int i=1;i<=nv;i++){
    if(allvals.Row(i).MaximumAbsoluteValue()==0)continue;
    float vm=0;
    for(int j=1;j<=nt;j++){
      if(j==1)vm=allvals(i,j);
      if(allvals(i,j)>=vm){values(i,1)=j;vm=allvals(i,j);}
    }
  }

  fout->setScalars(values);
  // set vectors so it can be loaded in fslview)
  Matrix vec(nv,3);
  vec=0;
  fout->setVectors(vec);
  
  fout->save(oname);

		      
  
}

int main ( int argc, char **argv ){
  if(argc<3){
    cerr<<" "<<endl;
    cerr<<"usage: find_the_biggest <lots of volumes> output"<<endl;
    cerr<<"output is index in order of inputs"<<endl;
    cerr<<"output will be of the same type as input, i.e. either volume or surface"<<endl;
    cerr<<" "<<endl;
    exit(1);
  }

  vector<string> innames;
  innames.clear();
  for(int i=1;i<=argc-2;i++){
    innames.push_back(argv[i]);
  }
  if(!fsl_imageexists(innames[0])){
    biggest_from_surfaces(innames,argv[argc-1]);
  }
  else{
    biggest_from_volumes(innames,argv[argc-1]);
  }

  return(0);

}









