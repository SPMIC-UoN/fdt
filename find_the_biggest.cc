/*  Copyright (C) 2004 University of Oxford  */

/*  CCOPYRIGHT  */

#include <iostream>
#include <fstream>
#include "newimage/newimageall.h"
#include <vector>
using namespace std;
using namespace NEWIMAGE;

void biggest_from_volumes(vector<string> innames,string oname){
  vector<volume<float> > tmpvec;
  tmpvec.reserve(innames.size());
  volume<float> tmp;
  cout<<"number of inputs "<<innames.size()<<endl;
  cout<<"Indices"<<endl;
  for(unsigned int i=0;i<innames.size();i++){
    cout<<i<<" "<<innames[i]<<endl;
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
	bool flag=true;
	if(index.AsScalar()==1){
	  // Check to see whether they're all zero.
	  flag=false;
	  for(unsigned int i=1;i<=innames.size();i++ ){
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
  save_volume(output,oname);

}

ReturnMatrix read_label(const string& labelfile){
  Matrix L;
  ifstream fs(labelfile.c_str());
  if (!fs) { 
    cerr << "Could not open label file " << labelfile << endl;
    L.Release();
    return L;
  }
  string cline;
  // skip header
  cline = skip_alpha(fs);
  // read number of vertices
  string ss="";
  fs >> ss;
  float nrows = atof(ss.c_str());
  L.ReSize(int(nrows),5);
  for(int r=1;r<=nrows;r++)
    for(int c=1;c<=5;c++){
      if(!fs.eof()){
	fs >> ss;
	while ( !isNumber(ss) && !fs.eof() ) {
	  fs >> ss;
	}
	L(r,c) = atof(ss.c_str());
      }
    }
  
  L.Release();
  return L;
}
void write_label(const Matrix& L,const string& filename){
  ofstream fs(filename.c_str());
  if (!fs) { 
    cerr << "Could not open file " << filename << " for writing" << endl;
    exit(1);
  }
  fs << "##!ascii label , written from FSL414" << endl;
  fs << L.Nrows() << endl;

#ifdef PPC64	
  int n=0;
#endif
  for (int i=1; i<=L.Nrows(); i++) {
    for (int j=1; j<=L.Ncols(); j++) {
      fs << L(i,j) << "  ";
#ifdef PPC64	
      if ((n++ % 50) == 0) fs.flush();
#endif
    }
    fs << endl;
  }
  
 
  fs.close();
}

void biggest_from_matrix(vector<string> innames,string oname){
  if(innames.size()!=2){
    cerr<<"usage: find_the_biggest <Matrix> <Label> <Output>"<<endl;
    exit(1);
  }

  cout << innames[0] << endl;
  Matrix M = read_ascii_matrix(innames[0]);
  ColumnVector maxMrow(M.Nrows()),coordMax(M.Nrows());
  maxMrow = sum(abs(M),2);

  cout << "number of vertices: " << M.Nrows() << endl;
  cout << "number of targets:  " << M.Ncols() << endl;

  cout << endl << "read label file" << endl;
  Matrix L = read_label(innames[1]);
  
  cout << "number of vertices: " << L.Nrows() << endl;
  cout << "number of columns:  " << L.Ncols() << endl;

  if(M.Nrows() != L.Nrows()){
    cerr<<"Error: matrix and label file do not have the same number of entries"<<endl;
    exit(1);
  }


  vector< vector<int> > Clusters(M.Ncols());
  float val;
  int cmax;
  for(int i=1;i<=M.Nrows();i++){
    val=(M.Row(i)).MaximumAbsoluteValue1(cmax);
    if(val!=0)
      Clusters[cmax-1].push_back(int(i));
  }
  // now store these into files
  for(unsigned int i=0;i<Clusters.size();i++){
    if(Clusters[i].size()>0){
      Matrix C(Clusters[i].size(),5);
      for(unsigned int j=0;j<Clusters[i].size();j++){
	C.Row(j+1) = L.Row(Clusters[i][j]); 
      }
      write_label(C,oname+"_"+num2str(i+1)+".label");
    }
  }
  
}

int main ( int argc, char **argv ){
  if(argc<3){
    cerr<<" "<<endl;
    cerr<<"usage: find_the_biggest <lots of volumes> output"<<endl;
    cerr<<"output is index in order of inputs"<<endl;
    cerr<<"Or: find_the_biggest <singleMatrixFile> <labelfile> output"<<endl;
    cerr<<" "<<endl;
    exit(1);
  }

  vector<string> innames;
  innames.clear();
  for(int i=1;i<=argc-2;i++){
    innames.push_back(argv[i]);
  }
  if(!fsl_imageexists(innames[0])){
    biggest_from_matrix(innames,argv[argc-1]);
  }
  else{
    biggest_from_volumes(innames,argv[argc-1]);
  }

  return(0);

}









