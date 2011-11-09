/*  Copyright (C) 2010 University of Oxford  */
/* Stam Sotiropoulos    */
/*  CCOPYRIGHT  */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include "newimage/newimageall.h"

using namespace std;
using namespace NEWIMAGE;
using namespace NEWMAT;

//Computes the model 1 predicted signal from the mode of the bedpostx samples
//An f0 compartment could also be included in the model   
int main ( int argc, char *argv[]){
  if(argc<6 || argc>7){
    cerr<<" "<<endl;
    cerr<<"usage: xfibres_pred data_dir samples_dir fibres_num f0(0 or 1) [brain_mask] ouput"<<endl;
    cerr<<" "<<endl;
    exit(1);
  }
  
  Matrix bvecs,bvals;   //Read Input 
  volume<float> mask; int num_fibres;
  volume<float> d,S0,f0, temp;
  volume4D<float> temp4D;
  vector< volume<float> > f;  
  vector< volume4D<float> > dyads;  
  
  string dir_name=argv[1], temp_name;
  temp_name=dir_name+"/bvals";
  bvals=read_ascii_matrix(temp_name);
  temp_name=dir_name+"/bvecs";
  bvecs=read_ascii_matrix(temp_name);
  temp_name=argv[3];
  num_fibres=atoi(temp_name.c_str());
  
  if(bvecs.Nrows()>3) bvecs=bvecs.t();   //Make sure the bvecs entries are normalized unit vectors
  if(bvals.Nrows()>1) bvals=bvals.t();
  for(int i=1;i<=bvecs.Ncols();i++){
    float tmpsum=sqrt(bvecs(1,i)*bvecs(1,i)+bvecs(2,i)*bvecs(2,i)+bvecs(3,i)*bvecs(3,i));
    if(tmpsum!=0){
      bvecs(1,i)=bvecs(1,i)/tmpsum;
      bvecs(2,i)=bvecs(2,i)/tmpsum;
      bvecs(3,i)=bvecs(3,i)/tmpsum;
    }  
  }
  
  dir_name=argv[2]; 
  temp_name=dir_name+"/mean_dsamples";   //Read bedpostx results 
  if (!fsl_imageexists(temp_name)){
    cout<<"No mean_dsamples file exists!"<<endl;
    exit(1); 
  }
  else read_volume(d,temp_name);
 
  temp_name=dir_name+"/mean_S0samples";   //In case no S0samples is saved, use the mean b=0 from the data
  if (!fsl_imageexists(temp_name)){
    string dir_name2=argv[1];
    temp_name=dir_name2+"/data"; 
    if (!fsl_imageexists(temp_name)){
      cerr<<"No mean_S0samples file exists or data file exists!"<<endl;
      exit(1);
    }
    else{
      cout<<"No mean_SOsamples found, getting the average of b=0's from data..."<<endl;
      volume4D<float> data;
      read_volume4D(data,temp_name);
      S0.reinitialize(data.xsize(),data.ysize(),data.zsize());
      S0=0;
      int b0count=0;
      for (int l=1; l<=bvals.Ncols(); l++)
	if (bvals(1,l)==0){
	  S0+=data[l-1];
	  b0count++;
	}
      S0/=b0count;
    }
  }
  else read_volume(S0,temp_name);
  
  for (int n=0; n<num_fibres; n++){   //Read dyads
    temp_name=dir_name+"/dyads"+num2str(n+1);
    if (!fsl_imageexists(temp_name)){
      cerr<<"No dyads"<<n+1<<" file exists!"<<endl;
      exit(1); }
    else{
      read_volume4D(temp4D,temp_name);
      dyads.push_back(temp4D);
    }
    temp_name=dir_name+"/mean_f"+num2str(n+1)+"samples";
    if (!fsl_imageexists(temp_name)){
      cerr<<"No mean_f"<<n+1<<"samples file exists!"<<endl;
      exit(1); }
    else{
      read_volume(temp,temp_name);
      f.push_back(temp);
    }
  }
  string argv4=argv[4];

  if (argv4=="1"){                     //Read f0
    temp_name=dir_name+"/mean_f0samples";    
    if (!fsl_imageexists(temp_name)){
      cout<<"No mean_f0samples file exists!"<<endl;
      exit(1); 
    }
    else read_volume(f0,temp_name);   
  }

  if (argc==7)   //Read Brain Mask
    read_volume(mask,argv[5]);
  else{ 
    mask.reinitialize(d.xsize(),d.ysize(),d.zsize());
    mask=1;
  }
  
  volume4D<float> output;
  copybasicproperties(d,output);
  output.reinitialize(d.xsize(),d.ysize(),d.zsize(),bvals.Ncols());
  output=0;

  for(int z=d.minz();z<=d.maxz();z++){   //Compute predicted signal for each voxel
    for(int y=d.miny();y<=d.maxy();y++){
      for(int x=d.minx();x<=d.maxx();x++){
	if (mask(x,y,z)!=0){
	  for (int l=0; l<bvals.Ncols(); l++){
	    float sumf=0;
	    for (int n=0; n<num_fibres; n++){
	      sumf+=f[n](x,y,z);
	      float angp=dyads[n](x,y,z,0)*bvecs(1,l+1)+dyads[n](x,y,z,1)*bvecs(2,l+1)+dyads[n](x,y,z,2)*bvecs(3,l+1);
	      output(x,y,z,l)+=f[n](x,y,z)*exp(-bvals(1,l+1)*d(x,y,z)*angp*angp);
	    }
	    if (argv4=="1")
	      output(x,y,z,l)+=(f0(x,y,z)+(1-sumf-f0(x,y,z))*exp(-bvals(1,l+1)*d(x,y,z)));
	    else  
	      output(x,y,z,l)+=(1-sumf)*exp(-bvals(1,l+1)*d(x,y,z));
	    output(x,y,z,l)*=S0(x,y,z); 
	  }
	}
      }
    }
    cout<<z+1<<" slices processed"<<endl;
  }
  save_volume4D(output,argv[argc-1]); 
  return 0;
}









