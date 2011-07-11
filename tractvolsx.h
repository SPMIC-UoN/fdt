/*  tractvolsx.h

    Tim Behrens, Saad Jbabdi, FMRIB Image Analysis Group

    Copyright (C) 2004 University of Oxford  */

/*  CCOPYRIGHT  */

#ifndef __TRACTVOLSX_H_
#define __TRACTVOLSX_H_

/////////////////////////////////////////////////////////
//         Class TractVolsx                             //
/////////////////////////////////////////////////////////

#include "newimage/newimageall.h"
#include <iostream>
#include "stdlib.h"
#include "probtrackxOptions.h"
using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;

namespace TRACTVOLSX{
  class Tractvolsx
    {
    private:
      probtrackxOptions& opts;
      Log&               logger;

      vector<Matrix> thsamples;
      vector<Matrix> phsamples;
      vector<Matrix> fsamples;

      volume<int>    lut_vol2mat;

      int            nfibres;
      int            nsamples;

      bool           init_sample;
      int            fibst;
      bool           usef;
      
    public:
      //constructors::
      Tractvolsx(const bool& usefin=false):opts(probtrackxOptions::getInstance()),
					   logger(LogSingleton::getInstance()),
					   init_sample(true),fibst(0),usef(usefin){}
      ~Tractvolsx(){}
      int get_nfibres()const{return nfibres;}
      int get_nsamples()const{return nsamples;}
      
      void reset(const int& fibst_in){
	init_sample=true;
	fibst=fibst_in;
      }

      int sample_fibre(int col,int samp,const int& mode=2){
	if(mode==0){
	  return 0;
	}
	if(mode==3){//sample all
	  return int(round((float)(nfibres-1)*(float)rand()/float(RAND_MAX)));
	}
	else{
	  if(mode==1){//sample all>thresh
	    vector<int> fibvec;
	    for(int fib=0;fib<nfibres;fib++){	    
	      float ft=fsamples[fib](samp,col);
	      if(ft>opts.fibthresh.value()){
		fibvec.push_back(fib);
	      }
	    }
	    if(fibvec.size()==0){
	      return 0;
	    }
	    else{
	      float rtmp=(float)rand()/float(RAND_MAX) * float(fibvec.size()-1);
	      return (fibvec[ (int)round(rtmp) ]);
	    }
	  }
	  else if(mode==2){//sample all>thresh in proportion of f (default)
	    float fsumtmp=0;
	    for(int fib=0;fib<nfibres;fib++){	    
	      float ft=fsamples[fib](samp,col);
	      if(ft>opts.fibthresh.value()){
		fsumtmp+=ft;  //count total weight of f in this voxel. 
	      }
	    } 
	    if(fsumtmp==0){
	      return(0);
	    }
	    else{
	      float ft,fsumtmp2=0;
	      float rtmp=fsumtmp * (float)rand()/float(RAND_MAX);	      
	      for(int fib=0;fib<nfibres;fib++){
		ft=fsamples[fib](samp,col);
		if(ft>opts.fibthresh.value())
		  fsumtmp2 += ft;
		if(rtmp<=fsumtmp2){
		  return(fib); 
		}
	      }
	    }
	  }
	  else{
	    cerr<<"TRACTVOLSX::sample_fibre:Error - unknown mode = "<<mode<<endl;
	    exit(1);
	  }
	}
	return 0;
      }


      //Initialise
      void initialise(const string& basename,const volume<float>& mask){
	volume4D<float> tmpvol;
	Matrix          tmpmat;
	
	cout<<"Load bedpostx samples"<<endl;
	if(fsl_imageexists(basename+"_thsamples")){
	  cout<<"1"<<endl;
	  read_volume4D(tmpvol,basename+"_thsamples");
	  tmpmat=tmpvol.matrix(mask);
	  thsamples.push_back(tmpmat);
	  cout<<"2"<<endl;
	  read_volume4D(tmpvol,basename+"_phsamples");
	  tmpmat=tmpvol.matrix(mask);
	  phsamples.push_back(tmpmat);
	  cout<<"3"<<endl;
	  read_volume4D(tmpvol,basename+"_fsamples");
	  tmpmat=tmpvol.matrix(mask);
	  fsamples.push_back(tmpmat);

	  lut_vol2mat = tmpvol.vol2matrixkey(mask);
	  nsamples    = tmpmat.Nrows();
	  nfibres     = 1;
	}
	else{
	  int fib=1;
	  bool fib_existed=true;
	  while(fib_existed){
	    if(fsl_imageexists(basename+"_th"+num2str(fib)+"samples")){
	      cout<<fib<<"_1"<<endl;
	      read_volume4D(tmpvol,basename+"_th"+num2str(fib)+"samples");
	      tmpmat=tmpvol.matrix(mask);
	      thsamples.push_back(tmpmat);
	      cout<<fib<<"_2"<<endl;
	      read_volume4D(tmpvol,basename+"_ph"+num2str(fib)+"samples");
	      tmpmat=tmpvol.matrix(mask);
	      phsamples.push_back(tmpmat);
	      cout<<fib<<"_3"<<endl;
	      read_volume4D(tmpvol,basename+"_f"+num2str(fib)+"samples");
	      tmpmat=tmpvol.matrix(mask);
	      fsamples.push_back(tmpmat);
	      fib++;
	    }
	    else{
	      fib_existed=false;
	    }
	  }
	  if(fib==1){
	      cerr<<"Could not find samples to load. Exit without doing anything"<<endl;
	      exit(1);
	  }
	  lut_vol2mat = tmpvol.vol2matrixkey(mask);
	  nsamples = thsamples[0].Nrows();
	  nfibres  = (int)thsamples.size();
	}
	copybasicproperties(mask,lut_vol2mat);

	cout<<endl;
	cout<<"nfibres  : "<<nfibres<<endl;
	cout<<"nsamples : "<<nsamples<<endl;
	cout<<endl;
	cout<<"Done loading samples."<<endl;
      }
      
      
      ColumnVector sample(const float& x,const float& y,const float&z,
			  const float& r_x,const float& r_y,const float& r_z,
			  float& prefer_x,float& prefer_y,float& prefer_z,
			  const int& sample_fib,int& sampled_fib,
			  int& newx,int& newy,int& newz){

	////////Probabilistic interpolation
	int cx =(int) ceil(x),fx=(int) floor(x);
	int cy =(int) ceil(y),fy=(int) floor(y);
	int cz =(int) ceil(z),fz=(int) floor(z);
	
	float pcx = (cx==fx)?1:(x-fx)/(cx-fx);
	float pcy = (cy==fy)?1:(y-fy)/(cy-fy);
	float pcz = (cz==fz)?1:(z-fz)/(cz-fz);
	
	newx = ((float)rand()/(float)RAND_MAX)>pcx?fx:cx;
	newy = ((float)rand()/(float)RAND_MAX)>pcy?fy:cy;
	newz = ((float)rand()/(float)RAND_MAX)>pcz?fz:cz;
	////////////////////////////////////	

	ColumnVector th_ph_f(3);	

	int col = lut_vol2mat(newx,newy,newz);
	if(col==0){//outside brain mask
	  th_ph_f=0;
	  return th_ph_f;
	}

	int samp=(int)round((float)rand()/float(RAND_MAX)*(float)(nsamples-1))+1;

	float theta=0,phi=0;
	float dotmax=0,dottmp=0;
	int fibind=0;
	if(nfibres>1){//more than 1 fibre
	  if(init_sample){//go for the specified fibre on the first jump or generate at random
	    fibst=sample_fibre(col,samp,opts.randfib.value());
	    theta=thsamples[fibst](samp,col);
	    phi=phsamples[fibst](samp,col);
	    init_sample=false;
	  }
	  else{
	    if(sample_fib>0){
	      fibind=sample_fibre(col,samp,sample_fib);	      
	      theta=thsamples[fibind](samp,col);
	      phi=phsamples[fibind](samp,col);
	    }
	    else{
	      if((fabs(prefer_x)+fabs(prefer_y)+fabs(prefer_z))==0){
		prefer_x=r_x;prefer_y=r_y;prefer_z=r_z;
	      }
	      for(int fib=0;fib<nfibres;fib++){
		if(fsamples[fib](samp,col)>opts.fibthresh.value()){
		  float phtmp=phsamples[fib](samp,col);
		  float thtmp=thsamples[fib](samp,col);
		  dottmp=fabs(sin(thtmp)*cos(phtmp)*prefer_x + sin(thtmp)*sin(phtmp)*prefer_y + cos(thtmp)*prefer_z);
		  if(dottmp>dotmax){
		    dotmax=dottmp;
		    theta=thtmp;
		    phi=phtmp;
		    fibind=fib;
		  }
		}
	      }
	      if(dotmax==0){
		theta=thsamples[0](samp,col);
		phi=phsamples[0](samp,col);
	      }
	    }
	  }
	}
	else{
	  theta=thsamples[0](samp,col);
	  phi=phsamples[0](samp,col);
	}
	
	float f;	
	if(usef){
	  f = fsamples[fibind](samp,col);
	}
	else{
	  f=1;
	}

	sampled_fib = fibind+1;

	th_ph_f(1)=theta;
	th_ph_f(2)=phi;
	th_ph_f(3)=f;
	return th_ph_f;
      }

      ColumnVector dimensions() const{
	ColumnVector dims(3);
	dims << lut_vol2mat.xdim() <<lut_vol2mat.ydim() << lut_vol2mat.zdim();
	return dims;
      }
    };
}

#endif



