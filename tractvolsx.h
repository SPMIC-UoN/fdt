/*  tractvolsx.h

    Tim Behrens, FMRIB Image Analysis Group

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
      vector<volume4D<float>* > thsamples;
      vector<volume4D<float>* > phsamples;
      vector<volume4D<float>* > fsamples;
      bool init_sample;
      int fibst;
      bool usef;
      
    public:
      //constructors::
      Tractvolsx(const bool& usefin=false):opts(probtrackxOptions::getInstance()),init_sample(true),fibst(1),usef(usefin){}
      Tractvolsx():opts(probtrackxOptions::getInstance()){}
      ~Tractvolsx(){
	for(unsigned int m=0;m<thsamples.size();m++)
	  delete thsamples[m]; //ask flitney, do you just delete the ptr??
	for(unsigned int m=0;m<phsamples.size();m++)
	  delete phsamples[m];
	for(unsigned int m=0;m<fsamples.size();m++)
	  delete fsamples[m];
      }
      inline int nfibres()const{return (int)thsamples.size();}
      
      void reset(const int& fibst_in){
	init_sample=true;
	fibst=fibst_in;
      }
      //Initialise
      void initialise(const string& basename){
	

	if(fsl_imageexists(basename+"_thsamples")){
	  volume4D<float> *tmpthptr= new volume4D<float>;
	  volume4D<float> *tmpphptr= new volume4D<float>;
	  volume4D<float> *tmpfptr= new volume4D<float>;
	  cout<<"1"<<endl;
	  read_volume4D(*tmpthptr,basename+"_thsamples");
	  cout<<"2"<<endl;
	  thsamples.push_back(tmpthptr);
	  cout<<"3"<<endl;
	  read_volume4D(*tmpphptr,basename+"_phsamples");
	  cout<<"4"<<endl;
	  phsamples.push_back(tmpphptr);
	  cout<<"5"<<endl;
	  if(usef){
	    read_volume4D(*tmpfptr,basename+"_fsamples");
	    fsamples.push_back(tmpfptr);
	  }
	  cout<<"6"<<endl;
	}
	else{
	  int fib=1;
	  bool fib_existed=true;
	  while(fib_existed){
	    if(fsl_imageexists(basename+"_th"+num2str(fib)+"samples")){
	      volume4D<float> *tmpthptr= new volume4D<float>;
	      volume4D<float> *tmpphptr= new volume4D<float>;
	      volume4D<float> *tmpfptr= new volume4D<float>;
	      cout<<fib<<"_1"<<endl;
	      read_volume4D(*tmpthptr,basename+"_th"+num2str(fib)+"samples");
	      thsamples.push_back(tmpthptr);
	      cout<<fib<<"_2"<<endl;
	      read_volume4D(*tmpphptr,basename+"_ph"+num2str(fib)+"samples");
	      phsamples.push_back(tmpphptr);
	      cout<<fib<<"_3"<<endl;
	      read_volume4D(*tmpfptr,basename+"_f"+num2str(fib)+"samples");
	      fsamples.push_back(tmpfptr);
	      fib++;
	    }
	    else{
	      fib_existed=false;
	    }
	  }
	  
	}
	cout<<"7"<<endl;
      }
      
      
      ColumnVector sample(const float& x,const float& y,const float&z,const float& r_x,const float& r_y,const float& r_z,
			  float& prefer_x,float& prefer_y,float& prefer_z){

	////////Probabilistic interpolation
	int cx =(int) ceil(x),fx=(int) floor(x);
	int cy =(int) ceil(y),fy=(int) floor(y);
	int cz =(int) ceil(z),fz=(int) floor(z);
	
	//cerr<<x<<" "<<y<<" "<<z<<" "<<cx<<" "<<cy<<" "<<cz<<" "<<fx<<" "<<fy<<" "<<fz<<endl;
	float pcx,pcy,pcz;
	if(cx==fx)
	  pcx=1;
	else
	  pcx=(x-fx)/(cx-fx);
	
	if(cy==fy)
	  pcy=1;
	else
	  pcy=(y-fy)/(cy-fy);
	
	if(cz==fz)
	  pcz=1;
	else
	  pcz=(z-fz)/(cz-fz);
	
	///////new xyz values from probabilistic interpolation
	int newx,newy,newz; 
	float tmp=rand(); tmp/=RAND_MAX;
	if(tmp>pcx)
	  newx=fx;
	else
	  newx=cx;
	
	tmp=rand(); tmp/=RAND_MAX;
	if(tmp>pcy)
	  newy=fy;
	else
	  newy=cy;
	
	tmp=rand(); tmp/=RAND_MAX;
	if(tmp>pcz)
	  newz=fz;
	else
	  newz=cz;
 
	ColumnVector th_ph_f(3);	
	float samp=rand(); samp/=RAND_MAX;
	samp=round(samp*((*thsamples[0]).tsize()-1));
	float theta=0,phi=0;
	float dotmax=0,dottmp=0;
	int fibind=0;
	if(thsamples.size()>1){//more than 1 fibre
	  if(init_sample){//go for the specified option on the first jump
	    theta=(*thsamples[fibst])(int(newx),int(newy),int(newz),int(samp));
	    phi=(*phsamples[fibst])(int(newx),int(newy),int(newz),int(samp));
	    init_sample=false;
	  }
	  else{
	    if((fabs(prefer_x)+fabs(prefer_y)+fabs(prefer_z))==0){
	      prefer_x=r_x;prefer_y=r_y;prefer_z=r_z;
	    }
	    for(unsigned int fib=0;fib<thsamples.size();fib++){
	      if((*fsamples[fib])(int(newx),int(newy),int(newz),int(samp))>opts.fibthresh.value()){
		float phtmp=(*phsamples[fib])(int(newx),int(newy),int(newz),int(samp));
		float thtmp=(*thsamples[fib])(int(newx),int(newy),int(newz),int(samp));
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
	      theta=(*thsamples[0])(int(newx),int(newy),int(newz),int(samp));
	      phi=(*phsamples[0])(int(newx),int(newy),int(newz),int(samp));
	    }
	  }
	}
	else{
	  theta=(*thsamples[0])(int(newx),int(newy),int(newz),int(samp));
	  phi=(*phsamples[0])(int(newx),int(newy),int(newz),int(samp));
	}

	
	float f;
	
	if(usef){
	  f = (*fsamples[fibind])(int(newx),int(newy),int(newz),int(samp));
	}
	else
	  f=1;
	
	th_ph_f(1)=theta;
	th_ph_f(2)=phi;
	th_ph_f(3)=f;
	return th_ph_f;
      }

      ColumnVector dimensions() const{
	ColumnVector dims(3);
	dims << (*thsamples[0]).xdim() <<(*thsamples[0]).ydim() << (*thsamples[0]).zdim();
	return dims;
      }
    };
}

#endif



