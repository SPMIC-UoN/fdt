/*  Copyright (C) 2004 University of Oxford  */

/*  CCOPYRIGHT  */

#include "pt_simple.h"

using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;
using namespace Utilities;
using namespace PARTICLE;
using namespace TRACTVOLS;
using namespace mesh;


void track(){
  probtrackOptions& opts =probtrackOptions::getInstance();
 
  ////////////////////////////
  Log& logger = LogSingleton::getInstance();
  if(opts.verbose.value()>1){
    logger.makeDir("particles","particle0",true,false);
  }
  ////////////////////////////////////
  float xst,yst,zst,x,y,z;
  int nparticles=opts.nparticles.value();
  int nsteps=opts.nsteps.value();
  ///////////////////////////////////
  volume<int> mask;
  volume<int> RUBBISH;
  volume<int> skipmask;
  volume<int> prob,beenhere;
  read_volume(mask,opts.maskfile.value());
  Matrix Seeds_to_DTI;
  if(opts.seedref.value()!=""){
    read_volume(prob,opts.seedref.value());
    beenhere=prob*0;
    read_ascii_matrix(Seeds_to_DTI,opts.seeds_to_dti.value());
  }
  else{
    prob=mask;
    beenhere=prob*0;
    Seeds_to_DTI=Identity(4);
  }
  
  float lcrat=5;
  volume4D<float> loopcheck((int)ceil(mask.xsize()/lcrat)+1,(int)ceil(mask.ysize()/lcrat)+1,(int)ceil(mask.zsize()/lcrat)+1,3);
  loopcheck=0;
  
  if(opts.rubbishfile.value()!="") read_volume(RUBBISH,opts.rubbishfile.value());
  if(opts.skipmask.value()!="") read_volume(skipmask,opts.skipmask.value());
  
  TractVols vols(opts.usef.value());
  vols.initialise(opts.basename.value());
  
  Matrix path(nsteps,3);
  path=1;

  
  Matrix Seeds = read_ascii_matrix(opts.seedfile.value());
  float tmp2;
  float randtmp1,randtmp2,randtmp3;
  ColumnVector th_ph_f;
  
  
  for(int SN=1; SN<=Seeds.Nrows();SN++){
    prob=0;
    ColumnVector xyz_seeds(3),dim_seeds(3),xyz_dti;

    
    if(opts.seedref.value()==""){
      xst=Seeds(SN,1);
      yst=Seeds(SN,2);
      zst=Seeds(SN,3);
    }
    else{
      xyz_seeds << Seeds(SN,1) << Seeds(SN,2) << Seeds(SN,3);
      dim_seeds <<prob.xdim()<<prob.ydim()<<prob.zdim();
      xyz_dti=vox_to_vox(xyz_seeds,dim_seeds,vols.dimensions(),Seeds_to_DTI); 
      xst=xyz_dti(1);yst=xyz_dti(2);zst=xyz_dti(3);

    }
      
    Particle part(0,0,0,0,0,0,opts.steplength.value(),mask.xdim(),mask.ydim(),mask.zdim(),false);
    
    for( int p = 0; p < nparticles ; p++ ){
      if(opts.verbose.value()>0){
	cout<<"particle number "<<p<<endl;
      }
      
      if(opts.verbose.value()>1)
	logger.setLogFile("particle"+num2str(p));
      if(Seeds.Ncols()==6){
	randtmp1=rand(); randtmp1/=RAND_MAX;
	randtmp2=rand(); randtmp2/=RAND_MAX;
	randtmp3=rand(); randtmp3/=RAND_MAX;
	xst = xst+randtmp1*Seeds(SN,4)/prob.xdim(); 
	yst = yst+randtmp2*Seeds(SN,5)/prob.ydim(); 
	zst = zst+randtmp3*Seeds(SN,6)/prob.zdim();
      }
      
      for(int direc=1;direc<=2;direc++){
	x=xst;y=yst;z=zst;
	part.change_xyz(x,y,z);	    
	if(direc==2){
	  part.restart_reverse();  //go again in the opposite direction
       	}
	
	for( int it = 1 ; it <= nsteps/2; it++){
	  if( (mask( round(part.x()), round(part.y()), round(part.z())) == 1) ){
	    if(opts.loopcheck.value()){
	      float oldrx=loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),0);
	      float oldry=loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),1);
	      float oldrz=loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),2);
	      if(part.rx()*oldrx+part.ry()*oldry+part.rz()*oldrz<0)
		{
		  break;
		}
	   
	      loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),0)=part.rx();
	      loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),1)=part.ry();
	      loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),2)=part.rz();
	      
	    }
	    
	    
	    if(opts.verbose.value()>1){
	      logger << part;
	    } 
	    
	    
	    int x_s,y_s,z_s;
	    
	    if(opts.seedref.value()!=""){
	      x=part.x();y=part.y();z=part.z();
	      xyz_dti <<x<<y<<z;
	      xyz_seeds=vox_to_vox(xyz_dti,vols.dimensions(),dim_seeds,Seeds_to_DTI.i());
	      x_s =(int)round((float)xyz_seeds(1));
	      y_s =(int)round((float)xyz_seeds(2));
	      z_s =(int)round((float)xyz_seeds(3));
	    }
	    else{
	      x_s=(int)round(part.x());
	      y_s=(int)round(part.y());
	      z_s=(int)round(part.z());
	    }
	    
	    if(opts.rubbishfile.value()!="")
	      {
		if(RUBBISH(x_s,y_s,z_s)==1) break;
	      }

	    path(it+(direc-1)*nsteps/2,1)=x_s; 
	    path(it+(direc-1)*nsteps/2,2)=y_s;
	    path(it+(direc-1)*nsteps/2,3)=z_s;
	    if(beenhere(x_s,y_s,z_s)==0){
	      prob(x_s,y_s,z_s)+=1;
	      beenhere(x_s,y_s,z_s)=1;
	    }
	    
	    
	    if(opts.skipmask.value() == ""){
	      th_ph_f=vols.sample(part.x(),part.y(),part.z());
	    }
	    else{
	      if(skipmask(x_s,y_s,z_s)==0)
		th_ph_f=vols.sample(part.x(),part.y(),part.z());
	    }
	    
	    
	    tmp2=rand(); tmp2/=RAND_MAX;
	    if(th_ph_f(3)>tmp2){
	      if(!part.check_dir(th_ph_f(1),th_ph_f(2),opts.c_thr.value())){
		break;
	      }
	      
	      if((th_ph_f(1)!=0&&th_ph_f(2)!=0)){
		if(!opts.modeuler.value())
		  part.jump(th_ph_f(1),th_ph_f(2));
		else
		  {
		    ColumnVector test_th_ph_f;
		    part.testjump(th_ph_f(1),th_ph_f(2));
		    test_th_ph_f=vols.sample(part.testx(),part.testy(),part.testz());
		    part.jump(test_th_ph_f(1),test_th_ph_f(2));
		  }
		
	      }
	      

	    }
	  }
	} // Close Step Number Loop
	
	if(opts.loopcheck.value()){
	  loopcheck=0;}
      } //   close direction loop
      
      part.reset();
      indexset(beenhere,path,0);
      
    } // Close Particle Number Loop    
    string thisout=opts.outfile.value()+num2str(Seeds(SN,1))+(string)"_"+num2str(Seeds(SN,2))+(string)"_"+num2str(Seeds(SN,3));
    save_volume(prob,thisout);
  } //Close Seed number Loop
}
