/*  Copyright (C) 2004 University of Oxford  */

/*  CCOPYRIGHT  */

#include "pt_alltracts.h"

using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;
using namespace Utilities;
using namespace PARTICLE;
using namespace TRACTVOLS;
using namespace mesh;


void alltracts(){
  
  probtrackOptions& opts = probtrackOptions::getInstance();
  Log& logger = LogSingleton::getInstance();

  float xst,yst,zst,x,y,z;
  int nparticles=opts.nparticles.value();
  int nsteps=opts.nsteps.value();
  ///////////////////////////////////
  volume<int> mask;
  read_volume(mask,opts.maskfile.value());  
  volume<int> skipmask;
  if(opts.skipmask.value()!="") read_volume(skipmask,opts.skipmask.value());
  float lcrat=5;
  volume4D<float> loopcheck;
  if(opts.loopcheck.value()){
    loopcheck.reinitialize(int(ceil(mask.xsize()/lcrat)+1),int(ceil(mask.ysize()/lcrat)+1),int(ceil(mask.zsize()/lcrat)+1),3);
    loopcheck=0;
  }

  volume<int> Seeds;
  read_volume(Seeds,opts.seedfile.value());

  volume<int> RUBBISH;
  if(opts.rubbishfile.value()!=""){
    read_volume(RUBBISH,opts.rubbishfile.value());
  }
  
  volume<int> prob;
  volume<int> beenhere;
  beenhere=Seeds;
  beenhere=0;
  prob=beenhere;

  //   int numseeds=0;
  //   for(int Wz=Seeds.minz();Wz<=Seeds.maxz();Wz++){
  //     for(int Wy=Seeds.miny();Wy<=Seeds.maxy();Wy++){
  //       for(int Wx=Seeds.minx();Wx<=Seeds.maxx();Wx++){
  // 	if(Seeds.value(Wx,Wy,Wz)>0){
  // 	  numseeds++;
  // 	}
  //       }
  //     }
  //   }
  

  cout<<"Loading MCMC volumes"<<endl;
  TractVols vols(opts.usef.value());
  vols.initialise(opts.basename.value());
  
  Matrix Seeds_to_DTI;
  if(opts.seeds_to_dti.value()!=""){
    Seeds_to_DTI = read_ascii_matrix(opts.seeds_to_dti.value());
  }
  else{
    Seeds_to_DTI=Identity(4);
  }
  
  Matrix path(nsteps,3);
  path=1;
  
  float tmp2;
  ColumnVector th_ph_f;

  for(int Sz=Seeds.minz();Sz<=Seeds.maxz();Sz++){
    cout<<Sz<<endl;
    for(int Sy=Seeds.miny();Sy<=Seeds.maxy();Sy++){
      for(int Sx=Seeds.minx();Sx<=Seeds.maxx();Sx++){
	if(Seeds(Sx,Sy,Sz)!=0){
	  
	  ColumnVector xyz_seeds(3),dim_seeds(3),xyz_dti;
	  xyz_seeds << Sx << Sy << Sz;
	  dim_seeds <<Seeds.xdim()<<Seeds.ydim()<<Seeds.zdim();
	  
	  xyz_dti=vox_to_vox(xyz_seeds,dim_seeds,vols.dimensions(),Seeds_to_DTI);    
	  xst=xyz_dti(1);yst=xyz_dti(2);zst=xyz_dti(3);
	  Particle part(0,0,0,0,0,0,opts.steplength.value(),mask.xdim(),mask.ydim(),mask.zdim(),false);
	  
	  for( int p = 0; p < nparticles ; p++ ){
	    
	    for(int direc=1;direc<=2;direc++){
	      x=xst;y=yst;z=zst;
	      part.change_xyz(x,y,z);	    
	      if(direc==2){
		part.restart_reverse();  //go again in the opposite direction
	      }
	      for( int it = 1 ; it <= nsteps/2; it++){
		if( (mask( round(part.x()), round(part.y()), round(part.z())) != 0) ){
		
		  ///////////////////////////////////
		  //loopchecking
		  ///////////////////////////////////
		  if(opts.loopcheck.value()){
		    float oldrx=loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),0);
		    float oldry=loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),1);
		    float oldrz=loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),2);
		    if(part.rx()*oldrx+part.ry()*oldry+part.rz()*oldrz<0)
		      {
			// p--;
			break;
		      }
		    
		    loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),0)=part.rx();
		    loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),1)=part.ry();
		    loopcheck((int)round(part.x()/lcrat),(int)round(part.y()/lcrat),(int)round(part.z()/lcrat),2)=part.rz();  
		  
		  }
		
		
		  x=part.x();y=part.y();z=part.z();
		  xyz_dti <<x<<y<<z;
		  xyz_seeds=vox_to_vox(xyz_dti,vols.dimensions(),dim_seeds,Seeds_to_DTI.i());
		  int x_s =(int)round((float)xyz_seeds(1));
		  int y_s =(int)round((float)xyz_seeds(2));
		  int z_s =(int)round((float)xyz_seeds(3));
		  
		  if(opts.rubbishfile.value()!=""){
		    if(RUBBISH(x_s,y_s,z_s)>0) break;
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
		      if( (mask( round(part.x()), round(part.y()), round(part.z())) != 0) ){
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
		
		
		}

	      } // Close Step Number Loop
	      
	      if(opts.loopcheck.value()){
		loopcheck=0;
	      }
	    }
	    indexset(beenhere,path,0);
	    part.reset();

	  } // Close Particle Number Loop


	}
      }
    }
    


  } //Close Seed number Loop

  save_volume(prob,logger.appendDir(opts.outfile.value()));
}

