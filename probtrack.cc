#include <iostream>
#include <fstream>
#include "newimage/newimageall.h"
#include "utils/log.h"
#include "meshclass/meshclass.h"
#include "probtrackOptions.h"
#include "particle.h"
#include "tractvols.h"

using namespace std;
using namespace NEWIMAGE;
using namespace TRACT;
using namespace Utilities;
using namespace PARTICLE;
using namespace TRACTVOLS;
using namespace mesh;
//using namespace NEWMAT;
//////////////////////////
/////////////////////////




 void read_masks(vector<string>& masks, const string& filename){
   ifstream fs(filename.c_str());
   string tmp;
   if(fs){
     fs>>tmp;
     while(!fs.eof()){
       masks.push_back(tmp);
       fs>>tmp;
     }
   }
   else{
     cerr<<filename<<" does not exist"<<endl;
     exit(0);
   }
 }

inline const double pi(){return 3.1416;}

ColumnVector vox_to_vox(const ColumnVector& xyz1,const ColumnVector& dims1,const ColumnVector& dims2,const Matrix& xfm){
  ColumnVector xyz1_mm(4),xyz2_mm,xyz2(3);
  xyz1_mm<<xyz1(1)*dims1(1)<<xyz1(2)*dims1(2)<<xyz1(3)*dims1(3)<<1;
  xyz2_mm=xfm*xyz1_mm;
  xyz2_mm=xyz2_mm/xyz2_mm(4);
  xyz2<<xyz2_mm(1)/dims2(1)<<xyz2_mm(2)/dims2(2)<<xyz2_mm(3)/dims2(3);
  return xyz2;
}


ColumnVector mni_to_imgvox(const ColumnVector& mni,const ColumnVector& mni_origin,const Matrix& mni2img, const ColumnVector& img_dims){
  ColumnVector mni_new_origin(4),img_mm;//homogeneous
  ColumnVector img_vox(3);
  mni_new_origin<<mni(1)+mni_origin(1)<<mni(2)+mni_origin(2)<<mni(3)+mni_origin(3)<<1;
  img_mm=mni2img*mni_new_origin;
  img_vox<<img_mm(1)/img_dims(1)<<img_mm(2)/img_dims(2)<<img_mm(3)/img_dims(3);
  return img_vox;
}




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
    read_ascii_matrix(Seeds_to_DTI,opts.seeds_to_dti.value());
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
	if(Seeds(Sx,Sy,Sz)>0){
	  
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
		if( (mask( round(part.x()), round(part.y()), round(part.z())) > 0) ){
		
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
		    if(skipmask((int)round(part.x()),(int)round(part.y()),(int)round(part.z()))==0)
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


void twomasks_symm(){
  
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

  volume<int> Seeds,mask2;
  read_volume(Seeds,opts.seedfile.value());

  if(opts.mask2.value()!=""){
    read_volume(mask2,opts.mask2.value());
    Seeds=Seeds+(2*mask2);
  }


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
    read_ascii_matrix(Seeds_to_DTI,opts.seeds_to_dti.value());
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
	  if(Seeds(Sx,Sy,Sz)>0){
	    int	 maskno=Seeds(Sx,Sy,Sz);
	    ColumnVector xyz_seeds(3),dim_seeds(3),xyz_dti;
	    xyz_seeds << Sx << Sy << Sz;
	    dim_seeds <<Seeds.xdim()<<Seeds.ydim()<<Seeds.zdim();
	  
	    xyz_dti=vox_to_vox(xyz_seeds,dim_seeds,vols.dimensions(),Seeds_to_DTI);    
	    xst=xyz_dti(1);yst=xyz_dti(2);zst=xyz_dti(3);
	    Particle part(0,0,0,0,0,0,opts.steplength.value(),mask.xdim(),mask.ydim(),mask.zdim(),false);

	    for( int p = 0; p < nparticles ; p++ ){
	      bool keepflag=false;
	      ColumnVector partlength(2);
	      partlength=0;
	      for(int direc=1;direc<=2;direc++){
		x=xst;y=yst;z=zst;
		part.change_xyz(x,y,z);
		if(direc==2){
		  part.restart_reverse();  //go again in the opposite direction
		}
	      
		for( int it = 1 ; it <= nsteps/2; it++){
		  if( (mask( round(part.x()), round(part.y()), round(part.z())) > 0) ){
		    
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
		    partlength(direc)+=1;
		    if( Seeds(x_s,y_s,z_s)==(maskno*-1 + 3) ){ //mult by -1 and add 3 makes 1 into 2 and 2 into 1
		      keepflag=true;
		    }
		    if(opts.skipmask.value() == ""){
		      th_ph_f=vols.sample(part.x(),part.y(),part.z());
		    }
		    else{
		      if(skipmask((int)round(part.x()),(int)round(part.y()),(int)round(part.z()))==0)
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
	      // Do the counting after a particle has finished
	      // And only if keepflag is set
	      if(keepflag){

		for(int direc=1;direc<=2;direc++){
		  for(int s=1;s<=int(partlength(direc));s++){
		    int x_s=int(path(s+(direc-1)*nsteps/2,1));
		    int y_s=int(path(s+(direc-1)*nsteps/2,2));
		    int z_s=int(path(s+(direc-1)*nsteps/2,3));
		    if(beenhere(x_s,y_s,z_s)==0){
		      prob(x_s,y_s,z_s)+=1;
		      beenhere(x_s,y_s,z_s)=1;
		    }
		  }
		}
		
		indexset(beenhere,path,0);
	      }
	      part.reset();

	   	    
	    
	    } // Close Particle Number Loop


	  }

	}// close x loop
      }//close y loop
   

    }//close z loop
  save_volume(prob,logger.appendDir(opts.outfile.value()));
}


void twomasks_asymm(){
  
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


  volume<int> Seeds,mask2;
  read_volume(Seeds,opts.seedfile.value());
  
  if(opts.mask2.value()!=""){
    read_volume(mask2,opts.mask2.value());
    Seeds=Seeds+(2*mask2);
  } 


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
    read_ascii_matrix(Seeds_to_DTI,opts.seeds_to_dti.value());
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
	  if(Seeds(Sx,Sy,Sz)==1){
	    int	 maskno=Seeds(Sx,Sy,Sz);
	    ColumnVector xyz_seeds(3),dim_seeds(3),xyz_dti;
	    xyz_seeds << Sx << Sy << Sz;
	    dim_seeds <<Seeds.xdim()<<Seeds.ydim()<<Seeds.zdim();
	  
	    xyz_dti=vox_to_vox(xyz_seeds,dim_seeds,vols.dimensions(),Seeds_to_DTI);    
	    xst=xyz_dti(1);yst=xyz_dti(2);zst=xyz_dti(3);
	    Particle part(0,0,0,0,0,0,opts.steplength.value(),mask.xdim(),mask.ydim(),mask.zdim(),false);

	    for( int p = 0; p < nparticles ; p++ ){
	      bool keepflag=false;
	      ColumnVector partlength(2);
	      partlength=0;
	      for(int direc=1;direc<=2;direc++){
		x=xst;y=yst;z=zst;
		part.change_xyz(x,y,z);
		if(direc==2){
		  part.restart_reverse();  //go again in the opposite direction
		}
	      
		for( int it = 1 ; it <= nsteps/2; it++){
		  if( (mask( round(part.x()), round(part.y()), round(part.z())) > 0) ){
		    
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
		    partlength(direc)+=1;
		    if( Seeds(x_s,y_s,z_s)==(maskno*-1 + 3) ){ //mult by -1 and add 3 makes 1 into 2 and 2 into 1
		      keepflag=true;
		    }
		    
		    
		     
		    if(opts.skipmask.value() == ""){
		      th_ph_f=vols.sample(part.x(),part.y(),part.z());
		    }
		    else{
		      if(skipmask((int)round(part.x()),(int)round(part.y()),(int)round(part.z()))==0)
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
	      // Do the counting after a particle has finished
	      // And only if keepflag is set
	      if(keepflag){

		for(int direc=1;direc<=2;direc++){
		  for(int s=1;s<=int(partlength(direc));s++){
		    int x_s=int(path(s+(direc-1)*nsteps/2,1));
		    int y_s=int(path(s+(direc-1)*nsteps/2,2));
		    int z_s=int(path(s+(direc-1)*nsteps/2,3));
		    if(beenhere(x_s,y_s,z_s)==0){
		      prob(x_s,y_s,z_s)+=1;
		      beenhere(x_s,y_s,z_s)=1;
		    }
		  }
		}
		
		indexset(beenhere,path,0);
	      }
	      part.reset();

	    
	    
	    
	    } // Close Particle Number Loop


	  }

	}// close x loop
      }//close y loop
   

    }//close z loop
  save_volume(prob,logger.appendDir(opts.outfile.value()));
}





void matrix2(){
  
  probtrackOptions& opts =probtrackOptions::getInstance();
  Log& logger = LogSingleton::getInstance();

  ////////////////////////////
  //  Log& logger = LogSingleton::getInstance();
  //  logger.makeDir(opts.logdir.value(),"logfile",true,false);
  ////////////////////////////////////
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
  //////////////////////////////////////////////
  // Segmented Volumes                        //
  //////////////////////////////////////////////
 //   vector<string> masknames;
//    read_masks(masknames,opts.cortexfile.value());
//    vector<volume<int> > cortex_masks;
//    vector<volume<int> > thal_segs;
//    volume<int> tmpcort;
//    cerr<<"Number of masks "<<masknames.size()<<endl;
//    for( unsigned int m = 0; m < masknames.size(); m++ ){
//      cerr<<"Reading "<<masknames[m]<<endl;
//      read_volume(tmpcort,masknames[m]);
//      cortex_masks.push_back(tmpcort);
//      tmpcort=0;
//      thal_segs.push_back(tmpcort);
//  }


  volume<int> Seeds;
  read_volume(Seeds,opts.seedfile.value());

  volume<int> RUBBISH;
  if(opts.rubbishfile.value()!=""){
    read_volume(RUBBISH,opts.rubbishfile.value());
  }

  volume<int> prob;
  volume<int> beenhere;
  prob=Seeds;prob=0;

  int numseeds=0;
  for(int Wz=Seeds.minz();Wz<=Seeds.maxz();Wz++){
    for(int Wy=Seeds.miny();Wy<=Seeds.maxy();Wy++){
      for(int Wx=Seeds.minx();Wx<=Seeds.maxx();Wx++){
	if(Seeds.value(Wx,Wy,Wz)>0){
	  numseeds++;
	}
      }
    }
  }

  volume<int> ConMat;
  volume<int> CoordMat(numseeds,3,1);
  volume<int> CoordMat_tracts_om; //for storing tractspace coords for othermatrix
  volume<int> lrmask;

  
  Matrix tempy;
  volume4D<int> lookup;
  read_volume(lrmask,opts.lrmask.value());
  beenhere=lrmask;beenhere=0;
  int numnz=0;
  for(int Wz=lrmask.minz();Wz<=lrmask.maxz();Wz++){
    for(int Wy=lrmask.miny();Wy<=lrmask.maxy();Wy++){
      for(int Wx=lrmask.minx();Wx<=lrmask.maxx();Wx++){
	if(lrmask.value(Wx,Wy,Wz)>0){
	  numnz++;
	}
      }
    }
  }
  if(numnz> pow(2,(float)sizeof(short)*8-1)-1){
    cerr<<"Output matrix too big for AVW - stopping."<<endl;
    cerr<<" Remember - you can store your tracts in "<<endl;
    cerr<<" low res even if you want your seeds in high res"<<endl;
    cerr<<" Just subsample the structural space mask"<<endl;
    cerr<<" Although, it must stay in line with the seeds"<<endl;
    exit(-1);
  }
  
  //    cerr<<"WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
  //      cerr<<"You will need significantly more than "<<numseeds*numnz*2<<" bytes of free memory!!"<<endl;
  
  ConMat.reinitialize(numseeds,numnz,1);
  ConMat=0;
  tempy.ReSize(numnz,1);
  for(int i=1;i<=numnz;i++){tempy(i,1)=i-1;}
  lookup.addvolume(lrmask);
  lookup.setmatrix(tempy.t(),lrmask);
  
  CoordMat_tracts_om.reinitialize(numnz,3,1);//store the tract space coordmat.
  int mytrow=0;
  for(int Wz=lrmask.minz();Wz<=lrmask.maxz();Wz++){
    for(int Wy=lrmask.miny();Wy<=lrmask.maxy();Wy++){
      for(int Wx=lrmask.minx();Wx<=lrmask.maxx();Wx++){
	if(lrmask(Wx,Wy,Wz)>0){
	  CoordMat_tracts_om(mytrow,0,0)=Wx;
	  CoordMat_tracts_om(mytrow,1,0)=Wy;
	  CoordMat_tracts_om(mytrow,2,0)=Wz;
	  mytrow++;
	}
      }
    }
  }
  
  
  
  int myrow=0;
  for(int Wz=Seeds.minz();Wz<=Seeds.maxz();Wz++){
    for(int Wy=Seeds.miny();Wy<=Seeds.maxy();Wy++){
      for(int Wx=Seeds.minx();Wx<=Seeds.maxx();Wx++){
	if(Seeds(Wx,Wy,Wz)>0){
	  CoordMat(myrow,0,0)=Wx;
	  CoordMat(myrow,1,0)=Wy;
	  CoordMat(myrow,2,0)=Wz;
	  myrow++;
	}
      }
    }
  }
  
  cout<<"Loading MCMC volumes"<<endl;
  TractVols vols(opts.usef.value());
  vols.initialise(opts.basename.value());
  
  Matrix Seeds_to_DTI;
  if(opts.seeds_to_dti.value()!=""){
    read_ascii_matrix(Seeds_to_DTI,opts.seeds_to_dti.value());
  }
  else{
    Seeds_to_DTI=Identity(4);
  }
  
  Matrix path(nsteps,3);
  path=1;
    
  float tmp2;
  ColumnVector th_ph_f;
  int other_conrow=0;
  int other_concol=0;

  for(int Sz=Seeds.minz();Sz<=Seeds.maxz();Sz++){
    cout<<Sz<<endl;
    for(int Sy=Seeds.miny();Sy<=Seeds.maxy();Sy++){
      for(int Sx=Seeds.minx();Sx<=Seeds.maxx();Sx++){
	if(Seeds(Sx,Sy,Sz)>0){
	  
	  ColumnVector xyz_seeds(3),dim_seeds(3),xyz_dti,xyz_othermatrix_tracts(3),dim_othermatrix_tracts(3);
	  xyz_seeds << Sx << Sy << Sz;
	  dim_seeds <<Seeds.xdim()<<Seeds.ydim()<<Seeds.zdim();
	  dim_othermatrix_tracts <<lrmask.xdim()<<lrmask.ydim()<<lrmask.zdim();
	  
	  
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
		if( (mask( round(part.x()), round(part.y()), round(part.z())) > 0) ){
		
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
		
		  ////////////////////////////////////////

		  x=part.x();y=part.y();z=part.z();
		  xyz_dti <<x<<y<<z;
		  xyz_seeds=vox_to_vox(xyz_dti,vols.dimensions(),dim_seeds,Seeds_to_DTI.i());
		  int x_s =(int)round((float)xyz_seeds(1));
		  int y_s =(int)round((float)xyz_seeds(2));
		  int z_s =(int)round((float)xyz_seeds(3));

		  if(opts.rubbishfile.value()!=""){
		    if(RUBBISH(x_s,y_s,z_s)>0) break;
		  }
		  
		  //find out where we are in the space of "lrmask" (which is in alignment with Seeds, but not 
		  // necessarily the same resolution
		  xyz_othermatrix_tracts=vox_to_vox(xyz_dti,vols.dimensions(),dim_othermatrix_tracts,Seeds_to_DTI.i()); 
		  int x_omt=(int)round((float)xyz_othermatrix_tracts(1));
		  int y_omt=(int)round((float)xyz_othermatrix_tracts(2));
		  int z_omt=(int)round((float)xyz_othermatrix_tracts(3));
		  path(it+(direc-1)*nsteps/2,1)=x_omt; 
		  path(it+(direc-1)*nsteps/2,2)=y_omt;
		  path(it+(direc-1)*nsteps/2,3)=z_omt;
		  
		  // Find out where this lrmask voxel is in the unwrapped matrix
		  other_concol=lookup(x_omt,y_omt,z_omt,0);
		  //load up the matrix
		  if(other_concol!=0){
		    if(beenhere(x_omt,y_omt,z_omt)==0){
		      ConMat(other_conrow,other_concol,0)+=1;
		      beenhere(x_omt,y_omt,z_omt)=1;
		    }
		  }
		  
		    
		  if(opts.skipmask.value() == ""){
		    th_ph_f=vols.sample(part.x(),part.y(),part.z());
		  }
		  else{
		    if(skipmask((int)round(part.x()),(int)round(part.y()),(int)round(part.z()))==0)
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
	    }//close direc loop
	    indexset(beenhere,path,0);
	    part.reset();

	  } // Close Particle Number Loop
	  other_conrow++;
	}
      }
    }
    


  } //Close Seed number Loop
  
  save_volume(ConMat,logger.appendDir(opts.outfile.value()));
  save_volume(CoordMat,logger.appendDir("coords_for_"+opts.outfile.value()));
  save_volume(CoordMat_tracts_om,logger.appendDir("tract_space_coords_for_"+opts.outfile.value()));
  save_volume4D(lookup,logger.appendDir("lookup_tractspace_"+opts.outfile.value()));
  
}




void matrix1(){
  
  probtrackOptions& opts =probtrackOptions::getInstance();
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
  prob=Seeds;prob=0;
  beenhere=Seeds;beenhere=0;
 
//  prob=tmpcort;beenhere=tmpcort;
  int numseeds=0;
  for(int Wz=Seeds.minz();Wz<=Seeds.maxz();Wz++){
    for(int Wy=Seeds.miny();Wy<=Seeds.maxy();Wy++){
      for(int Wx=Seeds.minx();Wx<=Seeds.maxx();Wx++){
	if(Seeds.value(Wx,Wy,Wz)>0){
	  numseeds++;
	}
      }
    }
  }

  volume<int> ConMat;
  volume<int> CoordMat(numseeds,3,1);
  ConMat.reinitialize(numseeds,numseeds,1);
  ConMat=0;
  
  
  int myrow=0;
  for(int Wz=Seeds.minz();Wz<=Seeds.maxz();Wz++){
    for(int Wy=Seeds.miny();Wy<=Seeds.maxy();Wy++){
      for(int Wx=Seeds.minx();Wx<=Seeds.maxx();Wx++){
	if(Seeds(Wx,Wy,Wz)>0){
	  CoordMat(myrow,0,0)=Wx;
	  CoordMat(myrow,1,0)=Wy;
	  CoordMat(myrow,2,0)=Wz;
	  myrow++;
	}
      }
    }
  }
  
  cout<<"Loading MCMC volumes"<<endl;
  TractVols vols(opts.usef.value());
  vols.initialise(opts.basename.value());
  
  Matrix Seeds_to_DTI;
  if(opts.seeds_to_dti.value()!=""){
    read_ascii_matrix(Seeds_to_DTI,opts.seeds_to_dti.value());
  }
  else{
    Seeds_to_DTI=Identity(4);
  }
  
  Matrix path(nsteps,3);
  path=1;
  
    
  float tmp2;
  ColumnVector th_ph_f;
  int Conrow=0;
  int hitcount=0;

  for(int Sz=Seeds.minz();Sz<=Seeds.maxz();Sz++){
    cout<<Sz<<endl;
    for(int Sy=Seeds.miny();Sy<=Seeds.maxy();Sy++){
      for(int Sx=Seeds.minx();Sx<=Seeds.maxx();Sx++){
	if(Seeds(Sx,Sy,Sz)>0){
	  //	  other_conrow++;
	  //	  int repeatcount=0;
	  ColumnVector xyz_seeds(3),dim_seeds(3),xyz_dti;
	  xyz_seeds << Sx << Sy << Sz;
	  dim_seeds <<Seeds.xdim()<<Seeds.ydim()<<Seeds.zdim();	  
	  xyz_dti=vox_to_vox(xyz_seeds,dim_seeds,vols.dimensions(),Seeds_to_DTI);    
	  xst=xyz_dti(1);yst=xyz_dti(2);zst=xyz_dti(3);
	  Particle part(0,0,0,0,0,0,opts.steplength.value(),mask.xdim(),mask.ydim(),mask.zdim(),false);
	  hitcount=0;
	  for( int p = 0; p < nparticles ; p++ ){
	    
	    for(int direc=1;direc<=2;direc++){
	      x=xst;y=yst;z=zst;
	      part.change_xyz(x,y,z);	    
	      if(direc==2){
		part.restart_reverse();  //go again in the opposite direction
	      }
	      for( int it = 1 ; it <= nsteps/2; it++){
		if( (mask( round(part.x()), round(part.y()), round(part.z())) > 0) ){
		
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
		
		  ////////////////////////////////////////


		  //		if(opts.verbose.value()>2){
		  //		  logger << part;
		  //		} 

		 
		  x=part.x();y=part.y();z=part.z();
		  xyz_dti <<x<<y<<z;
		  xyz_seeds=vox_to_vox(xyz_dti,vols.dimensions(),dim_seeds,Seeds_to_DTI.i());
		  int x_s =(int)round((float)xyz_seeds(1));
		  int y_s =(int)round((float)xyz_seeds(2));
		  int z_s =(int)round((float)xyz_seeds(3));
		  
		  if(opts.rubbishfile.value()!=""){
		    if(RUBBISH(x_s,y_s,z_s)>0) break;
		  }
		  
		  if(beenhere(x_s,y_s,z_s)==0){
		    prob(x_s,y_s,z_s)+=1;
		    beenhere(x_s,y_s,z_s)=1;
		  }
		  
		  if(opts.skipmask.value() == ""){
		    th_ph_f=vols.sample(part.x(),part.y(),part.z());
		  }
		  else{
		    if(skipmask((int)round(part.x()),(int)round(part.y()),(int)round(part.z()))==0)
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
	
	  
	  ///Store connectivity information in massive matrix -  matrix mode.  
	  int Concol=0;
	  for(int Wz=prob.minz();Wz<=prob.maxz();Wz++){
	    for(int Wy=prob.miny();Wy<=prob.maxy();Wy++){
	      for(int Wx=prob.minx();Wx<=prob.maxx();Wx++){
		if(Seeds(Wx,Wy,Wz)>0){
		  if(prob(Wx,Wy,Wz)>0){
		    ConMat(Conrow,Concol,0)=prob(Wx,Wy,Wz);
		  }
		  Concol++;
		}
		prob(Wx,Wy,Wz)=0;
		
	      }
	    }
	  }
	  
	  Conrow++;	
	}
	
      }
    }
  } //Close Seed number Loop
  save_volume(ConMat,logger.appendDir(opts.outfile.value()));
  save_volume(CoordMat,logger.appendDir("coords_for_"+opts.outfile.value()));
}




void maskmatrix(){
  
  probtrackOptions& opts =probtrackOptions::getInstance();
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
  prob=Seeds;prob=0;
  beenhere=prob;
 

  int numseeds=0;
  int  maxclusternum=0;
  for(int Wz=Seeds.minz();Wz<=Seeds.maxz();Wz++){
    for(int Wy=Seeds.miny();Wy<=Seeds.maxy();Wy++){
      for(int Wx=Seeds.minx();Wx<=Seeds.maxx();Wx++){
	if(Seeds.value(Wx,Wy,Wz)>0){
	  numseeds++;
	  if(Seeds.value(Wx,Wy,Wz) > maxclusternum){
	    maxclusternum=Seeds.value(Wx,Wy,Wz);
	  }
	}
      }
    }
  }

  Matrix meanconmat(maxclusternum,maxclusternum),maxconmat(maxclusternum,maxclusternum);
  maxconmat=0;meanconmat=0;
    
  cout<<"Loading MCMC volumes"<<endl;
  TractVols vols(opts.usef.value());
  vols.initialise(opts.basename.value());

  Matrix Seeds_to_DTI;
  if(opts.seeds_to_dti.value()!=""){
    read_ascii_matrix(Seeds_to_DTI,opts.seeds_to_dti.value());
  }
  else{
    Seeds_to_DTI=Identity(4);
  }

  Matrix path(nsteps,3);
  path=1;
  
    
  float tmp2;
  ColumnVector th_ph_f;
  int hitcount=0;

  
  
  for(int seedclust=1;seedclust<=maxclusternum;seedclust++){
    int clustsize=0;
    ColumnVector clustsums(maxclusternum),clustmaxes(maxclusternum);
    clustsums=0;clustmaxes=0;
    for(int Sz=Seeds.minz();Sz<=Seeds.maxz();Sz++){
      cout<<Sz<<endl;
      for(int Sy=Seeds.miny();Sy<=Seeds.maxy();Sy++){
	for(int Sx=Seeds.minx();Sx<=Seeds.maxx();Sx++){
	  if(Seeds(Sx,Sy,Sz)>0){

	    clustsize++;
	    ColumnVector flags(maxclusternum);  flags=0;
	    ColumnVector clustcounts(maxclusternum);  clustcounts=0;



	    ColumnVector xyz_seeds(3),dim_seeds(3),xyz_dti;
	    xyz_seeds << Sx << Sy << Sz;
	    dim_seeds <<Seeds.xdim()<<Seeds.ydim()<<Seeds.zdim();	  
	    xyz_dti=vox_to_vox(xyz_seeds,dim_seeds,vols.dimensions(),Seeds_to_DTI);    
	    xst=xyz_dti(1);yst=xyz_dti(2);zst=xyz_dti(3);
	    Particle part(0,0,0,0,0,0,opts.steplength.value(),mask.xdim(),mask.ydim(),mask.zdim(),false);
	    hitcount=0;
	    for( int p = 0; p < nparticles ; p++ ){
	    
	      for(int direc=1;direc<=2;direc++){
		x=xst;y=yst;z=zst;
		part.change_xyz(x,y,z);	    
		if(direc==2){
		  part.restart_reverse();  //go again in the opposite direction
		}
		for( int it = 1 ; it <= nsteps/2; it++){
		  if( (mask( round(part.x()), round(part.y()), round(part.z())) > 0) ){
		
		    ///////////////////////////////////
		    //loopchecking
		    ///////////////////////////////////
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
		    
		   
		    x=part.x();y=part.y();z=part.z();
		    xyz_dti <<x<<y<<z;
		    xyz_seeds=vox_to_vox(xyz_dti,vols.dimensions(),dim_seeds,Seeds_to_DTI.i());
		    int x_s =(int)round((float)xyz_seeds(1));
		    int y_s =(int)round((float)xyz_seeds(2));
		    int z_s =(int)round((float)xyz_seeds(3));
		    if(opts.rubbishfile.value()!=""){
		      if(RUBBISH(x_s,y_s,z_s)>0) break;
		    }
		    
		    int cn=Seeds(x_s,y_s,z_s);
		    if(cn!=0){
		      if(flags(cn)==0){
			clustcounts(cn)=clustcounts(cn)+1;flags(cn)=1;
		      }
		    }
		    
		    if(beenhere(x_s,y_s,z_s)==0){
		      prob(x_s,y_s,z_s)+=1;
		      beenhere(x_s,y_s,z_s)=1;
		    }
		    
		    if(opts.skipmask.value() == ""){
		      th_ph_f=vols.sample(part.x(),part.y(),part.z());
		    }
		    else{
		      if(skipmask((int)round(part.x()),(int)round(part.y()),(int)round(part.z()))==0)
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
	
	  
	    clustsums+=clustcounts;
	    for(int cn=1;cn<=maxclusternum;cn++){
	      if(clustcounts(cn)>clustmaxes(cn)){
		clustmaxes(cn)=clustcounts(cn);
	      }
	    }
	    
	  } 
	  
	}//close x seed
      } //close y seed
    }  //close z seed
    for(int targclust=1;targclust<=maxclusternum;targclust++){
      meanconmat(seedclust,targclust)=clustsums(targclust)/clustsize;
      maxconmat(seedclust,targclust)=clustmaxes(targclust);
    }
    
  } //close cluster number
  write_ascii_matrix(meanconmat,logger.appendDir("mean_"+opts.outfile.value()));
  write_ascii_matrix(maxconmat,logger.appendDir("max_"+opts.outfile.value()));
}






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
    prob=0;beenhere=prob;
    read_ascii_matrix(Seeds_to_DTI,opts.seeds_to_dti.value());
  }
  else{
    prob=mask;prob=0;
    beenhere=prob;
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
    string thisout=opts.outfile.value()+num2str(xst)+(string)"_"+num2str(yst)+(string)"_"+num2str(zst);
    save_volume(prob,thisout);
  } //Close Seed number Loop
}


void seeds_to_targets()
{ 
  probtrackOptions& opts =probtrackOptions::getInstance();
  Log& logger = LogSingleton::getInstance();
  //  opts.parse_command_line(argc,argv,logger);
  //   if(opts.verbose.value()>0){
  //     opts.status();
  //   }
  ////////////////////////////
  //  Log& logger = LogSingleton::getInstance();
  //  logger.makeDir(opts.logdir.value(),"logfile",true,false);
  
  ////////////////////////////////////
  
  float xst,yst,zst,x,y,z;
  int nparticles=opts.nparticles.value();
  int nsteps=opts.nsteps.value();
  ///////////////////////////////////
  volume<char> mask;
  read_volume(mask,opts.maskfile.value());  
  float lcrat=5;
  volume4D<float> loopcheck;
  if(opts.loopcheck.value()){
    loopcheck.reinitialize(int(ceil(mask.xsize()/lcrat)+1),int(ceil(mask.ysize()/lcrat)+1),int(ceil(mask.zsize()/lcrat)+1),3);
    loopcheck=0;
  }
  
  //////////////////////////////////////////////
  // Segmented Volumes                        //
  //////////////////////////////////////////////
  vector<string> masknames;
  read_masks(masknames,opts.targetfile.value());
  vector<volume<int> > target_masks;
  vector<volume<int> > thal_segs;
  volume<int> tmpcort;
  cout<<"Number of masks "<<masknames.size()<<endl;
  for( unsigned int m = 0; m < masknames.size(); m++ ){
    cout<<"Reading "<<masknames[m]<<endl;
    read_volume(tmpcort,masknames[m]);
    target_masks.push_back(tmpcort);
    tmpcort=0;
    thal_segs.push_back(tmpcort);
  }

  volume<int> Seeds;
  read_volume(Seeds,opts.seedfile.value());
  volume<int> skipmask;
  if(opts.skipmask.value()!="") read_volume(skipmask,opts.skipmask.value());
  volume<int> RUBBISH;
  if(opts.rubbishfile.value()!=""){
    read_volume(RUBBISH,opts.rubbishfile.value());
  }
  cout<<"Loading MCMC volumes"<<endl;
  TractVols vols(opts.usef.value());
  vols.initialise(opts.basename.value());
  
  Matrix Seeds_to_DTI;
  if(opts.seeds_to_dti.value()!=""){
    read_ascii_matrix(Seeds_to_DTI,opts.seeds_to_dti.value());
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
	if(Seeds(Sx,Sy,Sz)>0){
	  ColumnVector xyz_seeds(3),dim_seeds(3),xyz_dti;
	  xyz_seeds << Sx << Sy << Sz;
	  dim_seeds <<Seeds.xdim()<<Seeds.ydim()<<Seeds.zdim();
	  xyz_dti=vox_to_vox(xyz_seeds,dim_seeds,vols.dimensions(),Seeds_to_DTI);
	  xst=xyz_dti(1);yst=xyz_dti(2);zst=xyz_dti(3);
	  Particle part(0,0,0,0,0,0,opts.steplength.value(),mask.xdim(),mask.ydim(),mask.zdim(),false);
	  for( int p = 0; p < nparticles ; p++ ){
	    vector<int> flags;
	    for(unsigned int m=0;m<masknames.size();m++){flags.push_back(0);}

	    for(int direc=1;direc<=2;direc++){
	      x=xst;y=yst;z=zst;
	      part.change_xyz(x,y,z);	    
	      if(direc==2){
		part.restart_reverse();  //go again in the opposite direction
	      }
	      
	      
	      for( int it = 1 ; it < nsteps/2; it++){
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
		
		  //////////////////////////////////////////////////////////
		  // passes through    //
		  //////////////////////////////////////////////////////////
	 
		 
		  x=part.x();y=part.y();z=part.z();
		  xyz_dti <<x<<y<<z;
		  xyz_seeds=vox_to_vox(xyz_dti,vols.dimensions(),dim_seeds,Seeds_to_DTI.i());
		  int x_s =(int)round(float(xyz_seeds(1)));
		  int y_s =(int)round(float(xyz_seeds(2)));
		  int z_s =(int)round(float(xyz_seeds(3)));
		  

		  if(opts.rubbishfile.value()!=""){
		    if(RUBBISH(x_s,y_s,z_s)>0) break;
		  }

		  path(it+(direc-1)*nsteps/2,1)=round(part.x()); 
		  path(it+(direc-1)*nsteps/2,2)=round(part.y());
		  path(it+(direc-1)*nsteps/2,3)=round(part.z()); //stopping path in DTI space here
		  
		  for(unsigned int m=0;m<masknames.size();m++){
		    if(target_masks[m](x_s,y_s,z_s)>0 && flags[m]==0){
		      thal_segs[m](Sx,Sy,Sz)=thal_segs[m](Sx,Sy,Sz)+1;flags[m]=1;
		    }
		  }
		  
		  if(opts.skipmask.value() == ""){
		    th_ph_f=vols.sample(part.x(),part.y(),part.z());
		  }
		  else{
		    if(skipmask((int)round(part.x()),(int)round(part.y()),(int)round(part.z()))==0)
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
		
		} //close mask loop
		
	      } // Close Step Number Loop
	      if(opts.loopcheck.value()){
		loopcheck=0;
	      }
	    }
	    part.reset();
	  } // Close Particle Number Loop
	}
      }
    }
  } //Close Seed number Loop
  
  string odir=logger.getDir();

  for(unsigned int m=0;m<masknames.size();m++){
    string tmpname=masknames[m];

    int pos=tmpname.find("/",0);
    int lastpos=pos;
    
    while(pos>=0){
      lastpos=pos;
      pos=tmpname.find("/",pos);
      // replace / with _
      tmpname[pos]='_';
    }
    
    //only take things after the last pos
    tmpname=tmpname.substr(lastpos+1,tmpname.length()-lastpos-1);
    
    save_volume(thal_segs[m],odir+"/seeds_to_"+tmpname);
  }

}






void meshtrack(){
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
  if(opts.seedref.value()!=""){
    read_volume(prob,opts.seedref.value());
    prob=0;beenhere=prob;
    
  }
  else{
    prob=mask;prob=0;
    beenhere=prob;
  }
  
  Matrix Seeds_to_DTI;
  read_ascii_matrix(Seeds_to_DTI,opts.seeds_to_dti.value()); // Here seeds_to_dti should take the standard volume to diff space 
  
  float lcrat=5;
  volume4D<float> loopcheck((int)ceil(mask.xsize()/lcrat)+1,(int)ceil(mask.ysize()/lcrat)+1,(int)ceil(mask.zsize()/lcrat)+1,3);
  loopcheck=0;
  
  if(opts.rubbishfile.value()!="") read_volume(RUBBISH,opts.rubbishfile.value());
  if(opts.skipmask.value()!="") read_volume(skipmask,opts.skipmask.value());
  
  TractVols vols(opts.usef.value());
  vols.initialise(opts.basename.value());
  
  Matrix path(nsteps,3);
  path=1;

  float tmp2;
  float randtmp1,randtmp2,randtmp3;
  ColumnVector th_ph_f;
  
  Mesh mseeds;
  int ftype=mseeds.load(opts.meshfile.value()); 
  mseeds.load_fs_label(opts.seedfile.value());
  ColumnVector mni_origin(3),fs_coord_mm(3),xyz_dti,xyz_seeds,dim_seeds(3);
  dim_seeds<<prob.xdim()<<prob.ydim()<<prob.zdim(); //In seedref space if exists. Else in dti space
  mni_origin << 92 << 128 << 37;
  
  
  
  for (vector<Mpoint*>::iterator i = mseeds._points.begin(); i!=mseeds._points.end(); i++ ){
    if((*i)->get_value() >0){
      
      fs_coord_mm<<(*i)->get_coord().X<<(*i)->get_coord().Y<<(*i)->get_coord().Z; 
      xyz_dti=mni_to_imgvox(fs_coord_mm,mni_origin, Seeds_to_DTI,vols.dimensions()); //xyz_dti in voxels, not mm
      xst=xyz_dti(1);yst=xyz_dti(2);zst=xyz_dti(3); //xyz_dti in voxels,not mm
      
      
      
      Particle part(0,0,0,0,0,0,opts.steplength.value(),mask.xdim(),mask.ydim(),mask.zdim(),false);
      
      int length=0;
      for( int p = 0; p < nparticles ; p++ ){
	if(opts.verbose.value()>0){
	  cout<<"particle number "<<p<<endl;
	}
	
	if(opts.verbose.value()>1)
	  logger.setLogFile("particle"+num2str(p));
	
	//Don't have a direction loop as in other cases, as always want to track in from cortex.
	
	x=xst;y=yst;z=zst;
	part.change_xyz(x,y,z);	    
	part.set_dir((*i)->local_normal().X,(*i)->local_normal().Y,(*i)->local_normal().Z);//Set the start dir so that we track inwards from cortex
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
	      
	    path(it,1)=x_s; 
	    path(it,2)=y_s;
	    path(it,3)=z_s;

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
	      if(!part.check_dir(th_ph_f(1),th_ph_f(2),opts.c_thr.value()) && it!=1){ 
		//Don't do curvature checking on the first step as we have set the old direction to the surface normal 
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
	    length++;
	      
	  }
	    
	} // Close Step Number Loop
	  
	if(opts.loopcheck.value()){
	  loopcheck=0;}
	  
	part.reset();
	indexset(beenhere,path,0);
	
      } // Close Particle Number Loop
      (*i)->set_value(length/nparticles);

      //      string thisout=opts.outfile.value()+num2str(xst)+(string)"_"+num2str(yst)+(string)"_"+num2str(zst);
      // save_volume(prob,thisout);
    }
  } //Close Seed number Loop
  mseeds.save_fs_label(logger.AppendDir(opts.outfile.value()));
}




int main ( int argc, char **argv ){
  probtrackOptions& opts =probtrackOptions::getInstance();
  Log& logger = LogSingleton::getInstance();
  opts.parse_command_line(argc,argv,logger);
  if(opts.verbose.value()>0){
    opts.status();
  }
  if(opts.mode.value()=="simple")
    track();
  else if(opts.mode.value()=="seeds_to_targets")
    seeds_to_targets();
  else if(opts.mode.value()=="seedmask")
    alltracts();
  else if(opts.mode.value()=="twomasks_symm")
    twomasks_symm();
  else if(opts.mode.value()=="twomasks_asymm")
    twomasks_asymm();
  else if(opts.mode.value()=="matrix1")
    matrix1();
  else if(opts.mode.value()=="matrix2")
    matrix2();
  else if(opts.mode.value()=="maskmatrix")
    maskmatrix();
  else if(opts.mode.value()=="meshtrack")
    meshtrack();
  else{
    cout <<"Invalid setting for option  mode -- try setting mode=help"<<endl;
  }
      
  return 0;
}















