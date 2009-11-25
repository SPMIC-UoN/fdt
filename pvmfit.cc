/*  Copyright (C) 2009 University of Oxford  */



/*  CCOPYRIGHT  */

#include <iostream>
#include <cmath>
#include "miscmaths/miscmaths.h"
#include "miscmaths/nonlin.h"
#include "newmat.h"
#include "pvmfitOptions.h"
#include "newimage/newimageall.h"
#include "diffmodels.h"

using namespace std;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace PVMFIT;
using namespace NEWIMAGE;




int main(int argc, char** argv)
{
  //parse command line
  pvmfitOptions& opts = pvmfitOptions::getInstance();
  int success=opts.parse_command_line(argc,argv);
  if(!success) return 1;
   if(opts.verbose.value()){
    cout<<"data file "<<opts.datafile.value()<<endl;
    cout<<"mask file "<<opts.maskfile.value()<<endl;
    cout<<"bvecs     "<<opts.bvecsfile.value()<<endl;
    cout<<"bvals     "<<opts.bvalsfile.value()<<endl;
  }
  
  // Set random seed:
  Matrix bvecs = read_ascii_matrix(opts.bvecsfile.value());
  if(bvecs.Nrows()>3) bvecs=bvecs.t();
  for(int i=1;i<=bvecs.Ncols();i++){
    float tmpsum=sqrt(bvecs(1,i)*bvecs(1,i)+bvecs(2,i)*bvecs(2,i)+bvecs(3,i)*bvecs(3,i));
    if(tmpsum!=0){
      bvecs(1,i)=bvecs(1,i)/tmpsum;
      bvecs(2,i)=bvecs(2,i)/tmpsum;
      bvecs(3,i)=bvecs(3,i)/tmpsum;
    }  
  }
  Matrix bvals = read_ascii_matrix(opts.bvalsfile.value());
  if(bvals.Nrows()>1) bvals=bvals.t();


  volume4D<float> data;
  volume<int> mask;

  if(opts.verbose.value()) cout<<"reading data"<<endl;
  read_volume4D(data,opts.datafile.value());

  if(opts.verbose.value()) cout<<"reading mask"<<endl;
  read_volume(mask,opts.maskfile.value());

  if(opts.verbose.value()) cout<<"ok"<<endl;
  int minx=0;
  int maxx=mask.xsize();
  int miny=0;
  int maxy=mask.ysize();
  int minz=0;
  int maxz=mask.zsize();
  cout<<minx<<" "<<maxx<<" "<<miny<<" "<<maxy<<" "<<minz<<" "<<maxz<<endl;

  if(opts.verbose.value()) cout<<"setting up vols"<<endl;
  volume<float> S0(maxx-minx,maxy-miny,maxz-minz);
  volume<float> dvol(maxx-minx,maxy-miny,maxz-minz);
  volume<float> tmpvol(maxx-minx,maxy-miny,maxz-minz);
  volume4D<float> tmpvol4D(maxx-minx,maxy-miny,maxz-minz,3);

  vector< volume<float> > fvol,thvol,phvol;
  vector< volume4D<float> > dyads;

  if(opts.verbose.value()) cout<<"copying input properties to output volumes"<<endl;
  copybasicproperties(data[0],S0);
  copybasicproperties(data[0],dvol);
  copybasicproperties(data[0],tmpvol);
  copybasicproperties(data[0],tmpvol4D);


  tmpvol = 0;
  tmpvol4D = 0;
  for(int i=0;i<opts.nfibres.value();i++){
    fvol.push_back(tmpvol);
    thvol.push_back(tmpvol);
    phvol.push_back(tmpvol);
    dyads.push_back(tmpvol4D);
  }

  if(opts.verbose.value()) cout<<"zeroing output volumes"<<endl;
  S0=0;dvol=0;
  volume<float> dvol_std;
  if(opts.modelnum.value()==2){
    dvol_std.reinitialize(maxx-minx,maxy-miny,maxz-minz);
    dvol_std=0;
  }
  

  if(opts.verbose.value()) cout<<"ok"<<endl;

  //int counter=0;
  ColumnVector S(bvals.Ncols());
  if(opts.verbose.value()) cout<<"starting the fits"<<endl;
  for(int k = minz; k < maxz; k++){
    cout<<k<<" slices processed"<<endl;
    for(int j=miny; j < maxy; j++){
      for(int i =minx; i< maxx; i++){
	if(mask(i,j,k)==0)continue;
	for(int t=0;t < data.tsize();t++)
	  S(t+1)=data(i,j,k,t);

	if(opts.modelnum.value()==1){
	  PVM_single pvm(S,bvecs,bvals,opts.nfibres.value());
	  pvm.fit();
	  
	  S0(i-minx,j-miny,k-minz)   = pvm.get_s0();
	  dvol(i-minx,j-miny,k-minz) = pvm.get_d();
	  for(int f=0;f<opts.nfibres.value();f++){
	    fvol[f](i-minx,j-miny,k-minz)  = pvm.get_f(f+1);
	    thvol[f](i-minx,j-miny,k-minz) = pvm.get_th(f+1);
	    phvol[f](i-minx,j-miny,k-minz) = pvm.get_ph(f+1);
	  }
	}
	else{
	  PVM_multi pvm(S,bvecs,bvals,opts.nfibres.value());
	  pvm.fit();
	  
	  S0(i-minx,j-miny,k-minz)   = pvm.get_s0();
	  dvol(i-minx,j-miny,k-minz) = pvm.get_d();
	  dvol_std(i-minx,j-miny,k-minz) = pvm.get_d_std();
	  for(int f=0;f<opts.nfibres.value();f++){
	    fvol[f](i-minx,j-miny,k-minz)  = pvm.get_f(f+1);
	    thvol[f](i-minx,j-miny,k-minz) = pvm.get_th(f+1);
	    phvol[f](i-minx,j-miny,k-minz) = pvm.get_ph(f+1);
	  }
	}
	
	for(int f=0;f<opts.nfibres.value();f++){
	  dyads[f](i-minx,j-miny,k-minz,0) = sin(thvol[f](i-minx,j-miny,k-minz)) * cos(phvol[f](i-minx,j-miny,k-minz));
	  dyads[f](i-minx,j-miny,k-minz,1) = sin(thvol[f](i-minx,j-miny,k-minz)) * sin(phvol[f](i-minx,j-miny,k-minz));
	  dyads[f](i-minx,j-miny,k-minz,2) = cos(thvol[f](i-minx,j-miny,k-minz));
	}
      }
    }
  }
  
  
  if(opts.verbose.value())
    cout << "saving results" << endl;

  S0.setDisplayMaximumMinimum(S0.max(),0);
  save_volume(S0,opts.ofile.value()+"_S0");

  dvol.setDisplayMaximumMinimum(dvol.max(),0);
  save_volume(dvol,opts.ofile.value()+"_D");

  if(opts.modelnum.value()==2){
    dvol_std.setDisplayMaximumMinimum(dvol_std.max(),0);
    save_volume(dvol_std,opts.ofile.value()+"_D_STD");
  }

  for(int f=1;f<=opts.nfibres.value();f++){
    fvol[f-1].setDisplayMaximumMinimum(1,0);
    save_volume(fvol[f-1],opts.ofile.value()+"_f"+num2str(f));
    thvol[f-1].setDisplayMaximumMinimum(thvol[f-1].max(),thvol[f-1].min());
    save_volume(thvol[f-1],opts.ofile.value()+"_th"+num2str(f));
    phvol[f-1].setDisplayMaximumMinimum(phvol[f-1].max(),phvol[f-1].min());
    save_volume(phvol[f-1],opts.ofile.value()+"_ph"+num2str(f));
    dyads[f-1].setDisplayMaximumMinimum(-1,1);
    save_volume4D(dyads[f-1],opts.ofile.value()+"_dyads"+num2str(f));

  }

  return 0;
}













