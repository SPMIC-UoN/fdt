/*  Copyright (C) 2009 University of Oxford  */

/*  CCOPYRIGHT  */

#include <iostream>
#include <cmath>
#include "utils/options.h"
#include "miscmaths/miscmaths.h"
#include "newmat.h"
#include "newimage/newimageall.h"

using namespace std;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;



string title="dtigen - generate diffusion data using tensor model";
string examples="dtigen -t <input4Dtensor> -o <output4Ddata> -b <bvals> -r <bvecs> -m <brainmask> --s0=<s0file>";

Option<bool> help(string("-h,--help"),false,
		       string("display this message"),
		       false,no_argument);
Option<string> itensor(string("-t,--tensor"),string(""),
		       string("input tensor"),
		       true,requires_argument);
Option<string> s0file(string("--s0"),string(""),
		       string("input S0"),
		       true,requires_argument);
Option<string> odata(string("-o,--output"),string(""),
		       string("output data"),
		       true,requires_argument);
Option<string> bvecsfile(string("-r,--bvecs"),string(""),
		       string("bvecs ASCII text file"),
		       true,requires_argument);
Option<string> bvalsfile(string("-b,--bvals"),string(""),
		       string("bvals ASCII text file"),
		       true,requires_argument);
Option<string> maskfile(string("-m,--mask"),string(""),
		       string("brain mask"),
		       true,requires_argument);

int do_dtigen(){
  volume<float> mask,S0;
  volume4D<float> data,tensor;

  read_volume(mask,maskfile.value());
  read_volume(S0,s0file.value());
  read_volume4D(tensor,itensor.value());
  
  Matrix r = read_ascii_matrix(bvecsfile.value());
  if(r.Nrows()>3) r=r.t();
  for(int i=1;i<=r.Ncols();i++){
    float tmpsum=sqrt(r(1,i)*r(1,i)+r(2,i)*r(2,i)+r(3,i)*r(3,i));
    if(tmpsum!=0){
      r(1,i)=r(1,i)/tmpsum;
      r(2,i)=r(2,i)/tmpsum;
      r(3,i)=r(3,i)/tmpsum;
    }  
  }
  Matrix b = read_ascii_matrix(bvalsfile.value());
  if(b.Nrows()>1) b=b.t();
  if( b.Ncols() != r.Ncols() ){ cerr << "Error: bvecs and bvals don't have the same number of entries" << endl; return(-1);}
  if( r.Nrows() !=3 ){cerr << "Error: bvecs must be either 3xN or Nx3" << endl; return(-1);}
 
  data.reinitialize(mask.xsize(),mask.ysize(),mask.zsize(),b.Ncols());
  copybasicproperties(tensor[0],data);

  cout << "generate data" << endl << endl;;
  for(int z=mask.minz();z<=mask.maxz();z++){
    cout << "processing slice" << z << endl;
    for(int y=mask.miny();y<=mask.maxy();y++)
      for(int x=mask.minx();x<=mask.maxx();x++){
	if(mask(x,y,z)==0)continue;
	
	for(int t=1;t<=data.tsize();t++){
	  data(x,y,z,t-1) = S0(x,y,z);
	  data(x,y,z,t-1) *= exp( -b(1,t) * r(1,t) * r(1,t) * tensor(x,y,z,0) );
	  data(x,y,z,t-1) *= exp( -2*b(1,t) * r(1,t) * r(2,t) * tensor(x,y,z,1) );
	  data(x,y,z,t-1) *= exp( -2*b(1,t) * r(1,t) * r(3,t) * tensor(x,y,z,2) );
	  data(x,y,z,t-1) *= exp( -b(1,t) * r(2,t) * r(2,t) * tensor(x,y,z,3) );
	  data(x,y,z,t-1) *= exp( -2*b(1,t) * r(2,t) * r(3,t) * tensor(x,y,z,4) );
	  data(x,y,z,t-1) *= exp( -b(1,t) * r(3,t) * r(3,t) * tensor(x,y,z,5) );
	  
	}

      }
  }
  cout<<"saving results" << endl;
  data.setDisplayMaximumMinimum(1000,0);
  save_volume4D(data,odata.value());


  return 0;
}

int main(int argc,char *argv[]){

  Tracer tr("main");
  OptionParser options(title,examples);

  try{
    options.add(help);
    options.add(itensor);
    options.add(s0file);
    options.add(odata);
    options.add(bvecsfile);
    options.add(bvalsfile);
    options.add(maskfile);

    options.parse_command_line(argc,argv);

    
    if ( (help.value()) || (!options.check_compulsory_arguments(true)) ){
      options.usage();
      exit(EXIT_FAILURE);
    }
  }
  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } 
  catch(std::exception &e) {
    cerr << e.what() << endl;
  } 
  
  return do_dtigen();
  
  
}
