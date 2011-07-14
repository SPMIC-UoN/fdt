#include "utils/options.h"
#include "newimage/newimageall.h"
#include "csv.h"
#include "stdlib.h"
#include "string.h"
#include "miscmaths/miscmaths.h"



using namespace Utilities;
using namespace std;
using namespace NEWIMAGE;
using namespace MISCMATHS;


string title="";
string examples="";


Option<string> surfin(string("--surfin"),string(""),
		      string("input surface(s) - ascii list"),
		      true,requires_argument);
Option<string> surfout(string("--surfout"),string(""),
		       string("output surface(s) - ascii list"),
		       true,requires_argument);
Option<string> convin(string("--convin"),string(""),
		      string("input convention"),
		      true,requires_argument);
Option<string> convout(string("--convout"),string(""),
		       string("output convention [default=same as input]"),
		       false,requires_argument);
Option<string> volin(string("--volin"),string(""),
		     string("input volume"),
		     true,requires_argument);
Option<string> volout(string("--volout"),string(""),
		      string("output volume [default=same as input]"),
		      false,requires_argument);
Option<string> xfm(string("--xfm"),"",
		   string("in-to-out transformation (ascii matrix) [default=identity]"),
		   false,requires_argument);


int main(int argc,char *argv[]){

  OptionParser options(title,examples);
  
  
  options.add(surfin);
  options.add(surfout);
  options.add(convin);
  options.add(convout);
  options.add(volin);
  options.add(volout);
  options.add(xfm);

  
  options.parse_command_line(argc,argv);
  if(!options.check_compulsory_arguments(true)){
    options.usage();
    return(1);
  }
  

  volume<short int> refvolin,refvolout;
  read_volume(refvolin,volin.value());
  if(volout.set()){
    read_volume(refvolout,volout.value());
  }
  else{
    refvolout.reinitialize(refvolin.xsize(),refvolin.ysize(),refvolin.zsize());
    refvolout=refvolin;
  }

  CSV csv(refvolin);
  csv.set_convention(convin.value());
  csv.load_rois(surfin.value());

  Matrix vox2vox;
  if(xfm.set())
    vox2vox=read_ascii_matrix(xfm.value());
  else
    vox2vox=IdentityMatrix(4);

  csv.set_refvol(refvolout);
  if(convout.set())
    csv.switch_convention(convout.value(),vox2vox);
  else
    csv.switch_convention(convin.value(),vox2vox);

  vector<string> filenames;
  ifstream fs(surfout.value().c_str());
  string tmp;
  if(fs){
    fs>>tmp;
    do{
      filenames.push_back(tmp);
      fs>>tmp;
    }while(!fs.eof());
  }
  else{
    cerr<<surfout.value()<<" does not exist"<<endl;
    exit(1);
  }
  
  for(int i=0;i<csv.nSurfs();i++){
    csv.save_roi(i,filenames[i]);
  }

  return 0;
}



