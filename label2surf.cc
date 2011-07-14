/*  label2surf.cc

    Transform a set of label files into a surface

    Saad Jbabdi, FMRIB Image Analysis Group

    Copyright (C) 2010 University of Oxford */

/*  CCOPYRIGHT */


#include "meshclass/meshclass.h"
#include "utils/options.h"

using namespace Utilities;
using namespace mesh;



string title="label2surf \n\t Transform a group of labels into a surface";
string examples="label2surf -s <surface> -o <outputsurface> -l <labels>";

Option<bool> verbose(string("-v,--verbose"),false,
		       string("switch on diagnostic messages"),
		       false,no_argument);
Option<bool> help(string("-h,--help"),false,
		       string("display this message"),
		       false,no_argument);
Option<string> isurf(string("-s,--surf"),"",
		     string("input surface (freesurfer ascii)"),
		     true,requires_argument);
Option<string> osurf(string("-o,--out"),"",
		     string("output surface (freesurfer ascii)"),
		    true,requires_argument);
Option<string> labels(string("-l,--labels"),"",
		      string("ascii list of label files"),
		      true,requires_argument);


void read_fnames(vector<string>& fnames,const string& filename){
  fnames.clear();
  ifstream fs(filename.c_str());
  string tmp;
  if(fs){
    fs>>tmp;
    do{
      fnames.push_back(tmp);
      fs>>tmp;
    }while(!fs.eof());
  }
  else{
    cerr<<filename<<" does not exist"<<endl;
    exit(0);
  }
}


int main(int argc,char *argv[]){
  OptionParser options(title,examples);
 
  options.add(verbose);
  options.add(help);
  options.add(isurf);
  options.add(osurf);
  options.add(labels);

  options.parse_command_line(argc,argv);
  if ( (help.value()) || (!options.check_compulsory_arguments(true)) ){
    options.usage();
    return 1;
  }
  
  ////////
  if(verbose.value())
    cout<<"read input surface"<<endl;
  
  Mesh m;
  m.load(isurf.value());
  
  if(verbose.value())
    cout<<"read input labels"<<endl;
  
  vector<string> labs;
  read_fnames(labs,labels.value());
  for(unsigned int i=0;i<labs.size();i++)
    m.load_fs_label(labs[i]);

  m.save_fs(osurf.value());

  return 0;
}
