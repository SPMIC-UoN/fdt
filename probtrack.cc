#include <iostream>
#include <fstream>
#include "newimage/newimageall.h"
#include "utils/log.h"
#include "meshclass/meshclass.h"
#include "probtrack.h"


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
  else if(opts.mode.value()=="matrix2"){
    if(opts.meshfile.value()=="")
      matrix2();
    else
      mesh_matrix2();
  }
  else if(opts.mode.value()=="maskmatrix")
    maskmatrix();
  else if(opts.mode.value()=="meshlengths")
    mesh_lengths();
  else{
    cout <<"Invalid setting for option  mode -- try setting mode=help"<<endl;
  }
      
  return 0;
}















