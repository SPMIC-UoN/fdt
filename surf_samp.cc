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


Option<string> dir(string("--dir"),string(""),
		   string("bpx directory"),
		   true,requires_argument);
Option<string> surf(string("--surf"),string(""),
		    string("surface file"),
		    true,requires_argument);
Option<string> surfvol(string("--meshref"),string(""),
		    string("surface volume ref"),
		    true,requires_argument);
Option<string> xfm(string("--xfm"),string(""),
		   string("surf2diff transform"),
		   true,requires_argument);
Option<string> meshspace(string("--meshspace"),string("freesurfer"),
		   string("meshspace (default='freesurfer')"),
		   true,requires_argument);
Option<string> out(string("--out"),string(""),
		   string("output file"),
		   true,requires_argument);
Option<float> dist(string("--step"),1,
		   string("step (mm - default=1)"),
		   false,requires_argument);
Option<float> thr(string("--thr"),0,
		  string("f-threshold (default=0)"),
		  false,requires_argument);


int main(int argc,char *argv[]){

  OptionParser options(title,examples);
  
  
  options.add(dir);
  options.add(surf);
  options.add(surfvol);
  options.add(xfm);
  options.add(meshspace);
  options.add(out);
  options.add(dist);
  options.add(thr);
  
  options.parse_command_line(argc,argv);
  if(!options.check_compulsory_arguments(true)){
    options.usage();
    return(1);
  }
  
  cout<<"read"<<endl;
  volume<short int> refvol;
  read_volume(refvol,surfvol.value());

  CSV csv(refvol);
  csv.set_convention(meshspace.value());
  csv.load_rois(surf.value());

  vector< volume4D<float> > dyads;
  vector< volume<float> > fs;
  volume4D<float> tmpdyad;volume<float> tmpf;
  vector< CSV* > osurfs;CSV* tmpsurf;
  int f=1;
  while( fsl_imageexists(dir.value()+"/dyads"+num2str(f)) ){
    read_volume4D(tmpdyad,dir.value()+"/dyads"+num2str(f));
    dyads.push_back(tmpdyad);
    read_volume(tmpf,dir.value()+"/mean_f"+num2str(f)+"samples");
    fs.push_back(tmpf);
    tmpsurf = new CSV(csv);
    tmpsurf->reset_values();
    osurfs.push_back(tmpsurf);
    f++;
  }  

  Matrix surf2dti = read_ascii_matrix(xfm.value());
  Matrix ixfm(3,3);
  ixfm=surf2dti.SubMatrix(1,3,1,3);
  ixfm=ixfm.i();

  ColumnVector surfdim(3),dtidim(3);
  surfdim << refvol.xdim() << refvol.ydim() << refvol.zdim();
  dtidim << dyads[0].xdim()<<dyads[0].ydim()<<dyads[0].zdim();

  float step = dist.value(); // mm

  cout<<"loop"<<endl;
  // loop over surface vertices
  // for each vertex, go along normal for a fixed distance, and sample orientations
  // save the orientation, the f1/f2, etc.
  ColumnVector x(3),n(3),dti_vox(3),nn(3);
  int ix,iy,iz;
  for(int i=0;i<csv.get_mesh(0).nvertices();i++){
    x=csv.get_vertex_as_vox(0,i);
    n=csv.get_normal_as_vox(0,i);

    // flip normal in caret
    if(meshspace.value()=="caret")
      n*=-1;
    if(n.MaximumAbsoluteValue()>0)
      n/=sqrt(n.SumSquare());
    
    x<<x(1)+step/surfdim(1)*n(1)
     <<x(2)+step/surfdim(2)*n(2)
     <<x(3)+step/surfdim(3)*n(3);

    dti_vox = vox_to_vox(x,surfdim,dtidim,surf2dti);
    ix = (int)round((float)dti_vox(1));
    iy = (int)round((float)dti_vox(2));
    iz = (int)round((float)dti_vox(3));


    // dyad*n
    ColumnVector d(3);
    for(unsigned int f=0;f<dyads.size();f++){
      if(fs[f](ix,iy,iz)<thr.value())
	continue;
      d << dyads[f](ix,iy,iz,0)
	<< dyads[f](ix,iy,iz,1)
	<< dyads[f](ix,iy,iz,2);
      d = ixfm*d;
      if(d.MaximumAbsoluteValue()>0)
	d/=sqrt(d.SumSquare());

      float dp = abs(dot(n,d));
      osurfs[f]->set_value("surface",0,i,dp);
    }
    
  }
  cout<<"save"<<endl;

  for(unsigned int f=0;f<dyads.size();f++){
    osurfs[f]->save_roi(0,out.value()+"_dyads"+num2str(f+1));
  }
  


  return 0;
}



