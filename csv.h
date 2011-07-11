/*  Combined Surfaces and Volumes Class

    Saad Jbabdi  - FMRIB Image Analysis Group
 
    Copyright (C) 2010 University of Oxford  */

/*  CCOPYRIGHT  */


#if !defined (CSV_H)
#define CSV_H

#ifndef MAX
#define MAX(x,y) ((x)>(y)?(x):(y))
#endif

#define VOLUME  0
#define SURFACE 1

#ifndef EXPOSE_TREACHEROUS
#define EXPOSE_TREACHEROUS
#endif

#include "newimage/newimageio.h"
#include "meshclass/meshclass.h"
#include "miscmaths/miscmaths.h"
#include "fslvtkio/fslvtkio.h"
#include "utils/tracer_plus.h"

using namespace NEWIMAGE;
using namespace MISCMATHS;
using namespace mesh;
using namespace fslvtkio;
using namespace Utilities;

// useful routines for collision detection
bool triBoxOverlap(float boxcenter[3],float boxhalfsize[3],float triverts[3][3]);
bool rayBoxIntersection(float origin[3],float direction[3],float vminmax[2][3]);
bool segTriangleIntersection(float seg[2][3],float tri[3][3]);
float triDistPoint(const Triangle& t,const ColumnVector pos);

// This is the main class that can handle an ROI made of a bunch of voxels and surface vertices
// Its main role is to tell you whether a point is in the ROI (voxels) or whether a line segment
// crossed its surface
// This class is quite tailored to fit with probtrackx requirements :)
// For now: it can only detect collisions with triangle meshes
// For now: it uses meshclass for I/O. Eventually port into Brian's I/O to support GIFTI
class CSV {
private:
  volume<float>                          roivol;     // volume-like ROIs (all in the same volume for memory efficiency)
  vector<Mesh>                           roimesh;    // surface-like ROIs

  vector<Mesh>                           lowresmesh; // low res surfaces (useful for counting)
  vector<Matrix>                         high2lowres;// correspondence between high and low res surfaces

  vector<string>                         surfnames;  // file names
  vector<string>                         volnames;   // useful for seed_to_target storage in probtrackx

  Matrix                                 hitvol;     // hit counter (voxels)
  Matrix                                 isinroi;    // which voxel in which ROI 
  Matrix                                 mat2vol;    // lookup volume coordinates
  volume<int>                            vol2mat;    // lookup matrix row
  volume4D<int>                          vol2bigmat; // lookup in big matrix (with meshes etc)


  volume<int>                            surfvol;    // this is the main thing that speeds up surface collision detection
                                                     // where each voxel knows whether the surface goes through or not
                                                     // the values inicate which triangle list crosses the voxel
                                                     // this volume can be higher resolution than refvol
  vector< vector< pair<int,int> > >      triangles;  // triangles crossed by voxels as pair<mesh,triangle>

  Matrix                                 mm2vox;     // for surfaces, transform coords into voxel space
                                                     // this will depend on various software conventions (argh)
  Matrix                                 vox2mm;     // the inverse transform
  string                                 convention; // FREESURFER/FIRST/CARET/etc.?
  

  float                                  _xdim;      // voxel sizes
  float                                  _ydim;
  float                                  _zdim;
  ColumnVector                           _dims;

  Matrix                                 _identity;  // useful for highres<->lowres

  int                                    nvols;      // # volume-like ROIs
  int                                    nsurfs;     // # surface-like ROIs
  int                                    nvoxels;    // # voxels (total across ROIs - unrepeated)
  int                                    nlocs;      // # ROI locations (eventually repeated)
  int                                    nrois;      // # ROIs in total (volumes+surfaces)

  vector< vector<int> >                  mesh2loc;   // mesh->locs lookup
  vector< pair<int,int> >                loc2mesh;   // locs->mesh lookup
  vector<ColumnVector>                   loccoords;
  vector<int>                            locrois;
  vector<int>                            locindex;
  vector<int>                            volind;
  vector<int>                            surfind;

  volume<short int>                      refvol;     // reference volume (in case no volume-like ROI is used)


public:
  CSV(){
    init_dims();
  }
  CSV(const volume<short int>& ref):refvol(ref){
    init_dims();
    _xdim=refvol.xdim();
    _ydim=refvol.ydim();
    _zdim=refvol.zdim();
    _dims.ReSize(3);_dims<<_xdim<<_ydim<<_zdim;
    set_convention("freesurfer");
    roivol.reinitialize(refvol.xsize(),refvol.ysize(),refvol.zsize());
    copybasicproperties(refvol,roivol);
    roivol=0;
    vol2mat.reinitialize(refvol.xsize(),refvol.ysize(),refvol.zsize());
    vol2mat=0;
    copybasicproperties(refvol,surfvol);
    surfvol.reinitialize(refvol.xsize(),
			 refvol.ysize(),
			 refvol.zsize());
    surfvol=0;
  }
  CSV(const CSV& csv){*this=csv;}
  ~CSV(){}

  void init_dims(){
    nvols=0;
    nsurfs=0;
    nvoxels=0;
    nlocs=0;
    nrois=0;
    _identity=IdentityMatrix(4);
  }

  // get/set
  const volume<short int>& get_refvol()const{return refvol;}
  void set_refvol(const volume<short int> vol){
    refvol=vol;
    _xdim=refvol.xdim();
    _ydim=refvol.ydim();
    _zdim=refvol.zdim();
    roivol.reinitialize(refvol.xsize(),refvol.ysize(),refvol.zsize());
    copybasicproperties(refvol,roivol);
    roivol=0;
    vol2mat.reinitialize(refvol.xsize(),refvol.ysize(),refvol.zsize());
    vol2mat=0;    

  }
  void change_refvol(const vector<string>& fnames){
    for(unsigned int i=0;i<fnames.size();i++){
      if(fsl_imageexists(fnames[i])){
	volume<short int> tmpvol;
	read_volume(tmpvol,fnames[i]);
	set_refvol(tmpvol);
	return;
      }
    }
  }
  void init_vol2bigmat(const vector<string>& fnames){
    int cnt=0;
    for(unsigned int i=0;i<fnames.size();i++){
      if(fsl_imageexists(fnames[i])){
	cnt++;
      }
    }
    vol2bigmat.reinitialize(vol2mat.xsize(),vol2mat.ysize(),vol2mat.zsize(),cnt);
    vol2bigmat=0;
  }
  inline const float xdim()const{return _xdim;}
  inline const float ydim()const{return _ydim;}
  inline const float zdim()const{return _zdim;}
  inline const int   xsize()const{return refvol.xsize();}
  inline const int   ysize()const{return refvol.ysize();}
  inline const int   zsize()const{return refvol.zsize();}

  int nRois ()const{return nrois;}
  int nSurfs()const{return nsurfs;}
  int nVols ()const{return nvols;}
  int nLocs ()const{return nlocs;}

  int nActVertices(const int& i)const{return (int)mesh2loc[i].size();}
  int nVertices(const int& i)const{return (int)roimesh[i]._points.size();}

  string get_name(const int& i){
    if(volind[i]>=0){return volnames[volind[i]];}
    else{return surfnames[surfind[i]];}
  }
  Mesh& get_mesh(const int i){return roimesh[i];}
  int get_loc(int roi,int x,int y,int z)const{
    return vol2bigmat(x,y,z,roi);
  }
  int get_loc(int roi,int vert)const{
    return mesh2loc[roi][vert];
  }
  bool isVol(const int& i)const{
    return (volind[i]>=0);
  }
  
  ColumnVector get_loc_coord(const int& ind){
    return loccoords[ind];
  }
  ColumnVector get_loc_coord_as_vox(const int& ind){
    if(locindex[ind]<0){
      return loccoords[ind];
    }
    else{
      return get_vertex_as_vox(loc2mesh[ind].first,locindex[ind]);
    }
  }

  // indices start at 0
  ColumnVector get_vertex(const int& surfind,const int& vertind){
    Mpoint *P;
    P=roimesh[surfind].get_point(vertind);
    ColumnVector p(3);
    p << P->get_coord().X
      << P->get_coord().Y
      << P->get_coord().Z;
    return p;
  }
  ColumnVector get_vertex_as_vox(const int& surfind,const int& vertind){
    Mpoint *P;
    P=roimesh[surfind].get_point(vertind);
    ColumnVector p(4);
    p << P->get_coord().X
      << P->get_coord().Y
      << P->get_coord().Z
      << 1;
    p = mm2vox*p;
    p = p.SubMatrix(1,3,1,1);
    return p;
  }
  ColumnVector get_normal(const int& surfind,const int& vertind){
    Vec n = roimesh[surfind].get_point(vertind)->local_normal();
    ColumnVector ret(3);
    ret<<n.X<<n.Y<<n.Z;
    return ret;
  }
  ColumnVector get_normal_as_vox(const int& surfind,const int& vertind){
    Vec n = roimesh[surfind].get_point(vertind)->local_normal();
    ColumnVector ret(3);
    ret<<n.X<<n.Y<<n.Z;
    ret=mm2vox.SubMatrix(1,3,1,3)*ret;
    
    return ret;
  }
  
  const vector<ColumnVector>& get_locs_coords()const{return loccoords;}
  const vector<int>& get_locs_roi_index()const{return locrois;}
  const vector<int>& get_locs_coord_index()const{return locindex;}
  int get_loc_roi(const int& i)const{return locrois[i];}

  // routines
  void load_volume  (const string& filename);
  void load_surface (const string& filename);
  void load_rois    (const string& filename,bool do_change_refvol=true);
  void reload_rois  (const string& filename);
  void save_roi     (const int& roiind,const string& prefix);
  void fill_volume  (volume<float>& vol,const int& ind);
  void cleanup();

  void set_convention(const string& conv);
  void switch_convention(const string& new_convention,const Matrix& vox2vox);

  void init_surfvol();
  void update_surfvol(const vector<ColumnVector>& v,const int& id,const int& meshid);
  void save_surfvol(const string& filename,const bool& binarise=true)const;
  void init_hitvol(const vector<string>& fnames);

  void  add_value(const string& type,const int& roiind,const int& locind,const float& val);
  void  set_value(const string& type,const int& roiind,const int& locind,const float& val);
  float get_value(const string& type,const int& roiind,const int& locind);
  void  reset_values();
  void  save_values(const int& roi);

  bool isInRoi(int x,int y,int z)const{return (roivol(x,y,z)!=0);}
  bool isInRoi(int x,int y,int z,int roi)const{
    if(!isInRoi(x,y,z))return false;
    return (isinroi(vol2mat(x,y,z),roi));
  }

  bool has_crossed(const ColumnVector& x1,const ColumnVector& x2,
		   bool docount=false,bool docontinue=false,const float& val=1.0);
  bool has_crossed_roi(const ColumnVector& x1,const ColumnVector& x2,
		       vector<int>* crossedrois=NULL,vector<int>* crossedlocs=NULL)const;


  // "near-collision" detection
  // returns true if one of the surfaces (meshind) is closer than dist. 
  // also returns (in dir) the direction towards the nearest surface vertex that is within dist
  // as well as the index for that closest surface
  bool is_near_surface(const ColumnVector& pos,const float& dist,ColumnVector& dir,float& mindist);

  void find_crossed_voxels(const ColumnVector& x1,const ColumnVector& x2,vector<ColumnVector>& crossed);
  void tri_crossed_voxels(float tri[3][3],vector<ColumnVector>& crossed);
  void line_crossed_voxels(float line[2][3],vector<ColumnVector>& crossed)const;


  // COPY
  CSV& operator=(const CSV& rhs){
    roivol=rhs.roivol;
    roimesh=rhs.roimesh;
    surfnames=rhs.surfnames;
    volnames=rhs.volnames;
    hitvol=rhs.hitvol;
    isinroi=rhs.isinroi;
    mat2vol=rhs.mat2vol;
    vol2mat=rhs.vol2mat;
    vol2bigmat=rhs.vol2bigmat;
    surfvol=rhs.surfvol;
    triangles=rhs.triangles;
    mm2vox=rhs.mm2vox;
    vox2mm=rhs.vox2mm;
    convention=rhs.convention;
    _xdim=rhs._xdim;
    _ydim=rhs._ydim;
    _zdim=rhs._zdim;
    _dims=rhs._dims;
    _identity=rhs._identity;
    nvols=rhs.nvols;
    nsurfs=rhs.nsurfs;
    nvoxels=rhs.nvoxels;
    nlocs=rhs.nlocs;
    nrois=rhs.nrois;
    mesh2loc=rhs.mesh2loc;
    loc2mesh=rhs.loc2mesh;
    loccoords=rhs.loccoords;
    locrois=rhs.locrois;
    locindex=rhs.locindex;
    volind=rhs.volind;
    surfind=rhs.surfind;
    refvol=rhs.refvol;
    return *this;
  }

};

#endif
