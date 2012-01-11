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

#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include "newimage/newimageio.h"
#include "miscmaths/miscmaths.h"
#include "fslvtkio/fslvtkio.h"
#include "utils/tracer_plus.h"
#include "csv_mesh.h"
//#include "fslsurface/fslsurface.h"
//#include "fslsurface/fslsurfaceio.h"
//#include "fslsurface/fslsurface_dataconv.h"

using namespace std;
using namespace NEWIMAGE;
using namespace MISCMATHS;
using namespace fslvtkio;
using namespace Utilities;
//using namespace fslsurface_name;

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
  vector<CsvMesh>                        roimesh;    // surface-like ROIs

  vector<string>                         roinames;   // file names
                                                     // useful for seed_to_target storage in probtrackx

  Matrix                                 hitvol;     // hit counter (voxels)
  Matrix                                 isinroi;    // which voxel in which ROI 
  Matrix                                 mat2vol;    // lookup volume coordinates
  volume<int>                            vol2mat;    // lookup matrix row

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

  // Sort out the transformations
  // from a given loc, I need to be able to tell 
  // whether it is a surface or volume
  // + which surface or volume
  // + which sublocation within the surface or volume
  vector<int>           loctype;            // SURFACE OR VOLUME
  vector<int>           locroi;             // 0:NROIS-1
  vector<int>           locsubroi;          // 0:NSURFS or 0:NVOLS
  vector<int>           locsubloc;          // 0:n
  vector<ColumnVector>  loccoords;
  vector< vector<int> > mesh2loc;
  volume4D<int>         vol2loc;
  vector<int>           volind;             // vol is which roi index?
  vector<int>           surfind;            // surf is which roi index?
  vector<int>           roitype;
  vector<int>           roisubind;          // roi is which subindex?

  volume<short int>     refvol;     // reference volume (in case no volume-like ROI is used)


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
  ~CSV(){clear_all();}

  void init_dims(){
    nvols=0;
    nsurfs=0;
    nvoxels=0;
    nlocs=0;
    nrois=0;
    _identity=IdentityMatrix(4);
  }
  void reinitialize(const volume<short int>& vol){
    init_dims();
    _xdim=vol.xdim();
    _ydim=vol.ydim();
    _zdim=vol.zdim();
    _dims.ReSize(3);_dims<<_xdim<<_ydim<<_zdim;
    set_convention("freesurfer");
    set_refvol(vol);
    clear_all();
  }
  
  void clear_all(){
    roimesh.clear();roinames.clear();
    loctype.clear();locroi.clear();locsubroi.clear();locsubloc.clear();
    mesh2loc.clear();
  }

  // get/set
  const volume<short int>& get_refvol()const{return refvol;}
  void set_refvol(const volume<short int>& vol){
    refvol=vol;
    _xdim=refvol.xdim();
    _ydim=refvol.ydim();
    _zdim=refvol.zdim();
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
  void init_vol2loc(const vector<string>& fnames){
    int cnt=0;
    for(unsigned int i=0;i<fnames.size();i++){
      if(fsl_imageexists(fnames[i])){
	cnt++;
      }
    }
    vol2loc.reinitialize(vol2mat.xsize(),vol2mat.ysize(),vol2mat.zsize(),cnt);
    vol2loc=0;
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
    return roinames[i];
  }
  int get_roitype(const int& i){return roitype[i];}
  CsvMesh& get_mesh(const int i){return roimesh[i];}
  int get_volloc(int subroi,int x,int y,int z)const{
    return vol2loc(x,y,z,subroi);
  }
  int get_surfloc(int subroi,int subloc)const{
    return mesh2loc[subroi][subloc];
  }
  
  ColumnVector get_loc_coord(const int& loc){
    return loccoords[loc];
  }
  ColumnVector get_loc_coord_as_vox(const int& loc){
    if(loctype[loc]==VOLUME){
      return loccoords[loc];
    }
    else{
      return get_vertex_as_vox(locsubroi[loc],locsubloc[loc]);
    }
  }
  vector<ColumnVector> get_locs_coords()const{return loccoords;}
  vector<int> get_locs_roi_index()const{
    vector<int> ret;
    for(unsigned int i=0;i<locroi.size();i++){
      ret.push_back(locroi[i]);
    }
    return ret;
  }
  vector<int> get_locs_coord_index()const{
    vector<int> ret;
    for(unsigned int i=0;i<locsubloc.size();i++){
      ret.push_back(locsubloc[i]);
    }
    return ret;
  }
  

  // indices start at 0
  ColumnVector get_vertex(const int& surfind,const int& vertind){
    CsvMpoint P(roimesh[surfind].get_point(vertind));
    ColumnVector p(3);
    p << P.get_coord().X
      << P.get_coord().Y
      << P.get_coord().Z;
    return p;
  }
  ColumnVector get_vertex_as_vox(const int& surfind,const int& vertind){
    CsvMpoint P(roimesh[surfind].get_point(vertind));
    ColumnVector p(4);
    p << P.get_coord().X
      << P.get_coord().Y
      << P.get_coord().Z
      << 1;
    p = mm2vox*p;
    p = p.SubMatrix(1,3,1,1);
    return p;
  }
   ColumnVector get_normal(const int& surfind,const int& vertind){
     Vec n = roimesh[surfind].local_normal(vertind);
     ColumnVector ret(3);
     ret<<n.X<<n.Y<<n.Z;
     return ret;
   }
   ColumnVector get_normal_as_vox(const int& surfind,const int& vertind){
     Vec n = roimesh[surfind].local_normal(vertind);
     ColumnVector ret(3);
     ret<<n.X<<n.Y<<n.Z;
     ret=mm2vox.SubMatrix(1,3,1,3)*ret;
     if(ret.MaximumAbsoluteValue()>0)
       ret/=std::sqrt(ret.SumSquare());
     return ret;
   }
  

  // routines
  void load_volume  (const string& filename);
  void load_surface (const string& filename);
  void load_rois    (const string& filename,bool do_change_refvol=true);
  void reload_rois  (const string& filename);
  void save_roi     (const int& roiind,const string& prefix);
  void save_rois    (const string& prefix);
  void fill_volume  (volume<float>& vol,const int& ind);
  void cleanup();

  void set_convention(const string& conv);
  void switch_convention(const string& new_convention,const Matrix& vox2vox);

  void init_surfvol();
  void update_surfvol(const vector<ColumnVector>& v,const int& id,const int& meshid);
  void save_surfvol(const string& filename,const bool& binarise=true)const;
  void save_normalsAsVol(const int& meshind,const string& filename)const;
  void init_hitvol(const vector<string>& fnames);


  void  add_value(const int& loc,const float& val);
  void  set_value(const int& loc,const float& val);
  void  set_vol_values(const float& val);
  float get_value(const int& loc);
  void  reset_values();
  void  reset_values(const vector<int>& locs);
  void  save_values(const int& roi);
  void  loc_info(const int& loc)const;

  bool isInRoi(int x,int y,int z)const{return (roivol(x,y,z)!=0);}
  bool isInRoi(int x,int y,int z,int roi)const{
    if(!isInRoi(x,y,z))return false;
    return (isinroi(vol2mat(x,y,z),roi));
  }

  bool has_crossed(const ColumnVector& x1,const ColumnVector& x2,
		   bool docount=false,bool docontinue=false,const float& val=1.0);
  bool has_crossed_roi(const ColumnVector& x1,const ColumnVector& x2,
		       vector<int>& crossedrois)const;
  bool has_crossed_roi(const ColumnVector& x1,const ColumnVector& x2,
		       vector<int>& crossedrois,vector<int>& crossedlocs)const;


  int step_sign(const int& loc,const Vec& step)const;
  int coord_sign(const int& loc,const ColumnVector& x2)const;

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
    roinames=rhs.roinames;
    hitvol=rhs.hitvol;
    isinroi=rhs.isinroi;
    mat2vol=rhs.mat2vol;
    vol2mat=rhs.vol2mat;
    vol2loc=rhs.vol2loc;
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
    loctype=rhs.loctype;
    locroi=rhs.locroi;
    locsubroi=rhs.locsubroi;
    locsubloc=rhs.locsubloc;
    loccoords=rhs.loccoords;
    mesh2loc=rhs.mesh2loc;
    vol2loc=rhs.vol2loc;
    volind=rhs.volind;
    surfind=rhs.surfind;
    roitype=rhs.roitype;
    roisubind=rhs.roisubind;
    refvol=rhs.refvol;
    return *this;
  }

};

 
#endif
