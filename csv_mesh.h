#if !defined (CSV_MESH_H)
#define CSV_MESH_H


#include <stdio.h>
#include "meshclass/point.h"
#include "miscmaths/miscmaths.h"

using namespace std;
using namespace MISCMATHS;
using namespace mesh;

int  meshFileType(const string& filename);
bool meshExists(const string& filename);

class CsvMpoint {

 public:

  CsvMpoint(double x, double y, double z, int counter):_no(counter){
    _coord=Pt(x, y, z);
  }
  CsvMpoint(const Pt p, int counter,float val=0):_no(counter){
    _coord=p;
  }
  ~CsvMpoint(){}
  CsvMpoint(const CsvMpoint &p):_coord(p._coord),_no(p._no){
    *this=p;
  }
  CsvMpoint& operator=(const CsvMpoint& p){
    _coord=p._coord;
    _no=p._no;
    return(*this);
  }

  const Pt&    get_coord() const         {return _coord;};
  void         set_coord(const Pt& coord){_coord=coord;}
  const int&   get_no() const            {return _no;};


 private:
  Pt        _coord;
  int       _no;
};

class CsvTriangle {
 private:
  vector<CsvMpoint> _vertice;
  int               _no;

 public:
  CsvTriangle(){}
  CsvTriangle(const CsvMpoint& p1,const CsvMpoint& p2,const CsvMpoint& p3,int no):_no(no){
    _vertice.push_back(p1);
    _vertice.push_back(p2);
    _vertice.push_back(p3);
  }
  ~CsvTriangle(){}
  CsvTriangle(const CsvTriangle &t):_vertice(t._vertice),_no(t._no){
    *this=t;
  }
  CsvTriangle& operator=(const CsvTriangle& t){
    _vertice=t._vertice;_no=t._no;return *this;
  }

  Pt  centroid() const{
    Pt p ((_vertice[0].get_coord().X +_vertice[1].get_coord().X +_vertice[2].get_coord().X)/3,
	  (_vertice[0].get_coord().Y +_vertice[1].get_coord().Y +_vertice[2].get_coord().Y)/3,
	  (_vertice[0].get_coord().Z +_vertice[1].get_coord().Z +_vertice[2].get_coord().Z)/3);

    return p;
  }
  Vec normal() const{
    Vec result = (_vertice[2].get_coord() - _vertice[0].get_coord()) * (_vertice[1].get_coord() - _vertice[0].get_coord());
    result.normalize();
    return result;     
  }

  const CsvMpoint& get_vertice(const int& i) const{return _vertice[i];}
  const bool intersect(const vector<Pt> & p)const;          // checks if a segment intersects the triangle
  const bool intersect(const vector<Pt> & p,int& ind)const; // checks if a segment intersects the triangle+gives the index of the closest vertex
  int get_no(){return _no;}
};


const bool operator ==(const CsvMpoint &p2, const CsvMpoint &p1);
const bool operator ==(const CsvMpoint &p2, const Pt &p1);

const Vec operator -(const CsvMpoint&p1, const CsvMpoint &p2);
const Vec operator -(const Pt&p1, const CsvMpoint &p2);
const Vec operator -(const CsvMpoint&p1, const Pt &p2);

class CsvMesh {

 public:
  vector<CsvMpoint>   _points;
  vector<CsvTriangle> _triangles;
  vector<float>       _pvalues;
  vector<float>       _tvalues;
  

  CsvMesh(){}
  ~CsvMesh(){}
  CsvMesh(const CsvMesh &m):_points(m._points),_triangles(m._triangles),_pvalues(m._pvalues),_tvalues(m._tvalues){
    *this=m;
  }
  CsvMesh& operator=(const CsvMesh& m){
    _points=m._points;
    _triangles=m._triangles;
    _pvalues=m._pvalues;
    _tvalues=m._tvalues;
    return(*this);
  }
  
  void        clear(){
    _points.clear();_triangles.clear();_pvalues.clear();_tvalues.clear();
  }

  int   nvertices() const{return (int)_points.size();}
  int   ntriangles() const{return (int)_triangles.size();}
  const CsvMpoint&    get_point(int n)const{return _points[n];}
  const CsvTriangle&  get_triangle(int n)const{return _triangles[n];}

  float get_pvalue(const int& i)const{return _pvalues[i];}
  float get_tvalue(const int& i)const{return _tvalues[i];}
  void set_pvalue(const int& i,const float& val){_pvalues[i]=val;}
  void set_tvalue(const int& i,const float& val){_tvalues[i]=val;}

  void reset_pvalues(){
    _pvalues.clear();
    for(unsigned int i=0;i<_points.size();i++)
      _pvalues.push_back(0);
  }
  void reset_tvalues(){
    _tvalues.clear();
    for(unsigned int i=0;i<_triangles.size();i++)
      _tvalues.push_back(0);  
  }

  void load(const string& filename); 
  void load_ascii(const string& filename);
  void load_vtk(const string& filename);

  void save_ascii(const string& s);


  void print(){
    cout<<_points.size()<<" vertices"<<endl;
    for(unsigned int i=0;i<_points.size();i++){
      cout<<_points[i].get_coord().X<<" "
	  <<_points[i].get_coord().Y<<" "
	  <<_points[i].get_coord().Z<<endl;
    }
    cout<<_triangles.size()<<" triangles"<<endl;
    for(unsigned int i=0;i<_triangles.size();i++){
      cout<<_triangles[i].get_vertice(0).get_no()<<" "
	  <<_triangles[i].get_vertice(1).get_no()<<" "
	  <<_triangles[i].get_vertice(2).get_no()<<endl;
    }
  }

  //  ostream& operator <<(ostream& flot);
};



#endif
