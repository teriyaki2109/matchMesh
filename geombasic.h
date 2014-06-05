#ifndef __GEOMBASIC_H_
#define __GEOMBASIC_H_

#include <stdio.h>
#include <math.h>
#include "list.h"
#include <list>
#include <vector>


class attribute
{
  public:
  std::vector<float> value;
  int uniqueness;
  
  attribute(){value.reserve(300); uniqueness=0; };
  ~attribute(){value.clear();};
};

class point_3
{
  float value; // valeur de la fonction en ce point utilis'ee pour le graphe de Reeb

  /* Pour la fonction computeMuValues */
  float distanceSum; // somme de toutes les distances calcul'ees jusqu'ici
                     // permet de choisir la nouvelle base 
  bool infiniteValue, // vrai si value=oo
    base; // vrai si ce point a d'eja` 'et'e choisi comme base
  point_3* colored; // =NULL si ce point n'est pas dans la ``zone color'ee'' (g(v)<r, cf article) d'une base
                   // sinon, 'egal a` la base pivot de cette zone color'ee
                   // quand tous les points sont color'es, l'algorithme de Dikstra est fini
  /* Fin de la fonction computeMuValues */ 

  void pointInit(void) ;

public:
  void pointDeInit();
  float* distanceTab; // si ce point est choisi comme base, ce tableau sera affecte et contiendra les distances a partir de ce point

  class attribute field;

  std::list<class edge *> edges;
  float x, y, z;
  
  float r, g ,b;

  //ttung: courbure gaussienne et moyenne au point
  float curv_g,curv_m,curv_i,val1,val2,mu;
  

  ~point_3() { pointDeInit(); }

  point_3 (float a, float b, float c) :
    x(a), y(b), z(c), r(0.0), g(0.0), b(0.0) { pointInit(); }

  point_3 (float a, float b, float c, float d, float e, float f, float g) :
    x(a), y(b), z(c), r(d), g(e), b(f) { pointInit(); value=g; }

  point_3 (float a) :
    x(a), y(a), z(a), r(0.0), g(0.0), b(0.0) { pointInit(); }

  point_3 () :
    x(0.0f), y(0.0f), z(0.0f), r(0.0f), g(0.0f), b(0.0f) { pointInit(); }

  point_3 operator+(const point_3 &p)const 
  {
    return point_3(x + p.x, y + p.y, z + p.z);
  }
  point_3 operator-(const point_3 &p)const
  {
    return point_3(x - p.x, y - p.y, z - p.z);
  }
  point_3 operator*(float d)const
  {
    return point_3(d*x,d*y,d*z);
  }
  point_3 operator/(float d)const
  {
    return point_3(x/d,y/d,z/d);
  }
 
  float operator*(const point_3& p) const
  {
    return (x * p.x + y * p.y + z * p.z);
  }


  /*
   * Pour les operateurs >,<.==, on fait le test sur le champ value
   * Pour cette raison, on fait tous les calculs intermediaires de l'algorithme de Dijkstra dans value
   * (car le BST utilise ces operatuers)
   */
  bool operator<(point_3 &p) const {
    return value < p.getValue();
  }
  bool operator>(point_3 &p) const {
    return value > p.getValue();
  }
  bool operator==(point_3 &p) const {
    return value == p.getValue();
  }

  bool operator==(const point_3 p) const {
    if((x==p.x) && (y==p.y) && (z==p.z))
		return true;
	else return false;
  }

  bool operator!=(const point_3 p) const {
    if((x!=p.x) || (y!=p.y) || (z!=p.z))
		return true;
	else return false;
  }


  void operator/=(float c) { x /= c; y /= c; z /= c; };
  void operator*=(float c) { x *= c; y *= c; z *= c; };
  void operator+=(float c) { x += c; y += c; z += c; };

  void operator+=(point_3 p) { x += p.x; y += p.y; z += p.z; };
  void operator-=(point_3 p) { x -= p.x; y -= p.y; z -= p.z; };
  void operator/=(point_3 p) { x /= p.x; y /= p.y; z /= p.z; };
  void operator*=(point_3 p) { x *= p.x; y *= p.y; z *= p.z; };


  /*
   * methodes
   */
  void setValue(float f) { value = f; }
  float getValue(void) const { return value; }
  float getDistanceSum(void) const { return distanceSum; }
  void addDistanceToSum(float dist) { distanceSum += dist; }
  void markInfiniteValue(void) { infiniteValue = true; }
  void markFiniteValue(void) { infiniteValue = false; }
  bool isInfinite(void) const { return infiniteValue; }
  void markBase(void) { base = true; }
  bool isBase(void) { return base; }
  void setColored(point_3* const &base) { colored = base; }
  point_3* getColored(void) const { return colored; }
  bool isColored(void) const { return colored != NULL; }
  void createDistanceTab(int size);
  void setDistanceTo(int n,float dist) { distanceTab[n] = dist; }
  float getDistanceTo(int n) const { return distanceTab[n]; }
  std::list<class edge *> getEdges(void) const { return edges; }
  void setEdges(std::list<class edge *> const &_edges) { edges = _edges; }
  void addEdgeToList(class edge * const &e) { edges.push_front(e); }
  class edge * findEdgeTo(point_3 * const &p);

  point_3 vec(const point_3 &p) const
  {
    return (point_3(y * p.z - z * p.y,z * p.x - x * p.z,x * p.y - y * p.x));
  }
  float norm() const 
  {
    return (float)sqrt(x*x + y*y +z*z);
  }

  point_3& normalize();

  point_3 cartesian() const;    
  point_3 spherical() const;

  //tung
  void setPoint_3(float a, float b, float c) { x = a; y = b; z = c; }
};
typedef class point_3 Point_3;


class face 
{
  /* Pour la fonction subdivide */
  bool obsolete; //une face est obsolete si elle a 'et'e subidivs'ee en plusieurs faces
  float lastLevel; //dernier niveau de la fonction Mu pour lequel la face a 'et'e trait'ee 
  /* Fin de la fonction subdivide */

  public :
  
  class point_3* p[3];
  class edge* e[3];
  
  //int i[3];       //indice des points
  short int FaceNumber; //numero de la face dans l'objet 3d

  face(void);
  face(Point_3* const &p0,Point_3* const &p1,Point_3*const & p2,class edge* const &e0,class edge* const &e1,class edge* const &e2);

  void createPoints(void);

  bool isObsolete(void) const { return obsolete; }
  void markObsolete(void) { obsolete = true; }
  float getLastLevel(void) const { return lastLevel; }
  void setLastLevel(float newLevel) { lastLevel = newLevel; }
  float area(void) const;
  point_3 center(void) const;
  int nbBaseColored(Point_3* const &base) ;
  int edgeIndex(class edge* const &e);

  // ttung
  void setFaceNumber(int FaceNb) { FaceNumber = FaceNb; }
  int  getFaceNumber(void) const { return FaceNumber; }
  float angle(unsigned int f) const;
  point_3 normal(void) const;
  void setFacePoints(float x0,float y0,float z0,float x1,float y1,float z1,float x2,float y2,float z2);

};
typedef class face Face;



class edge 
{
  public :

  class point_3 *p1,*p2;
  class face *f1,*f2;

  /*Pour la fonction subdivide*/
  class edge* newEdge1,*newEdge2;
  class point_3* newPoint;
  bool obsolete;
  /*Fin de la fonction subdivide*/

  edge (void) :
    p1(NULL),p2(NULL),
    f1(NULL),f2(NULL),
    obsolete(false)
  {}

  edge (Point_3 * const &_p1, Point_3 * const &_p2) :
    p1(_p1), p2(_p2),
    f1(NULL), f2(NULL),
    obsolete(false)
  {}

  float length() const { return ((*p1) - (*p2)).norm(); }
  Point_3* otherPoint(Point_3* const &p) const;
  void setFaces(Face * const &_f1,Face * const &_f2) { f1 = _f1; f2 = _f2; }
  Face* otherFace(Face * const &f) const ;

  void replaceFace(Face* const &old,Face* const &_new);
  void markObsolete(class edge* const &_newEdge1,class edge* const &_newEdge2,Point_3* const &_newPoint);

  bool isObsolete(void) const { return obsolete; }
  class edge* otherSubEdge(class edge* const &subEdge) const;
  bool isIsoLine(float value) const;
};
typedef class edge Edge;


#endif




