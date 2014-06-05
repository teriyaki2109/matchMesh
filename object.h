#ifndef __OBJECT_H_
#define __OBJECT_H_

#include <stdlib.h>
#include <math.h>
#include "bstindex.h"
#include "list.h"
#include "geombasic.h"
#include "mrgraph.h"
#include "dataItem.h"
//#include "statusBar.h"


class object 
{
public:
  float maxLength;
  void computeMaxLength(void);

  int pointArraySize; // taille d'un tableau de points
  int faceArraySize;
  int edgeArraySize;
  
  Point_3** pointArrays;
  int nbPointArrays; //nombre de tableaux de points
  int indexLastPoint; //index du dernier point dans le dernier tableau
  // on doit toujours avoir :
  // nbPoints=max(nbPointArrays-1,0)*arraySize+(indexLastPoint+1)

  Face** faceArrays;
  int nbFaceArrays;
  int indexLastFace;

  Edge** edgeArrays;
  int nbEdgeArrays;
  int indexLastEdge;

  public :
    
  Point_3 *t3d; //ttung:triangles de texture, il y en a NbTriangles et ils représentent les coordonnées dans la carte

  std::list<int> bases;
  float totalArea;
  float totalVolume;  // ttung: total volume
  
  bool muComputed; // indique si les valeurs (fonction mu) des points ont ete calcule

  bool acp;        // ttung: indique si on applique l'acp (true) ou non (false)
  char *name;
  bool texture;    // ttung: tell if there is a texture associated to the object

  void objectInit();
  object(int _pointArraySize,int _faceArraySize,int _edgeArraySize);  
  ~object();

  // Gestion des points
  int getNbPoints(void) const;
  Point_3* addNewPoint(const Point_3 &p);
  int getPointIndex(Point_3* const &p) const;
  Point_3* pointNb(int index) const;
  void resetPoints(void);
  void pointCoord(int index,int &tableau,int &offset) const;
  void setPointNb(int i,Point_3* p);
  void lessPoints(int howMany);
  
  // Gestion des faces
  int getNbFaces(void) const;
  Face* addNewFace(const Face &f);
  int getFaceIndex(Face* const &f) const;
  Face* faceNb(int index) const;
  int getNbNonObsoleteFaces(void) const;
  Face* nonObsoleteFaceNb(int index) const;
  Face** getAllNonObsoleteFaces(void) const;
  void resetFaces(void);
  void faceCoord(int index,int &tableau,int &offset) const;
  void setFaceNb(int i,Face* f);
  void lessFaces(int howMany);

  // Gestion des aretes
  int getNbEdges(void) const;
  Edge* addNewEdge(const Edge &e);
  int getEdgeIndex(Edge* const &e) const;
  Edge* edgeNb(int index) const;
  int getNbNonObsoleteEdges(void) const;
  void resetEdges(void);
  void edgeCoord(int index,int &tableau,int &offset) const;
  void setEdgeNb(int i,Edge* e);
  void lessEdges(int howMany);

  void setTotalArea(float area) { totalArea = area; }
  float getTotalArea(void) const { return totalArea; }
  float* getValueArray(int resolution) const;

  float* getValueArrayOnTri(int resolution) const;

  float* getGaussCurvArray(void) const; //ttung : get gaussian curvature attribute
  float* getMeanCurvArray(void) const;  //ttung : get mean curvature attribute
  float* getIndexCurvArray(void) const;
  
  float computeTotalArea(void);
  float computeTotalVolume(void); // ttung: total volume calculation for relative volume attribute

  float* computeMean(float* K);   //  ttung: diffusion

  void  computeGaussCurvature();  //  ttung: local gaussian curvature + 10 x diffusion
  void  computeMeanCurvature();   //  ttung: local mean curvature + 10 x diffusion*
  void  computeCurvatureIndex();
  
  void  resetGaussCurvature();
  void  resetMeanCurvature();
 
  
  /* * * * * * * * * * * * * *
   *                         *
   * Algorithme de Dijkstra  *
   *                         *
   * * * * * * * * * * * * * */
  
  void markAllInfinite(void);
  float computeThreshold(void);
  Point_3* selectNewBase(void);

  float distance(const Point_3 &base,Point_3* const &p) const
  {return base.getDistanceTo(getPointIndex(p));}
 
  void setDistance(int index,float f) 
  {pointNb(index)->setValue(f);}
 
  float getDistance(int index) const
  {return pointNb(index)->getValue();}

  void getConnectedComponent(int index,std::list<int> &indices,bool* const &connected,int &card);
  bool* connectedComponents(int& nbComponents);
  void copyToDistanceTab(Point_3* const &base);
  void dijkstra(Point_3* const&base,float threshold) ;
  float baseArea(Point_3* const &base);
  void oneDijkstraIter(Point_3* &base,float threshold);
  void updateEdgesList(List<Edge*>* const & list,Edge** newEdges);
  //List<int>* getEdgeList(List<int>* const &pointList);
  //List<int>* getFaceList(List<int>* const &pointList);
  //void deleteElements(List<int>* const &pointList);

  void computeMuValues(void);
  
  void computeMuValues2(void);  //  ttung: heigth function
  void computeMuValues3(void);  //  ttung: center of mass function

  DataItem* basesItem(void);
  DataItem* oneBaseItem(void);

  /* * * * * * * * * * * * * * * * * * * * * *
   *                                         *
   * Calcul du Multiresolutionnal Reeb Graph *
   *                                         *
   * * * * * * * * * * * * * * * * * * * * * */

  void subdivideFaceTwoCross(float level,Face* const &f);
  void oneCross2TwoCross(float level,int res,Face* const & f);
  void subdivideFace(float level,int res,Face* const &f,int nbCross);
  void subdivide(int level);

  std::list<MRGNode*> nodesWithRange(muRange range) const;
  MRG* computeMRG(int level);
};
typedef class object Object;

#endif

