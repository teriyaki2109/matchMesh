#ifndef __MRGRAPH_H_
#define __MRGRAPH_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "list.h"
#include "geombasic.h"
#include "byteOrder.h"
//#include "statusBar.h"
#include <string.h>

enum Towards {DOWN,UP,NONE};

// Pour rajouter un attribut, il faut modifier le type attrNames et la constante howManyAttributes
// et bien sur son traitement lors du calcul du MRG

enum attrNames {attrA,attrL};
const int howManyAttributes=2;
const int defaultResMax=7; //see computeSimilarityOneTurn(...) in similarity.cpp line 668
const int AttrSize=64;
const int Ncurv=16;
const int NcordL=16;
const int Ncord1=16;
const int Ncord2=16;



const int NbAttributes=10;

const int Ns=4;
const int Ntheta=8;
const int Nphi=4;
const float T=0.7;

void changeDirection(Towards &direction);

/*>>>>>>>>>> classe TriangleSet <<<<<<<<<<
 *
 * Associe a un MRGNode
 * Decrit un ensemble de triangles connexe
 * Eventuellement, on stocke les faces voisines de cet ensemble pour faciliter les calculs d'adjacences
 *
 */
class TriangleSet
{
public:
  std::list<Face*> innerFaces;
  std::list<Face*> minBorderFaces;
  std::list<Face*> maxBorderFaces;
  
  std::list<Face*> minNeighbours;
  std::list<Face*> maxNeighbours;

  TriangleSet() {}
  TriangleSet(Face* &from,class muRange);
  ~TriangleSet();

  void getTriangleSetFrom(Face* &f,class muRange range);
  void getTriangleSetFrom(std::list<Face*> &theConnectedFaceSet,muRange range);
  int getNbTriangles(void) const;
  void writeTRI(char* name);

  Point_3 computeBarycenter(void);
  float computeVolume(Point_3 &bary);
  float computeClosedArea(Point_3 &bary);
};

float getArea(const TriangleSet &TSet);
float getMuLength(const TriangleSet &TSet);
void TSetConcat(TriangleSet* TS1,TriangleSet* TS2);


/*>>>>>>>>>> classe MRGNodeAttribute <<<<<<<<<< 
 * 
 * Associe a un noeud MRGNode
 * Contient les attributs de ce noeud
 * Les attributs sont utilises au moment du calcul de la similitude entre deux MRG
 *
 */

class MRGNodeAttribute
{
private:
  int nbAttr;
  float* attrArray;

// ttung
public:
  // attribut geometrique
  float Gdistance;
  float Gangle;
  float Gangle2;
  float Volume;
  float ratio;

	float mulength;

  // cordes
  float *cords;
  float *cordsX;
  float *cordsY;


  // hough 3D;
  float *hough3D;

  // courbures
  //float *curv_gauss;
  //float *curv_mean;
  float *curv_index;


  // attribut texture
  //float *red;
  //float *green;
  //float *blue;

  // resolution du noeud
  int resNode; 

  // resolution max (les attributs dépendent de la résolution totale)
  int resMax;

  // niveau du noeud. la position du noeud est [l*(1/2^resNode);(l+1)*(1/2^resNode)[
  int level;

  // description topologique
  int  nbUpNeighborNodes;
  int  nbDownNeighborNodes;
  int  nbUpEndNodes;
  int  nbDownEndNodes;  

public:
  MRGNodeAttribute(int _nbAttr,float* _attrArray);
  MRGNodeAttribute();
  MRGNodeAttribute(int _nbAttr,const TriangleSet &TSet,float totalArea,float totalVolume,int resolution, int nbPoints);
  MRGNodeAttribute(std::list<class MRGNode*> const &childs, int cas);
  MRGNodeAttribute(class MRGNode** childs,int sz);
  ~MRGNodeAttribute();
  MRGNodeAttribute operator+ (MRGNodeAttribute attr);
  

  int getNbAttr(void) const { return nbAttr; }
  float getAttr(int nb) const { return attrArray[nb]; }
  void setAttr(int nb,float _a) { attrArray[nb] = _a; }

  // non utilise
  float getGdist(void) const { return Gdistance; }
  float getGang(void) const { return Gangle; }
};
MRGNodeAttribute* sumOfChilds(List<MRGNode*>* const &childs);

/*>>>>>>>>>> classe MuRange <<<<<<<<<<
 *
 * D'ecrit un intervalle de la fonction mu auquel un MRGNode est associ'e 
 * mu ne peut prendre que des valeurs entre 0 et 1 : tout intervalle est donc inclus dans [0,1]
 * La resolution d'un intervalle est le nombre de fois ou il faut subdiviser egalement [0,1] pour l'obtenir
 */
class muRange
{
private:
  float muMin,muMax;

 public: 
  muRange(void) :
    muMin(0.0f), muMax(0.0f) {}

  muRange(float min,float max) :
    muMin(min), muMax(max){}

  muRange(int res,int order);

  bool operator==(muRange range) const {
    return ((range.muMin == muMin)
         && (range.muMax == muMax));
  }

  // Attention : les deux operateurs <= et << n'ont pas du tout le meme sens
  bool operator<=(muRange range) const { //au sens de l'inclusion
    return ((muMin >= range.muMin)
         && (muMax <= range.muMax));
  }

  bool operator<<(muRange range) const { //avant sur la droite des reels
    return (muMax <= range.muMin);
  }

  int getResolution(void) const;
  int getNumber(void) const;
  float getMin(void) const { return muMin; }
  float getMax(void) const { return muMax; }
  bool isInRange(Point_3* &p) const;
  bool isInRange(Face* &f) const;
  Towards directionOfSubRange(muRange subRange) const; 
};

/*>>>>>>>>>> classe MRGNode <<<<<<<<<< 
 * 
 * Contient un noeud d'un MRG
 *
 */
class MRGNode 
{
private:
  MRGNodeAttribute* attribute;
  muRange range;
  bool coarserLinked;
  std::list<int> MList; //liste des ''matching labels''
  int number; //utilise pour la numerotation des noeuds dans le calcul de la matrice d'adjacence

public:
  TriangleSet* TSet;
  std::list<class MRGEdge*> upwardsIsoEdges;
  std::list<class MRGEdge*> downwardsIsoEdges;
  std::list<class MRGEdge*> coarserNonIsoEdges;
  std::list<class MRGEdge*> finerNonIsoEdges;

  bool toFusion;

  Point_3 coord;
  class matchingEl* matchedEl; // egal au matchingEl dont ce noeud fait partie SEULEMENT quand un matching contenant ce noeud a eu lieu; NULL sinon
 
  MRGNode(void) :
    attribute(NULL),
    range(muRange()),
    coarserLinked(false),
    number(-1),
    TSet(NULL),
    toFusion(false),
    coord(Point_3(0.0f, 0.0f, 0.0f)),
    matchedEl(NULL)
  {}

  MRGNode(MRGNodeAttribute * const &_attribute,muRange _range,TriangleSet* const &_TSet);

  MRGNode(muRange _range,int nbAttr,TriangleSet *&_TSet,float totalArea,float totalVolume, int nbPoints);
  MRGNode(muRange _range,std::list<MRGNode*> &childs);
  MRGNode(std::list<MRGNode*> subnodes);
  ~MRGNode(void);

  muRange getRange(void) const { return range; }
  void setNumber(int nb) { number = nb; }
  int getNumber(void) const { return number; }
  bool ofThisRange(muRange _range) const { return range == _range; }
  void insertEdge(class MRGEdge* &e,Towards direction);
  void markCoarserLinked(void) { coarserLinked = true; }
  bool isCoarserLinked(void) const { return coarserLinked; }
  TriangleSet* getTSet(void) const { return TSet; }
  MRGNodeAttribute* getAttribute(void) const { return attribute; }
  std::list<int> getMList(void) const { return MList; }
  void appendToMList(int label) { MList.push_back(label); }
  void printMList(void);
  void propagateLabel(int label);
  void propagateLabelTowards(int label,Towards direction);
  void setCoord(float _x,float _y,float _z);
  void writeInMRG(FILE* f);
  void cleanAfterMatching(void);
};
void computeFinestAttributes(MRGNode** nodes,int nbNodes,int resolution);

void setNodeLevel(std::list<MRGNode*> &newNodesList,int level);
TriangleSet *TSetConcat(List<class MRGNode*>* childs);
TriangleSet *TSetConcat(class MRGNode** childs,int sz);


/*>>>>>>>>>> classe MRGEdge <<<<<<<<<<
 *
 * Contient une arete de MRG
 * Relie deux MRGNode
 *
 */
class MRGEdge
{
private:
  MRGNode *node1,*node2;
  bool isoLevel;

public:

  MRGEdge(void) :
    node1(NULL),
    node2(NULL)
  {}

  MRGEdge(MRGNode* const &_node1,MRGNode* const &_node2,bool iso) :
    node1(_node1),
    node2(_node2),
    isoLevel(iso) 
  {}
  MRGNode* getNode(int i) const {
    if (i == 0) return node1;
    else if (i == 1) return node2;
    else return NULL;
  }
  void	   setNode(int i, MRGNode* node) {
    if (i == 0) node1 = node;
    else if (i == 1) node2 = node;
  }

  MRGNode* otherNode(MRGNode* const &node) {
    if (node == node1) return node2;
    else if (node == node2) return node1;
    else return NULL;
  }

  bool isIso(void) const { return isoLevel; }

  std::list<int> relatedNodes;

//  Towards getDirection(void);
};
bool nodeIsAmongEdges(MRGNode* const & test,std::list<MRGEdge*> const &l);
std::list<MRGNode*> getOtherNodesList(const std::list<MRGEdge*> &lEdge,MRGNode* node);  

/*>>>>>>>>>> classe MRG <<<<<<<<<< 
 * 
 * Contient un MRG
 *
 */
class MRG
{
private:
  int resolution;
  int nbAttr;

public:
  int *nbNodes;
  MRGNode*** nodes;
  int nbEdges;
  MRGEdge** edges;
  
  MRG() :
    resolution(-1),
    nbAttr(0),
    nbNodes(NULL),
    nodes(NULL),
    nbEdges(0),
    edges(NULL)
  {}

  MRG(int _resolution,int _nbAttr);
  ~MRG();

  int getNbEdges(void) const { return nbEdges; }
  MRGNode** getNodes(int res) const { return nodes[res]; }
  MRGNode* getNode(int res,int index) const { return nodes[res][index]; }
  int getNbNodes(int res) const { return nbNodes[res]; }
  int getTotalNbNodes(void) const;
  void addEdge(MRGEdge* &edge);      
  void addNode(MRGNode* &node,int res);
  MRGEdge** createEdges(int _nbEdges);
  void addEdgeList(std::list<MRGEdge*> &newEdges);
  void addNodeList(std::list<MRGNode*> &newNodes,int res);
  std::list<MRGNode*> nodeArrayToList(int res) const;
  float getASum(int res) const;
  float getLSum(int res) const;
  int getResolution(void) const { return resolution; }
  int getNbAttr(void) const { return nbAttr; }
  int getNodeIndex(MRGNode* const &node,int res) const;
  int getNbIsoEdges(int res) const;
  void writeGraph(char* name,int resLvl);
  void writeMRG(char* name);
  void getBoundingBox(Point_3 &pMin,Point_3 &pMax,int res);
  void moveAllNodes(Point_3 const& offset,int res);
  void setAllNodesCoord(void); //ttung
  void printNodeCoord(MRGNode* node) const;
  void getNodeCoord(MRGNode* node,int &res,int &pos) const;
  void setNodeNumbers(void);
  void extractSubMatrixRec(MRGNode * const & root,bool* &indicesToKeep);
  void extractSubMatrix(MRGNode * const & root,bool* &indicesToKeep);
  void cleanAfterMatching(void);

  friend MRG* readMRG(char* name);
  friend void compare(MRG* const&mrg1,MRG* const &mrg2);

  //ttung
  MRG* createSubMRGfromNode(int res, int index);
  void simplify();

};

MRG*readMRG(char *name);
int findElementInArray(MRGNode** const &array,int size,MRGNode* const &el);
#endif
