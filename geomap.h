// *** PLEASE REFRAIN FROM DISTRIBUTING THIS BETA VERSION ***
// *** USE IT FOR RESEARCH PURPOSE ***
// Tony TUNG (tonytung.org)- 2014 (c)
// Surface Alignment by Geodesic Mapping [Tung et al., PAMI14] & [Tung et al., CVPR10]


#ifndef GEOMAP_H
#define GEOMAP_H

#include <vector>
#include <list>
#include <string>

#include "3dfiles.h" 
#include "object.h" 
#include "mrgraph.h" 
#include "byteOrder.h"
#include "point3.h"
#include "pointsi.h"
#include "colormap.h"

#define matchThres 0.5 //threshold to match first refPoints (0.5 for Horse16-Cat009)
#define DIST 20.0F // ext. points matching: 0.3 or 0.8 for close pose, higher 20.0 for different pose
#define UNIQUE 0.001F // threshold for unicity of a point wrt global geodesic coordinates (square distance)
#define EPSILON 0.1F // threshold for geodesic consistency
#define ambThres 0.2 // ambiguity threshold (percentage) wrt the ambiguity map: refPoints are selected in regions below that threshold (lower=accurate)
#define ggRefNb 6 // nb of refPoints to take into account when computing "local" geodesic distances

using namespace std;

class valRank
{ 
public:
    int index; 
    float value; 
    float distance;
    valRank(int aIndex, float aValue):index(aIndex), value(aValue) {}; 
    valRank(int aIndex, float aValue, float aDistance):index(aIndex), value(aValue), distance(aDistance) {};
    bool operator<(const valRank& b) { return (value < b.value); } 
};

class match
{
	Object* obj1; 
	Object* obj2;
	int nbPoints1;
	int nbFaces1;
	int nbEdges1;
	int nbPoints2;
	int nbFaces2;
	int nbEdges2;
	int nbPoints1_0; //before subdivision
	int nbFaces1_0; //before subdivision
	int nbEdges1_0; //before subdivision
	int nbPoints2_0; //before subdivision
	int nbFaces2_0; //before subdivision
	int nbEdges2_0; //before subdivision
	int res; //MRG matching resolution
	float valThres;
	vector<int> refPoints1;
	vector<int> refPoints2;
    
public:
	match():obj1(NULL),obj2(NULL),nbPoints1(0),nbFaces1(0),nbEdges1(0),nbPoints2(0),nbFaces2(0),nbEdges2(0),nbPoints1_0(0),nbFaces1_0(0),nbEdges1_0(0),nbPoints2_0(0),nbFaces2_0(0),nbEdges2_0(0),res(4),valThres(0.2){};
	~match(){ delete obj1; delete obj2; refPoints1.clear(); refPoints2.clear(); };
	void load(char* file1, char* file2);
	void computeExtremalVertices(Object* obj, vector<int> &extPoints);
	void computeGeodesicToRefPoints(Object* obj, vector<int> &refPoints);
	void computeExtremalPoints();
	void randomReference(Object *obj, vector<int> &refPoints, int N, float aThres, float refSpacing);
	void computeUnicityOnMesh(string outMesh, Object* obj, vector<int> &refPoints, float &thres);
	void matchAllReferences(float aThres, bool checkSpacing=false);
	void printList(list<valRank> &listRef) const;
	void assignReferences(list<valRank> *list1, list<valRank> *list2, int k, int maxSize, int assignedRef[], float &minCost, int &maxMatch, vector<int> &bestAssignment);
	void reorderReferences();
	void retrieve1stRing(Object* obj, Point_3* p, list<int> &candidateList) const;
	void retrieveNeighbors(Object* obj, Point_3* p, list<int> &candidateList, uint k) const; //!suspended!
	float energy(std::vector<int> &refPoints1, std::vector<int> &refPoints2, std::vector<int> &refIndex) const;
	void minEnergy();
	void getClosestRefPoints(point_3* point, vector<int> &refPoints, int nbRefIndex, int firstRef, vector<int> &refIndex);
	float ggDistance(Point_3* pt1, Point_3* pt2, vector<int> &refIndex) const;
	float sqMuDistance(Point_3* pt1, Point_3* pt2, int refSize) const;
	float sqMuConsDistance(Point_3* pt1, Point_3* pt2, int refSize) const;
	void updateGeodesicCoordinate(Object* obj, int cIndex, int pIndex);
	void updateAmbiguity(Object* obj, vector<int> &refPoints, float &thres);
	void outputRefPointsOnMesh(string outMesh, Object* obj, vector<int> &refPoints) const;
	int getMidPathVertex(Object *obj, vector<int> refPoints, int i, int j) const;
	void addMidPathReferences();
	void addReferenceIteratively();
	void SKELdump(string outSKEL) const;
	void getPath(Object* obj, int index, int refindex, std::vector<int> &path) const;
	void outputIsoLinesOnMesh(string outMesh, Object* obj, std::vector<int> &refPoints, int index) const;
	void matchMesh();
};

#endif
