// *** PLEASE REFRAIN FROM DISTRIBUTING THIS BETA VERSION ***
// *** USE IT FOR RESEARCH PURPOSE ***
// Tony TUNG (tonytung.org)- 2014 (c)
// Surface Alignment by Geodesic Mapping [Tung et al., PAMI14] & [Tung et al., CVPR10]

#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime> 
#include <cmath>
#include "geomap.h"

inline float min( float  x, float  y ) { return x < y ? x : y; }
inline float max( float  x, float  y ) { return x > y ? x : y; }

typedef unsigned int uint;
using namespace std;

// Load object meshes from file
void match::load(char* file1, char* file2)
{
	// Load the two input meshes:
	//Load first input mesh (and keep the biggest connected component)
	read3D(obj1, file1, NULL, 0, true); //see 3dfiles.cpp: no texture, normalize coordinates but no PCA 
	cout << ">>" << file1 << " loaded"<< endl;
	nbPoints1_0=obj1->getNbPoints();
	nbFaces1_0=obj1->getNbFaces();
	nbEdges1_0=obj1->getNbEdges();
	cout << ">> Points: " << nbPoints1_0 << " - Triangles: " << nbFaces1_0 << " - Edges: " << nbEdges1_0 << endl;

	//le nb de composantes connexes est calcule a partir du parcourt des aretes 
	int cpt=0;
	bool* pointsToKeep1=obj1->connectedComponents(cpt);     
	if(cpt>1)
	{
		cout << ">> " << cpt << " connected components" << endl; 
		cout << ">> Reloading the mesh keeping only the biggest connected component..." << endl; 
		reloadObject(obj1,(char*)(const char*) file1, pointsToKeep1, 0, true); //see 3dfiles.cpp: no texture, normalize coordinates but no PCA
		nbPoints1_0=obj1->getNbPoints(); 
		nbFaces1_0=obj1->getNbFaces();
		nbEdges1_0=obj1->getNbEdges();                     
		cout << ">> Points: " << nbPoints1_0 << " - Triangles: " << nbFaces1_0 << " - Edges: " << nbEdges1_0 << endl;
	} 
	else
	{ cout << ">> 1 connected component" << endl; } 
	delete [] pointsToKeep1;

	//Load second input mesh (and keep the biggest connected component)
	read3D(obj2, file2, NULL, 0, true); //see 3dfiles.cpp: no texture, normalize coordinates but no PCA
	cout << ">>" << file2 << " loaded"<< endl;
	nbPoints2_0=obj2->getNbPoints(); 
	nbFaces2_0=obj2->getNbFaces();
	nbEdges2_0=obj2->getNbEdges();
	cout << ">> Points: " << nbPoints2_0 << " - Triangles: " << nbFaces2_0 << " - Edges: " << nbEdges2_0 << endl;

	cpt=0;
	bool* pointsToKeep2=obj2->connectedComponents(cpt);     
	if(cpt>1)
	{ 
		cout << ">> " << cpt << " connected components" << endl; 
		cout << ">> Reloading the mesh keeping only the biggest connected component..." << endl; 
		reloadObject(obj2,(char*)(const char*) file2, pointsToKeep2, 0, true); //true=normalize coordinates
		nbPoints2_0=obj2->getNbPoints(); 
		nbFaces2_0=obj2->getNbFaces();
		nbEdges2_0=obj2->getNbEdges();
		cout << ">> Points: " << nbPoints2_0 << " - Triangles: " << nbFaces2_0 << " - Edges: " << nbEdges2_0 << endl;
	}
	else
	{ cout << ">> 1 connected component" << endl; } 
	delete [] pointsToKeep2;
}


void match::computeExtremalPoints()
{
	//Compute Geodesic function on mesh 1 and extract extremal points (store in extPoints1)
	cout << "\n>> Retrieve extremal vertices at res=" << res << "..." << endl;
	obj1->computeMuValues();
	obj1->subdivide(res);
	vector<int> extPoints1; computeExtremalVertices(obj1, extPoints1);

	//Compute Geodesic function on mesh 2 and extract extremal points (store in extPoints2)
	obj2->computeMuValues();
	obj2->subdivide(res);
	vector<int> extPoints2; computeExtremalVertices(obj2, extPoints2);
	
	for(uint i=0; i<extPoints1.size(); i++) { if(extPoints1[i]<nbPoints1_0) refPoints1.push_back(extPoints1[i]); } extPoints1.clear(); //copy extPoints1 in RefPoints1
	for(uint i=0; i<extPoints2.size(); i++) { if(extPoints2[i]<nbPoints2_0) refPoints2.push_back(extPoints2[i]); } extPoints2.clear(); //copy extPoints2 in RefPoints2

	nbPoints1=obj1->getNbPoints();
	nbFaces1=obj1->getNbFaces();
	nbEdges1=obj1->getNbEdges();   
	cout << ">> Points: " << nbPoints1 << " - Triangles: " << nbFaces1 << " - Edges: " << nbEdges1 << " after mesh subdivision" << endl;

	nbPoints2=obj2->getNbPoints(); 
	nbFaces2=obj2->getNbFaces();
	nbEdges2=obj2->getNbEdges();
	cout << ">> Points: " << nbPoints2 << " - Triangles: " << nbFaces2 << " - Edges: " << nbEdges2 << " after mesh subdivision" << endl;
}

void match::printList(list<valRank> &listRef) const
{
	list<valRank>::iterator it=listRef.begin(), itEnd=listRef.end();
	for(; it!=itEnd; it++)
	{
		cout << it->index << " " << it->value << endl;
	}
	cout << endl;
}

void match::assignReferences(list<valRank> *list1, list<valRank> *list2, int k, int maxSize, int assignedRef[], float &minCost, int &maxMatch, vector<int> &bestAssignment)
{
	//stop criteria
	if(k == maxSize)
	{
		float cost=0;
		int nbMatch=0;
		for(int i=0; i<maxSize; i++)
		{
			if(assignedRef[i]>=0) //*it==-1 if refPoints1[i] is unassigned
			{
				point_3* pointRef1=obj1->pointNb(refPoints1[i]);
				point_3* pointRef2=obj2->pointNb(refPoints2[assignedRef[i]]);
				float tmpDist= (*pointRef1-*pointRef2).norm();
				//if(tmpDist>0.1) continue; //!! this discards the matching between distant nodes!
				cost += tmpDist;
				nbMatch++;
			}
			//cout << assignedRef[i] << " ";
		}
		//cout << endl; cout << "cost=" << cost << endl;

		if(nbMatch > maxMatch || (cost<minCost && nbMatch >= maxMatch)) //keep min matching cost with max number of matches
		{
			minCost=cost;
			maxMatch=nbMatch;
			bestAssignment.clear();
			for(int j=0; j<maxSize; j++)
			{
				bestAssignment.push_back(assignedRef[j]); 
			}
		}

		return;
	}

	//find good candidates at each step (ranked as closest RefPoints2 and unmatched yet)
	if(!list1[k].empty())
	{
		list<valRank>::iterator itC=list1[k].begin(), itCEnd=list1[k].end();
		for(; itC!=itCEnd; itC++)
		{
			bool notMatched=true;
			for(int i=0; i<k; i++)
			{
				if(itC->index == assignedRef[i])
				{
					notMatched=false;
					break;	
				}
			}

			if(notMatched)
			{
				assignedRef[k]=itC->index;
			}
			else
			{
				assignedRef[k]=-1;
			}
			assignReferences(list1, list2, k+1, maxSize, assignedRef, minCost, maxMatch, bestAssignment);
		}
	}
	else
	{
			assignedRef[k]=-1;
			assignReferences(list1, list2, k+1, maxSize, assignedRef, minCost, maxMatch, bestAssignment);
	}
}


// reorder the reference points in refPoints1 so that they corresond to points in refPoints2
// !!! refPoints1 and refPoints2 can be modified !!!
// extremal vertices which are not matched are removed from the list
// usually refPoints1 and refPoints 2 have small size, so this is fast
void match::reorderReferences()
{
	int refSize1=refPoints1.size();
	int refSize2=refPoints2.size();
	std::vector<int> tmpRef1;
	std::vector<int> tmpRef2;
	list<valRank> *list1 = new list<valRank>[refSize1];
	list<valRank> *list2 = new list<valRank>[refSize2];

	cout << "Get closest refPoints2 to each refPoints1..." << endl;

	// we go through refPoints1 and compute the distances to refPoints2
	for(int i=0; i<refSize1; i++)
	{
		point_3* pointRef1=obj1->pointNb(refPoints1[i]);
		for(int j=0; j<refSize2; j++)
		{
			point_3* pointRef2=obj2->pointNb(refPoints2[j]);
			float tmpDist=(*pointRef1-*pointRef2).norm(); //cout << j << " " << tmpDist << endl;
			if(tmpDist<matchThres)
			{
				list1[i].push_back(valRank(j,tmpDist));
			}
		}
		list1[i].sort();
		//printList(list1[i]);
		tmpRef1.push_back(refPoints1[i]);
	}

	cout << "Get closest refPoints1 to each refPoints2..." << endl;
	for(int i=0; i<refSize2; i++)
	{
		/*
		point_3* pointRef2=obj2->pointNb(refPoints2[i]);
		for(int j=0; j<refSize1; j++)
		{
			point_3* pointRef1=obj1->pointNb(tmpRef1[j]);
			float tmpDist=(*pointRef1-*pointRef2).norm();
			if(tmpDist<0.5)
			{
				list2[i].push_back(valRank(j,tmpDist));
			}
		}
		list2[i].sort();
		printList(list2[i]);
		*/

		tmpRef2.push_back(refPoints2[i]);

		/*
		list<valRank>::iterator it=list2[i].begin(), itEnd=list2[i].end();		
		//try with unique closests
		if(i == list1[list2[i].front().index].front().index)
		{
			point_3* tmpPoint1=obj1->pointNb(tmpRef1[list2[i].front().index]);
			float tmpDist2=(*pointRef2-*tmpPoint1).norm(); //cout << j << " " << tmpDist2 << endl;
			if(tmpDist2<DIST)
			{
				tmpRef2.push_back(refPoints2[i]);
				refPoints1.push_back(tmpRef1[list2[i].front().index]);
				cout << "+";
			}
		}
		*/
	}

	cout << "Match refPoints1 to refPoints2...\n" << endl;
	cout << "refSize1=" << refSize1 << endl;
	cout << "refSize2=" << refSize2 << endl;

	int assignedRef[refSize1]; for (int i=0; i<refSize1; i++) assignedRef[i]=-1;
	float minCost=1e9;
	int maxMatch=0;
	vector<int> bestAssignment;
	assignReferences(list1, list2, 0, refSize1, assignedRef, minCost, maxMatch, bestAssignment);

	cout << "Min cost=" << minCost << endl;

	// update refPoints2
	refPoints2.clear();
	refPoints1.clear();
	for(int i=0; i<refSize1; i++)
	{
		if(bestAssignment[i]!=-1)
		{
			refPoints1.push_back(tmpRef1[i]);
			refPoints2.push_back(tmpRef2[bestAssignment[i]]);
			cout << ">> ref1=" << refPoints1[i] << " - ref2=" << refPoints2[i] << endl;
		}
	}
}

// get vertex index located between refPoint i and j
int match::getMidPathVertex(Object *obj, vector<int> refPoints, int i, int j) const
{
	float pathLength=obj->pointNb(refPoints[j])->field.value[i]; //distance of refPoints[j] to refPoints[i]
	cout << "length=" << pathLength << endl;
	float midLength=pathLength*0.5;
	int index=refPoints[j];
	
	while(pathLength>midLength)
	{
		list<int> ring;
		retrieve1stRing(obj, obj->pointNb(index), ring);
		list<int>::iterator it=ring.begin(), itEnd=ring.end();
		float minPathLength=1e9, tmpLength;
		int closest=*it;
		for(; it!=itEnd; it++)
		{
			tmpLength=obj->pointNb(*it)->field.value[i];
			if(tmpLength<minPathLength)
			{
				minPathLength=tmpLength;
				closest=*it;
			}
		}
		
		pathLength=obj->pointNb(closest)->field.value[i];
		index=closest;
	}

	cout << "midlength=" << pathLength << endl;

	return index;
}

// add reference points between each pair of existing pairs of references
void match::addMidPathReferences()
{
	int nbRef=refPoints1.size();	
	for(int i=0; i<nbRef; i++)
	{
		for(int j=i+1; j<nbRef; j++)
		{
			int index1= getMidPathVertex(obj1, refPoints1, i, j);
			int index2= getMidPathVertex(obj2, refPoints2, i, j);
			refPoints1.push_back(index1);
			refPoints2.push_back(index2);
		}
	}
}

void match::addReferenceIteratively()
{
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ADD HERE IN refPoints1 AND refPoints2 SOME SPARSE REFERENCE POINTS
	// FOUND FROM COLOR INFORMATION SUCH THAT (refPoints1[i], refPoints2[i])
	// IS A PAIR ~ IF YOU WANT~
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	cout << "Compute geodesics..." << endl;
	computeGeodesicToRefPoints(obj1, refPoints1); cout << endl; //compute all geodesics to refPoints1 from mesh1 vertices
	computeGeodesicToRefPoints(obj2, refPoints2); //compute all geodesics to refPoints1 from mesh1 vertices

//	addMidPathReferences(); //add reference points between every pair of references: unstable? symmetry issue?
//	computeGeodesicToRefPoints(obj1, refPoints1); cout << endl;
//	computeGeodesicToRefPoints(obj2, refPoints2);

	//output meshes with highlighted reference points (refPoints are diffused to neighbors)
	outputRefPointsOnMesh("output/outMesh1.off", obj1, refPoints1);
	outputRefPointsOnMesh("output/outMesh2.off", obj2, refPoints2);

	SKELdump("output/match_0.off");

	int sumRef=refPoints1.size();
	for(int k=1; k<4; k++)
	{
		//update and output uniqueness degree of each vertex wrt refPoints (update obj->pointNb()->field.uniqueness)
		stringstream outUnique1; outUnique1 << "output/outUnique1_" << sumRef << ".off";
		stringstream outUnique2; outUnique2 << "output/outUnique2_" << sumRef << ".off";
		float aThres1, aThres2;
		computeUnicityOnMesh(outUnique1.str(), obj1, refPoints1, aThres1);
		computeUnicityOnMesh(outUnique2.str(), obj2, refPoints2, aThres2);
		float aThres=min(aThres1,aThres2);

		//add random reference points on mesh2 (unicity is checked to avoid ambiguous location, refPoints2 is updated)
		randomReference(obj2, refPoints2, k*20, aThres, 0.3); // points are well spaced on the surface
		sumRef+=k*20;

		matchAllReferences(aThres, false); cout << endl; //this is not accurate when deformations are too strong
		minEnergy(); //optimize refPoints1 positions

		//output meshes with highlighted reference points (refPoints are diffused to neighbors)
		stringstream outMesh1; outMesh1 << "output/outMesh1_" << sumRef << ".off";
		stringstream outMesh2; outMesh2 << "output/outMesh2_" << sumRef << ".off";
		outputRefPointsOnMesh(outMesh1.str(), obj1, refPoints1);
		outputRefPointsOnMesh(outMesh2.str(), obj2, refPoints2);

		stringstream outMatch; outMatch << "output/match_" << k << ".off";
		SKELdump(outMatch.str());
	}
}


// get the extremal points using Morse function values computed on mesh surface
// require pre-computation of computeMu() and subdivide() at resolution res
// we use the "Reeb graph" structure at level of resolution res and keep extremal maximal points (not minimal) of interest!
void match::computeExtremalVertices(Object* obj, vector<int> &extPoints)
{
    MRG* mrg=obj->computeMRG(res);
    if(!mrg) { cout << "Error computing mrg!" << endl; return; }
    int nbNodes=mrg->nbNodes[res]; //mrg->writeMRG("output/toto.mrg");
    int nbPoints=obj->getNbPoints();

    // go through all the nodes at level of resolution res
    for(int cpt=0; cpt<nbNodes; ++cpt) 
    { 
        MRGNode* oneNode=mrg->nodes[res][cpt]; if (oneNode==NULL) continue; 
        /*if(oneNode->getAttribute()->nbDownEndNodes==1) //extrema (minima) of Mu
        {
            cout << ">> Extremal node " << cpt << " (minima)--> ";
            Point_3* ptMin=NULL; float tmpMu, muMin=888;
            //look for the extremal vertex in the region (submesh) corresponding to the node
            list<Face*> faces=oneNode->TSet->innerFaces;
            for(faces=oneNode->TSet->innerFaces;!faces.empty();faces.pop_front()) 
            {
                for(i=0; i<3; i++)
                {
                    Point_3* tmpPt=faces.front()->p[i];
                    tmpMu=tmpPt->getValue();                            
                    if(tmpMu<muMin)
                    {
                        muMin=tmpMu;
                        ptMin=tmpPt;
                    }
                }
            }
            
            int indMin=obj->getPointIndex(ptMin);
            cout << "index Min " << indMin << endl;
            //refPoint.push_back(indMin);
        }
        else */                   
        if(oneNode->getAttribute()->nbUpEndNodes==1) // extrema (maxima) of Mu //if(oneNode->upwardsIsoEdges.size()==0) //oneNode->getAttribute()->getAttr(0) > 1e-5
        {
	    // look for the vertex in the region (submesh) corresponding to the node having max value
            list<Face*> faces;
	    if(!oneNode->TSet->innerFaces.empty()) { faces=oneNode->TSet->innerFaces; }
	    else if(!oneNode->TSet->minBorderFaces.empty()) { faces=oneNode->TSet->minBorderFaces; }
	    else { continue; }

	    list<valRank> candidates;

            while(!faces.empty()) 
            {
                for(int i=0; i<3; i++)
                {
                    Point_3* tmpPt=faces.front()->p[i];
		    if(tmpPt)
		    {
			float tmpMu=tmpPt->getValue();
			int index=obj->getPointIndex(tmpPt);
			if(index>=0 && index < nbPoints ) { candidates.push_back(valRank(index, tmpMu)); }
		    }
                }
		faces.pop_front();
            }
            
	    if(!candidates.empty())
	    {
		cout << ">> Extremal node " << cpt << " (maxima)--> ";
		candidates.sort();
		int indMax=candidates.back().index;
		Point_3* pt=obj->pointNb(indMax);
                cout << "vert. index " << indMax<< " @ P(" << pt->x << ", " << pt->y << ", " << pt->z <<"): Mu=" << pt->getValue()
		     << ", Area=" << oneNode->getAttribute()->getAttr(0) << endl;
		extPoints.push_back(indMax);		
	    }
        } 
    }// end forall MRG nodes at res
}


// Compute geodesic distances to reference points for every vertex on the mesh
// distances are stored in the attributes of vertices (vector field)
void match::computeGeodesicToRefPoints(Object* obj, vector<int> &refPoints)
{
	int nbPoints=obj->getNbPoints();
	float threshold=obj->computeThreshold();
	int refSize=refPoints.size();

	// store reference points in vector reference[i] with geodesic distances to all vertices
	vector<Point_3*> reference;
	vector<int>::iterator refPtIt=refPoints.begin();
	for(; refPtIt!=refPoints.end(); refPtIt++)
	{
		Point_3* tmpRefPoint=obj->pointNb(*refPtIt);
		tmpRefPoint->markBase();
		tmpRefPoint->createDistanceTab(nbPoints);

		// compute the geodesic distances from all vertices to each reference point
		// all distances are stored in a (big) distanceTab table for each reference point
		// cout << "## Computing all geodesic distances to point " << *refPtIt << "..." << endl;
		obj->dijkstra(tmpRefPoint, threshold);
		reference.push_back(tmpRefPoint);
	}

	//obj->markAllInfinite(); //pour que tous les champs value soient a` zero
	// copy the geodesic distances from the refPoints.DistanceTab to each vertex of the mesh
	// we use the attribute fields p->field.value[i] to store the geodesic distances to each reference point i
	// see geombasic.h: max value[] size is 300!! --> 300 refPoints is the max
	vector<float> min; min.reserve(refSize);
	vector<float> max; max.reserve(refSize);
	vector<float> value; value.reserve(refSize);
	for(int i=0; i<refSize; i++) { min[i]=1e9; max[i]=0.0; }

	for(int cpt=0; cpt<nbPoints; cpt++)
	for(int i=0; i<refSize; i++)
	{
		value[i]=reference[i]->getDistanceTo(cpt);
		obj->pointNb(cpt)->field.value[i]=value[i];

		// compute min and max of geodesic distances in order to normalize!
		if (value[i]<min[i]) { min[i]=value[i]; }
		else if (value[i]>max[i]) { max[i]=value[i]; }
	}

	// normalize geodesics to [0,1] for each refPoints
	Point_3* p;
	for (int cpt=0; cpt<nbPoints; cpt++)
	{
		p=obj->pointNb(cpt); 
		for(int i=0; i<refSize; i++) { p->field.value[i]=((p->field.value[i]-min[i])/(max[i]-min[i])); }
	}
}


// compute uniqueness degree of each vertex wrt refPoints (to check position ambiguity)
// update member variable, and output mesh
void match::computeUnicityOnMesh(string outMesh, Object* obj, vector<int> &refPoints, float &thres)
{
	cout << "Writing ambiguity map mesh " << outMesh << endl;
	ofstream outMeshFile(outMesh.c_str(), ios::out);
	if (!outMeshFile.is_open()) { cout << "Error writing COFF file " << outMesh << endl; return; }
	else
	{
		int nbPoints=obj->getNbPoints();
        	int nbFaces=obj->getNbFaces();
       		int nbEdges=obj->getNbEdges();
		int refSize=refPoints.size();
		int uMax=0, uMin=1e9;

		// Using refIndex (instead of refPoints): if nbRefPoints > ggRefNb, get the ggRefNb closest refPoints to pointi
		int nbRefIndex=refSize; if(refSize>ggRefNb) { nbRefIndex=ggRefNb; }

		for(int i=0; i<nbPoints; ++i)
		{
			// get the visited point i and go through all vertices and check unicity
			point_3* pointi=obj->pointNb(i); pointi->field.uniqueness=0;			
			vector<int> refIndex; refIndex.reserve(nbRefIndex); getClosestRefPoints(pointi, refPoints, nbRefIndex, -1, refIndex);

		    	for(int j=0; j<nbPoints; ++j)
			{
				// check if there exists a point within geodesic coordinate distance (UNIQUE=0.001F)
				// sqMuDistance uses ALL the refPoints => we select only the very unambiguous points!
				// ggDistance uses a list of (closest) refPoints => better to handle surface deformations
				point_3* pointj=obj->pointNb(j);
				if(nbRefIndex!=refSize) { if(ggDistance(pointi, pointj, refIndex)<UNIQUE) { obj->pointNb(i)->field.uniqueness++; } } else
				{ if(sqMuDistance(pointi, pointj, refSize)<UNIQUE) { obj->pointNb(i)->field.uniqueness++; } }
			}
			if(obj->pointNb(i)->field.uniqueness > uMax) { uMax=obj->pointNb(i)->field.uniqueness; }
			else if(obj->pointNb(i)->field.uniqueness < uMin) { uMin=obj->pointNb(i)->field.uniqueness; }
		}

		thres=(float)uMin+ambThres*(float)(uMax-uMin); //give the ambiguity threshold for this mesh
//		printf("Min=%d, Max=%d -->", uMin, uMax); printf("Ambiguity threshold (%2.1f\%)=%2.3f\n", ambThres*100, thres);
//		int counter=0; for(int i=0; i<nbPoints; ++i) { if(obj->pointNb(i)->field.uniqueness<thres) counter++; } //count number of ambiguous points
//		printf("Non-ambiguous vertex ratio: %2.2f\% (out of %d)\n", (float) counter *100.0/(float) nbPoints, nbPoints);
		
       		// write meshes with ambiguity degrees
		outMeshFile << "COFF" << endl;
		outMeshFile << nbPoints << " " << nbFaces << " " << nbEdges << endl;
		outMeshFile << "#vertex" << endl;

		// go through all the vertices again and color the mesh with uniqueness degree
		colorMapFunc cmf=ColorMap::selectColorMap(1); //1: hot colormap
		for(int i=0; i<nbPoints; ++i)
		{
			int uDegree=obj->pointNb(i)->field.uniqueness;
			unsigned char color[3]; cmf(color, uDegree, uMin, uMax);		

			// write 3D coordinates and rgb color corresponding to uniqueness degree 
			outMeshFile << obj->pointNb(i)->x << " " << obj->pointNb(i)->y << " " << obj->pointNb(i)->z << " "
			<< float(color[0])/255.0 << " " <<  float(color[1])/255.0 << " " << float(color[2])/255.0 << " " << uDegree << endl;
		}

		outMeshFile << "#face" << endl;
		for(int j=0; j<nbFaces; ++j)
		{
			outMeshFile << "3 " << obj->getPointIndex(obj->faceNb(j)->p[0]) << " " << obj->getPointIndex(obj->faceNb(j)->p[1]) << " "
			<< obj->getPointIndex(obj->faceNb(j)->p[2]) << endl;
		}

	}

	outMeshFile.close();
}


// choose randomly N points on the mesh obj (actually Mesh2 only!)
// unicity degree is checked to avoid "poor" reference points which would lead to ambiguous localization
// refSpacing (=0.3): min distance between refPoints in the search process BB in [-1,1]
void match::randomReference(Object *obj, vector<int> &refPoints, int N, float aThres, float refSpacing)
{
	int nbPoints=obj->getNbPoints(), searchLoop=0, j=0;
	int nbPointMax=nbPoints;
	if(obj==obj2) { nbPointMax=nbPoints2_0; } else if(obj==obj1) { nbPointMax=nbPoints1_0; } 
	float threshold=obj->computeThreshold();
	cout << ">> Trying to add " << N << " reference points randomly on Mesh2... ";
	//Choose randomly N points amongst dstMesh
	vector<int> randomRef;
	srand ( time(NULL) );
	while(j<N)
	{
		// we give it 200 trials to find a good refpoint randomly else the search spacing is reduced (refPoints are getting closer)
		if(searchLoop > 200) { searchLoop=0; refSpacing-=0.025; if (refSpacing <= 0.05) break; }		
		int r=rand() % nbPointMax; // generate ramdom point index

		// check uniqueness degree of the point (uThreshold=30 with 5 refPoints)
		// same values as in matchAllReferences()
		if(obj->pointNb(r)->field.uniqueness > aThres) { searchLoop++; continue; }

		//the new references have to be at least 0.3 distant to other ref points
		bool isOk=true;
		for(uint i=0; i<refPoints.size(); i++) { if(obj->pointNb(refPoints[i])->getDistanceTo(r) < refSpacing) { isOk=false; break; } }

		// if the point is far enough from all others, insert it in the refPoints list and update the geodesics of all vertices
		if(isOk)
		{
			refPoints.push_back(r);

			//update value field with normalized geodesics
			Point_3* tmpRefPoint=obj->pointNb(r);
			tmpRefPoint->markBase();
			tmpRefPoint->createDistanceTab(nbPoints);
			obj->dijkstra(tmpRefPoint, threshold);

			float gmin=1e9, gmax=0;
			int lastRefIndex=refPoints.size()-1;
			for (int k=0; k<nbPoints; k++)
			{
				float tmpValue=tmpRefPoint->getDistanceTo(k);
				obj->pointNb(k)->field.value[lastRefIndex]=tmpValue;
				if (tmpValue<gmin) { gmin=tmpValue; }
				else if (tmpValue>gmax) { gmax=tmpValue; }
			}

			for (int k=0; k<nbPoints; k++)
			{
				Point_3* p=obj->pointNb(k); 
				p->field.value[lastRefIndex]=(p->field.value[lastRefIndex]-gmin)/(gmax-gmin);
			}

			j++;
			searchLoop=0;
		}
		else { searchLoop++; }
	}

	// update geodesics of all vertices to the refPoints
	cout << "There are now "<< refPoints.size() << " refpoints on Mesh2!" << endl;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------


// compute geodesic coordinate square distance only on certain (geodesically consistent) components
float match::sqMuConsDistance(Point_3* pt1, Point_3* pt2, int refSize) const
{
	float dist=0.0;
	for(int k=0; k<refSize; k++)
	{
		float tmpDist=(pt1->field.value[k]-pt2->field.value[k])*(pt1->field.value[k]-pt2->field.value[k]);
		if(tmpDist > EPSILON) continue;
		else dist+=tmpDist; 
	}
	return dist;
}


// Update geodesic coordinates: cIndex-th refPoint is replaced by pIndex-th vertex
void match::updateGeodesicCoordinate(Object* obj, int cIndex, int pIndex)
{
	int nbPoints=obj->getNbPoints();
	float threshold=obj->computeThreshold();

	// compute all geodesics to newRefpoint
	Point_3* newRefPoint=obj->pointNb(pIndex);
	newRefPoint->markBase();
	newRefPoint->createDistanceTab(nbPoints);
	obj->dijkstra(newRefPoint, threshold);

	float min=1e9, max=0.0, value;
	for (int i=0; i<nbPoints; i++)
	{
		value=newRefPoint->getDistanceTo(i);
		obj->pointNb(i)->field.value[cIndex]=value;
		if (value<min) { min=value; } else if (value>max) { max=value; }
//		if (value==0.0f) { cout << "!!! value of index " << i << " =0" << endl;}
	}

	// normalize geodesics to [0,1]
	for (int i=0;i<nbPoints;i++)
	{
		value=obj->pointNb(i)->field.value[cIndex];
		obj->pointNb(i)->field.value[cIndex]=(value-min)/(max-min);
	}

}


// then, there are additional points in the refPoints list that have to be matched
// match all reference points in refPoints2 to the closest point in obj1 (using global geodesic distance)
// points that fall on an ambiguous position in Mesh1 are not kept
void match::matchAllReferences(float aThres, bool checkSpacing)
{
	cout << "CheckSpacing="; if(checkSpacing) { cout << "TRUE" << endl; } else { cout << "FALSE" << endl; }
	float threshold1=obj1->computeThreshold();
	int nbPoints1=obj1->getNbPoints();
	int refSize1=refPoints1.size();
	int refSize2=refPoints2.size();
	vector<int> tmpRef1;
	vector<int> tmpRef2;

	//copy already matched reference points in tmplists
	for(int i=0; i<refSize1; i++) { tmpRef1.push_back(refPoints1[i]); tmpRef2.push_back(refPoints2[i]); }

	//same as in computeUnicityOnMesh(): if refSize1>ggRefNb, get the ggRefNb closest refPoints2 in Mesh2 to the visited point
	int nbRefIndex=refSize1; if(refSize1>ggRefNb) { nbRefIndex=ggRefNb; }

	//go through the newly added reference points in list2
	for(int i=refSize1; i<refSize2; i++)
	{
		int closestMuPointIndex=-1;
	    	point_3* pointRef2=obj2->pointNb(refPoints2[i]); // get the pointRef2 coordinates to be matched

		// sort closest refPoints to the visited pointRef2
		int tmpSize1=tmpRef1.size(); //update size
		vector<int> refIndex2; refIndex2.reserve(tmpSize1);
		getClosestRefPoints(pointRef2, tmpRef2, nbRefIndex, -1, refIndex2); //get the 6 closest refpoints from {tmpRef2} to pointRef2

		float minDist2=1e9; for(int j=0; j<tmpSize1; j++) { if(pointRef2->field.value[j]<minDist2) { minDist2=pointRef2->field.value[j]; } }
		
		// look for the closest vertex in whole Mesh1 using global geodesic coordinates wrt pointRef2
		// i.e., get candidates (index, nbhits) in Mesh1 which are geodesically consistent to pointRef2
		list<valRank> candidates1; 
		int maxCons=0;
				
		for(int j=0; j<nbPoints1_0; j++) //we take candidates in the non-subdivided mesh1 (refPoints2 are in the non-subdivided mesh2 anyway)
		{
			point_3* pointj=obj1->pointNb(j); // get visited point j
			if(pointj->field.uniqueness > aThres) { continue; } // check ambiguity value on Mesh1
			//if((pointRef2->norm() > 0.1) && (pointj->norm() > 0.1) && (*pointRef2 * *pointj < 0)) continue; // fix symmetry issue!
			
			//count the number of geodesically consistent coordinate components (this should help to overcome topology changes and large surface variations!)
			int nbConsistentCoordComponents=0;
			for(int k=0; k<tmpSize1; k++) { if(fabs(pointRef2->field.value[k]-pointj->field.value[k]) < EPSILON ) nbConsistentCoordComponents++; } // EPSILON=0.1
			if(nbConsistentCoordComponents < tmpSize1*0.6) continue; // visited point is a candidate if 70% of coordinate components are consistent
			
			if(checkSpacing) //check minimum spacing consistency between the candidate and existing {refPoints1} : not very good criteria
			{
				bool tooClose=false; for(int k=0; k<tmpSize1; k++) { if(pointj->field.value[k]<minDist2*0.01) { tooClose=true; break; } } if(tooClose) continue;
			}

			//float tmpMuDist=sqMuDistance(pointRef2, pointj, tmpSize1); // global geodesic distance between pointRef2 and pointj
			//float tmpMuDist=sqMuConsDistance(pointRef2, pointj, tmpSize1); // global geodesic distance between pointRef2 and pointj wrt closest refpoints
			float tmpMuDist=ggDistance(pointRef2, pointj, refIndex2); // global geodesic distance considering nbRefIndex closest

			if(tmpMuDist<EPSILON) //check overall gg distance consistency
			{
				valRank cand(j, nbConsistentCoordComponents, tmpMuDist);
				candidates1.push_back(cand); //visited point is a valid candidate
				if(nbConsistentCoordComponents>maxCons) { maxCons=nbConsistentCoordComponents; }
			}
		}
		
		// Manage candidates: they all have consistent coordinates with pointRef2!
		int candSize=candidates1.size(); //cout << "There are "<< candSize << " candidates" << endl;
		float minMuDist=1e9;
		if(candSize)
		{
			list<valRank>::iterator valIt=candidates1.begin(), valItE=candidates1.end();
			for(; valIt!=valItE; valIt++)
			{
				//HERE: ADD MORE SELECTION CRITERIA??
				if((*valIt).value==maxCons)
				{
					//choose the candidate with lowest gg distance and max number of consistent coordinate components
					float tmpDist=(*valIt).distance;
					if (tmpDist<minMuDist) { minMuDist=tmpDist; closestMuPointIndex=(*valIt).index; } 
				}
			}
		}

		if(minMuDist==1e9 || candSize==0) { /* cout << "No candidate found for " << i << endl; */ continue; }
		
		tmpRef1.push_back(closestMuPointIndex); tmpRef2.push_back(refPoints2[i]); //keep the pair; closestMuPointIndex=-1 if no candidate was found!

		//update value field with normalized geodesics
		Point_3* tmpRefPoint=obj1->pointNb(closestMuPointIndex);
		tmpRefPoint->markBase();
		tmpRefPoint->createDistanceTab(nbPoints1);
		obj1->dijkstra(tmpRefPoint, threshold1);

		float gmin=1e9, gmax=0;
		int lastRefIndex=tmpRef1.size()-1;
		for (int k=0; k<nbPoints1; k++)
		{
			float tmpValue=tmpRefPoint->getDistanceTo(k);
			obj1->pointNb(k)->field.value[lastRefIndex]=tmpValue;
			if (tmpValue<gmin) { gmin=tmpValue; }
			else if (tmpValue>gmax) { gmax=tmpValue; }
		}

		for (int k=0; k<nbPoints1; k++)
		{
			Point_3* p=obj1->pointNb(k); 
			p->field.value[lastRefIndex]=(p->field.value[lastRefIndex]-gmin)/(gmax-gmin);
		}

		updateAmbiguity(obj1, tmpRef1, aThres); //update ambiguity values	
	}

	// update the lists of reference points without the ambiguous ones
	refPoints1.clear(); refPoints2.clear();
	for(uint i=0; i<tmpRef1.size(); i++) { if(tmpRef1[i]!=-1) { refPoints2.push_back(tmpRef2[i]); refPoints1.push_back(tmpRef1[i]); } }
	tmpRef1.clear(); tmpRef2.clear();
}


// get geodesic path from a point index in mesh obj, to the refPoint identified by refindex
void match::getPath(Object* obj, int index, int refindex, std::vector<int> &path) const
{
	float pathLength=obj->pointNb(index)->field.value[refindex];
//	cout << "pathLength= " << pathLength << endl;
	while(pathLength>0.0)
	{
		list<int> ring;
		retrieve1stRing(obj, obj->pointNb(index), ring);
//		int ringSize=ring.size(); cout << "ringSize= " << ringSize << endl;

		list<int>::iterator it=ring.begin(), itEnd=ring.end();
		float minPathLength=1e10, tmpLength; int closest=*it;
		for(; it!=itEnd; it++)
		{
			tmpLength=obj->pointNb(*it)->field.value[refindex];
			if(tmpLength<minPathLength)
			{
				minPathLength=tmpLength;
				closest=*it;
			}
		}
		
		path.push_back(closest);
		pathLength=obj->pointNb(closest)->field.value[refindex];
		index=closest;
//		cout << "pathLength: " << pathLength << endl;
	}
//	stringstream outPath; outPath << "output/path.off"; outputRefPointsOnMesh(outPath.str(), obj, path);
}



// compute geodesic coordinate square distance
float match::sqMuDistance(Point_3* pt1, Point_3* pt2, int refSize) const
{
    int k; float tmpDist=0.0;
    for(k=0; k<refSize; k++) { tmpDist+=(pt1->field.value[k]-pt2->field.value[k])*(pt1->field.value[k]-pt2->field.value[k]); }
    return tmpDist;
}

// compute global geodesic distance using the indices given in ref
float match::ggDistance(Point_3* pt1, Point_3* pt2, std::vector<int> &refIndex) const
{
	float tmpDist=0.0;
	vector<int>::iterator it=refIndex.begin(), itEnd=refIndex.end();
//	list<float> geoDist;
	for(; it!=itEnd; it++)
	{
//		geoDist.push_back((pt1->field.value[*it]-pt2->field.value[*it])*(pt1->field.value[*it]-pt2->field.value[*it]));
		tmpDist+=(pt1->field.value[*it]-pt2->field.value[*it])*(pt1->field.value[*it]-pt2->field.value[*it]);
	}
/*
	geoDist.sort();

	// compute the distance from the 5 smallest differences
	int nbRefPoints=5;
	list<float>::iterator gdIt=geoDist.begin();
	list<float>::iterator gdItEnd=geoDist.end();
	for(int i=0; i<5, gdIt!=gdItEnd; i++, gdIt++) { tmpDist+=*gdIt; }
*/
	return tmpDist; //beware: it's a square distance!
}

// update Ambiguity on mesh, update member variable (field.uniqueness)
void match::updateAmbiguity(Object* obj, vector<int> &refPoints, float &thres)
{
	int nbPoints=obj->getNbPoints();
	int refSize=refPoints.size();
	int nbRefIndex=refSize; if(refSize>ggRefNb) { nbRefIndex=ggRefNb; }

	int uMax=0, uMin=1e9;
	for(int i=0; i<nbPoints; ++i)
	{
		point_3* pointi=obj->pointNb(i); // get the visited point i and reset ambiguity value
		if(pointi->field.uniqueness >= thres) //  to save time, update ONLY ambiguous points!
		{
			vector<int> refIndex; refIndex.reserve(nbRefIndex); getClosestRefPoints(pointi, refPoints, nbRefIndex, -1, refIndex);
			obj->pointNb(i)->field.uniqueness=0;
			for(int j=0; j<nbPoints; ++j)
			{
				point_3* pointj=obj->pointNb(j);
				if(nbRefIndex!=refSize) { if(ggDistance(pointi, pointj, refIndex)<UNIQUE) { obj->pointNb(i)->field.uniqueness++; } } else
				{ if(sqMuDistance(pointi, pointj, refSize)<UNIQUE) { obj->pointNb(i)->field.uniqueness++; } }

			}
		}
		if(obj->pointNb(i)->field.uniqueness > uMax) { uMax=obj->pointNb(i)->field.uniqueness; }
		if(obj->pointNb(i)->field.uniqueness < uMin) { uMin=obj->pointNb(i)->field.uniqueness; }
	}

	thres=(float)uMin+ambThres*(float)(uMax-uMin); //update ambiguity threshold for the mesh
}



//--------------------------------------------------------------------------------------------

//retrieve 1-ring neighbors of p in obj, and output in candidateList (first element is p index)
void match::retrieve1stRing(Object* obj, Point_3* p, list<int> &candidateList) const
{
    int indexP=obj->getPointIndex(p);
    candidateList.push_back(indexP);
    
    //retrieve 1-ring neighbor vertices to p
    list<edge *> neighbor=p->edges;
    list<edge *>::iterator it=neighbor.begin();
    list<edge *>::iterator itEnd=neighbor.end();
    for(; it!=itEnd; it++)
    {
        Point_3* neighborPoint=(*it)->otherPoint(p);
        if(neighborPoint != NULL)
        {
            candidateList.push_back(obj->getPointIndex(neighborPoint));
        }
    }

}

//retrieve kNN of p in obj, and output in candidateList (first element is p index)
void match::retrieveNeighbors(Object* obj, Point_3* p, list<int> &candidateList, uint k) const
{
    int nbPoints=obj->getNbPoints();
    int indexP=obj->getPointIndex(p);
    
    //candidateList is the output
    candidateList.push_back(indexP);

    list<int> neighborList;
    neighborList.push_back(indexP);
                            
    Point_3* visitedPoint=NULL;
    Point_3* neighborTmp=NULL;
    int neighborIndexTmp;
    
//    cout << ">> Retrieve the k nearest neighbors of " << indexP << endl;
    while(candidateList.size()<k)
    {
        visitedPoint=obj->pointNb(neighborList.front());
        
        if(visitedPoint!=NULL && !visitedPoint->edges.empty())
        {
            list<edge *> adjacentEdges=visitedPoint->edges;
            list<edge *>::iterator it=adjacentEdges.begin();
            list<edge *>::iterator itEnd=adjacentEdges.end();
            for(; it!=itEnd; it++)
            {
                //adjacent vertex to visitedPoint
                neighborTmp=(*it)->otherPoint(visitedPoint);
                if(neighborTmp!=NULL)
                {
                    //insert the neighbor point in the candidate list if it has not been inserted yet
                    neighborIndexTmp=obj->getPointIndex(neighborTmp);
                    if(neighborIndexTmp<nbPoints)
                    {
                        list<int>::iterator itCand=candidateList.begin();
                        bool contain=false;
                        for(; itCand!=candidateList.end(); itCand++)
                        {
                            if((*itCand) == neighborIndexTmp)
                            {
                                contain=true;
                                break;
                            }
                        }
                    
                        if(contain==false)
                        {
                            neighborList.push_back(neighborIndexTmp);
                            candidateList.push_back(neighborIndexTmp);
                            //cout << " " << neighborIndexTmp;
                        }
                    }
                }
            }
        }

    	candidateList.unique();    
        neighborList.pop_front();
    }

}

float match::energy(std::vector<int> &refPoints1, std::vector<int> &refPoints2, std::vector<int> &refIndex) const
{
	float energy=0.0;
	vector<int>::iterator it=refIndex.begin(), itEnd=refIndex.end();
	for(; it!=itEnd; it++)
	{
		vector<int>::iterator it2=refIndex.begin(), itEnd2=refIndex.end();
		for(; it2!=itEnd2; it2++)
		{
			energy+=fabs(obj1->pointNb(refPoints1[*it])->field.value[*it2]-obj2->pointNb(refPoints2[*it])->field.value[*it2])
				*fabs(obj1->pointNb(refPoints1[*it])->field.value[*it2]-obj2->pointNb(refPoints2[*it])->field.value[*it2]);

		}
	}
	return(energy);
}

//optimize refPoints1 positions
void match::minEnergy()
{
	int size1=refPoints1.size();
	vector<int> newRefPoints1; newRefPoints1.reserve(size1);
	for (int i=0; i<size1; i++) { newRefPoints1.push_back(refPoints1[i]); }

	int nbRefIndex=size1; if(size1>ggRefNb) { nbRefIndex=ggRefNb; }	
	// optimize the refPoints1 positions sequentially (not the first five ones)
	for (int i=6; i<size1; i++)
	{

		vector<int> refIndex2; refIndex2.reserve(nbRefIndex);
		getClosestRefPoints(obj2->pointNb(refPoints2[i]), refPoints2, nbRefIndex, -1, refIndex2);

		// actual energy before optimization
		float minEnergy=energy(newRefPoints1, refPoints2, refIndex2);
		int minPoint=newRefPoints1[i];

		// retrieve the neighborhood of the visited refPoint1
		list<int> neighbor;
		retrieveNeighbors(obj1, obj1->pointNb(newRefPoints1[i]), neighbor, 25); //12,25
		//retrieve1stRing(obj1, obj1->pointNb(newRefPoints1[i]), neighbor);

		list<int>::iterator it=neighbor.begin(), itEnd=neighbor.end();
		for(; it!=itEnd; it++)
		{
			if(*it < nbPoints1_0)
			{
				// have to recalculate the i-th geodesic coordinate wrt *it
				updateGeodesicCoordinate(obj1, i, *it);
				newRefPoints1[i]=*it;
				float tmpE=energy(newRefPoints1, refPoints2, refIndex2);
				if(tmpE<minEnergy) { minEnergy=tmpE; minPoint=*it; }
			}
		}
		
		// update the i-th refPoints with the neighbor that minimize the most the energy
		newRefPoints1[i]=minPoint;
		updateGeodesicCoordinate(obj1, i, minPoint);
	}
	
	refPoints1.clear();
	for(int i=0; i<size1; i++) { refPoints1.push_back(newRefPoints1[i]); }
	newRefPoints1.clear();
}

// get the nbRefIndex closest refPoints among the firstRef refPoints to a point: return a list of refPoints indices
void match::getClosestRefPoints(point_3* point, vector<int> &refPoints, int nbRefIndex, int firstRef, vector<int> &refIndex)
{
	int refSize=refPoints.size();
	if (firstRef<0) firstRef=refSize;

	// retrieve and sort geodesic distances to all refPoints 
	list<valRank> refInd;
	for(int i=0; i<firstRef; ++i)
	{
		valRank tmpRefInd(i, point->field.value[i]);
		refInd.push_back(tmpRefInd);
	}
	refInd.sort();

	// keep the nbRefIndex closest ones and store indices in refIndex	
	list<valRank>::iterator it=refInd.begin(), itEnd=refInd.end();
	for(int i=0; it!=itEnd, i<nbRefIndex; it++, i++) { refIndex.push_back(it->index); }
}


// output reference points on the mesh
void match::outputRefPointsOnMesh(string outMesh, Object* obj, vector<int> &refPoints) const
{
	cout << "Writing " << outMesh << " ..." << endl;
	ofstream outMeshFile(outMesh.c_str(), ios::out);
	if (!outMeshFile.is_open())  { cout << "Error writing COFF file " << outMesh << endl; return; }
	else
	{
		int nbPoints=obj->getNbPoints();
        	int nbFaces=obj->getNbFaces();
        	int nbEdges=obj->getNbEdges();
		int nbRefPoints=refPoints.size();

		outMeshFile << "COFF" << endl;
		outMeshFile << nbPoints << " " << nbFaces << " " << nbEdges << endl;
		outMeshFile << "#vertex" << endl;
		
		colorMapFunc cmf=ColorMap::selectColorMap(4); //4: cyclic colormap
		for(int i=0; i<nbPoints; ++i)
		{
			// if the visited vertex is in the neighborhood of a reference point, give it a color
			unsigned char color[3]={255, 255, 255}; int k=0; bool found=false;

			
			list<int> candidateList;
			retrieve1stRing(obj, obj->pointNb(i), candidateList); //retrieve the 1-ring of point(i). 1st ring = 6 vertices	
			//retrieveNeighbors(obj, obj->pointNb(i), candidateList, 10);

			list<int>::iterator it=candidateList.begin(), itEnd=candidateList.end();
			for(; it!=itEnd; it++)
			{
				for(int j=0; j<nbRefPoints; ++j) { if(refPoints[j]==(*it)) { cmf(color, j, 0, refPoints.size()); k=j; found=true; break;} }
				if(found) break;	
			}

			// write 3D coordinates and rgb color 
			outMeshFile << obj->pointNb(i)->x << " " << obj->pointNb(i)->y << " " << obj->pointNb(i)->z << " "
			<< float(color[0])/255 << " " << float(color[1])/255 << " " << float(color[2])/255 << " " << k << endl;
		}

		outMeshFile << "#face" << endl;
		for(int j=0; j<nbFaces; ++j)
		{
			outMeshFile << "3 " <<
			obj->getPointIndex(obj->faceNb(j)->p[0]) << " " << obj->getPointIndex(obj->faceNb(j)->p[1]) << " " << obj->getPointIndex(obj->faceNb(j)->p[2]) << endl;
		}
	}

	outMeshFile.close();
}


// output isovalue lines on mesh1 that correspond to a reference point index
void match::outputIsoLinesOnMesh(string outMesh, Object* obj, std::vector<int> &refPoints, int index) const
{
	cout << "Writing "<< outMesh << "..." << endl;
	ofstream outMeshFile(outMesh.c_str(), ios::out);
	if(!outMeshFile.is_open()) { cout << "Error writing COFF file " << outMesh << endl; return; }
	else
	{
		//get the geodesic coordinates of point index
		int nbRefPoints1=refPoints1.size();
		vector<float> value2; value2.reserve(nbRefPoints1);
	    	for(int j=0; j<nbRefPoints1; j++) { value2[j]=obj1->pointNb(index)->field.value[j]; printf("value1[%d]=%2.3f\n", j, value2[j]); }

		//write outMesh1
		int nbPoints1=obj1->getNbPoints();
	        int nbFaces1=obj1->getNbFaces();
	        int nbEdges1=obj1->getNbEdges();

		outMeshFile << "COFF" << endl;
		outMeshFile << nbPoints1 << " " << nbFaces1 << " " << nbEdges1 << endl;
		outMeshFile << "#vertex" << endl;
		
		colorMapFunc cmf=ColorMap::selectColorMap(4); //4: cyclic colormap
		for(int i=0, k=0; i<nbPoints1; ++i)
		{
			// if the visited vertex has a coordinate component that is on an isoline value
			// give the color of the corresponding reference point to the isoline
			unsigned char color[3]={255, 255, 255}; k=0;
			for(int j=0; j<nbRefPoints1; ++j)
			{
				if((refPoints1[j]==i) || (fabs(obj1->pointNb(i)->field.value[j]-value2[j])<0.01))
				{
					cmf(color, j, 0, refPoints1.size()); /*color[0]=255; color[1]=0; color[2]=0;*/
					k=j; break;
				}
			}

			// write 3D coordinates and rgb color 
			outMeshFile << obj1->pointNb(i)->x << " " << obj1->pointNb(i)->y << " " << obj1->pointNb(i)->z << " "
			<< float(color[0])/255 << " " << float(color[1])/255 << " " << float(color[2])/255 << k << endl;
		}

		outMeshFile << "#face" << endl;
		for(int j=0; j<nbFaces1; ++j)
		{
			outMeshFile << "3 " <<
			obj1->getPointIndex(obj1->faceNb(j)->p[0]) << " " << obj1->getPointIndex(obj1->faceNb(j)->p[1]) << " " << obj1->getPointIndex(obj1->faceNb(j)->p[2]) << endl;
		}
	}

	outMeshFile.close();
}


void match::SKELdump(string outSKEL) const
{
	int refSize1=refPoints1.size();
	int refSize2=refPoints2.size();
	cout << "Writing " << outSKEL << " ..." << endl;
	ofstream outFile(outSKEL.c_str(), ios::out);
	if(!outFile.is_open()) { cout << "Error creating " << outSKEL << endl; }
	outFile << "SKEL" << endl;
	outFile << refSize1+refSize2 << " " << refSize1 << endl;
	for(int i=0; i<refSize1; i++)
	{
		outFile << obj1->pointNb(refPoints1[i])->x << " "<< obj1->pointNb(refPoints1[i])->y << " " << obj1->pointNb(refPoints1[i])->z << endl;
	}
	for(int i=0; i<refSize2; i++)
	{
		outFile << obj2->pointNb(refPoints2[i])->x << " "<< obj2->pointNb(refPoints2[i])->y << " " << obj2->pointNb(refPoints2[i])->z << endl;
	}
	for(int i=0; i<refSize1; i++)
	{
		outFile << "2 " << i << " " << i+refSize1 << endl;
	}
	outFile.close();
}

void match::matchMesh()
{
	cout << "Let's align the surfaces!" << endl;
	int refSize2=refPoints2.size();
	int nbPoints1=nbPoints1_0; //obj1->getNbPoints(); // use non-subdivided mesh1
	int nbPoints2=nbPoints2_0; //obj2->getNbPoints(); // use non-subdivided mesh2
	vector<int> matchVec; matchVec.reserve(nbPoints2);
	for(int i=0; i<nbPoints2; i++){ matchVec[i]=-1; }
	for(int i=0; i<refSize2; i++) { matchVec[refPoints2[i]]=refPoints1[i]; }

	// visit all vertices of Mesh2
	int unmatched=nbPoints2-refSize2;
	int prev_unmatched=0;
	int unmatched_s=-1;
	float flowThres=0.1;
	float alpha=0.8, beta=1-alpha;
	float mult=0.8;
	int iter=3;
	int bs=sizeof(int); //bin size
	vector<int> bitvector(nbPoints1/bs+1, 0);

	while(unmatched>0)
	{
		while(unmatched>0)
		{
			for(int i=0; i<nbPoints2; i++)
			{
				if(matchVec[i]>=0) continue; //already matched
				point_3* point2=obj2->pointNb(i);
				vector<int> refIndex2; refIndex2.reserve(5);
				getClosestRefPoints(point2, refPoints2, 5, -1, refIndex2);

				list<int> candidateList;
				list<int> visitList;
				retrieveNeighbors(obj2, point2, candidateList, 12); //12, 25, 49
				//retrieve1stRing(obj2, point2, candidateList);

				//get closest matched points and compute average of flows from existing matches
				Point_3 meanRefFlowDst=Point_3(0.0, 0.0, 0.0);
				Point_3 meanRefFlowSrc=Point_3(0.0, 0.0, 0.0);
				float flowNb=0;
				if(!candidateList.empty())
				{
				    list<int>::iterator candIt=candidateList.begin(), candItEnd=candidateList.end();
				    for(; candIt!=candItEnd; candIt++)
				    {
					if(*candIt<nbPoints2)
					{
						Point_3* neighborPoint=obj2->pointNb(*candIt);
						if(neighborPoint != NULL && matchVec[*candIt] != -1) // visit a neighbor point in Mesh2 (that has been matched)
						{
						    if(matchVec[*candIt]<nbPoints1) // check if the indices of the pair is valid
						    {
							Point_3* neighborPoint1=obj1->pointNb(matchVec[*candIt]); // get the corresponding neighbor point in Mesh1
							if(neighborPoint1 != NULL)
							{
							    list<int> visitListTmp;
							    retrieveNeighbors(obj1, neighborPoint1, visitListTmp, 12); //12,25
							    visitListTmp.sort();
							    visitList.merge(visitListTmp); //merge items from visitlistTmp to visitList
							    visitList.unique(); //keep only unique elements

							    meanRefFlowDst.x+=neighborPoint->x;
							    meanRefFlowDst.y+=neighborPoint->y;
							    meanRefFlowDst.z+=neighborPoint->z;
							    meanRefFlowSrc.x+=neighborPoint1->x;
							    meanRefFlowSrc.y+=neighborPoint1->y;
							    meanRefFlowSrc.z+=neighborPoint1->z;
							    flowNb++;
							}
						    }
/*						    else
						    {
							cout << "The pair: *candIt=" << *candIt << "/" << nbPoints2 <<
								", matchVec[*candIt]=" << matchVec[*candIt] << "/" << nbPoints1 << " does not exist!" << endl;
						    }
*/						}
					}
				    }
				}

				int closest=-1;
				if(flowNb>0)
				{
					//find best matching point in Mesh1 among the visitList			
					float minDist=1e9, ggDist=1e9;
					//for(int j=0; j<nbPoints1; j++)
				
					list<int>::iterator visitIt=visitList.begin(), visitItEnd=visitList.end();
					for(; visitIt!=visitItEnd; visitIt++)
					{
						if(*visitIt<nbPoints1)
						{
							point_3* point1=obj1->pointNb(*visitIt); //j

							ggDist=ggDistance(point1, point2, refIndex2); //global geodesic distance
							if(ggDist > 10) continue;
			    
			    				Point_3 meanFlow=(meanRefFlowSrc-meanRefFlowDst)/flowNb;
							float flowDiff=((*point1-*point2)-meanFlow).norm();
							if(flowDiff > flowThres) continue; //check candidate consistency with the meanFlow

							// will choose the most geodesically consistent
							float tmpDist=alpha*ggDist+beta*flowDiff;
							if(tmpDist<minDist)
							{
								minDist=tmpDist;
								closest=*visitIt; //j
							}
						}
					}
				}

				if(closest != -1)
				{
					matchVec[i]=closest;
					bitvector[closest/bs] |= (1 << (closest % bs));
					unmatched--;
					//cout << "ggDist= " << ggDist << endl;
				}
			}

			cout << "Unmatched=" << unmatched << endl;

			if(unmatched == prev_unmatched)
			{
				flowThres*=2; //relax the constraint if there is no improvement
			}

			prev_unmatched=unmatched;
		}

		//ensure "complete" surface matching (of mesh1) by going through the unmatched vertices in mesh1
		long int unmatched1=0;
		flowThres=0.1;
		for(int j=0; j<nbPoints1; j++)
		{
			if((bitvector[j/bs] & (1 << (j % bs))) == 0) //not matched
			{
				point_3* point1=obj1->pointNb(j); //get the unmatched point
				vector<int> refIndex1; refIndex1.reserve(5);
				getClosestRefPoints(point1, refPoints1, 5, -1, refIndex1);

				//find best matching point in Mesh2
				float minDist=1e9;
				int closest=-1;
				for(int i=0; i<nbPoints2; i++)
				{
					point_3* point2=obj2->pointNb(i);

					// will choose the most geodesically consistent
					float ggDist=ggDistance(point1, point2, refIndex1);
					if(ggDist>10) continue;

					list<int> candidateList;
					retrieveNeighbors(obj2, point2, candidateList, 12); //12, 25, 49
					//retrieve1stRing(obj2, point2, candidateList);
	
					//get closest matched points and compute average of flows from existing matches
					Point_3 meanRefFlowDst=Point_3(0.0, 0.0, 0.0);
					Point_3 meanRefFlowSrc=Point_3(0.0, 0.0, 0.0);
					float flowNb=0;
	    
					if(!candidateList.empty())
					{
					    list<int>::iterator candIt=candidateList.begin(), candItEnd=candidateList.end();
					    for(; candIt!=candItEnd; candIt++)
					    {
						if(*candIt<nbPoints2 )
						{
							Point_3* neighborPoint=obj2->pointNb(*candIt);
							if(neighborPoint != NULL && matchVec[*candIt] != -1) // visit a neighbor point in Mesh2 (that has been matched)
							{
							    if(matchVec[*candIt]<nbPoints1) // check if the indices of the pair is valid
							    {
								Point_3* neighborPoint1=obj1->pointNb(matchVec[*candIt]); // get the corresponding neighbor point in Mesh1
								if(neighborPoint1 != NULL)
								{
								    meanRefFlowDst.x+=neighborPoint->x;
								    meanRefFlowDst.y+=neighborPoint->y;
								    meanRefFlowDst.z+=neighborPoint->z;
								    meanRefFlowSrc.x+=neighborPoint1->x;
								    meanRefFlowSrc.y+=neighborPoint1->y;
								    meanRefFlowSrc.z+=neighborPoint1->z;
								    flowNb++;
								}
							    }
/*							    else
							    {
								cout << "The pair: *candIt=" << *candIt << "/" << nbPoints2 <<
									", matchVec[*candIt]=" << matchVec[*candIt] << "/" << nbPoints1 << " does not exist!" << endl;
							    }
*/							}
						}
					    }
					}

					if(flowNb>0)
					{	    
		    				Point_3 meanFlow=(meanRefFlowSrc-meanRefFlowDst)/flowNb;
						float flowDiff=((*point1-*point2)-meanFlow).norm();
						//cout << ggDist << ":" << flowDiff << endl; //same range
						//if(flowDiff > flowThres) continue; //check candidate consistency with the meanFlow
				
						float tmpDist=alpha*ggDist+beta*flowDiff;
						if(tmpDist<minDist)
						{
							minDist=tmpDist;
							closest=i;
						}
					}
					else
					{
						float tmpDist=ggDist; //ggDist and flowDiff are in the same range [0,3]
						if(tmpDist<minDist)
						{
							minDist=tmpDist;
							closest=i;
						}
					}
				}			

				if(closest !=-1)
				{				
					matchVec[closest]=j; //update
				}
				unmatched1++;
			}
		}

		bitvector.clear();
		cout << unmatched1 << " unmatched vertices in Mesh1 were handled" << endl;

		iter--;
		if(iter==0) break;

		//smoothness constraint: each vertex should comply with its neighbors. if not, we remove both for re-iteration! ultimately it should converge to a smooth mapping~
		for(int i=0; i<nbPoints2; i++)
		{
			point_3* point2=obj2->pointNb(i);
			if(matchVec[i] == -1) continue;
			point_3* point1=obj1->pointNb(matchVec[i]);
			list<int> neighList;
			retrieve1stRing(obj2, point2, neighList);
			if(!neighList.empty())
			{
				list<int>::iterator neighIt=neighList.begin(), neighItEnd=neighList.end();
				neighIt++; //skip first as it is point2
				for(; neighIt!=neighItEnd; neighIt++)
				{
					if(*neighIt<nbPoints2 && matchVec[*neighIt] != -1 && matchVec[*neighIt]<nbPoints1) //visit a neighbor point in Mesh2 that has been matched
					{
						Point_3* neighborPoint2=obj2->pointNb(*neighIt); // get the corresponding neighbor point in Mesh1
						Point_3* neighborPoint1=obj1->pointNb(matchVec[*neighIt]); // get the corresponding neighbor point in Mesh1
						if(neighborPoint1 != NULL)
						{
							Point_3 tmpFlow=*neighborPoint1-*neighborPoint2;
							float flowDiff=((*point1-*point2)-tmpFlow).norm();
							//cout << flowDiff << endl; //~0.05
							if(flowDiff > flowThres*mult) //check consistency with neighbor matched pair: reduce threshold to add more iterations!
							{
								if (matchVec[i] != -1) { matchVec[i] = -1; unmatched++; } //reset!
								if (matchVec[*neighIt] != -1) { matchVec[*neighIt] = -1; unmatched++; } //reset!
							}
						}
					}
				}
			}
		}
		cout << "Unmatched=" << unmatched << " after smoothness check"<< endl; //number of unmatched vertices after smoothing check

		if(unmatched == unmatched_s)
		{
			mult*=1.5;  //relax the constraint if there is no improvement
		}
		unmatched_s=unmatched;
	}


	ofstream lOutFile("output/reMesh.off", ios::out);
	if(lOutFile.is_open() == 0) { cout << "Error creating reMesh.off" << endl; return; }

	ofstream lOutFile2("output/matchVec.txt", ios::out);
	if(lOutFile2.is_open() == 0) { cout << "Error creating matchVec.off" << endl; return; }
	
	lOutFile << "OFF" << endl;
	lOutFile << nbPoints2_0 << " " << nbFaces2_0 << " " << nbEdges2_0 << endl;
	
	for(int i=0; i<nbPoints2_0; i++)
	{
		point_3* aPoint=obj1->pointNb(matchVec[i]); //matched vertex on Mesh1
		lOutFile << aPoint->x << " " << aPoint->y << " " << aPoint->z << endl;
		lOutFile2 << i << " " << matchVec[i] << endl;
	}

	for(int i=0; i<nbFaces2_0; i++)
	{
		face* lFace=obj2->faceNb(i); //connectivity of Mesh2
		lOutFile << 3 << " " << obj2->getPointIndex(lFace->p[0]) << " " << obj2->getPointIndex(lFace->p[1]) << " " << obj2->getPointIndex(lFace->p[2]) << endl;
	}
	
	lOutFile.close();
	lOutFile2.close();
	
	cout << "reMesh.off created!" << endl;
	cout << "matchVec.txt created!" << endl;
}



