#ifndef __SIMILARITY_H_
#define __SIMILARITY_H_

#include "mrgraph.h"
#include "matching.h"
#include "list.h"

float computeSimilarity(MRGNodeAttribute* const &attr0,MRGNodeAttribute* const &attr1, int attribut, int ResX, float factor[NbAttributes], int cumul, int dist);
float computeSimilarity(matchingEl* const &mat0,matchingEl* const &mat1, int attribut, int ResX, float factor[NbAttributes], int cumul, int dist);
//ttung
float computeImp(MRGNodeAttribute* const &attr0);
float computeImp(matchingEl* const &mat0);

matchingEl* findMaxSimEl(std::list<matchingEl*> NLists[2],int &index);
bool sameMuRange(matchingEl* const &m,matchingEl* const &n);
bool sameResolution(matchingEl* const &m,matchingEl* const &n);
bool parentsMatch(matchingEl* m,matchingEl* n);
bool sameMList(matchingEl* const &m,matchingEl* const &n);
bool parentsListMatch(const std::list<MRGNode*> &nodes1, const std::list<MRGNode*> &nodes2);
bool neighboursParentsMatch(matchingEl* const &m,matchingEl* const &n);
bool topoConsistent(matchingEl* const &m,matchingEl* const &n,const bool mask[4]);
std::list<matchingEl*> findCandidates(matchingEl* const &m,const std::list<matchingEl*> &NList,const bool mask[4]);

float computeLoss(MRGNodeAttribute* const &m,MRGNodeAttribute* const & n, float factor[NbAttributes], int res, int cumul, int dist);
float computeLossAdjacent(matchingEl* const &m,matchingEl* const &n,Towards direction, float factor[NbAttributes], int res, int cumul, int dist);
float computeMat(matchingEl* const &m,matchingEl* const &n, float factor[NbAttributes], int res, int cumul, int dist);

void testNodeSetFromSource(matchingEl* const& m,MRGNode* const & source,std::list<matchingEl*> const & candidates,float &bestMat,matchingEl* &result,Towards direction, float factor[NbAttributes], int res, int cumul, int dist);
matchingEl* selectBestCandidate(matchingEl* const &m,std::list<matchingEl*> const &candidates, float factor[NbAttributes], int res, int cumul, int dist);
bool findMatchingPair(std::list<matchingEl*> NLists[2],matchingEl* &mat0,matchingEl* &mat1,const bool mask[4], float factor[NbAttributes], int res, int cumul, int dist);
void insertFinerChilds(matchingEl* const &n,std::list<matchingEl*> &l);

//ttung
void computeSimilarityOneTurn(MRG* mrgs[2],PairList &MPair,float &result,int &cptLabel,const bool mask[4],std::list<matchingEl*> activeLists[2],std::list<matchingEl*> failureLists[2], float factor[NbAttributes], int res, int attribut, int cumul, int dist);
float *computeSimilarity(MRG* mrgs[2],PairList &MPair, float factor[NbAttributes], int res, int cumul, int dist);
float *computeSimilarity(MRG* const &mrg1,MRG* const &mrg2,PairList &MPair, float factor[NbAttributes], int res, int cumul, int dist);

PairList findMatchingPairs(PairList const& MPair,int res);
void writeMatchingGraph(char* file,MRG* const &mrg1,MRG* const &mrg2,PairList const &MPair,int resLvl);

//ttung:
float computeSimilaritySubMRG(MRGNodeAttribute* const &attr0,MRGNodeAttribute* const &attr1, float factor[NbAttributes], int res, int cumul, int dist);
float computeSimilarityOneNode(MRGNodeAttribute* const &attr0,MRGNodeAttribute* const &attr1, float factor[NbAttributes], int res, int cumul, int dist);
void computeSimilarityOneTurnSubMRG(MRG* mrgs[2],PairList &MPair,float &result,int &cptLabel,const bool mask[4],std::list<matchingEl*> activeLists[2],std::list<matchingEl*> failureLists[2], float factor[NbAttributes], int res, int cumul, int dist);

void additionalTopologyMatching(MRG* mrg1, MRG* mrg2, PairList &MPair, int res);
void computeHumanSimilarity(MRG* mrg1, MRG* mrg2, PairList &MPair, int res);
float computeSIMSubMRG(MRG* mrgs[2],PairList &MPair, float factor[NbAttributes], int res, int cumul, int dist);     
float computeSIMSubMRG(MRG* const &mrg1,MRG* const &mrg2,PairList &MPair, float factor[NbAttributes], int res, int cumul, int dist); 


List<MRGNode*>* getSimilarNodes(MRG* mrgs[2], int res, int index, float threshold, float factor[NbAttributes], int resX, int cumul, int dist);
List<MRGNode*>* getSimilarNodes(MRG* const &mrg1,MRG* const &mrg2, int res, int index, float threshold, float factor[NbAttributes], int resX, int cumul, int dist);

List<MRGNode*>* getSimilarNodesFromSubMRG(MRG* mrgs[2], int res0, int res1, float threshold, float factor[NbAttributes], int resX, int cumul, int dist);
List<MRGNode*>* getSimilarNodesFromSubMRG(MRG* const &mrg0,MRG* const &mrg1, int res0, int res1, float threshold, float factor[NbAttributes], int resX, int cumul, int dist);


#endif
