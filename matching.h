#ifndef __MATCHING_H_
#define __MATCHING_H_

#include "mrgraph.h"
#include "list.h"
#include <string.h>

const float weight=0.75f; // pour a et l
const float weight2=0.5f; // pour Gdistance et Ndistance

/*>>>>>>>>>> classe matchingEl <<<<<<<<<<
 *
 * matchingEl = matching element
 * element qui est matche dans le calcul de similarite
 * ie un noeud ou plusieurs noeuds
 *
 */
enum matType {OneNode,SeveralNodes};

void countSeveral(int &severalNo);

class matchingEl
{
public:
     matType type;     
     
     // utilis'e si type==OneNode
     MRGNode* oneNode;

     // utilis'e si type==SeveralNodes
     MRGNode** severalNodes;
     int cptNodes;
     int severalNo; // permet de faire des comparaisons sans comparer les tableaux element par element

     // matchingEl parent de celui-ci
     class matchingEl* parentEl;
       
     // matchingEls qui sont a l'origine de celui-ci (utile si on cree un nouveau matchingEl en en concatenant plusieurs)
     class matchingEl** originalElArray;
     int nbOriginalEl;
     
     // matchingEl qui matche celui-ci
     matchingEl* coMatchingEl; 
     
     // attributs de ce matchingEl
     MRGNodeAttribute* attr;

     float importance; // =sim(this,this)

     void commonInit(void);
     matchingEl(MRGNode* _oneNode,matchingEl* parent);
     matchingEl(const std::list<MRGNode*> &listNodes,matchingEl* parent);
     matchingEl(const std::list<matchingEl*> &listEl,int sz); 
     ~matchingEl(void);

	 std::list<MRGNode*> getAdjacentNodes(Towards direction);
     MRGNode** getChilds(int &sz);
     int nbNodes(void);
     void printCoords(MRG* graph);
     MRGNode* getNodeNb(int nb);
     void removeOriginalElFromList(std::list<matchingEl*> &l);
     void putOriginalElInList(std::list<matchingEl*> &l);
     void deleteOriginalEl(void);
     void propagateInAllNodes(int label);
     MRGNodeAttribute* getAttribute(void) const { return attr; }
     muRange getRange(void) const;
     void markNodesWithMatching(void);
};
bool sameMatchingEl(matchingEl* m1,matchingEl* m2);
std::list<matchingEl*> NodeArray2MatchingElList(MRGNode** array,int sz,matchingEl* parent);
void theseMatchingElMatch(matchingEl* m0,matchingEl* m1);

class matchingPair
{
private:
  matchingEl* mat[2];
public:
  matchingPair(void) {
    mat[0] = NULL;
    mat[1] = NULL;
  }
  matchingPair(matchingEl* const &mat0,matchingEl* const &mat1) {
    mat[0] = mat0;
    mat[1] = mat1;
  }

  matchingEl* getMatchingEl(int index) {
    return mat[index];
  }
};
typedef std::list<matchingPair*> PairList;

void pairListInsert(matchingEl* const &node0,matchingEl* const &node1,PairList &l);
int nbEdgesInList(PairList MPair);
void deleteThis(PairList &MPair);

#endif
