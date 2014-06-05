#ifndef __LIST_H_
#define __LIST_H_

#include <stdlib.h>
#include <list>
#include <algorithm>

template <class O>
class List;

template<class O>
int indexBiggestElement(List<O>* const & list)
{
  if (list==NULL)
    return -1;
  else
    {
      List<O>* l;
      O element=list->getHead();
      int result=0;
      int cpt=1;
      for(l=list->getTail();l!=NULL;l=l->getTail())
	{
	  if(element<l->getHead())
	    {
	      element=l->getHead();
	      result=cpt;
	    }
	  ++cpt;
	}
      return result;
    }
}

template<class O>
int listLength(List<O>* const &l)
{
  if (l==NULL)
    return 0;
  else 
    return 1+listLength(l->tail);
};

template<class O>
void listInsert(const O hd,List<O>* &tl)
{
  List<O>* newList=new List<O>;
  newList->head=hd;
  newList->tail=tl;
  if (tl) newList->size=tl->size+1;
  else newList->size=1;
  tl=newList;
}

template<class O>
List<O> * listAppend(const O &hd,List<O> * const &tl) 
{
  if(tl==NULL)
    return new List<O>(hd,NULL);
  else
    return new List<O>(tl->head,listAppend(hd,tl->getTail()));
}

/*
 * listReplace: Attention, ne remplace que la premiere occurence trouvee
 */
template<class O>
void listReplace(const O old,const O _new,List<O>* &l)
{
  if (l==NULL)
      printf("!!! Attention : Aucun remplacement effectue dans la liste\n");
  else
    {
      if (l->head==old)
	l->head=_new;
      else
	listReplace(old,_new,l->tail);
    }
}

template<class O>
void listRemove(const O old,List<O>* &l)
{
  if(l==NULL)
    printf("!!! Attention : Aucune suppression effectuee dans la liste\n");
  else
    {
      if(l->getHead()==old)
	{
       List<O>* sauv=l->getTail();
      l->tail = NULL;
      delete l;
//      free (l);
      l=sauv;
    }
    else
      listRemove(old,l->tail);
    }
}

template<class O>
void listConcat(List<O>* &l,List<O>* const &ll)
{
  List<O>* plus=ll;
  while(plus!=NULL)
  {
    listInsert(plus->getHead(),l);
    plus=plus->getTail();
  }
}

template<class O>
O* getElementNb(int index,List<O*>* const &l)
{  
  if (l==NULL)
    return NULL;
  if (index==0)
    return l->getHead();
  return getElementNb(index-1,l->getTail());
}

template<class O>
bool isMember(const O thing,List<O>* const &l)
{
  if(l==NULL)
    return false;
  if(thing==l->getHead())
    return true;
  return isMember(thing,l->getTail());
}

template<class O>
bool commonElement(List<O>* const &l1,List<O>* const &l2)
{
  if(l1==NULL)
    return false;
  if (isMember(l1->getHead(),l2))
    return true;
  return commonElement(l1->getTail(),l2);
}

template<class O>
bool sameLists(List<O>* const &l1,List<O>* const &l2)
{
  if ((l1==NULL)&&(l2==NULL))
    return true;
  if ((l1==NULL)||(l2==NULL))
    return false;
  return (l1->getHead()==l2->getHead())&&(sameLists(l1->getTail(),l2->getTail()));
}

template<class O>
int getPosition(List<O>* const &l,O el)
{
  if(l==NULL)
    {
      printf("!!! Attention : element not found in the list\n");
      return -1;
    }
  if (el==l->head)
    return 0;
  else
    return (1+getPosition(l->tail,el));
}

/*
 * allSubSets :
 * Retourne une liste de tous les sous-ensembles (sauf le sous-ensemble NULL)
 */
template<class O>
List<List<O>*>* allSubSets(List<O>* l)
{
  if(l==NULL)
    return NULL;
  else {
    List<List<O>*>* suite,* res=NULL;
    for(suite=allSubSets(l->getTail());suite!=NULL;suite=suite->getTail()) {
	  listInsert(new List<O>(l->getHead(),suite->getHead()),res);
	  listInsert(suite->getHead(),res);
	}
    listInsert(new List<O>(l->getHead(),NULL),res);

    return res;
  }
}

template<class O>
void listConcatNoDoublon(List<O>* &l,List<O>* plus)
{
     O truc;
     while(plus!=NULL)
     {
          truc=plus->getHead();
          if(isMember(truc,l)==false)
               listInsert(truc,l);
          plus=plus->getTail();
     }
}

/*
// ttung
template<class O>
void replaceListByThisOne(List<O>* &l1, List<O>* &l2)
{
	l1->setHead(l2->getHead());
	l1->setTail(l2->getTail());
	l1->setSize(l2->size);
}
*/

template<class O>
List<List<O>*>* allSubSets(std::list<O> l)
{
  if(l.empty())
    return NULL;
  else {
    List<List<O>*>* suite,* res=NULL;
	std::list<O> l_temp=l;
	l_temp.pop_front();
    for(suite=allSubSets(l_temp);suite!=NULL;suite=suite->getTail()) {
	  listInsert(new List<O>(l.front(),suite->getHead()),res);
	  listInsert(suite->getHead(),res);
	}
    listInsert(new List<O>(l.front(),NULL),res);

    return res;
  }
}


template<class O>
std::list<std::list<O> > allSubSets2 (std::list<O> l)
{
	if(l.empty()) {
		std::list<std::list<O> > res;
		return res;
	}    
  else {
	std::list<std::list<O> > suite, res;
	std::list<O> l_temp=l;

	int sub_size=5;
	if (l.size()> 6) sub_size=1;

	l_temp.pop_front();
    for(suite=allSubSets2(l_temp);!suite.empty();suite.pop_front()) {
	  std::list<O> a_temp;
	  a_temp=suite.front();
	  a_temp.push_front(l.front());

	  if (a_temp.size() > sub_size) continue;

	  res.push_front(a_temp);
	  res.push_front(suite.front());
	}
	std::list<O> b_temp;
	b_temp.push_front(l.front());
    res.push_front(b_temp);
    return res;
  }
}
/*
template<class O>
std::list<List<O>*> allSubSets2 (std::list<O> l)
{
  if(l.empty())
    return NULL;
  else {
    std::list<List<O>*> suite, res;
	std::list<O> l_temp=l;
	l_temp.pop_front();
    for(suite=allSubSets2(l_temp);!suite.empty();suite.pop_front()) {
	  res.push_front(new List<O>(l.front(),suite.front()));
	  res.push_front(suite.front());
	}
    res.push_front(new List<O>(l.front(),NULL));

    return res;
  }
}
*/
/*
template<class O>
List<std::list<O>>* allSubSets2(std::list<O> l)
{
  if(l.empty())
    return NULL;
  else {
	  List<std::list<O>>* suite,* res=NULL;
	std::list<O> l_temp=l;
	l_temp.pop_front();
    for(suite=allSubSets2(l_temp);suite!=NULL;suite=suite->getTail()) {
	  std::list<O> a_temp;
	  a_temp=suite->getHead();
	  a_temp.push_front(l.front());
	  listInsert(a_temp,res);
	  listInsert(suite->getHead(),res);
	}
	std::list<O> b_temp;
	b_temp.push_front(l.front());
    listInsert(b_temp,res);

    return res;
  }
}
*/


template<class O>
bool sameLists(std::list<O> const &l1, std::list<O> const &l2)
{
  if ((l1.empty())&&(l2.empty()))
    return true;
  if ((l1.empty())||(l2.empty()))
    return false;
  if (l1.size()!=l2.size())
    return false;

  std::list<O> l1_temp=l1, l2_temp=l2;
  for (unsigned int i=0; i<l1_temp.size(); i++) {  
	  if (l1_temp.front() != l2_temp.front()) return false;
	  else {
		  l1_temp.pop_front();
		  l2_temp.pop_front();	  
	  }
  }
  return true;
}


template<class O>
void listConcat(std::list<O> &l,std::list<O> const &ll)
{
  l.insert(l.begin(), ll.begin(), ll.end());
//  std::list<O> plus=ll;
//  while(!plus.empty())
//  {
//    l.push_front(plus.front());
//    plus.pop_front();
//  }
}

template<class O>
void listConcatNoDoublon(std::list<O> &l,std::list<O> plus)
{
     O truc;
     while(!plus.empty())
     {
          truc=plus.front();
          if(isMember(truc,l)==false)
               l.push_front(truc);
          plus.pop_front();
     }
}

template<class O>
int getPosition(std::list<O> l,const O el)
{
  if(l.empty())
    {
      printf("!!! Attention : element not found in the list\n");
      return -1;
    }
  int i=0;
  while (el != l.front())
  {
	++i;
	l.pop_front();
  }
  return i;
}


template<class O>
bool isMember(const O thing,std::list<O> const &list)
{
  return (std::find(list.begin(), list.end(), thing) != list.end());
//	if(list.empty()) return false;
//	std::list<O> l;
//	for(l=list;!l.empty();l.pop_front()) if(thing==l.front()) return true;
//	return false;
}

template<class O>
void listRemove(const O old,std::list<O> &l)
{
/*
	if (isMember(old,l)) {
		std::list<O> l_temp;
		while (old != l.front()) {
			l_temp.push_back(l.front());
			l.pop_front();
		}
		l.pop_front();
		while (!l.empty()) {
			l_temp.push_back(l.front());
			l.pop_front();
		}
		l=l_temp;
    }
	else printf("!!! Attention : Aucune suppression effectuee dans la liste\n");
*/
  typename std::list<O>::iterator I = std::find(l.begin(), l.end(), old);
  if (I != l.end()) {
    l.erase(I);
  }
  else printf("!!! Attention : Aucune suppression effectuee dans la liste\n");
}

template<class O>
bool commonElement(std::list<O> const &l1,std::list<O> const &l2)
{
	if(l1.empty()) return false;
	std::list<O> l;
	for(l=l1;!l.empty();l.pop_front()) if(isMember(l.front(),l2)) return true;
	return false;
}



template<class O>
O getElementNb(int nb, std::list<O> l)
{
	int i;
	std::list<O> ltemp=l;
	if (l.empty()) return 0;
	else {
		for (i=0; i<nb; i++) {
			ltemp.pop_front();
		}
		return (ltemp.front());
	}
}


template<class O>
int indexBiggestElement(std::list<O> const & list)
{
	if (list.empty()) return -1;
	else {
		std::list<O> l;
		std::list<O> ltemp=list;
		O element=list.front();
		int result=0;
		int cpt=1;
		ltemp.pop_front();
		for(l=ltemp;!l.empty();l.pop_front()) {
			if(element<l.front()) {
				element=l.front();
				result=cpt;
			}
		++cpt;
		}
		return result;
    }
}

template <class O>
class List 
{  
private:
  O head;
  List<O> * tail;

public:
  int size;
  List<O>(void) 
    {tail=NULL; size=0;};
  
  List<O>(const List<O> &l)
    {
      head=l.head;
      if(l.tail!=NULL)
	tail=new List<O>(*l.tail);
      else
	tail=NULL;
	  size=l.size;
    }

  List<O>(O hd, List<O> * tl) 
    {
      head=hd;
      tail=tl;
	  if (tl) {size=1+tl->size;}
	  else size=1;
    };
 
  ~List<O>(void)
    {
      if (tail!=NULL)
	    delete tail;
      tail=NULL;	  
    };

  void setHead(O hd)
    {head=hd;}

  O getHead(void)
    {return head;};
  
  List<O> * getTail(void)
    {return tail;};

  O getElementNb(int nb)
    {
      if(nb==0)
	return this->getHead();
      if(this->tail==NULL)
	{
	  printf("!!! Warning : Error, too short list");
	  return this->getHead(); // Ce resultat est faux !!!
	}
      else
	return this->tail->getElementNb(nb-1);
    }

#ifdef WIN32
  friend int listLength(List<O>* const &l);
  friend void listInsert(const O hd,List<O>* &tl);
  friend List<O> * listAppend(const O &hd,List<O> * const&tl);
  friend void listReplace(const O old,const O _new,List<O>* &l);
  friend void listRemove(const O old,List<O>* &l);
  friend int getPosition(List<O>* const &l,O el);
#else
  friend int listLength<O>(List<O>* const &l);
  friend void listInsert<O>(const O hd,List<O>* &tl);
  friend List<O> * listAppend<O>(const O &hd,List<O> * const &tl);
  friend void listReplace<O>(const O old,const O _new,List<O>* &l);
  friend void listRemove<O>(const O old,List<O>* &l);
  friend int getPosition<O>(List<O>* const &l,O el);
#endif  
};

#endif
