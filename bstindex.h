#ifndef __BSTINDEX_H_
#define _BSTINDEX_H_

#include "geombasic.h"

const int indexNULL=-1;

struct Coord
{
  int page;
  int offset;
  Coord(void) :
    page(indexNULL),
    offset(indexNULL)
  {}

  Coord(int p,int o) :
    page(p),
    offset(o)
  {}
  bool isNull(void) const {
    return (page == indexNULL);
  }
};

class oneNode
{
 public:
  int indexPoint;
  Coord left,right;
  oneNode(void) :
    indexPoint(indexNULL),
    left(Coord()),
    right(Coord())
  {}
  void set(int index,Coord l,Coord r) {
    indexPoint = index;
    left = l;
    right = r;
  }
};

const int pageSize=16384;

class BSTIndex
{
  Point_3* point0;
  Coord rootCoord;
  oneNode** arrays;
  int lastPage;
  int lastOffset;

  inline oneNode* getNode(Coord coord) const;
  void getNewCoord(void);
  void insert_int(int index);
  int takeSmallest_int(void);
  
 public:
  BSTIndex();
  BSTIndex(Point_3* const  &_point0) :
    point0(_point0),
    rootCoord(Coord()),
    arrays(NULL),
    lastPage(-1),
    lastOffset(pageSize - 1)
  {}
  ~BSTIndex(void);

  void printFrom(Coord coord) const;
  void print(void) const;
  void insert(Point_3* p) { insert_int(p - point0); }
  Point_3* takeSmallest(void) { return &point0[takeSmallest_int()]; }
  bool isEmpty(void) const { return rootCoord.isNull(); }
};

#endif
