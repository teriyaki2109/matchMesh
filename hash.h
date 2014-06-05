#ifndef _HASH_H_
#define _HASH_H_

#define NO_COLLISION (unsigned int)(-1)

class hashCell
{
public:
  int key2;
  int pos;
  int next;
};

class Hash
{
  int maxSize;
  int colSize;
  int nCollisions;
  int nElements;
  int freeCell;
  hashCell* table;
 public:
  Hash();
  Hash(int mSize);
  ~Hash();
  void del(int key1,int key2);
  int addE(int key1,int key2);
  int add(int key1,int key2);
  int get(int key1,int key2);
  int getNCollisions(){return nCollisions;};
  int getNElements();
  void clear();
};
#endif
