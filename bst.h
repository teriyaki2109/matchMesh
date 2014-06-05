#ifndef __BST_H_
#define __BST_H_

template <class O>
class BSTNode {
 public:
  O* node;
  BSTNode<O>* leftChild, * rightChild;
  
  BSTNode<O>(O* _node,BSTNode<O> * left,BSTNode<O> * right)
    {
      node=_node;
      leftChild=left;
      rightChild=right;
    };
  
  ~BSTNode<O>(void)
    {
      if(leftChild!=NULL)
	delete leftChild;
      if(rightChild!=NULL)
	delete rightChild;
    };
};

template <class O>
class BST {
  
 public:
  
  BSTNode<O> * root;
  
  BST<O>(void)
    {root=NULL;};
  
  BST<O>(BSTNode<O> * _root)
    {root=_root;};
  
  BST<O>(O* _node,BSTNode<O> * left,BSTNode<O> * right)
    {root=new BSTNode<O>(_node,left,right);};
  
  ~BST<O>(void)
    {delete root;};
};

template<class O>
O* BSTNodeTakeSmallest(BSTNode<O>* &root)
{
  O* result;
  if (root->leftChild==NULL)
    {
      result=root->node;
      delete root->leftChild;
      root=root->rightChild;
      return result;
    }
  else
    return BSTNodeTakeSmallest(root->leftChild);
};

template<class O>
O* BSTTakeSmallest(BST<O>* &tree)
{
  if(tree==NULL)
    {
      printf(">>> Attention ! BST=NULL <<<");
      return NULL;
    }
  else
    return BSTNodeTakeSmallest(tree->root);
};

template <class O>
void BSTNodeInsert(O* const &obj,BSTNode<O>* &root)
{
  if(root==NULL)
    root=new BSTNode<O>(obj,NULL,NULL);
  else if ((*obj<*(root->node))||(*obj==*(root->node)))
    BSTNodeInsert(obj,root->leftChild);
  else
    BSTNodeInsert(obj,root->rightChild);
};

template <class O>
void BSTInsert(O* const &obj,BST<O>* &tree)
{
  if (tree==NULL)
    printf(">>> Attention ! BST=NULL <<<");
  else
    BSTNodeInsert(obj,tree->root);
};

inline int max(int x,int y)
{
  if(x>y)
    return x;
  else
    return y;
}

template<class O>
int BSTNodeDepth(BSTNode<O>* const &node)
{
  if (node==NULL)
    return 0;
  else
    return 1+max(BSTNodeDepth(node->leftChild),BSTNodeDepth(node->rightChild));
}

template<class O>
int BSTNodeDepthDiff(BSTNode<O>* const &node)
{
  if(node==NULL)
    return 0;
  else 
    return (BSTNodeDepth(node->leftChild)-BSTNodeDepth(node->rightChild));
}

#endif

