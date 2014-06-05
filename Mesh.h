#ifndef _MESH_H_
#define _MESH_H_
#include "point3.h"
#include "pointsi.h"

class Edge_M
{
 public:
  unsigned int v0;
  unsigned int v1;
  bool isInvalid()const{if(v0==(unsigned int) -1) return true;else return false;};
  void setInvalid(){v0=(unsigned int) -1;v1=(unsigned int) -1;};
};

class Mesh {
  public:
  Point3 *p3d;
  int    *p3d_index;

  Point3 *t3d; //triangles de texture, il y en a NbTriangles et ils représentent les coordonnées dans la carte

  int    NbPoints;
  int    NbTriangles;
  int    NbEdges;

  Edge_M *edge;

  Mesh(char *file, int texture);
  ~Mesh();

  void createEdges();
  float facetArea(unsigned int f);
  float facetAngle(unsigned int f);
  Point3 facetNormal(unsigned int f);
  void computeMean(float *K);
  void computeMeanCurvature(float *K);
  void computeGaussianCurvature(float *K);
  void computeCurvature(float *K);

  void CurvatureDescriptor(int N, FILE *scuf);
  void CurvatureIndex(float *K);

  void ArrangeTri();
  void ComputePCA(int bbsize);
  void CanonicalFrame();
  void CEGISphere(int Ntheta, int Nphi, int V, FILE *scuf);
  void CEGI(int V, FILE *scuf);

  void CordHistogram1(int B, int V,FILE *scuf);
  void CordHistogram2(int B, int V,FILE *scuf);
  void CordHistogram3(int B, int V,FILE *scuf);
  
  void ComputeMoments(int q, int r, int s, FILE *scuf);
  void MomentsInvariants(int q, int r, int s, FILE *scuf);
  void MomentsZernike3D(int q, int r, int s, FILE *scuf);

  void ShapeDistributionD2(int N, int B, int V,FILE *scuf);
  void ShapeSpectrumDescriptor(int N, FILE *scuf);
  
  void AreaVolume(FILE *scuf);
  void HoughTransformSphere(int Ns, int Ntheta, int Nphi, float T, int V, FILE *scuf);
  void HoughTransform(int Ns, float T, int V, FILE *scuf);
  void HoughTransform32(int Ns, float T, int V, FILE *scuf);
};
#endif
