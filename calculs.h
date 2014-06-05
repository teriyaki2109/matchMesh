#include "point3.h"
#include "geombasic.h"

#define PI 3.141592653589
#define DEG2RAD( x )     ((x) * PI / 180.0)
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	                      a[k][l]=h+s*(g-h*tau);

Point_3 MultMatPoint_3(float **A, Point_3 *X);

Point3 MultMatPoint3(float **A, Point3 *X);

void MultMatVect(float **A, float *X, float *Y ,int m, int n);

float **MultMatrix(float **A, float **B);

void InverseMatrixGreville(float **A, float **B, int NbL);

void InverseMatrice3x3(float **m, float **n);

void jacobi(float **a, int n, float *d, float **v, int *nrot);

void eigsrt(float *d, float **v, int n);

int Sphere2Octaedre128(Point3 SphericalPoint);

Point_3 Octaedre128Sph(int pos);

int Sphere2Octaedre32(Point_3 SphericalPoint);
int Sphere2Octaedre32(Point3 SphericalPoint);

Point_3 Octaedre32Sph(int pos);


