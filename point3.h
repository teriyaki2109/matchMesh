#ifndef   __POINT3_H__
#define   __POINT3_H__

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
class Matrix;
class Point3
{
  friend Point3 operator/(float f,Point3 p) 
    {
      return (Point3(f/p.x,f/p.y,f/p.z));
    };
 public:  
  float x, y, z;

  Point3 (float a, float b, float c)
     :x(a), y(b), z(c) {};
  Point3 (float a)
     :x(a), y(a), z(a) {};
  Point3 ()
     :x(0), y(0), z(0) {};
  Point3 operator+(Point3 p)const
    {
      return Point3(x + p.x, y + p.y, z + p.z);
    };
  Point3 operator+(float c) const
    {
      return Point3(x + c, y + c, z + c);
    };
  Point3 operator-(Point3 p) const
    {
      return Point3(x - p.x, y - p.y, z - p.z);
    };

  Point3 operator-()const
    {
      return Point3(-x, -y, -z);
    };
  
  /*
   ** operateurs
   */
  float operator*(const Point3& p)const { return (x * p.x + y * p.y + z * p.z); };
  Point3 operator/(const Point3& p)const { return Point3(x/p.x,y/p.y,z/p.z);};
  Point3 operator*(float c)const { return (Point3(x * c, y * c, z * c)); };
  Point3 operator/(float c)const { return (Point3(x / c, y / c, z / c)); };

  void operator/=(float c) { x /= c; y /= c; z /= c; };
  void operator*=(float c) { x *= c; y *= c; z *= c; };
  void operator+=(float c) { x += c; y += c; z += c; };

  void operator+=(Point3 p) { x += p.x; y += p.y; z += p.z; };
  void operator-=(Point3 p) { x -= p.x; y -= p.y; z -= p.z; };
  void operator/=(Point3 p) { x /= p.x; y /= p.y; z /= p.z; };
  void operator*=(Point3 p) { x *= p.x; y *= p.y; z *= p.z; };
  Matrix operator^(Point3 p);
  /*
   ** methodes
   */
  Point3 gr(Point3 p) const
    {
      return Point3(x * p.x, y * p.y, z * p.z); 
    };
  Point3 vec(Point3 p) const
    {
      return (Point3(y * p.z - z * p.y,z * p.x - x * p.z,x * p.y - y * p.x)); 
    };
  float sum() const
    {
      return (x + y + z); 
    };

  float norm() const
    {
      return sqrt(x*x + y*y +z*z); 
    };
  float norm2() const
    {
      return x*x + y*y +z*z; 
    };  
  Point3& normalize() 
    { 
      float l = norm(); 
      x /= l;
      y /= l;
      z /= l;
      return *this;
    };  

  Point3 cartesian() const
    { 
      Point3 car; 
      car.x=x*sin(y)*cos(z);
      car.y=x*sin(y)*sin(z);
      car.z=x*cos(y);
      return car;
    };    
  
  Point3 spherical() const
    { 
      Point3 sph;
      sph.x=norm();
      sph.y=acos(z/sph.x);
      sph.z=atan2(y,x);
      return sph;
    };
  float maxi() const
    {
      if(x>y)
	if(x>z) return x;
	else return z;
      else
	if(y>z) return y;
	else return z;
    }
  void print() const
    {
      printf("[%f,%f,%f]\n",x,y,z);
    }
};

#endif /* __POINT3_H__ */



