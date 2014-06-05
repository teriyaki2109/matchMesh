/*
 ** pointsi.h
 ** 
 ** Made by (Geoffroy Fouquier)
 ** Login   <fouquier@epita.fr>
 ** 
 ** Started on  Sat Oct 30 13:28:56 1999 Geoffroy Fouquier
 ** Last update Sat Oct 30 13:36:04 1999 Geoffroy Fouquier
 */

#ifndef   __POINTSI_H__
#define   __POINTSI_H__

class Pointsi
{
 public:
  /*
   ** attributs
   */
  short x;
  short	y;
  short	z;

  /*
   ** constructeurs
   */
  Pointsi (short a,short b,short c)
    :x(a), y(b), z(c) {};
  Pointsi ()
    :x(0), y(0), z(0) {};
  
  Pointsi (const Point3& p){
    x=(short)p.x;y=(short)p.y;z=(short)p.z;
  };
  /*
   ** operateurs
   */
  Pointsi	operator+(Pointsi p)
    {
      return Pointsi(x + p.x, y + p.y, z + p.z);
    };
  Pointsi	operator-(Pointsi p)
    {
      return Pointsi(x - p.x, y - p.y, z - p.z);
    };
  Point3 	operator+(Point3 p)
    {
      return Point3(x + p.x, y + p.y, z + p.z);
    };
  Point3	operator-(Point3 p)
    {
      return Point3(x - p.x, y - p.y, z - p.z);
    };
  float		operator*(Point3 p)
    {
      return (x * p.x + y * p.y + z * p.z);
    };
  float		operator*(Pointsi p)
    {
      return (x * p.x + y * p.y + z * p.z);
    };
  Pointsi	operator*(short c) const
    {
      return Pointsi(x * c, y * c, z * c);
    };
  void	operator*=(short c) 
    {
      x*=c;
      y*=c;
      z*=c;
    };

  Point3	operator*(float c)
    {
      return Point3(x * c, y * c, z * c);
    };
  void print()
    {
      printf("%d %d %d\n",x,y,z);
    };
  Pointsi operator/(short c)
    {
      return Pointsi(x / c, y / c, z / c);
    };
  /*
   ** methodes
   */
  Point3	normalize()
    { 
      float l = sqrt(pow(x, 2.0f) + pow(y, 2.0f) + pow(z, 2.0f));
      return Point3(x / l, y / l, z / l);
    };
  Point3	gr(Point3 p)
    {
      return Point3(x * p.x, y * p.y, z * p.z);
    };
  Point3	cp(Pointsi p)
    {
      return Point3(y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x);
    };
  
  Pointsi	cpp(Pointsi p)
    {
      return Pointsi(y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x);
    };
  float		dist(Pointsi p)
    {
      return ((float)(p.x - x) * (p.x - x) + (p.y - y) * (p.y - y) + (p.z - z)
	      * (p.z - z));
    };
};


#endif /* __POINTSI_H__ */
