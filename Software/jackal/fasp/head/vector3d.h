#ifndef _VECTOR
#define _VECTOR


#include <iostream.h>

#define vplus(w,u,v) { w.x=u.x+v.x; w.y=u.y+v.y; w.z=u.z+v.z; }
#define vminus(w,u,v) { w.x=u.x-v.x; w.y=u.y-v.y; w.z=u.z-v.z; }
#define dotp(u,v) ( u.x*v.x+u.y*v.y+u.z*v.z )
#define crossp(w,u,v)  { w.x=u.y*v.z-u.z*v.y; w.y=u.z*v.x-u.x*v.z; w.z=u.x*v.y-u.y*v.x; } 
#define scalarp(u,t) { u.x*=t; u.y*=t; u.z*=t; }

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);
void diag(float a[5][5],int n, float d[],float v[5][5],int *nrot);
void sorteig(float d[],float v[5][5], int n);


namespace Troll {


class Vector {
public:
  Vector() { }
  Vector(float nx,float ny,float nz) { x=nx; y=ny; z=nz; }
  float operator*(Vector v);
  Vector operator*(float t);
  Vector operator&(Vector v);
  Vector operator-(Vector v);
  Vector operator+(Vector v);
  Vector operator/(float t);
  void print();
  float norm();
  float normsq();
  void normalize();
  void normalize(float h);
  void operator+=(Vector);
  int operator==(Vector);
  double x,y,z;

} ;


class Matrix {
public:
  Matrix();
  Matrix(float,Vector);
  Matrix operator*(Matrix m);
  void operator=(Matrix m);
  Vector operator*(Vector v);
  float determinant();
  Matrix inverse();
  void print();
  Vector row[3];
} ;

Matrix RotationMatrix(Vector,float);

class Transformation {
public:
  Matrix m;
  Vector t;
} ;


extern Vector xaxis,yaxis,zaxis,zerovector;

float dihedral(Vector,Vector,Vector,Vector);
float angle(Vector,Vector,Vector);
float angle(Vector,Vector);


ostream& operator<<(ostream&,Vector&);

}


#endif


