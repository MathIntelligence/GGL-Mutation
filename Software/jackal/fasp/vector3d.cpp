#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "vector3d.h"


using namespace Troll;




namespace Troll {

Vector xaxis(1.0,0.0,0.0);
Vector yaxis(0.0,1.0,0.0);
Vector zaxis(0.0,0.0,1.0);
Vector zerovector(0.0,0.0,0.0);

}

void Vector::operator+=(Vector v)

{

x=x+v.x;
y=y+v.y;
z=z+v.z;

}


float Vector::operator*(Vector v)

{
float u=x*v.x+y*v.y+z*v.z;
return u;
}

Vector Vector::operator*(float t)

{
Vector u;
u.x=x*t;
u.y=y*t;
u.z=z*t;
return u;
}


Vector Vector::operator&(Vector v)
{
Vector t;

t.x=y*v.z-z*v.y;
t.y=z*v.x-x*v.z;
t.z=x*v.y-y*v.x;
    
return t;
}


Vector Vector::operator-(Vector v)

{

Vector t;
t.x=x-v.x;
t.y=y-v.y;
t.z=z-v.z;
return t;
}


Vector Vector::operator+(Vector v)

{

Vector t;
t.x=x+v.x;
t.y=y+v.y;
t.z=z+v.z;
return t;

}


int Vector::operator==(Vector v)

{

if(x!=v.x) return 0;
if(y!=v.y) return 0;
if(z!=v.z) return 0;

return 1;

}


Vector Vector::operator/(float t)

{

Vector u;
u.x=x/t;
u.y=y/t;
u.z=z/t;
return u;

}

void Vector::print() { printf("x: %f  y: %f  z: %f\n",x,y,z); }

float Vector::norm() { return sqrt(x*x+y*y+z*z); }

float Vector::normsq() { return x*x+y*y+z*z; }

void Vector::normalize() { float d=sqrt(x*x+y*y+z*z); x/=d; y/=d; z/=d; }

void Vector::normalize(float h) { float d=sqrt(x*x+y*y+z*z); x=x*h/d; y=y*h/d; z=z*h/d; }



Matrix::Matrix()

{

row[0].x=1.0;
row[0].y=0.0;
row[0].z=0.0;
row[1].x=0.0;
row[1].y=1.0;
row[1].z=0.0;
row[2].x=0.0;
row[2].y=0.0;
row[2].z=1.0;

}


Matrix::Matrix(float angle,Vector axis)

{

Vector v=axis;
v.normalize();

float s=sin(angle);
float c=cos(angle);

float uxx=v.x*v.x;
float uxy=v.x*v.y;
float uxz=v.x*v.z;
float uyy=v.y*v.y;
float uyz=v.y*v.z;
float uzz=v.z*v.z;

row[0].x=uxx+c*(1.0-uxx);
row[0].y=uxy-c*uxy-v.z*s;
row[0].z=uxz-c*uxz+v.y*s;
row[1].x=uxy-c*uxy+v.z*s;
row[1].y=uyy+c*(1.0-uyy);
row[1].z=uyz-c*uyz-v.x*s;
row[2].x=uxz-c*uxz-v.y*s;
row[2].y=uyz-c*uyz+v.x*s;
row[2].z=uzz+c*(1.0-uzz);

}


Matrix Matrix::operator*(Matrix m)

{
Matrix t;
t.row[0].x=row[0].x*m.row[0].x+row[0].y*m.row[1].x+
  row[0].z*m.row[2].x;
t.row[0].y=row[0].x*m.row[0].y+row[0].y*m.row[1].y+
  row[0].z*m.row[2].y;
t.row[0].z=row[0].x*m.row[0].z+row[0].y*m.row[1].z+
  row[0].z*m.row[2].z;
t.row[1].x=row[1].x*m.row[0].x+row[1].y*m.row[1].x+
  row[1].z*m.row[2].x;
t.row[1].y=row[1].x*m.row[0].y+row[1].y*m.row[1].y+
  row[1].z*m.row[2].y;
t.row[1].z=row[1].x*m.row[0].z+row[1].y*m.row[1].z+
  row[1].z*m.row[2].z;
t.row[2].x=row[2].x*m.row[0].x+row[2].y*m.row[1].x+
  row[2].z*m.row[2].x;
t.row[2].y=row[2].x*m.row[0].y+row[2].y*m.row[1].y+
  row[2].z*m.row[2].y;
t.row[2].z=row[2].x*m.row[0].z+row[2].y*m.row[1].z+
  row[2].z*m.row[2].z;
return t;
}


void Matrix::operator=(Matrix m)

{
for(int i=0;i<3;i++)
  {
  row[i].x=m.row[i].x;
  row[i].y=m.row[i].y;
  row[i].z=m.row[i].z;
  }
}


Vector Matrix::operator*(Vector v)

{
Vector t;
t.x=row[0].x*v.x+row[0].y*v.y+row[0].z*v.z;
t.y=row[1].x*v.x+row[1].y*v.y+row[1].z*v.z;
t.z=row[2].x*v.x+row[2].y*v.y+row[2].z*v.z;
return t;
}


float Matrix::determinant()

{
float result;

result=row[0].x*(row[1].y*row[2].z-row[1].z*row[2].y);
result+=row[0].y*(row[1].z*row[2].x-row[1].x*row[2].z);
result+=row[0].z*(row[1].x*row[2].y-row[1].y*row[2].x);

return(result);

}


Matrix Matrix::inverse()

{
Matrix t;
float d;

d=determinant();

t.row[0].x=(row[1].y*row[2].z-row[1].z*row[2].y)/d;
t.row[0].y=(row[0].z*row[2].y-row[0].y*row[2].z)/d;
t.row[0].z=(row[0].y*row[1].z-row[0].z*row[1].y)/d;
t.row[1].x=(row[1].z*row[2].x-row[1].x*row[2].z)/d;
t.row[1].y=(row[0].x*row[2].z-row[0].z*row[2].x)/d;
t.row[1].z=(row[0].z*row[1].x-row[0].x*row[1].z)/d;
t.row[2].x=(row[1].x*row[2].y-row[1].y*row[2].x)/d;
t.row[2].y=(row[0].y*row[2].x-row[0].x*row[2].y)/d;
t.row[2].z=(row[0].x*row[1].y-row[0].y*row[1].x)/d;

return t;

}



void Matrix::print()

{

printf("[ %6.2f  %6.2f  %6.2f ]\n",row[0].x,row[0].y,row[0].z); 
printf("| %6.2f  %6.2f  %6.2f |\n",row[1].x,row[1].y,row[1].z); 
printf("[ %6.2f  %6.2f  %6.2f ]\n",row[2].x,row[2].y,row[2].z);
 
}


Matrix RotationMatrix(Vector v,float angle)

{
Matrix m,n;
Vector u,w;


v=v/v.norm();

m.row[0].z=v.x;
m.row[1].z=v.y;
m.row[2].z=v.z;

u.x=rand();
u.y=rand();
u.z=rand();
u=u&v;
u=u/u.norm();

m.row[0].y=u.x;
m.row[1].y=u.y;
m.row[2].y=u.z;

w=u&v;
w=w/w.norm();

m.row[0].x=w.x;
m.row[1].x=w.y;
m.row[2].x=w.z;

n.row[0].x=cos(angle);
n.row[0].y=sin(angle);
n.row[0].z=0.0;

n.row[1].x=-sin(angle);
n.row[1].y=cos(angle);
n.row[1].z=0.0;

n.row[2].x=0.0;
n.row[2].y=0.0;
n.row[2].z=1.0;

m=m.inverse()*n*m;

return m;

}



float Troll::dihedral(Vector u1,Vector u2,Vector u3,Vector u4)

{
Vector v1,v2,v3,v4;
float sign,temp;

v1=u1-u2;
v2=u2-u3;
v3=u3-u4;

v4=v2&v3;

if((v1*v4)<0) sign=1.0;
else sign=-1.0;

v4=v1&v2;
v1=v2&v3;

temp=(v4*v1)/(v4.norm()*v1.norm());
if(temp>1.0) return 0.0;
if(temp<-1.0) return PI;
temp=sign*acos(temp);

return temp;
}



float angle(Vector v1,Vector v2,Vector v3)

{
Vector u1, u2;

u1=v1-v2;
u2=v3-v2;
  
return(acos((u1*u2)/(u1.norm()*u2.norm())));

}


float angle(Vector v1,Vector v2)

{

return(acos((v1*v2)/(v1.norm()*v2.norm())));

}


ostream& Troll::operator<<(ostream& os,Vector& v)

{

os << v.x << " " << v.y << " " << v.z ;

return os;

}

istream& operator>>(istream& is,Vector& v)

{

is >> v.x;
is >> v.y;
is >> v.z;

return is;

}




void diag(float t_a[5][5],int n, float d[],float v[5][5],int *nrot)

{
int j,iq,ip,i;
float tresh,theta,tau,t,sm,s,h,g,c,*vector();
static float a[5][5],b[5], z[5];

for(i=n;i>=0;i--) for(j=n;j>=0;j--) a[i+1][j+1]=t_a[i][j];

for(ip=1;ip<=n;ip++)
  {
  for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
  v[ip][ip]=1.0;
  }

for(ip=1;ip<=n;ip++) 
  {
  b[ip]=d[ip]=a[ip][ip];
  z[ip]=0.0;
  }

*nrot=0;

for (i=1;i<=50;i++)
  {
  sm=0.0;
  for(ip=1;ip<=n-1;ip++) for(iq=ip+1;iq<=n;iq++) sm+=fabs(a[ip][iq]);
 
  if(sm == 0.0) return;

  if (i < 4) tresh=0.2*sm/(n*n);
  else tresh=0.0;

  for(ip=1;ip<=n-1;ip++)
    {
	for(iq=ip+1;iq<=n;iq++)
	  {
	  g=100.0*fabs(a[ip][iq]);
	  if(i>4 && fabs(d[ip])+g==fabs(d[ip]) && fabs(d[iq])+g == fabs(d[iq])) a[ip][iq]=0.0;
      else if (fabs(a[ip][iq]) > tresh) 
	    {
		h=d[iq]-d[ip];
		if (fabs(h)+g == fabs(h)) t=(a[ip][iq])/h;
		else 
		  {
		  theta=0.5*h/(a[ip][iq]);
		  t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
		  if (theta < 0.0) t=-t;
		  }
        c=1.0/sqrt(1+t*t);
		s=t*c;
		tau=s/(1.0+c);
		h=t*a[ip][iq];
		z[ip] -= h;
		z[iq] += h;
		d[ip] -= h;
		d[iq] += h;
		a[ip][iq]=0.0;

		for(j=1;j<=ip-1;j++) 
		  {
		  ROTATE(a,j,ip,j,iq)
		  }

		for(j=ip+1;j<=iq-1;j++)
		  {
          ROTATE(a,ip,j,j,iq)
		  }
         
		for(j=iq+1;j<=n;j++)
		  {
		  ROTATE(a,ip,j,iq,j)
		  }
		for(j=1;j<=n;j++) 
		  {
		  ROTATE(v,j,ip,j,iq)
		  }
        
		++(*nrot);
		}
	  }
	}

  for (ip=1;ip<=n;ip++)
    {
	b[ip] += z[ip];
	d[ip]=b[ip];
	z[ip]=0.0;
	}
  }

}

#undef ROTATE

void sorteig(float d[],float v[5][5], int n)
{
	int k,j,i;
	float p;

	for (i=1;i<n;i++) {
		p=d[k=i];
		for (j=i+1;j<=n;j++)
			if (d[j] >= p) p=d[k=j];
			//if (d[j] <= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=1;j<=n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}

