#ifndef GRID
#define GRID

#include <assert.h>
#include <stdio.h>
#include <string.h>
//#include "structure.h"



#ifndef min
#define min(a, b)  (((a) < (b)) ? (a) : (b)) 
#endif

#ifndef max
#define max(a, b)  (((a) > (b)) ? (a) : (b)) 
#endif



//namespace Troll {


class ObjectBox {

public:
  int ne,max;
  int *element;
} ;


template <class type> 
class Grid {

public:
  Grid(std::vector<type>& data,float sz,int l=1);
  ~Grid();
  ObjectBox& operator[](Vector);
  ObjectBox *gridvalue,defaultvalue;  
  Vector origin,extent;
  float gridsize,resolution;
  int npoints;
  int level;
  int nx, ny, nz;
  int slicesize;
    
};



template <class type>
Grid<type>::Grid(std::vector<type>& data,float sz,int le) : gridvalue(NULL), npoints(0)

{
int i,j,k,l,o;
int x_ll, x_ul, y_ll, y_ul, z_ll, z_ul;
int nelements;
int *temp;
type *element;

nelements=data.size();
if(elements==0) return;

defaultvalue.ne=0;
gridsize=sz;
element=data.begin();

level=le;

origin=data[0]->pos;
extent=data[0]->pos;

origin=origin-Vector(sz,sz,sz);
extent=extent+Vector(sz,sz,sz);


for(i=1;i<nelements;i++) {
  if(data[i]->pos.x<origin.x) origin.x=data[i]->pos.x;
  if(data[i]->pos.y<origin.y) origin.y=data[i]->pos.y;
  if(data[i]->pos.z<origin.z) origin.z=data[i]->pos.z;
  if(data[i]->pos.x>extent.x) extent.x=data[i]->pos.x;
  if(data[i]->pos.y>extent.y) extent.y=data[i]->pos.y;
  if(data[i]->pos.z>extent.z) extent.z=data[i]->pos.z;
  }

resolution=1.0/gridsize;


// Determine number of grid boxes.

origin.x-=gridsize;
origin.y-=gridsize;
origin.z-=gridsize;
extent.x+=gridsize;
extent.y+=gridsize;
extent.z+=gridsize;

nx=(int)(1+(extent.x-origin.x)/gridsize);
ny=(int)(1+(extent.y-origin.y)/gridsize);
nz=(int)(1+(extent.z-origin.z)/gridsize);


// Allocate memory for grid boxes;

npoints=nx*ny*nz+1;
gridvalue=new ObjectBox[npoints];
if(gridvalue==NULL) throw(TrollError(TE_OutOfMemory,"Grid::Grid(), too many grid boxes"));


for(i=0;i<npoints;i++) {
  gridvalue[i].ne=0;
  gridvalue[i].max=30;
  gridvalue[i].element=NULL;
  }


// Fill the boxes with objects.

for(i=0;i<nelements;i++) {
  Vector v=data[i]->pos;
  j=(int)((v.x-origin.x)*resolution);
  k=(int)((v.y-origin.y)*resolution);
  l=(int)((v.z-origin.z)*resolution);
  x_ll=max(0,j-level);
  x_ul=min(nx-1,j+level);
  y_ll=max(0,k-level);
  y_ul=min(ny-1,k+level);
  z_ll=max(0,l-level);
  z_ul=min(nz-1,l+level);
  
  for(j=x_ll;j<=x_ul; j++) for(k=y_ll;k<=y_ul;k++) for(l=z_ll;l<=z_ul;l++) {
    o=l*nx*ny + k*nx + j;
    if(gridvalue[o].ne==0) gridvalue[o].element=new int[30];
    if(gridvalue[o].ne==gridvalue[o].max) {
      gridvalue[o].max+=10;
      temp=new int[gridvalue[o].max];
      assert(temp!=NULL);
      memmove(temp,gridvalue[o].element,(gridvalue[o].max-10)*sizeof(int));
      delete gridvalue[o].element;
      gridvalue[o].element=temp;
      }
    gridvalue[o].element[gridvalue[o].ne++]=i;
    }
  }

}



template <class type>
ObjectBox& Grid<type>::operator[](Vector v)

{
int i, j, k;

i=(int)((v.x-origin.x)*resolution);
j=(int)((v.y-origin.y)*resolution);
k=(int)((v.z-origin.z)*resolution);
if(i<0 || i>=nx) return defaultvalue;
if(j<0 || j>=ny) return defaultvalue;
if(k<0 || k>=nz) return defaultvalue;

i=k*nx*ny + j*nx + i;

return gridvalue[i];

}


template <class type>
Grid<type>::~Grid()

{
int i;

if(npoints) for(i=0;i<npoints;i++) if(gridvalue[i].ne) if(gridvalue[i].element!=NULL) 
  delete gridvalue[i].element;
if(gridvalue!=NULL) delete gridvalue;

}


//} // End namespace


#endif
