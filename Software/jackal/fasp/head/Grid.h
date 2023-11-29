#ifndef _GRID
#define _GRID

#include <assert.h>
#include <stdio.h>
#include <string.h>

#ifndef _min
#define min(a, b)  (((a) < (b)) ? (a) : (b)) 
#endif

#ifndef _max
#define max(a, b)  (((a) > (b)) ? (a) : (b)) 
#endif


class Vector;
class ObjectBox {

public:
  int ne,max;
  int *element;
  ObjectBox();
  void remove(int);
  void add(int);
  int exist(int);
};

class Grid {

public:
  Grid(Atom **data,float sz,int l=1);
  ~Grid();
  ObjectBox &operator[](Vector);
  ObjectBox *gridvalue,defaultvalue;  
  Vector origin,extent;
  float gridsize,resolution;
  int npoints;
  int level;
  int nx, ny, nz;
  int slicesize;
};

#endif
