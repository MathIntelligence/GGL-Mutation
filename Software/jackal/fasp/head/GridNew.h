#ifndef _GRIDNEW
#define _GRIDNEW

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
class ObjectBox;
class GridNew {

public:
  GridNew();
  void ready(Atom **data);
  void ready(Pdb *);
  ~GridNew();
  ObjectBox &operator[](Vector);
  void add(int);
  void initial();
  void remove(int,int);
  void puton(Atm *);
  void puton(Res *);
  void putoff();
  void putoff(Atm *);
  void putoff(Res *);
  void putoff(Res *,int);
  void putoff(Res *,int,int);
  void puton(Res *,int);
  void puton(Res *,int,int);
  void puton(Chn *);
  void puton(Pdb *);
  void putoff(Chn *);
  void putoff(Pdb *);
  void setGridNew(Atom **);
  int gettotal();

  //data
  ObjectBox *gridvalue,defaultvalue;  
  Vector origin,extent;
  float gridsize,resolution;
  float leg;
  int npoints;
  int level;
  int nx, ny, nz;
  int slicesize;
  Atom **data;
};

#endif
