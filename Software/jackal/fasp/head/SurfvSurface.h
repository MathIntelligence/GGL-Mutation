#ifndef SURFACE
#define SURFACE


//#include "grid.h"


//namespace Troll {


class Vector;
class Structure;
class Grid;
class Atom;
class Icosahedron;


class SurfvSurface {

public:
  SurfvSurface();
  SurfvSurface(Atom**,int ,float);
  ~SurfvSurface();
  void ready(Atom**,int ,float);
  void Build(unsigned long flag=0);
  float Surface(Vector,float);

  Atom** atom;
  float probe_radius;
  Grid *gb;
  float *rp2;
  Icosahedron *sp;

};


//} // End namespace Troll


#endif
