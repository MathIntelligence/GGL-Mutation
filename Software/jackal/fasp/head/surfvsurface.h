#ifndef SURFACE
#define SURFACE


#include "grid.h"


namespace Troll {


class Vector;
class Structure;
class AtomGrid;
class Atom;
class Icosahedron;


class SurfvSurface {

public:
  SurfvSurface(std::vector<Atom*>,int spherelevel=2);
  ~SurfvSurface();
  void Build(unsigned long flag=0);
  float Surface(Vector,float,Grid<Atom*> *gridbox=NULL,Icosahedron *sphere=NULL,float *rplusp2=NULL);

protected:
  std::vector<Atom*> atom;
  float probe_radius;
  Grid<Atom*> *gb;
  float *rp2;
  Icosahedron *sp;

};


} // End namespace Troll


#endif
