#ifndef _LatSURFACE
#define _LatSURFACE


//#include "grid.h"


//namespace Troll {


class Vector;
class Structure;
class GridNew;
class Atom;
class Icosahedron;


class LatSurfv {

public:
  LatSurfv();
  ~LatSurfv();
  void ready(Atom**,int ,float);
  void assignarea(Pdb *);
  void ready(Pdb*,int ,float);
  float Surface(Vector,float);
  float Surface(float *,float);
  float Surface(float,float,float,float);
  float calcarea(Chn*);
  float calcarea(Atm*);
  float calcarea(Res*);
  float calcarea(Res *,int);
  float calcarea(Pdb *);
  void update(Atm *);
  void update(Res *);
  void update(Chn *);
  void update(Pdb *);
  void update(Res *,int);
  void clearsurfvmask(Atm*);
  void clearsurfvmask(Res *);
  void clearsurfvmask(Pdb *);
  void clearsurfvmask(Chn *);
  void clearsurfvmask(Res *,int);
  Atom** atom;
  float probe_radius;
  GridNew *gb;
  float *rp2;
  Icosahedron *sp;
};


//} // End namespace Troll


#endif
